"""Multiscale image cache: build an OME-Zarr style pyramid from big TIFFs and
serve 8-bit tiles out of it.

Why a cache: the lab's Akoya exports are single-strip uncompressed TIFFs
(one channel per file, ~1 GB each, no pyramid) and the Visium microscope
images are LZW RGB TIFFs with one row per strip. Neither can be tiled at
arbitrary zoom on the fly, so `prepare.py` streams them ONCE into a chunked
zarr pyramid (chunks = 1 channel x tile x tile, zstd) and the server reads
one chunk per tile after that.

Layout of <dataset>/image.zarr (OME-NGFF 0.4 multiscales, axes c,y,x):
    0/   (C, H, W)        full resolution
    1/   (C, H/2, W/2)
    ...  until the coarsest level fits in one tile
.zattrs carries `multiscales` + `omero.channels`, so napari / vizarr can open
the cache directly too.

Tile addressing (same as the FIVE-Atlas editor / deck.gl TileLayer):
    server z in [0 .. nlevels-1];  pyramid level = (nlevels-1) - z
    z = nlevels-1  -> level 0 (full res);  z = 0 -> coarsest
"""
from __future__ import annotations

import io
import json
import re
import threading
import time
from functools import lru_cache
from pathlib import Path
from typing import Iterator, Optional

import numpy as np
import tifffile
import zarr
from numcodecs import Blosc
from PIL import Image

TILE = 1024

# Default colors assigned by channel NAME (case-insensitive substring), then
# by index for anything unnamed. DAPI/Hoechst blue is the one everyone expects.
_BY_NAME = [
    ("dapi", "#4040ff"), ("hoechst", "#4040ff"),
    ("red", "#ff0000"), ("green", "#00ff00"), ("blue", "#0000ff"), ("hoechst", "#4040ff"),
]
_PALETTE = ["#00ff00", "#ff00ff", "#ffff00", "#00ffff", "#ff8000",
            "#ff4080", "#80ff80", "#8080ff", "#ffffff", "#ff0000"]


def color_for(name: str, idx: int) -> str:
    n = (name or "").lower()
    for key, col in _BY_NAME:
        if key in n:
            return col
    return _PALETTE[idx % len(_PALETTE)]


# --------------------------------------------------------------------------
# Sources: anything that can hand back one channel, one horizontal band at a
# time, without loading the whole image.
# --------------------------------------------------------------------------
class Source:
    names: list[str]
    height: int
    width: int
    dtype: np.dtype

    @property
    def nchannels(self) -> int:
        return len(self.names)

    def read_band(self, c: int, y0: int, y1: int) -> np.ndarray:
        raise NotImplementedError

    def close(self):
        pass


class CropSource(Source):
    """A rectangle of another source, in full-res pixels of the original.

    Whole slides carry several samples with a lot of empty glass between them;
    cropping before the pyramid is built means the empty part is never tiled,
    never exported and never uploaded. Everything downstream (levels, tiles,
    the scale bar) just sees a smaller image.
    """

    def __init__(self, inner: Source, box):
        x0, y0, x1, y1 = (int(round(v)) for v in box)
        x0, y0 = max(0, x0), max(0, y0)
        x1, y1 = min(inner.width, x1), min(inner.height, y1)
        if x1 - x0 < 2 or y1 - y0 < 2:
            raise ValueError(f"crop {box} is empty against a {inner.width}x{inner.height} image")
        self._inner, self._x0, self._y0 = inner, x0, y0
        self.box = (x0, y0, x1, y1)
        self.names, self.dtype = inner.names, inner.dtype
        self.width, self.height = x1 - x0, y1 - y0
        if getattr(inner, "pixel_size_um", None):
            self.pixel_size_um = inner.pixel_size_um          # cropping does not rescale

    def read_band(self, c, y0, y1):
        band = self._inner.read_band(c, y0 + self._y0, y1 + self._y0)
        return band[:, self._x0:self._x0 + self.width]

    def close(self):
        self._inner.close()


class FlipSource(Source):
    """Mirror another source horizontally, vertically or both.

    Applied after any crop, so the flip is about the cropped image's own centre.
    Bin positions must be transformed to match - see _flip_xy in prepare.py.
    """

    def __init__(self, inner: Source, flip: str):
        self._inner = inner
        self.h = "h" in flip
        self.v = "v" in flip
        self.names, self.dtype = inner.names, inner.dtype
        self.width, self.height = inner.width, inner.height
        if getattr(inner, "pixel_size_um", None):
            self.pixel_size_um = inner.pixel_size_um

    def read_band(self, c, y0, y1):
        if self.v:
            band = self._inner.read_band(c, self.height - y1, self.height - y0)[::-1]
        else:
            band = self._inner.read_band(c, y0, y1)
        return band[:, ::-1] if self.h else band

    def close(self):
        self._inner.close()


def _page_reader(page: tifffile.TiffPage):
    """Return (get_band(y0,y1)->2D/3D array, close) for one TIFF page, picking
    the cheapest access path for how the page is laid out on disk."""
    tf_lock = threading.Lock()
    H = page.imagelength
    if page.is_memmappable:
        # uncompressed + contiguous: the Akoya exports. Slicing a memmap is a
        # sequential read of just those rows.
        mm = page.asarray(out="memmap")
        return (lambda y0, y1: np.asarray(mm[y0:y1])), (lambda: None)
    store = page.aszarr()
    z = zarr.open(store, mode="r")
    nchunks_y = int(np.ceil(H / z.chunks[0])) if z.chunks else 1
    if nchunks_y > 1:
        # strips / tiles: zarr decodes only the strips a band touches
        def get(y0, y1):
            with tf_lock:
                return np.asarray(z[y0:y1])
        return get, store.close
    # single compressed strip: no way around decoding the whole thing once
    store.close()
    full = page.asarray()
    return (lambda y0, y1: full[y0:y1]), (lambda: None)


class TiffFileSource(Source):
    """A list of single-channel (YX) or RGB (YXS) TIFF files, one per channel
    (the Akoya `NN_Marker_Sample_gray.tif` layout). RGB files become one
    channel by max over RGB unless `rgb_split` is set, in which case they
    contribute three channels (R, G, B)."""

    def __init__(self, files: list[tuple[str, Path]], rgb_split: bool = False):
        self.names = []
        self._plan = []          # (name, path, mode) mode: 'gray' | 'max' | ('rgb', k)
        self._open: dict[Path, tuple] = {}
        self.height = self.width = None
        self.dtype = None
        for name, path in files:
            # is_ome=False: a per-channel OME file (Xenium morphology_focus) is read as
            # its own pyramid instead of tifffile trying to assemble the multi-file set
            with tifffile.TiffFile(str(path), is_ome=False) as tf:
                s = tf.series[0]
                ax = s.axes
                shp = s.shape
                H, W = shp[ax.index("Y")], shp[ax.index("X")]
                if self.height is None:
                    self.height, self.width, self.dtype = H, W, np.dtype(s.dtype)
                elif (H, W) != (self.height, self.width):
                    raise ValueError(f"{path.name}: {H}x{W} does not match {self.height}x{self.width}")
                if "S" in ax and shp[ax.index("S")] >= 3:
                    if rgb_split:
                        for k, suf in enumerate(("Red", "Green", "Blue")):
                            self.names.append(f"{name} {suf}" if name else suf)
                            self._plan.append((path, ("rgb", k)))
                    else:
                        self.names.append(name)
                        self._plan.append((path, "max"))
                else:
                    self.names.append(name)
                    self._plan.append((path, "gray"))

    def _reader(self, path: Path):
        if path not in self._open:
            tf = tifffile.TiffFile(str(path), is_ome=False)
            get, close = _page_reader(tf.series[0].pages[0])
            self._open[path] = (tf, get, close)
        return self._open[path][1]

    def read_band(self, c, y0, y1):
        path, mode = self._plan[c]
        band = self._reader(path)(y0, y1)
        if mode == "gray":
            return band if band.ndim == 2 else band[..., 0]
        if mode == "max":
            return band.max(axis=-1)
        return band[..., mode[1]]

    def release(self, c):
        """Drop the handle once a channel is done (keeps SMB handles low)."""
        path, _ = self._plan[c]
        if path in self._open and all(p != path for p, _m in self._plan[c + 1:]):
            tf, _g, close = self._open.pop(path)
            close()
            tf.close()

    def close(self):
        for tf, _g, close in self._open.values():
            close()
            tf.close()
        self._open = {}


class TiffStackSource(Source):
    """One TIFF holding all channels: CYX / IYX / ZYX stacks (ImageJ, OME,
    qptiff level 0) or a plain YX / YXS image."""

    def __init__(self, path: Path, names: Optional[list[str]] = None, rgb_split: bool = True):
        self.path = Path(path)
        self._tf = tifffile.TiffFile(str(path))
        s = self._tf.series[0]
        self._series = s
        ax = s.axes
        shp = s.shape
        self.height, self.width = shp[ax.index("Y")], shp[ax.index("X")]
        self.dtype = np.dtype(s.dtype)
        self._lock = threading.Lock()
        # channel axis = whatever non-YXS axis is left (C, I, Z, Q ...)
        other = [i for i, a in enumerate(ax) if a not in "YXS"]
        self._cax = other[0] if other else None
        nC = shp[self._cax] if self._cax is not None else 1
        self._nS = shp[ax.index("S")] if "S" in ax else 1
        self._rgb_split = rgb_split and self._nS >= 3
        found = names or self._channel_names(nC)
        self.is_rgb = self._rgb_split
        if self._rgb_split:
            self.names = []
            for i in range(nC):
                base = found[i] if i < len(found) else f"ch{i}"
                for suf in ("Red", "Green", "Blue"):
                    self.names.append(f"{base} {suf}" if nC > 1 else suf)
        else:
            self.names = [found[i] if i < len(found) else f"ch{i}" for i in range(nC)]
        self._z = zarr.open(s.aszarr(level=0), mode="r") if len(s.levels) > 1 else zarr.open(s.aszarr(), mode="r")
        self._ax = ax
        self.pixel_size_um = self._pixel_size()

    def _pixel_size(self):
        """µm per pixel from the TIFF resolution tags (qptiff/OME carry them)."""
        try:
            pg = self._series.pages[0]
            xres = pg.resolution[0]
            unit = getattr(pg.resolutionunit, "value", pg.resolutionunit)
            if not xres or xres <= 1:
                return None
            if unit == 3:      # centimetre
                return 10000.0 / float(xres)
            if unit == 2:      # inch
                return 25400.0 / float(xres)
        except Exception:
            pass
        return None

    def _channel_names(self, nC) -> list[str]:
        tf = self._tf
        try:
            ij = tf.imagej_metadata or {}
            labels = ij.get("Labels")
            if labels and len(labels) >= nC:
                return [str(x) for x in labels[:nC]]
        except Exception:
            pass
        try:
            ome = tf.ome_metadata
            if ome:
                names = re.findall(r'<Channel[^>]*?\sName="([^"]*)"', ome)
                if len(names) >= nC:
                    return names[:nC]
        except Exception:
            pass
        # qptiff / PerkinElmer: per-page XML with <Name> or <Biomarker>
        names = []
        try:
            for pg in self._series.pages[:nC]:
                d = pg.description or ""
                m = re.search(r"<Biomarker>(.*?)</Biomarker>", d) or re.search(r"<Name>(.*?)</Name>", d)
                names.append(m.group(1).strip() if m else "")
        except Exception:
            names = []
        if names and any(names):
            return [n or f"ch{i}" for i, n in enumerate(names)]
        return [f"ch{i}" for i in range(nC)]

    def read_band(self, c, y0, y1):
        if self._rgb_split:
            ci, k = divmod(c, 3)
        else:
            ci, k = c, None
        idx = []
        for a in self._ax:
            if a == "Y":
                idx.append(slice(y0, y1))
            elif a == "X":
                idx.append(slice(None))
            elif a == "S":
                idx.append(k if k is not None else slice(None))
            else:
                idx.append(ci)
        with self._lock:
            band = np.asarray(self._z[tuple(idx)])
        if band.ndim == 3:            # YXS with rgb_split off -> max over S
            band = band.max(axis=-1)
        return band

    def close(self):
        try:
            self._tf.close()
        except Exception:
            pass


# --------------------------------------------------------------------------
# Builder
# --------------------------------------------------------------------------
def _downsample2(a: np.ndarray, how: str = "mean") -> np.ndarray:
    """2x2 reduce, edge-padded to even dims. Integer maths for uint8/uint16.

    `how` matters more than it looks. Fluorescence is sparse: a structure a few pixels
    wide sits in near-black background, so a 2x2 MEAN quarters an isolated bright pixel
    at every level - about 1000x over five levels - and the signal averages away into
    black when zoomed out, then reappears on zoom-in. MAX keeps the brightest pixel of
    each block, so a marker looks the same at every zoom. Brightfield and H&E are dense
    and continuous, and want the mean: max would make them look blown out and noisy.
    """
    h, w = a.shape
    if h % 2 or w % 2:
        a = np.pad(a, ((0, h % 2), (0, w % 2)), mode="edge")
    q = (a[0::2, 0::2], a[1::2, 0::2], a[0::2, 1::2], a[1::2, 1::2])
    if how == "max":
        return np.maximum(np.maximum(q[0], q[1]), np.maximum(q[2], q[3]))
    if a.dtype == np.uint8:
        s = q[0].astype(np.uint16) + q[1] + q[2] + q[3]
        return ((s + 2) >> 2).astype(np.uint8)
    s = q[0].astype(np.uint32) + q[1] + q[2] + q[3]
    return ((s + 2) >> 2).astype(a.dtype)


def _to_storage_dtype(band: np.ndarray, dtype: np.dtype) -> np.ndarray:
    if band.dtype == dtype:
        return band
    if dtype == np.uint16 and band.dtype.kind == "f":
        return np.clip(band, 0, 65535).astype(np.uint16)
    return band.astype(dtype)


def _reduce_for(names: list[str]) -> str:
    """Photos average, fluorescence takes the max - see _downsample2."""
    rgb = {"red", "green", "blue"}
    is_photo = len(names) == 3 and {n.strip().lower() for n in names} == rgb
    return "mean" if is_photo else "max"


def build_image(src: Source, out_dir: Path, tile: int = TILE, colors: Optional[list[str]] = None,
                log=print) -> dict:
    """Stream `src` into <out_dir>/image.zarr and return the channel manifest."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    zpath = out_dir / "image.zarr"
    dtype = np.dtype(src.dtype)
    if dtype not in (np.dtype(np.uint8), np.dtype(np.uint16)):
        dtype = np.dtype(np.uint16)
    C, H, W = src.nchannels, src.height, src.width
    reduce_how = _reduce_for(list(src.names))
    log(f"[image] pyramid reduce: {reduce_how}")
    comp = Blosc(cname="zstd", clevel=3, shuffle=Blosc.BITSHUFFLE)
    g = zarr.open_group(str(zpath), mode="w")

    # level shapes: halve until the coarsest level fits in one tile
    shapes = [(H, W)]
    while max(shapes[-1]) > tile:
        h, w = shapes[-1]
        shapes.append(((h + 1) // 2, (w + 1) // 2))
    arrays = []
    for L, (h, w) in enumerate(shapes):
        arrays.append(g.create_dataset(str(L), shape=(C, h, w), chunks=(1, tile, tile), dtype=dtype,
                                       compressor=comp, dimension_separator="/", overwrite=True))
    log(f"[image] {C} channels, {W}x{H} {dtype.name}, {len(shapes)} levels -> {zpath}")

    t_all = time.time()
    for c in range(C):
        t0 = time.time()
        # level 0: band by band (tile-aligned rows => each chunk written once)
        a0 = arrays[0]
        for y0 in range(0, H, tile):
            y1 = min(y0 + tile, H)
            a0[c, y0:y1, :] = _to_storage_dtype(src.read_band(c, y0, y1), dtype)
        if hasattr(src, "release"):
            src.release(c)
        # coarser levels from the level just written, 2*tile rows at a time
        for L in range(1, len(shapes)):
            prev, cur = arrays[L - 1], arrays[L]
            ph = shapes[L - 1][0]
            for y0 in range(0, ph, 2 * tile):
                y1 = min(y0 + 2 * tile, ph)
                cur[c, y0 // 2:(y1 + 1) // 2, :] = _downsample2(np.asarray(prev[c, y0:y1, :]), reduce_how)
        log(f"[image]   ch {c:2d} {src.names[c]:<20s} {time.time() - t0:6.1f}s")
    log(f"[image] done in {time.time() - t_all:.0f}s")

    # per-channel default display window from a mid-res level (<= 4096 px)
    stats_level = next((i for i, (h, w) in enumerate(shapes) if max(h, w) <= 4096), len(shapes) - 1)
    chans = []
    for c in range(C):
        data = np.asarray(arrays[stats_level][c])
        nz = data[data > 0]
        dmax = int(data.max()) if data.size else 0
        if nz.size:
            lo, hi = float(np.percentile(nz, 1.0)), float(np.percentile(nz, 99.8))
        else:
            lo, hi = 0.0, 1.0
        if hi <= lo:
            hi = lo + 1.0
        name = src.names[c]
        col = (colors[c] if colors and c < len(colors) and colors[c] else color_for(name, c))
        if getattr(src, "is_rgb", False):
            # a color photo (H&E / brightfield): show it as-is, full range
            lo, hi = 0.0, float(dmax or 255)
        # `bits`: raw -> 8-bit mapping the server applies (identity for uint8).
        # `window`: the slider defaults the client starts with, in 8-bit units.
        if dtype == np.uint8:
            bits = [0, 255]
            window = [lo, hi]
        else:
            # 16-bit: map raw 0..p99.99 (not 0..dataMax) onto the 8-bit tile so a
            # dim channel keeps its low end; the rare pixels above saturate.
            top = float(np.percentile(nz, 99.99)) if nz.size else float(max(dmax, 1))
            top = float(max(top, hi, 1))
            bits = [0, top]
            window = [lo * 255.0 / top, hi * 255.0 / top]
        chans.append({"index": c, "name": name, "color": col, "bits": bits,
                      "window": [round(window[0], 2), round(window[1], 2)], "dataMax": dmax})

    # OME-NGFF metadata so other tools can read the cache
    g.attrs["multiscales"] = [{
        "version": "0.4", "name": out_dir.name,
        "axes": [{"name": "c", "type": "channel"}, {"name": "y", "type": "space", "unit": "pixel"},
                 {"name": "x", "type": "space", "unit": "pixel"}],
        "datasets": [{"path": str(L), "coordinateTransformations": [{"type": "scale", "scale": [1, 2 ** L, 2 ** L]}]}
                     for L in range(len(shapes))],
        "type": "mean 2x2",
    }]
    g.attrs["omero"] = {"channels": [
        {"label": ch["name"], "color": ch["color"].lstrip("#"), "active": i < 4,
         "window": {"min": 0, "max": ch["bits"][1], "start": ch["window"][0] * ch["bits"][1] / 255,
                    "end": ch["window"][1] * ch["bits"][1] / 255}}
        for i, ch in enumerate(chans)]}

    return {"width": W, "height": H, "levels": len(shapes), "tileSize": tile, "dtype": dtype.name,
            "levelShapes": [[h, w] for h, w in shapes], "channels": chans}


# --------------------------------------------------------------------------
# Reader (server side)
# --------------------------------------------------------------------------
class PyramidReader:
    def __init__(self, ds_dir: Path, manifest: dict):
        self.dir = Path(ds_dir)
        self.m = manifest
        img = manifest["image"]
        self.W, self.H = img["width"], img["height"]
        self.tile = img["tileSize"]
        self.nlevels = img["levels"]
        self.shapes = [tuple(s) for s in img["levelShapes"]]
        self.channels = img["channels"]
        self._g = zarr.open_group(str(self.dir / "image.zarr"), mode="r")
        self._arrays = [self._g[str(L)] for L in range(self.nlevels)]
        self._lock = threading.Lock()

    def tile_u8(self, c: int, level: int, tx: int, ty: int, lo: Optional[float] = None,
                hi: Optional[float] = None) -> Optional[np.ndarray]:
        if not (0 <= level < self.nlevels) or not (0 <= c < len(self.channels)):
            return None
        H, W = self.shapes[level]
        T = self.tile
        y0, x0 = ty * T, tx * T
        if y0 >= H or x0 >= W or y0 < 0 or x0 < 0:
            return None
        y1, x1 = min(y0 + T, H), min(x0 + T, W)
        a = self._arrays[level]
        # zarr DirectoryStore reads are safe to run concurrently
        win = np.asarray(a[c, y0:y1, x0:x1])
        bits = self.channels[c]["bits"]
        blo = float(bits[0] if lo is None else lo)
        bhi = float(bits[1] if hi is None else hi)
        if win.dtype == np.uint8 and blo == 0 and bhi == 255:
            u8 = win
        else:
            if bhi <= blo:
                bhi = blo + 1
            f = (win.astype(np.float32) - blo) * (255.0 / (bhi - blo))
            u8 = np.clip(f, 0, 255).astype(np.uint8)
        if u8.shape != (T, T):
            out = np.zeros((T, T), np.uint8)
            out[: u8.shape[0], : u8.shape[1]] = u8
            u8 = out
        return u8

    def tile_bytes(self, c, level, tx, ty, fmt="png", lo=None, hi=None) -> Optional[bytes]:
        u8 = self.tile_u8(c, level, tx, ty, lo, hi)
        if u8 is None:
            return None
        buf = io.BytesIO()
        img = Image.fromarray(u8, mode="L")
        if fmt == "jpg":
            img.save(buf, format="JPEG", quality=88)
        else:
            img.save(buf, format="PNG", compress_level=1)
        return buf.getvalue()

    def thumbnail(self, c: int, max_side: int = 512) -> bytes:
        a = self._arrays[self.nlevels - 1]
        data = np.asarray(a[c])
        bits = self.channels[c]["bits"]
        f = (data.astype(np.float32) - bits[0]) * (255.0 / max(bits[1] - bits[0], 1))
        img = Image.fromarray(np.clip(f, 0, 255).astype(np.uint8), mode="L")
        img.thumbnail((max_side, max_side))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
