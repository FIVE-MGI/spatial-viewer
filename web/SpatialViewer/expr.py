"""Visium / Visium HD expression -> per-gene point layers.

Space Ranger writes the count matrix CSC by barcode (matrix/indptr runs over
barcodes), so pulling ONE gene means scanning every non-zero. `build_expr`
transposes it once into a gene-major CSR inside the cache so the server can
answer "all bins where gene X > 0" with two array slices.

Cache layout  <dataset>/expr/<key>.zarr
    xy        (B, 2) float32   bin/spot centre in full-res image pixels (x, y)
    indptr    (G+1,) int64     CSR over genes
    indices   (nnz,) int32     -> row of xy
    data      (nnz,) float32   UMI count
    totals    (G,)   float32   per-gene sum (for ranking the picker)
    attrs: genes [..], binSizeUm, spotDiameterPx, micronsPerPixel, ...

Positions are whatever frame Space Ranger's `pxl_*_in_fullres` uses, i.e. the
pixel grid of the image given to `spaceranger count --image`. The pyramid must
be built from THAT image (or `--image-scale` must map between the two).
"""
from __future__ import annotations

import csv
import json
import re
import threading
import time
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import zarr
from numcodecs import Blosc


def find_visium_layers(outs: Path) -> list[dict]:
    """Every count layer under a Space Ranger `outs` folder.
    Visium HD: outs/binned_outputs/square_XXXum/;  classic Visium: outs itself."""
    outs = Path(outs)
    layers = []
    binned = outs / "binned_outputs"
    if binned.is_dir():
        for d in sorted(binned.iterdir()):
            m = re.match(r"square_(\d+)um$", d.name)
            if not m:
                continue
            h5 = d / "filtered_feature_bc_matrix.h5"
            pos = d / "spatial" / "tissue_positions.parquet"
            sf = d / "spatial" / "scalefactors_json.json"
            if h5.exists() and pos.exists():
                layers.append({"key": f"{int(m.group(1)):03d}um", "label": f"{int(m.group(1))} µm bins",
                               "h5": h5, "positions": pos, "scalefactors": sf if sf.exists() else None,
                               "binSizeUm": float(m.group(1))})
    h5 = outs / "filtered_feature_bc_matrix.h5"
    if h5.exists():
        sp = outs / "spatial"
        pos = next((p for p in (sp / "tissue_positions.parquet", sp / "tissue_positions.csv",
                                sp / "tissue_positions_list.csv") if p.exists()), None)
        sf = sp / "scalefactors_json.json"
        if pos is not None:
            layers.append({"key": "spots", "label": "Visium spots (55 µm)", "h5": h5, "positions": pos,
                           "scalefactors": sf if sf.exists() else None, "binSizeUm": None})
    return layers


def _read_positions(path: Path) -> dict[str, tuple[float, float]]:
    """barcode -> (x, y) in full-res pixels."""
    path = Path(path)
    if path.suffix == ".parquet":
        import pyarrow.parquet as pq
        t = pq.read_table(str(path), columns=["barcode", "pxl_col_in_fullres", "pxl_row_in_fullres"])
        bc = t.column("barcode").to_pylist()
        x = t.column("pxl_col_in_fullres").to_numpy().astype(np.float64)
        y = t.column("pxl_row_in_fullres").to_numpy().astype(np.float64)
        return dict(zip(bc, zip(x.tolist(), y.tolist())))
    out = {}
    with open(path, newline="") as f:
        rd = csv.reader(f)
        first = next(rd)
        if first and first[0].lower() == "barcode":       # header (spaceranger >= 2.0)
            cols = {n: i for i, n in enumerate(first)}
            ix, iy = cols["pxl_col_in_fullres"], cols["pxl_row_in_fullres"]
        else:                                             # tissue_positions_list.csv, no header
            ix, iy = 5, 4
            out[first[0]] = (float(first[ix]), float(first[iy]))
        for row in rd:
            out[row[0]] = (float(row[ix]), float(row[iy]))
    return out


def _move_xy(xy, crop_box, flip, log, image_size=None):
    """Put bin positions through the same geometry as the image.

    Order matters and mirrors prepare.py: crop first, then flip about the cropped
    image's own centre. Bins outside the crop become NaN, which is how this file
    already marks a barcode with no position, and the viewer skips them.
    """
    if crop_box:
        x0, y0, x1, y1 = crop_box
        outside = (xy[:, 0] < x0) | (xy[:, 0] > x1) | (xy[:, 1] < y0) | (xy[:, 1] > y1)
        xy[:, 0] -= x0
        xy[:, 1] -= y0
        xy[outside] = np.nan
        w, h = x1 - x0, y1 - y0
        log(f"[expr]   crop: {int((~outside).sum()):,} of {len(xy):,} bins kept")
    else:
        w = h = None
    if flip:
        if w is None:
            # Must be the IMAGE size, never the extent of the bins. The capture area is a
            # 6.5 mm square and a section can be taller than it, so the bins often cover
            # only part of the image; mirroring about their own extent shifts every bin.
            if not image_size:
                raise ValueError("flip without a crop needs image_size")
            w, h = image_size
        if "h" in flip:
            xy[:, 0] = w - xy[:, 0]
        if "v" in flip:
            xy[:, 1] = h - xy[:, 1]
        log(f"[expr]   flip {flip} applied to bin positions")
    return xy


def build_expr(layer: dict, out_dir: Path, image_scale: float = 1.0, log=print,
               crop_box=None, flip: str = "", image_size=None) -> dict:
    """Transpose one count layer into the cache. Returns its manifest entry."""
    t0 = time.time()
    from scipy import sparse
    key = layer["key"]
    with h5py.File(layer["h5"], "r") as f:
        m = f["matrix"]
        G, B = (int(v) for v in m["shape"][()])
        data = m["data"][()]
        indices = m["indices"][()]
        indptr = m["indptr"][()]
        barcodes = [b.decode() for b in m["barcodes"][()]]
        names = [n.decode() for n in m["features/name"][()]]
        ftype = [t.decode() for t in m["features/feature_type"][()]]
    keep = np.array([t == "Gene Expression" for t in ftype])
    csc = sparse.csc_matrix((data.astype(np.float32), indices, indptr), shape=(G, B))
    if not keep.all():
        csc = csc[keep]
        names = [n for n, k in zip(names, keep) if k]
    csr = csc.tocsr()
    csr.sort_indices()
    del csc
    log(f"[expr] {key}: {csr.shape[0]} genes x {csr.shape[1]} bins, nnz {csr.nnz:,} ({time.time() - t0:.1f}s)")

    pos = _read_positions(layer["positions"])
    xy = np.full((B, 2), np.nan, np.float32)
    miss = 0
    for i, bc in enumerate(barcodes):
        p = pos.get(bc)
        if p is None:
            miss += 1
            continue
        xy[i, 0] = p[0] * image_scale
        xy[i, 1] = p[1] * image_scale
    if miss:
        log(f"[expr]   WARNING {miss} barcodes have no position")

    xy = _move_xy(xy, crop_box, flip, log, image_size)

    sf = {}
    if layer.get("scalefactors"):
        with open(layer["scalefactors"]) as f:
            sf = json.load(f)
    spot_px = float(sf.get("spot_diameter_fullres", 0) or 0) * image_scale
    if not spot_px and layer.get("binSizeUm") and sf.get("microns_per_pixel"):
        spot_px = layer["binSizeUm"] / float(sf["microns_per_pixel"]) * image_scale

    zdir = Path(out_dir) / "expr" / f"{key}.zarr"
    g = zarr.open_group(str(zdir), mode="w")
    comp = Blosc(cname="zstd", clevel=3, shuffle=Blosc.SHUFFLE)
    g.create_dataset("xy", data=xy, chunks=(1 << 20, 2), compressor=comp, overwrite=True)
    g.create_dataset("indptr", data=csr.indptr.astype(np.int64), compressor=comp, overwrite=True)
    g.create_dataset("indices", data=csr.indices.astype(np.int32), chunks=(1 << 20,), compressor=comp, overwrite=True)
    g.create_dataset("data", data=csr.data.astype(np.float32), chunks=(1 << 20,), compressor=comp, overwrite=True)
    totals = np.asarray(csr.sum(axis=1)).ravel().astype(np.float32)
    g.create_dataset("totals", data=totals, compressor=comp, overwrite=True)
    g.attrs["genes"] = names
    meta = {"key": key, "label": layer["label"], "nBins": int(B), "nGenes": int(len(names)),
            "binSizeUm": layer.get("binSizeUm"), "spotDiameterPx": spot_px,
            "micronsPerPixel": (float(sf["microns_per_pixel"]) / image_scale) if sf.get("microns_per_pixel") else None,
            "imageScale": image_scale}
    g.attrs.update(meta)
    log(f"[expr] {key}: written ({time.time() - t0:.1f}s)")
    return meta


class ExprReader:
    def __init__(self, zdir: Path):
        self.g = zarr.open_group(str(zdir), mode="r")
        self.genes: list[str] = list(self.g.attrs["genes"])
        self.gidx = {n: i for i, n in enumerate(self.genes)}
        self.meta = {k: self.g.attrs.get(k) for k in
                     ("key", "label", "nBins", "nGenes", "binSizeUm", "spotDiameterPx", "micronsPerPixel")}
        self.xy = np.asarray(self.g["xy"])              # small enough to keep (8 B / bin)
        self.indptr = np.asarray(self.g["indptr"])
        self.totals = np.asarray(self.g["totals"])
        self._lock = threading.Lock()

    def gene_list(self) -> list[dict]:
        return [{"n": n, "t": float(t)} for n, t in zip(self.genes, self.totals)]

    def gene_xyv(self, name: str) -> Optional[np.ndarray]:
        i = self.gidx.get(name)
        if i is None:
            return None
        s, e = int(self.indptr[i]), int(self.indptr[i + 1])
        if e <= s:
            return np.zeros((0, 3), np.float32)
        with self._lock:
            ids = np.asarray(self.g["indices"][s:e])
            vals = np.asarray(self.g["data"][s:e])
        out = np.empty((e - s, 3), np.float32)
        out[:, :2] = self.xy[ids]
        out[:, 2] = vals
        return out[~np.isnan(out[:, 0])]

    def positions(self) -> np.ndarray:
        return self.xy[~np.isnan(self.xy[:, 0])]
