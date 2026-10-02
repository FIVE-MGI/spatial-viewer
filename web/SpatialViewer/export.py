"""Export a prepared dataset as a plain static file tree (no server needed).

    <dest>/<id>/manifest.json
    <dest>/<id>/tiles/<c>/<z>/<x>/<y>.png       one channel, 8-bit grey (or .jpg)
    <dest>/<id>/expr/<key>/genes.json           [{n, t}]
    <dest>/<id>/expr/<key>/positions.bin        float32 [x,y]*
    <dest>/<id>/expr/<key>/gene/<NAME>.bin      float32 [x,y,v]*   (NAME percent-encoded like encodeURIComponent)
    <dest>/<id>/overlays/<file>.geojson

Upload the <id> folder to any static host (Azure blob, S3, a web server) and
point a `<id>.link.txt` in web/SpatialViewer at its URL - the same pattern the
Nerve page uses for its CSVs. Re-running skips files that already exist, so an
interrupted export can be resumed.
"""
from __future__ import annotations

import json
import math
import shutil
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import quote

import numpy as np

from expr import ExprReader
from pyramid import PyramidReader


def _encode(name: str) -> str:
    # identical to JavaScript encodeURIComponent
    return quote(name, safe="-_.!~*'()")


def verify_tiles(r: PyramidReader, out: Path, ext: str, log=print) -> list[int]:
    """Re-derive one tile per channel and compare it with what is on disk.

    s49_nonpainful_v2 shipped with twelve of its thirty-two channels wrong - eight
    written as solid black, four saturated to white - while the pyramid they came
    from was perfectly intact. Nothing noticed until someone opened the dataset and
    saw markers missing, days later. The tile loop has no way to tell a blank tile
    from a blank marker, so the only honest check is to recompute a tile and look.

    One tile per channel at the coarsest level: that is a single small read each,
    it covers the whole slide, and it catches any channel that came out uniformly
    wrong - which is the shape every failure of this kind has taken.
    """
    from PIL import Image

    bad = []
    z = 0                                            # coarsest level is z 0 on disk
    for c in range(len(r.channels)):
        want = r.tile_u8(c, r.nlevels - 1, 0, 0)
        p = out / "tiles" / str(c) / str(z) / "0" / f"0.{ext}"
        if want is None:
            continue
        if not p.exists():
            bad.append(c); continue
        got = np.asarray(Image.open(p).convert("L"), dtype=np.float32)
        w = float(want.mean())
        # jpeg moves the mean of a near-black tile by a grey level or so; a channel
        # that is blank-when-it-should-not-be, or saturated, is off by far more
        if abs(w - float(got.mean())) > max(1.5, 0.15 * w):
            bad.append(c)
    if bad:
        names = ", ".join(f"{c} ({r.channels[c]['name']})" for c in bad)
        log(f"[export] *** {len(bad)} channel(s) do not match the pyramid: {names}")
        log(f"[export] *** delete {out / 'tiles'} for those channels and re-run - "
            f"an export skips tiles that already exist, so re-running alone will not fix it")
    else:
        log(f"[export] verified {len(r.channels)} channel(s) against the pyramid")
    return bad


def export_static(ds_dir: Path, dest: Path, fmt: str = "png", workers: int = 8, log=print,
                  verify: bool = True) -> Path:
    ds_dir = Path(ds_dir)
    m = json.load(open(ds_dir / "manifest.json", encoding="utf-8"))
    out = Path(dest) / m["id"]
    out.mkdir(parents=True, exist_ok=True)
    r = PyramidReader(ds_dir, m)
    T = r.tile
    ext = "jpg" if fmt == "jpg" else "png"

    # ---- tiles -------------------------------------------------------------
    jobs = []
    for L, (H, W) in enumerate(r.shapes):
        z = r.nlevels - 1 - L
        nx, ny = math.ceil(W / T), math.ceil(H / T)
        for c in range(len(r.channels)):
            for ty in range(ny):
                for tx in range(nx):
                    p = out / "tiles" / str(c) / str(z) / str(tx) / f"{ty}.{ext}"
                    if not p.exists():
                        jobs.append((c, L, tx, ty, p))
    log(f"[export] {m['id']}: {len(jobs):,} tiles to write ({len(r.channels)} ch x {r.nlevels} levels) -> {out}")
    t0 = time.time()
    done = [0]

    def work(job):
        c, L, tx, ty, p = job
        data = r.tile_bytes(c, L, tx, ty, fmt=ext)
        if data is None:
            return
        p.parent.mkdir(parents=True, exist_ok=True)
        tmp = p.with_suffix(p.suffix + ".part")
        tmp.write_bytes(data)
        tmp.replace(p)
        done[0] += 1
        if done[0] % 2000 == 0:
            log(f"[export]   {done[0]:,}/{len(jobs):,} tiles  {time.time() - t0:.0f}s")

    with ThreadPoolExecutor(max_workers=workers) as ex:
        list(ex.map(work, jobs))
    log(f"[export] tiles done in {time.time() - t0:.0f}s")
    if verify:
        verify_tiles(r, out, ext, log=log)

    # ---- expression --------------------------------------------------------
    for layer in m.get("expr", []):
        key = layer["key"]
        er = ExprReader(ds_dir / "expr" / f"{key}.zarr")
        d = out / "expr" / key
        (d / "gene").mkdir(parents=True, exist_ok=True)
        (d / "genes.json").write_text(json.dumps(er.gene_list(), separators=(",", ":")), encoding="utf-8")
        (d / "positions.bin").write_bytes(np.ascontiguousarray(er.positions(), dtype="<f4").tobytes())
        # slice in memory: reading each gene through zarr decompresses whole chunks
        indices = np.asarray(er.g["indices"])
        data = np.asarray(er.g["data"], dtype=np.float32)
        indptr = er.indptr
        xy = er.xy
        n = 0
        for gi, name in enumerate(er.genes):
            p = d / "gene" / f"{_encode(name)}.bin"
            if p.exists():
                continue
            a, b = int(indptr[gi]), int(indptr[gi + 1])
            if b <= a:
                continue
            ids = indices[a:b]
            xyv = np.empty((b - a, 3), np.float32)
            xyv[:, :2] = xy[ids]
            xyv[:, 2] = data[a:b]
            xyv = xyv[~np.isnan(xyv[:, 0])]
            if len(xyv) == 0:
                continue
            p.write_bytes(np.ascontiguousarray(xyv, dtype="<f4").tobytes())
            n += 1
        log(f"[export] expr {key}: {n:,} gene files")

    # ---- overlays + manifest -----------------------------------------------
    for ov in m.get("overlays", []):
        src = ds_dir / "overlays" / ov["file"]
        if src.exists():
            (out / "overlays").mkdir(exist_ok=True)
            shutil.copyfile(src, out / "overlays" / ov["file"])
    sm = dict(m)
    sm["static"] = True
    sm["tileFormat"] = ext
    sm.pop("source", None)          # local paths are nobody else's business
    (out / "manifest.json").write_text(json.dumps(sm, indent=1), encoding="utf-8")
    log(f"[export] done -> {out}")
    return out
