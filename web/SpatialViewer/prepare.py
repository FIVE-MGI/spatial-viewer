r"""Build the viewer cache for one dataset.

    python prepare.py akoya  --src <folder of NN_Marker_Sample_gray.tif> --id sample1 --out <cache folder>
    python prepare.py tiff   --src <stack.tif|qptiff|ome.tif>            --id x        --out <export root>
    python prepare.py visium --image <hires.tif> --outs <spaceranger outs> --id sample1_visiumhd --out <export root>

Common options: --label, --pixel-size (µm/px), --tile, --geojson "Label=path" (repeatable),
--rgb-split (keep RGB files as three channels instead of max-projecting).
The result is <out>/<id>/ with image.zarr, expr/, overlays/ and manifest.json;
the server just scans <out> for manifest.json files.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

import pyramid
from expr import build_expr, find_visium_layers


def log(*a):
    print(*a, flush=True)


def akoya_files(folder: Path):
    """Marker files in an Akoya export folder, ordered by their leading index.
    `NN_Marker_<Sample>_gray.tif` is a plain grey channel; `_blue` etc. are
    RGB renderings and only used when no grey twin exists."""
    folder = Path(folder)
    tifs = [p for p in folder.iterdir() if p.suffix.lower() in (".tif", ".tiff")]
    sample = folder.name
    entries = {}
    for p in tifs:
        stem = p.stem
        m = re.match(r"^(.*)_(gray|grey|blue|red|green|rgb)$", stem, re.I)
        variant = m.group(2).lower() if m else "gray"
        base = m.group(1) if m else stem
        if base.endswith("_" + sample):
            base = base[: -len(sample) - 1]
        m2 = re.match(r"^(\d+)[_ ]+(.*)$", base)
        order = int(m2.group(1)) if m2 else 10 ** 6
        name = (m2.group(2) if m2 else base).strip() or stem
        cur = entries.get(name)
        pref = 0 if variant in ("gray", "grey") else 1
        if cur is None or pref < cur[0]:
            entries[name] = (pref, order, p)
    return [(n, e[2]) for n, e in sorted(entries.items(), key=lambda kv: (kv[1][1], kv[0]))]


def parse_geojson_args(items, ds_dir: Path):
    out = []
    for i, it in enumerate(items or []):
        if "=" in it:
            label, path = it.split("=", 1)
        else:
            label, path = Path(it).stem, it
        src = Path(path)
        if not src.exists():
            log(f"[overlay] missing: {src}")
            continue
        dst = ds_dir / "overlays" / f"{i:02d}_{re.sub(r'[^A-Za-z0-9_.-]+', '_', label)}.geojson"
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        out.append({"key": dst.stem, "label": label, "file": dst.name})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="kind", required=True)

    def common(p):
        p.add_argument("--id", required=True, help="dataset id (folder name in the cache, URL-safe)")
        p.add_argument("--out", required=True, help="cache root folder")
        p.add_argument("--label")
        p.add_argument("--pixel-size", type=float, help="µm per full-res pixel (scale bar)")
        p.add_argument("--tile", type=int, default=pyramid.TILE)
        p.add_argument("--geojson", action="append", help='"Label=path.geojson" overlay in image pixel coords')
        p.add_argument("--rgb-split", action="store_true", help="RGB files -> three channels")

    a = sub.add_parser("akoya", help="folder of per-marker TIFFs")
    a.add_argument("--src", required=True)
    a.add_argument("--flip", help="h, v or hv - mirror the image")
    a.add_argument("--crop", help="x0,y0,x1,y1 in full-res pixels: build only that rectangle. Run once per region with a different --id to split a whole slide into samples.")
    common(a)
    t = sub.add_parser("tiff", help="one multi-channel TIFF / qptiff / OME-TIFF")
    t.add_argument("--src", required=True)
    t.add_argument("--channels", help="comma-separated channel names (override)")
    t.add_argument("--flip", help="h, v or hv - mirror the image")
    t.add_argument("--crop", help="x0,y0,x1,y1 in full-res pixels: build only that rectangle. Run once per region with a different --id to split a whole slide into samples.")
    common(t)
    e = sub.add_parser("export", help="write a prepared dataset as static files for blob/web hosting (no server)")
    e.add_argument("--id", required=True)
    e.add_argument("--out", required=True, help="cache root that prepare wrote into")
    e.add_argument("--dest", required=True, help="folder to write <dest>/<id>/ into (then upload it)")
    e.add_argument("--fmt", choices=["png", "jpg"], default="png")
    e.add_argument("--workers", type=int, default=8)
    r = sub.add_parser("rnaseq", help="single-cell / nucleus RNA-seq -> static UMAP dataset (written straight to --dest)")
    r.add_argument("--id", required=True)
    r.add_argument("--label")
    r.add_argument("--dest", required=True, help="static export root, e.g. D:\\SpatialVizStatic")
    r.add_argument("--h5ad", help="AnnData file with an embedding in obsm")
    r.add_argument("--embedding", default="X_umap")
    r.add_argument("--mtx", help="10x matrix folder (matrix.mtx, features.tsv, barcodes.tsv)")
    r.add_argument("--umap", help="CSV: barcode, UMAP1, UMAP2 (with --mtx)")
    r.add_argument("--obs", action="append", default=[], help="CSV: barcode, column... (repeatable, with --mtx)")
    r.add_argument("--tmp", help="local folder for the temporary gene-major arrays of big files (needs ~8 bytes per non-zero)")
    i = sub.add_parser("inspect", help="print channels / levels / dtype of a TIFF or qptiff without building")
    i.add_argument("--src", required=True)
    x = sub.add_parser("xenium", help="Xenium output folder: morphology_focus stains (+ Custom_annotation.geojson if present)")
    x.add_argument("--src", required=True)
    common(x)
    v = sub.add_parser("visium", help="Space Ranger outs + the full-res image it was run on")
    v.add_argument("--image", required=True, help="the --image TIFF given to spaceranger (RGB or grey)")
    v.add_argument("--outs", required=True)
    v.add_argument("--bins", help="comma list of layer keys to keep, e.g. 008um,016um (default all)")
    v.add_argument("--crop", help="x0,y0,x1,y1 in full-res pixels; bin positions move with it")
    v.add_argument("--flip", help="h, v or hv - mirror the image and the bins together")
    v.add_argument("--image-scale", type=float, default=1.0,
                   help="multiply spaceranger positions by this if --image differs in scale")
    common(v)
    args = ap.parse_args(argv)

    if args.kind == "rnaseq":
        import rnaseq
        if args.h5ad and Path(args.h5ad).stat().st_size > 2e9:          # large files: stream, never load
            rnaseq.export_h5ad_streaming(Path(args.h5ad), Path(args.dest), args.id, args.label or args.id,
                                         embedding=args.embedding, tmpdir=args.tmp, log=log)
            return
        if args.h5ad:
            X, genes, emb, obs = rnaseq.read_h5ad(Path(args.h5ad), args.embedding)
        elif args.mtx and args.umap:
            X, genes, emb, obs = rnaseq.read_10x(Path(args.mtx), Path(args.umap), [Path(p) for p in args.obs])
        else:
            sys.exit("rnaseq needs --h5ad, or --mtx together with --umap")
        rnaseq.export_rnaseq(X, genes, emb, obs, Path(args.dest), args.id, args.label or args.id, log=log)
        return

    if args.kind == "export":
        from export import export_static
        export_static(Path(args.out) / args.id, Path(args.dest), fmt=args.fmt, workers=args.workers, log=log)
        return

    if args.kind == "inspect":
        import tifffile
        with tifffile.TiffFile(args.src) as tf:
            s0 = tf.series[0]
            print(f"file: {args.src}")
            print(f"  bigtiff={tf.is_bigtiff} qpi={getattr(tf, 'is_qpi', None)} ome={tf.is_ome} series={len(tf.series)}")
            print(f"  series[0]: axes={s0.axes} shape={s0.shape} dtype={s0.dtype} levels={len(s0.levels)}")
            for L, lv in enumerate(s0.levels[:10]):
                print(f"    level {L}: {lv.shape}")
            pg = s0.pages[0]
            print(f"  page0: tiled={pg.is_tiled} tile={pg.tilewidth}x{pg.tilelength} compression={pg.compression.name}")
        src = pyramid.TiffStackSource(Path(args.src), rgb_split=True)
        print(f"  channels ({src.nchannels}): " + ", ".join(src.names))
        src.close()
        return

    ds_dir = Path(args.out) / args.id
    ds_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    manifest = {"id": args.id, "label": args.label or args.id, "kind": args.kind,
                "pixelSizeUm": args.pixel_size, "created": datetime.now().isoformat(timespec="seconds"),
                "source": {}, "expr": [], "overlays": []}

    if args.kind == "akoya":
        files = akoya_files(Path(args.src))
        if not files:
            sys.exit(f"no TIFFs in {args.src}")
        log(f"[akoya] {len(files)} markers: " + ", ".join(n for n, _ in files))
        src = pyramid.TiffFileSource(files, rgb_split=args.rgb_split)
        manifest["source"] = {"folder": str(args.src), "files": [str(p) for _, p in files]}
    elif args.kind == "xenium":
        import re as _re
        import tifffile as _tf
        folder = Path(args.src)
        files = sorted((folder / "morphology_focus").glob("*.ome.tif*"))
        if not files:
            sys.exit(f"no morphology_focus/*.ome.tif in {folder}")
        with _tf.TiffFile(str(files[0])) as tf0:           # the OME block lists every channel
            ome = tf0.ome_metadata or ""
        names = _re.findall(r'<Channel[^>]*?\sName="([^"]*)"', ome)
        px = _re.findall(r'PhysicalSizeX="([^"]*)"', ome)
        pairs = []
        for f in files:
            m = _re.match(r"ch(\d+)_", f.name)
            idx = int(m.group(1)) if m else len(pairs)
            pairs.append((names[idx] if idx < len(names) else f.stem.split("_", 1)[-1], f))
        src = pyramid.TiffFileSource(pairs)
        if args.pixel_size is None and px:
            manifest["pixelSizeUm"] = float(px[0])
        ann = folder / "Custom_annotation.geojson"
        if ann.exists() and not args.geojson:
            args.geojson = [f"Annotations={ann}"]
        manifest["source"] = {"folder": str(folder), "files": [str(p) for _, p in pairs]}
    elif args.kind == "tiff":
        names = args.channels.split(",") if args.channels else None
        src = pyramid.TiffStackSource(Path(args.src), names=names, rgb_split=True)
        manifest["source"] = {"file": str(args.src)}
    else:
        src = pyramid.TiffStackSource(Path(args.image), rgb_split=True)
        manifest["source"] = {"image": str(args.image), "outs": str(args.outs)}

    if getattr(args, "crop", None):
        try:
            box = [int(v) for v in str(args.crop).replace(" ", "").split(",")]
            if len(box) != 4:
                raise ValueError
        except ValueError:
            raise SystemExit(f"--crop wants four numbers x0,y0,x1,y1 (got {args.crop!r})")
        full_w, full_h = src.width, src.height
        src = pyramid.CropSource(src, box)
        manifest["crop"] = {"box": list(src.box), "ofFullSize": [full_w, full_h]}
        log(f"[image] crop {src.box} of {full_w}x{full_h} -> {src.width}x{src.height}")

    flip = (getattr(args, "flip", "") or "").lower()
    if flip:
        if set(flip) - set("hv"):
            raise SystemExit(f"--flip takes h, v or hv (got {args.flip!r})")
        src = pyramid.FlipSource(src, flip)
        manifest["flip"] = flip
        log(f"[image] flip {flip} ({'horizontal ' if 'h' in flip else ''}{'vertical' if 'v' in flip else ''})")

    if manifest["pixelSizeUm"] is None and getattr(src, "pixel_size_um", None):
        manifest["pixelSizeUm"] = float(src.pixel_size_um)
        log(f"[image] pixel size from file: {manifest['pixelSizeUm']:.4f} µm/px")
    try:
        manifest["image"] = pyramid.build_image(src, ds_dir, tile=args.tile, log=log)
    finally:
        src.close()

    if args.kind == "visium":
        layers = find_visium_layers(Path(args.outs))
        keep = set(args.bins.split(",")) if args.bins else None
        for layer in layers:
            if keep and layer["key"] not in keep:
                continue
            meta = build_expr(layer, ds_dir, image_scale=args.image_scale, log=log,
                              crop_box=(manifest.get("crop") or {}).get("box"),
                              flip=manifest.get("flip", ""),
                              image_size=(manifest["image"]["width"], manifest["image"]["height"]))
            manifest["expr"].append(meta)
            if manifest["pixelSizeUm"] is None and meta.get("micronsPerPixel"):
                manifest["pixelSizeUm"] = meta["micronsPerPixel"]
        if not layers:
            log("[visium] WARNING: no count layers found under", args.outs)

    manifest["overlays"] = parse_geojson_args(args.geojson, ds_dir)
    with open(ds_dir / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=1)
    log(f"[done] {ds_dir}\\manifest.json  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
