r"""RNA-seq (single-cell / single-nucleus) -> static UMAP dataset for the viewer.

Inputs (pick one):
    --h5ad file.h5ad                       AnnData from scanpy, or Seurat via SeuratDisk/sceasy.
                                           Needs an embedding in obsm (X_umap by default).
    --mtx <10x folder> --umap umap.csv     matrix.mtx(.gz) + features.tsv(.gz) + barcodes.tsv(.gz),
                                           a CSV of barcode,UMAP1,UMAP2, and optional --obs CSVs
                                           (barcode,<column>...) such as clusters or sample.

Output  <dest>/<id>/   (upload like the image datasets)
    manifest.json      kind "rnaseq", cell/gene counts, UMAP bounds, color-by columns
    umap.bin           float32 [x, y] per cell
    obs/<n>.bin        uint16 category code per cell   (or float32 for numeric columns)
    genes.json         [{n: name, p: fraction of cells expressing, m: max value}]
    gene/<NAME>.bin    uint32 k, uint32[k] cell index, float32[k] value   (expressing cells only)

Values: if the matrix holds integer counts they are normalised to log1p(counts per 10k);
otherwise they are exported as provided (already normalised).

HIPAA: every sample ID in column names, category values and the label is replaced
through deid.py before writing, cell barcodes are never exported, and the manifest is
checked for leftover IDs - the export aborts rather than publish one.
"""
from __future__ import annotations

import csv
import gzip
import json
import time
from pathlib import Path
from urllib.parse import quote

import h5py
import numpy as np
from scipy import sparse

import deid

MAX_CATEGORIES = 200          # columns with more distinct values (barcodes, ids) are not exported


def _encode(name: str) -> str:
    return quote(name, safe="-_.!~*'()")


def _s(x) -> str:
    return x.decode("utf-8", "replace") if isinstance(x, (bytes, bytearray)) else str(x)


# ---------------------------------------------------------------------------------------
# readers -> (X csr [cells x genes], genes[list], umap (n,2), obs{name: ndarray}, source-note)
# ---------------------------------------------------------------------------------------
def _h5_sparse(g) -> sparse.spmatrix:
    enc = _s(g.attrs.get("encoding-type", g.attrs.get("h5sparse_format", "csr")))
    shape = tuple(int(v) for v in (g.attrs["shape"] if "shape" in g.attrs else g.attrs["h5sparse_shape"]))
    data, indices, indptr = g["data"][()], g["indices"][()], g["indptr"][()]
    cls = sparse.csc_matrix if "csc" in enc else sparse.csr_matrix
    return cls((data, indices, indptr), shape=shape)


def _h5_column(grp, name):
    obj = grp[name]
    if isinstance(obj, h5py.Group):                               # categorical (anndata >= 0.8)
        cats = np.array([_s(c) for c in obj["categories"][()]], dtype=object)
        codes = obj["codes"][()]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = "NA"
        return out
    arr = obj[()]
    if arr.dtype.kind in "OSU":
        return np.array([_s(v) for v in arr], dtype=object)
    legacy = grp.file.get(f"uns/{name}_categories")               # anndata < 0.7
    if legacy is not None and arr.dtype.kind in "iu":
        cats = np.array([_s(c) for c in legacy[()]], dtype=object)
        return cats[arr]
    return arr


def read_h5ad(path: Path, embedding: str = "X_umap"):
    f = h5py.File(path, "r")
    X = f["X"]
    X = _h5_sparse(X) if isinstance(X, h5py.Group) else sparse.csr_matrix(X[()])
    var = f["var"]
    genes = None
    for col in ("gene_symbols", "feature_name", "gene_name", "gene_names", "features", "name"):
        if col in var and not isinstance(var[col], h5py.Group) and var[col].dtype.kind in "OSU":
            genes = [_s(v) for v in var[col][()]]
            break
    if genes is None:
        idx = _s(var.attrs.get("_index", "_index"))
        genes = [_s(v) for v in var[idx][()]]
    obsm = f.get("obsm")
    emb = None
    if obsm is not None:
        for k in (embedding, embedding.upper(), "X_UMAP", "umap", "UMAP"):
            if k in obsm:
                emb = np.asarray(obsm[k][()], dtype=np.float32)[:, :2]
                break
    if emb is None:
        raise SystemExit(f"{path.name}: no '{embedding}' in obsm - export the UMAP from Seurat/scanpy first "
                         f"(found: {list(obsm.keys()) if obsm is not None else 'no obsm'})")
    obs = {}
    og = f["obs"]
    idx = _s(og.attrs.get("_index", "_index"))
    for name in og.keys():
        if name == idx:
            continue
        try:
            obs[name] = _h5_column(og, name)
        except Exception:
            pass
    return X.tocsr(), genes, emb, obs


def _open_text(p: Path):
    return gzip.open(p, "rt", encoding="utf-8") if p.suffix == ".gz" else open(p, encoding="utf-8")


def _find(folder: Path, stem: str) -> Path:
    for cand in (folder / stem, folder / f"{stem}.gz"):
        if cand.exists():
            return cand
    raise SystemExit(f"missing {stem}(.gz) in {folder}")


def read_10x(folder: Path, umap_csv: Path, obs_csvs: list[Path]):
    from scipy.io import mmread
    with _open_text(_find(folder, "barcodes.tsv")) as fh:
        barcodes = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
    feats = []
    with _open_text(_find(folder, "features.tsv") if (folder / "features.tsv").exists() or (folder / "features.tsv.gz").exists()
                    else _find(folder, "genes.tsv")) as fh:
        for ln in fh:
            parts = ln.rstrip("\n").split("\t")
            if parts and parts[0]:
                feats.append((parts[1] if len(parts) > 1 else parts[0], parts[2] if len(parts) > 2 else "Gene Expression"))
    mtx = _find(folder, "matrix.mtx")
    M = mmread(gzip.open(mtx) if mtx.suffix == ".gz" else str(mtx))          # genes x cells
    X = sparse.csr_matrix(M).T.tocsr()                                        # cells x genes
    keep = np.array([t == "Gene Expression" for _, t in feats])
    genes = [n for (n, _), k in zip(feats, keep) if k]
    X = X[:, np.where(keep)[0]]

    pos = {b: i for i, b in enumerate(barcodes)}
    emb = np.full((len(barcodes), 2), np.nan, np.float32)
    with open(umap_csv, newline="", encoding="utf-8") as fh:
        rd = csv.reader(fh)
        next(rd)
        for row in rd:
            i = pos.get(row[0])
            if i is not None:
                emb[i] = (float(row[1]), float(row[2]))
    obs = {}
    for oc in obs_csvs:
        with open(oc, newline="", encoding="utf-8") as fh:
            rd = csv.reader(fh)
            header = next(rd)
            cols = {h: np.full(len(barcodes), "NA", dtype=object) for h in header[1:]}
            for row in rd:
                i = pos.get(row[0])
                if i is None:
                    continue
                for h, v in zip(header[1:], row[1:]):
                    cols[h][i] = v
        for h, v in cols.items():
            name = h if h.lower() not in ("cluster", "clusters") else oc.parent.name.replace("gene_expression_", "")
            obs[name] = v
    ok = ~np.isnan(emb[:, 0])                                                 # cells without a UMAP are dropped
    return X[ok], genes, emb[ok], {k: v[ok] for k, v in obs.items()}


# ---------------------------------------------------------------------------------------
def write_obs(out: Path, obs: dict, ds_id: str, colors: dict | None = None, log=print) -> list:
    """Color-by columns, de-identified. Sample-identity columns are replaced value-for-value
    by numbers from the private key; every other text value goes through the ID pattern.
    `colors` = {column: {original_value: "#rrggbb"}} (from uns/<col>_colors) keeps the
    analyst's palette."""
    (out / "obs").mkdir(parents=True, exist_ok=True)
    obs_meta = []
    for i, (name, values) in enumerate(obs.items()):
        safe_name = deid.anonymize(str(name), form="text", context=f"obs column {ds_id}")
        v = np.asarray(values)
        if v.dtype.kind in "fiu":
            vals = v.astype(np.float32)
            uniq = np.unique(vals[~np.isnan(vals)])
            if len(uniq) <= 30 and np.all(np.mod(uniq, 1) == 0):             # small integer codes (cluster ids)
                v = np.array([str(int(x)) if not np.isnan(x) else "NA" for x in vals], dtype=object)
            else:
                fn = f"{i:02d}.bin"
                (out / "obs" / fn).write_bytes(np.ascontiguousarray(vals, dtype="<f4").tobytes())
                finite = vals[np.isfinite(vals)]
                obs_meta.append({"name": safe_name, "type": "numeric", "file": fn,
                                 "min": float(np.percentile(finite, 1)) if finite.size else 0.0,
                                 "max": float(np.percentile(finite, 99)) if finite.size else 1.0})
                continue
        raw_cats, raw_codes = np.unique(np.array([_s(x) for x in v], dtype=object), return_inverse=True)
        if len(raw_cats) > MAX_CATEGORIES or len(raw_cats) < 2:
            continue
        if deid.is_sample_column(name, raw_cats):
            anon = [deid.sample_value(c, context=f"sample column '{name}' {ds_id}") for c in raw_cats]
            log(f"[deid] column '{safe_name}': {len(raw_cats)} sample values replaced by numbers")
        else:
            anon = [deid.anonymize(c, form="value", context=f"obs {name} {ds_id}") for c in raw_cats]
        cats = sorted(set(anon), key=lambda c: (not c.lstrip("-").isdigit(), int(c) if c.lstrip("-").isdigit() else 0, c))
        pos = {c: k for k, c in enumerate(cats)}
        remap = np.array([pos[a] for a in anon], dtype=np.uint16)
        codes = remap[raw_codes]
        fn = f"{i:02d}.bin"
        (out / "obs" / fn).write_bytes(np.ascontiguousarray(codes, dtype="<u2").tobytes())
        meta = {"name": safe_name, "type": "category", "file": fn, "categories": cats,
                "counts": [int(c) for c in np.bincount(codes, minlength=len(cats))]}
        pal = (colors or {}).get(name)
        if pal:
            first = {}
            for rc, a in zip(raw_cats, anon):
                first.setdefault(a, pal.get(rc))
            if all(first.get(c) for c in cats):
                meta["colors"] = [first[c] for c in cats]
        obs_meta.append(meta)
    log(f"[rnaseq] {len(obs_meta)} color-by column(s): {[o['name'] for o in obs_meta]}")
    return obs_meta


def export_rnaseq(X, genes, emb, obs, dest: Path, ds_id: str, label: str, log=print) -> Path:
    t0 = time.time()
    ds_id = deid.anonymize(ds_id, form="slug", context=f"rnaseq id")
    label = deid.anonymize(label, form="text", context=f"rnaseq label {ds_id}")
    out = Path(dest) / ds_id
    (out / "gene").mkdir(parents=True, exist_ok=True)
    (out / "obs").mkdir(exist_ok=True)
    n_cells, n_genes = X.shape
    log(f"[rnaseq] {n_cells:,} cells x {n_genes:,} genes -> {out}")

    X = X.astype(np.float32)
    data = X.data
    is_counts = data.size > 0 and np.all(np.mod(data[: min(200000, data.size)], 1) == 0) and data.max() > 20
    if is_counts:
        totals = np.asarray(X.sum(axis=1)).ravel()
        scale = np.where(totals > 0, 1e4 / totals, 0).astype(np.float32)
        X = sparse.diags(scale) @ X
        X.data = np.log1p(X.data)
        value_label = "log1p(counts per 10k)"
    else:
        value_label = "expression (as provided)"
    X = X.tocsc()
    X.sort_indices()
    log(f"[rnaseq] values: {value_label}  ({time.time() - t0:.0f}s)")

    (out / "umap.bin").write_bytes(np.ascontiguousarray(emb, dtype="<f4").tobytes())

    obs_meta = write_obs(out, obs, ds_id, colors=None, log=log)

    # ---- genes -----------------------------------------------------------------------
    listing = []
    for gi, name in enumerate(genes):
        a, b = X.indptr[gi], X.indptr[gi + 1]
        k = int(b - a)
        if k == 0:
            continue
        idx = X.indices[a:b].astype("<u4")
        val = X.data[a:b].astype("<f4")
        (out / "gene" / f"{_encode(name)}.bin").write_bytes(np.uint32(k).tobytes() + idx.tobytes() + val.tobytes())
        listing.append({"n": name, "p": round(k / n_cells, 4), "m": _safe_max(val)})
    listing.sort(key=lambda d: -d["p"])
    (out / "genes.json").write_text(json.dumps(listing, separators=(",", ":"), allow_nan=False), encoding="utf-8")
    log(f"[rnaseq] {len(listing):,} gene files  ({time.time() - t0:.0f}s)")

    lo, hi = np.nanmin(emb, axis=0), np.nanmax(emb, axis=0)
    manifest = {"id": ds_id, "label": label, "kind": "rnaseq", "static": True,
                "nCells": int(n_cells), "nGenes": len(listing), "valueLabel": value_label,
                "bounds": [float(lo[0]), float(lo[1]), float(hi[0]), float(hi[1])], "obs": obs_meta}
    deid.assert_clean({k: v for k, v in manifest.items()}, f"manifest of {ds_id}")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    log(f"[rnaseq] done -> {out}  ({time.time() - t0:.0f}s)")
    return out


# ---------------------------------------------------------------------------------------
# Very large .h5ad (tens of GB): stream the sparse counts, never load the matrix.
# ---------------------------------------------------------------------------------------
QC_COLUMNS = {"n_genes", "n_genes_by_counts", "total_counts", "pct_counts_mt", "doublet_score"}


def export_h5ad_streaming(path: Path, dest: Path, ds_id: str, label: str, embedding: str = "X_umap",
                          tmpdir: Path | None = None, min_cells: int = 3, chunk_nnz: int = 20_000_000,
                          log=print) -> Path:
    """Two passes over the sparse matrix (raw/X if present, else X):
    1) count non-zeros per gene;  2) normalise each cell (log1p counts-per-10k when the values
    are integer counts) and scatter every entry into a gene-major memmap on local disk.
    Then one small file per gene.  Peak RAM ~1 GB per 20M non-zeros."""
    import tempfile
    t0 = time.time()
    f = h5py.File(path, "r")
    src = f["raw"]["X"] if "raw" in f and "X" in f["raw"] else f["X"]
    if not isinstance(src, h5py.Group):
        raise SystemExit("streaming needs a sparse matrix (raw/X or X as csr_matrix)")
    enc = _s(src.attrs.get("encoding-type", "csr_matrix"))
    if "csr" not in enc:
        raise SystemExit(f"streaming expects csr_matrix, found {enc}")
    var = f["raw"]["var"] if src.name.startswith("/raw") else f["var"]
    genes = [_s(g) for g in var[_s(var.attrs.get("_index", "_index"))][()]]
    n_genes = len(genes)
    indptr = src["indptr"][()].astype(np.int64)
    n_cells, nnz = len(indptr) - 1, int(indptr[-1])
    emb = None
    for k in (embedding, "X_umap", "X_UMAP"):
        if k in f["obsm"]:
            emb = np.asarray(f["obsm"][k][()], dtype=np.float32)[:, :2]
            break
    if emb is None:
        raise SystemExit(f"no {embedding} in obsm")

    ds_id = deid.anonymize(ds_id, form="slug", context="rnaseq id")
    label = deid.anonymize(label, form="text", context=f"rnaseq label {ds_id}")
    out = Path(dest) / ds_id
    (out / "gene").mkdir(parents=True, exist_ok=True)
    log(f"[rnaseq] {n_cells:,} cells x {n_genes:,} genes, {nnz:,} non-zeros (source {src.name}) -> {out}")

    # ---- obs + palettes (small) -----------------------------------------------------------
    og = f["obs"]
    idx_name = _s(og.attrs.get("_index", "_index"))
    obs, colors = {}, {}
    for name in og.keys():
        if name == idx_name:
            continue
        obj = og[name]
        if isinstance(obj, h5py.Dataset) and obj.dtype.kind in "fiu" and name not in QC_COLUMNS:
            continue                                   # skip unknown numeric columns (coordinates etc.)
        try:
            obs[name] = _h5_column(og, name)
        except Exception:
            continue
        if isinstance(obj, h5py.Group) and f.get(f"uns/{name}_colors") is not None:
            cats = [_s(c) for c in obj["categories"][()]]
            pal = [_s(c) for c in f[f"uns/{name}_colors"][()]]
            if len(pal) == len(cats):
                colors[name] = {c: (p if p.startswith("#") else "#" + p)[:7] for c, p in zip(cats, pal)}
    (out / "umap.bin").write_bytes(np.ascontiguousarray(emb, dtype="<f4").tobytes())
    obs_meta = write_obs(out, obs, ds_id, colors=colors, log=log)

    # ---- values: counts -> log1p CP10k -------------------------------------------------
    probe = src["data"][: min(nnz, 2_000_000)]
    is_counts = bool(probe.size) and bool(np.all(np.mod(probe, 1) == 0)) and float(probe.max()) > 20
    if is_counts:
        if "total_counts" in og and not isinstance(og["total_counts"], h5py.Group) and src.name.startswith("/raw"):
            totals = og["total_counts"][()].astype(np.float64)
        else:
            totals = np.zeros(n_cells)
            for r0, r1 in _row_chunks(indptr, chunk_nnz):
                a, b = indptr[r0], indptr[r1]
                d = src["data"][a:b]
                totals[r0:r1] = np.add.reduceat(d, indptr[r0:r1] - a) if b > a else 0
        scale = np.where(totals > 0, 1e4 / np.maximum(totals, 1e-9), 0.0).astype(np.float32)
        value_label = "log1p(counts per 10k)"
    else:
        scale = None
        value_label = "expression (as provided)"
    log(f"[rnaseq] values: {value_label}")

    # ---- pass 1: non-zeros per gene ------------------------------------------------------
    per_gene = np.zeros(n_genes, np.int64)
    for r0, r1 in _row_chunks(indptr, chunk_nnz):
        per_gene += np.bincount(src["indices"][indptr[r0]:indptr[r1]], minlength=n_genes)
    offsets = np.concatenate([[0], np.cumsum(per_gene)])
    log(f"[rnaseq] pass 1 done ({time.time() - t0:.0f}s)")

    # ---- pass 2: scatter into gene-major memmaps -------------------------------------------
    tmp = Path(tempfile.mkdtemp(prefix="rnaseq_", dir=str(tmpdir) if tmpdir else None))
    try:
        mm_idx = np.lib.format.open_memmap(tmp / "idx.npy", mode="w+", dtype=np.uint32, shape=(nnz,))
        mm_val = np.lib.format.open_memmap(tmp / "val.npy", mode="w+", dtype=np.float32, shape=(nnz,))
        ptr = offsets[:-1].copy()
        done = 0
        for r0, r1 in _row_chunks(indptr, chunk_nnz):
            a, b = indptr[r0], indptr[r1]
            if b == a:
                continue
            ind = src["indices"][a:b]
            val = src["data"][a:b].astype(np.float32)
            rows = np.repeat(np.arange(r0, r1, dtype=np.uint32), np.diff(indptr[r0:r1 + 1]))
            if scale is not None:
                val = np.log1p(val * scale[rows])
            order = np.argsort(ind, kind="stable")
            g = ind[order]
            cc = np.bincount(g, minlength=n_genes)
            first = np.cumsum(cc) - cc
            pos = ptr[g] + (np.arange(len(g), dtype=np.int64) - first[g])
            mm_idx[pos] = rows[order]
            mm_val[pos] = val[order]
            ptr += cc
            done += b - a
            log(f"[rnaseq]   pass 2 {done / nnz:5.1%}  ({time.time() - t0:.0f}s)")
        mm_idx.flush(); mm_val.flush()

        # ---- one file per gene --------------------------------------------------------------
        listing = []
        for gi, name in enumerate(genes):
            a, b = int(offsets[gi]), int(offsets[gi + 1])
            k = b - a
            if k < min_cells:
                continue
            v = np.array(mm_val[a:b], dtype="<f4")          # copies: no view may keep the memmap open
            (out / "gene" / f"{_encode(name)}.bin").write_bytes(
                np.uint32(k).tobytes() + np.array(mm_idx[a:b], dtype="<u4").tobytes() + v.tobytes())
            listing.append({"n": name, "p": round(k / n_cells, 4), "m": _safe_max(v)})
            if len(listing) % 5000 == 0:
                log(f"[rnaseq]   {len(listing):,} gene files ({time.time() - t0:.0f}s)")
    finally:
        import gc
        import shutil
        for mm in ("mm_idx", "mm_val"):                    # Windows cannot delete a mapped file
            m = locals().get(mm)
            if m is not None and getattr(m, "_mmap", None) is not None:
                m._mmap.close()
        mm_idx = mm_val = None
        gc.collect()
        shutil.rmtree(tmp, ignore_errors=True)
        if tmp.exists():
            log(f"[rnaseq] WARNING: could not remove temporary folder {tmp} - delete it by hand")
    listing.sort(key=lambda d: -d["p"])
    (out / "genes.json").write_text(json.dumps(listing, separators=(",", ":"), allow_nan=False), encoding="utf-8")

    lo, hi = np.nanmin(emb, axis=0), np.nanmax(emb, axis=0)
    manifest = {"id": ds_id, "label": label, "kind": "rnaseq", "static": True,
                "nCells": int(n_cells), "nGenes": len(listing), "valueLabel": value_label,
                "bounds": [float(lo[0]), float(lo[1]), float(hi[0]), float(hi[1])], "obs": obs_meta}
    deid.assert_clean(manifest, f"manifest of {ds_id}")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    f.close()
    log(f"[rnaseq] done -> {out}  ({time.time() - t0:.0f}s)")
    return out


def _safe_max(v) -> float:
    """Largest finite value (JSON cannot hold NaN/Infinity)."""
    v = np.asarray(v, dtype=np.float64)
    v = v[np.isfinite(v)]
    return round(float(v.max()), 3) if v.size else 0.0


def repair_gene_listing(ds_dir: Path, log=print) -> int:
    """Recompute any non-finite "m" in genes.json from the gene files themselves."""
    p = Path(ds_dir) / "genes.json"
    listing = json.loads(p.read_text(encoding="utf-8").replace(":NaN", ":null").replace(":Infinity", ":null"))
    fixed = 0
    for e in listing:
        if not isinstance(e.get("m"), (int, float)):
            b = (Path(ds_dir) / "gene" / f"{_encode(e['n'])}.bin").read_bytes()
            k = int(np.frombuffer(b, np.uint32, 1)[0])
            e["m"] = _safe_max(np.frombuffer(b, np.float32, k, 4 + 4 * k))
            fixed += 1
    p.write_text(json.dumps(listing, separators=(",", ":"), allow_nan=False), encoding="utf-8")
    log(f"[rnaseq] genes.json: {fixed} entr{'y' if fixed == 1 else 'ies'} repaired")
    return fixed


def _row_chunks(indptr: np.ndarray, chunk_nnz: int):
    """Row ranges holding ~chunk_nnz non-zeros each."""
    n = len(indptr) - 1
    r0 = 0
    while r0 < n:
        target = indptr[r0] + chunk_nnz
        r1 = int(np.searchsorted(indptr, target, side="right")) - 1
        r1 = max(r1, r0 + 1)
        r1 = min(r1, n)
        yield r0, r1
        r0 = r1
