r"""Tile + transcript server for the SpatialViewer page in this folder.

    python server.py --cache <export root> [--port 8760] [--web ..]

API (all GET):
    /api/datasets                                   list
    /api/datasets/{id}                              manifest
    /api/datasets/{id}/tiles/{c}/{z}/{x}/{y}.png    one channel, 8-bit grey  (?fmt=jpg, ?lo=&hi= raw-unit override)
    /api/datasets/{id}/thumb/{c}.png
    /api/datasets/{id}/expr/{key}/genes             [{n, t}]   name + total count
    /api/datasets/{id}/expr/{key}/gene/{name}       float32 [x,y,v]* for bins with v>0
    /api/datasets/{id}/expr/{key}/positions         float32 [x,y]*  every in-tissue bin
    /api/datasets/{id}/overlays/{key}               GeoJSON

The repo's web/ folder is served at / so http://host:8760/SpatialViewer/ works
as a single process; the Azure-hosted copy of the page can point here with
?api=http://host:8760 instead (CORS is open). Everything is behind one shared
password (HTTP Basic) on /api/*: --password / SPATIALVIZ_PASSWORD / a .password file.
The page files are served openly; the server's own sources are never served.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import secrets
import threading
from functools import lru_cache
from pathlib import Path
from typing import Optional

import numpy as np
from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from expr import ExprReader
from pyramid import PyramidReader

CACHE_ROOT = Path(os.environ.get("SPATIALVIZ_CACHE", "cache"))
PASSWORD: Optional[str] = os.environ.get("SPATIALVIZ_PASSWORD") or None
app = FastAPI(title="SpatialViz tiles")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["GET"], allow_headers=["*"])


# Files in this folder that must never be served (the page's own html/js/css are fine).
_BLOCKED_SUFFIXES = (".py", ".pyc", ".ps1", ".txt", ".md", ".password")
_BLOCKED_PARTS = ("/.venv/", "/__pycache__/", "/.password")


def _blocked_path(path: str) -> bool:
    low = path.lower()
    if not low.startswith("/spatialviewer/"):
        return False
    if low.endswith(".link.txt"):          # dataset link files are part of the page
        return False
    return low.endswith(_BLOCKED_SUFFIXES) or any(p in low for p in _BLOCKED_PARTS)


@app.middleware("http")
async def require_password(request: Request, call_next):
    """The data API (/api/*) sits behind one shared password (HTTP Basic, any
    username). The viewer page itself is served freely - it is public on the
    Azure site anyway - and asks for the password in its own login box, which
    works the same whether it is served from here or from Azure with ?api=.
    Server sources next to the page are never served."""
    path = request.url.path
    if _blocked_path(path):
        return Response("not found", status_code=404)
    if PASSWORD is None or request.method == "OPTIONS" or not path.startswith("/api/"):
        return await call_next(request)
    ok = False
    auth = request.headers.get("authorization", "")
    if auth.lower().startswith("basic "):
        try:
            _user, _, pw = base64.b64decode(auth[6:]).decode("utf-8", "replace").partition(":")
            ok = secrets.compare_digest(pw.encode(), PASSWORD.encode())
        except Exception:
            ok = False
    if ok:
        return await call_next(request)
    # The 401 is produced outside CORSMiddleware, so add the CORS headers by
    # hand or a cross-origin page cannot even read the status code.
    # Deliberately no WWW-Authenticate header: the page has its own login box,
    # and the browser's native dialog on top of it only confuses people.
    return Response("password required", status_code=401,
                    headers={"Access-Control-Allow-Origin": "*",
                             "Access-Control-Expose-Headers": "*"})

_lock = threading.Lock()
_readers: dict[str, PyramidReader] = {}
_expr: dict[tuple[str, str], ExprReader] = {}
_manifests: dict[str, tuple[float, dict]] = {}


def _scan() -> dict[str, Path]:
    if not CACHE_ROOT.is_dir():
        return {}
    return {d.name: d for d in sorted(CACHE_ROOT.iterdir()) if (d / "manifest.json").exists()}


def manifest(ds_id: str) -> dict:
    d = _scan().get(ds_id)
    if d is None:
        raise HTTPException(404, "unknown dataset")
    mp = d / "manifest.json"
    mt = mp.stat().st_mtime
    cached = _manifests.get(ds_id)
    if cached and cached[0] == mt:
        return cached[1]
    with open(mp) as f:
        m = json.load(f)
    m["_dir"] = str(d)
    with _lock:
        _manifests[ds_id] = (mt, m)
        _readers.pop(ds_id, None)          # re-prepared: drop stale handles
        for k in [k for k in _expr if k[0] == ds_id]:
            _expr.pop(k)
    return m


def reader(ds_id: str) -> PyramidReader:
    m = manifest(ds_id)
    r = _readers.get(ds_id)
    if r is None:
        with _lock:
            r = _readers.get(ds_id)
            if r is None:
                r = _readers[ds_id] = PyramidReader(Path(m["_dir"]), m)
    return r


def expr_reader(ds_id: str, key: str) -> ExprReader:
    m = manifest(ds_id)
    if not any(e["key"] == key for e in m.get("expr", [])):
        raise HTTPException(404, "unknown expression layer")
    k = (ds_id, key)
    r = _expr.get(k)
    if r is None:
        with _lock:
            r = _expr.get(k)
            if r is None:
                r = _expr[k] = ExprReader(Path(m["_dir"]) / "expr" / f"{key}.zarr")
    return r


@app.get("/api/health")
def health():
    return {"cache": str(CACHE_ROOT.resolve()), "exists": CACHE_ROOT.is_dir(), "datasets": list(_scan())}


@app.get("/api/datasets")
def list_datasets():
    out = []
    for ds_id in _scan():
        try:
            m = manifest(ds_id)
        except Exception:
            continue
        img = m.get("image", {})
        out.append({"id": ds_id, "label": m.get("label", ds_id), "kind": m.get("kind"),
                    "width": img.get("width"), "height": img.get("height"),
                    "channels": len(img.get("channels", [])), "expr": [e["key"] for e in m.get("expr", [])]})
    return out


@app.get("/api/datasets/{ds_id}")
def get_manifest(ds_id: str):
    m = dict(manifest(ds_id))
    m.pop("_dir", None)
    return m


@lru_cache(maxsize=3000)
def _tile_cached(ds_id: str, mtime: float, c: int, z: int, x: int, y: int, fmt: str,
                 lo: Optional[float], hi: Optional[float]) -> Optional[bytes]:
    r = reader(ds_id)
    level = (r.nlevels - 1) - z
    return r.tile_bytes(c, level, x, y, fmt=fmt, lo=lo, hi=hi)


@app.get("/api/datasets/{ds_id}/tiles/{c}/{z}/{x}/{y}.{ext}")
def tile(ds_id: str, c: int, z: int, x: int, y: int, ext: str,
         lo: Optional[float] = None, hi: Optional[float] = None):
    m = manifest(ds_id)
    fmt = "jpg" if ext.lower() in ("jpg", "jpeg") else "png"
    data = _tile_cached(ds_id, _manifests[ds_id][0], c, z, x, y, fmt, lo, hi)
    if data is None:
        return Response(status_code=204)
    return Response(content=data, media_type="image/jpeg" if fmt == "jpg" else "image/png",
                    headers={"Cache-Control": "public, max-age=86400"})


@app.get("/api/datasets/{ds_id}/thumb/{c}.png")
def thumb(ds_id: str, c: int):
    r = reader(ds_id)
    if not (0 <= c < len(r.channels)):
        raise HTTPException(404, "bad channel")
    return Response(content=r.thumbnail(c), media_type="image/png")


@app.get("/api/datasets/{ds_id}/expr/{key}/genes")
def genes(ds_id: str, key: str):
    return JSONResponse(expr_reader(ds_id, key).gene_list(),
                        headers={"Cache-Control": "public, max-age=3600"})


@app.get("/api/datasets/{ds_id}/expr/{key}/gene/{name}")
def gene(ds_id: str, key: str, name: str):
    xyv = expr_reader(ds_id, key).gene_xyv(name)
    if xyv is None:
        raise HTTPException(404, "unknown gene")
    return Response(content=np.ascontiguousarray(xyv, dtype="<f4").tobytes(),
                    media_type="application/octet-stream",
                    headers={"X-Count": str(len(xyv)), "Cache-Control": "public, max-age=3600"})


@app.get("/api/datasets/{ds_id}/expr/{key}/positions")
def positions(ds_id: str, key: str):
    xy = expr_reader(ds_id, key).positions()
    return Response(content=np.ascontiguousarray(xy, dtype="<f4").tobytes(),
                    media_type="application/octet-stream",
                    headers={"X-Count": str(len(xy)), "Cache-Control": "public, max-age=3600"})


@app.get("/api/datasets/{ds_id}/overlays/{key}")
def overlay(ds_id: str, key: str):
    m = manifest(ds_id)
    ov = next((o for o in m.get("overlays", []) if o["key"] == key), None)
    if ov is None:
        raise HTTPException(404, "unknown overlay")
    return FileResponse(Path(m["_dir"]) / "overlays" / ov["file"], media_type="application/geo+json")


def mount_web(web_dir: Path):
    if web_dir.is_dir():
        app.mount("/", StaticFiles(directory=str(web_dir), html=True), name="web")


if __name__ == "__main__":
    import uvicorn
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True, help="folder that prepare.py wrote into")
    ap.add_argument("--port", type=int, default=8760)
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--web", default=str(Path(__file__).resolve().parents[1]),
                    help="static folder to serve at / (the repo's web/)")
    ap.add_argument("--password", help="shared password (or set SPATIALVIZ_PASSWORD, or put it in .password next to this file)")
    ap.add_argument("--no-password", action="store_true", help="serve without a password (local testing only)")
    args = ap.parse_args()
    CACHE_ROOT = Path(args.cache).resolve()
    pw_file = Path(__file__).resolve().parent / ".password"
    if args.no_password:
        PASSWORD = None
    else:
        PASSWORD = args.password or PASSWORD or (pw_file.read_text(encoding="utf-8").strip() if pw_file.exists() else None)
        if not PASSWORD:
            raise SystemExit("Refusing to start without a password: pass --password, set SPATIALVIZ_PASSWORD, "
                             f"write it to {pw_file}, or use --no-password for local testing.")
    mount_web(Path(args.web))
    print(f"cache: {CACHE_ROOT}  datasets: {list(_scan())}", flush=True)
    print(f"open  http://localhost:{args.port}/SpatialViewer/   ({'password protected' if PASSWORD else 'NO PASSWORD'})", flush=True)
    uvicorn.run(app, host=args.host, port=args.port, log_level="warning")
