# Spatial Viewer

A browser viewer for spatial biology data that is too large to open: multiplexed
immunofluorescence, spatial transcriptomics, and single-cell RNA-seq UMAPs.

Slides are converted once into a pyramid of small tiles. The browser fetches only the
tiles for the region and channels on screen, so a 30-channel slide opens immediately
and pans without downloading the rest.

## Quick start

```bash
pip install -r web/SpatialViewer/requirements.txt
cd web/SpatialViewer

python prepare.py tiff --src slide.qptiff --id sample1 --out cache/
python prepare.py export --id sample1 --out cache/ --dest exports/
python static_serve.py exports/ 8765
```

Copy `config.example.js` to `config.js`, list `sample1` in it, serve `web/`, and open
`/SpatialViewer/index.html?static=http://localhost:8765`.

## Ingest

| Input | Command |
|---|---|
| OME-TIFF, QPTIFF, TIFF stack | `prepare.py tiff` |
| Folder of per-marker greyscale TIFFs | `prepare.py akoya` |
| Visium / Visium HD `outs/` | `prepare.py visium` |
| Xenium output bundle | `prepare.py xenium` |
| AnnData `.h5ad`, 10x matrix | `prepare.py rnaseq` |

`prepare.py inspect` reports what is in a file; `prepare.py export` turns a cache into
the static tree.

## Private data

The viewer implements no authentication. Point `tokenEndpoint` in config at anything
satisfying:

```
POST <tokenEndpoint>  {"password": "..."}   -> 200 {base, sas, expires}  |  401
GET  <tokenEndpoint>                        -> same, via a session cookie
```

`base` is the prefix the datasets sit under, `sas` a query string appended to requests
for them. Empty `tokenEndpoint` reads plain URLs.

## Configuration

`config.js` holds everything specific to a study — datasets, logos, byline — and is
gitignored. `config.example.js` documents every option.

`tools/label_editor.html` places UMAP cluster labels by hand when the automatic
placement is not good enough.

## License

Not yet chosen. Until one is added, no rights are granted.
