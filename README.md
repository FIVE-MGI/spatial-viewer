# Spatial Viewer

A browser viewer for large multiplexed microscopy and spatial transcriptomics, with
the pipeline that prepares the data for it.

It handles three kinds of dataset:

- **Multiplexed immunofluorescence** — tens of antibody channels on one section, each
  with its own color and display window, composited on the GPU.
- **Spatial transcriptomics** — whole-transcriptome counts on a bin grid, registered to
  the section image, drawn per gene or as a combined density.
- **Single-cell RNA-seq** — a UMAP with named clusters and per-gene expression panels.

Images are never loaded whole. Each visible channel is a tile layer reading an 8-bit
grey tile pyramid, windowed and tinted in a shader, so a 30-channel slide of tens of
thousands of pixels square opens immediately and stays responsive.

## Layout

```
web/SpatialViewer/     the viewer, and the pipeline that builds its inputs
  viewer.js            all of the client
  index.html
  styles.css
  config.example.js    copy to config.js and edit
  prepare.py           CLI: ingest a dataset
  pyramid.py           image -> zarr pyramid + per-channel display windows
  rnaseq.py            h5ad -> UMAP, cluster columns, per-gene vectors
  expr.py              transcript/bin expression storage
  export.py            zarr -> a flat tile tree for static hosting
  deid.py              de-identification of names and labels
  server.py            optional live tile server
  static_serve.py      serves an export folder locally, with CORS
api/token/             password check -> short-lived read-only SAS (Azure Function)
tools/label_editor.html  place UMAP cluster labels by hand
```

## Running it

Prepare a dataset, export it, serve it:

```bash
pip install -r web/SpatialViewer/requirements.txt
python web/SpatialViewer/prepare.py --help
python web/SpatialViewer/export.py          # zarr -> static tiles
python web/SpatialViewer/static_serve.py /path/to/exports 8765
```

Then copy `config.example.js` to `config.js`, list your datasets in it, serve the
`web/` folder, and open:

```
/SpatialViewer/index.html?static=http://localhost:8765
```

`?static=` reads an export folder directly and skips the password, which is how you
check a dataset before uploading it.

## Hosting with a password

For data that should not be public, put the exports in a private container and deploy
`api/token/` alongside the site. The browser posts a shared password to `/api/token`;
the function verifies it server-side and returns a container-scoped, read-only SAS that
expires. The storage key never reaches the browser.

Settings the function reads:

| Variable | Meaning |
|---|---|
| `VIEWER_PASSWORD` | the shared password |
| `STORAGE_KEY` | account key; mints a fresh SAS per session. Preferred. |
| `BLOB_SAS` | a ready-made read-only SAS, used only when no key is available |
| `STORAGE_ACCOUNT`, `STORAGE_CONTAINER` | where the exports live |
| `TOKEN_HOURS`, `SESSION_HOURS` | lifetimes, both one hour by default |

The password is part of the session cookie's signing secret, so changing it ends every
existing session immediately.

## Cluster labels

Automatic placement keeps a label off the thick of its own cluster, out of its
neighbours, inside the frame, and on the outward side of the plot. It is decent and it
is not a person looking at a figure, so `tools/label_editor.html` lets you drag each
label where it belongs and emits coordinates for `labelPositions` in config. Anything
placed by hand wins; anything missing is still placed automatically.

## Configuration

`config.js` is the only file that knows about a particular study — its datasets,
institution, logos and byline. It is gitignored. `config.example.js` documents every
option.

## License

Not yet chosen. Until a license is added, no rights are granted: default copyright
applies and the code may not be copied, used or redistributed.
