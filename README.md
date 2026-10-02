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

## Hosting data that should not be public

The viewer does not implement authentication; it consumes a token endpoint, so the
backend is yours to choose. Point `tokenEndpoint` in config at anything that satisfies
this contract:

```
POST <tokenEndpoint>   {"password": "..."}      -> 200 {base, sas, expires}
                                                -> 401 wrong or missing
GET  <tokenEndpoint>                            -> same, authorised by a session cookie
```

`base` is the URL prefix the datasets sit under and `sas` is a query string appended to
every request for them; `expires` is epoch milliseconds, and the viewer renews shortly
before it. With `tokenEndpoint` empty the viewer reads plain URLs instead, which is what
`?static=` does for local testing.

A sound implementation keeps the storage credential server-side and returns something
scoped and short-lived rather than the credential itself, ties its session cookie to the
password so that changing the password ends existing sessions, and rate-limits attempts.

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
