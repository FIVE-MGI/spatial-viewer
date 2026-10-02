# Spatial Viewer

A browser viewer for spatial biology data that is too large to open, together with the
pipeline that prepares it.

A multiplexed immunofluorescence slide is routinely 30 channels of 16-bit pixels across
tens of thousands of pixels square — tens of gigabytes before anyone has looked at
anything. The usual answers are to send someone a handful of exported PNGs, or to ask
them to install desktop software and download the raw file. This does neither: it
converts a slide once into a pyramid of small tiles, and the browser fetches only the
tiles for the region and channels on screen. Opening a 30-channel slide costs about as
much as opening a web page, and panning it costs the tiles you pan onto.

The same viewer handles three kinds of data, because the questions people ask of them
overlap:

- **Multiplexed immunofluorescence** — tens of antibody channels on one section. Each
  gets its own colour and display window, composited on the GPU.
- **Spatial transcriptomics** — whole-transcriptome counts on a bin grid, registered to
  the section image. Draw one gene, or several combined as a density.
- **Single-cell RNA-seq** — a UMAP with named clusters and per-gene expression panels,
  so expression in tissue can be read against the populations it comes from.

---

## Quick start

No server, no account, no cloud. Prepare a dataset, export it, serve the folder:

```bash
pip install -r web/SpatialViewer/requirements.txt

cd web/SpatialViewer

# ingest: image -> pyramid, counts -> per-gene vectors
python prepare.py tiff --src slide.qptiff --id sample1 --out cache/

# export: pyramid -> a flat tree of tiles any static host can serve
python prepare.py export --id sample1 --out cache/ --dest exports/

# serve the exports with CORS
python static_serve.py exports/ 8765
```

Copy `web/SpatialViewer/config.example.js` to `config.js`, list `sample1` in it, serve
the `web/` folder on another port, and open:

```
/SpatialViewer/index.html?static=http://localhost:8765
```

`?static=` points the viewer straight at an export folder. It is how you check a
dataset before it goes anywhere.

## Preparing data

`prepare.py` has a subcommand per input:

| Input | Command | Produces |
|---|---|---|
| OME-TIFF, QPTIFF, TIFF stack | `prepare.py tiff` | channel pyramid + per-channel display windows |
| A folder of per-marker greyscale TIFFs | `prepare.py akoya` | the same, with markers read from the filenames |
| Visium / Visium HD `outs/` | `prepare.py visium` | pyramid of the section image + per-gene bin vectors |
| Xenium output bundle | `prepare.py xenium` | morphology pyramid + per-gene transcript positions |
| AnnData `.h5ad`, 10x matrix | `prepare.py rnaseq` | UMAP, cluster columns, per-gene vectors |

`prepare.py inspect` reports what is in a file before you commit to ingesting it, and
`prepare.py export` turns a finished cache into the static tree.

Images become a [Zarr](https://zarr.dev) pyramid; `export.py` flattens that into
`tiles/<channel>/<level>/<x>/<y>.png`, which any static host serves without code.
Expression is stored per gene as a sparse `float32` vector, so adding a gene to the
view is one small request rather than a slice of a large matrix.

Cropping is supported for the image ingests (`--crop x0,y0,x1,y1`), and crop and flip
for Visium, with bin coordinates moved to match.

## Serving data that should not be public

The viewer implements no authentication. It consumes a token endpoint, so the backend
is yours. Point `tokenEndpoint` in config at anything satisfying:

```
POST <tokenEndpoint>  {"password": "..."}   -> 200 {base, sas, expires}
                                            -> 401 wrong or missing
GET  <tokenEndpoint>                        -> same, authorised by a session cookie
```

`base` is the URL prefix the datasets sit under, `sas` a query string appended to every
request for them, `expires` epoch milliseconds — the viewer renews shortly before it.
Leave `tokenEndpoint` empty and the viewer reads plain URLs instead.

A sound implementation keeps the storage credential server-side and hands back
something scoped and short-lived; ties the session cookie to the password, so changing
the password ends existing sessions; and rate-limits attempts.

## Things that are easy to get wrong

Most of these cost real time to find, and none of them announce themselves — the
picture looks plausible either way.

**Build pyramids with max, not mean, for fluorescence.** Averaging 2×2 is right for
brightfield and wrong for sparse signal. A 4×4 object at intensity 60000, four levels
down: mean leaves 938 — under 2% — while max leaves 60000. Markers quietly fade out as
you zoom away, which reads as biology rather than as resampling.

**A density map needs counts over an area.** On a 2 µm grid almost every occupied bin
holds exactly one transcript, so colouring by per-bin value produces a presence mask
with no variation in it. Sum counts into squares large enough to hold several, and take
the scale from a high percentile of the occupied squares rather than the maximum.

**Composite channels with MAX, not addition.** Additive blending makes two markers on
one spot brighter than either, confusing intensity with overlap, and it makes a coarse
tile drawn under its own children flash as they refine. Under MAX an overlapping copy
of the same tile is a no-op, so a parent can stay on screen while children load.

**Take histograms from the coarsest pyramid level.** The obvious index is often the
wrong end. Reading tile (0,0) of the *full-resolution* level samples the top-left corner
of the slide, which is usually blank glass — the histogram then describes background,
and on sparse channels it measures nothing at all.

**Fluorescence sliders need a non-linear axis.** Intensities pile up near zero, so on a
linear 0–255 track every useful window lands in the first few percent of travel. Map
positions through a power curve and draw the histogram on the same axis, so a peak sits
under the handle that windows it.

**Watch the units when you combine scoring terms.** The automatic label placer weighs
local point density against distance from the cluster. Density counts cells per grid
square — hundreds, on 175,000 cells — while distance was a fraction of the plot width,
about three. Distance had no vote at all, and no amount of reweighting helped until both
were scaled to their own range.

## Cluster labels

Labels are placed automatically: off the thick of their own cluster, clear of their
neighbours, inside the frame, and on the outward face of the plot. It is decent, and it
is not a person looking at a figure. `tools/label_editor.html` draws the UMAP, puts each
label on its centroid, and lets you drag them where they belong; it emits coordinates
for `labelPositions` in config. Hand-placed positions win, anything missing is still
placed automatically, and the coordinates are in UMAP units so they hold at any zoom.

## Configuration

`config.js` is the only file that knows about a particular study — its datasets,
institution, logos and byline — and it is gitignored. `config.example.js` documents
every option. The viewer itself names no institution and carries no data.

## Limits

- Everything is in the browser, so a very large gene list or a very large point cloud
  is bounded by the tab's memory. The limits are configurable and conservative.
- Label placement is solved for the fitted view. Zooming keeps labels anchored to their
  coordinates, which is usually what you want, but it is not re-solved as you go.
- Transcript density is rendered as instanced squares, which is fast to tens of
  thousands of bins and heavy beyond that; coarser bins are the answer.
- Built against [deck.gl](https://deck.gl) 9.1, loaded from a CDN. No build step, no
  bundler, no `npm install` — the client is three files.

## License

Not yet chosen. Until a license is added no rights are granted: default copyright
applies and the code may not be copied, used or redistributed.
