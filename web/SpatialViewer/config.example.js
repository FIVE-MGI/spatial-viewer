// Copy to config.js and edit. config.js is gitignored: it names your datasets,
// your institution and your storage, none of which belong in this repository.
//
// The viewer reads nothing else. Everything specific to a study lives here.
window.SpatialViewerConfig = {

  // ---- datasets -----------------------------------------------------------------
  // One entry per folder in the storage container. The key is the folder name.
  //   group         which heading it appears under (see `groups` below)
  //   label         what the picker shows; the tissue and its condition, nothing else
  //   hideChannels  exact channel names to leave out of this dataset, e.g. stain
  //                 tests that are not panel markers. Matched case-insensitively.
  datasets: {
    // example_visium:    { group: "transcriptomics", label: "Normal tissue" },
    // example_proteomics: { group: "proteomics", label: "Condition A",
    //                       hideChannels: ["CD11c FFPE"] },
    // example_rnaseq:    { group: "scrnaseq", label: "All cells" },
  },

  // Headings for the picker, in the order they should appear. A dataset whose
  // `group` is missing from here falls under "Other"; with no groups at all the
  // picker is a flat list.
  groups: {
    transcriptomics: "Spatial transcriptomics",
    proteomics: "Spatial proteomics",
    scrnaseq: "Single-cell RNA-seq",
  },

  defaultDataset: "",

  // ---- whose instance this is ----------------------------------------------------
  // The viewer has no institution of its own. Logo paths are relative to this folder.
  brand: {
    title: "Spatial Viewer",
    logos: [
      // { src: "logos/institution.png", alt: "Institution", className: "logo-wide" },
    ],
    byline: "",
  },

  // Shown as "Link to paper". Empty renders it greyed and inert.
  paperUrl: "",

  // ---- data access ----------------------------------------------------------------
  // The managed function that checks the shared password and returns a short-lived,
  // read-only SAS for the private container. Empty disables the password entirely,
  // which only makes sense when the data is public or served locally.
  tokenEndpoint: "/api/token",

  // A live tile server (server.py) instead of, or as well as, static exports.
  // Empty means none; ?api=http://host at runtime also works.
  apiBase: "",
  useLocalServer: false,

  // ---- UMAP cluster labels --------------------------------------------------------
  // Hand-placed label positions in UMAP coordinates, keyed by dataset then cluster.
  // Anything listed wins over the automatic placement; anything missing is placed by
  // the solver, so a partial list is fine. Produced by tools/label_editor.html.
  labelPositions: {},

  // Clustering columns to color the UMAP by, first match wins. The picker was removed
  // because each dataset here carries exactly one; add names as needed.
  rnaColorBy: ["leiden"],

  // ---- display --------------------------------------------------------------------
  // Most genes at once, so one dataset cannot exhaust the tab.
  rnaseqMaxGenes: 12,   // UMAP panels, one per gene
  maxGenes: 8,          // spatial transcript layers

  // Smallest bin layer offered. Below roughly 8 um a bin holds a single transcript
  // nearly everywhere, which makes a density map of it a presence map.
  exprMinBinUm: 8,

  // Channels never shown, whatever dataset they are in: unassigned detector slots.
  // A JavaScript regex source string; omit to use the built-in default.
  // hideChannels: "^\\s*(open|blank|empty|unused)([\\s_-]|$)",

  tileFormat: "png",    // "jpg" is ~4x smaller per tile (server mode only)
  tileCacheSize: 400,   // decoded tiles held per page, all channels together

  // Window-slider axis. Fluorescence piles up near zero, so a linear 0-255 slider
  // crams every useful window into the first few percent of travel. Positions map
  // through v = max * (p ** gamma); 1 is linear.
  sliderGamma: 3.2,

  // Histogram shape: count -> bar height, then normalised to this percentile rather
  // than the maximum so one background spike cannot flatten the rest.
  histGamma: 0.35,
  histClip: 0.985,

  // Colors cycled through as genes are added.
  geneColors: ["#ff3355", "#33ddff", "#ffdd33", "#66ff66", "#ff88ff", "#ff9933", "#99aaff", "#ffffff"],

  // Floors for categorical colors, as relative luminance, against a dark background.
  // tab20 was drawn for white paper and several of its entries vanish otherwise.
  catMinLuminance: 0.21,    // the points
  labelMinLuminance: 0.5,   // cluster label text, which sits over its own points
};
