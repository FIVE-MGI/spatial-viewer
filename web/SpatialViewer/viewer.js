// SpatialViewer - a thin viewing client over server.py in this folder.
// Images never come to the browser whole: each visible channel is a deck.gl
// TileLayer pulling 8-bit grey tiles, tinted + windowed on the GPU and
// composited with MAX blending like Xenium Explorer / Viv. Transcripts are per-gene
// point clouds fetched as raw float32 and drawn with ScatterplotLayer.
(function () {
  const cfg = window.SpatialViewerConfig || {};
  const params = new URLSearchParams(location.search);
  const API = (params.get("api") || cfg.apiBase || "").replace(/\/$/, "");
  const TILE_FMT = params.get("fmt") || cfg.tileFormat || "png";
  const GENE_COLORS = cfg.geneColors || ["#ff3355", "#33ddff", "#ffdd33", "#66ff66"];
  const NUCLEAR_CHANNEL = /^\s*(dapi|hoechst)\b/i;   // pinned to the top of the list
  const HIDE_CHANNEL = cfg.hideChannels ? new RegExp(cfg.hideChannels, "i")
                                        : /^\s*(open|blank|empty|unused)([\s_-]|$)/i;
  const $ = (id) => document.getElementById(id);

  if (!window.deck) {
    $("status").textContent = "deck.gl failed to load (no internet?)";
    $("status").hidden = false; $("status").classList.add("err");
    return;
  }
  const { Deck, OrthographicView, TileLayer, BitmapLayer, ScatterplotLayer, GeoJsonLayer, TextLayer, ColumnLayer } = deck;

  // ---- channel tile layer: grey texture -> window -> color, on the GPU ------
  const chanModule = {
    name: "chan",
    vs: "uniform chanUniforms { float lo; float hi; vec3 color; } chan;",
    fs: "uniform chanUniforms { float lo; float hi; vec3 color; } chan;",
    uniformTypes: { lo: "f32", hi: "f32", color: "vec3<f32>" },
  };
  class ChannelBitmapLayer extends BitmapLayer {
    getShaders() {
      const s = super.getShaders();
      return {
        ...s,
        modules: [...(s.modules || []), chanModule],
        inject: {
          ...(s.inject || {}),
          "fs:DECKGL_FILTER_COLOR": `
            float g = color.r;
            float v = clamp((g - chan.lo) / max(chan.hi - chan.lo, 1e-4), 0.0, 1.0);
            color = vec4(chan.color * v * layer.opacity, v * layer.opacity);`,
        },
      };
    }
    draw(opts) {
      const { lo, hi, color } = this.props;
      this.setShaderModuleProps({ chan: { lo, hi, color } });
      super.draw(opts);
    }
  }
  ChannelBitmapLayer.layerName = "ChannelBitmapLayer";
  ChannelBitmapLayer.defaultProps = {
    ...BitmapLayer.defaultProps,
    lo: { type: "number", value: 0 },
    hi: { type: "number", value: 1 },
    color: { type: "array", value: [1, 1, 1] },
  };
  // Additive sums whatever overlaps. That is why a coarse tile drawn under its own
  // children flashes bright mid-zoom, and it is also why two markers on the same spot
  // blow out. MAX takes the brighter of the two instead: overlapping copies of the same
  // tile are then a no-op, so the parent can stay on screen while children load - sharp
  // refinement and no flash at once.
  const MAXBLEND = {
    blend: true, depthTest: false,
    blendColorOperation: "max", blendColorSrcFactor: "one", blendColorDstFactor: "one",
    blendAlphaOperation: "max", blendAlphaSrcFactor: "one", blendAlphaDstFactor: "one",
  };
  const NORMAL = {
    blend: true, depthTest: false,
    blendColorOperation: "add", blendColorSrcFactor: "src-alpha", blendColorDstFactor: "one-minus-src-alpha",
    blendAlphaOperation: "add", blendAlphaSrcFactor: "one", blendAlphaDstFactor: "one-minus-src-alpha",
  };

  // ---- state ------------------------------------------------------------------
  const S = {
    datasets: [], dsId: null, m: null, src: null, mode: "image", rna: null,
    channels: [],          // {index,name,color:[r,g,b] 0-1, hex, lo, hi, visible}
    imgOpacity: 1, imgOn: true, geneOpacity: 1,
    exprKey: null, geneIndex: [], genes: [],   // genes: {name,hex,color,max,p99,n,positions,values,colors}
    geneBlend: "heatmap", pointScale: 1, genePalette: "genes", geneBin: 1, heat: null,
    tissue: { on: false, positions: null, key: null },
    overlays: [],          // {key,label,on,data}
    viewState: null,
  };
  const tileCache = new Map();   // "ds/c/z/x/y" -> ImageBitmap (LRU)
  const TILE_CACHE_MAX = cfg.tileCacheSize || 1500;
  // How many genes may be shown at once, so no dataset can exhaust the browser tab.
  // UMAP panels are light (one small panel per gene); a spatial transcript layer draws every
  // bin of that gene, so it is capped more tightly.
  const MAX_RNA_GENES = Math.max(1, Math.min(cfg.rnaseqMaxGenes || 12, 24));
  const MAX_SPATIAL_GENES = Math.max(1, Math.min(cfg.maxGenes || 8, 16));
  let deckgl = null;

  // ---- helpers ----------------------------------------------------------------
  const hexToRgb01 = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);
  const hexToRgb255 = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
  function status(msg, err) {
    const el = $("status");
    if (!msg) { el.hidden = true; return; }
    el.textContent = msg; el.hidden = false; el.classList.toggle("err", !!err);
  }
  function apiUrl(p) { return `${API}/api/${p}`; }

  // ---- data sources -------------------------------------------------------------
  // "static": a folder of files exported by prepare.py export, hosted anywhere.
  //           Normally a PRIVATE blob container: `qs` appends the access token
  //           from /api/token to every request (see below).
  // "api":    the tile server (server.py), password protected.
  function staticSource(base, qs) {
    base = base.replace(/\/?$/, "/");
    const u = (p) => base + p + (qs ? qs() : "");
    return {
      kind: "static",
      manifest: () => getJSON(u("manifest.json"), true),
      tileUrl: (c, z, x, y, fmt) => u(`tiles/${c}/${z}/${x}/${y}.${fmt}`),
      genes: (key) => getJSON(u(`expr/${key}/genes.json`), true),
      gene: (key, name) => getF32(u(`expr/${key}/gene/${encodeURIComponent(name)}.bin`), true),
      positions: (key) => getF32(u(`expr/${key}/positions.bin`), true),
      overlay: (file) => getJSON(u(`overlays/${file}`), true),
      json: (p) => getJSON(u(p), true),
      buf: async (p) => { const r = await dataFetch(u(p)); if (!r.ok) throw new Error(`${r.status} ${p}`); return r.arrayBuffer(); },
    };
  }
  function apiSource(dsId) {
    return {
      kind: "api",
      manifest: () => getJSON(apiUrl(`datasets/${dsId}`)),
      tileUrl: (c, z, x, y, fmt) => apiUrl(`datasets/${dsId}/tiles/${c}/${z}/${x}/${y}.${fmt}`),
      genes: (key) => getJSON(apiUrl(`datasets/${dsId}/expr/${key}/genes`)),
      gene: (key, name) => getF32(apiUrl(`datasets/${dsId}/expr/${key}/gene/${encodeURIComponent(name)}`)),
      positions: (key) => getF32(apiUrl(`datasets/${dsId}/expr/${key}/positions`)),
      overlay: (key) => getJSON(apiUrl(`datasets/${dsId}/overlays/${key}`)),
    };
  }
  // ---- private container: short-lived read-only access token --------------------
  // The exported datasets live in a private blob container, so a URL on its own
  // fetches nothing. /api/token checks the shared password on the server and
  // returns a read-only SAS that expires after a few hours; the storage key never
  // reaches the browser. A session cookie set by that call lets an open tab renew
  // silently, so nobody has to retype the password in the middle of a session.
  const TOKEN_URL = cfg.tokenEndpoint || "";
  let token = null, tokenWait = null;
  const tokenQuery = () => (token ? "?" + token.sas : "");
  async function requestToken(password) {
    const opts = password == null
      ? { method: "GET", cache: "no-store", credentials: "same-origin" }
      : { method: "POST", cache: "no-store", credentials: "same-origin",
          headers: { "Content-Type": "application/json" }, body: JSON.stringify({ password }) };
    const r = await fetch(TOKEN_URL, opts);
    if (r.status === 401) return null;                     // no session yet, or wrong password
    if (r.status === 429) throw new Error("too many password attempts - wait a few minutes");
    if (!r.ok) throw new Error(`${r.status} from ${TOKEN_URL}`);
    return r.json();
  }
  async function ensureToken(force) {
    if (!force && token && Date.now() < token.expires - 60000) return token;
    if (tokenWait) return tokenWait;                       // one request, however many tiles wait on it
    tokenWait = (async () => {
      try {
        let t = await requestToken(null);                  // cookie still valid? then no prompt
        for (let tries = 0; !t && tries < 6; tries++) {
          const pw = await askPassword(tries > 0, "This viewer is password protected.");
          t = await requestToken(pw);
        }
        if (!t) throw new Error("password not accepted");
        token = t;
        return t;
      } finally { tokenWait = null; }
    })();
    return tokenWait;
  }
  // A URL built with a token that expired while it was queued: sign it again.
  const reSign = (url) => (token ? url.split("?")[0] + tokenQuery() : url);
  async function dataFetch(url, opts = {}) {
    let r = await fetch(url, opts);
    if ((r.status === 403 || r.status === 401) && TOKEN_URL && token) {
      await ensureToken(true);
      r = await fetch(reSign(url), opts);
    }
    return r;
  }
  // Renew a little before expiry so a long-open tab never sees a failed tile.
  if (TOKEN_URL) setInterval(() => {
    if (token && Date.now() > token.expires - 5 * 60000) ensureToken(true).catch(() => {});
  }, 60000);

  async function resolveLink(path) {
    // "<id>.link.txt" holds the real URL (same convention as the Nerve page),
    // either plain or encrypted with the shared password ("enc:v1:...", see links.py).
    // Legacy: only used for datasets still configured with `linkPath`.
    if (!/\.link\.txt$/i.test(path)) return path;
    const r = await fetch(path, { cache: "no-cache" });
    if (!r.ok) throw new Error(`${r.status} ${path}`);
    const line = (await r.text()).trim().split(/\r?\n/)[0].trim();
    if (!line.startsWith("enc:v1:")) return line;
    return decryptLink(line);
  }
  // ---- encrypted links: PBKDF2-SHA256 (200k) -> AES-GCM, all in the browser ----
  let linkPassword = null;
  try { linkPassword = sessionStorage.getItem("spatialviewer.pw"); } catch (e) { /* ignore */ }
  const b64 = (str) => Uint8Array.from(atob(str), (ch) => ch.charCodeAt(0));
  async function tryDecrypt(payload, password) {
    const [, , salt, iv, ct] = payload.split(":");
    const base = await crypto.subtle.importKey("raw", new TextEncoder().encode(password), "PBKDF2", false, ["deriveKey"]);
    const key = await crypto.subtle.deriveKey({ name: "PBKDF2", salt: b64(salt), iterations: 200000, hash: "SHA-256" },
      base, { name: "AES-GCM", length: 256 }, false, ["decrypt"]);
    const plain = await crypto.subtle.decrypt({ name: "AES-GCM", iv: b64(iv) }, key, b64(ct));
    return new TextDecoder().decode(plain);
  }
  async function decryptLink(payload) {
    if (!window.isSecureContext || !crypto.subtle) throw new Error("this page needs https (or localhost) to decrypt links");
    for (let tries = 0; tries < 6; tries++) {
      if (linkPassword) {
        try { return await tryDecrypt(payload, linkPassword); } catch (e) { linkPassword = null; }
      }
      const pw = await askPassword(tries > 0, "This dataset is password protected.");
      linkPassword = pw;
      try { sessionStorage.setItem("spatialviewer.pw", pw); } catch (e) { /* ignore */ }
    }
    throw new Error("password not accepted");
  }

  // ---- password (HTTP Basic) ----------------------------------------------------
  // Same-origin: the browser already asked once when the page loaded and sends
  // the credentials itself. Cross-origin (?api=...): browsers never show the
  // login box for fetch(), so we ask here and send the header ourselves.
  let authHeader = null;
  try { const s = sessionStorage.getItem("spatialviewer.auth"); if (s) authHeader = s; } catch (e) { /* ignore */ }
  // One inline login box shared by every request that hits a 401; concurrent
  // requests (tiles) all wait on the same submission.
  let loginWait = null, loginResolve = null;
  // The form must never actually submit (that would put the password in the URL).
  $("login").addEventListener("submit", (e) => {
    e.preventDefault();
    const pw = $("loginPw").value;
    if (!pw) return;
    authHeader = "Basic " + btoa("viewer:" + pw);
    try { sessionStorage.setItem("spatialviewer.auth", authHeader); } catch (err) { /* ignore */ }
    $("login").hidden = true;
    const r = loginResolve; loginWait = null; loginResolve = null;
    if (r) r(pw);
  });
  function askPassword(wrong, text) {
    if (loginWait) return loginWait;
    const form = $("login"), msg = $("loginMsg"), pw = $("loginPw");
    msg.textContent = wrong ? "Wrong password - try again." : (text || `This server (${API || location.origin}) needs a password.`);
    msg.classList.toggle("err", !!wrong);
    form.hidden = false; pw.value = ""; pw.focus();
    loginWait = new Promise((resolve) => { loginResolve = resolve; });
    return loginWait;
  }
  async function apiFetch(url, opts = {}) {
    let r = null;
    for (let tries = 0; tries < 6; tries++) {
      const sent = authHeader;
      r = await fetch(url, { ...opts, headers: { ...(opts.headers || {}), ...(sent ? { Authorization: sent } : {}) } });
      if (r.status !== 401) return r;
      if (authHeader !== sent) continue;           // a password arrived while this was in flight: just retry
      if (sent) {                                  // this request carried the current password, so it is wrong
        authHeader = null;
        try { sessionStorage.removeItem("spatialviewer.auth"); } catch (e) { /* ignore */ }
      }
      await askPassword(!!sent);
    }
    return r;
  }
  async function getJSON(url, plain) {
    const r = plain ? await dataFetch(url) : await apiFetch(url);
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return r.json();
  }
  async function getF32(url, plain) {
    const r = plain ? await dataFetch(url) : await apiFetch(url);
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return new Float32Array(await r.arrayBuffer());
  }
  function percentile(sortedCopy, q) {
    if (!sortedCopy.length) return 1;
    return sortedCopy[Math.min(sortedCopy.length - 1, Math.floor(q * sortedCopy.length))];
  }

  // ---- layers -----------------------------------------------------------------
  function tileUrl(c, z, x, y) {
    const fmt = S.src.kind === "static" ? (S.m.tileFormat || "png") : TILE_FMT;
    return S.src.tileUrl(c, z, x, y, fmt);
  }
  function fetchTile(c, backendZ, x, y, signal) {
    const key = `${S.dsId}/${c}/${backendZ}/${x}/${y}`;
    const hit = tileCache.get(key);
    if (hit) { tileCache.delete(key); tileCache.set(key, hit); return Promise.resolve(hit); }
    const req = S.src.kind === "static" ? dataFetch(tileUrl(c, backendZ, x, y), { signal }) : apiFetch(tileUrl(c, backendZ, x, y), { signal });
    return req
      .then((r) => (r.ok && r.status !== 204 ? r.blob() : null))
      .then((b) => (b ? createImageBitmap(b, { premultiplyAlpha: "none", colorSpaceConversion: "none" }) : null))
      .then((bmp) => {
        if (bmp) {
          tileCache.set(key, bmp);
          // never close() evicted bitmaps: deck may still be drawing them; GC reclaims them
          if (tileCache.size > TILE_CACHE_MAX) tileCache.delete(tileCache.keys().next().value);
        }
        return bmp;
      });
  }

  function channelLayer(ch) {
    const img = S.m.image;
    const maxLevel = img.levels - 1;
    return new TileLayer({
      id: `ch-${S.dsId}-${ch.index}`,
      tileSize: img.tileSize,
      minZoom: -maxLevel,
      maxZoom: 0,
      extent: [0, 0, img.width, img.height],
      maxRequests: 12,
      // Safe under MAX: an overlapping parent tile cannot brighten anything, so a coarse
      // tile can show immediately while its children load.
      refinementStrategy: "best-available",
      // renderSubLayers closes over ch.lo / ch.hi / ch.hex and the global image
      // settings. Without them here deck keeps the sublayers it already built, so
      // moving a window slider changed nothing on screen until you panned.
      updateTriggers: {
        getTileData: [S.dsId],
        renderSubLayers: [ch.lo, ch.hi, ch.hex, S.imgOpacity],
      },
      getTileData: ({ index, signal }) => {
        const z = maxLevel + index.z;
        if (z < 0 || z > maxLevel) return Promise.resolve(null);
        return fetchTile(ch.index, z, index.x, index.y, signal).catch((e) => {
          // An abort is "we moved on", not "this tile is empty". Swallowing it as
          // null marked the tile loaded-with-no-content, so deck stopped drawing the
          // parent tile in its place and never asked for it again: holes that only a
          // pan would fill, and detail that never sharpened past the coarse level.
          // Re-throw so deck discards the tile and re-requests it.
          if (e && (e.name === "AbortError" || signal?.aborted)) throw e;
          return null;                                  // a real error: leave it blank
        });
      },
      onTileError: () => {},
      // Force a draw once the visible tiles are in. On the big slides (8 levels, first
      // open, cold cache) the sublayers were built and loaded but nothing was painted
      // until something else forced a redraw - the dataset opened black until a manual
      // Fit. Costs one redraw per viewport load.
      onViewportLoad: () => { if (deckgl) deckgl.redraw(true); },
      renderSubLayers: (props) => {
        if (!props.data) return null;
        const bb = props.tile.boundingBox;
        return new ChannelBitmapLayer(props, {
          data: null,
          image: props.data,
          bounds: [bb[0][0], bb[1][1], bb[1][0], bb[0][1]],
          lo: ch.lo / 255, hi: ch.hi / 255, color: ch.color,
          opacity: S.imgOpacity,
          parameters: MAXBLEND,
          // Always interpolated. This was a checkbox, but BitmapLayer declares
          // textureParameters as `ignore: true`, so the sampler is fixed when the
          // texture is built and nothing revisits it - the control could only ever
          // be made to work by rebuilding every tile layer to change one sampler.
          textureParameters: { minFilter: "linear", magFilter: "linear" },
        });
      },
    });
  }

  // ---- transcript heat map, the way FIVE Atlas draws one --------------------------
  // One square per bin, colored by the COMBINED density of every gene on show,
  // through a single scientific ramp. Two differences from coloring each gene on
  // its own ramp, which is what this used to do:
  //   - genes add up instead of hiding one another. With a layer per gene the last
  //     one drawn simply covered the others, so a second gene could make the first
  //     disappear rather than reading as more signal.
  //   - a bin holding transcripts never falls off the bottom of the ramp. The
  //     window's low end lifts to 6% rather than 0, so "a little here" stays
  //     visible and only genuinely empty bins are blank.
  // Squares, not discs, and sized to the bin pitch: the bins tile the slide, so
  // circles left gaps that read as structure that is not there.
  function geneHeat() {
    const e = currentExpr();
    const pitch = (e && e.spotDiameterPx) ? e.spotDiameterPx : 8;
    const cell = pitch * Math.max(1, S.geneBin);
    const ink = S.genePalette === "genes";
    const on = S.genes.filter((g) => g.visible !== false);
    const key = on.map((g) => `${g.name}:${g.hex}:${g.min}`).join("|") +
                `#${S.genePalette}#${S.geneBin}#${S.exprKey}#${cell}`;
    if (S.heat && S.heat.key === key) return S.heat;

    const idx = new Map(), KX = 1 << 21;
    const binOf = (g, i) => Math.floor(g.positions[2 * i + 1] / cell) * KX + Math.floor(g.positions[2 * i] / cell);
    for (const g of on) for (let i = 0; i < g.n; i++) {
      const k = binOf(g, i);
      if (!idx.has(k)) idx.set(k, idx.size);
    }
    const n = idx.size;
    const positions = new Float32Array(n * 3), colors = new Uint8Array(n * 4);
    for (const [k, j] of idx) {
      const iy = Math.floor(k / KX), ix = k - iy * KX;
      positions[j * 3] = (ix + 0.5) * cell; positions[j * 3 + 1] = (iy + 0.5) * cell;
    }

    // DENSITY, not presence. Sum the counts falling in each square and scale by how
    // busy the busiest squares are. Normalising each spot against the gene's own
    // per-spot window instead was the bug that made every bin the same color: at
    // 2 um almost every occupied bin holds exactly one UMI, so a 0-1 window sent all
    // of them to the top of the scale. A density has to count over an area, and the
    // scale has to come from the spread of those counts.
    const norm = on.map((g) => {
      const sum = new Float32Array(n);
      const lo = g.min || 0;
      for (let i = 0; i < g.n; i++) {
        const v = g.values[i];
        if (v >= lo) sum[idx.get(binOf(g, i))] += v;       // the window gates spots, not bins
      }
      // top of scale = a high percentile of the occupied squares, so one hot square
      // cannot flatten everything else, and empty squares are excluded
      const nz = Array.from(sum).filter((v) => v > 0).sort((x, y) => x - y);
      const hi = nz.length ? Math.max(nz[Math.floor((nz.length - 1) * 0.99)], 1e-6) : 1;
      for (let j = 0; j < n; j++) sum[j] = Math.min(1, sum[j] / hi);
      return sum;
    });

    if (ink) {
      const paper = new Float32Array(n * 3).fill(1), miss = new Float32Array(n).fill(1);
      on.forEach((g, gi) => {
        const d = norm[gi], [tr, tg, tb] = hexToRgb255(g.hex);
        for (let j = 0; j < n; j++) {
          if (d[j] <= 0) continue;
          // a square with anything in it is at least a fifth of full ink
          const sEff = 0.2 + 0.8 * d[j];
          paper[j * 3] *= 1 - sEff * (1 - tr / 255);
          paper[j * 3 + 1] *= 1 - sEff * (1 - tg / 255);
          paper[j * 3 + 2] *= 1 - sEff * (1 - tb / 255);
          miss[j] *= 1 - sEff;
        }
      });
      for (let j = 0; j < n; j++) {
        const c = j * 4;
        colors[c] = Math.round(paper[j * 3] * 255); colors[c + 1] = Math.round(paper[j * 3 + 1] * 255);
        colors[c + 2] = Math.round(paper[j * 3 + 2] * 255);
        colors[c + 3] = Math.round(Math.sqrt(Math.max(0, 1 - miss[j])) * 255);
      }
    } else {
      const lut = PALETTES[S.genePalette] || PALETTES.viridis;
      // The AVERAGE of the genes on show, not the sum. Summing meant two genes at
      // half density each already hit the top of the ramp, so adding a second gene
      // flooded the map and every dense area looked identical. An average keeps the
      // ramp meaning the same thing however many genes are on.
      const inv = 1 / Math.max(1, norm.length);
      for (let j = 0; j < n; j++) {
        let total = 0, any = false;
        for (const d of norm) { total += d[j]; if (d[j] > 0) any = true; }
        const t = Math.min(255, Math.round((0.06 + 0.94 * Math.min(1, total * inv)) * 255)) * 3;
        const c = j * 4;
        colors[c] = lut[t]; colors[c + 1] = lut[t + 1]; colors[c + 2] = lut[t + 2];
        colors[c + 3] = any ? 235 : 0;
      }
    }
    S.heat = { key, n, positions, colors, cell, index: Array.from({ length: n }, (_, j) => j) };
    return S.heat;
  }


  // With one gene the picture explains itself. With three it does not: nothing on
  // screen says whether a dark square means "dense" or "two genes overlapping", or
  // what the ramp is averaged over. This says so in a sentence.
  function drawDensityKey() {
    const el = $("densityKey");
    const on = S.genes.filter((g) => g.visible !== false);
    if (S.mode !== "image" || S.geneBlend !== "heatmap" || !on.length) { el.hidden = true; return; }
    const e = currentExpr();
    const um = (e && e.binSizeUm) ? `${e.binSizeUm.toLocaleString()} µm` : "each";
    if (S.genePalette === "genes") {
      el.innerHTML =
        `<div class="genes">${on.map((g) => `<span class="g"><i style="background:${inkCss(g.hex)}"></i>${g.name}</span>`).join("")}</div>` +
        `<div class="note">pale = few transcripts, saturated = many, per ${um} square.` +
        (on.length > 1 ? " Where genes overlap the colors mix and darken." : "") + `</div>`;
    } else {
      el.innerHTML =
        `<div class="ramp"><span>low</span><i style="background:${rampCss(S.genePalette)}"></i><span>high</span></div>` +
        `<div class="note">` +
        (on.length > 1
          ? `transcripts per ${um} square, averaged over ${on.length} genes (${on.map((g) => g.name).join(", ")})`
          : `${on[0].name} transcripts per ${um} square`) + `</div>`;
    }
    el.hidden = false;
  }

  function geneHeatLayer() {
    const h = geneHeat();
    if (!h.n) return null;
    // A 2 um bin is under a tenth of a screen pixel when the whole slide is in
    // view, so the squares vanished entirely at the fit zoom. ScatterplotLayer
    // has radiusMinPixels for exactly this; ColumnLayer has no equivalent, so
    // hold a floor in screen pixels by hand. zoom is log2(pixels per unit).
    const ppu = Math.pow(2, (S.viewState && S.viewState.zoom) || 0);
    const side = Math.max(h.cell * S.pointScale, HEAT_MIN_PX / ppu);
    return new ColumnLayer({
      id: `geneheat-${S.dsId}-${S.exprKey}`,
      // ColumnLayer does not take the binary-attribute form ScatterplotLayer does;
      // fed that way it draws nothing at all, silently. An index array with plain
      // accessors reading out of the same typed arrays costs no extra memory.
      data: h.index,
      getPosition: (i) => [h.positions[i * 3], h.positions[i * 3 + 1], 0],
      getFillColor: (i) => [h.colors[i * 4], h.colors[i * 4 + 1], h.colors[i * 4 + 2], h.colors[i * 4 + 3]],
      diskResolution: 4, angle: 45,            // four sides, corner-up: an axis-aligned square
      radius: side / Math.SQRT2,               // circumradius of a square of side `side`
      radiusUnits: "common", extruded: false, filled: true, stroked: false,
      opacity: S.geneOpacity, parameters: NORMAL,
      updateTriggers: { getFillColor: [h.key], getPosition: [h.key] },
    });
  }

  function geneLayer(g) {
    const e = currentExpr();
    const r = Math.max(0.5, (e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4) * S.pointScale);
    return new ScatterplotLayer({
      id: `gene-${S.dsId}-${S.exprKey}-${g.name}`,
      data: { length: g.n, attributes: { getPosition: { value: g.positions, size: 2 }, getFillColor: { value: g.colors, size: 4, normalized: true } } },
      dataComparator: () => false,
      radiusUnits: "common", getRadius: r, radiusMinPixels: 0.6,
      stroked: false, filled: true,
      opacity: S.geneOpacity,
      parameters: NORMAL,
      updateTriggers: { getFillColor: [g.colorKey] },
    });
  }

  function tissueLayer() {
    const t = S.tissue;
    if (!t.on || !t.positions) return null;
    const e = currentExpr();
    const r = Math.max(0.5, (e && e.spotDiameterPx ? e.spotDiameterPx / 2 : 4) * S.pointScale * 0.9);
    return new ScatterplotLayer({
      id: `tissue-${S.dsId}-${t.key}`,
      data: { length: t.positions.length / 2, attributes: { getPosition: { value: t.positions, size: 2 } } },
      radiusUnits: "common", getRadius: r, radiusMinPixels: 0.5,
      getFillColor: [120, 130, 150, 90], stroked: false,
      parameters: NORMAL,
    });
  }

  function overlayLayer(o) {
    if (!o.on || !o.data) return null;
    return new GeoJsonLayer({
      id: `ov-${S.dsId}-${o.key}`, data: o.data,
      stroked: true, filled: true, getFillColor: [255, 255, 255, 18],
      getLineColor: (f) => {
        const c = f.properties && f.properties.classification && f.properties.classification.color;
        return Array.isArray(c) ? [c[0], c[1], c[2], 220] : [255, 220, 80, 220];
      },
      lineWidthUnits: "pixels", getLineWidth: 1.5, pointRadiusMinPixels: 2,
      pickable: true, parameters: NORMAL,
    });
  }

  let renderPending = false;  // asked to draw while the pane had no size (hidden tab)
  // Dragging a slider fires `input` continuously, and render() rebuilds layers - on a
  // multi-channel slide that is a tile layer and all its per-tile sublayers, per pixel of
  // travel. Coalesce to at most one render per animation frame: the work still happens,
  // just once per frame the screen can actually show.
  let renderQueued = false;
  function requestRender() {
    if (renderQueued) return;
    renderQueued = true;
    requestAnimationFrame(() => { renderQueued = false; render(); });
  }

  function render() {
    if (!deckgl || !S.m) return;
    const el = $("viewport");
    if (!el.clientWidth || !el.clientHeight) { renderPending = true; return; }   // deck can't build a viewport at 0 px
    renderPending = false;
    if (S.mode === "rna") { renderRna(); return; }
    if (!S.viewState || !Number.isFinite(S.viewState.zoom) || fitPending) fitView();
    const layers = [];
    if (S.imgOn) for (const ch of S.channels) if (ch.visible) layers.push(channelLayer(ch));
    const tl = tissueLayer(); if (tl) layers.push(tl);
    if (S.geneBlend === "heatmap") {
      const hl = geneHeatLayer();
      if (hl) layers.push(hl);
    } else {
      for (const g of S.genes) if (g.visible !== false) layers.push(geneLayer(g));
    }
    drawDensityKey();
    for (const o of S.overlays) { const l = overlayLayer(o); if (l) layers.push(l); }
    deckgl.setProps({ layers });
  }

  // ---- gene colors -------------------------------------------------------------
  function recolorGene(g) {
    const n = g.n, out = new Uint8Array(n * 4);
    const [r, gg, b] = hexToRgb255(g.hex);
    const lo = g.min || 0, inv = 1 / Math.max(g.max - lo, 1e-6);
    const mode = S.geneBlend;
    for (let i = 0; i < n; i++) {
      const raw = g.values[i];
      const k = i * 4;
      if (raw < lo) { out[k + 3] = 0; continue; }                       // below min display: hidden
      const v = Math.min(1, (raw - lo) * inv);
      if (mode === "heatmap") {
        // Viridis, opaque - the same ramp the RNA-seq panels use, so a gene reads the
        // same way in both views. Low is dark blue, high is yellow, and unlike the glow
        // the brightness of a bin is not confounded with how many bins overlap it.
        const t = Math.min(255, Math.round(v * 255)) * 3;
        out[k] = VIRIDIS[t]; out[k + 1] = VIRIDIS[t + 1]; out[k + 2] = VIRIDIS[t + 2]; out[k + 3] = 255;
      } else {
        out[k] = r; out[k + 1] = gg; out[k + 2] = b; out[k + 3] = 40 + 215 * v;
      }
    }
    g.colors = out;
    g.colorKey = `${g.hex}|${lo}|${g.max}|${S.geneBlend}`;
    // Only the low gate decides which spots the density counts, so only a change to
    // it need throw the bins away. Clearing on every window change meant dragging
    // the upper handle recomputed 300k+ squares, and a percentile sort, per pixel of
    // drag - to redraw exactly the same picture.
    if (g.heatKey !== `${g.hex}|${lo}`) { g.heatKey = `${g.hex}|${lo}`; S.heat = null; }
  }

  // ---- view -----------------------------------------------------------------------
  let fitPending = false;   // the pane was 0 px when we fitted (hidden tab): refit on resize
  function fitView() {
    if (!S.m) return;
    if (S.mode === "rna") { fitRna(); return; }
    const img = S.m.image;
    const el = $("viewport");
    const w = el.clientWidth || 0, h = el.clientHeight || 0;
    let zoom = Math.log2(Math.min(w / img.width, h / img.height)) - 0.03;
    if (!Number.isFinite(zoom)) { zoom = Math.log2(1000 / Math.max(img.width, img.height)); fitPending = true; }
    else fitPending = false;
    S.viewState = { target: [img.width / 2, img.height / 2, 0], zoom, minZoom: zoom - 2, maxZoom: 6 };
    deckgl.setProps({ viewState: S.viewState });
    updateHud();
  }
  function updateHud(xy) {
    if (S.mode === "rna") return;
    const vs = S.viewState; if (!vs) return;
    $("hudZoom").textContent = `zoom ${vs.zoom.toFixed(2)}  (${(Math.pow(2, vs.zoom)).toFixed(3)} px/px)`;
    if (xy) {
      const um = S.m.pixelSizeUm;
      $("hudXY").textContent = `x ${xy[0].toFixed(0)}  y ${xy[1].toFixed(0)} px` +
        (um ? `   (${(xy[0] * um).toFixed(1)}, ${(xy[1] * um).toFixed(1)} µm)` : "");
    }
    // scale bar
    const sb = $("scalebar");
    const um = S.m.pixelSizeUm;
    const unitsPerPx = 1 / Math.pow(2, vs.zoom);              // image px per screen px
    const worldPerScreen = um ? unitsPerPx * um : unitsPerPx; // µm (or px) per screen px
    const nice = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000];
    let len = nice.find((n) => n / worldPerScreen >= 70) || nice[nice.length - 1];
    sb.querySelector(".bar").style.width = `${len / worldPerScreen}px`;
    sb.querySelector(".lbl").textContent = um ? (len >= 1000 ? `${len / 1000} mm` : `${len} µm`) : `${len} px`;
  }

  function imageViews() {
    return [new OrthographicView({ id: "ortho", flipY: true, controller: { inertia: 200, scrollZoom: { speed: 0.01, smooth: false }, doubleClickZoom: true, keyboard: false } })];
  }
  function initDeck() {
    deckgl = new Deck({
      parent: $("viewport"),
      views: imageViews(),
      // UMAP mode: "umap-base" is drawn in every panel (except under the color-by panel's own
      // colors); "umap-<i>-..." layers only in panel i. Image-mode layers are unaffected.
      layerFilter: ({ layer, viewport }) => {
        const id = layer.id;
        if (!id.startsWith("umap-")) return true;
        const m = id.match(/^umap-(\d+)-/);
        if (m) return viewport.id === `umap-${m[1]}`;
        return !(viewport.id === "umap-0" && S.rna && S.rna.obsData);
      },
      viewState: { target: [0, 0, 0], zoom: 0 },
      layers: [],
      useDevicePixels: true,
      // keep the last frame in the buffer so screenshots / "save image" see it
      deviceProps: { webgl: { preserveDrawingBuffer: true } },
      style: { background: "#000" },
      getCursor: ({ isDragging }) => (isDragging ? "grabbing" : "crosshair"),
      onViewStateChange: ({ viewState }) => {
        if (S.mode === "rna" && S.rna) { S.rna.vs = viewState; deckgl.setProps({ viewState: rnaViewStates() }); return; }
        const zoomed = !S.viewState || Math.abs(S.viewState.zoom - viewState.zoom) > 1e-6;
        S.viewState = viewState; deckgl.setProps({ viewState }); updateHud();
        // the heat map's squares have a screen-space floor, so their size in world
        // units changes with the zoom
        if (zoomed && S.geneBlend === "heatmap" && S.genes.length) requestRender();
      },
      onHover: (info) => { if (info.coordinate) updateHud(info.coordinate); },
      getTooltip: (info) => {
        if (!info.object || !info.object.properties) return null;
        const p = info.object.properties;
        const name = p.name || (p.classification && p.classification.name) || p.label || p.id;
        return name ? { text: String(name) } : null;
      },
    });
    window.addEventListener("keydown", (e) => { if (e.key === "r" || e.key === "R") fitView(); });
    window.__spatialViewer = { deck: deckgl, state: S, render, fitView };   // for debugging from the console
  }

  // ---- dual-handle range slider (QuPath "min display / max display") ------------------
  // Two native range inputs stacked on one track; the band between the handles is drawn
  // in the layer's color; an optional histogram sits behind the track. `fmt` turns the
  // 0..max slider units into the label the user should read (raw intensity, counts...).
  // Fluorescence intensities pile up near zero: on a linear 0-255 axis every useful
  // window lands in the first few percent of travel, so the handles are unusable
  // exactly where the detail is. The slider therefore works in position units and
  // maps them through a power curve, stretching the dark end. The histogram is drawn
  // on the same axis, so a peak sits directly under the handle that windows it.
  const SLIDER_GAMMA = cfg.sliderGamma || 2.5;
  const HIST_GAMMA = cfg.histGamma || 0.35;        // count -> bar height; 1 = linear, lower = more shape
  const HIST_CLIP = cfg.histClip || 0.985;         // normalise to this percentile, not the max
  const EXPR_MIN_UM = cfg.exprMinBinUm || 8;        // smallest bin layer offered
  const HEAT_MIN_PX = cfg.heatMinPixels || 1.5;    // smallest a heat-map square may draw
  const HIST_BAR = cfg.histBarWidth || 3;          // px, on the 256-wide histogram canvas
  const HIST_GAP = cfg.histBarGap || 1;
  const POS_MAX = 1000;
  const posToVal = (p, min, span) => min + span * Math.pow(Math.min(1, Math.max(0, p / POS_MAX)), SLIDER_GAMMA);
  const valToPos = (v, min, span) => POS_MAX * Math.pow(Math.min(1, Math.max(0, (v - min) / span)), 1 / SLIDER_GAMMA);

  function rangeSlider(host, o) {
    host.classList.add("rs");
    const span = o.max - o.min || 1;
    const toVal = (p) => Math.round(posToVal(+p, o.min, span));
    const toPos = (v) => Math.round(valToPos(v, o.min, span));
    host.innerHTML = `<canvas class="rs-hist"></canvas><div class="rs-track"></div><div class="rs-band"></div>
      <input class="rs-lo" type="range" min="0" max="${POS_MAX}" step="1" value="${toPos(o.lo)}" title="min display" />
      <input class="rs-hi" type="range" min="0" max="${POS_MAX}" step="1" value="${toPos(o.hi)}" title="max display" />
      <span class="rs-lbl rs-lbl-lo"></span><span class="rs-lbl rs-lbl-hi"></span>`;
    const lo = host.querySelector(".rs-lo"), hi = host.querySelector(".rs-hi"), band = host.querySelector(".rs-band");
    const llo = host.querySelector(".rs-lbl-lo"), lhi = host.querySelector(".rs-lbl-hi");
    const hist = host.querySelector(".rs-hist");

    const paintHist = () => {
      const cols = hist.__cols;
      if (!cols) return;
      const W = hist.width, H = hist.height;
      const g = hist.getContext("2d");
      g.clearRect(0, 0, W, H);
      const a = (+lo.value / POS_MAX) * W, b = (+hi.value / POS_MAX) * W;
      // Drawn as separated bars rather than a filled curve: 1-px columns butted
      // together read as one smooth shape, which is what made a histogram with
      // plenty of structure in it still look like a flat line.
      const step = HIST_BAR + HIST_GAP;
      for (let x0 = 0; x0 < W; x0 += step) {
        let v = 0;
        for (let x = x0; x < Math.min(W, x0 + step); x++) if (cols[x] > v) v = cols[x];
        const bh = Math.round(v * (H - 1));
        if (bh <= 0) continue;
        const mid = x0 + step / 2;
        g.fillStyle = (mid >= a && mid <= b) ? o.color() : "rgba(255,255,255,0.20)";
        g.fillRect(x0, H - bh, HIST_BAR, bh);
      }
    };
    const paint = () => {
      const a = +lo.value / POS_MAX, b = +hi.value / POS_MAX;
      band.style.left = `${a * 100}%`; band.style.width = `${Math.max(0, b - a) * 100}%`;
      band.style.background = o.color();
      llo.textContent = o.fmt(toVal(lo.value)); lhi.textContent = o.fmt(toVal(hi.value));
      paintHist();
    };
    // keep at least one raw unit between the handles, in value space not position space
    lo.oninput = () => { if (toVal(lo.value) >= toVal(hi.value)) lo.value = Math.max(0, +hi.value - 1); paint(); o.onChange(toVal(lo.value), toVal(hi.value)); };
    hi.oninput = () => { if (toVal(hi.value) <= toVal(lo.value)) hi.value = Math.min(POS_MAX, +lo.value + 1); paint(); o.onChange(toVal(lo.value), toVal(hi.value)); };
    // handle on top = the one nearer the pointer, so both are reachable when they overlap
    host.onpointermove = (e) => {
      const r = host.getBoundingClientRect(); const f = (e.clientX - r.left) / r.width;
      const a = +lo.value / POS_MAX, b = +hi.value / POS_MAX;
      const nearLo = Math.abs(f - a) < Math.abs(f - b);
      lo.style.zIndex = nearLo ? 4 : 3; hi.style.zIndex = nearLo ? 3 : 4;
    };
    paint();
    return {
      set: (a, b) => { lo.value = toPos(a); hi.value = toPos(b); paint(); },
      repaint: paint, hist, toPos,
    };
  }
  // 64-bin histogram of one 8-bit channel from its coarsest tile (background pixels = 0 dropped)
  async function drawHistogram(canvas, ch, slider) {
    try {
      // The coarsest level, not level 0. Level 0 is full resolution, so this used to
      // pull one full-res tile per channel on open - 33 of them on the big Akoya
      // slides - which saturated the connection pool and starved the tiles actually
      // being drawn. Worse, tile (0,0) at full res is the top-left corner, which is
      // usually blank glass, so the histogram described background rather than tissue.
      // The coarsest level is one small tile covering the whole slide.
      // z 0 is the coarsest level here: getTileData maps deck's index.z through
      // `maxLevel + index.z`, so `levels - 1` is FULL resolution. Asking for it by
      // that name fetched tile (0,0) of the largest level - the top-left corner of
      // the slide, which is blank glass - so most channels measured nothing and drew
      // a flat line. One small tile covering the whole slide is what was wanted.
      const bmp = await fetchTile(ch.index, 0, 0, 0);
      if (!bmp) return;
      // Sample the tile at its own resolution. Squeezing a 1024-px tile into a
      // 256x128 thumbnail averaged neighbouring pixels together before counting
      // them, which is precisely what destroys a distribution: ~32k blurred samples
      // with the tails averaged off, so every channel drew the same smooth hump.
      const w = Math.min(bmp.width || 256, 1024), h = Math.min(bmp.height || 256, 1024);
      const oc = document.createElement("canvas"); oc.width = w; oc.height = h;
      const cx = oc.getContext("2d", { willReadFrequently: true });
      cx.imageSmoothingEnabled = false; cx.drawImage(bmp, 0, 0, w, h);
      const d = cx.getImageData(0, 0, w, h).data;
      const raw = new Float32Array(256);                      // one bin per 8-bit value
      for (let i = 0; i < d.length; i += 4) { const v = d[i]; if (v > 0) raw[v] += 1; }

      // Spread each value bin across the columns it occupies on the stretched axis and
      // divide by that width, so the bars show density rather than piling every dark
      // value into one spike. Without this the curve would misrepresent the low end
      // exactly where the axis gives it the most room.
      const W = 256, cols = new Float32Array(W);
      const colOf = (v) => (valToPos(v, 0, 255) / POS_MAX) * W;
      for (let v = 1; v < 256; v++) {
        if (!raw[v]) continue;
        const c0 = Math.max(0, Math.floor(colOf(v)));
        const c1 = Math.min(W, Math.max(c0 + 1, Math.ceil(colOf(v + 1))));
        const share = raw[v] / (c1 - c0);
        for (let c = c0; c < c1; c++) cols[c] += share;
      }
      // Counts run over orders of magnitude - a background bin can hold 1e5 pixels where
      // real signal holds 1e2. log1p squashed that to within a factor of two, so every
      // channel drew the same lazy hump. A gentler power compression keeps the shape,
      // and normalising to a high percentile rather than the maximum stops one
      // background spike from flattening everything else against the floor.
      for (let x = 0; x < W; x++) cols[x] = Math.pow(cols[x], HIST_GAMMA);
      const sorted = Array.from(cols).filter((v) => v > 0).sort((a, b) => a - b);
      let top = sorted.length ? sorted[Math.floor((sorted.length - 1) * HIST_CLIP)] : 0;
      if (!(top > 0)) for (const v of cols) if (v > top) top = v;
      if (top > 0) for (let x = 0; x < W; x++) cols[x] = Math.min(1, cols[x] / top);

      canvas.width = W; canvas.height = 34;
      canvas.__cols = cols;
      if (slider) slider.repaint();                            // draws it, tinted by the window
    } catch (e) { /* no histogram is fine */ }
  }
  function rawLabel(ch, v255) {
    // slider units are 0..255 display; show the intensity in the file's own units
    const bits = (S.m.image.channels[ch.index] || {}).bits || [0, 255];
    const raw = bits[0] + (v255 / 255) * (bits[1] - bits[0]);
    return bits[1] > 255 ? String(Math.round(raw)) : String(Math.round(v255));
  }

  // ---- UI: channels ---------------------------------------------------------------
  function syncImageButton() {
    const btn = $("btnImage"); if (!btn) return;
    btn.textContent = S.imgOn ? "Image off" : "Image on";
    btn.title = S.imgOn ? "Hide the image, keeping your channel settings" : "Show the image again";
    $("imgOpacity").disabled = !S.imgOn;
  }

  function filterChannelUI() {
    const q = ($("chanSearch").value || "").trim().toLowerCase();
    let shown = 0;
    for (const row of $("chanList").children) {
      const hit = !q || (row.dataset.name || "").includes(q);
      row.hidden = !hit;
      if (hit) shown += 1;
    }
    $("chanCount").textContent = q ? `(${shown} / ${S.channels.length})` : `(${S.channels.length})`;
  }
  function buildChannelUI() {
    const list = $("chanList"); list.innerHTML = "";
    $("chanCount").textContent = `(${S.channels.length})`;
    for (const ch of S.channels) {
      const row = document.createElement("div");
      row.className = "chan" + (ch.visible ? "" : " off");
      row.innerHTML = `
        <input class="cb" type="checkbox" ${ch.visible ? "checked" : ""} />
        <input class="col" type="color" value="${ch.hex}" />
        <span class="name" title="${ch.name}">${ch.name}</span>
        <button class="x reset" title="reset to the default display range">reset</button>
        <div class="range"></div>`;
      row.dataset.name = ch.name.toLowerCase();
      const cb = row.querySelector(".cb"), col = row.querySelector(".col");
      const slider = rangeSlider(row.querySelector(".range"), {
        min: 0, max: 255, step: 1, lo: ch.lo, hi: ch.hi,
        color: () => ch.hex, fmt: (v) => rawLabel(ch, v),
        onChange: (a, b) => { ch.lo = a; ch.hi = b; requestRender(); },
      });
      drawHistogram(slider.hist, ch, slider);
      cb.onchange = () => { ch.visible = cb.checked; row.classList.toggle("off", !ch.visible); render(); };
      col.oninput = () => { ch.hex = col.value; ch.color = hexToRgb01(col.value); slider.repaint(); requestRender(); };
      row.querySelector(".reset").onclick = () => {
        const w = S.m.image.channels[ch.index].window; ch.lo = Math.round(w[0]); ch.hi = Math.round(w[1]);
        slider.set(ch.lo, ch.hi); render();
      };
      row.querySelector(".name").onclick = () => { cb.checked = !cb.checked; cb.onchange(); };
      list.appendChild(row);
    }
    filterChannelUI();
  }

  // ---- UI: genes ------------------------------------------------------------------
  function currentExpr() {
    return (S.m.expr || []).find((e) => e.key === S.exprKey) || null;
  }
  async function loadGeneIndex() {
    S.geneIndex = [];
    if (!S.exprKey) return;
    status("Loading gene list…");
    try {
      const list = await S.src.genes(S.exprKey);
      list.sort((a, b) => b.t - a.t);
      S.geneIndex = list;
    } finally { status(null); }
  }
  function showGeneResults(q) {
    const box = $("geneResults");
    q = q.trim().toLowerCase();
    if (!q) { box.hidden = true; return; }
    const hits = [];
    const exact = S.geneIndex.find((g) => g.n.toLowerCase() === q);
    if (exact) hits.push(exact);
    for (const g of S.geneIndex) {
      const n = g.n.toLowerCase();
      if (g !== exact && n.startsWith(q)) hits.push(g);
      if (hits.length >= 40) break;
    }
    if (hits.length < 40) for (const g of S.geneIndex) {
      const n = g.n.toLowerCase();
      if (!n.startsWith(q) && n.includes(q)) hits.push(g);
      if (hits.length >= 40) break;
    }
    box.innerHTML = "";
    hits.forEach((g, i) => {
      const d = document.createElement("div");
      d.innerHTML = `<span>${g.n}</span><small>${Math.round(g.t).toLocaleString()} UMI</small>`;
      if (i === 0) d.classList.add("hot");
      d.onmousedown = (e) => { e.preventDefault(); addGene(g.n); };
      box.appendChild(d);
    });
    box.hidden = hits.length === 0;
  }
  async function addGene(name) {
    if (!S.exprKey || S.genes.some((g) => g.name === name)) return;
    $("geneResults").hidden = true; $("geneSearch").value = "";
    if (S.genes.length >= MAX_SPATIAL_GENES) {
      status(`Up to ${MAX_SPATIAL_GENES} genes at once. Remove one to add ${name}.`, true);
      setTimeout(() => status(null), 4000);
      return;
    }
    status(`Loading ${name}…`);
    try {
      const xyv = await S.src.gene(S.exprKey, name);
      const n = xyv.length / 3;
      const positions = new Float32Array(n * 2), values = new Float32Array(n);
      for (let i = 0; i < n; i++) { positions[2 * i] = xyv[3 * i]; positions[2 * i + 1] = xyv[3 * i + 1]; values[i] = xyv[3 * i + 2]; }
      const sorted = Float32Array.from(values).sort();
      const p99 = Math.max(1, percentile(sorted, 0.99));
      const g = { name, hex: GENE_COLORS[S.genes.length % GENE_COLORS.length], n, positions, values,
                  min: 0, max: p99, p99, vmax: sorted[sorted.length - 1] || 1, visible: true };
      recolorGene(g);
      S.genes.push(g);
      buildGeneUI(); render();
    } catch (e) { status(`Could not load ${name}: ${e.message}`, true); return; }
    status(null);
  }
  function buildGeneUI() {
    $("geneCount").textContent = `(${S.genes.length} / ${MAX_SPATIAL_GENES})`;
    const list = $("geneList"); list.innerHTML = "";
    for (const g of S.genes) {
      const row = document.createElement("div");
      row.className = "chan";
      const top = Math.max(2, Math.ceil(Math.min(g.vmax, g.p99 * 3)));
      row.innerHTML = `
        <input class="cb" type="checkbox" ${g.visible === false ? "" : "checked"} />
        <input class="col" type="color" value="${g.hex}" />
        <span class="name" title="${g.name}">${g.name}</span>
        <button class="rm x" title="remove">×</button>
        <div class="range"></div>`;
      row.querySelector(".x").onclick = () => { S.genes = S.genes.filter((x) => x !== g); buildGeneUI(); render(); };
      const slider = rangeSlider(row.querySelector(".range"), {
        min: 0, max: top, step: top > 50 ? 1 : 0.5, lo: g.min, hi: g.max,
        color: () => g.hex, fmt: (v) => `${v} UMI`,
        onChange: (a, b) => { g.min = a; g.max = b; recolorGene(g); requestRender(); },
      });
      row.querySelector(".col").oninput = (e) => { g.hex = e.target.value; slider.repaint(); recolorGene(g); requestRender(); };
      const cb = row.querySelector(".cb");
      cb.title = `${g.n.toLocaleString()} bins > 0`;
      cb.onchange = () => {
        g.visible = cb.checked;
        row.classList.toggle("off", !cb.checked);
        S.heat = null;                  // the density view sums only what is on
        render();
      };
      row.classList.toggle("off", g.visible === false);
      list.appendChild(row);
    }
  }
  async function setExprKey(key) {
    S.exprKey = key || null;
    S.genes = []; S.tissue.positions = null; S.heat = null; S.tissue.key = key;
    // One square per bin, always: the layer already sets the resolution, and a
    // second multiplier over it only invited the picker and the squares to disagree.
    S.geneBin = 1;
    buildGeneUI();
    await loadGeneIndex();
    if (S.tissue.on) await loadTissue();
    render();
  }
  async function loadTissue() {
    if (!S.exprKey) return;
    status("Loading bin positions…");
    try { S.tissue.positions = await S.src.positions(S.exprKey); }
    finally { status(null); }
  }

  // ---- UI: overlays -----------------------------------------------------------------
  function buildOverlayUI() {
    const list = $("overlayList"); list.innerHTML = "";
    $("panelOverlays").hidden = S.overlays.length === 0;
    for (const o of S.overlays) {
      const row = document.createElement("label");
      row.className = "row check";
      row.innerHTML = `<input type="checkbox" ${o.on ? "checked" : ""}/><span>${o.label}</span>`;
      row.querySelector("input").onchange = async (e) => {
        o.on = e.target.checked;
        if (o.on && !o.data) {
          status(`Loading ${o.label}…`);
          try { o.data = await S.src.overlay(S.src.kind === "static" ? o.file : o.key); } finally { status(null); }
        }
        render();
      };
      list.appendChild(row);
    }
  }

  // ---- dataset ------------------------------------------------------------------------

  // ---- RNA-seq: one linked UMAP panel per gene --------------------------------------
  // All panels share one WebGL context (one Deck, several views) and one pan/zoom state,
  // and the number of gene panels is capped so large datasets cannot exhaust memory.
  const CAT_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#bcbd22",
    "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5",
    "#7f7f7f", "#393b79"];
  const catColor = (i) => (i < CAT_COLORS.length ? hexToRgb255(CAT_COLORS[i]) : hslRgb((i * 137.508) % 360, 0.65, 0.55));
  function hslRgb(h, s, l) {
    const a = s * Math.min(l, 1 - l), f = (n) => { const k = (n + h / 30) % 12; return Math.round(255 * (l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1)))); };
    return [f(0), f(8), f(4)];
  }
  // tab20 was drawn for white paper. On this background #1f77b4 (blue), #8c564b
  // (brown) and #393b79 (navy) sink into it, and a cluster you cannot see is one
  // whose label you cannot read either. Lift anything under the floor toward its
  // own lighter self, keeping the hue and most of the saturation so the legend
  // swatch, the points and the label stay recognisably the same color.
  const CAT_MIN_Y = cfg.catMinLuminance || 0.21;
  const LABEL_MIN_Y = cfg.labelMinLuminance || 0.5;
  const _lin = (c) => { c /= 255; return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4); };
  const _lum = ([r, g, b]) => 0.2126 * _lin(r) + 0.7152 * _lin(g) + 0.0722 * _lin(b);
  function liftTo(rgb, target) {
    // Relative luminance, not HSL lightness: by lightness pure orange scores 0.53 and
    // would get "brightened" into a duller, darker orange. Anything already bright
    // enough is left exactly alone; the rest is mixed toward white by the smallest
    // amount that clears the floor, which keeps the hue.
    if (_lum(rgb) >= target) return rgb;
    let lo = 0, hi = 0.95;
    for (let i = 0; i < 14; i++) {
      const t = (lo + hi) / 2;
      if (_lum(rgb.map((v) => v + (255 - v) * t)) >= target) hi = t; else lo = t;
    }
    return rgb.map((v) => Math.round(v + (255 - v) * hi));
  }
  const liftForDark = (rgb) => liftTo(rgb, CAT_MIN_Y);
  // Label text sits over its own cluster's points, so clearing the background is not
  // enough - it has to clear the color it is sitting on. A higher floor than the
  // points themselves, and the halo below does the rest.
  const labelColor = (rgb) => liftTo(rgb, LABEL_MIN_Y);


  // Viridis, 16 anchors interpolated to 256. FIVE Atlas offers five maps; one is
  // enough here, and viridis is the one to keep - perceptually uniform and readable
  // to color-blind viewers, which turbo and the rest are not.
  const RAMP = (stops) => {
    const lut = new Uint8Array(256 * 3);
    for (let i = 0; i < 256; i++) {
      const t = (i / 255) * (stops.length - 1), a = Math.floor(t), b = Math.min(stops.length - 1, a + 1), f = t - a;
      for (let c = 0; c < 3; c++) lut[i * 3 + c] = Math.round(stops[a][c] * (1 - f) + stops[b][c] * f);
    }
    return lut;
  };
  const PALETTES = {
    viridis: RAMP([[68,1,84],[72,26,108],[71,47,125],[65,68,135],[57,86,140],[49,104,142],[42,120,142],
      [35,136,142],[31,152,139],[34,168,132],[53,183,121],[84,197,104],[122,209,81],[165,219,54],[210,226,27],[253,231,37]]),
  };
  const VIRIDIS = PALETTES.viridis;
  const VIRIDIS_CSS = "linear-gradient(90deg, #440154, #3b528b, #21918c, #5ec962, #fde725)";
  // A gene in ink mode spans from its lightest printable tint - the 0.2 floor every
  // occupied square gets - to the full color. Show that span rather than one chip,
  // so the legend reads on the same scale as the map.
  const inkCss = (hex) => {
    const [r, g, b] = hexToRgb255(hex);
    const pale = [r, g, b].map((v) => Math.round((1 - 0.2 * (1 - v / 255)) * 255));
    return `linear-gradient(90deg, rgb(${pale.join(",")}), ${hex})`;
  };
  const rampCss = (name) => {
    const lut = PALETTES[name] || PALETTES.viridis;
    const stops = [];
    for (let i = 0; i <= 6; i++) {
      const t = Math.round((i / 6) * 255) * 3;
      stops.push(`rgb(${lut[t]},${lut[t + 1]},${lut[t + 2]})`);
    }
    return `linear-gradient(90deg, ${stops.join(", ")})`;
  };

  // The only clustering columns we publish. First match wins; a dataset carries one.
  const RNA_COLOR_BY = cfg.rnaColorBy || ["leiden_heli", "leiden_bigGroup"];

  // Put a name on every cluster, outside its own cloud of points and out of anyone
  // else's. Placement is solved once, in UMAP coordinates rather than screen pixels,
  // so a label keeps its bearing relative to its cluster at every zoom instead of
  // swimming around as you scroll.
  //
  // For each cluster: the centroid, then a ring of eight candidate positions just
  // beyond the cluster's own median radius. Each candidate is scored on how many
  // cells - of any cluster - sit under it, plus a penalty for crowding a label that
  // is already placed. Big clusters are placed first, so the crowded small ones
  // inherit the leftover room rather than the other way round.
  function placeClusterLabels(R, o, codes, lut) {
    const K = o.categories.length, n = R.n, pos = R.pos;
    const b = S.m.bounds, spanX = Math.max(b[2] - b[0], 1e-6), spanY = Math.max(b[3] - b[1], 1e-6);
    const span = Math.max(spanX, spanY);

    // occupancy grid over the whole embedding: one O(1) density lookup per candidate
    const GN = 160, gx = (x) => Math.min(GN - 1, Math.max(0, ((x - b[0]) / spanX * GN) | 0));
    const gy = (y) => Math.min(GN - 1, Math.max(0, ((y - b[1]) / spanY * GN) | 0));
    const grid = new Float32Array(GN * GN);
    for (let i = 0; i < n; i++) grid[gy(pos[i * 2 + 1]) * GN + gx(pos[i * 2])] += 1;

    let mx0 = 0, my0 = 0;
    for (let i = 0; i < n; i++) { mx0 += pos[i * 2]; my0 += pos[i * 2 + 1]; }
    mx0 /= Math.max(n, 1); my0 /= Math.max(n, 1);   // centre of mass of the whole plot

    const sum = new Float64Array(K * 2), cnt = new Float64Array(K);
    for (let i = 0; i < n; i++) {
      const k = codes[i]; if (k >= K) continue;
      sum[k * 2] += pos[i * 2]; sum[k * 2 + 1] += pos[i * 2 + 1]; cnt[k] += 1;
    }
    // median distance to the centroid: robust to the stragglers every cluster has
    const rad = new Float64Array(K), bucket = [];
    for (let k = 0; k < K; k++) bucket.push([]);
    const stride = Math.max(1, Math.floor(n / 60000));        // sampled: exact radius is not needed
    for (let i = 0; i < n; i += stride) {
      const k = codes[i]; if (k >= K || !cnt[k]) continue;
      const dx = pos[i * 2] - sum[k * 2] / cnt[k], dy = pos[i * 2 + 1] - sum[k * 2 + 1] / cnt[k];
      bucket[k].push(Math.hypot(dx, dy));
    }
    for (let k = 0; k < K; k++) {
      const a = bucket[k]; if (!a.length) { rad[k] = span * 0.02; continue; }
      a.sort((x, y) => x - y); rad[k] = Math.max(a[(a.length * 0.6) | 0], span * 0.01);
    }

    // east first, so a label lands to the right of its cluster whenever that is free
    const DIRS = [[1, 0], [0.71, 0.71], [0.71, -0.71], [0, 1], [0, -1],
                  [-0.71, 0.71], [-0.71, -0.71], [-1, 0]];
    const probe = span * 0.035;                               // radius the label text occupies
    // The text is drawn at a fixed pixel size, so its footprint in UMAP units depends
    // on the zoom. Guessing it as a fraction of `span` was out by half on a narrow
    // panel, which is why two labels could overlap while the solver thought they were
    // clear. Derive it from the same fit the renderer uses instead.
    const L0 = rnaLayout();
    const ppu = Math.max(1e-9, Math.min((L0.W / L0.cols) / spanX, (L0.H / L0.rows) / spanY) * 0.88);
    const labelBox = (name, x, y, anchor) => {
      const hw = (3.7 * name.length + 7) / ppu, hh = 9 / ppu;
      return { x: anchor === "start" ? x + hw : anchor === "end" ? x - hw : x, y, hw, hh };
    };
    const density = (x, y) => {
      let t = 0;
      const r = Math.max(1, Math.round(probe / spanX * GN)), cx = gx(x), cy = gy(y);
      for (let j = -r; j <= r; j++) for (let i2 = -r; i2 <= r; i2++) {
        const a = cx + i2, c = cy + j;
        if (a >= 0 && a < GN && c >= 0 && c < GN) t += grid[c * GN + a];
      }
      return t / ((2 * r + 1) * (2 * r + 1));
    };

    // Positions placed by hand in content/label_editor.html win outright. No solver
    // gets this as right as someone looking at the plot, and once a figure is settled
    // it should not move because the heuristic was retuned.
    const fixed = (cfg.labelPositions && cfg.labelPositions[S.dsId]) || {};
    const order = Array.from({ length: K }, (_, k) => k).filter((k) => cnt[k] > 0)
      .sort((a, c) => cnt[c] - cnt[a]);
    const placed = [], out = [];
    const STEPS = [0.15, 0.3, 0.45, 0.65, 0.9, 1.2, 1.6];
    // the fit leaves a small margin round the data, so a label pushed past the
    // bounds is a label nobody can read - keep every one inside the frame
    const mx = spanX * 0.04, my = spanY * 0.04;
    const clampX = (v) => Math.min(b[2] - mx, Math.max(b[0] + mx, v));
    const clampY = (v) => Math.min(b[3] - my, Math.max(b[1] + my, v));
    for (const k of order) {
      const name = o.categories[k];
      const hand = fixed[name];
      if (hand && hand.length === 2 && Number.isFinite(hand[0]) && Number.isFinite(hand[1])) {
        // the editor stores the centre of the label, so it anchors dead centre
        const box = labelBox(name, hand[0], hand[1], "middle");
        placed.push(box);
        out.push({
          name, x: hand[0], y: hand[1], count: cnt[k],
          color: [...labelColor(lut[k] || [200, 200, 200]), 255],
          anchor: "middle", baseline: "center",
        });
        continue;
      }
      const cx = sum[k * 2] / cnt[k], cy = sum[k * 2 + 1] / cnt[k];
      // Scale both scoring terms to 0..1 before weighing them against each other.
      // density() counts cells per grid square, which on 175k cells runs into the
      // hundreds, while the distance term was a fraction of the plot width - about
      // three. So distance counted for nothing and every label went wherever the
      // plot was emptiest, which is a long way from the cluster it names.
      const core = Math.max(density(cx, cy), 1);
      const reach = Math.max(rad[k] * 2, 1e-6);
      // Which way is "out" for this cluster: away from the middle of the plot. A
      // label on an outlying cluster belongs on its far side, in the empty space,
      // rather than tucked back toward the crowd - that is where there is room, and
      // it is the side a reader's eye is already on. Fixed east said nothing about
      // the plot and put Perineurial's label between it and everything else.
      const ox0 = cx - mx0, oy0 = cy - my0, om = Math.hypot(ox0, oy0) || 1;
      const outX = ox0 / om, outY = oy0 / om;
      let best = null;
      // Every direction at every distance is a candidate, not just the first gap
      // along each ray. With one candidate per direction a label boxed in by its
      // neighbours had nowhere to go and simply sat on top of one of them; now it
      // can step further out, and the distance term keeps it as close as it can be.
      DIRS.forEach(([ux, uy], di) => {
        for (const t of STEPS) {
          const d = rad[k] * t + probe * 0.5;
          const x = clampX(cx + ux * d), y = clampY(cy + uy * d);
          const anchor = ux > 0.3 ? "start" : ux < -0.3 ? "end" : "middle";
          const box = labelBox(o.categories[k], x, y, anchor);
          // clamping the anchor is not enough: the text grows away from it, so a
          // label anchored just inside the frame can still run off the edge
          const over = Math.max(0, (box.x + box.hw) - (b[2] - mx)) + Math.max(0, (b[0] + mx) - (box.x - box.hw))
                     + Math.max(0, (box.y + box.hh) - (b[3] - my)) + Math.max(0, (b[1] + my) - (box.y - box.hh));
          // Sit on the edge of the cloud: far enough out that the label is not lost
          // in the thick of the points, no further. Both terms are now fractions of
          // the same size, so the trade between them is the one intended.
          // Emptiness decides, the rest only breaks ties. These were close enough in
          // magnitude that "do not point back toward the middle" could outvote "do
          // not sit on top of the points", which is how a label ended up over its own
          // cluster with clear space just the other side of it.
          const inward = 1 - (ux * outX + uy * outY);   // 0 straight out, 2 straight in
          let sc = (density(x, y) / core) * 6 + (d / reach) * 3 + inward * 0.3
                 + (over / span) * 150;
          for (const q of placed) {
            const ox = (box.hw + q.hw) - Math.abs(box.x - q.x);
            const oy = (box.hh + q.hh) - Math.abs(box.y - q.y);
            if (ox > 0 && oy > 0) sc += Math.min(ox / (box.hw + q.hw), oy / (box.hh + q.hh)) * 200;
          }
          if (!best || sc < best.sc) best = { sc, x, y, ux, uy, box };
        }
      });
      placed.push(best.box);
      out.push({
        name, x: best.x, y: best.y, count: cnt[k],
        color: [...labelColor(lut[k] || [200, 200, 200]), 255],
        anchor: best.ux > 0.3 ? "start" : best.ux < -0.3 ? "end" : "middle",
        // world +y is up here (flipY is false), so a label above its cluster sits on
        // its own bottom edge
        baseline: best.uy > 0.3 ? "bottom" : best.uy < -0.3 ? "top" : "center",
      });
    }
    return out;
  }

  function leaveRna() {
    if (S.mode !== "rna") return;
    S.mode = "image"; S.rna = null; S.viewState = null;
    deckgl.setProps({ views: imageViews(), viewState: { target: [0, 0, 0], zoom: 0 }, layers: [] });
    $("umapLabels").hidden = true; $("umapLabels").innerHTML = "";
    $("densityKey").hidden = true;
    $("panelRna").hidden = true; $("panelChannels").hidden = false;
    $("hud").hidden = false; $("scalebar").hidden = false;
  }

  async function openRna() {
    S.mode = "rna";
    const m = S.m;
    for (const id of ["panelChannels", "panelGenes", "panelOverlays", "btnAllOn", "btnAllOff", "btnImage", "hud", "scalebar", "densityKey"]) $(id).hidden = true;
    $("panelRna").hidden = false; $("umapLabels").hidden = false;
    $("dsInfo").textContent = `${m.nCells.toLocaleString()} cells · ${m.nGenes.toLocaleString()} genes · ${m.valueLabel}`;
    status("Loading UMAP…");
    let pos, index;
    try {
      [pos, index] = await Promise.all([S.src.buf("umap.bin").then((b) => new Float32Array(b)),
        S.src.buf("genes.json").then((b) => JSON.parse(new TextDecoder().decode(b).replace(/:(NaN|-?Infinity)\b/g, ":null")))]);
    } catch (e) { status(`Cannot load UMAP: ${e.message}`, true); return; }
    const n = pos.length / 2;
    S.rna = { pos, n, index, genes: [], obs: m.obs || [], colorBy: null, obsData: null, grid: "", vs: null,
              labels: null, labelsKey: null, cat: null, showLabels: $("rnaLabels").checked,
              view: "normal", palette: "viridis", bin: 1, heat: null,
              pointSize: n > 200000 ? 1 : n > 50000 ? 1.5 : 2.5,
              baseData: { length: n, attributes: { getPosition: { value: pos, size: 2 } } } };
    $("rnaPointSize").value = S.rna.pointSize;
    $("rnaView").value = S.rna.view; $("rnaPalette").value = S.rna.palette;
    $("rnaBin").value = String(S.rna.bin); $("rowRnaBin").hidden = S.rna.view !== "density";
    // One clustering per dataset, so there is nothing to choose between and no picker.
    // RNA_COLOR_BY is tried in order; the fallback only matters if a future dataset
    // arrives without either column.
    const cats = S.rna.obs.filter((o) => o.type === "category" && !/sample|doublet/i.test(o.name));
    const pick = RNA_COLOR_BY.map((want) => cats.find((o) => o.name === want)).find(Boolean)
      || cats.find((o) => /cluster|leiden|louvain/i.test(o.name)) || cats[0] || null;
    await setColorBy(pick ? pick.name : "");
    buildRnaGeneUI();
    status(null);
  }

  async function setColorBy(name) {
    const R = S.rna; if (!R) return;
    R.colorBy = name || null; R.obsData = null; R.labels = null; R.labelsKey = null; R.cat = null;
    const o = R.obs.find((x) => x.name === name);
    const leg = $("rnaLegend"); leg.innerHTML = "";
    if (o) {
      const buf = await S.src.buf(`obs/${o.file}`);
      const colors = new Uint8Array(R.n * 4);
      if (o.type === "category") {
        const codes = new Uint16Array(buf);
        const lut = o.categories.map((_, i) => liftForDark(o.colors && o.colors[i] ? hexToRgb255(o.colors[i]) : catColor(i)));
        for (let i = 0; i < R.n; i++) { const c = lut[codes[i]] || [128, 128, 128]; const k = i * 4; colors[k] = c[0]; colors[k + 1] = c[1]; colors[k + 2] = c[2]; colors[k + 3] = 255; }
        leg.innerHTML = o.categories.map((c, i) => `<div class="leg"><i style="background:rgb(${lut[i].join(",")})"></i><span>${c}</span><small>${(o.counts[i] || 0).toLocaleString()}</small></div>`).join("");
        // Placement depends on how big a label is in UMAP units, which depends on
        // the fitted zoom - and at this point the panel has not been laid out yet
        // (clientWidth is still 0 on a first open). Keep the inputs and let
        // renderRna solve it once it knows the real geometry.
        R.cat = { o, codes, lut };
      } else {
        const vals = new Float32Array(buf), lo = o.min, span = Math.max(o.max - o.min, 1e-6);
        for (let i = 0; i < R.n; i++) {
          const t = Math.max(0, Math.min(255, Math.round(((vals[i] - lo) / span) * 255))), k = i * 4;
          colors[k] = VIRIDIS[t * 3]; colors[k + 1] = VIRIDIS[t * 3 + 1]; colors[k + 2] = VIRIDIS[t * 3 + 2]; colors[k + 3] = Number.isFinite(vals[i]) ? 255 : 0;
        }
        leg.innerHTML = `<div class="cbar"><span>${o.min.toFixed(2)}</span><i style="background:${VIRIDIS_CSS}"></i><span>${o.max.toFixed(2)}</span></div>`;
      }
      R.obsData = { length: R.n, attributes: { getPosition: { value: R.pos, size: 2 }, getFillColor: { value: colors, size: 4, normalized: true } } };
    }
    renderRna();
  }

  function showRnaResults(q) {
    const box = $("rnaResults"), R = S.rna;
    q = q.trim().toLowerCase();
    if (!q || !R) { box.hidden = true; return; }
    const hits = [];
    const exact = R.index.find((g) => g.n.toLowerCase() === q);          // an exact name always comes first
    if (exact) hits.push(exact);
    for (const g of R.index) if (g !== exact && g.n.toLowerCase().startsWith(q)) { hits.push(g); if (hits.length >= 40) break; }
    if (hits.length < 40) for (const g of R.index) { const n = g.n.toLowerCase(); if (!n.startsWith(q) && n.includes(q)) { hits.push(g); if (hits.length >= 40) break; } }
    box.innerHTML = "";
    hits.forEach((g, i) => {
      const d = document.createElement("div");
      d.innerHTML = `<span>${g.n}</span><small>${(g.p * 100).toFixed(1)}% of cells</small>`;
      if (i === 0) d.classList.add("hot");
      d.onmousedown = (e) => { e.preventDefault(); addRnaGene(g.n); };
      box.appendChild(d);
    });
    box.hidden = hits.length === 0;
  }

  async function addRnaGene(name) {
    const R = S.rna; if (!R) return;
    $("rnaResults").hidden = true; $("rnaSearch").value = "";
    if (R.genes.some((g) => g.name === name)) return;
    if (R.genes.length >= MAX_RNA_GENES) {
      status(`Up to ${MAX_RNA_GENES} genes at once. Remove one to add ${name}.`, true);
      setTimeout(() => status(null), 4000);
      return;
    }
    status(`Loading ${name}…`);
    let buf;
    try { buf = await S.src.buf(`gene/${encodeURIComponent(name)}.bin`); }
    catch (e) { status(`${name}: not expressed in this dataset`, true); setTimeout(() => status(null), 4000); return; }
    const k = new Uint32Array(buf, 0, 1)[0];
    const idx = new Uint32Array(buf, 4, k), val = new Float32Array(buf, 4 + 4 * k, k);
    const order = Array.from({ length: k }, (_, i) => i).sort((a, b) => val[a] - val[b]);   // brightest drawn last
    const pos = new Float32Array(k * 2), v = new Float32Array(k);
    order.forEach((o, j) => { pos[2 * j] = R.pos[2 * idx[o]]; pos[2 * j + 1] = R.pos[2 * idx[o] + 1]; v[j] = val[o]; });
    // default display range = 5th..99th percentile of the expressing cells (quantile cutoffs, as Seurat FeaturePlot)
    const q = (f) => (k ? v[Math.min(k - 1, Math.floor(f * k))] : 0);
    const vmax = k ? v[k - 1] : 1, p99 = q(0.99) || vmax;
    let p05 = q(0.05);
    if (p05 >= p99) p05 = 0;
    const g = { name, k, pos, val: v, min: +p05.toFixed(2), max: +p99.toFixed(2), vmax: Math.max(vmax, 0.01),
                frac: k / R.n, hex: GENE_COLORS[R.genes.length % GENE_COLORS.length], visible: true };
    recolorRnaGene(g);
    if (!S.rna || S.rna !== R) return;                 // dataset switched while loading
    R.genes.push(g); R.heat = null; R.grid = "";
    buildRnaGeneUI();
    renderRna();
    status(null);
  }

  function recolorRnaGene(g) {
    const R = S.rna;
    const colors = new Uint8Array(g.k * 4), lo = g.min, span = Math.max(g.max - lo, 1e-6);
    const useInk = R && R.palette === "genes";
    const lut = (R && PALETTES[R.palette]) || VIRIDIS;
    const [hr, hg, hb] = hexToRgb255(g.hex || "#33ddff");
    for (let i = 0; i < g.k; i++) {
      const x = g.val[i], k = i * 4;
      if (x < lo || x <= 0) { colors[k + 3] = 0; continue; }
      const f = Math.min(1, (x - lo) / span);
      if (useInk) {
        // the gene's own color, faint where it is barely expressed
        colors[k] = hr; colors[k + 1] = hg; colors[k + 2] = hb; colors[k + 3] = 40 + Math.round(215 * f);
      } else {
        const t = Math.min(255, Math.round(f * 255)) * 3;
        colors[k] = lut[t]; colors[k + 1] = lut[t + 1]; colors[k + 2] = lut[t + 2]; colors[k + 3] = 255;
      }
    }
    g.data = { length: g.k, attributes: { getPosition: { value: g.pos, size: 2 }, getFillColor: { value: colors, size: 4, normalized: true } } };
  }

  function buildRnaGeneUI() {
    const R = S.rna; if (!R) return;
    $("rnaCount").textContent = `(${R.genes.length} / ${MAX_RNA_GENES})`;
    const list = $("rnaGeneList"); list.innerHTML = "";
    R.genes.forEach((g) => {
      const row = document.createElement("div");
      row.className = "chan";
      const sw = R.palette === "genes" ? inkCss(g.hex) : rampCss(R.palette);
      row.innerHTML = `<input class="cb" type="checkbox" ${g.visible === false ? "" : "checked"} /><span class="col sw" style="background:${sw}"></span>
        <span class="name" title="${g.name}">${g.name}</span><button class="rm x" title="remove">×</button><div class="range"></div>`;
      row.querySelector(".x").onclick = () => { R.genes = R.genes.filter((x) => x !== g); buildRnaGeneUI(); renderRna(); };
      const step = g.vmax > 20 ? 0.5 : 0.05;
      rangeSlider(row.querySelector(".range"), {
        min: 0, max: +g.vmax.toFixed(2), step, lo: g.min, hi: Math.min(g.max, +g.vmax.toFixed(2)),
        color: () => "#5ec962", fmt: (x) => x.toFixed(2),
        onChange: (a, b) => {
          g.min = a; g.max = b; recolorRnaGene(g);
          const hk = `${g.hex}|${a}`;
          if (g.heatKey !== hk) { g.heatKey = hk; R.heat = null; }
          renderRna();
        },
      });
      const cb = row.querySelector(".cb");
      cb.title = `${g.k.toLocaleString()} cells expressing`;
      cb.onchange = () => {
        g.visible = cb.checked;
        row.classList.toggle("off", !cb.checked);
        // a hidden gene loses its panel in Normal view, so the grid is re-fitted
        R.heat = null; R.grid = "";
        renderRna();
      };
      row.classList.toggle("off", g.visible === false);
      list.appendChild(row);
    });
  }

  // The UMAP's density view: the same square-bin model as the slide, in embedding
  // space. Density collapses the small multiples into ONE panel on purpose - the
  // point of it is that genes combine, either as inks that multiply or as one ramp
  // over their total, and neither says anything if each gene keeps its own panel.
  function rnaHeat() {
    const R = S.rna;
    const b = S.m.bounds, span = Math.max(b[2] - b[0], b[3] - b[1], 1e-6);
    const cell = (span / 200) * Math.max(1, R.bin);
    const ink = R.palette === "genes";
    const on = R.genes.filter((g) => g.visible !== false);
    const key = on.map((g) => `${g.name}:${g.hex}:${g.min}`).join("|") + `#${R.palette}#${R.bin}`;
    if (R.heat && R.heat.key === key) return R.heat;

    const idx = new Map(), KX = 1 << 21;
    const binOf = (x, y) => (Math.floor((y - b[1]) / cell)) * KX + Math.floor((x - b[0]) / cell);
    for (const g of on) for (let i = 0; i < g.k; i++) {
      const k = binOf(g.pos[2 * i], g.pos[2 * i + 1]);
      if (!idx.has(k)) idx.set(k, idx.size);
    }
    const n = idx.size;
    const positions = new Float32Array(n * 3), colors = new Uint8Array(n * 4);
    for (const [k, j] of idx) {
      const iy = Math.floor(k / KX), ix = k - iy * KX;
      positions[j * 3] = b[0] + (ix + 0.5) * cell; positions[j * 3 + 1] = b[1] + (iy + 0.5) * cell;
    }
    // expression summed over the cells in each square, scaled by how much the busiest
    // squares hold - the same correction the slide view needed
    const norm = on.map((g) => {
      const sum = new Float32Array(n), lo = g.min;
      for (let i = 0; i < g.k; i++) {
        const x = g.val[i];
        if (x >= lo && x > 0) sum[idx.get(binOf(g.pos[2 * i], g.pos[2 * i + 1]))] += x;
      }
      const nz = Array.from(sum).filter((v) => v > 0).sort((x, y) => x - y);
      const hi = nz.length ? Math.max(nz[Math.floor((nz.length - 1) * 0.99)], 1e-6) : 1;
      for (let j = 0; j < n; j++) sum[j] = Math.min(1, sum[j] / hi);
      return sum;
    });

    if (ink) {
      const paper = new Float32Array(n * 3).fill(1), miss = new Float32Array(n).fill(1);
      on.forEach((g, gi) => {
        const d = norm[gi], [tr, tg, tb] = hexToRgb255(g.hex);
        for (let j = 0; j < n; j++) {
          if (d[j] <= 0) continue;
          const sEff = 0.2 + 0.8 * d[j];
          paper[j * 3] *= 1 - sEff * (1 - tr / 255);
          paper[j * 3 + 1] *= 1 - sEff * (1 - tg / 255);
          paper[j * 3 + 2] *= 1 - sEff * (1 - tb / 255);
          miss[j] *= 1 - sEff;
        }
      });
      for (let j = 0; j < n; j++) {
        const c = j * 4;
        colors[c] = Math.round(paper[j * 3] * 255); colors[c + 1] = Math.round(paper[j * 3 + 1] * 255);
        colors[c + 2] = Math.round(paper[j * 3 + 2] * 255);
        colors[c + 3] = Math.round(Math.sqrt(Math.max(0, 1 - miss[j])) * 255);
      }
    } else {
      const lut = PALETTES[R.palette] || VIRIDIS;
      const inv = 1 / Math.max(1, norm.length);      // average, not sum - see geneHeat
      for (let j = 0; j < n; j++) {
        let total = 0, any = false;
        for (const d of norm) { total += d[j]; if (d[j] > 0) any = true; }
        const t = Math.min(255, Math.round((0.06 + 0.94 * Math.min(1, total * inv)) * 255)) * 3;
        const c = j * 4;
        colors[c] = lut[t]; colors[c + 1] = lut[t + 1]; colors[c + 2] = lut[t + 2];
        colors[c + 3] = any ? 235 : 0;
      }
    }
    R.heat = { key, n, positions, colors, cell, index: Array.from({ length: n }, (_, j) => j) };
    return R.heat;
  }


  // A gene switched off keeps its row and its window, but loses its panel and
  // stops contributing to the combined density.
  const rnaShownGenes = () => (S.rna ? S.rna.genes.filter((g) => g.visible !== false) : []);

  function rnaLayout() {
    const R = S.rna, el = $("viewport");
    const shown = rnaShownGenes();
    const dens = R.view === "density" && shown.length > 0;
    const n = dens ? 2 : 1 + shown.length;
    const W = el.clientWidth, H = el.clientHeight;
    let best = { cols: 1, rows: n, size: 0 };
    for (let c = 1; c <= n; c++) { const r = Math.ceil(n / c), size = Math.min(W / c, H / r); if (size > best.size) best = { cols: c, rows: r, size }; }
    return { n, W, H, ...best };
  }

  function fitRna() {
    const R = S.rna; if (!R) return;
    R.vs = null;                                              // renderRna refits when there is no view state
    renderRna();
  }

  function rnaViewStates() {
    const R = S.rna, out = {};
    for (let i = 0; i <= R.genes.length; i++) out[`umap-${i}`] = R.vs;
    return out;
  }

  function renderRna() {
    const R = S.rna; if (!R || !deckgl) return;
    const L = rnaLayout();
    if (!L.W || !L.H) { renderPending = true; return; }
    const gridKey = `${L.cols}x${L.rows}@${L.W}x${L.H}`;   // L.n already follows what is shown
    if (!R.vs || R.grid !== gridKey) {                        // first draw, or panel size changed: refit
      const b = S.m.bounds, pw = L.W / L.cols, ph = L.H / L.rows;
      const zoom = Math.log2(Math.min(pw / Math.max(b[2] - b[0], 1e-6), ph / Math.max(b[3] - b[1], 1e-6)) * 0.88);
      R.vs = { target: [(b[0] + b[2]) / 2, (b[1] + b[3]) / 2, 0], zoom, minZoom: zoom - 3, maxZoom: zoom + 8 };
    }
    R.grid = gridKey;
    // Re-solved whenever the fit changes - a resize changes how much of the UMAP a
    // label covers, and so which placement is actually free.
    const key = `${R.colorBy}|${gridKey}`;
    if (R.cat && R.labelsKey !== key) {
      R.labels = placeClusterLabels(R, R.cat.o, R.cat.codes, R.cat.lut);
      R.labelsKey = key;
    }
    const ctl = { scrollZoom: { speed: 0.01, smooth: false }, inertia: false, doubleClickZoom: true, keyboard: false };
    const views = [], labels = [];
    const shown = rnaShownGenes();
    const dens = R.view === "density" && shown.length > 0;
    for (let i = 0; i < L.n; i++) {
      const c = i % L.cols, r = Math.floor(i / L.cols);
      const box = { x: `${(c * 100) / L.cols}%`, y: `${(r * 100) / L.rows}%`, width: `${100 / L.cols}%`, height: `${100 / L.rows}%` };
      views.push(new OrthographicView({ id: `umap-${i}`, ...box, flipY: false, controller: ctl }));
      const style = `left:${box.x};top:${box.y};width:${box.width};height:${box.height}`;
      if (i === 0) {
        labels.push(`<div class="umap-cell" style="${style}"><div class="t">${R.colorBy ? `Color: ${R.colorBy}` : "All cells"}</div></div>`);
      } else if (dens) {
        const names = shown.map((g) => g.name).join(" + ");
        // say what the colors mean, same as the slide view's key
        const key = R.palette === "genes"
          ? `<div class="genes">${shown.map((g) => `<span class="g"><i style="background:${inkCss(g.hex)}"></i>${g.name}</span>`).join("")}</div>` +
            `<div class="note">pale = few cells expressing, saturated = many` +
            (shown.length > 1 ? ". Where genes overlap the colors mix and darken." : "") + `</div>`
          : `<div class="ramp"><span>low</span><i style="background:${rampCss(R.palette)}"></i><span>high</span></div>` +
            `<div class="note">expression per square` +
            (shown.length > 1 ? `, averaged over ${shown.length} genes` : ``) + `</div>`;
        labels.push(`<div class="umap-cell" style="${style}"><div class="t">${names} <small>combined density</small></div>` +
          `<div class="dkey inpanel">${key}</div></div>`);
      } else {
        const g = shown[i - 1];
        const sw = R.palette === "genes" ? inkCss(g.hex) : rampCss(R.palette);
        labels.push(`<div class="umap-cell" style="${style}"><div class="t">${g.name} <small>${(g.frac * 100).toFixed(1)}% expressing</small></div>` +
          `<div class="cbar"><span>${g.min.toFixed(2)}</span><i style="background:${sw}"></i><span>${g.max.toFixed(2)}</span></div></div>`);
      }
    }
    $("umapLabels").innerHTML = labels.join("");
    const common = { radiusUnits: "pixels", getRadius: R.pointSize, radiusMinPixels: 0.5, stroked: false, filled: true, parameters: NORMAL, antialiasing: true };
    const layers = [new ScatterplotLayer({ id: "umap-base", data: R.baseData, getFillColor: [80, 86, 100, 255], ...common,
      updateTriggers: { getRadius: [R.pointSize] } })];
    if (R.obsData) layers.push(new ScatterplotLayer({ id: `umap-0-obs-${R.colorBy}`, data: R.obsData, ...common }));
    if (R.obsData && R.showLabels && R.labels && R.labels.length) {
      layers.push(new TextLayer({
        id: `umap-0-labels-${R.colorBy}`, data: R.labels,
        getPosition: (d) => [d.x, d.y], getText: (d) => d.name, getColor: (d) => d.color,
        getTextAnchor: (d) => d.anchor, getAlignmentBaseline: (d) => d.baseline,
        sizeUnits: "pixels", getSize: 14, sizeMinPixels: 12, sizeMaxPixels: 15,
        characterSet: "auto", fontFamily: "system-ui, sans-serif", fontWeight: 700,
        // An SDF outline was tried here and made the glyphs thinner without reading as
        // a halo at this size. A near-opaque plate with generous padding does the job:
        // it is the plate, not the letterform, that lifts the label off the cloud.
        background: true, getBackgroundColor: [8, 11, 18, 232],
        backgroundPadding: [7, 4, 7, 4],
        getBorderColor: (d) => [...d.color.slice(0, 3), 110], getBorderWidth: 1,
        parameters: NORMAL, pickable: false,
        updateTriggers: { getPosition: [R.colorBy], getText: [R.colorBy] },
      }));
    }
    if (dens) {
      const h = rnaHeat();
      if (h.n) layers.push(new ColumnLayer({
        id: `umap-1-heat-${R.palette}`, data: h.index,
        getPosition: (j) => [h.positions[j * 3], h.positions[j * 3 + 1], 0],
        getFillColor: (j) => [h.colors[j * 4], h.colors[j * 4 + 1], h.colors[j * 4 + 2], h.colors[j * 4 + 3]],
        diskResolution: 4, angle: 45, radius: (h.cell * 1.02) / Math.SQRT2,
        radiusUnits: "common", extruded: false, filled: true, stroked: false,
        parameters: NORMAL, pickable: false,
        updateTriggers: { getFillColor: [h.key], getPosition: [h.key] },
      }));
    } else {
      shown.forEach((g, i) => layers.push(new ScatterplotLayer({ id: `umap-${i + 1}-g-${g.name}`, data: g.data, ...common })));
    }
    deckgl.setProps({ views, viewState: rnaViewStates(), layers });
  }

  async function openDataset(id) {
    status("Opening dataset…");
    S.dsId = id;
    tileCache.clear();
    S.channels = []; S.overlays = []; S.m = null;
    deckgl.setProps({ layers: [] });           // drop the previous dataset's layers right away
    S.genes = []; S.tissue.positions = null; S.heat = null;
    leaveRna();
    try {
      const d = S.datasets.find((x) => x.id === id);
      const staticRoot = params.get("static");   // ?static=http://host/folder -> <folder>/<id>/ (test exports before upload)
      if (staticRoot) S.src = staticSource(`${staticRoot.replace(/\/$/, "")}/${id}/`);
      else if (TOKEN_URL && !(d && (d.linkPath || d.base))) {
        const t = await ensureToken();           // asks for the password once per session
        S.src = staticSource(`${t.base}${id}/`, tokenQuery);
      }
      else if (d && d.linkPath) S.src = staticSource(await resolveLink(d.linkPath));
      else if (d && d.base) S.src = staticSource(d.base);
      else S.src = apiSource(id);
      S.m = await S.src.manifest();
    } catch (e) { status(`Cannot open ${id}: ${e.message}`, true); return; }
    if (S.m.kind === "rnaseq") {
      await openRna();
      const u = new URL(location.href); u.searchParams.set("ds", id); history.replaceState(null, "", u);
      return;
    }
    const img = S.m.image;
    const wanted = (cfg.defaultChannels || []).map((s) => s.toLowerCase());
    // "open-atto550", "open AF488", "open-cy5": detector slots with no antibody in
    // them. They carry autofluorescence and nothing a viewer should reason about,
    // and on the Akoya panels they are the three biggest files in the set.
    // Panel order out of the scanner is acquisition order, which is nothing to the
    // reader: finding CD8 among thirty markers meant reading the whole list. Sort
    // them, with the nuclear counterstain pinned at the top because it is the one
    // channel people leave on as a reference for all the others. `index` keeps
    // pointing at the right tiles, so only the row order moves.
    // Left alone on an RGB photo, where the channels are not markers at all.
    // Per-dataset exclusions, by exact name, on top of the blanket `open-*` rule.
    // Some panels carry an antibody that was only ever there to test the stain; it is
    // real data but not part of the panel, and leaving it in the list invites someone
    // to read it as one more marker.
    const dsCfg = (cfg.datasets && cfg.datasets[S.dsId]) || {};
    const dropped = new Set((dsCfg.hideChannels || []).map((x) => String(x).trim().toLowerCase()));
    const named = img.channels.filter((c) =>
      !HIDE_CHANNEL.test(c.name) && !dropped.has(c.name.trim().toLowerCase()));
    const photo = named.length === 3 && named.every((c) => /^(red|green|blue)$/i.test(c.name));
    if (!photo && named.length > 4) {
      const nuclear = (c) => (NUCLEAR_CHANNEL.test(c.name) ? 0 : 1);
      named.sort((a, b) => nuclear(a) - nuclear(b) ||
        a.name.localeCompare(b.name, undefined, { numeric: true, sensitivity: "base" }));
    }
    S.channels = named.map((c, i) => {
      const on = wanted.length ? wanted.some((w) => c.name.toLowerCase().includes(w)) : i < 3;
      return { index: c.index, name: c.name, hex: c.color, color: hexToRgb01(c.color),
               lo: Math.round(c.window[0]), hi: Math.round(c.window[1]), visible: on };
    });
    if (!S.channels.some((c) => c.visible)) S.channels.slice(0, 3).forEach((c) => (c.visible = true));
    buildChannelUI();
    const um = S.m.pixelSizeUm;
    $("dsInfo").textContent = `${img.width.toLocaleString()} × ${img.height.toLocaleString()} px · ${img.levels} levels · ${img.dtype}` +
      (um ? ` · ${um.toFixed(4)} µm/px` : "");

    const ex = S.m.expr || [];
    // Panels follow the platform: sequencing datasets (Visium) get transcripts and just an
    // image opacity control for the H&E; imaging datasets (Akoya, Xenium stains, TIFF) get
    // the per-channel marker list and no transcript panel.
    const isVisium = S.m.kind === "visium" || ex.length > 0;
    $("panelGenes").hidden = !isVisium;
    $("chanList").hidden = isVisium;
    $("chanCount").hidden = isVisium;
    $("btnAllOn").hidden = isVisium;
    $("btnAllOff").hidden = isVisium;
    // Only where there is something underneath to see. On a Visium slide the image
    // is an H&E the transcripts are drawn over, so hiding it is useful; on a marker
    // panel the image IS the data, and the button only ever blanked the view.
    $("btnImage").hidden = !isVisium;
    if (!isVisium) S.imgOn = true;
    $("chanTitle").textContent = isVisium ? "Image" : "Markers / stains";
    if (isVisium) S.channels.forEach((c) => (c.visible = true));
    // Only the resolutions worth offering. The 2 um grid stays in the exports but is
    // not listed: at that size a bin holds a single transcript nearly everywhere, so
    // a density map of it is a presence map. The counts that give the view meaning
    // only appear once a bin is 8 um or larger.
    const offered = ex.filter((e) => !e.binSizeUm || e.binSizeUm >= EXPR_MIN_UM);
    const shownEx = offered.length ? offered : ex;
    const sel = $("exprSelect"); sel.innerHTML = "";
    for (const e of shownEx) {
      const o = document.createElement("option"); o.value = e.key;
      o.textContent = `${e.label} · ${e.nBins.toLocaleString()} bins · ${e.nGenes.toLocaleString()} genes`;
      sel.appendChild(o);
    }
    const isPhoto = img.channels.length === 3 && img.channels.every((c) => /^(red|green|blue)$/i.test(c.name));
    S.geneBlend = "heatmap";        // density: glow conflates intensity with overlap
    S.genePalette = "genes"; S.geneBin = 1; S.heat = null;
    $("geneBlend").value = S.geneBlend;
    $("genePalette").value = S.genePalette;
    S.imgOn = true; S.imgOpacity = 1; S.geneOpacity = 1;
    $("imgOpacity").value = 1; $("geneOpacity").value = 1;
    syncImageButton();
    S.overlays = (S.m.overlays || []).map((o) => ({ ...o, on: false, data: null }));
    buildOverlayUI();
    fitView();
    render();
    status(null);
    if (shownEx.length) { sel.value = shownEx[0].key; await setExprKey(shownEx[0].key); }
    const u = new URL(location.href); u.searchParams.set("ds", id); history.replaceState(null, "", u);
  }

  async function init() {
    initDeck();
    const cfgSets = cfg.datasets || {};
    S.datasets = Object.keys(cfgSets).map((id) => ({ id, label: cfgSets[id].label || id, linkPath: cfgSets[id].linkPath,
      group: cfgSets[id].group, base: cfgSets[id].base }));
    // ?static=<root>&ds=<id> can open an exported dataset that is not in config.js yet
    if (params.get("static") && params.get("ds") && !S.datasets.some((d) => d.id === params.get("ds")))
      S.datasets.push({ id: params.get("ds"), label: `${params.get("ds")} (local export)` });
    const useServer = !!(params.get("api") || cfg.apiBase || cfg.useLocalServer);
    $("apiInfo").textContent = useServer ? `server: ${API || location.origin}` : "data: static links (blob)";
    if (useServer) {
      try {
        const fromServer = await getJSON(apiUrl("datasets"));
        for (const d of fromServer) if (!S.datasets.some((x) => x.id === d.id)) S.datasets.push(d);
      } catch (e) {
        if (!S.datasets.length) {
          status(/^401/.test(e.message) ? "Password not accepted." :
            `No tile server at ${API || location.origin} — run run.ps1 or add ?api=http://host:8760`, true);
          return;
        }
      }
    }
    // Grouped by what was measured. The channel count used to be appended here, but
    // it counted raw channels including the `open-*` slots the viewer hides, so it
    // disagreed with both the label and the stain list; the labels carry it now.
    const sel = $("dsSelect"); sel.innerHTML = "";
    const groups = cfg.groups || {};
    const order = [...Object.keys(groups), undefined];
    const bucket = new Map(order.map((g) => [g, []]));
    for (const d of S.datasets) bucket.get(bucket.has(d.group) ? d.group : undefined).push(d);
    for (const g of order) {
      const items = bucket.get(g);
      if (!items || !items.length) continue;
      // an ungrouped dataset (a ?ds= local export, or a server listing) goes in flat
      const into = g === undefined && order.length === 1 ? sel
        : g === undefined ? Object.assign(document.createElement("optgroup"), { label: "Other" })
        : Object.assign(document.createElement("optgroup"), { label: groups[g] });
      for (const d of items) {
        const o = document.createElement("option"); o.value = d.id; o.textContent = d.label;
        into.appendChild(o);
      }
      if (into !== sel) sel.appendChild(into);
    }
    sel.onchange = () => openDataset(sel.value);
    $("btnFit").onclick = fitView;
    // One button for both directions: if anything is on it turns everything off,
    // otherwise it turns everything on. Rebuilding the whole channel list would refetch
    // every histogram, so just flip the checkboxes that are already there.
    // Two buttons rather than one toggle: "turn everything on" was otherwise two clicks
    // (off, then on) and you could not tell which way the toggle would go.
    const setAllChannels = (on) => {
      S.channels.forEach((c) => (c.visible = on));
      const rows = $("chanList").querySelectorAll(".chan");
      S.channels.forEach((c, i) => {
        const row = rows[i]; if (!row) return;
        row.querySelector(".cb").checked = on;
        row.classList.toggle("off", !on);
      });
      render();
    };
    $("btnAllOn").onclick = () => setAllChannels(true);
    $("btnAllOff").onclick = () => setAllChannels(false);
    $("imgOpacity").oninput = (e) => { S.imgOpacity = +e.target.value; requestRender(); };
    $("btnImage").onclick = () => { S.imgOn = !S.imgOn; syncImageButton(); render(); };
    $("geneOpacity").oninput = (e) => { S.geneOpacity = +e.target.value; requestRender(); };
    $("exprSelect").onchange = (e) => setExprKey(e.target.value);
    $("geneSearch").oninput = (e) => showGeneResults(e.target.value);
    $("geneSearch").onkeydown = (e) => {
      if (e.key === "Enter" || e.key === "Return" || e.keyCode === 13) {
        e.preventDefault();
        const hot = $("geneResults").querySelector(".hot span") || $("geneResults").querySelector("div span");
        if (hot) addGene(hot.textContent);
        else if (e.target.value.trim()) { status(`No gene named "${e.target.value.trim()}" in this layer`, true); setTimeout(() => status(null), 3500); }
      }
      if (e.key === "Escape") $("geneResults").hidden = true;
    };
    $("geneSearch").onblur = () => setTimeout(() => ($("geneResults").hidden = true), 150);
    const syncGeneMode = () => {
      const heat = S.geneBlend === "heatmap";
      $("rowGenePalette").hidden = !heat;
    };
    $("geneBlend").onchange = (e) => {
      S.geneBlend = e.target.value; S.heat = null;
      S.genes.forEach(recolorGene); syncGeneMode(); buildGeneUI(); render();
    };
    $("genePalette").onchange = (e) => { S.genePalette = e.target.value; S.heat = null; render(); };
    syncGeneMode();
    $("pointScale").oninput = (e) => { S.pointScale = +e.target.value; requestRender(); };
    $("tissueBins").onchange = async (e) => { S.tissue.on = e.target.checked; if (S.tissue.on && !S.tissue.positions) await loadTissue(); render(); };
    // The description is held with the data, not in this page. Anything written
    // into index.html is readable by anyone who opens the page - including the
    // parts of the project that are not public yet - so the body is fetched from
    // the private container, behind the same password and the same expiring
    // credential the imaging data uses.
    let aboutLoaded = false;
    async function loadAbout() {
      if (aboutLoaded) return;
      const body = $("aboutBody");
      body.innerHTML = '<p class="pending">Loading…</p>';
      try {
        const staticRoot = params.get("static");
        let url;
        if (staticRoot) url = `${staticRoot.replace(/\/$/, "")}/about.html`;
        else if (TOKEN_URL) { const t = await ensureToken(); url = `${t.base}about.html${tokenQuery()}`; }
        else throw new Error("no source configured");
        const r = await dataFetch(url, { cache: "no-cache" });
        if (r.status === 404) {
          body.innerHTML = '<p class="pending">No description has been published yet.</p>';
          return;                                   // not an error: it just is not written
        }
        if (!r.ok) throw new Error(`${r.status}`);
        body.innerHTML = await r.text();            // our own authored copy, from our own container
        aboutLoaded = true;
      } catch (e) {
        body.innerHTML = `<p class="pending">Could not load the description: ${e.message}</p>`;
      }
    }
    // Everything naming this particular project lives in config, so the viewer
    // itself carries no dataset, no institution and no byline. A different study
    // points it at different data and supplies its own three lines.
    const brand = cfg.brand || {};
    if (brand.title) { $("brandTitle").textContent = brand.title; document.title = brand.title; }
    $("brandLogos").innerHTML = (brand.logos || []).map((l) =>
      `<img class="${l.className || ""}" src="${l.src}" alt="${l.alt || ""}" />`).join("");
    $("colophon").textContent = brand.byline || "";

    const paper = $("lnkPaper"), paperUrl = (cfg.paperUrl || "").trim();
    if (paperUrl) {
      paper.href = paperUrl;
      paper.title = paperUrl;
    } else {
      paper.classList.add("pending");
      paper.removeAttribute("href");                 // not a link until there is one
      paper.title = "Added once the paper is published";
      paper.onclick = (e) => e.preventDefault();
    }
    const about = (open) => { $("aboutBack").hidden = !open; if (open) loadAbout(); };
    $("btnAbout").onclick = () => about(true);
    $("aboutClose").onclick = () => about(false);
    // clicking the backdrop closes; clicking inside the dialog must not
    $("aboutBack").onclick = (e) => { if (e.target === $("aboutBack")) about(false); };
    document.addEventListener("keydown", (e) => { if (e.key === "Escape") about(false); });
    $("chanSearch").oninput = () => filterChannelUI();
    const syncRnaMode = () => { $("rowRnaBin").hidden = !(S.rna && S.rna.view === "density"); };
    $("rnaView").onchange = (e) => {
      if (!S.rna) return;
      S.rna.view = e.target.value; S.rna.heat = null; S.rna.grid = "";
      syncRnaMode(); renderRna();
    };
    $("rnaPalette").onchange = (e) => {
      if (!S.rna) return;
      S.rna.palette = e.target.value; S.rna.heat = null;
      S.rna.genes.forEach(recolorRnaGene); buildRnaGeneUI(); renderRna();
    };
    $("rnaBin").onchange = (e) => { if (!S.rna) return; S.rna.bin = +e.target.value; S.rna.heat = null; renderRna(); };
    $("rnaLabels").onchange = (e) => { if (S.rna) { S.rna.showLabels = e.target.checked; renderRna(); } };
    $("rnaSearch").oninput = (e) => showRnaResults(e.target.value);
    $("rnaSearch").onkeydown = (e) => {
      if (e.key === "Enter") {
        e.preventDefault();
        const hot = $("rnaResults").querySelector(".hot span") || $("rnaResults").querySelector("div span");
        if (hot) addRnaGene(hot.textContent);
        else if (e.target.value.trim()) { status(`No gene named "${e.target.value.trim()}" in this dataset`, true); setTimeout(() => status(null), 3500); }
      }
      if (e.key === "Escape") $("rnaResults").hidden = true;
    };
    $("rnaSearch").onblur = () => setTimeout(() => ($("rnaResults").hidden = true), 150);
    $("rnaPointSize").oninput = (e) => { if (S.rna) { S.rna.pointSize = +e.target.value; renderRna(); } };
    $("rnaClear").onclick = () => {
      if (!S.rna) return;
      S.rna.genes = []; S.rna.heat = null; S.rna.grid = "";
      buildRnaGeneUI(); renderRna();
    };
    window.addEventListener("resize", () => updateHud());
    if (window.ResizeObserver) new ResizeObserver(() => {
      const el = $("viewport");
      if (!el.clientWidth || !el.clientHeight) return;
      if (fitPending) fitView();
      if (renderPending) render(); else updateHud();
    }).observe($("viewport"));

    if (!S.datasets.length) { status("No datasets: add entries to config.js (static links) or set apiBase", true); return; }
    const want = params.get("ds") || cfg.defaultDataset;
    const first = S.datasets.find((d) => d.id === want) || S.datasets[0];
    sel.value = first.id;
    await openDataset(first.id);
  }

  init();
})();
