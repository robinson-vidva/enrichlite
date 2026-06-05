// app.js - UI controller. Loads manifest, lazy-loads collections, drives worker.
(function () {
  "use strict";

  var state = {
    manifest: null,
    species: "human",
    collectionKey: null,
    cache: {},          // "species/collection" -> {symbols, aliases, terms, ...}
    lastResult: null,
    viewRows: [],       // full filtered+sorted set (paginated for display)
    sortKey: "p",
    sortDir: 1,
    page: 1,
    pageSize: 25,       // number, or "all"
    chartType: "dot",   // "dot" | "bar"
    topN: 25,
    currentSvg: null
  };

  var GENE_PREVIEW = 4;  // genes shown before "+N more" in the table

  // Coherent cell-cycle example set, per species (HGNC upper / MGI title-case).
  var EXAMPLE = {
    human: "CDK1 CDK2 CCNB1 CCNB2 CCNA2 CDC20 BUB1 BUB1B MAD2L1 AURKA AURKB PLK1 CCNE1 CDC25A CDC25B CDC25C ESPL1 NDC80 CENPA",
    mouse: "Cdk1 Cdk2 Ccnb1 Ccnb2 Ccna2 Cdc20 Bub1 Bub1b Mad2l1 Aurka Aurkb Plk1 Ccne1 Cdc25a Cdc25b Cdc25c Espl1 Ndc80 Cenpa"
  };

  // NOTE: only one collection runs per analysis (single background universe,
  // single BH/Bonferroni family). Combining collections is a deliberate future
  // decision: it would redefine the background universe AND make the multiple-
  // testing correction span every collection at once (inflating the family
  // size and the correction). Decide that model before implementing it.

  var worker = null;
  var workerFailed = false;
  var VER = "0";  // cache-bust token from manifest.version

  var el = {
    demoBanner: document.getElementById("demoBanner"),
    genes: document.getElementById("genes"),
    speciesToggle: document.getElementById("speciesToggle"),
    collection: document.getElementById("collection"),
    background: document.getElementById("background"),
    customBgWrap: document.getElementById("customBgWrap"),
    customBg: document.getElementById("customBg"),
    loadExample: document.getElementById("loadExample"),
    clearGenes: document.getElementById("clearGenes"),
    fdr: document.getElementById("fdr"),
    adjust: document.getElementById("adjust"),
    run: document.getElementById("run"),
    dlCsv: document.getElementById("dlCsv"),
    dlJson: document.getElementById("dlJson"),
    report: document.getElementById("report"),
    results: document.getElementById("results"),
    tbody: document.querySelector("#results tbody"),
    noResults: document.getElementById("noResults"),
    attrib: document.getElementById("attrib"),
    sizeMin: document.getElementById("sizeMin"),
    sizeMax: document.getElementById("sizeMax"),
    pageSize: document.getElementById("pageSize"),
    prevPage: document.getElementById("prevPage"),
    nextPage: document.getElementById("nextPage"),
    showingInfo: document.getElementById("showingInfo"),
    filterChips: document.getElementById("filterChips"),
    clearFilters: document.getElementById("clearFilters"),
    chartToggle: document.getElementById("chartToggle"),
    topN: document.getElementById("topN"),
    vizCaption: document.getElementById("vizCaption"),
    chartWrap: document.getElementById("chartWrap"),
    dlSvg: document.getElementById("dlSvg"),
    dlPng: document.getElementById("dlPng")
  };

  function fetchJson(path) {
    return fetch(path).then(function (r) {
      if (!r.ok) throw new Error("Failed to load " + path + " (" + r.status + ")");
      return r.json();
    });
  }

  function dataUrl(p) {
    return p + (p.indexOf("?") < 0 ? "?" : "&") + "v=" + VER;
  }

  function tokenize(text) {
    // Split on whitespace/comma/semicolon, then strip leading/trailing
    // punctuation/quotes/brackets only - internal characters (e.g. the hyphen
    // in HLA-DRB1 or NKX2-1) are preserved.
    return text.split(/[\s,;]+/)
      .map(function (t) { return t.replace(/^[^A-Za-z0-9]+/, "").replace(/[^A-Za-z0-9]+$/, ""); })
      .filter(Boolean);
  }

  function collectionsForSpecies(sp) {
    return state.manifest.collections.filter(function (c) { return c.species === sp; });
  }

  function populateCollections() {
    var cols = collectionsForSpecies(state.species);
    el.collection.innerHTML = "";
    cols.forEach(function (c) {
      var o = document.createElement("option");
      o.value = c.key;
      o.textContent = c.label + (c.available === false ? " (Stage 2)" : "");
      if (c.available === false) o.disabled = true;
      el.collection.appendChild(o);
    });
    var firstAvail = cols.find(function (c) { return c.available !== false; });
    if (firstAvail) {
      el.collection.value = firstAvail.key;
      state.collectionKey = firstAvail.key;
    }
  }

  function currentCollectionMeta() {
    return state.manifest.collections.find(function (c) {
      return c.species === state.species && c.key === state.collectionKey;
    });
  }

  function loadCollection() {
    var meta = currentCollectionMeta();
    var cacheKey = state.species + "/" + state.collectionKey;
    if (state.cache[cacheKey]) return Promise.resolve(state.cache[cacheKey]);
    var symPath = state.manifest.symbols[state.species];
    return Promise.all([fetchJson(dataUrl(symPath)), fetchJson(dataUrl(meta.path))]).then(function (res) {
      var bundle = {
        symbols: res[0].symbols,
        aliases: res[0].aliases || null,
        codingN: res[0].codingN || 0,
        terms: res[1].terms
      };
      state.cache[cacheKey] = bundle;
      return bundle;
    });
  }

  function setReport(html) { el.report.innerHTML = html; }

  // Fill the textarea with the species-appropriate example set (replaces any
  // current content). Does not auto-run.
  function loadExample() {
    el.genes.value = EXAMPLE[state.species] || EXAMPLE.human;
    setReport("");
    el.genes.focus();
  }

  // Reset the input and clear prior results, status, and charts (no reload).
  function clearGenes() {
    el.genes.value = "";
    state.lastResult = null;
    state.viewRows = [];
    state.page = 1;
    el.tbody.innerHTML = "";
    el.showingInfo.textContent = "";
    el.filterChips.innerHTML = "";
    el.noResults.classList.add("hidden");
    el.dlCsv.disabled = true;
    el.dlJson.disabled = true;
    setReport("");
    renderChart();           // lastResult null -> placeholder, disables chart export
    el.genes.focus();
  }

  // Lock the controls that define what ran (collection + species) plus Run,
  // so the dropdown/toggle can never disagree with the displayed results.
  function setRunning(on) {
    el.run.disabled = on;
    el.collection.disabled = on;
    Array.prototype.forEach.call(el.speciesToggle.children, function (b) { b.disabled = on; });
  }

  function run() {
    if (workerFailed) { onWorkerError(); return; }
    var tokens = tokenize(el.genes.value);
    if (tokens.length === 0) { setReport('<span class="warn">Paste at least one gene symbol.</span>'); return; }
    setRunning(true);
    var rm = currentCollectionMeta();
    state.lastRunMeta = { species: state.species, collection: rm ? rm.label : state.collectionKey };
    setReport("Loading collection and computing...");
    loadCollection().then(function (bundle) {
      var msg = {
        symbols: bundle.symbols,
        aliases: bundle.aliases,
        terms: bundle.terms,
        tokens: tokens,
        bgMode: el.background.value,
        codingN: bundle.codingN,
        customBg: el.background.value === "custom" ? tokenize(el.customBg.value) : null
      };
      worker.postMessage(msg);
    }).catch(function (err) {
      setRunning(false);
      setReport('<span class="warn">' + err.message + '</span>');
    });
  }

  // Worker script error or undeliverable message: turn a silent hang into a
  // visible, recoverable error and unlock the controls.
  function onWorkerError() {
    workerFailed = true;
    setRunning(false);
    setReport('<span class="warn">Analysis failed to run - please reload the page.</span>');
  }

  function onWorkerMessage(e) {
    setRunning(false);
    if (!e.data.ok) { setReport('<span class="warn">Error: ' + e.data.error + '</span>'); return; }
    state.lastResult = e.data.result;
    state.page = 1;
    renderReport(e.data.result);
    renderTable();
    el.dlCsv.disabled = false;
    el.dlJson.disabled = false;
  }

  function renderReport(res) {
    var bgLabel = { annotated: "annotated in collection", coding: "all protein-coding", custom: "custom" }[res.bgMode];
    var m = state.lastRunMeta || {};
    var sp = m.species ? m.species.charAt(0).toUpperCase() + m.species.slice(1) : "";
    var head = "<strong>" + sp + " | " + (m.collection || "") + "</strong> | ";
    var recog = '<span class="ok">recognized ' + res.recognized + " of " + res.uniqueGenes +
      " unique gene" + (res.uniqueGenes === 1 ? "" : "s") + "</span>";
    if (res.duplicates) {
      recog += '<span class="muted"> (' + res.duplicates + " duplicate token" +
        (res.duplicates === 1 ? "" : "s") + " ignored)</span>";
    }

    // Custom background that is empty or fully unrecognized -> N=0; surface it
    // instead of presenting a contradictory "n=0 | N=0" as a valid run.
    if (res.bgMode === "custom" && res.N === 0) {
      setReport(head + recog + '<br><span class="warn">Custom background has N=0: ' +
        (res.bgUnique ? "none of its " + res.bgUnique + " genes were recognized" : "it is empty") +
        ". No test was run. Paste a valid background gene list for the selected species, or " +
        "switch the Background option.</span>");
      return;
    }

    var html = head + recog + " | n=" + res.n + " in universe | background N=" + res.N + " (" + bgLabel + ")";
    if (res.bgMode === "custom") {
      html += '<span class="muted"> (custom background: ' + res.bgRecognized + " of " +
        res.bgUnique + " recognized)</span>";
    }
    if (res.dropped.length) {
      html += '<br><span class="warn">Unrecognized ' + res.dropped.length + ": " +
        res.dropped.slice(0, 25).join(", ") + (res.dropped.length > 25 ? " ..." : "") + "</span>";
    }
    setReport(html);
  }

  function passesThreshold(row) {
    var thr = parseFloat(el.fdr.value);
    var key = el.adjust.value; // "fdr" | "bonferroni"
    return row[key] <= thr;
  }

  function passesSize(row) {
    var mn = parseInt(el.sizeMin.value, 10);
    var mx = parseInt(el.sizeMax.value, 10);
    if (!isNaN(mn) && row.K < mn) return false;
    if (!isNaN(mx) && row.K > mx) return false;
    return true;
  }

  function fmtP(x) {
    if (x === 0) return "0";
    if (x < 1e-4) return x.toExponential(2);
    return x.toFixed(4);
  }

  // Full filtered + sorted result set (paginated only at render time).
  function buildView() {
    if (!state.lastResult) return [];
    var rows = state.lastResult.rows.filter(function (r) {
      return passesThreshold(r) && passesSize(r);
    });
    var k = state.sortKey, dir = state.sortDir;
    rows.sort(function (a, b) {
      var av = a[k], bv = b[k];
      if (typeof av === "string") return dir * av.localeCompare(bv);
      return dir * (av - bv);
    });
    return rows;
  }

  function genesCell(genes) {
    var full = genes.join(", ");
    if (genes.length <= GENE_PREVIEW) {
      return '<td class="genes">' + esc(full) + "</td>";
    }
    var preview = genes.slice(0, GENE_PREVIEW).join(", ");
    var extra = genes.length - GENE_PREVIEW;
    // Full list in title (hover) and in a hidden span revealed on expand;
    // export still uses the full r.genes array, never this truncation.
    return '<td class="genes" title="' + esc(full) + '">' +
      '<span class="gene-preview">' + esc(preview) + "</span>" +
      '<span class="gene-full">' + esc(full) + "</span>" +
      ' <button type="button" class="more" data-extra="' + extra + '">+' + extra + " more</button></td>";
  }

  function rowHtml(r) {
    return "<tr>" +
      '<td class="term" title="' + esc(r.name) + '">' + esc(r.name) + "</td>" +
      '<td class="ns">' + esc(r.namespace) + "</td>" +
      '<td class="num">' + r.K + "</td>" +
      '<td class="num">' + r.k + "</td>" +
      '<td class="num">' + r.fold.toFixed(2) + "</td>" +
      '<td class="num">' + fmtP(r.p) + "</td>" +
      '<td class="num">' + fmtP(r.fdr) + "</td>" +
      '<td class="num">' + fmtP(r.bonferroni) + "</td>" +
      genesCell(r.genes) +
      "</tr>";
  }

  function renderTable() {
    if (!state.lastResult) return;
    var rows = buildView();
    state.viewRows = rows;
    var total = rows.length;
    el.noResults.classList.toggle("hidden", total > 0);

    var size = state.pageSize === "all" ? total : state.pageSize;
    var pages = size > 0 ? Math.ceil(total / size) : 1;
    if (state.page > pages) state.page = pages || 1;
    if (state.page < 1) state.page = 1;
    var start = size > 0 ? (state.page - 1) * size : 0;
    var end = size > 0 ? Math.min(start + size, total) : total;

    var html = "";
    for (var i = start; i < end; i++) html += rowHtml(rows[i]);
    el.tbody.innerHTML = html;

    var shown = total === 0 ? "0 of 0 terms" :
      (start + 1) + "-" + end + " of " + total + " terms";
    el.showingInfo.textContent = "Showing " + shown +
      (pages > 1 ? "  (page " + state.page + " of " + pages + ")" : "");
    el.prevPage.disabled = state.page <= 1;
    el.nextPage.disabled = state.page >= pages;

    renderSortIndicators();
    renderActiveFilters();
    renderChart();
  }

  function renderSortIndicators() {
    Array.prototype.forEach.call(el.results.querySelectorAll("th[data-sort]"), function (th) {
      th.classList.remove("sort-asc", "sort-desc");
      if (th.dataset.sort === state.sortKey) {
        th.classList.add(state.sortDir === 1 ? "sort-asc" : "sort-desc");
      }
    });
  }

  var SORT_LABELS = {
    name: "Term", namespace: "Namespace", K: "Size", k: "Overlap", fold: "Fold",
    p: "P", fdr: "FDR", bonferroni: "Bonferroni"
  };

  function renderActiveFilters() {
    var chips = [];
    var sigLabel = el.adjust.value === "bonferroni" ? "Bonferroni" : "FDR";
    chips.push(sigLabel + " < " + (parseFloat(el.fdr.value) || 0));

    var mn = parseInt(el.sizeMin.value, 10);
    var mx = parseInt(el.sizeMax.value, 10);
    if (!isNaN(mn) && !isNaN(mx)) chips.push("size " + mn + "-" + mx);
    else if (!isNaN(mn)) chips.push("size >= " + mn);
    else if (!isNaN(mx)) chips.push("size <= " + mx);

    chips.push("sorted by " + (SORT_LABELS[state.sortKey] || state.sortKey) +
      " " + (state.sortDir === 1 ? "asc" : "desc"));
    chips.push((state.pageSize === "all" ? "all" : state.pageSize) + " / page");

    el.filterChips.innerHTML = chips.map(function (c) {
      return '<span class="chip">' + esc(c) + "</span>";
    }).join("");
  }

  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
    });
  }

  // ---- charts (hand-rolled SVG, no dependency) ----
  var SVGNS = "http://www.w3.org/2000/svg";

  function svgEl(tag, attrs, text) {
    var e = document.createElementNS(SVGNS, tag);
    for (var k in attrs) if (Object.prototype.hasOwnProperty.call(attrs, k)) e.setAttribute(k, attrs[k]);
    if (text != null) e.appendChild(document.createTextNode(text));
    return e;
  }

  function withTitle(node, full) {
    node.appendChild(svgEl("title", {}, full));
    return node;
  }

  function nlog10(x) { return -Math.log10(Math.max(x, 1e-300)); }

  // Word-wrap a term name into up to maxLines lines of ~maxChars each; the last
  // line is ellipsized only if it still overflows (full name stays in tooltip).
  function wrapLabel(name, maxChars, maxLines) {
    var words = name.split(" ");
    var lines = [], cur = "";
    for (var i = 0; i < words.length; i++) {
      var test = cur ? cur + " " + words[i] : words[i];
      if (test.length > maxChars && cur) {
        lines.push(cur);
        if (lines.length === maxLines - 1) {
          var rest = words.slice(i).join(" ");
          lines.push(rest.length > maxChars ? rest.slice(0, maxChars - 1) + "..." : rest);
          return lines;
        }
        cur = words[i];
      } else {
        cur = test;
      }
    }
    if (cur) lines.push(cur);
    return lines;
  }

  function sigLabelText() {
    var thr = parseFloat(el.fdr.value) || 0;
    var name = el.adjust.value === "bonferroni" ? "Bonferroni" : "FDR";
    return name + " < " + thr;
  }

  // Filtered (threshold + size) rows, sorted by significance ascending, for
  // charts. Independent of the table's current sort column.
  function chartRows() {
    if (!state.lastResult) return [];
    var rows = state.lastResult.rows.filter(function (r) {
      return passesThreshold(r) && passesSize(r);
    });
    rows.sort(function (a, b) { return a.p - b.p; });
    return rows;
  }

  function newSvg(w, h) {
    var svg = svgEl("svg", {
      xmlns: SVGNS, viewBox: "0 0 " + w + " " + h, width: w, height: h,
      "font-family": "-apple-system, Segoe UI, Roboto, sans-serif"
    });
    svg.appendChild(svgEl("rect", { x: 0, y: 0, width: w, height: h, fill: "#ffffff" }));
    return svg;
  }

  function axisLabels(rows, svg, marginLeft, marginTop, rowH) {
    // Chars that fit the gutter at ~11px; wrap to at most 2 lines.
    var maxChars = Math.max(14, Math.floor((marginLeft - 16) / 5.8));
    rows.forEach(function (r, i) {
      var cy = marginTop + i * rowH + rowH / 2;
      var lines = wrapLabel(r.name, maxChars, 2);
      var t = svgEl("text", {
        x: marginLeft - 8, "text-anchor": "end", "dominant-baseline": "middle",
        "font-size": 11, fill: "#1c2330"
      });
      var startY = cy - (lines.length - 1) * 6;
      lines.forEach(function (ln, j) {
        t.appendChild(svgEl("tspan", { x: marginLeft - 8, y: startY + j * 12 }, ln));
      });
      withTitle(t, r.name + " (" + r.namespace + ")");
      svg.appendChild(t);
    });
  }

  function buildBarSvg(rows, W) {
    var mL = 250, mR = 56, mT = 18, mB = 34, rowH = 26, barH = 14;
    var H = mT + rows.length * rowH + mB;
    var plotW = W - mL - mR;
    var svg = newSvg(W, H);
    var xmax = Math.max.apply(null, rows.map(function (r) { return nlog10(r.fdr); })) || 1;

    // x gridlines + ticks
    var ticks = 4;
    for (var t = 0; t <= ticks; t++) {
      var xv = xmax * t / ticks;
      var x = mL + (xv / xmax) * plotW;
      svg.appendChild(svgEl("line", { x1: x, y1: mT, x2: x, y2: mT + rows.length * rowH, stroke: "#eceff5" }));
      svg.appendChild(svgEl("text", { x: x, y: H - mB + 16, "text-anchor": "middle", "font-size": 10, fill: "#66718a" }, xv.toFixed(1)));
    }
    svg.appendChild(svgEl("text", { x: mL + plotW / 2, y: H - 4, "text-anchor": "middle", "font-size": 11, fill: "#66718a" }, "-log10(FDR)"));

    rows.forEach(function (r, i) {
      var y = mT + i * rowH;
      var val = nlog10(r.fdr);
      var len = (val / xmax) * plotW;
      var bar = svgEl("rect", { x: mL, y: y + (rowH - barH) / 2, width: Math.max(len, 0.5), height: barH, fill: "#2d6cdf", rx: 2 });
      withTitle(bar, r.name + "  FDR=" + r.fdr.toExponential(2) + "  overlap=" + r.k + "/" + r.K);
      svg.appendChild(bar);
    });
    axisLabels(rows, svg, mL, mT, rowH);
    return svg;
  }

  // 3-stop YlOrRd-ish scale; t in [0,1], 1 = most significant.
  function fdrColor(t) {
    var stops = [[255, 237, 160], [254, 178, 76], [227, 26, 28]];
    var seg = t >= 1 ? 1 : t * (stops.length - 1);
    var i = Math.min(Math.floor(seg), stops.length - 2);
    var f = seg - i, a = stops[i], b = stops[i + 1];
    var c = a.map(function (av, j) { return Math.round(av + f * (b[j] - av)); });
    return "rgb(" + c[0] + "," + c[1] + "," + c[2] + ")";
  }

  function buildDotSvg(rows, W) {
    var mL = 250, mR = 150, mT = 18, mB = 40, rowH = 26;
    var H = mT + rows.length * rowH + mB;
    var plotW = W - mL - mR;
    var svg = newSvg(W, H);

    var xs = rows.map(function (r) { return Math.log10(Math.max(r.fold, 0.001)); });
    var xmin = Math.min.apply(null, xs), xmax = Math.max.apply(null, xs);
    if (xmin === xmax) { xmin -= 0.5; xmax += 0.5; }
    var pad = (xmax - xmin) * 0.08; xmin -= pad; xmax += pad;
    var xpix = function (lv) { return mL + (lv - xmin) / (xmax - xmin) * plotW; };

    var ks = rows.map(function (r) { return r.k; });
    var kmin = Math.min.apply(null, ks), kmax = Math.max.apply(null, ks);
    var rad = function (k) { return kmax === kmin ? 7 : 4 + (k - kmin) / (kmax - kmin) * 8; };

    var vs = rows.map(function (r) { return nlog10(r.fdr); });
    var vmin = Math.min.apply(null, vs), vmax = Math.max.apply(null, vs);
    var colorT = function (v) { return vmax === vmin ? 1 : (v - vmin) / (vmax - vmin); };

    // x gridlines + ticks (labelled in fold-enrichment units)
    var ticks = 4;
    for (var t = 0; t <= ticks; t++) {
      var lv = xmin + (xmax - xmin) * t / ticks;
      var x = xpix(lv);
      svg.appendChild(svgEl("line", { x1: x, y1: mT, x2: x, y2: mT + rows.length * rowH, stroke: "#eceff5" }));
      var fold = Math.pow(10, lv);
      svg.appendChild(svgEl("text", { x: x, y: H - mB + 16, "text-anchor": "middle", "font-size": 10, fill: "#66718a" },
        fold >= 10 ? String(Math.round(fold)) : fold.toFixed(1)));
    }
    svg.appendChild(svgEl("text", { x: mL + plotW / 2, y: H - 6, "text-anchor": "middle", "font-size": 11, fill: "#66718a" }, "fold enrichment (log scale)"));

    rows.forEach(function (r, i) {
      var y = mT + i * rowH + rowH / 2;
      var dot = svgEl("circle", {
        cx: xpix(Math.log10(Math.max(r.fold, 0.001))), cy: y, r: rad(r.k),
        fill: fdrColor(colorT(nlog10(r.fdr))), stroke: "#7a1d1d", "stroke-width": 0.5
      });
      withTitle(dot, r.name + "  fold=" + r.fold.toFixed(1) + "  overlap=" + r.k + "/" + r.K + "  FDR=" + r.fdr.toExponential(2));
      svg.appendChild(dot);
    });
    axisLabels(rows, svg, mL, mT, rowH);

    // legends in the right margin: FDR color gradient + overlap size
    var lx = W - mR + 24, ly = mT + 6;
    svg.appendChild(svgEl("text", { x: lx, y: ly, "font-size": 10, fill: "#66718a" }, "FDR"));
    var gradH = 70;
    for (var g = 0; g < gradH; g++) {
      var tt = 1 - g / gradH;
      svg.appendChild(svgEl("rect", { x: lx, y: ly + 6 + g, width: 12, height: 1, fill: fdrColor(tt) }));
    }
    svg.appendChild(svgEl("text", { x: lx + 16, y: ly + 12, "font-size": 9, fill: "#66718a" }, "most sig"));
    svg.appendChild(svgEl("text", { x: lx + 16, y: ly + 6 + gradH, "font-size": 9, fill: "#66718a" }, "least"));
    var sy = ly + 6 + gradH + 24;
    svg.appendChild(svgEl("text", { x: lx, y: sy - 8, "font-size": 10, fill: "#66718a" }, "overlap"));
    [kmin, kmax].forEach(function (k, j) {
      if (j === 1 && kmax === kmin) return;
      var cy = sy + 8 + j * 26;
      svg.appendChild(svgEl("circle", { cx: lx + 8, cy: cy, r: rad(k), fill: "#ccd3e0", stroke: "#7a1d1d", "stroke-width": 0.5 }));
      svg.appendChild(svgEl("text", { x: lx + 24, y: cy + 3, "font-size": 9, fill: "#66718a" }, String(k)));
    });
    return svg;
  }

  function renderChart() {
    if (!state.lastResult) {
      el.chartWrap.innerHTML = "";
      el.vizCaption.textContent = "Run an analysis to see charts.";
      state.currentSvg = null;
      el.dlSvg.disabled = true; el.dlPng.disabled = true;
      return;
    }
    var rows = chartRows();
    var M = rows.length;
    if (M === 0) {
      el.chartWrap.innerHTML = '<p class="empty-chart">No terms pass ' + esc(sigLabelText()) +
        ". Adjust the threshold or filters.</p>";
      el.vizCaption.textContent = "";
      state.currentSvg = null;
      el.dlSvg.disabled = true; el.dlPng.disabled = true;
      return;
    }
    var top = rows.slice(0, state.topN);
    // Size the chart to the panel width (clamped) so it fills the space on wide
    // screens and scrolls inside its container on narrow ones; high enough
    // intrinsic width keeps SVG/PNG export crisp.
    var cw = el.chartWrap.clientWidth || 680;
    var W = Math.max(620, Math.min(cw, 1100));
    var svg = state.chartType === "bar" ? buildBarSvg(top, W) : buildDotSvg(top, W);
    el.chartWrap.innerHTML = "";
    el.chartWrap.appendChild(svg);
    state.currentSvg = svg;
    el.vizCaption.textContent = "Showing top " + top.length + " of " + M +
      " significant terms (" + sigLabelText() + "), ordered by significance.";
    el.dlSvg.disabled = false; el.dlPng.disabled = false;
  }

  function serializeSvg(svg) {
    var s = new XMLSerializer().serializeToString(svg);
    if (s.indexOf("xmlns=") === -1) s = s.replace("<svg", '<svg xmlns="' + SVGNS + '"');
    return s;
  }

  function downloadSvg() {
    if (!state.currentSvg) return;
    saveBlob(serializeSvg(state.currentSvg), "image/svg+xml", "enrichlite_" + state.chartType + ".svg");
  }

  function downloadPng() {
    if (!state.currentSvg) return;
    var svg = state.currentSvg;
    var W = +svg.getAttribute("width"), H = +svg.getAttribute("height"), scale = 2;
    var data = "data:image/svg+xml;base64," + btoa(unescape(encodeURIComponent(serializeSvg(svg))));
    var img = new Image();
    img.onload = function () {
      var c = document.createElement("canvas");
      c.width = W * scale; c.height = H * scale;
      var ctx = c.getContext("2d");
      ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, c.width, c.height);
      ctx.drawImage(img, 0, 0, c.width, c.height);
      c.toBlob(function (b) {
        var a = document.createElement("a");
        a.href = URL.createObjectURL(b);
        a.download = "enrichlite_" + state.chartType + ".png";
        a.click();
        URL.revokeObjectURL(a.href);
      }, "image/png");
    };
    img.src = data;
  }

  // Export uses the full filtered + sorted view (all pages, not just the
  // current page).
  function visibleRows() {
    return buildView();
  }

  function downloadCsv() {
    var head = ["term", "namespace", "set_size_K", "overlap_k", "query_n", "background_N", "fold", "p", "fdr", "bonferroni", "genes"];
    var lines = [head.join(",")];
    visibleRows().forEach(function (r) {
      var cells = [r.name, r.namespace, r.K, r.k, r.n, r.N, r.fold, r.p, r.fdr, r.bonferroni, r.genes.join(" ")];
      lines.push(cells.map(csvCell).join(","));
    });
    saveBlob(lines.join("\n"), "text/csv", "enrichlite_results.csv");
  }

  function csvCell(v) {
    var s = String(v);
    if (/[",\n]/.test(s)) return '"' + s.replace(/"/g, '""') + '"';
    return s;
  }

  function downloadJson() {
    var payload = {
      species: state.species,
      collection: state.collectionKey,
      background: state.lastResult.bgMode,
      N: state.lastResult.N,
      n: state.lastResult.n,
      results: visibleRows()
    };
    saveBlob(JSON.stringify(payload, null, 2), "application/json", "enrichlite_results.json");
  }

  function saveBlob(text, type, name) {
    var blob = new Blob([text], { type: type });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = name;
    a.click();
    URL.revokeObjectURL(a.href);
  }

  function renderAttribution() {
    var s = state.manifest.sources || {};
    var parts = [];
    if (s.go) parts.push("GO release " + (s.go.release || "n/a") + (s.go.doi ? " (DOI: " + s.go.doi + ")" : ""));
    if (s.msigdb) parts.push("MSigDB " + (s.msigdb.version || "n/a"));
    if (s.reactome) parts.push("Reactome " + (s.reactome.version || "n/a"));
    el.attrib.textContent = parts.length ? "Data: " + parts.join(" | ") : "";
  }

  function wireEvents() {
    el.speciesToggle.addEventListener("click", function (e) {
      var b = e.target.closest("button[data-species]");
      if (!b) return;
      Array.prototype.forEach.call(el.speciesToggle.children, function (c) { c.classList.remove("active"); });
      b.classList.add("active");
      state.species = b.dataset.species;
      populateCollections();
    });
    el.collection.addEventListener("change", function () { state.collectionKey = el.collection.value; });
    el.background.addEventListener("change", function () {
      el.customBgWrap.classList.toggle("hidden", el.background.value !== "custom");
    });
    el.run.addEventListener("click", run);
    el.loadExample.addEventListener("click", loadExample);
    el.clearGenes.addEventListener("click", clearGenes);
    // Filter/sort changes reset to page 1; pagination keeps the page.
    var repage = function () { state.page = 1; renderTable(); };
    el.fdr.addEventListener("change", repage);
    el.adjust.addEventListener("change", repage);
    el.sizeMin.addEventListener("input", repage);
    el.sizeMax.addEventListener("input", repage);
    el.pageSize.addEventListener("change", function () {
      var v = el.pageSize.value;
      state.pageSize = v === "all" ? "all" : parseInt(v, 10);
      repage();
    });
    el.prevPage.addEventListener("click", function () { state.page -= 1; renderTable(); });
    el.nextPage.addEventListener("click", function () { state.page += 1; renderTable(); });
    el.clearFilters.addEventListener("click", function () {
      // Reset size range, sort, and page size to defaults. Does not touch the
      // FDR threshold/correction and does not re-run the worker.
      el.sizeMin.value = "";
      el.sizeMax.value = "";
      el.pageSize.value = "25";
      state.pageSize = 25;
      state.sortKey = "p";
      state.sortDir = 1;
      state.page = 1;
      renderTable();
    });
    // Expand/collapse a truncated genes cell.
    el.tbody.addEventListener("click", function (e) {
      var btn = e.target.closest(".more");
      if (!btn) return;
      var td = btn.closest("td.genes");
      var expanded = td.classList.toggle("expanded");
      btn.textContent = expanded ? "show less" : "+" + btn.dataset.extra + " more";
    });
    el.dlCsv.addEventListener("click", downloadCsv);
    el.dlJson.addEventListener("click", downloadJson);
    el.chartToggle.addEventListener("click", function (e) {
      var b = e.target.closest("button[data-chart]");
      if (!b) return;
      Array.prototype.forEach.call(el.chartToggle.children, function (c) { c.classList.remove("active"); });
      b.classList.add("active");
      state.chartType = b.dataset.chart;
      renderChart();
    });
    el.topN.addEventListener("change", function () {
      state.topN = parseInt(el.topN.value, 10);
      renderChart();
    });
    el.dlSvg.addEventListener("click", downloadSvg);
    el.dlPng.addEventListener("click", downloadPng);
    // Re-fit the chart to the panel width on resize (debounced).
    var rzT;
    window.addEventListener("resize", function () {
      clearTimeout(rzT);
      rzT = setTimeout(function () { if (state.lastResult) renderChart(); }, 150);
    });
    Array.prototype.forEach.call(el.results.querySelectorAll("th[data-sort]"), function (th) {
      th.addEventListener("click", function () {
        var key = th.dataset.sort;
        if (state.sortKey === key) state.sortDir *= -1;
        else { state.sortKey = key; state.sortDir = 1; }
        state.page = 1;
        renderTable();
      });
    });
  }

  // Make the "?" helper badges reachable by tap and keyboard (not just hover).
  // The native title is kept for desktop hover; a popover carries the same text
  // for click/tap/focus, with screen-reader association and viewport clamping.
  function wireHelpBadges() {
    var pop = document.createElement("div");
    pop.className = "help-pop hidden";
    pop.id = "helpPopover";
    pop.setAttribute("role", "tooltip");
    document.body.appendChild(pop);
    var openFor = null;

    function place(badge) {
      var r = badge.getBoundingClientRect();
      var margin = 8;
      var docW = document.documentElement.clientWidth;
      var docH = document.documentElement.clientHeight;
      var pw = pop.offsetWidth, ph = pop.offsetHeight;
      var left = r.left;
      left = Math.min(left, docW - pw - margin);
      left = Math.max(margin, left);
      var top = r.bottom + 6;
      if (docH - r.bottom < ph + 14) top = r.top - ph - 6;  // flip above near bottom
      pop.style.left = (left + window.scrollX) + "px";
      pop.style.top = (top + window.scrollY) + "px";
    }
    function close() {
      if (!openFor) return;
      pop.classList.add("hidden");
      openFor.setAttribute("aria-expanded", "false");
      openFor.removeAttribute("aria-describedby");
      openFor = null;
    }
    function open(badge) {
      if (openFor && openFor !== badge) {
        openFor.setAttribute("aria-expanded", "false");
        openFor.removeAttribute("aria-describedby");
      }
      pop.textContent = badge.getAttribute("data-help") || "";
      pop.classList.remove("hidden");
      badge.setAttribute("aria-expanded", "true");
      badge.setAttribute("aria-describedby", "helpPopover");
      openFor = badge;
      place(badge);
    }

    Array.prototype.forEach.call(document.querySelectorAll(".helpq"), function (b) {
      b.setAttribute("data-help", b.getAttribute("title") || "");
      b.setAttribute("tabindex", "0");
      b.setAttribute("role", "button");
      b.setAttribute("aria-label", "Help");
      b.setAttribute("aria-expanded", "false");
      b.addEventListener("click", function (e) {
        e.preventDefault();        // don't activate the wrapping label/control
        e.stopPropagation();
        if (openFor === b) close(); else open(b);
      });
      b.addEventListener("keydown", function (e) {
        if (e.key === "Enter" || e.key === " " || e.key === "Spacebar") {
          e.preventDefault();
          if (openFor === b) close(); else open(b);
        } else if (e.key === "Escape" && openFor === b) {
          close(); b.focus();
        }
      });
    });

    document.addEventListener("click", function (e) {
      if (openFor && !pop.contains(e.target) && !openFor.contains(e.target)) close();
    });
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") close(); });
    window.addEventListener("resize", function () { if (openFor) place(openFor); });
    window.addEventListener("scroll", function () { if (openFor) place(openFor); }, true);
  }

  function init(m) {
    state.manifest = m;
    VER = encodeURIComponent(String(m.version || "0"));
    worker = new Worker("js/worker.js?v=" + VER);
    worker.onmessage = onWorkerMessage;
    worker.onerror = onWorkerError;
    worker.onmessageerror = onWorkerError;
    el.demoBanner.classList.toggle("hidden", !m.demo);
    populateCollections();
    renderAttribution();
    wireEvents();
    wireHelpBadges();
    renderChart();
  }

  // Manifest is preloaded by the inline bootstrap in index.html (fetched with
  // no-cache). Fall back to fetching it if app.js was loaded standalone.
  if (window.__ENRICH_MANIFEST__) {
    init(window.__ENRICH_MANIFEST__);
  } else {
    fetchJson("data/manifest.json").then(init).catch(function (err) {
      setReport('<span class="warn">Could not load data/manifest.json: ' + err.message + "</span>");
    });
  }
})();
