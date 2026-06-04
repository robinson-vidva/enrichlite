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
    pageSize: 50        // number, or "all"
  };

  var worker = null;
  var VER = "0";  // cache-bust token from manifest.version

  var el = {
    demoBanner: document.getElementById("demoBanner"),
    genes: document.getElementById("genes"),
    speciesToggle: document.getElementById("speciesToggle"),
    collection: document.getElementById("collection"),
    background: document.getElementById("background"),
    customBgWrap: document.getElementById("customBgWrap"),
    customBg: document.getElementById("customBg"),
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
    showingInfo: document.getElementById("showingInfo")
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
    return text.split(/[\s,;]+/).map(function (t) { return t.trim(); }).filter(Boolean);
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

  function run() {
    var tokens = tokenize(el.genes.value);
    if (tokens.length === 0) { setReport('<span class="warn">Paste at least one gene symbol.</span>'); return; }
    el.run.disabled = true;
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
      el.run.disabled = false;
      setReport('<span class="warn">' + err.message + '</span>');
    });
  }

  function onWorkerMessage(e) {
    el.run.disabled = false;
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
    var html = head + '<span class="ok">recognized ' + res.recognized + "/" + res.queryTotal + "</span>";
    html += " | n=" + res.n + " in universe | background N=" + res.N + " (" + bgLabel + ")";
    if (res.dropped.length) {
      html += '<br><span class="warn">Dropped ' + res.dropped.length + ": " +
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

  function rowHtml(r) {
    var genes = r.genes.join(", ");
    return "<tr>" +
      '<td class="term" title="' + esc(r.name) + '">' + esc(r.name) + "</td>" +
      '<td class="ns">' + esc(r.namespace) + "</td>" +
      '<td class="num">' + r.K + "</td>" +
      '<td class="num">' + r.k + "</td>" +
      '<td class="num">' + r.fold.toFixed(2) + "</td>" +
      '<td class="num">' + fmtP(r.p) + "</td>" +
      '<td class="num">' + fmtP(r.fdr) + "</td>" +
      '<td class="num">' + fmtP(r.bonferroni) + "</td>" +
      '<td class="genes" title="' + esc(genes) + '">' + esc(genes) + "</td>" +
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
  }

  function renderSortIndicators() {
    Array.prototype.forEach.call(el.results.querySelectorAll("th[data-sort]"), function (th) {
      th.classList.remove("sort-asc", "sort-desc");
      if (th.dataset.sort === state.sortKey) {
        th.classList.add(state.sortDir === 1 ? "sort-asc" : "sort-desc");
      }
    });
  }

  function esc(s) {
    return String(s).replace(/[&<>"]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c];
    });
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
    el.dlCsv.addEventListener("click", downloadCsv);
    el.dlJson.addEventListener("click", downloadJson);
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

  function init(m) {
    state.manifest = m;
    VER = encodeURIComponent(String(m.version || "0"));
    worker = new Worker("js/worker.js?v=" + VER);
    worker.onmessage = onWorkerMessage;
    el.demoBanner.classList.toggle("hidden", !m.demo);
    populateCollections();
    renderAttribution();
    wireEvents();
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
