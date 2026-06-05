// worker.js - runs ORA off the main thread.
var _v = new URLSearchParams(self.location.search).get("v");
importScripts("stats.js" + (_v ? "?v=" + encodeURIComponent(_v) : ""));

var S = self.EnrichStats;

// Normalize a symbol for within-species matching (case-fold only).
function norm(sym) {
  return sym.toUpperCase();
}

// Build a lookup from normalized symbol -> gene index, including aliases.
function buildIndex(symbols, aliases) {
  var map = new Map();
  for (var i = 0; i < symbols.length; i++) map.set(norm(symbols[i]), i);
  if (aliases) {
    for (var a in aliases) {
      if (Object.prototype.hasOwnProperty.call(aliases, a)) {
        var target = aliases[a];
        if (!map.has(norm(a)) && map.has(norm(symbols[target]))) {
          map.set(norm(a), target);
        }
      }
    }
  }
  return map;
}

// Resolve pasted tokens to gene indices within this species universe.
// De-duplicates by normalized symbol so repeated tokens are not counted as
// separate genes or conflated with unrecognized ones.
function resolveQuery(tokens, symIndex) {
  var seen = new Set();        // normalized tokens already counted
  var recognized = new Set();  // matched gene indices
  var dropped = [];            // unique unrecognized tokens (original form)
  var duplicates = 0;
  for (var i = 0; i < tokens.length; i++) {
    var t = tokens[i];
    if (!t) continue;
    var key = norm(t);
    if (seen.has(key)) { duplicates++; continue; }
    seen.add(key);
    var gi = symIndex.get(key);
    if (gi === undefined) dropped.push(t);
    else recognized.add(gi);
  }
  return { recognized: recognized, dropped: dropped, duplicates: duplicates, unique: seen.size };
}

function run(msg) {
  var symbols = msg.symbols;
  var aliases = msg.aliases;
  var terms = msg.terms;
  var tokens = msg.tokens;
  var bgMode = msg.bgMode;            // "annotated" | "coding" | "custom"
  var codingN = msg.codingN || 0;     // protein-coding count for species
  var customBg = msg.customBg;        // array of tokens for custom background

  var symIndex = buildIndex(symbols, aliases);

  // Universe = set of gene indices that count toward N and toward matching.
  var universe = null;
  if (bgMode === "annotated") {
    universe = new Set();
    for (var ti = 0; ti < terms.length; ti++) {
      var g = terms[ti].genes;
      for (var gi = 0; gi < g.length; gi++) universe.add(g[gi]);
    }
  }
  var bgRecognized = null, bgUnique = null;
  if (bgMode === "custom") {
    var rc = resolveQuery(customBg || [], symIndex);
    universe = rc.recognized;
    bgRecognized = rc.unique - rc.dropped.length;
    bgUnique = rc.unique;
  }
  // "coding" uses codingN directly; universe stays null (no membership filter
  // beyond the symbol table).

  var q = resolveQuery(tokens, symIndex);
  var queryGenes = q.recognized;

  // Restrict query to the universe when one is defined.
  var queryInUni = new Set();
  queryGenes.forEach(function (gi) {
    if (!universe || universe.has(gi)) queryInUni.add(gi);
  });

  var N;
  if (bgMode === "coding") N = codingN;
  else N = universe.size;

  var n = queryInUni.size;

  var lnFact = S.makeLnFactorial(N + 1);

  var rows = [];
  for (var t2 = 0; t2 < terms.length; t2++) {
    var term = terms[t2];
    var tg = term.genes;
    // K = term genes within universe.
    var K = 0;
    var overlapIdx = [];
    var inUni;
    if (universe) {
      K = 0;
      for (var j = 0; j < tg.length; j++) {
        if (universe.has(tg[j])) {
          K++;
          if (queryInUni.has(tg[j])) overlapIdx.push(tg[j]);
        }
      }
    } else {
      // coding background: term genes all count, overlap vs query.
      K = tg.length;
      var tgSet = new Set(tg);
      queryInUni.forEach(function (gi) { if (tgSet.has(gi)) overlapIdx.push(gi); });
    }
    var k = overlapIdx.length;
    if (k === 0) continue;

    var p = S.hyperUpperTail(k, n, K, N, lnFact);
    var fold = S.foldEnrichment(k, n, K, N);
    rows.push({
      id: term.id,
      name: term.name,
      namespace: term.namespace || "",
      K: K,
      k: k,
      n: n,
      N: N,
      fold: fold,
      p: p,
      genes: overlapIdx.map(function (gi) { return symbols[gi]; })
    });
  }

  var pvals = rows.map(function (r) { return r.p; });
  var qvals = S.bhFDR(pvals);
  var bvals = S.bonferroni(pvals);
  for (var r = 0; r < rows.length; r++) {
    rows[r].fdr = qvals[r];
    rows[r].bonferroni = bvals[r];
  }

  rows.sort(function (a, b) { return a.p - b.p; });

  return {
    rows: rows,
    N: N,
    n: n,
    uniqueGenes: q.unique,
    recognized: q.unique - q.dropped.length,
    duplicates: q.duplicates,
    dropped: q.dropped,
    bgMode: bgMode,
    bgRecognized: bgRecognized,
    bgUnique: bgUnique
  };
}

self.onmessage = function (e) {
  try {
    var res = run(e.data);
    self.postMessage({ ok: true, result: res });
  } catch (err) {
    self.postMessage({ ok: false, error: String(err && err.message || err) });
  }
};
