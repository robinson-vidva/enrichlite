// test_data.js - guards on the built data + per-collection universe isolation.
// Run: node tests/test_data.js
var fs = require("fs");
var path = require("path");
var S = require("../js/stats.js");

var ROOT = path.join(__dirname, "..");
function L(p) { return JSON.parse(fs.readFileSync(path.join(ROOT, p), "utf8")); }

var fails = 0;
function ok(name, cond, extra) {
  console.log((cond ? "PASS " : "FAIL ") + name + (extra ? " (" + extra + ")" : ""));
  if (!cond) fails++;
}

function universe(coll) {
  var u = new Set();
  coll.terms.forEach(function (t) { t.genes.forEach(function (g) { u.add(g); }); });
  return u;
}

// Skip if still on demo data.
var man = L("data/manifest.json");
if (man.demo) { console.log("SKIP: manifest is demo data"); process.exit(0); }

var hsym = L("data/human/symbols.json").symbols;
var hidx = {}; hsym.forEach(function (s, i) { hidx[s] = i; });

var hReact = L("data/human/reactome.json");
var hHall = L("data/human/hallmark.json");
var rUni = universe(hReact);
var hUni = universe(hHall);

// 1) Reactome universe is in a sane range, well above Hallmark's ~4.4k.
ok("human Reactome universe > 9000", rUni.size > 9000, "size=" + rUni.size);
ok("human Hallmark universe ~4.4k (3000-6000)", hUni.size > 3000 && hUni.size < 6000, "size=" + hUni.size);
ok("Reactome universe >> Hallmark (>2x)", rUni.size > 2 * hUni.size,
   rUni.size + " vs " + hUni.size);

// 2) Query genes present where expected (documents the n=2 vs n=3 symptom).
["BRAF", "KRAS", "MAPK1"].forEach(function (g) {
  ok("Reactome universe contains " + g, hidx[g] !== undefined && rUni.has(hidx[g]));
});
ok("Hallmark universe MISSING KRAS (explains n=2 fingerprint)", !hUni.has(hidx["KRAS"]));

// 3) Manifest N matches each collection file's own computed universe (per species).
man.collections.filter(function (c) { return c.available; }).forEach(function (c) {
  var coll = L(c.path);
  var n = universe(coll).size;
  ok("manifest N matches universe for " + c.species + "/" + c.key,
     c.N === n, "manifest=" + c.N + " computed=" + n);
});

// 4) Background switches with the collection: same query, different N and n.
function runOra(coll, querySyms, symbols) {
  var idx = {}; symbols.forEach(function (s, i) { idx[s] = i; });
  var uni = universe(coll);
  var N = uni.size;
  var q = new Set();
  querySyms.forEach(function (s) { var gi = idx[s.toUpperCase()]; if (gi !== undefined && uni.has(gi)) q.add(gi); });
  var n = q.size;
  var lnFact = S.makeLnFactorial(N + 1);
  var rows = [];
  coll.terms.forEach(function (t) {
    var K = t.genes.length, k = 0;
    t.genes.forEach(function (g) { if (q.has(g)) k++; });
    if (k === 0) return;
    rows.push({ name: t.name, K: K, k: k, p: S.hyperUpperTail(k, n, K, N, lnFact) });
  });
  var fdr = S.bhFDR(rows.map(function (r) { return r.p; }));
  rows.forEach(function (r, i) { r.fdr = fdr[i]; });
  return { N: N, n: n, rows: rows };
}

var query = ["BRAF", "KRAS", "MAPK1"];
var rRes = runOra(hReact, query, hsym);
var hRes = runOra(hHall, query, hsym);
ok("Reactome run n=3", rRes.n === 3, "n=" + rRes.n);
ok("Reactome run N matches manifest", rRes.N === 11963, "N=" + rRes.N);
ok("Hallmark run n=2 (KRAS not in its universe)", hRes.n === 2, "n=" + hRes.n);
ok("collection isolation: Reactome N != Hallmark N", rRes.N !== hRes.N, rRes.N + " vs " + hRes.N);

// 5) Acceptance: a RAF/MAPK/ERK term is significant at FDR < 0.05 on Reactome.
var sig = rRes.rows.filter(function (r) { return r.fdr < 0.05; });
var mapk = sig.filter(function (r) { return /MAPK|RAF|ERK/i.test(r.name); });
ok("Reactome: >=1 significant term at FDR<0.05", sig.length > 0, sig.length + " significant");
ok("Reactome: a RAF/MAPK/ERK term is significant", mapk.length > 0,
   mapk.slice(0, 3).map(function (r) { return r.name; }).join("; "));

// 6) GO checks (only when the GO data is present; guarded so this stays green
// if a build omits GO).
if (fs.existsSync(path.join(ROOT, "data/human/go_bp.json"))) {
  var goBp = L("data/human/go_bp.json");
  var goUni = universe(goBp);
  ok("human GO-BP has thousands of terms", goBp.terms.length > 2000, "terms=" + goBp.terms.length);
  ok("human GO-BP universe sane (>5000 genes)", goUni.size > 5000, "size=" + goUni.size);
  var gs = goBp.terms.map(function (t) { return t.genes.length; });
  ok("GO-BP terms within size bounds [5,500]",
     Math.min.apply(null, gs) >= 5 && Math.max.apply(null, gs) <= 500,
     "min=" + Math.min.apply(null, gs) + " max=" + Math.max.apply(null, gs));
  var goMeta = man.collections.find(function (c) {
    return c.species === "human" && c.key === "go_bp" && c.available;
  });
  ok("manifest go_bp available", !!goMeta);
  // namespace tag is BP only in the BP file
  ok("GO-BP file holds only BP terms",
     goBp.terms.every(function (t) { return t.namespace === "GO-BP"; }));
  // manifest carries a GO release date for the footer
  ok("manifest GO release populated", !!(man.sources.go && man.sources.go.release),
     String(man.sources.go && man.sources.go.release));
}

console.log(fails === 0 ? "\nALL DATA TESTS PASSED" : "\n" + fails + " TEST(S) FAILED");
process.exit(fails === 0 ? 0 : 1);
