// test_stats.js - sanity checks for the ORA numerics. Run: node tests/test_stats.js
var S = require("../js/stats.js");

var fails = 0;
function approx(name, got, want, tol) {
  tol = tol || 1e-6;
  var ok = Math.abs(got - want) <= tol * Math.max(1, Math.abs(want));
  console.log((ok ? "PASS " : "FAIL ") + name + " got=" + got + " want=" + want);
  if (!ok) fails++;
}
function eq(name, got, want) {
  var ok = got === want;
  console.log((ok ? "PASS " : "FAIL ") + name + " got=" + got + " want=" + want);
  if (!ok) fails++;
}

// lnGamma against known factorials: lnGamma(n+1) = ln(n!)
approx("lnGamma(6)=ln(120)", S.lnGamma(6), Math.log(120), 1e-9);
approx("lnGamma(11)=ln(10!)", S.lnGamma(11), Math.log(3628800), 1e-9);

// Independent hypergeometric tail via exact binomials for a small case.
// N=50, K=10, n=8, k>=3. Compute reference with exact choose.
function chooseExact(n, k) {
  if (k < 0 || k > n) return 0;
  k = Math.min(k, n - k);
  var num = 1, den = 1;
  for (var i = 0; i < k; i++) { num *= (n - i); den *= (i + 1); }
  return num / den;
}
function hyperRef(k, n, K, N) {
  var kMax = Math.min(n, K);
  var denom = chooseExact(N, n);
  var s = 0;
  for (var i = k; i <= kMax; i++) s += chooseExact(K, i) * chooseExact(N - K, n - i) / denom;
  return s;
}
var lnFact = S.makeLnFactorial(60);
for (var k = 1; k <= 5; k++) {
  approx("hyper tail N50 K10 n8 k>=" + k,
    S.hyperUpperTail(k, 8, 10, 50, lnFact), hyperRef(k, 8, 10, 50), 1e-9);
}

// Large-N stability: must stay in [0,1] and be non-increasing in k.
var bigFact = S.makeLnFactorial(20001);
var prev = 1.0;
for (var kk = 1; kk <= 12; kk++) {
  var p = S.hyperUpperTail(kk, 200, 100, 20000, bigFact);
  if (p < 0 || p > 1) { console.log("FAIL bounds k=" + kk + " p=" + p); fails++; }
  if (p > prev + 1e-12) { console.log("FAIL monotonic k=" + kk); fails++; }
  prev = p;
}
console.log("PASS large-N bounds+monotonic (k=1..12), p(k=10)=" +
  S.hyperUpperTail(10, 200, 100, 20000, bigFact).toExponential(3));

// Fold enrichment
approx("fold", S.foldEnrichment(10, 200, 100, 20000), (10 / 200) / (100 / 20000), 1e-12);

// BH-FDR: known small example. p=[0.01,0.02,0.03,0.04,0.05], m=5
// adjusted (step-up): [0.05,0.05,0.05,0.05,0.05]
var q = S.bhFDR([0.01, 0.02, 0.03, 0.04, 0.05]);
approx("bh[0]", q[0], 0.05, 1e-9);
approx("bh[4]", q[4], 0.05, 1e-9);
// Bonferroni cap at 1
var b = S.bonferroni([0.3, 0.5]);
eq("bonf cap", b[0], 0.6);
eq("bonf cap2", b[1], 1);

console.log(fails === 0 ? "\nALL TESTS PASSED" : "\n" + fails + " TEST(S) FAILED");
process.exit(fails === 0 ? 0 : 1);
