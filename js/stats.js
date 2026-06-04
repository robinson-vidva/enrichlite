// stats.js - pure numerics for ORA. Usable in a worker via importScripts
// and in Node-less unit tests by attaching to globalThis.
(function (root) {
  "use strict";

  // Lanczos approximation for ln(Gamma(z)), z > 0.
  var G = 7;
  var C = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
  ];

  function lnGamma(z) {
    if (z < 0.5) {
      // reflection formula
      return Math.log(Math.PI / Math.sin(Math.PI * z)) - lnGamma(1 - z);
    }
    z -= 1;
    var x = C[0];
    for (var i = 1; i < G + 2; i++) x += C[i] / (z + i);
    var t = z + G + 0.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  }

  // Precomputed ln(n!) table to avoid repeated lnGamma calls.
  function makeLnFactorial(nMax) {
    var tbl = new Float64Array(nMax + 1);
    tbl[0] = 0;
    for (var i = 1; i <= nMax; i++) tbl[i] = tbl[i - 1] + Math.log(i);
    return tbl;
  }

  // lnChoose using a factorial table when in range, else lnGamma.
  function lnChoose(n, k, lnFact) {
    if (k < 0 || k > n) return -Infinity;
    if (lnFact && n < lnFact.length) {
      return lnFact[n] - lnFact[k] - lnFact[n - k];
    }
    return lnGamma(n + 1) - lnGamma(k + 1) - lnGamma(n - k + 1);
  }

  // Hypergeometric upper tail P(X >= k).
  // N background, K term genes in universe, n query genes in universe, k overlap.
  // Summed via pmf recurrence in log space for stability at large N.
  function hyperUpperTail(k, n, K, N, lnFact) {
    if (k <= 0) return 1.0;
    var kMax = Math.min(n, K);
    if (k > kMax) return 0.0;
    var lnDenom = lnChoose(N, n, lnFact);
    // pmf(i) = C(K,i) C(N-K,n-i) / C(N,n)
    var lnP = lnChoose(K, k, lnFact) + lnChoose(N - K, n - k, lnFact) - lnDenom;
    var p = Math.exp(lnP);
    var sum = p;
    for (var i = k; i < kMax; i++) {
      // ratio pmf(i+1)/pmf(i)
      var ratio = ((K - i) * (n - i)) / ((i + 1) * (N - K - n + i + 1));
      p *= ratio;
      sum += p;
    }
    if (sum > 1) sum = 1;
    if (sum < 0) sum = 0;
    return sum;
  }

  function foldEnrichment(k, n, K, N) {
    if (n === 0 || K === 0) return 0;
    return (k / n) / (K / N);
  }

  // Benjamini-Hochberg step-up. Returns adjusted q in original order.
  function bhFDR(pvals) {
    var m = pvals.length;
    var idx = pvals.map(function (p, i) { return i; });
    idx.sort(function (a, b) { return pvals[a] - pvals[b]; });
    var q = new Array(m);
    var prev = 1;
    for (var rank = m - 1; rank >= 0; rank--) {
      var i = idx[rank];
      var val = pvals[i] * m / (rank + 1);
      if (val > prev) val = prev;
      if (val > 1) val = 1;
      prev = val;
      q[i] = val;
    }
    return q;
  }

  function bonferroni(pvals) {
    var m = pvals.length;
    return pvals.map(function (p) {
      var v = p * m;
      return v > 1 ? 1 : v;
    });
  }

  var api = {
    lnGamma: lnGamma,
    makeLnFactorial: makeLnFactorial,
    lnChoose: lnChoose,
    hyperUpperTail: hyperUpperTail,
    foldEnrichment: foldEnrichment,
    bhFDR: bhFDR,
    bonferroni: bonferroni
  };

  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.EnrichStats = api;
})(typeof self !== "undefined" ? self : this);
