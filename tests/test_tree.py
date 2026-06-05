#!/usr/bin/env python3
# test_tree.py - GO tree-view reduced-hierarchy file invariants.
# Run: python3 tests/test_tree.py
# File-only checks always run when tree files exist. The ontology-backed checks
# (exact reduction + skipE==0 + roots-are-root) run only when go-basic.obo and
# goatools are present (they are, in a full local build environment).
import glob
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import build_genesets as B

ROOT = os.path.join(os.path.dirname(__file__), "..")
fails = 0


def check(name, cond, extra=""):
    global fails
    print(("PASS " if cond else "FAIL ") + name + ((" (" + extra + ")") if extra else ""))
    if not cond:
        fails += 1


def has_cycle(parents):
    # parents[i] = list of parent indices. DFS three-color cycle detection.
    color = [0] * len(parents)  # 0 white, 1 gray, 2 black
    for s in range(len(parents)):
        if color[s]:
            continue
        stack = [(s, 0)]
        while stack:
            node, k = stack.pop()
            if k == 0:
                color[node] = 1
            if k < len(parents[node]):
                stack.append((node, k + 1))
                nxt = parents[node][k]
                if color[nxt] == 1:
                    return True
                if color[nxt] == 0:
                    stack.append((nxt, 0))
            else:
                color[node] = 2
    return False


tree_files = sorted(glob.glob(os.path.join(ROOT, "data", "*", "go_*_tree.json")))
if not tree_files:
    print("SKIP: no GO tree files present")
    sys.exit(0)

cases = []  # (label, terms, parents)
for tf in tree_files:
    parents = json.load(open(tf))["parents"]
    coll_path = tf[:-len("_tree.json")] + ".json"
    terms = json.load(open(coll_path))["terms"]
    label = os.path.relpath(tf, ROOT)
    cases.append((label, terms, parents))

    n = len(terms)
    check(label + ": parents length matches terms", len(parents) == n,
          "%d vs %d" % (len(parents), n))
    valid = all(all(0 <= p < n for p in row) for row in parents)
    check(label + ": all parent indices in range", valid)
    no_self = all(i not in parents[i] for i in range(len(parents)))
    check(label + ": no self-parent", no_self)
    sorted_uniq = all(row == sorted(set(row)) for row in parents)
    check(label + ": parent lists sorted and unique", sorted_uniq)
    check(label + ": acyclic", not has_cycle(parents))
    roots = sum(1 for row in parents if not row)
    check(label + ": has at least one root", roots > 0, "roots=%d" % roots)

# Ontology-backed exact check: re-derive the reduction and compare. build_tree
# also asserts skipE==0 internally, so a passing recompute proves the build-time
# honesty guarantee and that empty rows are genuinely root (no shipped ancestor).
obo = os.path.join(ROOT, "data", "raw", "go-basic.obo")
try:
    from goatools.obo_parser import GODag
    have_obo = os.path.exists(obo)
except Exception:
    have_obo = False

if have_obo:
    godag = GODag(obo, optional_attrs={"relationship"}, prt=None)
    pm, ns = B.build_parent_map(godag)
    for label, terms, parents in cases:
        try:
            recomputed = B.build_tree(terms, pm, ns)  # asserts skipE==0
            check(label + ": reduction matches build_tree (skipE==0, roots genuine)",
                  recomputed == parents)
        except AssertionError as e:
            check(label + ": reduction matches build_tree (skipE==0, roots genuine)",
                  False, str(e))
        # explicit roots-are-root spot check via ancestors()
        ids = [t["id"] for t in terms]
        shipped = set(ids)
        cache = {}
        bad_root = next((ids[i] for i, row in enumerate(parents)
                         if not row and any(a in shipped
                         for a in B.ancestors(ids[i], pm, ns, cache))), None)
        check(label + ": every root has no shipped ancestor", bad_root is None,
              "offender=" + str(bad_root) if bad_root else "")
else:
    print("SKIP: ontology checks (go-basic.obo / goatools not available)")

print("\n" + ("ALL TREE TESTS PASSED" if fails == 0 else str(fails) + " TEST(S) FAILED"))
sys.exit(0 if fails == 0 else 1)
