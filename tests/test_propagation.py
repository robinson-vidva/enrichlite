#!/usr/bin/env python3
# test_propagation.py - GO true-path propagation, namespace isolation, and
# GAF annotation filtering. Run: python3 tests/test_propagation.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import build_genesets as B

fails = 0


def check(name, cond):
    global fails
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails += 1


# True-path: A is_a B, B part_of C (all BP). A gene on A must reach B and C.
pm = {"A": {"B"}, "B": {"C"}, "C": set()}
ns = {"A": "biological_process", "B": "biological_process", "C": "biological_process"}
t2g = B.propagate([("g1", "A")], pm, ns)
check("true-path: gene under annotated term A", "g1" in t2g.get("A", set()))
check("true-path: gene under is_a ancestor B", "g1" in t2g.get("B", set()))
check("true-path: gene under part_of ancestor C", "g1" in t2g.get("C", set()))

# Namespace isolation: B (BP) has a parent C in a different namespace (MF).
# A gene annotated to A (BP) must never reach the MF term C.
ns2 = {"A": "biological_process", "B": "biological_process", "C": "molecular_function"}
t2 = B.propagate([("g2", "A")], pm, ns2)
check("ns isolation: gene under A (BP)", "g2" in t2.get("A", set()))
check("ns isolation: gene under B (BP)", "g2" in t2.get("B", set()))
check("ns isolation: gene NOT under C (MF)", "g2" not in t2.get("C", set()))

# MF/CC annotated genes must not surface under unrelated namespaces.
pm3 = {"bp1": set(), "mf1": set(), "cc1": set()}
ns3 = {"bp1": "biological_process", "mf1": "molecular_function", "cc1": "cellular_component"}
t3 = B.propagate([("mfg", "mf1"), ("ccg", "cc1")], pm3, ns3)
check("mf gene stays in MF, not BP", t3.get("mf1") == {"mfg"} and "mfg" not in t3.get("bp1", set()))
check("cc gene stays in CC, not BP", t3.get("cc1") == {"ccg"} and "ccg" not in t3.get("bp1", set()))

# Annotation filtering.
check("ND dropped", B.keep_annotation("involved_in", "ND", False) is False)
check("IEA dropped by default", B.keep_annotation("involved_in", "IEA", False) is False)
check("IEA kept with --with-iea", B.keep_annotation("involved_in", "IEA", True) is True)
check("normal annotation kept", B.keep_annotation("involved_in", "IDA", False) is True)
check("empty qualifier kept", B.keep_annotation("", "IDA", False) is True)
# NOT must be a pipe-delimited token, not a substring.
check("NOT token dropped", B.keep_annotation("NOT|involved_in", "IDA", False) is False)
check("involved_in kept", B.keep_annotation("involved_in", "IDA", False) is True)
check("token-based NOT (NOTE|involved_in kept)",
      B.keep_annotation("NOTE|involved_in", "IDA", False) is True)

print("\n" + ("ALL PROPAGATION TESTS PASSED" if fails == 0 else str(fails) + " TEST(S) FAILED"))
sys.exit(0 if fails == 0 else 1)
