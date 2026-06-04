#!/usr/bin/env python3
# build_genesets.py - download sources and emit compact index-based JSON.
# Stage 1: Hallmark + Reactome for human and mouse. GO is Stage 2 (stubbed).

import argparse
import gzip
import io
import json
import os
import re
import sys
import zipfile

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RAW = os.path.join(DATA, "raw")

REACTOME_GMT_ZIP = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
REACTOME_PATHWAYS = "https://reactome.org/download/current/ReactomePathways.txt"
REACTOME_NCBI = "https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt"
GENEINFO = {
    "human": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    "mouse": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
}
MSIGDB_BASE = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"
SPECIES_TAX = {"human": "9606", "mouse": "10090"}
REACTOME_PREFIX = {"human": "R-HSA-", "mouse": "R-MMU-"}


def log(msg):
    print("[build] " + msg, flush=True)


def fetch_bytes(url, timeout=120):
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.content


def fetch_text(url, timeout=120):
    return fetch_bytes(url, timeout).decode("utf-8", "replace")


# --- gene_info: symbol table, protein-coding count, alias map ---

def load_gene_info(species):
    log("downloading gene_info for " + species)
    raw = fetch_bytes(GENEINFO[species])
    txt = gzip.decompress(raw).decode("utf-8", "replace")
    geneid_to_symbol = {}     # NCBI GeneID -> authoritative symbol
    coding_symbols = []       # protein-coding authoritative symbols
    synonyms = {}             # synonym -> symbol (protein-coding only)
    tax = SPECIES_TAX[species]
    for line in txt.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 12 or f[0] != tax:
            continue
        gid = f[1]
        sym = f[2]
        syn = f[4]
        type_of_gene = f[9]
        nom_sym = f[10]
        symbol = nom_sym if nom_sym and nom_sym != "-" else sym
        geneid_to_symbol[gid] = symbol
        if type_of_gene == "protein-coding":
            coding_symbols.append(symbol)
            if syn and syn != "-":
                for s in syn.split("|"):
                    if s and s != "-":
                        synonyms.setdefault(s, symbol)
    return geneid_to_symbol, coding_symbols, synonyms


# --- Reactome ---

def load_reactome_names():
    # stable id -> name, and name -> stable id per species prefix
    txt = fetch_text(REACTOME_PATHWAYS)
    id_to_name = {}
    name_to_id = {}
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) < 3:
            continue
        sid, name = f[0], f[1]
        id_to_name[sid] = name
        name_to_id.setdefault(name, sid)
    return id_to_name, name_to_id


def reactome_human_from_gmt(name_to_id):
    log("downloading Reactome human GMT")
    raw = fetch_bytes(REACTOME_GMT_ZIP)
    zf = zipfile.ZipFile(io.BytesIO(raw))
    gmt_name = [n for n in zf.namelist() if n.endswith(".gmt")][0]
    text = zf.read(gmt_name).decode("utf-8", "replace")
    terms = []
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        name = f[0]
        genes = [g for g in f[2:] if g]
        sid = name_to_id.get(name, name)
        terms.append({"id": sid, "name": name, "namespace": "Reactome Pathway", "symbols": genes})
    return terms


def reactome_mouse_from_ncbi(geneid_to_symbol, id_to_name):
    log("downloading NCBI2Reactome mapping (mouse)")
    txt = fetch_text(REACTOME_NCBI)
    by_path = {}
    for line in txt.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        gid, sid, name = f[0], f[1], f[3]
        if not sid.startswith(REACTOME_PREFIX["mouse"]):
            continue
        sym = geneid_to_symbol.get(gid)
        if not sym:
            continue
        entry = by_path.setdefault(sid, {"name": id_to_name.get(sid, name), "symbols": set()})
        entry["symbols"].add(sym)
    terms = []
    for sid, e in by_path.items():
        terms.append({"id": sid, "name": e["name"], "namespace": "Reactome Pathway",
                      "symbols": sorted(e["symbols"])})
    return terms


# --- MSigDB Hallmark ---

def msigdb_detect_release(session_get):
    # Probe the release dir listing for the newest x.y.Hs directory.
    try:
        html = session_get(MSIGDB_BASE + "/")
    except Exception as e:
        log("could not list MSigDB releases: " + str(e))
        return None
    rels = sorted(set(re.findall(r"(\d+\.\d+)\.Hs/", html)),
                  key=lambda v: [int(x) for x in v.split(".")])
    return rels[-1] if rels else None


def looks_like_gmt(text):
    first = text.lstrip().split("\n", 1)[0]
    return "\t" in first and "<html" not in text[:200].lower()


def load_hallmark(species, raw_dir):
    # Returns list of {id,name,namespace,symbols} or None if unavailable.
    fname = {"human": "h.all", "mouse": "mh.all"}[species]
    hs = {"human": "Hs", "mouse": "Mm"}[species]

    # 1) manual fallback: any matching file dropped in data/raw/
    if os.path.isdir(raw_dir):
        for n in sorted(os.listdir(raw_dir)):
            if n.startswith(fname + ".") and n.endswith("." + hs + ".symbols.gmt"):
                log("using manual Hallmark file data/raw/" + n)
                with open(os.path.join(raw_dir, n), "r", encoding="utf-8") as fh:
                    return parse_gmt(fh.read()), n
        # also accept exact simple name
        simple = os.path.join(raw_dir, fname + "." + hs + ".symbols.gmt")
        if os.path.isfile(simple):
            with open(simple, "r", encoding="utf-8") as fh:
                return parse_gmt(fh.read()), os.path.basename(simple)

    # 2) mirror auto-detect
    rel = msigdb_detect_release(fetch_text)
    if not rel:
        return None, None
    url = MSIGDB_BASE + "/" + rel + "." + hs + "/" + fname + ".v" + rel + "." + hs + ".symbols.gmt"
    log("trying MSigDB mirror " + url)
    try:
        text = fetch_text(url)
    except Exception as e:
        log("MSigDB download failed: " + str(e))
        return None, None
    if not looks_like_gmt(text):
        log("MSigDB response was not a GMT (likely auth wall)")
        return None, None
    return parse_gmt(text), os.path.basename(url)


def parse_gmt(text):
    terms = []
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        terms.append({"id": f[0], "name": f[0], "namespace": "hallmark",
                      "symbols": [g for g in f[2:] if g]})
    return terms


# --- assembly ---

def index_terms(terms, sym_to_idx, symbols, casefold):
    out = []
    for t in terms:
        idxs = []
        for s in t["symbols"]:
            key = s.upper() if casefold else s
            gi = sym_to_idx.get(key)
            if gi is None:
                gi = len(symbols)
                symbols.append(s)
                sym_to_idx[key] = gi
            idxs.append(gi)
        out.append({"id": t["id"], "name": t["name"], "namespace": t["namespace"],
                    "genes": sorted(set(idxs))})
    return out


def write_json(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, separators=(",", ":"))
    log("wrote " + os.path.relpath(path, ROOT))


def build_species(species, want, raw_dir, sources):
    log("=== building " + species + " ===")
    geneid_to_symbol, coding_symbols, synonyms = load_gene_info(species)

    symbols = []
    sym_to_idx = {}
    for s in coding_symbols:
        if s not in sym_to_idx:
            sym_to_idx[s] = len(symbols)
            symbols.append(s)
    coding_n = len(symbols)

    id_to_name, name_to_id = load_reactome_names()

    collections = {}
    if "reactome" in want:
        if species == "human":
            rt = reactome_human_from_gmt(name_to_id)
        else:
            rt = reactome_mouse_from_ncbi(geneid_to_symbol, id_to_name)
        collections["reactome"] = index_terms(rt, sym_to_idx, symbols, casefold=False)

    if "hallmark" in want:
        ht, src = load_hallmark(species, raw_dir)
        if ht is None:
            log("HALLMARK UNAVAILABLE for " + species + ".")
            log("Mirror appears to need auth. Download the GMT manually and drop it in data/raw/,")
            log("e.g. data/raw/h.all.vXXXX.Hs.symbols.gmt (human) or mh.all.vXXXX.Mm.symbols.gmt (mouse),")
            log("then re-run. Stopping so nothing is guessed.")
            sys.exit(2)
        collections["hallmark"] = index_terms(ht, sym_to_idx, symbols, casefold=False)
        sources["msigdb"]["version"] = src

    # alias map limited to symbols present in the table
    aliases = {}
    for syn, sym in synonyms.items():
        gi = sym_to_idx.get(sym)
        if gi is not None and syn not in sym_to_idx:
            aliases[syn] = gi

    write_json(os.path.join(DATA, species, "symbols.json"), {
        "species": species, "codingN": coding_n, "symbols": symbols, "aliases": aliases
    })
    for key, terms in collections.items():
        write_json(os.path.join(DATA, species, key + ".json"), {
            "collection": key, "species": species, "terms": terms
        })
    return collections, coding_n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", nargs="+", default=["human", "mouse"])
    ap.add_argument("--collections", nargs="+", default=["hallmark", "reactome"])
    ap.add_argument("--raw", default=RAW)
    ap.add_argument("--min", type=int, default=5, help="GO min term size (Stage 2)")
    ap.add_argument("--max", type=int, default=500, help="GO max term size (Stage 2)")
    ap.add_argument("--with-iea", action="store_true", help="GO include IEA (Stage 2)")
    ap.add_argument("--go", action="store_true", help="build GO collections (Stage 2, not yet implemented)")
    args = ap.parse_args()

    if args.go:
        log("GO build is Stage 2 and not yet implemented. Exiting.")
        sys.exit(2)

    sources = {"reactome": {"version": "current"}, "msigdb": {"version": None},
               "go": {"release": None, "doi": None}}

    manifest_collections = []
    symbols_paths = {}
    for sp in args.species:
        cols, coding_n = build_species(sp, set(args.collections), args.raw, sources)
        symbols_paths[sp] = "data/" + sp + "/symbols.json"
        labels = {"hallmark": "Hallmark", "reactome": "Reactome"}
        for key, terms in cols.items():
            uni = set()
            for t in terms:
                uni.update(t["genes"])
            manifest_collections.append({
                "key": key, "label": labels.get(key, key), "species": sp,
                "path": "data/" + sp + "/" + key + ".json", "available": True, "N": len(uni)
            })
    # keep GO placeholders visible but unavailable
    for sp in args.species:
        for key, label in [("go_bp", "GO-BP"), ("go_mf", "GO-MF"), ("go_cc", "GO-CC")]:
            manifest_collections.append({
                "key": key, "label": label, "species": sp,
                "path": "data/" + sp + "/" + key + ".json", "available": False, "N": 0
            })

    manifest = {"demo": False, "symbols": symbols_paths,
                "collections": manifest_collections, "sources": sources}
    write_json(os.path.join(DATA, "manifest.json"), manifest)
    log("done.")


if __name__ == "__main__":
    main()
