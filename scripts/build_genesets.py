#!/usr/bin/env python3
# build_genesets.py - download sources and emit compact index-based JSON.
# Hallmark + Reactome + GO (BP/MF/CC) for human and mouse.

import argparse
import gzip
import io
import json
import os
import re
import shutil
import sys
import tempfile
import zipfile
from datetime import datetime, timezone

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
GO_OBO = "https://current.geneontology.org/ontology/go-basic.obo"
GAF = {
    "human": "https://current.geneontology.org/annotations/goa_human.gaf.gz",
    "mouse": "https://current.geneontology.org/annotations/mgi.gaf.gz",
}
# GO concept DOI (all-versions, not release-specific). The OBO header carries
# no DOI, so this is intentionally set rather than parsed.
GO_CONCEPT_DOI = "10.5281/zenodo.1205166"
SPECIES_TAX = {"human": "9606", "mouse": "10090"}
REACTOME_PREFIX = {"human": "R-HSA-", "mouse": "R-MMU-"}
GO_NS = {
    "biological_process": ("go_bp", "GO-BP"),
    "molecular_function": ("go_mf", "GO-MF"),
    "cellular_component": ("go_cc", "GO-CC"),
}
LABELS = {
    "hallmark": "Hallmark", "reactome": "Reactome",
    "go_bp": "GO-BP", "go_mf": "GO-MF", "go_cc": "GO-CC",
    "go_bp_iea": "GO-BP (+IEA)", "go_mf_iea": "GO-MF (+IEA)", "go_cc_iea": "GO-CC (+IEA)",
}


def log(msg):
    print("[build] " + msg, flush=True)


def fetch_bytes(url, timeout=120):
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.content


def fetch_text(url, timeout=120):
    return fetch_bytes(url, timeout).decode("utf-8", "replace")


def human_size(nbytes):
    if nbytes >= 1024 * 1024:
        return "%.1f MB" % (nbytes / (1024.0 * 1024.0))
    if nbytes >= 1024:
        return "%.1f KB" % (nbytes / 1024.0)
    return str(nbytes) + " B"


def atomic_write(path, data, mode="wb"):
    # Write to a temp name and rename, so a failure never leaves a partial file.
    tmp = path + ".tmp"
    with open(tmp, mode) as fh:
        fh.write(data)
    os.replace(tmp, path)


class RawCache:
    # When enabled, save downloads under raw_dir and reuse them on later runs.
    def __init__(self, raw_dir, enabled):
        self.raw = raw_dir
        self.enabled = enabled
        self._tmp = None

    def get_bytes(self, url, fname, label):
        path = os.path.join(self.raw, fname)
        if self.enabled and os.path.isfile(path):
            with open(path, "rb") as fh:
                data = fh.read()
            log(label + ": cached (" + os.path.relpath(path, ROOT) + ")")
            return data
        data = fetch_bytes(url)  # raises before any write on failure
        log(label + ": downloaded " + human_size(len(data)))
        if self.enabled:
            os.makedirs(self.raw, exist_ok=True)
            atomic_write(path, data)
        return data

    def get_text(self, url, fname, label):
        return self.get_bytes(url, fname, label).decode("utf-8", "replace")

    # Returns a filesystem path (for tools like goatools that need a file).
    def get_path(self, url, fname, label):
        if self.enabled:
            os.makedirs(self.raw, exist_ok=True)
            path = os.path.join(self.raw, fname)
            if os.path.isfile(path):
                log(label + ": cached (" + os.path.relpath(path, ROOT) + ")")
                return path
            data = fetch_bytes(url)
            log(label + ": downloaded " + human_size(len(data)))
            atomic_write(path, data)
            return path
        if self._tmp is None:
            self._tmp = tempfile.mkdtemp(prefix="enrichlite_")
        data = fetch_bytes(url)
        log(label + ": downloaded " + human_size(len(data)))
        path = os.path.join(self._tmp, fname)
        atomic_write(path, data)
        return path

    def cleanup(self):
        if self._tmp:
            shutil.rmtree(self._tmp, ignore_errors=True)
            self._tmp = None


# --- gene_info: symbol table, protein-coding count, alias map ---

def load_gene_info(species, dl):
    raw = dl.get_bytes(GENEINFO[species], os.path.basename(GENEINFO[species]),
                       "gene_info " + species)
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

def load_reactome_names(dl):
    # stable id -> name, and name -> stable id per species prefix
    txt = dl.get_text(REACTOME_PATHWAYS, "ReactomePathways.txt", "Reactome pathway names")
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


def reactome_human_from_gmt(name_to_id, dl):
    raw = dl.get_bytes(REACTOME_GMT_ZIP, "ReactomePathways.gmt.zip", "Reactome human GMT")
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


def reactome_mouse_from_ncbi(geneid_to_symbol, id_to_name, dl):
    txt = dl.get_text(REACTOME_NCBI, "NCBI2Reactome_All_Levels.txt", "NCBI2Reactome mapping")
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


def parse_release(name):
    # Extract the MSigDB release (e.g. 2026.1) from a GMT filename.
    m = re.search(r"\.v(\d+\.\d+(?:\.\d+)?)\.", name)
    return m.group(1) if m else None


def load_hallmark(species, dl):
    # Returns (terms, filename, release); terms is None if unavailable.
    fname = {"human": "h.all", "mouse": "mh.all"}[species]
    hs = {"human": "Hs", "mouse": "Mm"}[species]
    label = "MSigDB Hallmark " + species

    # 1) cached or manually dropped file in data/raw/ (treated the same)
    raw_dir = dl.raw
    if os.path.isdir(raw_dir):
        for n in sorted(os.listdir(raw_dir)):
            if n.startswith(fname + ".") and n.endswith("." + hs + ".symbols.gmt"):
                log(label + ": cached (data/raw/" + n + ")")
                with open(os.path.join(raw_dir, n), "r", encoding="utf-8") as fh:
                    return parse_gmt(fh.read()), n, parse_release(n)

    # 2) mirror auto-detect + download (auth-wall sniff kept exactly as-is)
    rel = msigdb_detect_release(fetch_text)
    if not rel:
        return None, None, None
    url = MSIGDB_BASE + "/" + rel + "." + hs + "/" + fname + ".v" + rel + "." + hs + ".symbols.gmt"
    try:
        text = fetch_text(url)
    except Exception as e:
        log(label + ": download failed " + str(e))
        return None, None, None
    if not looks_like_gmt(text):
        log(label + ": response was not a GMT (likely auth wall)")
        return None, None, None
    log(label + ": downloaded " + human_size(len(text.encode("utf-8"))))
    name = os.path.basename(url)
    if dl.enabled:
        os.makedirs(raw_dir, exist_ok=True)
        atomic_write(os.path.join(raw_dir, name), text, mode="w")
    return parse_gmt(text), name, rel


def parse_gmt(text):
    terms = []
    for line in text.splitlines():
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        terms.append({"id": f[0], "name": f[0], "namespace": "hallmark",
                      "symbols": [g for g in f[2:] if g]})
    return terms


# --- Gene Ontology ---

def keep_annotation(qualifier, evidence, with_iea):
    # Drop ND always; drop NOT-qualified (pipe-delimited token, not substring);
    # drop IEA unless explicitly included.
    if evidence == "ND":
        return False
    if "NOT" in (qualifier.split("|") if qualifier else []):
        return False
    if not with_iea and evidence == "IEA":
        return False
    return True


def build_parent_map(godag):
    # parent_map: term -> set of direct parents over is_a + part_of.
    # namespaces: term -> namespace. Both keyed by primary GO id.
    parent_map = {}
    namespaces = {}
    for term in godag.values():
        pid = term.item_id
        if pid in parent_map:
            continue
        namespaces[pid] = term.namespace
        parents = set(p.item_id for p in term.parents)
        rel = getattr(term, "relationship", {}) or {}
        for p in rel.get("part_of", set()):
            parents.add(p.item_id)
        parent_map[pid] = parents
    return parent_map, namespaces


def ancestors(term, parent_map, namespaces, cache):
    # All ancestors reachable while staying within the term's own namespace.
    if term in cache:
        return cache[term]
    base = namespaces.get(term)
    res = set()
    for p in parent_map.get(term, ()):
        if namespaces.get(p) != base:
            continue
        res.add(p)
        res |= ancestors(p, parent_map, namespaces, cache)
    cache[term] = res
    return res


def propagate(annotations, parent_map, namespaces):
    # True-path rule: a gene annotated to a term counts for that term and all
    # of its (same-namespace) ancestors.
    cache = {}
    term2genes = {}
    for symbol, go_id in annotations:
        if go_id not in namespaces:
            continue
        targets = ancestors(go_id, parent_map, namespaces, cache) | {go_id}
        for t in targets:
            term2genes.setdefault(t, set()).add(symbol)
    return term2genes


def parse_gaf(path, godag, with_iea):
    # GAF 2.2: symbol col3, qualifier col4, GO ID col5, evidence col7, aspect col9.
    ann = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line[0] == "!":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            symbol, qualifier, go_raw, evidence = f[2], f[3], f[4], f[6]
            if not keep_annotation(qualifier, evidence, with_iea):
                continue
            term = godag.get(go_raw)
            if term is None:
                continue
            ann.append((symbol, term.item_id))
    return ann


def build_go_collections(godag, parent_map, namespaces, gaf_path, with_iea, mn, mx):
    ann = parse_gaf(gaf_path, godag, with_iea)
    term2genes = propagate(ann, parent_map, namespaces)
    out = {"go_bp": [], "go_mf": [], "go_cc": []}
    # Emit terms in sorted GO-id order so output (and the downstream symbol
    # append order) is deterministic regardless of set/dict iteration order.
    for go_id in sorted(term2genes):
        genes = term2genes[go_id]
        if len(genes) < mn or len(genes) > mx:
            continue
        ns = namespaces.get(go_id)
        if ns not in GO_NS:
            continue
        key = GO_NS[ns][0]
        term = godag.get(go_id)
        out[key].append({"id": go_id, "name": term.name, "namespace": GO_NS[ns][1],
                         "symbols": sorted(genes)})
    return out


def parse_obo_meta(obo_path):
    # Release date from data-version; Zenodo DOI only if present (never guessed).
    release = None
    doi = None
    with open(obo_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("[Term]"):
                break
            if line.startswith("data-version:"):
                v = line.split(":", 1)[1].strip()
                m = re.search(r"(\d{4}-\d{2}-\d{2})", v)
                release = m.group(1) if m else v
            m = re.search(r"(10\.5281/zenodo\.\d+)", line)
            if m:
                doi = m.group(1)
    return release, doi


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


def build_species(species, want, dl, sources, go_ctx=None):
    log("=== building " + species + " ===")
    geneid_to_symbol, coding_symbols, synonyms = load_gene_info(species, dl)

    symbols = []
    sym_to_idx = {}
    for s in coding_symbols:
        if s not in sym_to_idx:
            sym_to_idx[s] = len(symbols)
            symbols.append(s)
    coding_n = len(symbols)

    id_to_name, name_to_id = load_reactome_names(dl)

    collections = {}
    if "reactome" in want:
        if species == "human":
            rt = reactome_human_from_gmt(name_to_id, dl)
        else:
            rt = reactome_mouse_from_ncbi(geneid_to_symbol, id_to_name, dl)
        collections["reactome"] = index_terms(rt, sym_to_idx, symbols, casefold=False)

    if "hallmark" in want:
        ht, src_file, src_rel = load_hallmark(species, dl)
        if ht is None:
            log("HALLMARK UNAVAILABLE for " + species + ".")
            log("Mirror appears to need auth. Download the GMT manually and drop it in data/raw/,")
            log("e.g. data/raw/h.all.vXXXX.Hs.symbols.gmt (human) or mh.all.vXXXX.Mm.symbols.gmt (mouse),")
            log("then re-run. Stopping so nothing is guessed.")
            sys.exit(2)
        collections["hallmark"] = index_terms(ht, sym_to_idx, symbols, casefold=False)
        if src_rel:
            sources["msigdb"]["version"] = src_rel
        if src_file:
            sources["msigdb"]["files"].append(src_file)

    if go_ctx:
        gaf_path = dl.get_path(GAF[species], os.path.basename(GAF[species]), "GAF " + species)
        go_cols = build_go_collections(go_ctx["godag"], go_ctx["pm"], go_ctx["ns"], gaf_path,
                                       False, go_ctx["min"], go_ctx["max"])
        for key in ("go_bp", "go_mf", "go_cc"):
            collections[key] = index_terms(go_cols[key], sym_to_idx, symbols, casefold=False)
            log(species + " " + key + ": " + str(len(go_cols[key])) + " terms")
        if go_ctx["with_iea"]:
            iea = build_go_collections(go_ctx["godag"], go_ctx["pm"], go_ctx["ns"], gaf_path,
                                       True, go_ctx["min"], go_ctx["max"])
            for key in ("go_bp", "go_mf", "go_cc"):
                collections[key + "_iea"] = index_terms(iea[key], sym_to_idx, symbols, casefold=False)

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
    ap.add_argument("--cache-raw", action="store_true",
                    help="save downloads into data/raw/ and reuse them on later runs")
    ap.add_argument("--min", type=int, default=5, help="GO min term size")
    ap.add_argument("--max", type=int, default=500, help="GO max term size")
    ap.add_argument("--with-iea", action="store_true", help="also emit IEA-inclusive GO variants")
    ap.add_argument("--go", action="store_true", help="build GO collections (full rebuild)")
    args = ap.parse_args()

    sources = {"reactome": {"version": "current"}, "msigdb": {"version": None, "files": []},
               "go": {"release": None, "doi": None}}

    dl = RawCache(args.raw, args.cache_raw)

    go_ctx = None
    if args.go:
        from goatools.obo_parser import GODag
        obo_path = dl.get_path(GO_OBO, "go-basic.obo", "GO ontology (go-basic.obo)")
        godag = GODag(obo_path, optional_attrs={"relationship"}, prt=None)
        pm, ns = build_parent_map(godag)
        release, doi = parse_obo_meta(obo_path)
        sources["go"]["release"] = release
        # OBO has no DOI; set the GO concept DOI intentionally.
        sources["go"]["doi"] = doi or GO_CONCEPT_DOI
        log("GO release=" + str(release) + " doi=" + str(sources["go"]["doi"]))
        go_ctx = {"godag": godag, "pm": pm, "ns": ns,
                  "min": args.min, "max": args.max, "with_iea": args.with_iea}

    manifest_collections = []
    symbols_paths = {}
    coding_counts = {}
    built = {}
    for sp in args.species:
        cols, coding_n = build_species(sp, set(args.collections), dl, sources, go_ctx)
        built[sp] = cols
        coding_counts[sp] = coding_n
        symbols_paths[sp] = "data/" + sp + "/symbols.json"

        def universe_n(terms):
            u = set()
            for t in terms:
                u.update(t["genes"])
            return len(u)

        # Base (non-IEA) collections become dropdown entries; an IEA variant is
        # attached to its base as `iea` (path + N), not a separate dropdown entry.
        for key, terms in cols.items():
            if key.endswith("_iea"):
                continue
            entry = {
                "key": key, "label": LABELS.get(key, key), "species": sp,
                "path": "data/" + sp + "/" + key + ".json", "available": True,
                "N": universe_n(terms)
            }
            iea_key = key + "_iea"
            if iea_key in cols:
                entry["iea"] = {"path": "data/" + sp + "/" + iea_key + ".json",
                                "N": universe_n(cols[iea_key])}
            manifest_collections.append(entry)
    # GO placeholders for namespaces not built this run (keep dropdown stable)
    for sp in args.species:
        for key, label in [("go_bp", "GO-BP"), ("go_mf", "GO-MF"), ("go_cc", "GO-CC")]:
            if key not in built[sp]:
                manifest_collections.append({
                    "key": key, "label": label, "species": sp,
                    "path": "data/" + sp + "/" + key + ".json", "available": False, "N": 0
                })

    version = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest = {"demo": False, "version": version, "symbols": symbols_paths,
                "collections": manifest_collections, "sources": sources}
    write_json(os.path.join(DATA, "manifest.json"), manifest)
    dl.cleanup()
    for sp in args.species:
        log("codingN " + sp + "=" + str(coding_counts[sp]))
    log("done. version=" + version)


if __name__ == "__main__":
    main()
