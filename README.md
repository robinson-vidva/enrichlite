# enrichlite

Client-side gene-set over-representation analysis (ORA). Everything runs in the
browser: paste a gene list, pick a species and collection, and get a
hypergeometric enrichment table with BH-FDR and Bonferroni correction. No
server and no live API calls. Data is prepared at build time in Python and
shipped as static JSON.

Human and mouse are separate universes with separate symbol tables and gene-set
files. Matching is case-insensitive within a species only.

## Run locally

```
python3 -m http.server 8000
```

Open http://localhost:8000/ . The repo ships tiny DEMO gene sets so the page
works before you build real data; a banner marks demo mode.

## Build real data

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python3 scripts/build_genesets.py
```

Stage 1 builds Hallmark and Reactome for human and mouse. If the MSigDB mirror
needs auth, download the Hallmark GMTs manually and drop them in `data/raw/`
(e.g. `h.all.v2025.1.Hs.symbols.gmt`, `mh.all.v2025.1.Mm.symbols.gmt`), then
re-run. GO collections (BP/MF/CC) are added in Stage 2.

Pass `--cache-raw` to save every downloaded source into `data/raw/` and reuse
it on later runs instead of re-downloading (useful for Stage 2 GO work). The
build prints whether each source was downloaded or served from cache. To force
fresh downloads, clear the cache with `rm -rf data/raw/`.

## GitHub Pages

Hosted from the repository root of the `main` branch. In the GitHub repo:
Settings > Pages > Build and deployment > Source: "Deploy from a branch",
Branch: `main`, folder: `/ (root)`. Live URL:
https://robinson-vidva.github.io/enrichlite/

## Licensing

Code: MIT (see `LICENSE`). Data: see `DATA_LICENSES.md` (Reactome CC0, MSigDB
Hallmark CC BY 4.0, GO CC BY 4.0, NCBI public domain).
