# Data licenses

The MIT LICENSE covers the enrichlite source code only. The gene-set data under
`data/` is derived from third-party sources and is governed by their licenses:

## Reactome
- License: CC0 1.0 (public domain dedication).
- Source: https://reactome.org/ (download/current).
- Files derived: human and mouse Reactome pathway gene sets.

## MSigDB Hallmark
- License: CC BY 4.0.
- Attribution: Broad Institute, MIT, and the Regents of the University of California.
- Source: https://www.gsea-msigdb.org/gsea/msigdb/
- Files derived: `h.all` (human) and `mh.all` (mouse) Hallmark gene sets.

## Gene Ontology (Stage 2)
- License: CC BY 4.0.
- Source: https://geneontology.org/
- The app footer must display the GO release date and the Zenodo DOI of the
  release used.
- Files derived: GO-BP, GO-MF, GO-CC gene sets (added in Stage 2).

## NCBI Gene
- NCBI gene_info is used at build time to map NCBI GeneIDs to symbols, to count
  protein-coding genes, and to build alias maps. NCBI data are in the public
  domain (U.S. Government work); see https://www.ncbi.nlm.nih.gov/home/about/policies/.
