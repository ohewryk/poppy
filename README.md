# POPPy — Phyto-Ontology Platform for Pharmacology

POPPy is a reproducible phytotherapy ontology / knowledge graph and a static website to
explore it. It links medicinal plants to their natural-product compounds, drug targets, and
clinical evidence — enriched with names and cross-references from PubChem, ChEBI, UniChem,
and COCONUT — and surfaces multi-omic data (gene / protein / UniProt, NCBI Taxonomy) live,
on demand, from UniProt/NCBI.

## Explore it
- Live site: https://romanolab.github.io/poppy/
- Full ontology (RDF/XML, ~2 GB): https://upenn.box.com/v/poppyontology

## Repository
- notebooks/Ontology_Work_clean.ipynb — the ontology build pipeline (authoritative build-of-record)
- website/ — the static site; website/data/ is the sharded JSON the Explore page loads
- data/ontology/poppystructure.rdf — hand-curated TBox scaffold
- data/SOURCES.md — input datasets, versions, and where to obtain them
- src/poppy/, scripts/, configs/ — secondary config-driven pipeline (the notebook is authoritative)

## Build and deploy
- Assemble / rebuild: see ASSEMBLY.md (repo map, build order, data-layer generation)
- Publish the site: see BRINGING_ONLINE.md (merge, then enable GitHub Pages)

## How it works
The site is fully static (no build step), published to GitHub Pages by
.github/workflows/pages.yml. Explore searches all ~59,700 organisms and lazy-loads each
plant's compounds and targets from website/data/. Gene/protein/UniProt and NCBI Taxonomy
details are fetched live from UniProt/NCBI when a target or plant is opened — no NCBI data is
stored in the ontology.

## Run locally

    cd website && python3 -m http.server 8000
    # then open http://localhost:8000/Home.html

## License
Code: MIT. Ontology data: per-source licenses listed in data/SOURCES.md.
Developed in the Romano Lab at the University of Pennsylvania.
