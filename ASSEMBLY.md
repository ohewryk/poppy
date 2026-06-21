# POPPy — Assembly & Build Guide

How the repository fits together and how to build/deploy the website. There are **two
paths**: (A) just publish the site that's already in the repo, and (B) rebuild the
ontology + data from sources. Most of the time you want **A**.

---

## Repository map

```
poppy/
├── notebooks/Ontology_Work_clean.ipynb   # the ontology build pipeline (build-of-record)
├── data/
│   ├── ontology/poppystructure.rdf        # hand-curated TBox scaffold (Protégé)
│   └── SOURCES.md                          # every input dataset + version + where to get it
├── website/                                # the static site (what gets published)
│   ├── Home.html  Explore.html  Download.html  … index.html
│   ├── ontology-data.js                    # demo graph (fallback)
│   ├── poppy-ontology-real.js              # curated 24-species subset (real ontology data)
│   ├── assets/                             # logo, plant cut-outs + plates
│   └── data/                               # ← the browse backend (sharded JSON)
│       ├── plants_index.json               #   all ~59.7k organisms (search)
│       ├── plant_edges/0..255.json         #   plant → compound links (lazy-loaded)
│       ├── compounds/0..255.json           #   compound details (lazy-loaded)
│       ├── targets.json                    #   target id → name
│       ├── compound_targets.json           #   compound → target links
│       └── meta.json                       #   counts
├── .github/workflows/pages.yml             # GitHub Pages deploy (publishes website/)
├── src/poppy/  scripts/  configs/          # secondary config-driven pipeline (NOT the build-of-record)
└── ASSEMBLY.md  BRINGING_ONLINE.md         # this guide + deploy notes
```

The **full ~2 GB enriched ontology RDF is NOT in the repo** — it's on Box:
<https://upenn.box.com/v/poppyontology> (linked from the Download page). Too large for Git.

---

## Path A — Publish the existing site (the usual case)

The site and its data layer are already committed, so there is **no build step** — just deploy.

1. Merge the open PR into `RomanoLab/poppy:main` (or work on `main` directly).
2. **Settings → Pages → Build and deployment → Source: GitHub Actions.**
3. On the next push to `main`, `.github/workflows/pages.yml` publishes `website/`. The live
   URL appears at the top of Settings → Pages (`https://romanolab.github.io/poppy/`).
4. Verify: open the URL → **Explore** → search a plant (e.g. *Panax ginseng*) → confirm
   compounds load; open a compound → click a **Target** → a live UniProt/NCBI panel appears.

That's it. Nothing is compiled; the HTML/JS is served as-is and reads the JSON in `website/data/`.

---

## Path B — Rebuild the ontology + data from sources

Only needed when the underlying data changes. Run on a machine with ~32 GB RAM (the graph is large).

**1. Environment**
```bash
git clone https://github.com/RomanoLab/poppy.git && cd poppy
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,chem]"     # rdflib, rdkit, pandas, requests, psycopg2, biopython, tqdm
export ROMANO_DB_PASSWORD=…       # PMACS opendata DB (lab access)
export ENTREZ_EMAIL=you@upenn.edu
```

**2. Inputs** — download per `data/SOURCES.md` (COCONUT, CMAUP, DrugCentral, ChEMBL, Dr. Duke's,
ChEBI, COCONUT 2.0 SDF). `data/ontology/poppystructure.rdf` is already in the repo.

**3. Build the ontology** — open `notebooks/Ontology_Work_clean.ipynb` and run top to bottom:
- Stages 1–13: base graph → COCONUT/CMAUP/DrugCentral/ChEMBL/Dr. Duke's → cleanup → ChEBI.
- Stage 14: external names (PubChem/ChEBI), cross-refs, NP classes.
- Stage 15: IUPAC names → `phytotherapies_named.rdf`.
- Stage 16: COCONUT 2.0 occurrences (SDF) → `phytotherapies_named_coconut.rdf`.
- Stage 17: final sanitize + dedupe.

The result, **`phytotherapies_named_coconut.rdf`**, is the canonical ontology → upload to Box.

**4. Generate the website data layer** from that graph (Appendix A) → writes `website/data/`.

**5. Commit & deploy**
```bash
git add website/data && git commit -m "Rebuild browse data" && git push origin main
```
Pages redeploys automatically. Upload the new RDF to Box and update the link if it changed.

---

## Live multi-omics (no build step)

Gene/protein/UniProt and plant taxonomy data are **not stored** in the ontology or the repo —
they're fetched **on demand from UniProt/NCBI** when a user opens a Target or Plant page
(`showTargetOmics` / `showPlantOmics` in `Explore.html`). Nothing to build or refresh.

---

## Notes & known limitations

- **Big files → Box, not Git.** The 2 GB RDF (and any bulk CSV exports) live on Box; the repo
  keeps only the compact sharded JSON the site needs (~50 MB).
- **Non-plant organisms.** The full index currently includes some bacteria/fungi/animals
  (e.g. *Streptomyces*, *Homo sapiens*) pulled in by COCONUT 2.0, which skipped the NCBI/POWO
  plant validation the base pipeline uses. A taxonomy filter is a planned follow-up.
- **Targets are sparse** (~3,200 compound→target links), so multi-omics appears mainly on
  drug-like compounds.
- `src/poppy/` + `scripts/` + `configs/` are a parallel config-driven pipeline kept for reuse;
  the **notebook is the authoritative build**.

---

## Appendix A — Data-layer generation (run after Stage 17, on the canonical `g`)

These cells turn the ontology into the sharded JSON the site reads. (They are also the cells
used to produce the committed `website/data/`.)

```python
import json, os, re
from collections import defaultdict
NS_ = NS  # phytotherapies namespace, defined in Setup
D = "website/data"; os.makedirs(f"{D}/plant_edges", exist_ok=True); os.makedirs(f"{D}/compounds", exist_ok=True)
IK = re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')
def djb2(s):
    h=5381
    for c in s: h=((h*33)^ord(c))&0xFFFFFFFF
    return h&255
def cname(u):                       # compound display name
    for _,_,o in g.triples((u,hasCommonName,None)):
        x=str(o).strip()
        if x and not IK.match(x): return x
    for _,_,o in g.triples((u,hasIUPACName,None)):
        x=str(o).strip()
        if x and x.lower()!="nan": return x
    for _,_,o in g.triples((u,hasMolecularFormula,None)): return str(o)
    return str(u).rsplit('#',1)[-1]
def pname(u):                       # plant display name = rdfs:label (scientific name)
    for _,_,o in g.triples((u,RDFS.label,None)):
        if str(o).strip(): return str(o).strip()
    return str(u).rsplit('#',1)[-1]
def n(u,p): return sum(1 for _ in g.triples((u,p,None)))

# compounds + plants + edges
compounds={}
for s,_,o in g.triples((None,hasInChIKey,None)):
    cid=str(s).rsplit('#',1)[-1]
    if cid in compounds: continue
    compounds[cid]={"name":cname(s),
        "formula":next((str(x) for _,_,x in g.triples((s,hasMolecularFormula,None))),""),
        "mw":next((str(x) for _,_,x in g.triples((s,hasMolecularWeight,None))),""),
        "inchikey":str(o).strip(),"trials":n(s,hasClinicalStudy)}
plants={}; shards=defaultdict(dict); ncomp=defaultdict(int)
for s in g.subjects(RDF.type, NS_.Plant):
    pid=str(s).rsplit('#',1)[-1]
    plants[pid]={"id":pid,"name":pname(s),"trials":n(s,hasClinicalStudy),"papers":n(s,hasPaper),"nc":0}
for s,_,o in g.triples((None,hasCompound,None)):
    pid=str(s).rsplit('#',1)[-1]; cid=str(o).rsplit('#',1)[-1]
    if pid in plants and cid in compounds:
        shards[djb2(pid)].setdefault(pid,[]).append(cid); ncomp[pid]+=1
for pid,c in ncomp.items():
    if pid in plants: plants[pid]["nc"]=c

# targets (names only — gene/UniProt fetched live on the site)
tP=NS_.targetsPathway; targets={}; comp_tgt={}
for c,_,t in g.triples((None,tP,None)):
    cid=str(c).rsplit('#',1)[-1]; tid=str(t).rsplit('#',1)[-1]
    comp_tgt.setdefault(cid,[]);
    if tid not in comp_tgt[cid]: comp_tgt[cid].append(tid)
    targets.setdefault(tid,{"name":tid.replace("TargetedPathway_","").replace("Pathway_","").replace("_"," ").strip()})

import json as J
for b in range(256):
    J.dump(shards.get(b,{}), open(f"{D}/plant_edges/{b}.json","w"), separators=(",",":"))
    J.dump({k:v for k,v in compounds.items() if djb2(k)==b}, open(f"{D}/compounds/{b}.json","w"), ensure_ascii=False, separators=(",",":"))
J.dump(sorted(plants.values(), key=lambda x:(-x["nc"],x["name"])), open(f"{D}/plants_index.json","w"), ensure_ascii=False, separators=(",",":"))
J.dump(targets, open(f"{D}/targets.json","w"), ensure_ascii=False, separators=(",",":"))
J.dump(comp_tgt, open(f"{D}/compound_targets.json","w"), ensure_ascii=False, separators=(",",":"))
J.dump({"plants":len(plants),"compounds":len(compounds),
        "plant_compound_links":sum(ncomp.values()),"shard":"djb2(id)&255"}, open(f"{D}/meta.json","w"), indent=1)
print("plants",len(plants),"compounds",len(compounds),"targets",len(targets))
```

Then `git add website/data && git commit && git push` — Pages redeploys.
