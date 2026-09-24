# mega27-08-phage-design

Phage-host interaction prediction and falsifiable phage-cocktail design for
antibiotic-resistant *E. coli* and *Klebsiella pneumoniae*, from real public
data (790 phage genomes, NCBI/UniProt/AlphaFold/PDB). CNN+GNN interaction
model, receptor-level discovery screen, exact escape analysis, and a held-out-
genus generalization benchmark.

## The phagedesign tool

Works fully offline from committed artifacts (trained checkpoint, receptor
panel, normalization stats). No network, no training data needed.

```bash
python3 scripts/phagedesign.py predict --genome path/to/phage.gb
python3 scripts/phagedesign.py screen  --genome path/to/phage.gb
python3 scripts/phagedesign.py cocktail [--design PD8-EC1]
```

- `predict` - per-species host-range logits, probabilities, and ranking from a
  GenBank genome's annotated receptor-binding proteins.
- `screen` - per-receptor scores (which host receptor the phage likely uses).
- `cocktail` - the committed cocktail designs (members, targeted receptors,
  escape-probability bounds).

## Reproduce

```bash
PYTHONPATH=src python3 -m pytest tests/ -q   # 48 hermetic tests, no network
```

All result numbers cited in paper/paper.md live in results/*.json and
results/*.csv (see the artifact map, section 33 of the paper). The 55+ page
PDF is paper/paper.pdf; figures in paper/figures/.

## Data

Corpus: data/processed/phage_rbp_table.csv (790 phages, two host species) +
data/processed/salmonella_rbp_table.csv (84 Salmonella enterica phages,
held-out genus). Panels: host_receptor_panel.csv, salmonella_receptor_panel.csv.
Raw GenBank caches: data/raw/genbank/, data/raw/genbank_salmonella/.
