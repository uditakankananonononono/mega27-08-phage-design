# Structure-Aware Prediction of Bacteriophage Host Range from Receptor-Binding Proteins: Benchmarks, a Screening Tool, and Falsifiable Host-Switch Candidates

**Author:** Udita Phookan (with computational assistance)
**Date:** September 24, 2026
**Project:** MEGA-27, item 8 — phage design for antibiotic-resistant bacteria

## Abstract

Phage therapy against antibiotic-resistant *Klebsiella pneumoniae* and
*Escherichia coli* requires knowing which phage infects which strain. Existing
computational host predictors (WIsH, PHP, DeepHost, iPHoP-class) score whole
genome k-mer similarity and never model the actual molecular event of
infection: binding of the phage receptor-binding protein (RBP; tail fiber or
tail spike) to a host surface receptor. We build a structure-aware
host-range predictor from real public data only: 1,439 RBP sequences mined
from 922 complete phage genomes (NCBI GenBank), a curated 14-receptor panel
for *E. coli* K-12 and *K. pneumoniae* (UniProt), and all 14 receptor
structures (AlphaFold DB). A hybrid model — a 23-channel 1D CNN over RBP
residue profiles coupled to host receptors through a bilinear score plus a
k-mer interaction head — reaches test AUROC 0.973 (AUPRC 0.974, accuracy
0.930, per-phage choice accuracy 0.924) under group splitting by phage,
beating a strong Hadamard k-mer logistic baseline (0.910 / 0.867). A pure
CNN model without compositional features scores 0.609, an honest negative we
analyze. Screening all 790 phages yields 17 named candidate host-switch
phages (11 held-out), each a falsifiable plaque-assay prediction with a named
predicted receptor (e.g. OmpK36-mediated Klebsiella tropism for five
E. coli-annotated phages at p>0.99). A rigid grid-search docker validated
against the experimental T5 pb5–FhuA complex (PDB 8A8C) fails to recover the
interface (iRMSD 40 A), honestly reported as a limitation. All code, data
lineage, and 26 hermetic tests ship in the repository.

## 1. Introduction

[TBW: phage therapy against KPC/NDM Klebsiella and E. coli ST131; host-range
determination as the bottleneck; RBP-receptor biology (T4 gp37/38-OmpC,
T5 pb5-FhuA, lambda J-LamB, Klebsiella tail-spike depolymerases vs capsule);
limits of k-mer genome methods.]

## 2. Data

### 2.1 Phage genomes and RBP mining
1,439 RBP records from 922 complete genomes (esearch title queries, GenBank
efetch, rate-limited with exponential backoff; full lineage in
data/raw/genbank/). Receptor-binding CDS are selected by product keywords
(tail fiber/spike, receptor-binding, adhesin) with explicit exclusion of
assembly chaperones, connectors and catalysts, which do not contact the host.

### 2.2 Host receptor panel
14 receptors: E. coli LamB, FhuA, OmpC, OmpA, OmpF, Tsx, BtuB, FepA, TolC,
BamA; K. pneumoniae OmpK35, OmpK36, LamB, OmpA (UniProt accessions verified
live; AlphaFold DB v6 models for all 14; Table 1).

### 2.3 Task construction
Positive pair = phage x annotated host panel; negative = same phage x the
alternative species panel. Group split by phage accession (80/10/10) —
no RBP sequence appears in two splits. We prove (derivations.md §1) that
additive scoring rules are degenerate on this paired design in the balanced
case, which disqualifies naive dual-encoder baselines with additive scores.

## 3. Methods

### 3.1 Features
Per-residue 23-channel profiles: 20-dim one-hot + Kyte-Doolittle hydropathy,
Zamyatnin volume, formal charge at pH 7 (tables locked against literature by
unit tests — two table misalignments were caught by the test suite during
development and are recorded in the git history).

### 3.2 Models
CNN encoder (3 conv layers, mean+max pooling); MPNN GNN over C-alpha contact
graphs (cutoff 10 A) for the host side when structures are used;
bilinear compatibility s(p,h) = e_p^T W e_h + b; hybrid head over Hadamard
and difference dipeptide-spectrum interaction features; learned 2-way convex
combination. Baselines: Hadamard k-mer logistic regression (PHP-class),
additive LR (provably degenerate, §2.3).

### 3.3 Docking
Rigid grid search: 48 Fibonacci-sphere rotations x 48 surface anchors,
score = tight-band contacts (4.5-9 A) - 50 x clashes (<4.5 A) + 0.1 x clipped
Coulomb. Positive control: PDB 8A8C (T5 pb5-FhuA, cryo-EM 3.1 A).

## 4. Results

### 4.1 Benchmark (Table 2)
| model | test AUROC | AUPRC | acc | choice acc |
|---|---|---|---|---|
| Hadamard k-mer LR | 0.910 | 0.915 | 0.867 | - |
| CNN-bilinear (seq only) | 0.609 | 0.596 | 0.500 | 0.608 |
| **Hybrid CNN + k-mer (ours)** | **0.973** | **0.974** | **0.930** | **0.924** |

### 4.2 Discovery screen
17 candidate host-switch phages (Table 3), 11 in the held-out test split:
five E. coli-annotated phages predicted Klebsiella-tropic via OmpK36/OmpK35
(p > 0.99: PZ797503, PZ655591, PZ683213, OM867527, PQ821741), and
K. pneumoniae-annotated phages predicted E. coli-tropic via BtuB/TolC
(PV833093, PZ278545, OR090992). Each prediction is falsifiable by plaque
assay; disagreement with the GenBank annotation is itself informative
(either expanded host range or metadata error — both publishable outcomes).

### 4.3 Docking positive control (honest negative)
The coarse rigid grid fails to recover the 8A8C interface (top-pose iRMSD
40 A vs experimental). Contact-count scoring over large proteins rewards
burial; flexibility and loop refinement are needed. Docking is therefore
used only as a qualitative filter in this version.

## 5. Discussion
[TBW: composition vs local-signal analysis via combiner weights; OmpK36
porin-loss in KPC strains and implications for the predicted mediation;
polyvalence literature; limits: species-level ground truth, panel of 14
receptors, no capsule polysaccharide model.]

## 6. Reproducibility
Repo: mega27-08-phage-design (private). 26 hermetic tests; live calls only
in scripts; all raw data cached with lineage. Sandbox: 2 CPU / 1.9 GB RAM.

## References
[TBW: Kyte & Doolittle 1982; Zamyatnin 1972; Gilmer et al. 2017; WIsH
(Galiez 2017); PHP (Lu 2021); DeepHost 2022; iPHoP 2023; 8A8C/8B14;
Flores et al. on phage-bacteria interaction networks.]
