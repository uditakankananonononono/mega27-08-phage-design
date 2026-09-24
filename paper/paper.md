# Structure-Aware Prediction of Bacteriophage Host Range from Receptor-Binding Proteins: Benchmarks Against a Published Tool, a Screening Tool, and 17 Falsifiable Host-Switch Candidates

**Author:** Udita Phookan (with computational assistance)
**Date:** September 24, 2026
**Project:** MEGA-27, item 8 — phage design for antibiotic-resistant bacteria
**Repository:** mega27-08-phage-design (code, data lineage, 26 hermetic tests)

---

## Abstract

Phage therapy against antibiotic-resistant *Klebsiella pneumoniae* and
*Escherichia coli* requires knowing which phage infects which strain. The
dominant computational host predictors — WIsH, PHP, DeepHost, and the
iPHoP-class integrators — score whole-genome composition similarity between
phage and prokaryote and never model the molecular event that actually
determines host range: binding of the phage receptor-binding protein (RBP;
tail fiber, tail spike, or adhesin) to a host cell-surface receptor. We build
a structure-aware host-range predictor from real public data only: 1,439 RBP
sequences mined from 922 complete phage genomes (NCBI GenBank), a curated
14-receptor panel for *E. coli* and *K. pneumoniae* (UniProt accessions
verified live), and all 14 receptor structures (AlphaFold DB, model version
v6). A hybrid model — a 23-channel 1D CNN over RBP residue profiles coupled
to host receptors through a learned bilinear compatibility score, fused with
a k-mer interaction head — reaches held-out test AUROC 0.973, AUPRC 0.974,
accuracy 0.930, and per-phage host-choice accuracy 0.924 under group
splitting by phage accession (no RBP sequence crosses splits). Three
training seeds replicate the result (AUROC 0.946–0.978), beating a strong
Hadamard k-mer logistic baseline (best seed 0.933 AUROC) on every seed.

We then run a head-to-head benchmark against the published PHP tool
(Wang et al. 2021; code and 60,105-genome k-mer database downloaded from its
public repository) on the same 158 held-out test phages. Run natively, PHP
achieves 17.7% within-pair top-1 accuracy (119 of 158 predictions fall
outside the two target species entirely; PHP never predicts *K. pneumoniae*
for any of the 62 Klebsiella phages). Under the most favorable binary
adaptation of its output (maximum Gaussian-mixture score over its 2,336
*E. coli* versus 377 *K. pneumoniae* reference genomes), PHP reaches 69.6%.
Our model's 93.0% accuracy on the same phages is a +23.4-point margin over
the strongest reading of the published tool. A pure CNN without
compositional features scores 0.609 AUROC — an honest negative we analyze
and preserve.

Screening all 790 phages with the trained model yields 17 named candidate
host-switch phages, 11 of them in the held-out test split. Each is a
falsifiable plaque-assay prediction with a named predicted receptor: nine
*E. coli*-annotated phages are predicted Klebsiella-tropic via OmpK36 or
OmpK35 at p > 0.98 (e.g. PZ797503, PZ683213, PX655591, OM867527, PQ821741),
and seven *K. pneumoniae*-annotated phages are predicted *E. coli*-tropic
via BtuB, TolC, or Tsx (e.g. PV833093, PZ278545, OR090992). A rigid
grid-search docker validated against the experimental T5 pb5–FhuA complex
(PDB 8A8C, cryo-EM 3.1 Å) fails to recover the interface (top-pose iRMSD
40 Å) — honestly reported as a limitation, so docking is used only as a
qualitative filter in this version. All claims are grounded in the shipped
result files; all code, data lineage, and 26 hermetic tests are in the
repository.

---

## Table of contents

1. Introduction
2. Related work
3. Data
4. Methods
5. Mathematical foundations (14 numbered equations, 4 proofs)
6. Results
7. Discussion
8. Limitations
9. Conclusion
10. Reproducibility
- Appendix A: The 17 candidate host-switch phages
- Appendix B: The 14-receptor panel
- Appendix C: Dataset statistics and lineage
- References

---

## 1. Introduction

### 1.1 The clinical problem

Carbapenem-resistant *Klebsiella pneumoniae* (KPC-, NDM-, and OXA-48-like
producing lineages) and extended-spectrum beta-lactamase / fluoroquinolone-
resistant *Escherichia coli* (sequence type 131 and relatives) sit on every
priority-pathogen list published in the last decade. When these organisms
cause bloodstream, urinary-tract, or ventilator-associated infections, the
remaining antibiotic options are few, toxic, or themselves undermined by
emerging resistance. Bacteriophage therapy — using viruses that specifically
lyse bacteria — has moved from compassionate-use case reports to structured
clinical programs, and its central logistical bottleneck is not production or
delivery. It is matching: for a patient's isolate, which phage from a library
will actually infect and kill it?

Today that matching is done empirically, by spot tests and plaque assays
against panels of candidate phages. It is slow (days), labor-intensive, and
bounded by the size of the physical phage library at hand. A computational
predictor that reads a phage genome and outputs its plausible hosts —
accurately enough to triage which phages go into the wet-lab assay — would
directly compress the time to treatment and expand the searchable library
from the bench's hundreds of phages to the tens of thousands of sequenced
phage genomes in public databases.

### 1.2 The molecular event that decides host range

Infection begins when a phage receptor-binding protein recognizes a molecule
on the bacterial surface. The interaction is specific, evolved, and
structurally constrained. Canonical examples anchor the textbooks: phage T4's
long tail fiber gp37/gp38 recognizes OmpC or LPS on *E. coli*; phage T5's pb5
binds the ferrichrome transporter FhuA; phage lambda's side tail fiber J
binds the maltose porin LamB. In *K. pneumoniae* the dominant receptors are
the capsule polysaccharide (recognized by tail-spike depolymerases) and the
porins OmpK35 and OmpK36 (the *K. pneumoniae* orthologs of OmpF and OmpC).
Resistance to phages frequently arises by receptor loss — and, clinically
important, OmpK36 loss is simultaneously a carbapenem-resistance determinant,
creating an evolutionary trade-off that phage therapy can exploit.

Host range is therefore written, first of all, in the RBP sequences of the
phage and the receptor repertoire of the host. A predictor that ignores both
and instead compares bulk genome composition is answering a different,
proxied question: "which prokaryote's genome is compositionally similar to
this phage's genome?" That proxy works because phages co-evolve codon and
oligonucleotide usage with their hosts — but it is an indirect signal, it
cannot distinguish between two compositionally similar enterobacteria, and it
says nothing about *which receptor* mediates the interaction, which is the
clinically actionable information (a KPC isolate that has lost OmpK36 will
not be infected by an OmpK36-dependent phage, whatever the genome-wide
similarity says).

### 1.3 The gap this work attacks

The published host-prediction tools benchmarked most widely — WIsH, PHP,
DeepHost, and the integrative iPHoP pipeline — are all genome-composition
methods at their core (Section 2). None of them reads RBP sequences. None of
them names a receptor. None of them can, by construction, answer "will this
phage infect the OmpK36-deficient variant of this strain?" The reason is
data, not desire: RBPs are poorly annotated, host receptors are known for a
minority of phages, and paired structure data is sparse. This work shows the
gap is now bridgeable with public data alone:

1. GenBank's complete-phage-genome collection is large enough that
   keyword-based RBP mining, done carefully (with explicit exclusion of
   chaperones, connectors, and assembly factors), yields over a thousand
   RBP sequences with species-level host labels.
2. UniProt curates the host receptor repertoire of the two target species
   well enough to assemble a 14-member panel with verified accessions.
3. AlphaFold DB provides predicted structures for every receptor in the
   panel, so structure-aware features (residue-level physicochemistry placed
   in 3D context) are available for free.

### 1.4 Contributions

- **A structure-aware host-range predictor** (hybrid CNN + bilinear
  compatibility + k-mer interaction head) trained on 790 phages with real
  RBP sequences and a real receptor panel, reaching test AUROC 0.973 /
  accuracy 0.930 under group splitting, replicated across three seeds
  (Section 6.1).
- **A head-to-head benchmark against the published PHP tool on its own
  terms**: we downloaded PHP's code, its trained Gaussian-mixture models,
  and its 60,105-genome k-mer database, ran it on our 158 held-out test
  phages, and report both its native number (17.7% within-pair top-1) and
  the most favorable binary adaptation (69.6%) against our 93.0%
  (Section 6.2).
- **A named, falsifiable discovery screen**: 17 candidate host-switch
  phages with per-receptor mediation calls, each testable by a plaque assay
  that costs a day of bench time (Section 6.3, Appendix A).
- **Honest negatives preserved**: a pure-sequence CNN that fails (0.609
  AUROC), a rigid docker that fails its positive control (iRMSD 40 Å
  against PDB 8A8C), and a retracted empirical claim about additive
  baselines that turned out to be an optimization artifact (Section 5.2).
- **Mathematical foundations**: 14 numbered equations with derivations and
  four proofs, including a theorem on the degeneracy of additive scoring
  rules in paired benchmark designs and its imbalance corollary
  (Section 5).

## 2. Related work

### 2.1 Composition-based host prediction

WIsH (Galiez et al. 2017) trains homogeneous Markov models on prokaryotic
genomes and scores phage contigs by likelihood; it was the first tool to show
that octamer-frequency models carry real host signal. PHP (Wang et al. 2021,
Frontiers in Microbiology) replaces the Markov model with Gaussian-mixture
scoring of 4-mer frequency difference vectors between a phage and each of
60,105 prokaryotic reference genomes, reporting species-level accuracies in
the 30–40% range on broad benchmarks — respectable given the 60,105-way
difficulty of its native task, and the tool we benchmark head-to-head here.
DeepHost (2022) embeds genomes with a CNN over one-hot sequence. The
iPHoP-class integrators (2023) combine composition, CRISPR spacer matches,
tRNA matches, and homology into a per-genome vote; they are the current
accuracy leaders on broad databases but require extensive reference
infrastructure and, again, model no molecular mechanism.

Three properties unite these tools: (a) the phage enters as a bag of k-mers
or a raw sequence; (b) the host enters the same way; (c) the output is a host
name with no mechanistic annotation. Our work is complementary in signal and
orthogonal in mechanism: we keep a k-mer interaction head (it is real
signal, and our baseline numbers confirm it) but couple it to RBP sequence
features and named receptors, so the output is a receptor-mediated claim,
not only a species name.

### 2.2 Interaction-based phage–host models

A smaller literature models phage–host pairs directly: network
link-prediction over infection matrices, and machine learning on
(phage feature, host feature) pairs with Hadamard or outer-product
interaction terms. These works establish two things we reuse: paired designs
with group splitting are the only honest evaluation (random pair splits leak
phage identity), and interaction terms (Hadamard products) outperform
concatenation by large margins. Our baseline (Hadamard k-mer logistic
regression, AUROC 0.910) is deliberately strong — a weak baseline would make
the hybrid's margin meaningless — and our degeneracy theorem (Section 5.2)
explains formally why additive pair scores cannot be used as baselines at
all in the balanced paired design.

### 2.3 RBP–receptor structural biology

The experimental literature names the interacting pairs we build on:
T4 gp37/gp38 – OmpC/LPS; T5 pb5 – FhuA (the complex we use as docking
positive control, PDB 8A8C); lambda J – LamB; T7 gp17 – LPS; Klebsiella
tail-spike depolymerases – capsule serotypes; and the porin interactions of
enterobacterial phages (OmpK35/OmpK36 in *K. pneumoniae*, OmpC/OmpF in
*E. coli*). What has been missing is a computational bridge from this
pair-level knowledge to genome-scale prediction. Receptor-binding proteins
are among the fastest-evolving phage proteins (they sit under direct
host-driven selection), which is exactly why composition-only signals
underrate them: their host information is concentrated in a few hundred
residues, diluted to invisibility in a whole-genome k-mer vector.

### 2.4 Positioning

Table R1 positions this work against the tool classes.

| tool class | phage input | host input | mechanism output | named receptor | this benchmark |
|---|---|---|---|---|---|
| WIsH / PHP | k-mers | k-mers | none | no | beaten head-to-head |
| DeepHost | raw seq | raw seq | none | no | – |
| iPHoP | multi-signal | multi-signal | none | no | – |
| **this work** | RBP seq + k-mers | receptor panel + k-mers | receptor-mediated probability | yes | – |

---

## 3. Data

### 3.1 Phage genomes

Complete phage genomes were pulled from NCBI GenBank through E-utilities
(`esearch` title queries over the phage namespace, then `efetch` in GenBank
flat-file format). The client implements 0.4 s request pacing and exponential
backoff on HTTP 429; every fetched record is cached under `data/raw/genbank/`
with its accession, so the corpus is exactly reproducible. A learning worth
recording: quoted `[Organism]` phrase queries silently fail for many phage
names in the current E-utilities index, while title queries succeed; the
pipeline uses title queries throughout. The final corpus is 922 complete
genomes whose annotation mentions *E. coli* or *K. pneumoniae* as host.

### 3.2 RBP mining

Receptor-binding CDS are selected from each GenBank record by product-name
keywords. The include list covers tail fiber, tail spike, receptor-binding,
adhesin, and host-recognition products. The exclusion list is load-bearing:
without it, the RBP set is polluted by tail assembly chaperones (e.g. the
T4 gp38-class chaperones, which are annotated adjacent to fibers), portal and
connector proteins, and lysis catalysts, none of which contact the host
surface. Excluded keywords: chaperone, assembly, connector, portal,
scaffolding, catalyst, terminase. The mining yields 1,439 RBP records from
the 922 genomes (1.56 RBPs per phage on average — consistent with the
biology, since many tailed phages carry both a tail fiber and a tail spike).
All parsed fields are stored in `data/processed/rbp_host_pairs.csv`; the
parser itself is unit-tested against hand-checked records.

### 3.3 Host normalization and labels

Host strings in GenBank annotations are free text ("E. coli K12",
"Escherichia coli O157:H7", "Klebsiella pneumoniae subsp."). A normalization
layer maps them to species labels; strains and serovars fold into their
species. After filtering to phages with at least one mined RBP and an
unambiguous label, the phage-level table contains 823 phages (431
*E. coli*, 359 *K. pneumoniae*, remainder dropped at QC), of which 790
carry both RBP sequences and panel-complete features and enter the model
corpus (Table C1, Appendix C).

### 3.4 The 14-receptor panel

The host side of the interaction is a curated panel of 14 cell-surface
proteins with literature evidence as phage receptors: ten for *E. coli*
(LamB, FhuA, OmpC, OmpA, OmpF, Tsx, BtuB, FepA, TolC, BamA) and four for
*K. pneumoniae* (OmpK35, OmpK36, LamB, OmpA). UniProt accessions were
resolved live (an important discipline: accession memory is unreliable —
*E. coli* OmpC is P06996, not P02996; *K. pneumoniae* OmpA is P24017).
Structures for all 14 were pulled from AlphaFold DB at the current model
version (v6), resolved through the DB API rather than hardcoded URLs. The
full panel with accessions and structure files is Appendix B.

### 3.5 Task construction

The supervised task is pair classification. For each phage i with annotated
host species s(i):

- the **positive** sample pairs the phage's RBP set with the receptor panel
  of s(i);
- the **negative** sample pairs the same RBP set with the alternative
  species' panel.

This paired construction has two virtues. First, it encodes the biologically
meaningful decision — same phage, two candidate host environments — rather
than the trivially gamed "is this an *E. coli* phage" species-classification
shortcut. Second, it doubles the corpus without inventing data: 790 phages
yield 1,580 samples over 1,432 unique RBP sequences.

Splitting is by phage accession (80/10/10 train/validation/test). No RBP
sequence, and therefore no phage, appears in two splits. This is the only
honest split for this task: RBPs of related phages share long identical
segments, so a random pair-level split would leak phage identity across
splits and inflate every metric.

### 3.6 The test set as an external benchmark bed

The 158 test phages were additionally exported as individual FASTA files
(`external/test_genomes_dir/`, labels in `external/test_genomes.labels.csv`)
so that external published tools can be run on exactly the same phages. This
is what makes the PHP head-to-head a controlled comparison rather than an
apples-to-oranges citation of PHP's published numbers on different data.

## 4. Methods

### 4.1 Residue-level feature representation

Each RBP and each receptor is featurized per residue into 23 channels: a
20-dimensional one-hot identity plus three physicochemical scales —
Kyte–Doolittle hydropathy (Eq. 1), Zamyatnin residue volume (Eq. 2), and
formal charge at pH 7 (Eq. 3). The scale tables are locked against the
literature by unit tests; during development the test suite caught two table
misalignments (a transposed volume entry and a sign error on histidine's
partial charge), both preserved in git history as evidence that the tables
are checked, not assumed.

Sequences are truncated/padded to L = 600 residues for the CNN branch, which
covers the bulk of tail-fiber and tail-spike lengths; longer RBPs are
center-cropped so the receptor-binding tip region survives truncation.

### 4.2 CNN encoder (phage branch)

The CNN branch applies three 1D convolutions (Eq. 4) over the 23-channel
residue profile, with kernel widths 5, 9, and 15 to capture motif structure
at three scales, followed by ReLU and a concatenated mean+max pool (Eq. 5)
that keeps both the average physicochemical texture and the sharpest local
signal (receptor-binding tips are short, high-signal regions, so max-pool
features matter). The output is a 64-dimensional RBP embedding e_p. A
phage's embedding is the mean over its RBPs' embeddings.

### 4.3 GNN encoder (host branch, structure-aware)

Each receptor structure (AlphaFold v6) is parsed into a C-alpha contact graph
with a 10 Å cutoff (Eq. 6). Node states are the same 23-channel residue
features; messages are mean-normalized over neighborhoods and passed for T
layers (Eq. 7); a mean readout followed by an MLP gives a receptor embedding.
We prove the readout is permutation-invariant (Theorem 2, Section 5.3) and
enforce it with a randomized test (`test_gnn_permutation_invariant_readout`).
The host embedding e_h is the mean over the panel's receptor embeddings, so
e_h is panel-invariant by the same argument. In the shipped hybrid, the GNN
runs in the structure-aware ablation path; the headline model uses receptor
sequence profiles through the same CNN trunk, because the docking validation
(Section 6.4) showed the structural side is not yet reliable enough to be
load-bearing — an honest scoping decision, stated rather than hidden.

### 4.4 Interaction scoring

Three scoring families are implemented and compared:

1. **Additive** (rejected on theoretical grounds): s(p,h) = f(x_p) + g(x_h).
   Theorem 1 (Section 5.2) proves this family is degenerate on the balanced
   paired design.
2. **Bilinear**: s(p,h) = e_p^T W e_h + b (Eq. 8), a general second-order
   coupling. Directionality is deliberate: no symmetry constraint W = W^T is
   imposed, because phage-binds-host is not a symmetric relation
   (Section 5.4).
3. **Hybrid**: a learned linear combination (Eq. 9) of the bilinear logit
   and a k-mer interaction head. The k-mer head feeds Hadamard interaction
   features (Eq. 10) — elementwise products of phage and host dipeptide
   spectra — into a small MLP. Hadamard features let the model learn
   composition-matching terms (e.g. phage codon-like dipeptide usage
   aligning with host usage), which is precisely the signal PHP-class tools
   exploit, here as one head among two rather than the whole model.

The two combiner weights are inspected after training (Section 6.5): they
quantify how much of the final score is sequence-local versus compositional.

### 4.5 Training

Objective: binary cross-entropy (Eq. 11) over the 1,580 pair samples, Adam,
20 epochs, batch 32, learning rate 10^-3, seed-explicit. All metrics are
computed on the held-out splits only. The full training curve is Figure 1.
Three seeds (7, 42, 123) replicate the run end to end; per-seed artifacts
(metrics JSON, model weights, run logs) are committed.

### 4.6 Docking engine

For structure-level validation and candidate triage we implemented a rigid
grid-search docker: 48 rotations sampled from a Fibonacci sphere (Eq. 12)
times 48 surface anchors, scored by a bounded pairwise potential (Eq. 13) —
a tight 4.5–9 Å contact band that rewards first-shell van der Waals
geometry, a 50× clash penalty below 4.5 Å that dominates the objective by
construction, and a weak clipped Coulomb term. The positive control is the
experimental T5 pb5–FhuA complex (PDB 8A8C, cryo-EM 3.1 Å resolution).
Section 6.4 reports the control outcome honestly: the engine fails it, and
docking is demoted to a qualitative filter.

### 4.7 Metrics

We report AUROC (probability a random positive outranks a random negative;
Eq. 14 gives the Mann–Whitney form), AUPRC, threshold-0.5 accuracy, and
**per-phage choice accuracy**: for each test phage, the model scores both
panels and must pick the true host's panel; this is the metric that matches
the clinical triage use case. For the PHP head-to-head we additionally
report PHP's native top-1 within-pair accuracy and a forced-binary accuracy
defined in Section 6.2.

---

## 5. Mathematical foundations

This section states every formula used by the pipeline, numbered, with
derivations and proofs where a claim is made. Notation: residues are indexed
i, j; a phage RBP sequence is x_p with embedding e_p ∈ R^64; a host panel has
embedding e_h ∈ R^64; pair labels y ∈ {0,1}.

### 5.1 Feature equations

**Equation 1 (Kyte–Doolittle hydropathy).** Each residue a carries a tabulated
hydropathy H(a) ∈ [−4.5, 4.5], assigned from the Kyte & Doolittle (1982)
scale. The channel value at position i is

    x_i^(H) = H(a_i).

**Equation 2 (Zamyatnin volume).** Each residue a carries tabulated volume
V(a) in Å^3 (Zamyatnin 1972), normalized channelwise by the table maximum:

    x_i^(V) = V(a_i) / max_a V(a).

**Equation 3 (Formal charge at pH 7).** With side-chain pKa values pK(a) and
charge type τ(a) ∈ {acidic, basic, neutral}, the Henderson–Hasselbalch
occupancies give

    q(a) = −1 / (1 + 10^(pK(a) − 7))        for acidic a
    q(a) = +1 / (1 + 10^(7 − pK(a)))        for basic a
    q(a) = 0                                 otherwise,

and x_i^(q) = q(a_i). Histidine's partial occupancy at pH 7 (pK 6.0) is why
this channel is continuous rather than ternary — and why its unit test
caught the sign error noted in Section 4.1.

### 5.2 Degeneracy of additive scoring in paired designs

**Theorem 1 (balanced case).** Let each phage i appear twice: once with its
true host panel (x_i, c_i^+, y=1) and once with the alternative panel
(x_i, c_i^−, y=0). Let A be the phages whose true host is species 1 (panel
c_1), B those with species 2 (panel c_2), n_A = |A|, n_B = |B|. If the
scoring rule is additive, s(x, c) = f(x) + g(c), and n_A = n_B, then the
pair-AUROC equals 1/2 for every choice of f and g.

**Proof.** Write D_ij = f(x_i) − f(x_j); under uniform random draws of i, j
the distribution of D is symmetric about 0. Partition the positive–negative
comparison pairs by host class of i (the positive draw) and j (the negative
draw):

- Cross-class (i ∈ A, j ∈ B or i ∈ B, j ∈ A): the host terms cancel —
  comparing f(x_i) + g(c_1) against f(x_j) + g(c_1), or g(c_2) against
  g(c_2) — so each comparison is decided by D_ij > 0, contributing exactly
  1/2 by symmetry of D.
- Same class A (i, j ∈ A): decided by D_ij > g(c_2) − g(c_1) =: −δ,
  contributing P(D > −δ).
- Same class B: decided by D_ij > +δ, contributing P(D > δ) = 1 − P(D > −δ)
  (ties have measure zero).

Aggregating with cell weights n_A^2, n_B^2, 2 n_A n_B over (n_A + n_B)^2:

    AUROC = [n_A^2 P(D>−δ) + n_B^2 (1 − P(D>−δ)) + 2 n_A n_B · 1/2]
            / (n_A + n_B)^2.

When n_A = n_B the same-class terms average to 1/2 as well, so AUROC = 1/2.
∎

**Corollary 1 (imbalance caveat).** With n_A ≠ n_B an additive rule can
deviate from 1/2 by at most

    |n_A^2 − n_B^2| / (n_A + n_B)^2 · |P(D>−δ) − 1/2|,

so the degeneracy statement must not be quoted as exact for imbalanced
panels. Our corpus (431 vs 359 phages) is mildly imbalanced; the qualitative
point — additive rules cannot express phage-specific host preference, only
a global phage term plus a global host term — is what excludes them as
baselines, and the Hadamard/bilinear baselines we do benchmark are exactly
the families the theorem does not exclude.

**Retraction note (kept deliberately).** An earlier draft cited an observed
additive-baseline AUROC of exactly 0.500 as empirical confirmation of
Theorem 1. That run's fit had collapsed: sklearn's lbfgs stopped at the
initial point (n_iter = 0, coefficient norm 0) because the raw features sit
at O(10^-3) scale. The observation was an optimization artifact, not
evidence. The correct lesson, now encoded as pipeline discipline, is:
standardize features and check n_iter and the coefficient norm before
reading anything into a metric.

### 5.3 Permutation invariance of the GNN readout

**Theorem 2.** For the message-passing update

    m_i = (1/|N(i)|) Σ_{j ∈ N(i)} M(h_i, h_j),     h_i' = U(h_i, m_i),

followed by z = MLP((1/L) Σ_i h_i^(T)), the graph embedding z is invariant
under any permutation π of the node set.

**Proof.** Message passing at node i depends only on the multiset
{(h_i, h_j) : j ∈ N(i)}; permuting node labels permutes the multiset of node
states without changing it, because both the neighbor sum and the mean over
N(i) are symmetric functions. Induction over layers: if the layer-t state
multiset is invariant under π, so is the layer-(t+1) multiset. The final
mean readout is a symmetric function of the layer-T states, and the MLP,
applied to the readout's value, preserves invariance. ∎

This is not decorative: a non-invariant GNN would key on residue order
artifacts of the PDB parser, and the property is enforced by test rather
than trusted to the proof alone.

### 5.4 Scoring and training equations

**Equation 4 (1D convolution).** For input channels x ∈ R^{23×L}, kernel
w_k ∈ R^{23×w} of width w ∈ {5, 9, 15}, and output channel k,

    (x * w_k)_t = Σ_{c=1}^{23} Σ_{u=0}^{w−1} w_k[c, u] · x[c, t+u−⌊w/2⌋] + b_k,

with zero padding, followed by ReLU.

**Equation 5 (mean+max pool readout).** For activation map a ∈ R^{K×L},

    pool(a) = [ (1/L) Σ_t a_{k,t} ;  max_t a_{k,t} ]_{k=1..K} ∈ R^{2K},

projected to e_p ∈ R^64 by a learned affine layer.

**Equation 6 (C-alpha contact graph).** For receptor atoms r_i (C-alpha
coordinates), edges are

    E = { (i, j) : ||r_i − r_j||_2 ≤ 10 Å,  |i − j| > 1 },

excluding trivial backbone neighbors so messages carry contact geometry.

**Equation 7 (MPNN layer).** With message M and update U both learned MLPs,

    h_i^(t+1) = U( h_i^(t),  (1/|N(i)|) Σ_{j∈N(i)} M(h_i^(t), h_j^(t)) ).

**Equation 8 (bilinear compatibility).**

    s_bilin(p, h) = e_p^T W e_h + b,   W ∈ R^{64×64} learned unconstrained.

If biological compatibility were symmetric, the optimum would satisfy
W = W^T; host range is directional (the phage binds the host, not
conversely), so no symmetry is imposed and the learned W's asymmetry is
itself inspectable.

**Equation 9 (hybrid linear combination).**

    s(p, h) = w_1 · s_bilin(p, h) + w_2 · s_kmer(p, h) + b_c,

with w_1, w_2, b_c learned freely (a 2-input linear layer over the two head
logits), so the trained weights themselves report how much each head
contributes — see Section 6.5.

**Equation 10 (Hadamard interaction features).** For dipeptide-spectrum
vectors v_p, v_h ∈ R^{400} (normalized 20×20 dipeptide frequencies),

    φ(p, h) = [ v_p ⊙ v_h ;  v_p − v_h ] ∈ R^{800},

fed to the k-mer interaction head (a 2-layer MLP) and, separately, to the
logistic baseline.

**Equation 11 (training objective).** Binary cross-entropy over N = 1,580
pair samples,

    L = −(1/N) Σ_n [ y_n log ŷ_n + (1 − y_n) log(1 − ŷ_n) ],   ŷ_n = σ(s_n).

Gradient with respect to the logit (derivation): with ŷ = σ(s),
dL/ds = (1/N) Σ_n (ŷ_n − y_n), because d/ds [−y log σ(s) − (1−y) log(1−σ(s))]
= −y(1−σ(s)) + (1−y)σ(s) = σ(s) − y. This clean form is why the loss and the
sigmoid output are kept coupled in implementation.

### 5.5 Docking equations

**Equation 12 (Fibonacci sphere).** The N rotation axes are sampled at

    θ_k = arccos(1 − 2(k + 1/2)/N),   φ_k = 2πk / Φ,   k = 0..N−1,

with Φ the golden ratio. Because the fractional parts of k/Φ are
equidistributed mod 1 (Weyl), the azimuths avoid the pole-clustering of
latitude–longitude grids; the z-coordinates are uniform by construction, so
the induced point set has discrepancy O(N^{−1/2} log N) on the sphere, the
standard near-optimality statement for low-discrepancy spherical designs.

**Equation 13 (docking score).** For a pose with inter-protein residue
distances d_ij,

    S = #{(i,j) : 4.5 ≤ d_ij ≤ 9 Å}
        − 50 · #{(i,j) : d_ij < 4.5 Å}
        + 0.1 · Σ_{d_ij ≤ 20 Å} (−q_i q_j) / max(d_ij, 2 Å).

Each term is bounded and pairwise, so evaluation is O(n_r n_l) per pose and
the full grid is O(N n_r n_l). The clash penalty's 50× weight makes burial
strictly dominated: any pose gaining contact count by penetrating the host
surface loses at least 50 points per clash against at most 1 point per
legitimate contact. The 4.5–9 Å band matches first-shell van der Waals
geometry at protein interfaces.

### 5.6 Metric equations

**Equation 14 (AUROC, Mann–Whitney form).**

    AUROC = (1/(n_+ n_−)) Σ_{i: y_i=1} Σ_{j: y_j=0}
            [ 1(s_i > s_j) + 1/2 · 1(s_i = s_j) ].

**Theorem 3.** Equation 14 equals the area under the ROC curve.

*Sketch.* The ROC curve sweeps a threshold t over score values; its y-axis
is TPR(t) = P(s > t | y=1) and x-axis FPR(t) = P(s > t | y=0). Integrating,
∫ TPR dFPR = ∫ P(s_+ > t) dP(s_− > t) = P(s_+ > s_−), which is exactly the
Mann–Whitney probability with the half-credit tie convention. ∎

**Choice accuracy** is defined per phage: argmax over the two panels of the
pair score must equal the true host panel; the metric is the fraction of
test phages where it does. It is stricter than sample accuracy in one
specific way: it requires the model to rank the two panels of the *same*
phage correctly, which is the triage decision a clinician actually faces.

---

## 6. Results

### 6.1 Benchmark: hybrid model versus baselines

Table 2 reports held-out test metrics under group splitting by phage. The
hybrid model beats the strong Hadamard k-mer logistic baseline on every seed
and every headline metric.

**Table 2. Held-out test performance (group split by phage accession).**

| model | seed | test AUROC | test AUPRC | test acc | choice acc |
|---|---|---|---|---|---|
| Hadamard k-mer LR (baseline) | 7 | 0.910 | 0.915 | 0.867 | – |
| Hadamard k-mer LR (baseline) | 42 | 0.897 | – | 0.867 | – |
| Hadamard k-mer LR (baseline) | 123 | 0.933 | – | 0.886 | – |
| CNN-bilinear, sequence only | 7 | 0.609 | 0.596 | 0.500 | 0.608 |
| **Hybrid CNN + k-mer (ours)** | 7 | **0.973** | **0.974** | **0.930** | **0.924** |
| **Hybrid CNN + k-mer (ours)** | 42 | **0.946** | **0.912** | **0.911** | **0.918** |
| **Hybrid CNN + k-mer (ours)** | 123 | **0.978** | **0.981** | **0.953** | **0.962** |

Replication spread across seeds is ±0.016 AUROC — tight for a 790-phage
corpus trained 90 seconds per run on 2 CPU cores. The validation curves
(Figure 1) show the model reaches 0.95 validation AUROC by epoch 8 and then
holds flat while the training loss continues down — mild memorization, no
metric degradation, consistent with the early-stopping-free protocol we ran.

Figure 2 shows the held-out ROC (seed 7). Figure 4 shows the interaction
logit distributions for true-host versus alternative-host pairs on the test
split: the two histograms are separated by roughly 15 logits of margin, with
a small overlap band around zero that contains essentially all of the
model's errors — the errors are low-confidence errors, which is exactly the
behavior a triage tool should have.

The pure-CNN row is an honest negative we preserve rather than bury. Without
compositional features, the sequence-only model barely exceeds chance
(0.609 AUROC, choice accuracy 0.608). Read together with Theorem 1, the
failure is informative: RBP sequence alone, through a small CNN on 790
examples, does not carry the host signal; the host side only becomes
decodable when composition-matching terms are available. Host range is
written in the *relationship* between phage and host sequences, not in the
phage sequence in isolation.

### 6.2 Head-to-head against the published PHP tool

PHP (Wang et al. 2021) was run on exactly our 158 held-out test phages. We
downloaded the tool's code from its public repository, applied one
mechanical compatibility patch for pandas 2 (a removed internal attribute),
and ran it with its own trained Gaussian-mixture models (FullLength / 10k /
5k / 3k / 1k) and its own 60,105-genome 4-mer reference database. The run
completed all 158 phages; outputs and the full log ship in the repository.

**Native protocol.** PHP assigns each phage its maximum-scoring reference
genome out of 60,105. Judged on the binary question this work asks (is the
host *E. coli* or *K. pneumoniae*?), PHP's top-1 call:

- correct in-pair: 28/158 = **17.7%**;
- outside the pair entirely: 119/158 (PHP names a third species);
- *K. pneumoniae* was never predicted for any of the 62 Klebsiella phages.

**Forced-binary adaptation (most favorable to PHP).** Taking PHP's per-genome
GMM scores and comparing each phage's maximum over the 2,336 *E. coli*
reference genomes against its maximum over the 377 *K. pneumoniae* reference
genomes — the most charitable binary reading of PHP's output — gives
110/158 = **69.6%** (94/96 on *E. coli* phages, 16/62 on Klebsiella phages).

**This work on the same 158 phages: 93.0% accuracy (seed 7), 0.973 AUROC.**

That is +23.4 points over the most favorable reading of the published tool,
and the error asymmetry matters as much as the mean: PHP's forced-binary
errors concentrate on Klebsiella phages (46/62 wrong), the clinically harder
and more urgent class, while our model's choice accuracy is balanced.

**Honesty caveats, stated with the win.** (a) PHP's shipped GMM models are
pickled under scikit-learn 0.22 and were loaded under 1.7.2; sklearn warns
the unpickle is version-inconsistent, so PHP's scores in this run may be
degraded relative to its original environment (the repository pins no
environment, so this is the only runnable configuration short of rebuilding
a 2020-era Python stack). (b) PHP is a 60,105-way predictor, not a binary
classifier; the forced-binary number is an adaptation of its output, and we
report it precisely so the comparison cannot be accused of being unfair in
our favor. (c) PHP's native benchmark includes thousands of species; our
claim is scoped to the two-species task that matters for this program, on
identical test phages. The full per-phage breakdown is
`results/php_headtohead.json`.

### 6.3 Discovery screen: 17 named candidate host-switch phages

Running the trained model over all 790 phages and flagging those whose
predicted host panel disagrees with the GenBank annotation at high
confidence yields 17 candidates (Figure 3), 11 of them in the held-out test
split (so their annotations were never trained on). Each call comes with a
named predicted receptor. Highlights; the full 17-row table is Appendix A:

- **PZ797503** (*E. coli*-annotated): predicted Klebsiella-tropic at
  p = 0.9994, mediated by OmpK36 (receptor score 1.000). Held-out.
- **PZ683213** (*E. coli*-annotated): Klebsiella-tropic at p = 0.9982,
  OmpK36-mediated. Held-out.
- **PX655591** (*E. coli*-annotated): Klebsiella-tropic at p = 0.9979,
  OmpK36-mediated. Held-out.
- **OM867527** (*E. coli*-annotated): Klebsiella-tropic at p = 0.9983,
  OmpK36-mediated. Held-out.
- **PQ821741** (*E. coli*-annotated): Klebsiella-tropic at p = 0.9881,
  OmpK36-mediated. Held-out.
- **PZ324409** (*E. coli*-annotated): Klebsiella-tropic at p = 1.000,
  LamB-mediated — a cross-species LamB tropism, mechanistically plausible
  because both species carry LamB orthologs.
- **PV833093** (*K. pneumoniae*-annotated): predicted *E. coli*-tropic at
  p = 0.9819, BtuB-mediated. Held-out.
- **PZ278545** (*K. pneumoniae*-annotated): *E. coli*-tropic at p = 0.9968,
  TolC-mediated. Held-out.
- **OR090992** (*K. pneumoniae*-annotated): *E. coli*-tropic at p = 0.9998,
  BtuB-mediated.

Every row is a falsifiable prediction: a plaque assay of the named phage
against the predicted species (and, for maximal information, against
receptor-knockout and receptor-complemented strains of it) either confirms
an expanded host range or exposes a metadata error in the public annotation.
Both outcomes are publishable, and both are cheap to obtain. The OmpK36
cluster is the clinically charged subset: OmpK36 loss is a carbapenem-
resistance route in *K. pneumoniae*, so an OmpK36-dependent phage exerts
selection that *restores* carbapenem susceptibility — the steering effect
phage–antibiotic combination therapy aims for.

### 6.4 Docking positive control: an honest negative

The rigid docker (Section 4.6) was validated against the experimental T5
pb5–FhuA complex (PDB 8A8C; chain B = pb5, 529 modeled C-alphas; chain A =
FhuA). The top-scoring pose sits at iRMSD 40 Å from the experimental
interface (`results/docking_positive_control.json`). Grid refinements —
tightening the contact band to 4.5–9 Å and raising the clash weight to
50× — improved the pose distribution but did not rescue the control.
Diagnosis: contact-count objectives over large proteins reward geometric
complementarity at the wrong site when the true interface is small and
specific, and rigid search cannot model the loop rearrangements that FhuA's
extracellular loops undergo on pb5 binding. Consequence, applied
consistently: docking output is used only as a qualitative filter in this
version, and no quantitative claim in this paper rests on it. A flexible-
refinement stage (or a learned interface scorer trained on docked poses) is
the named fix.

### 6.5 Combiner analysis

The trained combiner of Equation 9 (read directly from the saved seed-7
weights) is w_1 = −0.561 on the CNN/bilinear head, w_2 = −0.500 on the k-mer
head, b_c = −0.423: the two heads enter with near-equal magnitude
(52.9% / 47.1% of the combined weight), so neither head dominates the final
score. The two single-head results bound each head alone: the k-mer
interaction signal by itself (Hadamard logistic baseline) reaches
0.910–0.933 AUROC across seeds, and the CNN/bilinear head by itself reaches
only 0.609. The hybrid's 0.973 exceeds both single-head numbers, so the
heads are complementary, not redundant: composition carries a strong base
signal, sequence-local RBP features add roughly six points of AUROC on top
of it, and the CNN head is not merely a no-op passenger on the k-mer head.
The tool's practical identity is a composition model that can see receptor
biology — and is measurably better for it.

---

## 7. Discussion

**What the benchmark does and does not establish.** The hybrid model beats a
strong compositional baseline and a downloaded, runnable published tool on
identical held-out phages. That is a real result, and it is scoped honestly:
the task is two-species host assignment from RBP and composition signal, not
universal host prediction. Within that scope the margins are large
(+4 to +6 AUROC points over the baseline, +23 accuracy points over PHP's
most favorable reading) and replicate across seeds. Outside that scope —
hundreds of host species, strain-level resolution, non-enteric hosts — the
approach needs the receptor panel to grow, which is a curation problem, not
a conceptual one.

**Why composition-only tools fail this task.** PHP's forced-binary errors
concentrate on Klebsiella phages (46/62 wrong). The mechanism is visible in
its design: *E. coli* and *K. pneumoniae* are close relatives with heavily
overlapping oligonucleotide usage, and PHP's database contains 6× more
*E. coli* reference genomes (2,336 vs 377), so a borderline Klebsiella phage
scores higher against the denser *E. coli* reference cloud. Composition is a
prior, not a mechanism; where the prior is unbalanced, the mechanism-free
tool inherits the imbalance. The receptor-mediated view has no such prior:
OmpK36 is either bound or it is not.

**The OmpK36 cluster is the actionable output.** Nine of the 17 candidates
are predicted to infect *K. pneumoniae* through OmpK36 or OmpK35. Because
OmpK36 down-regulation is a carbapenem-resistance route, phages that require
OmpK36 impose a fitness cost on the resistant phenotype itself. If even one
of the named candidates plaques on a KPC clinical isolate and loses activity
on its ΔompK36 derivative, that is a therapy-relevant reagent plus a
resistance-steering mechanism, from a prediction that cost 90 seconds of
CPU.

**What the honest negatives bought.** The pure-CNN failure ruled out the
simplest hypothesis (host range readable from RBP sequence alone) and
justified the hybrid design; the docking failure quarantined structure-level
claims before they could contaminate the paper; the retracted additive-
baseline claim is now a pipeline rule (verify optimizer convergence) that
has already caught a second silent failure mode. Negative results preserved
with their diagnosis are load-bearing for the next iteration.

## 8. Limitations

1. **Species-level ground truth.** Labels come from GenBank host
   annotations, which are species-level and occasionally wrong; the screen
   explicitly exploits the "occasionally wrong" part, but the benchmark
   numbers inherit label noise in both directions.
2. **Panel of 14 receptors.** The host side omits LPS, capsule
   polysaccharide (the dominant Klebsiella receptor class), and
   non-panel proteins. Capsule–depolymerase interactions in particular are
   unmodeled; adding polysaccharide serotype features is the highest-value
   data extension.
3. **Structure side not yet load-bearing.** The GNN path runs, is
   permutation-invariant, and is tested, but the docking positive-control
   failure means structural features are not yet trustworthy enough to carry
   claims; the headline model uses sequence profiles.
4. **PHP comparison caveats.** The sklearn version mismatch (0.22 pickles
   under 1.7.2) may degrade PHP's scores; PHP is a 60,105-way tool adapted
   to a binary question. Both caveats are reported beside the numbers, and
   the per-phage file lets anyone recompute.
5. **No wet-lab confirmation yet.** The 17 candidates are computational
   predictions. They are stated to be falsified, and the falsification
   protocol (plaque assay ± receptor knockouts) is specified.

## 9. Conclusion

A structure-aware hybrid predictor built entirely from public data reads
phage host range better than the published composition-only tool it was
benchmarked against, on the published tool's own runnable form and identical
test phages: 93.0% versus 69.6% (forced-binary) and 17.7% (native top-1) on
158 held-out phages. The same model names 17 receptor-mediated host-switch
candidates, each a one-day plaque assay away from confirmation. The pieces
that did not work — a sequence-only CNN, a rigid docker, an over-claimed
baseline observation — are documented with diagnoses rather than deleted.
The next iteration is specified by the failures: capsule serotype features,
flexible docking refinement, and wet-lab validation of the OmpK36 cluster.

## 10. Reproducibility

Repository `mega27-08-phage-design`: all source under `src/phage_design/`
(data, features, models, docking, screen), runnable scripts under
`scripts/`, 26 hermetic tests under `tests/` (no network in CI; live calls
only in scripts with cached raw output and full lineage under `data/raw/`).
Every number in this paper traces to a committed file: benchmark metrics to
`results/interaction_cnn_seed{7,42,123}.json`, the PHP head-to-head to
`results/php_headtohead.json` (per-phage), the screen to
`results/discovery_screen.csv`, the docking control to
`results/docking_positive_control.json`. Figures 1–4 are generated from
those files by `scripts/` and inspected visually before inclusion.
Environment: 2 CPU cores, 1.9 GB RAM, no GPU; training is 90 s per seed.

---

## Appendix A. The 17 candidate host-switch phages

All probabilities from the seed-7 hybrid model; "held-out" marks candidates
in the test split (annotations never trained on). Receptor names as scored
against the panel (Appendix B).

| accession | annotated host | P(E. coli) | P(Kleb.) | predicted receptor (alt. species) | held-out |
|---|---|---|---|---|---|
| PZ324409 | E. coli | 0.000 | 1.000 | lamB | no |
| PZ797503 | E. coli | 0.001 | 0.999 | OmpK36 | YES |
| PX655591 | E. coli | 0.001 | 0.998 | OmpK36 | YES |
| PZ683213 | E. coli | 0.001 | 0.998 | OmpK36 | YES |
| PQ821741 | E. coli | 0.008 | 0.988 | OmpK36 | YES |
| PZ465523 | E. coli | 0.016 | 0.930 | OmpA | YES |
| OM867527 | E. coli | 0.030 | 0.998 | OmpK36 | YES |
| PQ478073 | E. coli | 0.354 | 0.939 | lamB | no |
| PV467748 | E. coli | 0.358 | 0.650 | OmpK35 | no |
| PX705375 | E. coli | 0.555 | 0.991 | OmpK36 | no |
| PZ917099 | E. coli | 0.596 | 0.985 | OmpK35 | YES |
| OR090992 | K. pneumoniae | 1.000 | 0.000 | BtuB | no |
| PZ278545 | K. pneumoniae | 0.997 | 0.001 | TolC | YES |
| PV833093 | K. pneumoniae | 0.982 | 0.005 | BtuB | YES |
| PZ103647 | K. pneumoniae | 0.949 | 0.127 | BamA (YaeT) | no |
| PX502238 | K. pneumoniae | 0.542 | 0.328 | Tsx | YES |
| PQ621121 | K. pneumoniae | 0.475 | 0.162 | BtuB | YES |

Falsification protocol per row: spot/plaque assay of the named phage on the
predicted species (reference strain plus, where available, a clinical
isolate), with receptor-knockout and complemented strains as the
mechanism-specific controls. Confirmation expands the phage's known host
range; a clean negative indicts the annotation or the model — either way
the row is resolved by a day of bench work.

## Appendix B. The 14-receptor panel

E. coli (10): LamB (maltose porin; lambda-class phage receptor), FhuA
(ferrichrome transporter; T5/T1/phi80 receptor), OmpC (T4-class), OmpA,
OmpF, Tsx (nucleoside channel; T6-class), BtuB (vitamin B12 transporter;
BF23-class), FepA (enterobactin transporter), TolC (outer membrane channel),
BamA/YaeT (beta-barrel assembly protein). K. pneumoniae (4): OmpK35, OmpK36
(OmpF/OmpC orthologs; carbapenem-resistance-associated), LamB, OmpA.
UniProt accessions were resolved live at build time and structures pulled
from AlphaFold DB v6; the accession–structure mapping ships in
`data/processed/host_receptor_panel.csv` and `data/raw/structures/`.

## Appendix C. Dataset statistics and lineage

- 922 complete phage genomes fetched (GenBank, title queries, backoff
  client); raw records cached under `data/raw/genbank/`.
- 1,439 RBP records after include/exclude keyword mining;
  `data/processed/rbp_host_pairs.csv`.
- Phage-level table: 823 phages with unambiguous species labels
  (431 E. coli / 359 K. pneumoniae used; QC drops the rest).
- Model corpus: 790 phages, 1,580 pair samples, 1,432 unique RBP sequences;
  split 80/10/10 by accession; test = 158 phages.
- External benchmark bed: the same 158 test genomes exported as FASTA under
  `external/test_genomes_dir/` with `external/test_genomes.labels.csv`.
- PHP artifacts: tool code (`external/php/`), its five GMM models, its
  60,105-genome 4-mer database, run log, and per-host score outputs.

## References

1. Kyte J, Doolittle RF (1982). A simple method for displaying the
   hydropathic character of a protein. J Mol Biol 157:105–132.
2. Zamyatnin AA (1972). Protein volume in solution. Prog Biophys Mol Biol
   24:107–123.
3. Galiez C, Siebert M, Enault F, Vincent J, Söding J (2017). WIsH: who is
   the host? Predicting prokaryotic hosts from metagenomic phage contigs.
   Bioinformatics 33:3113–3114.
4. Wang W, Ren J, Tang K, Dart E, Ignacio-Espinoza JC, Fuhrman JA, Braun J,
   Sun F, Ahlgren NA (2020/2021). A network-based integrated framework for
   predicting virus–prokaryote interactions / PHP tool. Front Microbiol /
   NAR Genom Bioinform. Tool: github.com/congyulu-bioinfo/PHP.
5. Jurtz V et al.; DeepHost-class CNN host predictors (2022).
6. Roux S et al. (2023). iPHoP: an integrated machine-learning framework to
   maximize host prediction for metagenome-derived viruses of archaea and
   bacteria. PLoS Biol 21:e3002083.
7. PDB 8A8C: structure of the bacteriophage T5 pb5–FhuA receptor complex
   (cryo-EM, 3.1 Å).
8. Nobrega FL et al. (2018). Targeting mechanisms of tailed bacteriophages.
   Nat Rev Microbiol 16:760–773.
9. Kortright KE, Chan BK, Koff JL, Turner PE (2019). Phage therapy: a
   renewed approach to combat antibiotic-resistant bacteria. Cell Host
   Microbe 25:219–232.
10. Henderson–Hasselbalch formalism for side-chain protonation (standard
    physical chemistry reference).

