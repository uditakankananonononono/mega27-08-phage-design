# mega27-08-phage-design

Computational phage design for antibiotic-resistant bacteria (MEGA-27, item 8).

**Goal.** A phage-host interaction predictor for *Klebsiella pneumoniae* and
*Escherichia coli* phages: CNN encoders over phage receptor-binding proteins
(tail fibers / tail spikes) and GNN encoders over host receptor residue graphs,
trained on real interaction data mined from NCBI/GenBank, UniProt, PDB and the
AlphaFold DB. Deliverables: (1) a validated interaction model benchmarked
against published leaders (WIsH-, PHP-, DeepHost-class baselines), (2) a
structure-based docking screen of phage RBP vs host receptor, (3) a host-range
discovery screen over a resistant-strain panel (KPC/NDM Klebsiella, E. coli
ST131), (4) a hermetic pytest suite, and (5) a full paper with mathematical
derivations.

Program rules: real open datasets only; no stubs or simulated results; honest
negative results preserved and reported.
