# Capsule-lane retrospective audit, 2026-10-08

This legacy repository is distinct from the capsule-specific lane at https://github.com/uditakankananonononono/phage-capsule-rbp-kpneumoniae . An audit of that lane reads existing code and data only, without new predictions, sequence designs, nominations or clinical recommendations. This note does not verify or certify this legacy repository's models, paper, cocktails, genus-holdout claims or source rights.

At capsule-lane source commit `32786465a1b072f2e4caf9c40babf2b03a228db1`, the prototype archives 131 held-out-host diagnostic metrics: 103/131 hit@5, macro host AUROC 0.81132348. Arithmetic and fixed-result bootstrap intervals were independently checked. Its original pair scores and individual fold files are not committed, so original model execution was not rerun or certified.

The newer author-parameter reimplementation has a complete prediction archive covering 133 locus groups, 200 hosts, 105 phages, 10,006 pairs and 333 positives. Every label and group assignment was checked against freshly downloaded author files from https://zenodo.org/records/11061100 . Archived predictions yield 95/131 descriptive host hit@5, macro host AUROC 0.71530035 and pooled pair AUROC 0.75132307. It differs from the original author runtime and uses a different split from the prototype; no fair paired superiority or independent validation claim follows.

The capsule lane's 105-genome inventory matches deposited bytes. Its cargo summary reports a limited protein-only search with zero hits, but original reports, logs and ORF files are absent from the committed tree. Search execution cannot be independently certified from hashes alone. All 105 remain uncleared; no nomination exists. Viral lifecycle, genome closure, complete hazard coverage and a commensal off-target panel remain unverified.

Full scoped evidence, limitations and the executed audit-only checker belong in the capsule lane's `docs/RETROSPECTIVE-AUDIT-2026-10-08.md` and `audit/retrospective-check-2026-10-08.json`. Author software: https://github.com/dimiboeckaerts/PhageHostLearn . No design or host-range optimization is part of the continuation.
