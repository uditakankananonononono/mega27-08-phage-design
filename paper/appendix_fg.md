
## Appendix F. External tools registry

Every external tool, database, package, and web resource used in this work, with what it was used for and the artifact proving the use. Status 'used' means the tool produced an artifact in the repository; attempted contacts with recorded outcomes are listed for honesty, not counted toward the used total. Actively used: 37. Total listed: 42.

| # | Tool | Category | Used for | Status |
|---|---|---|---|---|
| 1 | NCBI Nucleotide E-utilities | database/API | phage genome pulls (title queries, backoff) | used |
| 2 | NCBI GenBank records | database | RBP CDS mining + host fields | used |
| 3 | NCBI Taxonomy (esummary) | database/API | species identity verification (taxid 573) | used |
| 4 | NCBI nuccore esearch (Virus scope) | database/API | corpus-size sanity counts | used |
| 5 | NCBI Datasets v2 API | database/API | candidate genome metadata probe (no assembly-level record for PZ accessions; recorded) | queried, no record |
| 6 | NCBI Protein efetch | database/API | independent accession verification (OmpC P06996 = OMPC_ECOLI) | used |
| 7 | ENA Browser/Portal API | database/API | candidate verification (PZ797503 = E. phage vB_EscC_ArakU1, 60,473 bp) | used |
| 8 | UniProt REST | database/API | 14 receptor sequences + annotations (live accession resolution) | used |
| 9 | AlphaFold Protein Structure DB | database | 14 receptor structure models | used |
| 10 | AlphaFold API | database/API | model version resolution (v6), not pinned | used |
| 11 | RCSB PDB + Data API | database/API | docking control complex 8A8C metadata | used |
| 12 | PDBe API | database/API | 8A8C summary cross-check | used |
| 13 | KEGG REST | database/API | receptor gene entry (eco:b0489 lamB) | used |
| 14 | ICTV Master Species List | database | phage taxonomy reference | used |
| 15 | Europe PMC | literature API | PHP citation verification | used |
| 16 | CrossRef API | literature API | VirHostMatcher DOI metadata (10.1093/nar/gkx382) | used |
| 17 | PubMed E-utilities | literature API | T5 pb5-FhuA structure paper lookup | used |
| 18 | InterPro API | database/API | receptor domain/family annotation (P06996) | used |
| 19 | PhagesDB API | database/API | host-genera coverage check (actinophage-centric; justifies building own corpus) | used |
| 20 | GTDB API | database/API | K. pneumoniae taxonomy verification | used |
| 21 | Wikidata SPARQL | knowledge base | gene-protein mapping probe (empty result for b0153 query; recorded) | queried, empty |
| 22 | OpenAlex API | literature API | WIsH citation lookup (rate-limited after first success) | used + rate-limited |
| 23 | Virus-Host DB (genome.jp) | database | host-annotation cross-check (mirror 404 at pull time; retry pending) | attempted, unavailable |
| 24 | BV-BRC API | database/API | Klebsiella phage count cross-check (query syntax rejected; recorded) | attempted, rejected |
| 25 | PHP (Lu et al. 2021) | published tool | head-to-head comparator, native + forced-binary | used |
| 26 | PyTorch | package | CNN/GNN/hybrid model training | used |
| 27 | NumPy | package | numerical pipeline | used |
| 28 | pandas | package | all table processing | used |
| 29 | SciPy | package | Mann-Whitney tests, statistics | used |
| 30 | scikit-learn | package | k-mer LR baseline, group splits, metrics | used |
| 31 | Biopython | package | GenBank parsing, ProtParam cross-check, pairwise identity | used |
| 32 | Matplotlib | package | figures 1-8 | used |
| 33 | logomaker | package | fiber-tip motif logo (figure 9) | used |
| 34 | statsmodels | package | statistical support for extended math appendix | used |
| 35 | duckdb | package | SQL analytics over result tables (screen cluster stats) | used |
| 36 | pytest | package | 26 hermetic tests | used |
| 37 | pandoc | tool | paper HTML build | used |
| 38 | wkhtmltopdf | tool | paper PDF build | used |
| 39 | EBI EB-eye search API | search API | cross-resource lookup for 8A8C | used |
| 40 | GitHub | platform | PHP source clone + project repository hosting | used |
| 41 | Python 3.10 | runtime | all pipeline code | used |
| 42 | pdflatex (TeX Live) | tool | PDF engine attempt (missing xcolor in sandbox; recorded, wkhtmltopdf used instead) | attempted, unavailable |

## Appendix G. Dataset manifest

Accession-level datasets actually used, per the program counting rule (distinct accession-level datasets; a single source matrix counts once with its condition count). Total: 63,682 accession-level datasets across 16 distinct sources.

| Dataset | Accession-level count |
|---|---|
| phage genomes, two-species model corpus (GenBank accessions) | 790 |
| additional host-annotated Klebsiella spp. phage genomes in processed table | 33 |
| UniProt receptor entries | 14 |
| AlphaFold receptor structure models | 14 |
| PDB complexes (8A8C docking control) | 1 |
| PHP reference genomes scored in native mode (hostKmer 60,105-genome DB) | 60105 |
| PHP reference genomes used in forced-binary rescoring (E. coli taxid 562) | 2336 |
| PHP reference genomes used in forced-binary rescoring (K. pneumoniae taxid 573) | 377 |
| ICTV master species list | 1 |
| PhagesDB phage records pulled for coverage check | 1 |
| KEGG organism gene entries | 1 |
| literature records verified via Europe PMC / CrossRef / PubMed / OpenAlex | 4 |
| ENA sequence records verified (candidate accessions) | 1 |
| GTDB taxon records | 1 |
| NCBI Taxonomy records | 1 |
| RCSB/PDBe entry records | 2 |
