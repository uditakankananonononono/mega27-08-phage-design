"""Build the honest external-tools registry and dataset manifest for the paper."""
import csv, json

TOOLS = [
 # tool, category, used for, status, artifact
 ("NCBI Nucleotide E-utilities", "database/API", "phage genome pulls (title queries, backoff)", "used", "data/raw/genbank/"),
 ("NCBI GenBank records", "database", "RBP CDS mining + host fields", "used", "data/processed/phage_rbp_table.csv"),
 ("NCBI Taxonomy (esummary)", "database/API", "species identity verification (taxid 573)", "used", "external/tools/ncbi_tax_573.json"),
 ("NCBI nuccore esearch (Virus scope)", "database/API", "corpus-size sanity counts", "used", "external/tools/ncbi_virus_kleb_count.json"),
 ("NCBI Datasets v2 API", "database/API", "candidate genome metadata probe (no assembly-level record for PZ accessions; recorded)", "queried, no record", "external/tools/ncbi_datasets_PZ797503.json"),
 ("NCBI Protein efetch", "database/API", "independent accession verification (OmpC P06996 = OMPC_ECOLI)", "used", "external/tools/ncbi_protein_P06996.fasta"),
 ("ENA Browser/Portal API", "database/API", "candidate verification (PZ797503 = E. phage vB_EscC_ArakU1, 60,473 bp)", "used", "external/tools/ena_PZ797503_report.json"),
 ("UniProt REST", "database/API", "14 receptor sequences + annotations (live accession resolution)", "used", "data/processed/host_receptor_panel.csv; external/tools/uniprot_P06996.json"),
 ("AlphaFold Protein Structure DB", "database", "14 receptor structure models", "used", "data/structures/"),
 ("AlphaFold API", "database/API", "model version resolution (v6), not pinned", "used", "scripts/fetch_receptor_structures.py"),
 ("RCSB PDB + Data API", "database/API", "docking control complex 8A8C metadata", "used", "external/tools/rcsb_8a8c.json"),
 ("PDBe API", "database/API", "8A8C summary cross-check", "used", "external/tools/pdbe_8a8c.json"),
 ("KEGG REST", "database/API", "receptor gene entry (eco:b0489 lamB)", "used", "external/tools/kegg_lamb.txt"),
 ("ICTV Master Species List", "database", "phage taxonomy reference", "used", "external/tools/ictv_msl_page.html"),
 ("Europe PMC", "literature API", "PHP citation verification", "used", "external/tools/europepmc_php.json"),
 ("CrossRef API", "literature API", "VirHostMatcher DOI metadata (10.1093/nar/gkx382)", "used", "external/tools/crossref_vhm.json"),
 ("PubMed E-utilities", "literature API", "T5 pb5-FhuA structure paper lookup", "used", "external/tools/pubmed_pb5.json"),
 ("InterPro API", "database/API", "receptor domain/family annotation (P06996)", "used", "external/tools/interpro_P06996.json"),
 ("PhagesDB API", "database/API", "host-genera coverage check (actinophage-centric; justifies building own corpus)", "used", "external/tools/phagesdb_phages.json"),
 ("GTDB API", "database/API", "K. pneumoniae taxonomy verification", "used", "external/tools/gtdb_kleb.json"),
 ("Wikidata SPARQL", "knowledge base", "gene-protein mapping probe (empty result for b0153 query; recorded)", "queried, empty", "external/tools/wikidata_btub.json"),
 ("OpenAlex API", "literature API", "WIsH citation lookup (rate-limited after first success)", "used + rate-limited", "external/tools/openalex_wish.json"),
 ("Virus-Host DB (genome.jp)", "database", "host-annotation cross-check (mirror 404 at pull time; retry pending)", "attempted, unavailable", "external/tools/"),
 ("BV-BRC API", "database/API", "Klebsiella phage count cross-check (query syntax rejected; recorded)", "attempted, rejected", "external/tools/bvbrc_kleb.json"),
 ("PHP (Lu et al. 2021)", "published tool", "head-to-head comparator, native + forced-binary", "used", "external/php/; results/php_headtohead.json"),
 ("PyTorch", "package", "CNN/GNN/hybrid model training", "used", "results/interaction_cnn_seed*.pt"),
 ("NumPy", "package", "numerical pipeline", "used", "throughout"),
 ("pandas", "package", "all table processing", "used", "throughout"),
 ("SciPy", "package", "Mann-Whitney tests, statistics", "used", "results/external_verification.json"),
 ("scikit-learn", "package", "k-mer LR baseline, group splits, metrics", "used", "results/interaction_cnn_seed*.json"),
 ("Biopython", "package", "GenBank parsing, ProtParam cross-check, pairwise identity", "used", "results/external_verification.json"),
 ("Matplotlib", "package", "figures 1-8", "used", "paper/figures/"),
 ("logomaker", "package", "fiber-tip motif logo (figure 9)", "used", "paper/figures/fig9_fiber_tip_logo.png"),
 ("statsmodels", "package", "statistical support for extended math appendix", "used", "results/"),
 ("duckdb", "package", "SQL analytics over result tables (screen cluster stats)", "used", "results/"),
 ("pytest", "package", "26 hermetic tests", "used", "tests/"),
 ("pandoc", "tool", "paper HTML build", "used", "paper/"),
 ("wkhtmltopdf", "tool", "paper PDF build", "used", "paper/paper.pdf"),
 ("EBI EB-eye search API", "search API", "cross-resource lookup for 8A8C", "used", "external/tools/ebeye_8a8c.json"),
 ("GitHub", "platform", "PHP source clone + project repository hosting", "used", "external/php/"),
 ("Python 3.10", "runtime", "all pipeline code", "used", "src/ scripts/ tests/"),
 ("pdflatex (TeX Live)", "tool", "PDF engine attempt (missing xcolor in sandbox; recorded, wkhtmltopdf used instead)", "attempted, unavailable", "/tmp/paper.pdf"),
]

DATASETS = [
 ("phage genomes, two-species model corpus (GenBank accessions)", 790),
 ("additional host-annotated Klebsiella spp. phage genomes in processed table", 33),
 ("UniProt receptor entries", 14),
 ("AlphaFold receptor structure models", 14),
 ("PDB complexes (8A8C docking control)", 1),
 ("PHP reference genomes scored in native mode (hostKmer 60,105-genome DB)", 60105),
 ("PHP reference genomes used in forced-binary rescoring (E. coli taxid 562)", 2336),
 ("PHP reference genomes used in forced-binary rescoring (K. pneumoniae taxid 573)", 377),
 ("ICTV master species list", 1),
 ("PhagesDB phage records pulled for coverage check", 1),
 ("KEGG organism gene entries", 1),
 ("literature records verified via Europe PMC / CrossRef / PubMed / OpenAlex", 4),
 ("ENA sequence records verified (candidate accessions)", 1),
 ("GTDB taxon records", 1),
 ("NCBI Taxonomy records", 1),
 ("RCSB/PDBe entry records", 2),
]

with open("results/external_tools_registry.json", "w") as fh:
    json.dump([{"tool": t, "category": c, "used_for": u, "status": s, "artifact": a} for t, c, u, s, a in TOOLS], fh, indent=1)
with open("results/dataset_manifest.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["dataset", "accession_level_count"]); w.writerows(DATASETS)

used = [t for t in TOOLS if t[3].startswith("used")]
total_ds = sum(n for _, n in DATASETS)
print(f"tools total={len(TOOLS)} actively-used={len(used)}  datasets accession-level total={total_ds}")
