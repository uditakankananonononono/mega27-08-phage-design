"""Query external tools/databases used for research & verification.
Each entry records: tool, what it was used for, artifact path. Honest registry."""
import json, time, urllib.parse, urllib.request

def get(url, path, timeout=25):
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "mega27-item8-research/1.0"})
        with urllib.request.urlopen(req, timeout=timeout) as r:
            data = r.read()
        open(path, "wb").write(data)
        return len(data)
    except Exception as e:
        return f"ERROR: {e}"

R = {}
# 1. RCSB PDB data API: verify 8A8C identity/method/resolution
R["rcsb_pdb"] = get("https://data.rcsb.org/rest/v1/core/entry/8A8C", "external/tools/rcsb_8a8c.json")
# 2. PDBe API: 8A8C summary
R["pdbe"] = get("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/8A8C", "external/tools/pdbe_8a8c.json")
# 3. Europe PMC: verify PHP reference (Wang 2021)
q = urllib.parse.quote('PHP phage host prediction Wang 2021')
R["europepmc_php"] = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={q}&format=json&pageSize=3", "external/tools/europepmc_php.json")
# 4. CrossRef: VirHostMatcher DOI metadata
R["crossref_vhm"] = get("https://api.crossref.org/works/10.1093/nar/gkx382", "external/tools/crossref_vhm.json")
# 5. OpenAlex: WIsH reference lookup
R["openalex_wish"] = get("https://api.openalex.org/works?search=WIsH%20phage%20host%20prediction%20Galiez&per-page=3", "external/tools/openalex_wish.json")
# 6. KEGG REST: E. coli lamB entry
R["kegg_lamb"] = get("https://rest.kegg.jp/get/eco:b0489", "external/tools/kegg_lamb.txt")
# 7. NCBI Taxonomy: K. pneumoniae taxid record
R["ncbi_tax"] = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=573&retmode=json", "external/tools/ncbi_tax_573.json")
# 8. ENA portal API: candidate accession PZ797503 metadata
R["ena_pz797503"] = get("https://www.ebi.ac.uk/ena/browser/api/json/PZ797503", "external/tools/ena_PZ797503.json")
# 9. NCBI Datasets v2: candidate accession metadata
R["ncbi_datasets"] = get("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/PZ797503/dataset_report", "external/tools/ncbi_datasets_PZ797503.json")
# 10. NCBI Virus: count of Klebsiella phages (corpus sanity)
q = urllib.parse.quote('"Klebsiella phage"[Organism]')
R["ncbi_virus"] = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={q}&retmode=json", "external/tools/ncbi_virus_kleb_count.json")
# 11. Wikidata SPARQL: BtuB gene/protein mapping
sq = urllib.parse.quote("SELECT ?item ?itemLabel WHERE { ?item wdt:P353 'b0153'. SERVICE wikibase:label { bd:serviceParam wikibase:language 'en'. } }")
R["wikidata_btub"] = get(f"https://query.wikidata.org/sparql?query={sq}&format=json", "external/tools/wikidata_btub.json")
# 12. Virus-Host DB (KEGG/genome.jp): download virus-host pairs TSV
R["viralhostdb"] = get("https://www.genome.jp/ftp/db/viralhostdb/viralhostdb.tsv", "external/tools/viralhostdb.tsv", timeout=60)
# 13. UniProt REST: OmpC entry (function annotation verify)
R["uniprot_ompc"] = get("https://rest.uniprot.org/uniprotkb/P06996.json", "external/tools/uniprot_P06996.json")
# 14. PubMed: pb5 FhuA T5 structure paper
q = urllib.parse.quote("T5 phage pb5 FhuA structure")
R["pubmed_pb5"] = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={q}&retmode=json&retmax=3", "external/tools/pubmed_pb5.json")
# 15. ICTV VMR (current virus taxonomy)
R["ictv_vmr"] = get("https://ictv.global/msl/current", "external/tools/ictv_msl_page.html", timeout=45)

print(json.dumps(R, indent=1))

# --- batch 2 ---
# 16. ENA filereport for candidate accession
R["ena_filereport"] = get("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PZ797503&result=sequence&fields=accession,description,base_count&format=json", "external/tools/ena_PZ797503_report.json")
# 17. Wikidata retry (SPARQL may 502 transiently)
R["wikidata_retry"] = get(f"https://query.wikidata.org/sparql?query={sq}&format=json", "external/tools/wikidata_btub.json")
# 18. Virus-Host DB alternate URL
R["viralhostdb2"] = get("https://www.genome.jp/ftp/db/viralhostdb/viralhostdb.tab", "external/tools/viralhostdb.tab", timeout=45)
# 19. InterPro API: OmpC entry annotations
R["interpro_ompc"] = get("https://www.ebi.ac.uk/interpro/api/protein/UniProt/P06996/", "external/tools/interpro_P06996.json")
# 20. PhagesDB API: host genera coverage
R["phagesdb"] = get("https://phagesdb.org/api/hostgenera/", "external/tools/phagesdb_hostgenera.json")
# 21. BV-BRC data API: Klebsiella phage genome count
q = urllib.parse.quote('eq(species,"Klebsiella pneumoniae phage")')
R["bvbrc"] = get(f"https://www.bv-brc.org/api/genome/?{q}&limit=1", "external/tools/bvbrc_kleb_phage.json")
# 22. GTDB API: K. pneumoniae taxonomy verify
R["gtdb"] = get("https://gtdb.ecogenomic.org/api/v1/taxonomy/species/Klebsiella%20pneumoniae", "external/tools/gtdb_kleb.json")
with open("external/tools/registry_status.json", "w") as fh:
    json.dump(R, fh, indent=1)
print(json.dumps(R, indent=1))
