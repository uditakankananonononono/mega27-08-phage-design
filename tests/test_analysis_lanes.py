"""Hermetic checks that committed analysis-lane artifacts exist and are sane."""
import json, os
import pandas as pd

R = "results"

def test_accession_evidence():
    df = pd.read_csv(os.path.join(R, "genbank_accession_evidence.tsv"), sep="\t", header=None)
    assert len(df) >= 790

def test_lane_jsons_present():
    for f in ["pyrodigal_annotation_qc.json", "fastani_species_check.json",
              "skbio_rbp_permanova.json", "networkx_corpus_graph.json",
              "sourmash_genome_containment.json", "pyhmmer_fiber_scan.json",
              "statsmodels_rbp_dose_response.json", "viennarna_rbs.json",
              "diamond_rbp_receptor.json", "structure_qc.json"]:
        assert os.path.exists(os.path.join(R, f)), f

def test_fastani_signal():
    d = json.load(open(os.path.join(R, "fastani_species_check.json")))
    assert d["same_host_genus_mean_ani"] > d["diff_host_genus_mean_ani"]

def test_pyhmmer_hits():
    df = pd.read_csv(os.path.join(R, "pyhmmer_fiber_scan.csv"))
    assert len(df) > 100

def test_registry_size():
    reg = json.load(open(os.path.join(R, "external_tools_registry.json")))
    assert len(reg) >= 40
