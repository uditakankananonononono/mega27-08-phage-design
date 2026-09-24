"""Hermetic tests for the phagedesign CLI / predict library (random weights, tmp artifacts)."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from phage_design.models.cnn import ProteinCNNEncoder
from phage_design.models.grouped import GroupedInteractionModel
from phage_design.predict import (EMB, load_model, load_norm, load_panel,
                                  predict_host_range, rbps_from_genbank,
                                  receptor_screen)
import phagedesign

FIXTURE_GB = ROOT / "tests" / "fixtures" / "mini_phage.gb"


@pytest.fixture()
def artifacts(tmp_path):
    torch.manual_seed(0)
    model = GroupedInteractionModel(
        ProteinCNNEncoder(23, 48, EMB), ProteinCNNEncoder(23, 48, EMB), emb_dim=EMB)
    ckpt = tmp_path / "ckpt.pt"
    torch.save(model.state_dict(), ckpt)
    panel = tmp_path / "panel.csv"
    panel.write_text("species,receptor,accession,length,sequence\n"
                     "Species A,RecA,X1,32,ACDEFGHIKLMNPQRSTVWYACDEFGHIKL\n"
                     "Species A,RecB,X2,32,MNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n"
                     "Species B,RecC,Y1,32,GGGGACDEFHIKLMNPQRSTVWYAAAAAAA\n")
    norm = tmp_path / "norm.json"
    norm.write_text(json.dumps({"mu": [0.0] * 800, "sd": [1.0] * 800}))
    return ckpt, panel, norm


def test_predict_structure_and_consistency(artifacts):
    ckpt, panel, norm = artifacts
    model = load_model(ckpt)
    pan = load_panel(panel)
    mu, sd = load_norm(norm)
    rbps = [aa for _, aa in rbps_from_genbank(FIXTURE_GB)]
    out = predict_host_range(model, pan, rbps, mu, sd)
    assert set(out["probabilities"]) == {"Species A", "Species B"}
    assert all(0.0 <= p <= 1.0 for p in out["probabilities"].values())
    # ranking consistent with logits; predicted_host is the argmax
    by_logit = sorted(out["logits"], key=out["logits"].get, reverse=True)
    assert out["ranking"] == by_logit
    assert out["predicted_host"] == by_logit[0]


def test_predict_deterministic(artifacts):
    ckpt, panel, norm = artifacts
    model = load_model(ckpt)
    pan = load_panel(panel)
    mu, sd = load_norm(norm)
    rbps = [aa for _, aa in rbps_from_genbank(FIXTURE_GB)]
    a = predict_host_range(model, pan, rbps, mu, sd)
    b = predict_host_range(model, pan, rbps, mu, sd)
    assert a == b


def test_receptor_screen_rows_sorted(artifacts):
    ckpt, panel, norm = artifacts
    model = load_model(ckpt)
    pan = load_panel(panel)
    mu, sd = load_norm(norm)
    rbps = [aa for _, aa in rbps_from_genbank(FIXTURE_GB)]
    rows = receptor_screen(model, pan, rbps, mu, sd)
    assert len(rows) == 3  # one row per receptor
    logits = [r["logit"] for r in rows]
    assert logits == sorted(logits, reverse=True)


def test_rbps_from_genbank_raises_without_rbp(tmp_path):
    text = FIXTURE_GB.read_text().replace("long tail fiber protein", "hypothetical protein")
    gb = tmp_path / "norbp.gb"
    gb.write_text(text)
    with pytest.raises(ValueError, match="no receptor-binding proteins"):
        rbps_from_genbank(gb)


def test_cli_predict_end_to_end(artifacts, capsys):
    ckpt, panel, norm = artifacts
    rc = phagedesign.main(["predict", "--genome", str(FIXTURE_GB),
                           "--checkpoint", str(ckpt), "--panel", str(panel),
                           "--norm", str(norm)])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["predicted_host"] in ("Species A", "Species B")
    assert out["genome"] == str(FIXTURE_GB)


def test_cli_cocktail_reads_committed_format(tmp_path, capsys):
    dj = tmp_path / "d.json"
    dj.write_text(json.dumps({"parameters": {}, "species": {"ecoli": {
        "cocktail_id": "T-1", "strain_panel": "panel X",
        "leg_a_argmax_design": {"k_redundancy": 2, "targeted_receptors": ["R1"],
                                "escape_union_size": 1,
                                "escape_probability": {"mu=1e-06": {"union_bound": 1e-20}},
                                "members": []}}}}))
    mc = tmp_path / "m.csv"
    mc.write_text("cocktail_id,member_id,accession,phage_name,annotated_host,"
                  "discovery_candidate,p_species,argmax_receptor,argmax_score,in_heldout_test\n"
                  "T-1,T-1.1,ACC1,phage one,Host A,False,1.0,R1,1.0,False\n")
    rc = phagedesign.main(["cocktail", "--design-json", str(dj), "--members-csv", str(mc)])
    assert rc == 0
    out = json.loads(capsys.readouterr().out)
    assert out["designs"]["T-1"]["members"][0]["accession"] == "ACC1"
