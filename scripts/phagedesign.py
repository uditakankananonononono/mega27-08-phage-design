#!/usr/bin/env python3
"""phagedesign - phage host-range prediction and cocktail reporting CLI.

Usage:
  phagedesign predict --genome phage.gb [--checkpoint ...] [--panel ...] [--norm ...]
  phagedesign screen  --genome phage.gb   # per-receptor scores (discovery screen)
  phagedesign cocktail [--design PD8-EC1]  # committed cocktail design + falsification matrix

Offline: uses only committed artifacts (trained checkpoint, receptor panel,
k-mer norm stats, cocktail design). No network access.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from phage_design.predict import (load_model, load_norm, load_panel,
                                  predict_host_range, rbps_from_genbank,
                                  receptor_screen)

DEF_CKPT = ROOT / "results" / "interaction_cnn.pt"
DEF_PANEL = ROOT / "data" / "processed" / "host_receptor_panel.csv"
DEF_NORM = ROOT / "results" / "kmer_norm.json"
DEF_COCKTAIL = ROOT / "results" / "cocktail_design.json"
DEF_MEMBERS = ROOT / "results" / "cocktail_members.csv"


def _load_rbps(args):
    rbps = [aa for _, aa in rbps_from_genbank(args.genome)]
    print(f"[phagedesign] {len(rbps)} receptor-binding protein(s) from {args.genome}",
          file=sys.stderr)
    return rbps


def cmd_predict(args) -> int:
    model = load_model(args.checkpoint)
    panel = load_panel(args.panel)
    mu, sd = load_norm(args.norm)
    rbps = _load_rbps(args)
    out = predict_host_range(model, panel, rbps, mu, sd)
    out["genome"] = str(args.genome)
    print(json.dumps(out, indent=2))
    return 0


def cmd_screen(args) -> int:
    model = load_model(args.checkpoint)
    panel = load_panel(args.panel)
    mu, sd = load_norm(args.norm)
    rbps = _load_rbps(args)
    rows = receptor_screen(model, panel, rbps, mu, sd)
    print(json.dumps({"genome": str(args.genome), "receptor_scores": rows}, indent=2))
    return 0


def cmd_cocktail(args) -> int:
    design = json.loads(Path(args.design_json).read_text())
    members = list(csv.DictReader(open(args.members_csv)))
    out = {"designs": {}}
    for spkey, blk in design["species"].items():
        cid = blk["cocktail_id"]
        if args.design and cid != args.design:
            continue
        leg = blk["leg_a_argmax_design"]
        ms = [m for m in members if m["cocktail_id"] == cid]
        out["designs"][cid] = {
            "strain_panel": blk["strain_panel"],
            "k_redundancy": leg["k_redundancy"],
            "targeted_receptors": leg["targeted_receptors"],
            "escape_union_size": leg["escape_union_size"],
            "escape_probability": leg["escape_probability"],
            "members": [{"accession": m["accession"], "phage": m["phage_name"],
                         "argmax_receptor": m["argmax_receptor"],
                         "discovery_candidate": m["discovery_candidate"]} for m in ms],
        }
    print(json.dumps(out, indent=2))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="phagedesign", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in [("predict", cmd_predict), ("screen", cmd_screen)]:
        p = sub.add_parser(name)
        p.add_argument("--genome", required=True, help="GenBank file with RBP annotations")
        p.add_argument("--checkpoint", default=str(DEF_CKPT))
        p.add_argument("--panel", default=str(DEF_PANEL))
        p.add_argument("--norm", default=str(DEF_NORM))
        p.set_defaults(fn=fn)
    p = sub.add_parser("cocktail")
    p.add_argument("--design", default=None, help="design name (default: all)")
    p.add_argument("--design-json", default=str(DEF_COCKTAIL))
    p.add_argument("--members-csv", default=str(DEF_MEMBERS))
    p.set_defaults(fn=cmd_cocktail)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
