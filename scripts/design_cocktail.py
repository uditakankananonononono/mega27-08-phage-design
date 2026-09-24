"""Design named, falsifiable phage cocktails from the discovery-screen
per-receptor scores.

Two legs, both fully reported:

Leg A (PRIMARY, PD8-EC1 / PD8-KP1): argmax-usage k-redundant design.
A phage's usage claim is its top-scoring receptor (argmax) -- the only
usage claim with direct evidentiary support (per-receptor discrimination
analysis, results/per_receptor_analysis.json). Targeted receptors must
(a) discriminate own- vs other-species phages with AUC >= 0.9 and
(b) have >= k pool phages with that argmax. Members: top-k pool
phages per receptor ranked by p_species * argmax_score, with the last seat
per receptor reserved for the best discovery candidate (host-switch claim)
when one exists -- those seats are direct tests of the discovery claims. Cocktail escape depth
|U| = number of targeted receptors; depth k per receptor guards against
per-member model error (P(all k claims on r wrong) = eps^k).

Leg B (SENSITIVITY): absolute-threshold multicover (score >= theta) at
theta in {0.5, 0.7, 0.9} via greedy k-cover. Reported to show why
absolute thresholds are degenerate for this model (supports of size 6-10
inconsistent with the discrimination analysis).

Outputs:
  results/cocktail_design.json               full design record
  results/cocktail_members.csv               leg-A member table
  results/cocktail_falsification_matrix.csv  leg-A plaque predictions
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from phage_design import cocktail as ck

SPECIES_CFG = {
    "ecoli": {
        "cocktail_id": "PD8-EC1",
        "p_col": "p_ecoli",
        "top_col": "top_receptor_ecoli",
        "prefix": "E-",
        "strain": "Escherichia coli K-12 MG1655 (WT) and single-knockout panel",
    },
    "klebsiella": {
        "cocktail_id": "PD8-KP1",
        "p_col": "p_klebsiella",
        "top_col": "top_receptor_klebsiella",
        "prefix": "K-",
        "strain": "Klebsiella pneumoniae KPN2146 (WT) and single-knockout panel",
    },
}
P_MIN, K_RED, AUC_MIN = 0.9, 3, 0.9
THETAS = [0.5, 0.7, 0.9]
MU_GRID = [1e-7, 1e-6, 1e-5]
EPS_GRID = [0.05, 0.1, 0.2]
GENERATIONS, POPULATION = 20, 1e8


def escape_block(u, depth_vec, receptors):
    return {
        "escape_union": [r for r, d in zip(receptors, depth_vec) if d > 0],
        "escape_union_size": int(u),
        "escape_probability": {
            f"mu={mu:g}": ck.exact_escape_probability(u, mu, GENERATIONS, POPULATION)
            for mu in MU_GRID},
        "escape_model": {"generations": GENERATIONS, "population": POPULATION,
                         "assumption": "independent per-receptor LOF, no fitness cost"},
        "expected_true_escape_depth": {
            f"eps={eps:g}": ck.expected_true_escape_depth(np.asarray(depth_vec), eps)
            for eps in EPS_GRID},
        "single_class_escape": {
            f"mu={mu:g}": ck.exact_escape_probability(1, mu, GENERATIONS, POPULATION)
            for mu in MU_GRID},
    }


def leg_a(df, cfg, receptor_auc):
    rec_cols = [c for c in df.columns if c.startswith(cfg["prefix"])]
    receptors = [c[len(cfg["prefix"]):] for c in rec_cols]
    pool = df[df[cfg["p_col"]] >= P_MIN].reset_index(drop=True)
    argmax_rec = pool[cfg["top_col"]].str.split(":").str[0]
    argmax_score = pool[cfg["top_col"]].str.split(":").str[1].astype(float)
    rank_score = pool[cfg["p_col"]] * argmax_score

    validated = [r for r in receptors if receptor_auc.get(cfg["prefix"] + r, 0) >= AUC_MIN]
    excluded_auc = [r for r in receptors if r not in validated]
    targeted, members, depth_map = [], [], {}
    for r in validated:
        mask = (argmax_rec == r)
        n_avail = int(mask.sum())
        if n_avail < K_RED:
            depth_map[r] = n_avail
            continue
        # Seat rule: top-k by rank score, with the LAST seat reserved for the
        # best discovery candidate on this receptor when one exists -- each
        # such seat is a falsifiable test of the host-switch claim.
        idx_all = pool.index[mask][np.argsort(-rank_score[mask])]
        disc_idx = [i for i in idx_all
                    if pool.loc[i, "call"] == "candidate-host-switch"]
        non_idx = [i for i in idx_all
                   if pool.loc[i, "call"] != "candidate-host-switch"]
        if disc_idx:
            idx = list(non_idx[:K_RED - 1]) + [disc_idx[0]]
        else:
            idx = list(idx_all[:K_RED])
        targeted.append(r)
        depth_map[r] = K_RED
        for i in idx:
            members.append({
                "accession": pool.loc[i, "accession"],
                "phage_name": pool.loc[i, "phage_name"],
                "annotated_host": pool.loc[i, "annotated_host"],
                "discovery_candidate": bool(pool.loc[i, "call"] == "candidate-host-switch"),
                "p_species": float(pool.loc[i, cfg["p_col"]]),
                "argmax_receptor": r,
                "argmax_score": float(argmax_score[i]),
                "rank_score": float(rank_score[i]),
                "in_heldout_test": bool(pool.loc[i, "in_heldout_test"]),
            })
    for j, m in enumerate(members):
        m["member_id"] = f"{cfg['cocktail_id']}.{j+1}"
        m["receptor_support"] = [m["argmax_receptor"]]
        m["support_size"] = 1
    infeasible = [r for r in validated if r not in targeted]

    u = len(targeted)
    depth_vec = [depth_map.get(r, 0) for r in receptors]
    out = {"pool_size": int(len(pool)), "k_redundancy": K_RED,
           "validated_receptors": validated,
           "excluded_receptors_low_auc": excluded_auc,
           "infeasible_receptors": infeasible,
           "targeted_receptors": targeted,
           "depth": depth_map, "members": members,
           **escape_block(u, depth_vec, receptors)}
    # falsifiable matrix: member x {WT, single KO of each targeted receptor}
    fals = []
    for m in members:
        preds = {"WT": 1}
        for r in targeted:
            preds[f"delta_{r}"] = 0 if m["argmax_receptor"] == r else 1
        fals.append({"member_id": m["member_id"], **preds})
    out["falsification_matrix"] = fals
    out["weakest_member_predictions"] = [
        {"member_id": m["member_id"], "phage_name": m["phage_name"],
         "prediction": f"host knockout of {m['argmax_receptor']} alone escapes this member"}
        for m in members]
    return out, fals


def leg_b(df, cfg):
    rec_cols = [c for c in df.columns if c.startswith(cfg["prefix"])]
    receptors = [c[len(cfg["prefix"]):] for c in rec_cols]
    pool = df[df[cfg["p_col"]] >= P_MIN].reset_index(drop=True)
    X = pool[rec_cols].to_numpy(dtype=float)
    per_theta = {}
    for theta in THETAS:
        U = ck.usage_matrix(X, theta)
        sel = ck.greedy_k_cover(U, K_RED)
        U_sel = U[pool.index[sel]] if sel else np.zeros((0, len(receptors)), bool)
        depth = ck.coverage_depth(U_sel) if sel else np.zeros(len(receptors), int)
        union = ck.escape_union(U_sel) if sel else np.zeros(len(receptors), bool)
        per_theta[f"theta={theta}"] = {
            "n_members": len(sel),
            "feasible_k_cover": bool((depth >= K_RED).all()),
            "support_sizes": ck.member_support_sizes(U_sel) if sel else [],
            "depth": {r: int(d) for r, d in zip(receptors, depth)},
            "escape_union_size": int(union.sum()),
            "member_accessions": pool.loc[pool.index[sel], "accession"].tolist(),
        }
    return per_theta


def main():
    screen = pd.read_csv("results/discovery_screen.csv")
    scores = pd.read_csv("results/per_receptor_scores.csv")
    names = pd.read_csv("data/processed/phage_rbp_table.csv",
                        usecols=["accession", "organism"])
    auc = json.load(open("results/per_receptor_analysis.json"))
    receptor_auc = {k: v["auc"] for k, v in auc.items()}
    name_of = dict(zip(names.accession, names.organism))
    df = screen.merge(scores, on=["accession", "annotated_host"], how="inner",
                      validate="one_to_one")
    df["phage_name"] = df.accession.map(name_of).fillna("")
    assert len(df) == len(screen)

    design = {"parameters": {"p_min": P_MIN, "k_redundancy": K_RED,
                             "auc_min": AUC_MIN, "mu_grid": MU_GRID,
                             "eps_grid": EPS_GRID, "generations": GENERATIONS,
                             "population": POPULATION},
              "species": {}}
    member_rows, fals_rows = [], []
    for sp_key, cfg in SPECIES_CFG.items():
        a, fals = leg_a(df, cfg, receptor_auc)
        b = leg_b(df, cfg)
        design["species"][sp_key] = {
            "cocktail_id": cfg["cocktail_id"], "strain_panel": cfg["strain"],
            "leg_a_argmax_design": a, "leg_b_threshold_sensitivity": b}
        for m in a["members"]:
            member_rows.append({"cocktail_id": cfg["cocktail_id"], **{
                k2: m[k2] for k2 in ["member_id", "accession", "phage_name",
                                     "annotated_host", "discovery_candidate",
                                     "p_species", "argmax_receptor",
                                     "argmax_score", "in_heldout_test"]}})
        for fr in fals:
            fals_rows.append({"cocktail_id": cfg["cocktail_id"], **fr})
    json.dump(design, open("results/cocktail_design.json", "w"), indent=2)
    pd.DataFrame(member_rows).to_csv("results/cocktail_members.csv", index=False)
    pd.DataFrame(fals_rows).to_csv("results/cocktail_falsification_matrix.csv",
                                   index=False)
    for sp_key, blk in design["species"].items():
        a = blk["leg_a_argmax_design"]
        print(f"{blk['cocktail_id']}: {len(a['members'])} members over "
              f"{a['escape_union_size']} targeted receptors {a['targeted_receptors']}; "
              f"excluded(low AUC)={a['excluded_receptors_low_auc']}; "
              f"infeasible={a['infeasible_receptors']}")
        ndisc = sum(m['discovery_candidate'] for m in a['members'])
        print(f"  discovery candidates included: {ndisc}")
        for m in a["members"]:
            print(f"  {m['member_id']} {m['phage_name']} ({m['accession']}) "
                  f"-> {m['argmax_receptor']} score={m['argmax_score']:.3f} "
                  f"p={m['p_species']:.4f}{' DISCOVERY' if m['discovery_candidate'] else ''}")
    print("wrote results/cocktail_design.json, cocktail_members.csv, "
          "cocktail_falsification_matrix.csv")


if __name__ == "__main__":
    main()
