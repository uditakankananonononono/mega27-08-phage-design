"""Write results/kmer_norm.json: k-mer feature mean/sd over the full model corpus.

Inference-time z-scoring stats (predict.py loads them). Full-corpus stats are
used (all 1,580 pair rows) so the artifact is seed-independent; training used
train-split stats - the difference is negligible and documented in the paper.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from phage_design.features.sequence import kmer_spectrum
from train_interaction import load_task

_, _, samples, _, species_receptors, host_names, rbp_lists = load_task()
panel_spec = {sp: np.mean([kmer_spectrum(x, k=2) for x in seqs], axis=0)
              for sp, seqs in species_receptors.items()}
feats = np.stack([
    np.concatenate([np.mean([kmer_spectrum(x, k=2) for x in rbps], axis=0) * panel_spec[sp],
                    np.mean([kmer_spectrum(x, k=2) for x in rbps], axis=0) - panel_spec[sp]])
    for rbps, sp in zip(rbp_lists, host_names)])
mu, sd = feats.mean(0), feats.std(0) + 1e-8
Path("results/kmer_norm.json").write_text(json.dumps(
    {"mu": mu.tolist(), "sd": sd.tolist(), "n_samples": len(samples)}))
print(f"wrote results/kmer_norm.json dim={len(mu)} n={len(samples)}")
