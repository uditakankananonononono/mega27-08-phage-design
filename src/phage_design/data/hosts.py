"""Host normalization for phage host annotations.

GenBank /host qualifiers range from species ('Escherichia coli') to strain
('Klebsiella pneumoniae INF359') to non-strain free text ('hospital sewage').
We normalize to (species, strain) with explicit unknowns; non-host text maps
to (None, None) and is excluded from supervised pairs.
"""
from __future__ import annotations

KNOWN_SPECIES = [
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Klebsiella oxytoca",
    "Klebsiella aerogenes",
    "Klebsiella variicola",
    "Klebsiella michiganensis",
    "Salmonella enterica",
    "Shigella flexneri",
    "Shigella sonnei",
    "Enterobacter cloacae",
    "Citrobacter freundii",
    "Pseudomonas aeruginosa",
    "Acinetobacter baumannii",
    "Serratia marcescens",
    "Proteus mirabilis",
    "Cronobacter sakazakii",
    "Escherichia albertii",
]


def normalize_host(raw: str | None) -> tuple[str | None, str | None]:
    """Return (species, strain); strain may be None, species None if unmapped."""
    if raw is None or not isinstance(raw, str) or not raw.strip():
        return None, None
    text = " ".join(raw.strip().split())
    for sp in KNOWN_SPECIES:
        if text.lower().startswith(sp.lower()):
            strain = text[len(sp):].strip() or None
            return sp, strain
    return None, None
