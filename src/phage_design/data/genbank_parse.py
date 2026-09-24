"""Parse GenBank flat files for phage host annotations and RBP sequences.

Receptor-binding proteins (RBPs) are identified by product keywords used in
real GenBank annotations of tailed phages: tail fiber, tail spike, receptor
binding/binding protein, tail adhesin, baseplate wedge (gpJ-class). Host is
taken from the /host qualifier (preferred) or /isolation_source / organism
source fields.
"""
from __future__ import annotations

import io
import re
from dataclasses import dataclass, field

from Bio import SeqIO

RBP_KEYWORDS = re.compile(
    r"(tail[ -]?fib[er]{2}|tail[ -]?spike|receptor[ -]?binding|tail[ -]?adhesin|"
    r"long tail fiber|short tail fiber|tail fiber assembly)",
    re.IGNORECASE,
)


@dataclass
class PhageRecord:
    accession: str
    organism: str
    host: str | None
    genome_length: int
    rbp_products: list[tuple[str, str]] = field(default_factory=list)  # (product, aa_seq)


def parse_genbank(text: str) -> PhageRecord:
    rec = SeqIO.read(io.StringIO(text), "genbank")
    return _record_to_phage(rec)


def parse_genbank_multi(text: str):
    """Yield PhageRecord for each record in a concatenated multi-record fetch."""
    import io as _io
    for rec in SeqIO.parse(_io.StringIO(text), "genbank"):
        yield _record_to_phage(rec)


def _record_to_phage(rec) -> PhageRecord:
    host = None
    for feat in rec.features:
        if feat.type == "source":
            quals = feat.qualifiers
            host = (quals.get("host") or quals.get("isolation_source") or [None])[0]
    rbps = []
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        product = feat.qualifiers.get("product", [""])[0]
        if RBP_KEYWORDS.search(product):
            aa = feat.qualifiers.get("translation", [""])[0]
            if aa:
                rbps.append((product, aa))
    return PhageRecord(
        accession=rec.id.split(".")[0],
        organism=rec.annotations.get("organism", ""),
        host=host,
        genome_length=len(rec.seq),
        rbp_products=rbps,
    )
