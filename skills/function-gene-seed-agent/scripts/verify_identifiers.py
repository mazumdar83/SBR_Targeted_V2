#!/usr/bin/env python3
"""Standalone identifier resolver: PMIDs, gene symbols, organism names."""
import argparse, sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from bfgm.clients.kegg import KeggClient
from bfgm.clients.literature import LiteratureClient
from bfgm.clients.uniprot import UniProtClient

ap = argparse.ArgumentParser()
ap.add_argument("--pmids", nargs="*", default=[])
ap.add_argument("--genes", nargs="*", default=[])
ap.add_argument("--organisms", nargs="*", default=[])
a = ap.parse_args()
res = {}
if a.pmids:
    good = LiteratureClient().verify_pmids(a.pmids)
    res["pmids"] = {p: ("RESOLVES" if p in good else "NOT_FOUND") for p in a.pmids}
if a.genes:
    k, u = KeggClient(), UniProtClient()
    k.load_ko_catalogue()
    res["genes"] = {g: {"kegg_kos": k.match_symbol(g),
                        "uniprot_bacterial_reviewed": u.count(f"gene_exact:{g}")}
                    for g in a.genes}
if a.organisms:
    res["organisms"] = LiteratureClient().verify_taxa(a.organisms)
print(json.dumps(res, indent=1))
