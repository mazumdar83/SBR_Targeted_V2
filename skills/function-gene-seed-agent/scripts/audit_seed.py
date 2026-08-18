#!/usr/bin/env python3
"""Four-check hallucination audit for seed_genes.csv.

  1. PMID exists            -> quarantine PMID_UNRESOLVED
  2. gene symbol resolves   -> quarantine SYMBOL_UNRESOLVED
  3. symbol collision       -> flag COLLISION_RISK (not a rejection)
  4. organism is real       -> flag ORGANISM_FLAG / quarantine ORGANISM_UNRESOLVED

Usage:
    python audit_seed.py seed_genes.csv term_lexicon.json --out .
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from bfgm.clients.kegg import KeggClient          # noqa: E402
from bfgm.clients.literature import LiteratureClient  # noqa: E402
from bfgm.clients.uniprot import UniProtClient    # noqa: E402
from bfgm.lexicon import TermLexicon              # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("seed")
    ap.add_argument("lexicon")
    ap.add_argument("--out", default=".")
    ap.add_argument("--skip-taxa", action="store_true",
                    help="skip check 4 (slow: one E-utilities call per organism)")
    a = ap.parse_args()

    out = Path(a.out)
    df = pd.read_csv(a.seed)
    lex = TermLexicon.load(a.lexicon)
    lit, kegg, up = LiteratureClient(), KeggClient(), UniProtClient()

    # check 1
    all_pmids = {p for v in df.pmids.fillna("") for p in str(v).split(";") if p.strip()}
    good = lit.verify_pmids(all_pmids)
    bad_pmids = all_pmids - good
    print(f"[1] PMIDs: {len(good)}/{len(all_pmids)} resolve")

    # check 2 and 3
    kegg.load_ko_catalogue()
    status, reasons = [], []
    for _, r in df.iterrows():
        gene = str(r["gene"]).strip()
        cited = [p for p in str(r.get("pmids", "")).split(";") if p.strip()]
        if cited and all(p in bad_pmids for p in cited):
            status.append("QUARANTINE"); reasons.append("PMID_UNRESOLVED"); continue

        kos = kegg.match_symbol(gene)
        in_uniprot = up.count(f"gene_exact:{gene}") > 0 if not kos else True
        if not kos and not in_uniprot:
            status.append("QUARANTINE"); reasons.append("SYMBOL_UNRESOLVED"); continue

        defs = [kegg.ko_definition(k) for k in kos]
        on = [d for d in defs if lex.classify_text(d) == "ON_TERM"]
        if len(kos) > 1 and len(on) != 1:
            status.append("COLLISION_RISK")
            reasons.append(f"{len(kos)} KO hits, {len(on)} on-term")
        else:
            status.append("PASS"); reasons.append("")
    df["audit_status"] = status
    df["audit_note"] = reasons

    # check 4
    if not a.skip_taxa:
        orgs = [o for o in df.organism_characterised.fillna("").unique() if o.strip()]
        taxa = lit.verify_taxa(orgs)
        unresolved = {o for o, t in taxa.items() if not t}
        df.loc[df.organism_characterised.isin(unresolved) & (df.audit_status == "PASS"),
               "audit_status"] = "ORGANISM_FLAG"
        print(f"[4] organisms: {len(orgs) - len(unresolved)}/{len(orgs)} resolve")

    q = df[df.audit_status == "QUARANTINE"].copy()
    q["quarantine_reason"] = q.audit_note
    keep = df[df.audit_status != "QUARANTINE"]

    keep.to_csv(out / "seed_genes.audited.csv", index=False)
    q.to_csv(out / "quarantine.csv", index=False)
    print(f"\nkept {len(keep)}, quarantined {len(q)}")
    print(keep.audit_status.value_counts().to_string())
    print("\nQuarantine is a deliverable. Do not delete it.")
    print("Next: microbiologist critic review per references/critic-persona.md")


if __name__ == "__main__":
    main()
