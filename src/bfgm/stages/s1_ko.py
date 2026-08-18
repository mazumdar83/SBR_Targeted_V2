"""Stage 1: seed gene symbols -> KEGG KOs, with collision curation.

Exact token matching against the KO symbol field. Every multi-hit symbol is scored
against the term lexicon so the collision resolution is term-driven rather than
hardcoded, and every rejected KO is written to an audit trail.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd

from ..clients.kegg import KeggClient
from ..lexicon import TermLexicon

STATUS = {
    "UNIQUE": "KEGG_CONFIRMED (unique symbol match)",
    "CURATED": "KEGG_CONFIRMED (auto-disambiguated from symbol collision)",
    "AMBIGUOUS": "AMBIGUOUS (multiple on-term KOs; needs human decision)",
    "OFF_TERM": "OFF_TERM (only matches KOs unrelated to the term)",
    "NOT_IN_KO": "NOT_IN_KEGG_KO (no orthology entry exists)",
}


def run(seed_csv: str | Path, lexicon: TermLexicon, out_dir: str | Path,
        kegg: KeggClient | None = None) -> pd.DataFrame:
    kegg = kegg or KeggClient()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    seed = pd.read_csv(seed_csv)
    kegg.load_ko_catalogue()
    k2p, k2m = kegg.ko_to_pathway(), kegg.ko_to_module()
    pnames, mnames = kegg.pathway_names(), kegg.module_names()

    rows, rejected = [], []
    for _, r in seed.iterrows():
        gene = str(r["gene"]).strip()
        hits = kegg.match_symbol(gene)

        if not hits:
            rows.append(_row(r, gene, "", "NOT_IN_KO", kegg, k2p, k2m, pnames, mnames))
            continue

        scored = [(k, kegg.ko_definition(k), lexicon.classify_text(kegg.ko_definition(k)))
                  for k in hits]
        on = [s for s in scored if s[2] == "ON_TERM"]

        if len(hits) == 1:
            status = "UNIQUE" if on else "OFF_TERM"
            rows.append(_row(r, gene, hits[0], status, kegg, k2p, k2m, pnames, mnames))
        elif len(on) == 1:
            keep = on[0][0]
            rows.append(_row(r, gene, keep, "CURATED", kegg, k2p, k2m, pnames, mnames,
                             rejected_kos="; ".join(f"{k} ({d[:60]})"
                                                    for k, d, c in scored if k != keep)))
            for k, d, c in scored:
                if k != keep:
                    rejected.append({"gene": gene, "rejected_ko": k, "ko_definition": d,
                                     "reason": "off-term in a multi-hit symbol",
                                     "kept_ko": keep})
        elif len(on) > 1:
            rows.append(_row(r, gene, "; ".join(k for k, _, _ in on), "AMBIGUOUS",
                             kegg, k2p, k2m, pnames, mnames,
                             rejected_kos="; ".join(f"{k}" for k, d, c in scored if c != "ON_TERM")))
        else:
            rows.append(_row(r, gene, "; ".join(hits), "OFF_TERM", kegg, k2p, k2m,
                             pnames, mnames))
            for k, d, c in scored:
                rejected.append({"gene": gene, "rejected_ko": k, "ko_definition": d,
                                 "reason": "no on-term KO for this symbol", "kept_ko": ""})

    df = pd.DataFrame(rows)
    df.to_csv(out / "gene_ko_map.csv", index=False)
    pd.DataFrame(rejected).to_csv(out / "rejected_ko_collisions.csv", index=False)

    # bootstrap the lexicon from confirmed KO definitions
    confirmed = df[df.status_code.isin(["UNIQUE", "CURATED"])]["ko_definition"].tolist()
    added = lexicon.expand_from_ko_definitions(confirmed)
    lexicon.save(out / "term_lexicon.expanded.json")

    summary = df.status_code.value_counts().to_dict()
    summary["lexicon_terms_added"] = len(added)
    pd.Series(summary).to_csv(out / "stage1_summary.csv")
    return df


def _row(r, gene, ko, status, kegg, k2p, k2m, pnames, mnames, rejected_kos=""):
    first = ko.split(";")[0].strip() if ko else ""
    paths = k2p.get(first, [])
    mods = k2m.get(first, [])
    return {
        "gene": gene,
        "system": r.get("system", ""),
        "organism_characterised": r.get("organism_characterised", ""),
        "evidence_tier": r.get("evidence_tier", ""),
        "critic_verdict": r.get("critic_verdict", ""),
        "kegg_ko": ko if ko else "NONE",
        "ko_symbols": ", ".join(kegg.load_ko_catalogue()[0].get(first, [])),
        "ko_definition": kegg.ko_definition(first),
        "ec_number": kegg.ec_for(first),
        "kegg_pathways": "; ".join(f"{p} {pnames.get(p,'')}" for p in sorted(paths)),
        "kegg_modules": "; ".join(f"{m} {mnames.get(m,'')}" for m in sorted(mods)),
        "status_code": status,
        "evidence_status": STATUS[status],
        "rejected_ko_collisions": rejected_kos,
        "pmids": r.get("pmids", ""),
    }
