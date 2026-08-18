"""Stage 0: function term -> candidate bacterial genes from literature.

This module provides the *deterministic* half of stage 0: retrieval, candidate
extraction, and the four-check hallucination audit. The judgement half (term
expansion, careful extraction, evidence tiering, and the microbiologist critic review)
is carried out by the agent following ``skills/function-gene-seed-agent/SKILL.md``.

Run this to get a machine-extracted first pass, then have the agent review and correct
it. Do not ship the machine pass unreviewed: it has no way to tell a demonstrated
function from a hedged one.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List

import pandas as pd

from ..clients.literature import LiteratureClient
from ..lexicon import TermLexicon

SEED_COLUMNS = ["gene", "system", "function_role", "organism_characterised",
                "organism_scope", "evidence_tier", "tier_basis", "pmids",
                "audit_status", "critic_verdict", "critic_note", "term_match"]

# Bacterial gene symbols: 3-4 lowercase letters + optional uppercase letter/digit.
GENE_RE = re.compile(r"\b([a-z]{3}[A-Z][0-9]?|[a-z]{3,4}[0-9]?)\b")

# Words that fit the gene-symbol shape but never are one.
NOT_GENES = {
    "the", "and", "for", "with", "that", "this", "was", "were", "has", "have", "not",
    "але", "from", "these", "those", "into", "than", "then", "when", "which", "were",
    "also", "both", "such", "been", "more", "most", "some", "each", "other", "same",
    "gene", "genes", "cell", "cells", "host", "iron", "data", "used", "using", "show",
    "shown", "here", "well", "between", "during", "after", "over", "under", "high",
    "low", "were", "via", "our", "one", "two", "three", "may", "can", "are", "its",
}


def retrieve(lexicon: TermLexicon, max_papers: int = 80,
             out_dir: str | Path = ".") -> pd.DataFrame:
    lc = LiteratureClient()
    q = lc.build_query({"term": lexicon.term, "synonyms": lexicon.synonyms,
                        "mechanism_terms": lexicon.mechanism_terms,
                        "negative_terms": lexicon.negative_terms})
    papers = pd.DataFrame(lc.search_epmc(q, max_results=max_papers))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    papers.to_csv(out / "papers.csv", index=False)
    (out / "query.txt").write_text(q)
    return papers


def extract_candidates(papers: pd.DataFrame, lexicon: TermLexicon,
                       window: int = 220) -> pd.DataFrame:
    """Machine first pass. Pulls gene-shaped tokens that sit near lexicon vocabulary.

    Deliberately high-recall and low-precision. The audit and the critic cut it down.
    """
    pos = lexicon.positive_re
    rows: List[Dict] = []
    for _, p in papers.iterrows():
        text = f"{p.get('title','')} {p.get('abstract','')}"
        if not text.strip():
            continue
        for m in pos.finditer(text):
            lo, hi = max(0, m.start() - window), min(len(text), m.end() + window)
            for g in GENE_RE.finditer(text[lo:hi]):
                sym = g.group(1)
                if sym.lower() in NOT_GENES or len(sym) < 3:
                    continue
                if sym.islower() and len(sym) == 3:
                    continue  # too generic without a capital suffix
                rows.append({
                    "gene": sym,
                    "system": "",
                    "function_role": "",
                    "organism_characterised": "",
                    "organism_scope": "",
                    "evidence_tier": "T5",
                    "tier_basis": "machine extraction from abstract; NOT REVIEWED",
                    "pmids": str(p.get("pmid", "")),
                    "audit_status": "",
                    "critic_verdict": "",
                    "critic_note": "",
                    "term_match": m.group(0),
                })
    if not rows:
        return pd.DataFrame(columns=SEED_COLUMNS)
    df = pd.DataFrame(rows)
    agg = (df.groupby("gene")
             .agg(pmids=("pmids", lambda s: ";".join(sorted({x for x in s if x}))),
                  term_match=("term_match", lambda s: "; ".join(sorted(set(s))[:5])),
                  n_mentions=("gene", "size"))
             .reset_index()
             .sort_values("n_mentions", ascending=False))
    for c in SEED_COLUMNS:
        if c not in agg.columns:
            agg[c] = ""
    agg["evidence_tier"] = "T5"
    agg["tier_basis"] = "machine extraction from abstract; NOT REVIEWED"
    return agg[SEED_COLUMNS + ["n_mentions"]]
