"""Stage 2: metadata-only UniProt harvest (no sequences pulled).

Gene symbol only, via ``gene_exact``. High recall, imperfect precision -
``curate()`` below is what resolves the resulting symbol collisions.

Metadata first, sequences later, deliberately: the manifest is inspected and curated
before any bulk sequence download.
"""
from __future__ import annotations

from pathlib import Path
from typing import List

import pandas as pd

from ..clients.uniprot import UniProtClient
from ..lexicon import TermLexicon


def run(gene_ko_csv: str | Path, lexicon: TermLexicon, out_dir: str | Path,
        client: UniProtClient | None = None, progress=None) -> pd.DataFrame:
    up = client or UniProtClient()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    genes = sorted({str(g).strip() for g in pd.read_csv(gene_ko_csv)["gene"] if str(g).strip()})

    frames, manifest = [], []
    for i, g in enumerate(genes, 1):
        d = up.by_gene_symbol(g)
        if not d.empty:
            d = d.copy()
            d.insert(0, "query_term", g)
            d.insert(0, "query_axis", "gene_symbol")
            frames.append(d)
        manifest.append({"axis": "gene_symbol", "query_term": g,
                         "n_entries": 0 if d.empty else len(d)})
        if progress and i % 25 == 0:
            progress("gene_symbol", i, len(genes))
    allrows = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    pd.DataFrame(manifest).to_csv(out / "uniprot_manifest.csv", index=False)
    if not allrows.empty:
        allrows.to_csv(out / "uniprot_stage2_metadata.tsv", sep="\t", index=False)
    return allrows


def curate(meta: pd.DataFrame, lexicon: TermLexicon, out_dir: str | Path,
           keep_symbols: List[str] | None = None) -> pd.DataFrame:
    """Classify every metadata row; reject symbol collisions.

    A row survives if its protein name is on-term, or its symbol is on an explicit
    keep list (for genes whose protein name legitimately lacks the term vocabulary,
    e.g. a transporter named only for its family).
    """
    out = Path(out_dir)
    keep_symbols = {s.lower() for s in (keep_symbols or [])}

    def verdict(r):
        cls = lexicon.classify_text(r.get("Protein names", ""))
        if cls == "OFF_TERM":
            return "REJECTED_off_term"
        if cls == "ON_TERM":
            return "CONFIRMED_name"
        if str(r.get("query_term", "")).lower() in keep_symbols:
            return "CURATED_KEEP"
        return "REJECTED_collision"

    meta = meta.copy()
    meta["curation_verdict"] = meta.apply(verdict, axis=1)
    kept = meta[~meta.curation_verdict.str.startswith("REJECTED")]
    rej = meta[meta.curation_verdict.str.startswith("REJECTED")]

    meta.to_csv(out / "uniprot_curated.tsv", sep="\t", index=False)
    rej[["query_term", "Entry", "Protein names", "Organism", "curation_verdict"]] \
        .to_csv(out / "rejected_symbol_collisions.tsv", sep="\t", index=False)
    return kept
