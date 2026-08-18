"""Stage 2: metadata-only UniProt harvest (no sequences pulled).

Two axes, because neither alone is sound:
  A. gene symbol via ``gene_exact``  - high recall, poor precision
  B. Pfam domain                     - high precision, catches non-standard symbols

Axis B is **discovered, not hardcoded**. Pfam families are ranked by enrichment in the
axis-A confirmed set, so the domain axis adapts to whatever function the lexicon
describes. This is what makes the stage term-agnostic.

Metadata first, sequences later, deliberately: the manifest is inspected and curated
before any bulk sequence download.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Dict, List

import pandas as pd

from ..clients.uniprot import UniProtClient
from ..lexicon import TermLexicon


def discover_pfam(df: pd.DataFrame, lexicon: TermLexicon,
                  min_count: int = 3, top_n: int = 15) -> List[Dict]:
    """Rank Pfam families by frequency among on-term axis-A hits."""
    on = df[df["Protein names"].apply(lambda x: lexicon.classify_text(x) == "ON_TERM")]
    c = Counter()
    for v in on.get("Pfam", pd.Series(dtype=str)).dropna():
        for pf in str(v).split(";"):
            pf = pf.strip()
            if pf.startswith("PF"):
                c[pf] += 1
    return [{"pfam": pf, "n_on_term_hits": n}
            for pf, n in c.most_common(top_n) if n >= min_count]


def run(gene_ko_csv: str | Path, lexicon: TermLexicon, out_dir: str | Path,
        client: UniProtClient | None = None, use_pfam_axis: bool = True,
        progress=None) -> pd.DataFrame:
    up = client or UniProtClient()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    genes = sorted({str(g).strip() for g in pd.read_csv(gene_ko_csv)["gene"] if str(g).strip()})

    # ---- axis A: gene symbols ----
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
    axis_a = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    # ---- axis B: discovered Pfam families ----
    axis_b = pd.DataFrame()
    pfams: List[Dict] = []
    if use_pfam_axis and not axis_a.empty:
        pfams = discover_pfam(axis_a, lexicon)
        bframes = []
        for j, p in enumerate(pfams, 1):
            d = up.by_pfam(p["pfam"])
            if not d.empty:
                d = d.copy()
                d.insert(0, "query_term", p["pfam"])
                d.insert(0, "query_axis", "pfam_domain")
                bframes.append(d)
            manifest.append({"axis": "pfam_domain", "query_term": p["pfam"],
                             "n_entries": 0 if d.empty else len(d)})
            if progress:
                progress("pfam_domain", j, len(pfams))
        axis_b = pd.concat(bframes, ignore_index=True) if bframes else pd.DataFrame()

    both = [d for d in (axis_a, axis_b) if not d.empty]
    allrows = pd.concat(both, ignore_index=True) if both else pd.DataFrame()

    pd.DataFrame(manifest).to_csv(out / "uniprot_manifest.csv", index=False)
    pd.DataFrame(pfams).to_csv(out / "discovered_pfam.csv", index=False)
    if not allrows.empty:
        allrows.to_csv(out / "uniprot_stage2_metadata.tsv", sep="\t", index=False)
    return allrows


def curate(meta: pd.DataFrame, lexicon: TermLexicon, out_dir: str | Path,
           keep_symbols: List[str] | None = None) -> pd.DataFrame:
    """Classify every metadata row; reject symbol collisions.

    A row survives if it was found by the Pfam axis, OR its protein name is on-term,
    OR its symbol is on an explicit keep list (for genes whose protein name legitimately
    lacks the term vocabulary, e.g. a transporter named only for its family).
    """
    out = Path(out_dir)
    keep_symbols = {s.lower() for s in (keep_symbols or [])}
    pfam_accs = set(meta[meta.query_axis == "pfam_domain"]["Entry"]) if not meta.empty else set()

    def verdict(r):
        cls = lexicon.classify_text(r.get("Protein names", ""))
        if cls == "OFF_TERM":
            return "REJECTED_off_term"
        if r["Entry"] in pfam_accs:
            return "CONFIRMED_pfam"
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
