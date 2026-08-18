"""Stage 4: anchor UniProt accessions back to KEGG KOs.

UniProt has no KO cross-reference, so the link is built in two hops:

    accession -> KEGG gene ID (UniProt xref_kegg) -> KO (KEGG /link/ko/<org>)

The per-organism route is used; the global /link/genes/ko endpoint is unreliable.

The anchor is not just a lookup. It is the validation step: it catches domain-adjacent
false positives that the Pfam axis let through, and it surfaces on-term KOs that the
seed missed. Feed those back into the lexicon.
"""
from __future__ import annotations

from pathlib import Path
from typing import Set

import pandas as pd

from ..clients.kegg import KeggClient
from ..clients.uniprot import UniProtClient
from ..lexicon import TermLexicon


def run(proteins: pd.DataFrame, gene_ko_csv: str | Path, lexicon: TermLexicon,
        out_dir: str | Path, kegg: KeggClient | None = None,
        up: UniProtClient | None = None, progress=None) -> pd.DataFrame:
    kegg = kegg or KeggClient()
    up = up or UniProtClient()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # hop 1
    x = up.kegg_xrefs(sorted(proteins["Entry"].unique()))
    x = x[x.kegg_xref.notna()].copy()
    x["genes"] = x.kegg_xref.str.rstrip(";").str.split(";")
    ex = x.explode("genes")
    ex["genes"] = ex.genes.str.strip()
    ex = ex[ex.genes.str.contains(":", na=False)]
    ex["org"] = ex.genes.str.split(":").str[0]

    # hop 2
    g2k = kegg.bulk_organism_kos(sorted(ex.org.unique()), progress=progress)
    ex["KO"] = ex.genes.map(g2k)
    ex["ko_definition"] = ex.KO.map(lambda k: kegg.ko_definition(k) if pd.notna(k) else "")

    seed_kos: Set[str] = set()
    for v in pd.read_csv(gene_ko_csv)["kegg_ko"].dropna():
        for t in str(v).replace(";", " ").split():
            if t.startswith("K") and len(t) == 6:
                seed_kos.add(t)

    def classify(r):
        if pd.isna(r.KO):
            return "NO_KO_ASSIGNED"
        if r.KO in seed_kos:
            return "ANCHOR_SEED_KO"
        cls = lexicon.classify_text(r.ko_definition)
        if cls == "ON_TERM":
            return "ANCHOR_NEW_ON_TERM"      # KO the seed missed; feed back to lexicon
        if cls == "OFF_TERM":
            return "ANCHOR_ADJACENT_OFF_TERM"
        return "ANCHOR_OFFTARGET"

    ex["anchor_class"] = ex.apply(classify, axis=1)
    anchor = ex[["accession", "genes", "org", "KO", "ko_definition", "anchor_class"]] \
        .rename(columns={"genes": "kegg_gene_id", "org": "kegg_org_code"})
    anchor.to_csv(out / "kegg_anchor_table.tsv", sep="\t", index=False)

    # KOs the seed proposed that have no sequence support
    hit = set(anchor[anchor.anchor_class == "ANCHOR_SEED_KO"].KO)
    gaps = pd.DataFrame({"KO": sorted(seed_kos - hit)})
    gaps["ko_definition"] = gaps.KO.map(kegg.ko_definition) if not gaps.empty else pd.Series(dtype=str)
    gaps["status"] = "no reviewed bacterial sequence via this route" if not gaps.empty else pd.Series(dtype=str)
    gaps.to_csv(out / "ko_coverage_gaps.csv", index=False)

    # newly discovered on-term KOs, for lexicon feedback
    new = (anchor[anchor.anchor_class == "ANCHOR_NEW_ON_TERM"]
           [["KO", "ko_definition"]].drop_duplicates())
    new.to_csv(out / "discovered_kos.csv", index=False)

    roll = (anchor.groupby("accession")
            .agg(kegg_gene_ids=("kegg_gene_id", lambda s: "; ".join(sorted(set(s)))),
                 kegg_org_codes=("kegg_org_code", lambda s: "; ".join(sorted(set(s)))),
                 kegg_KOs=("KO", lambda s: "; ".join(sorted({v for v in s if pd.notna(v)}))),
                 anchor_class=("anchor_class", lambda s: "; ".join(sorted(set(s)))))
            .reset_index())
    merged = proteins.merge(roll, left_on="Entry", right_on="accession", how="left") \
                     .drop(columns=["accession"])
    merged["anchor_class"] = merged.anchor_class.fillna("NO_KEGG_XREF")
    merged.to_csv(out / "proteins_anchored.tsv", sep="\t", index=False)
    return anchor
