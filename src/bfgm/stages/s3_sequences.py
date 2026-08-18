"""Stage 3: sequences and taxonomic mapping for the curated accession set."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ..clients.uniprot import UniProtClient

RANKS = ["phylum", "class", "order", "family", "genus"]


def _rank(lineage: str, key: str) -> str:
    for part in str(lineage).split(","):
        part = part.strip()
        if part.endswith(f"({key})"):
            return part.rsplit("(", 1)[0].strip()
    return ""


def run(curated: pd.DataFrame, out_dir: str | Path,
        client: UniProtClient | None = None, progress=None) -> pd.DataFrame:
    up = client or UniProtClient()
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    accs = sorted(curated["Entry"].unique())
    df = up.by_accessions(accs, progress=progress)
    if df.empty:
        raise RuntimeError("stage 3 returned no records; check the curated accession set")

    prov = (curated.groupby("Entry")
            .agg(matched_query_terms=("query_term", lambda s: "; ".join(sorted(set(s)))),
                 curation_verdict=("curation_verdict", lambda s: "; ".join(sorted(set(s)))))
            .reset_index())
    df = df.merge(prov, on="Entry", how="left")
    df.to_csv(out / "proteins.tsv", sep="\t", index=False)

    (out / "proteins.fasta").write_text(up.fasta_for(accs))

    tax = df[["Organism (ID)", "Organism", "Taxonomic lineage"]] \
        .drop_duplicates("Organism (ID)").copy()
    tax.columns = ["taxid", "organism", "lineage"]
    for r in RANKS:
        tax[r] = tax["lineage"].apply(lambda x, r=r: _rank(x, r))
    tax = tax.merge(df.groupby("Organism (ID)").size().rename("n_proteins"),
                    left_on="taxid", right_index=True, how="left")
    tax.sort_values("n_proteins", ascending=False).to_csv(
        out / "taxonomy_mapping.tsv", sep="\t", index=False)

    # QA gates: these must hold or the run is not deliverable
    assert df["Sequence"].notna().all(), "empty sequences returned"
    leaked = tax[~tax.lineage.fillna("").str.contains("Bacteria")]
    if len(leaked):
        leaked.to_csv(out / "WARNING_non_bacterial.tsv", sep="\t", index=False)
    return df
