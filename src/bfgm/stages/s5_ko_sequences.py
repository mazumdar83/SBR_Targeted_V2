"""Stage 5: pull sequences *for the KOs themselves*, not for the genes that found them.

Why this stage exists
---------------------
Stages 2 and 3 are gene-first: they retrieve UniProt accessions by gene symbol and Pfam
domain, and stage 4 attaches KOs to whatever came back. A KO that no seed gene happened
to retrieve therefore ends up with zero sequence support, even though KEGG knows about
it. In the iron run that was 49 of 176 KOs, clustered in the NRPS modules and ABC
subunits.

This stage inverts the direction: KO -> KEGG genes -> amino acid sequences. It closes
the gap and does not depend on UniProt at all.

Endpoint note
-------------
``/link/genes/<KO>`` is a *targeted* query and is reliable. This is not the same as the
global ``/link/genes/ko`` dump, which is unreliable at scale and is used nowhere in this
package. ``/get/<gene>/aaseq`` accepts up to 10 gene IDs joined with ``+``.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
import requests

BASE = "https://rest.kegg.jp"
UA = {"User-Agent": "bfgm/1.0 (bacterial-function-gene-mapper)"}


def _get(path: str, retries: int = 3, throttle: float = 0.12) -> str:
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(f"{BASE}/{path}", headers=UA, timeout=90)
            if r.status_code == 200:
                time.sleep(throttle)
                return r.text
            if r.status_code == 404:
                return ""
            last = f"HTTP {r.status_code}"
        except Exception as e:
            last = str(e)
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"KEGG GET {path} failed: {last}")


def genes_for_ko(ko: str) -> List[str]:
    """KEGG gene IDs for one KO. Targeted endpoint, not the global dump."""
    txt = _get(f"link/genes/{ko}")
    out = []
    for line in txt.rstrip("\n").split("\n"):
        if "\t" in line:
            out.append(line.split("\t")[1].strip())
    return out


def pick_representatives(genes: List[str], per_ko: int,
                         prefer_orgs: Optional[Iterable[str]] = None) -> List[str]:
    """One gene per organism, so a KO with 4000 members does not return 4000 E. coli strains.

    Organisms in `prefer_orgs` are taken first; the rest fill the remaining slots in
    KEGG's own order, which is roughly taxonomically grouped.
    """
    prefer = list(prefer_orgs or [])
    seen_org, chosen = set(), []
    for want in (True, False):
        for g in genes:
            if len(chosen) >= per_ko:
                break
            org = g.split(":")[0]
            if org in seen_org:
                continue
            is_pref = org in prefer
            if is_pref is want:
                seen_org.add(org)
                chosen.append(g)
    return chosen


def fetch_aaseq(gene_ids: List[str], batch: int = 10) -> str:
    """Amino acid FASTA. KEGG /get accepts at most 10 entries per call."""
    out = []
    for i in range(0, len(gene_ids), batch):
        txt = _get("get/" + "+".join(gene_ids[i:i + batch]) + "/aaseq")
        if txt.strip():
            out.append(txt.rstrip("\n"))
    return "\n".join(out)


def run(run_dir: str | Path, kos: Optional[Iterable[str]] = None,
        per_ko: int = 5, gaps_only: bool = True,
        prefer_orgs: Optional[Iterable[str]] = None, progress=None) -> pd.DataFrame:
    """Retrieve sequences per KO.

    gaps_only=True  -> only KOs with no sequence from stage 3 (reads ko_coverage_gaps.csv)
    gaps_only=False -> every KO in gene_ko_map.csv
    """
    out = Path(run_dir)
    out.mkdir(parents=True, exist_ok=True)

    if kos is None:
        if gaps_only and (out / "ko_coverage_gaps.csv").exists():
            g = pd.read_csv(out / "ko_coverage_gaps.csv")
            kos = [k for k in g.get("KO", pd.Series(dtype=str)).dropna()]
        else:
            m = pd.read_csv(out / "gene_ko_map.csv")
            kos = sorted({t for v in m.kegg_ko.dropna()
                          for t in str(v).replace(";", " ").split()
                          if t.startswith("K") and len(t) == 6})
    kos = list(kos)
    if not kos:
        return pd.DataFrame(columns=["KO", "kegg_gene_id", "organism_code", "n_members"])

    rows, fasta = [], []
    for i, ko in enumerate(kos, 1):
        try:
            members = genes_for_ko(ko)
        except Exception:
            members = []
        reps = pick_representatives(members, per_ko, prefer_orgs)
        if reps:
            fasta.append(fetch_aaseq(reps))
        for g in reps:
            rows.append({"KO": ko, "kegg_gene_id": g,
                         "organism_code": g.split(":")[0],
                         "n_members_in_ko": len(members)})
        if progress and i % 10 == 0:
            progress(i, len(kos), len(rows))

    df = pd.DataFrame(rows)
    df.to_csv(out / "ko_sequence_manifest.csv", index=False)
    (out / "ko_sequences.fasta").write_text("\n".join(fasta) + "\n" if fasta else "")

    covered = set(df.KO) if not df.empty else set()
    still = sorted(set(kos) - covered)
    pd.DataFrame({"KO": still,
                  "reason": "no KEGG gene members; likely a KO with no sequenced representative"}
                 ).to_csv(out / "ko_still_uncovered.csv", index=False)
    return df
