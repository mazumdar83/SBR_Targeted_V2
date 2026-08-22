"""Stage 5: pull sequences *for the KOs themselves*, not for the genes that found them.

Why this stage exists
---------------------
Stages 2 and 3 are gene-first: they retrieve UniProt accessions by gene symbol, and
stage 4 attaches KOs to whatever came back. A KO that no seed gene happened
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

Domain filter
-------------
A KO is a function, not a lineage: KEGG assigns the same K number across bacteria,
archaea, and eukaryotes wherever the orthology holds. Stage 3 (UniProt) is scoped to
Bacteria with ``--taxonomy 2``; this stage had no equivalent and would silently hand
back whatever organism KEGG listed first, which is fungal, plant, or archaeal in some
cases (verified: K10531 in the iron run picked *Neurospora tetrasperma* over any
bacterium, and some coverage-gap KOs list hundreds of eukaryotic members before the
first bacterium -- a per-organism live check tried first made stage 5 impractically
slow for exactly that reason).

``/list/organism`` would answer this in one call but currently 400s, and the FTP bulk
taxonomy dump (``ftp.genome.jp/pub/kegg/genes/taxonomy``) is unreachable from here
(likely the same licensed distribution the KO columns need Pathway Solutions for
commercially, per NOTICE). What does work in one call is KEGG's own organism
classification tree, ``/get/br:br08601`` (~500KB: domain -> Bacteria/Archaea ->
... -> organism code). ``_load_domain_map`` fetches and parses it once per run;
``is_bacterial`` is then a plain dict lookup, not a network call. Organisms it
excludes are logged to ``rejected_non_bacterial_kegg.csv``, never silently dropped.
"""
from __future__ import annotations

import re
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
import requests

BASE = "https://rest.kegg.jp"
UA = {"User-Agent": "bfgm/1.0 (bacterial-function-gene-mapper)"}
_LEAF = re.compile(r"^(\S+)\s{2,}(.*)$")


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


def _load_domain_map(cache: Dict[str, bool]) -> None:
    """One bulk fetch of KEGG's organism tree; populates `cache` code -> is_bacterial.

    Tree shape (see module docstring for why this endpoint): a top-level ``A`` line is
    either ``Eukaryotes`` or ``Prokaryotes``; under ``Prokaryotes`` a ``B`` line is
    either ``Bacteria`` or ``Archaea``. Every deeper line that isn't a further category
    header is a leaf: ``<code>  <name>``, at whatever depth that lineage happens to
    nest to. Bacteria is the only branch worth naming; everything else defaults False.
    """
    txt = _get("get/br:br08601", throttle=0.0)
    cur_a, cur_b = None, None
    for line in txt.split("\n"):
        if not line or not line[0].isalpha() or not line[0].isupper():
            continue
        rest = line[1:].strip()
        if line[0] == "A":
            cur_a, cur_b = rest, None
            continue
        if line[0] == "B":
            cur_b = rest
            continue
        m = _LEAF.match(rest)
        if not m:
            continue
        code = m.group(1)
        cache[code] = bool(cur_a and cur_a.startswith("Prokaryotes")
                            and cur_b and cur_b.startswith("Bacteria"))


def is_bacterial(org: str, cache: Dict[str, bool]) -> bool:
    """True if a KEGG organism code is classified under Bacteria.

    `cache` doubles as the domain map: empty means not yet loaded, so the first call
    in a run triggers one `_load_domain_map` fetch and every call after is a dict
    lookup. An organism code absent from the map (obsolete, or a NCBI-only genome not
    yet organized into KEGG's tree) is treated as non-bacterial, not as bacterial by
    default -- absence of evidence is not evidence of Bacteria.
    """
    if not cache:
        _load_domain_map(cache)
    return cache.get(org, False)


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

    domain_cache: Dict[str, bool] = {}
    rows, fasta, rejected_rows = [], [], []
    members_by_ko: Dict[str, int] = {}
    for i, ko in enumerate(kos, 1):
        try:
            members = genes_for_ko(ko)
        except Exception:
            members = []
        members_by_ko[ko] = len(members)
        bacterial = [g for g in members if is_bacterial(g.split(":")[0], domain_cache)]
        rejected_orgs = {g.split(":")[0] for g in members} - {g.split(":")[0] for g in bacterial}
        for org in sorted(rejected_orgs):
            rejected_rows.append({"KO": ko, "organism_code": org,
                                   "reason": "non-bacterial lineage"})
        reps = pick_representatives(bacterial, per_ko, prefer_orgs)
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
    pd.DataFrame(rejected_rows, columns=["KO", "organism_code", "reason"]
                ).to_csv(out / "rejected_non_bacterial_kegg.csv", index=False)

    covered = set(df.KO) if not df.empty else set()
    still = sorted(set(kos) - covered)
    reasons = [
        "no KEGG gene members; likely a KO with no sequenced representative"
        if members_by_ko.get(ko, 0) == 0 else
        "all KEGG members are non-bacterial (see rejected_non_bacterial_kegg.csv)"
        for ko in still
    ]
    pd.DataFrame({"KO": still, "reason": reasons}
                 ).to_csv(out / "ko_still_uncovered.csv", index=False)
    return df
