"""UniProtKB REST client.

Empirically established constraints baked in here:

* UniProt carries **no KEGG Orthology cross-reference**. Verified against P06971
  (E. coli FhuA): the entry has KEGG gene IDs, eggNOG, InterPro and Pfam, but no KO
  line, and ``xref:ko-`` queries return zero. KO linkage must therefore be built via
  ``xref_kegg`` -> KEGG ``/link/ko/<org>``. See ``bfgm.stages.s4_kegg_anchor``.
* ``gene_exact`` is used rather than ``gene`` to avoid substring matches.
* Pagination is cursor-based via the ``Link`` response header.
* Reviewed-only is the default. TrEMBL returns hundreds of thousands of entries for
  large transporter families and carries propagated annotation that cannot be trusted
  at gene-symbol level.
"""
from __future__ import annotations

import io
import re
import time
import urllib.parse
from typing import Dict, List, Optional, Sequence

import pandas as pd
import requests

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UA = {"User-Agent": "bfgm/1.0 (bacterial-function-gene-mapper)"}

META_FIELDS = (
    "accession,id,reviewed,protein_name,gene_names,organism_name,organism_id,"
    "lineage,length,protein_existence,xref_pfam,cc_subcellular_location"
)
FULL_FIELDS = (
    "accession,id,reviewed,protein_name,gene_names,organism_name,organism_id,"
    "lineage,lineage_ids,length,mass,sequence,protein_existence,xref_pfam,"
    "xref_interpro,ec,cc_function"
)


class UniProtClient:
    def __init__(self, taxonomy_id: int = 2, reviewed_only: bool = True,
                 throttle: float = 0.2):
        self.tax = taxonomy_id          # 2 = Bacteria
        self.reviewed_only = reviewed_only
        self.throttle = throttle

    def _filters(self) -> str:
        f = f" AND taxonomy_id:{self.tax}"
        if self.reviewed_only:
            f += " AND reviewed:true"
        return f

    def _request(self, url: str, retries: int = 4):
        last = None
        for attempt in range(retries):
            try:
                r = requests.get(url, headers=UA, timeout=120)
                if r.status_code == 200:
                    return r.text, r.headers.get("Link", "")
                last = f"HTTP {r.status_code}"
            except Exception as e:
                last = str(e)
            time.sleep(3 * (attempt + 1))
        raise RuntimeError(f"UniProt request failed ({last}): {url[:160]}")

    def count(self, query: str) -> int:
        url = f"{SEARCH}?query={urllib.parse.quote(query + self._filters())}&size=1"
        for attempt in range(3):
            try:
                r = requests.get(url, headers=UA, timeout=60)
                return int(r.headers.get("x-total-results", 0))
            except Exception:
                time.sleep(2 * (attempt + 1))
        return -1

    def search(self, query: str, fields: str = META_FIELDS,
               page: int = 500, cap: int = 20000) -> pd.DataFrame:
        """Paginated TSV search. Filters (taxonomy, reviewed) are applied automatically."""
        url = (f"{SEARCH}?query={urllib.parse.quote(query + self._filters())}"
               f"&format=tsv&fields={fields}&size={page}")
        frames, n = [], 0
        while url and n < cap:
            body, link = self._request(url)
            if not body.strip():
                break
            df = pd.read_csv(io.StringIO(body), sep="\t", low_memory=False)
            if df.empty:
                break
            frames.append(df)
            n += len(df)
            m = re.search(r'<([^>]+)>;\s*rel="next"', link)
            url = m.group(1) if m else None
            time.sleep(self.throttle)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def by_gene_symbol(self, gene: str, fields: str = META_FIELDS) -> pd.DataFrame:
        return self.search(f"gene_exact:{gene}", fields=fields)

    def by_accessions(self, accessions: Sequence[str], fields: str = FULL_FIELDS,
                      batch: int = 100, progress=None) -> pd.DataFrame:
        """Fetch full records for a known accession set. Batches of 100 stay under URL limits."""
        accs = list(accessions)
        frames = []
        for i in range(0, len(accs), batch):
            q = " OR ".join(f"accession:{a}" for a in accs[i:i + batch])
            url = (f"{SEARCH}?query={urllib.parse.quote(q)}"
                   f"&format=tsv&fields={fields}&size=500")
            body, _ = self._request(url)
            if body.strip():
                frames.append(pd.read_csv(io.StringIO(body), sep="\t", low_memory=False))
            if progress:
                progress(min(i + batch, len(accs)), len(accs))
            time.sleep(self.throttle)
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True).drop_duplicates("Entry")

    def fasta_for(self, accessions: Sequence[str], batch: int = 100) -> str:
        accs, out = list(accessions), []
        for i in range(0, len(accs), batch):
            q = " OR ".join(f"accession:{a}" for a in accs[i:i + batch])
            url = f"{SEARCH}?query={urllib.parse.quote(q)}&format=fasta&size=500"
            body, _ = self._request(url)
            if body.strip():
                out.append(body.rstrip("\n"))
            time.sleep(self.throttle)
        return "\n".join(out) + "\n"

    def kegg_xrefs(self, accessions: Sequence[str], batch: int = 100) -> pd.DataFrame:
        """accession -> KEGG gene IDs. First hop of the KO anchor."""
        accs, frames = list(accessions), []
        for i in range(0, len(accs), batch):
            q = " OR ".join(f"accession:{a}" for a in accs[i:i + batch])
            url = (f"{SEARCH}?query={urllib.parse.quote(q)}"
                   f"&format=tsv&fields=accession,xref_kegg,organism_id&size=500")
            body, _ = self._request(url)
            if body.strip():
                frames.append(pd.read_csv(io.StringIO(body), sep="\t"))
            time.sleep(self.throttle)
        if not frames:
            return pd.DataFrame(columns=["accession", "kegg_xref", "taxid"])
        df = pd.concat(frames, ignore_index=True).drop_duplicates(
            df_col := "Entry") if False else pd.concat(frames, ignore_index=True)
        df = df.drop_duplicates(df.columns[0])
        df.columns = ["accession", "kegg_xref", "taxid"]
        return df
