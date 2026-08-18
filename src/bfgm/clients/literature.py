"""Europe PMC and NCBI E-utilities client. Used by stage 0 for retrieval and by the
audit for PMID verification.

No API key is required for the volumes this pipeline uses, but set NCBI_API_KEY in the
environment to raise the E-utilities rate limit if you are running many terms.
"""
from __future__ import annotations

import os
import time
from typing import Dict, Iterable, List, Set

import requests

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = {"User-Agent": "bfgm/1.0 (bacterial-function-gene-mapper)"}


class LiteratureClient:
    def __init__(self, throttle: float = 0.34):
        self.throttle = throttle
        self.api_key = os.environ.get("NCBI_API_KEY")

    def search_epmc(self, query: str, page_size: int = 100, max_results: int = 200) -> List[Dict]:
        out, cursor = [], "*"
        while len(out) < max_results:
            params = {"query": query, "format": "json", "pageSize": page_size,
                      "cursorMark": cursor, "resultType": "core"}
            r = requests.get(EPMC, params=params, headers=UA, timeout=60)
            r.raise_for_status()
            j = r.json()
            hits = j.get("resultList", {}).get("result", [])
            if not hits:
                break
            for h in hits:
                out.append({
                    "pmid": h.get("pmid", ""),
                    "pmcid": h.get("pmcid", ""),
                    "doi": h.get("doi", ""),
                    "title": h.get("title", ""),
                    "abstract": h.get("abstractText", ""),
                    "journal": h.get("journalTitle", ""),
                    "year": h.get("pubYear", ""),
                    "type": h.get("pubType", ""),
                    "is_open_access": h.get("isOpenAccess", "N"),
                })
            nxt = j.get("nextCursorMark")
            if not nxt or nxt == cursor:
                break
            cursor = nxt
            time.sleep(self.throttle)
        return out[:max_results]

    def build_query(self, lexicon: Dict) -> str:
        """Term lexicon -> Europe PMC boolean query, bacteria-restricted."""
        terms = [lexicon["term"]] + lexicon.get("synonyms", [])
        mech = lexicon.get("mechanism_terms", [])
        pos = " OR ".join(f'"{t}"' for t in terms + mech)
        q = f"({pos}) AND (bacteria OR bacterial OR microbial OR prokaryot*)"
        q += " AND (gene OR genes OR operon OR protein OR transporter OR enzyme)"
        for neg in lexicon.get("negative_terms", []):
            q += f' NOT "{neg}"'
        return q

    def verify_pmids(self, pmids: Iterable[str]) -> Set[str]:
        """Return the subset that actually resolves. Anything absent is a fabricated citation."""
        pmids = [str(p).strip() for p in pmids if str(p).strip()]
        good: Set[str] = set()
        for i in range(0, len(pmids), 200):
            batch = pmids[i:i + 200]
            params = {"db": "pubmed", "id": ",".join(batch), "retmode": "json"}
            if self.api_key:
                params["api_key"] = self.api_key
            try:
                r = requests.get(f"{EUTILS}/esummary.fcgi", params=params,
                                 headers=UA, timeout=60)
                j = r.json().get("result", {})
                for p in batch:
                    rec = j.get(p)
                    if rec and not rec.get("error"):
                        good.add(p)
            except Exception:
                pass
            time.sleep(self.throttle)
        return good

    def verify_taxa(self, names: Iterable[str]) -> Dict[str, str]:
        """Organism name -> NCBI taxid. Empty string means unresolved."""
        out = {}
        for n in {str(x).strip() for x in names if str(x).strip()}:
            params = {"db": "taxonomy", "term": n, "retmode": "json"}
            if self.api_key:
                params["api_key"] = self.api_key
            try:
                r = requests.get(f"{EUTILS}/esearch.fcgi", params=params,
                                 headers=UA, timeout=45)
                ids = r.json().get("esearchresult", {}).get("idlist", [])
                out[n] = ids[0] if ids else ""
            except Exception:
                out[n] = ""
            time.sleep(self.throttle)
        return out
