"""KEGG REST client.

Two behaviours here are load-bearing and were established empirically:

1. The global ``/link/genes/ko`` endpoint is unreliable at scale. The per-organism
   ``/link/ko/<org>`` loop is used instead everywhere in this package.
2. KO symbol matching must be **exact token** matching against the symbol field.
   Substring matching produces false positives at a high rate (``desA`` matching
   acyl-lipid desaturases, ``entA`` matching enterotoxin A, and so on).
"""
from __future__ import annotations

import os
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import requests

BASE = "https://rest.kegg.jp"
UA = {"User-Agent": "bfgm/1.0 (bacterial-function-gene-mapper)"}


class KeggClient:
    def __init__(self, cache_dir: str | Path = "data/cache/kegg", throttle: float = 0.12):
        self.cache = Path(cache_dir)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.throttle = throttle
        self._ko_sym: Optional[Dict[str, List[str]]] = None
        self._ko_def: Optional[Dict[str, str]] = None
        self._sym_index: Optional[Dict[str, List[str]]] = None

    # ---------- low level ----------

    def _get(self, path: str, cache_name: Optional[str] = None, retries: int = 3) -> str:
        if cache_name:
            p = self.cache / cache_name
            if p.exists() and p.stat().st_size > 0:
                return p.read_text()
        url = f"{BASE}/{path}"
        last = None
        for attempt in range(retries):
            try:
                r = requests.get(url, headers=UA, timeout=90)
                if r.status_code == 200:
                    if cache_name:
                        (self.cache / cache_name).write_text(r.text)
                    time.sleep(self.throttle)
                    return r.text
                last = f"HTTP {r.status_code}"
            except Exception as e:  # network flake
                last = str(e)
            time.sleep(2 * (attempt + 1))
        raise RuntimeError(f"KEGG GET {path} failed: {last}")

    # ---------- KO catalogue ----------

    def load_ko_catalogue(self) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
        """Full KO list. Returns (K -> [symbols], K -> definition)."""
        if self._ko_sym is not None:
            return self._ko_sym, self._ko_def  # type: ignore[return-value]
        txt = self._get("list/ko", "ko_list.tsv")
        sym: Dict[str, List[str]] = {}
        dfn: Dict[str, str] = {}
        for line in txt.rstrip("\n").split("\n"):
            if not line.strip():
                continue
            k, rest = line.split("\t", 1)
            if ";" in rest:
                syms, d = rest.split(";", 1)
            else:
                syms, d = "", rest
            sym[k] = [s.strip() for s in syms.split(",") if s.strip()]
            dfn[k] = d.strip()
        self._ko_sym, self._ko_def = sym, dfn
        return sym, dfn

    def symbol_index(self) -> Dict[str, List[str]]:
        """Lowercased exact symbol -> list of K numbers."""
        if self._sym_index is not None:
            return self._sym_index
        sym, _ = self.load_ko_catalogue()
        idx: Dict[str, List[str]] = defaultdict(list)
        for k, syms in sym.items():
            for s in syms:
                idx[s.lower()].append(k)
        self._sym_index = dict(idx)
        return self._sym_index

    def ko_definition(self, ko: str) -> str:
        _, dfn = self.load_ko_catalogue()
        return dfn.get(ko, "")

    def ec_for(self, ko: str) -> str:
        m = re.search(r"\[EC:([^\]]+)\]", self.ko_definition(ko))
        return m.group(1) if m else ""

    def match_symbol(self, gene: str) -> List[str]:
        """Exact token match only. Returns every KO whose symbol field contains `gene`."""
        return list(self.symbol_index().get(gene.strip().lower(), []))

    # ---------- pathway / module ----------

    def ko_to_pathway(self) -> Dict[str, List[str]]:
        txt = self._get("link/pathway/ko", "ko_pathway.tsv")
        out: Dict[str, List[str]] = defaultdict(list)
        for line in txt.rstrip("\n").split("\n"):
            if not line.strip():
                continue
            k, p = line.split("\t")
            p = p.replace("path:", "")
            if p.startswith("map"):  # drop ko-prefixed duplicates
                out[k.replace("ko:", "")].append(p)
        return dict(out)

    def ko_to_module(self) -> Dict[str, List[str]]:
        txt = self._get("link/module/ko", "ko_module.tsv")
        out: Dict[str, List[str]] = defaultdict(list)
        for line in txt.rstrip("\n").split("\n"):
            if not line.strip():
                continue
            k, m = line.split("\t")
            out[k.replace("ko:", "")].append(m.replace("md:", ""))
        return dict(out)

    def _name_map(self, path: str, cache: str, prefix: str) -> Dict[str, str]:
        txt = self._get(path, cache)
        out = {}
        for line in txt.rstrip("\n").split("\n"):
            if "\t" in line:
                a, b = line.split("\t", 1)
                out[a.replace(prefix, "")] = b
        return out

    def pathway_names(self) -> Dict[str, str]:
        return self._name_map("list/pathway", "pathway_list.tsv", "path:")

    def module_names(self) -> Dict[str, str]:
        return self._name_map("list/module", "module_list.tsv", "md:")

    # ---------- per-organism gene -> KO ----------

    def organism_kos(self, org: str) -> Dict[str, str]:
        """Gene-to-KO for one KEGG organism code.

        Uses ``/link/ko/<org>``. The global ``/link/genes/ko`` endpoint is NOT used;
        it is unreliable at scale.
        """
        txt = self._get(f"link/ko/{org}", f"ko_{org}.tsv")
        out = {}
        for line in txt.rstrip("\n").split("\n"):
            if "\t" in line:
                g, k = line.split("\t")
                out[g] = k.replace("ko:", "")
        return out

    def bulk_organism_kos(self, orgs: Iterable[str], progress=None) -> Dict[str, str]:
        """Union gene-to-KO map across many organism codes. Failures are skipped, not fatal."""
        out: Dict[str, str] = {}
        orgs = list(orgs)
        failed = []
        for i, o in enumerate(orgs, 1):
            try:
                out.update(self.organism_kos(o))
            except Exception:
                failed.append(o)
            if progress and i % 50 == 0:
                progress(i, len(orgs), len(out), len(failed))
        if failed:
            (self.cache / "failed_orgs.txt").write_text("\n".join(failed))
        return out
