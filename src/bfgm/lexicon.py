"""Term lexicon: the mechanism that makes this pipeline term-agnostic.

Everything downstream that used to be hardcoded iron regex is now derived from the
lexicon produced by stage 0. Swap the lexicon and the same code maps a different
function.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List


@dataclass
class TermLexicon:
    term: str
    synonyms: List[str] = field(default_factory=list)
    mechanism_terms: List[str] = field(default_factory=list)
    negative_terms: List[str] = field(default_factory=list)

    @classmethod
    def load(cls, path: str | Path) -> "TermLexicon":
        d = json.loads(Path(path).read_text())
        return cls(term=d["term"], synonyms=d.get("synonyms", []),
                   mechanism_terms=d.get("mechanism_terms", []),
                   negative_terms=d.get("negative_terms", []))

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps({
            "term": self.term, "synonyms": self.synonyms,
            "mechanism_terms": self.mechanism_terms,
            "negative_terms": self.negative_terms}, indent=1))

    def _pattern(self, terms: List[str]) -> re.Pattern:
        if not terms:
            return re.compile(r"(?!x)x")  # matches nothing
        parts = sorted({re.escape(t) for t in terms if t}, key=len, reverse=True)
        return re.compile("|".join(parts), re.IGNORECASE)

    @property
    def positive_re(self) -> re.Pattern:
        return self._pattern([self.term] + self.synonyms + self.mechanism_terms)

    @property
    def negative_re(self) -> re.Pattern:
        return self._pattern(self.negative_terms)

    def classify_text(self, text: str) -> str:
        """ON_TERM, OFF_TERM (matches a negative term), or NO_MATCH."""
        t = str(text or "")
        if self.negative_terms and self.negative_re.search(t):
            return "OFF_TERM"
        return "ON_TERM" if self.positive_re.search(t) else "NO_MATCH"

    def expand_from_ko_definitions(self, definitions: List[str], top_n: int = 40) -> List[str]:
        """Grow mechanism_terms from the KO definitions of confirmed hits.

        This is the bootstrap that lets the classifier recognise vocabulary the user
        never supplied. Called by stage 1 after KO mapping.
        """
        from collections import Counter
        stop = {"protein", "system", "family", "subunit", "transport", "putative",
                "component", "binding", "type", "domain", "containing", "like",
                "uncharacterized", "probable", "chain", "small", "large"}
        c = Counter()
        for d in definitions:
            for w in re.findall(r"[A-Za-z][A-Za-z-]{3,}", str(d).lower()):
                if w not in stop:
                    c[w] += 1
        new = [w for w, n in c.most_common(top_n * 3) if n >= 2]
        added = [w for w in new if w not in {m.lower() for m in self.mechanism_terms}]
        self.mechanism_terms.extend(added[:top_n])
        return added[:top_n]
