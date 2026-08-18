#!/usr/bin/env python3
"""Europe PMC retrieval for a term lexicon. Writes papers.csv and query.txt."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
from bfgm.lexicon import TermLexicon
from bfgm.stages import s0_seed

ap = argparse.ArgumentParser()
ap.add_argument("lexicon"); ap.add_argument("--out", default="."); ap.add_argument("--max", type=int, default=80)
a = ap.parse_args()
lex = TermLexicon.load(a.lexicon)
papers = s0_seed.retrieve(lex, max_papers=a.max, out_dir=a.out)
print(f"{len(papers)} papers -> {a.out}/papers.csv")
cand = s0_seed.extract_candidates(papers, lex)
cand.to_csv(Path(a.out) / "seed_candidates_MACHINE.csv", index=False)
print(f"{len(cand)} candidate genes (UNREVIEWED) -> seed_candidates_MACHINE.csv")
