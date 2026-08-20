"""bfgm command line interface.

    bfgm init-term    "iron sequestration"     scaffold a run and draft a lexicon
    bfgm seed         --term-dir runs/iron/    stage 0 machine pass (agent reviews after)
    bfgm ko-map       --seed seed_genes.csv    stage 1
    bfgm uniprot-meta --run runs/iron/         stage 2
    bfgm sequences    --run runs/iron/         stage 3
    bfgm kegg-anchor  --run runs/iron/         stage 4
    bfgm ko-sequences --run runs/iron/         stage 5: sequences per KO (closes gaps)
    bfgm report       --run runs/iron/         build the workbook
    bfgm all          --run runs/iron/         stages 1 to 4 plus report
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from .clients.kegg import KeggClient
from .clients.uniprot import UniProtClient
from .lexicon import TermLexicon
from .stages import (s0_seed, s1_ko, s2_uniprot_meta, s3_sequences,
                     s4_kegg_anchor, s5_ko_sequences)


def _bar(label):
    def p(*a):
        print(f"  [{label}] " + " ".join(str(x) for x in a), file=sys.stderr)
    return p


def _lex(run: Path) -> TermLexicon:
    exp = run / "term_lexicon.expanded.json"
    return TermLexicon.load(exp if exp.exists() else run / "term_lexicon.json")


def cmd_init_term(a):
    run = Path(a.out or f"runs/{a.term.replace(' ', '_')}")
    run.mkdir(parents=True, exist_ok=True)
    lex = TermLexicon(term=a.term, synonyms=[], mechanism_terms=[], negative_terms=[])
    lex.save(run / "term_lexicon.json")
    print(f"Scaffolded {run}")
    print(f"Next: edit {run}/term_lexicon.json, or better, run the")
    print("function-gene-seed-agent skill to expand the term and build the seed properly.")


def cmd_seed(a):
    run = Path(a.run)
    lex = _lex(run)
    papers = s0_seed.retrieve(lex, max_papers=a.max_papers, out_dir=run)
    print(f"retrieved {len(papers)} papers")
    cand = s0_seed.extract_candidates(papers, lex)
    cand.to_csv(run / "seed_candidates_MACHINE.csv", index=False)
    print(f"extracted {len(cand)} candidate genes -> seed_candidates_MACHINE.csv")
    print("\nThis is an UNREVIEWED machine pass. Hand it to the")
    print("function-gene-seed-agent skill for audit, tiering and critic review")
    print("before running ko-map. Do not ship it as-is.")


def cmd_ko_map(a):
    run = Path(a.run)
    seed = Path(a.seed) if a.seed else run / "seed_genes.csv"
    if not seed.exists():
        sys.exit(f"missing {seed}. Run the seed agent first, or pass --seed.")
    df = s1_ko.run(seed, _lex(run), run, KeggClient(cache_dir=a.cache))
    print(df.status_code.value_counts().to_string())


def cmd_uniprot_meta(a):
    run = Path(a.run)
    lex = _lex(run)
    up = UniProtClient(taxonomy_id=a.taxonomy, reviewed_only=not a.include_unreviewed)
    meta = s2_uniprot_meta.run(run / "gene_ko_map.csv", lex, run, up,
                               progress=lambda ax, i, n: print(f"  {ax} {i}/{n}", file=sys.stderr))
    if meta.empty:
        sys.exit("no UniProt hits; widen the lexicon or check the seed")
    keep = json.loads(Path(a.keep_symbols).read_text()) if a.keep_symbols else []
    kept = s2_uniprot_meta.curate(meta, lex, run, keep_symbols=keep)
    print(f"metadata rows {len(meta)}, kept {len(kept)}, "
          f"unique accessions {kept.Entry.nunique()}")
    print("Inspect uniprot_curated.tsv and rejected_symbol_collisions.tsv "
          "before running `sequences`.")


def cmd_sequences(a):
    run = Path(a.run)
    cur = pd.read_csv(run / "uniprot_curated.tsv", sep="\t", low_memory=False)
    kept = cur[~cur.curation_verdict.str.startswith("REJECTED")]
    up = UniProtClient(taxonomy_id=a.taxonomy, reviewed_only=not a.include_unreviewed)
    df = s3_sequences.run(kept, run, up, progress=_bar("seq"))
    print(f"proteins {len(df)}")


def cmd_kegg_anchor(a):
    run = Path(a.run)
    prot = pd.read_csv(run / "proteins.tsv", sep="\t", low_memory=False)
    anchor = s4_kegg_anchor.run(prot, run / "gene_ko_map.csv", _lex(run), run,
                                KeggClient(cache_dir=a.cache), progress=_bar("kegg"))
    print(anchor.anchor_class.value_counts().to_string())
    new = pd.read_csv(run / "discovered_kos.csv")
    if len(new):
        print(f"\n{len(new)} on-term KOs found that the seed missed. "
              "Add them to the lexicon and consider re-running:")
        print(new.head(10).to_string(index=False))


def cmd_ko_sequences(a):
    run = Path(a.run)
    df = s5_ko_sequences.run(run, per_ko=a.per_ko, gaps_only=not a.all_kos,
                             progress=_bar("ko-seq"))
    if df.empty:
        print("no KOs to fetch (no coverage gaps, or run --all-kos)")
        return
    print(f"{df.KO.nunique()} KOs, {len(df)} sequences -> ko_sequences.fasta")
    import pandas as _pd
    unc = _pd.read_csv(run / "ko_still_uncovered.csv")
    if len(unc):
        print(f"{len(unc)} KOs still have no sequenced representative in KEGG")


def cmd_report(a):
    from .report import build
    p = build(Path(a.run))
    print(f"wrote {p}")


def cmd_all(a):
    cmd_ko_map(a); cmd_uniprot_meta(a); cmd_sequences(a); cmd_kegg_anchor(a)
    cmd_ko_sequences(a); cmd_report(a)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="bfgm", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, run=True):
        if run:
            p.add_argument("--run", required=True, help="run directory")
        p.add_argument("--cache", default="data/cache/kegg")
        p.add_argument("--taxonomy", type=int, default=2, help="NCBI taxid filter (2=Bacteria)")
        p.add_argument("--include-unreviewed", action="store_true",
                       help="include TrEMBL. Expect very large result sets.")

    p = sub.add_parser("init-term"); p.add_argument("term"); p.add_argument("--out")
    p.set_defaults(func=cmd_init_term)

    p = sub.add_parser("seed"); common(p); p.add_argument("--max-papers", type=int, default=80)
    p.set_defaults(func=cmd_seed)

    p = sub.add_parser("ko-map"); common(p); p.add_argument("--seed")
    p.set_defaults(func=cmd_ko_map)

    p = sub.add_parser("uniprot-meta"); common(p)
    p.add_argument("--keep-symbols", help="JSON list of symbols to keep despite off-term names")
    p.set_defaults(func=cmd_uniprot_meta)

    p = sub.add_parser("sequences"); common(p); p.set_defaults(func=cmd_sequences)
    p = sub.add_parser("kegg-anchor"); common(p); p.set_defaults(func=cmd_kegg_anchor)
    p = sub.add_parser("ko-sequences"); common(p)
    p.add_argument("--per-ko", type=int, default=5,
                   help="representative sequences per KO, one per organism")
    p.add_argument("--all-kos", action="store_true",
                   help="fetch for every KO, not just those with no sequence from stage 3")
    p.set_defaults(func=cmd_ko_sequences)

    p = sub.add_parser("report"); common(p); p.set_defaults(func=cmd_report)

    p = sub.add_parser("all"); common(p)
    p.add_argument("--seed"); p.add_argument("--keep-symbols")
    p.add_argument("--per-ko", type=int, default=5)
    p.add_argument("--all-kos", action="store_true")
    p.set_defaults(func=cmd_all)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == "__main__":
    main()
