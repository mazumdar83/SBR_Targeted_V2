# CLAUDE.md

Operating instructions for Claude Code in this repo.

## What this is

`bfgm` maps a biological function term to the bacterial genes that perform it, their
KEGG Orthology assignments, their UniProt sequences, and a taxonomic map. The term is
the only input that changes between runs.

## Quick start

```
/map-function "iron sequestration"
```

That orchestrates the whole thing. `/seed-only "<term>"` stops after stage 0.

## Architecture

```
stage 0  function-seed-researcher (subagent)  literature -> genes, audited
         microbiologist-critic (subagent)     independent review
stage 1  bfgm ko-map                          symbols -> KEGG KOs
stage 2  bfgm uniprot-meta                    metadata, two axes
stage 3  bfgm sequences                       sequences + taxonomy
stage 4  bfgm kegg-anchor                     validates chain, finds missed KOs
stage 5  bfgm ko-sequences                    KO-first sequences, closes gaps
stage 6  bfgm build-db                        merged BLAST/DIAMOND reference
```

Stage 4's `discovered_kos.csv` feeds back into the stage 0 lexicon. That loop is the
intended use, not an afterthought.

## Setup

```bash
pip install -e ".[dev]"
make test          # 16 offline tests
make smoke         # 4 live API tests
```

No API keys needed. `NCBI_API_KEY` raises E-utilities rate limits if set.

## Rules for working in this repo

**Never skip the audit or the critic.** Stage 0 has two gates: the four-check
hallucination audit (`audit_seed.py`) and the microbiologist critic subagent. Both are
mandatory. A gene list that skipped them is not a bfgm output.

**Never let one context both extract and review.** The critic runs as a separate
subagent for a reason. An agent that produced a gene set is a poor judge of whether it
over-claimed.

**Never hardcode the function.** Everything term-specific lives in `term_lexicon.json`.
If you find yourself writing a regex for a particular biology, that belongs in the
lexicon instead. This was the whole point of the refactor.

**Never invent identifiers.** No PMID, gene symbol, or organism name goes into an
output unless it resolves against a live database. Use
`/verify-ids` or `verify_identifiers.py` when unsure.

**Quarantine is a deliverable.** Rejected rows are written with reasons, never silently
dropped. `quarantine.csv` and the `rejected_*` files are the audit trail.

**Never upgrade an evidence tier.** T1 to T5 reflects the strongest single study design.
Ten association studies do not make a T2.

## Sequence retrieval runs in two directions

Stage 3 is gene-first (UniProt, by symbol). Stage 5 is KO-first (KEGG aaseq).
Neither alone is complete: stage 3 misses KOs no seed gene retrieved, stage 5 only
covers what KEGG has. Stage 6 merges and dedupes them into the searchable reference.

Absence from the reference is not absence from biology. Systems with no KO and no
reviewed UniProt record produce no sequences, so a genome screen cannot detect them.
Check `ko_still_uncovered.csv` and the NOT_IN_KO rows before concluding anything is
missing from a genome.

## Two constraints the code depends on

**UniProt has no KEGG Orthology cross-reference.** Verified against P06971 (E. coli
FhuA): KEGG gene IDs, eggNOG, InterPro and Pfam are present, but no KO line, and
`xref:ko-` returns zero. Stage 4 therefore anchors in two hops, via `xref_kegg` and
KEGG's per-organism `/link/ko/<org>`. `tests/test_smoke_live.py` asserts this; if that
test starts failing, upstream changed and stage 4 can be simplified.

**The global KEGG `/link/genes/ko` endpoint is unreliable at scale.** Only the
per-organism route is used. Do not "optimise" this into a single call.

## Symbol collisions

The dominant failure mode. Matching is exact-token throughout, never substring. Real
examples from one iron run: `acsA` returned 145 acetyl-CoA synthetases, `hasA` returned
six hyaluronan synthases against one hemophore, `tbpA` returned thiamine-binding
proteins, `shr` matched a plant transcription factor, `entA` through `entH` matched
staphylococcal enterotoxins.

Collisions are resolved against the lexicon, never by taking the first hit, and every
rejection is logged with the definition that triggered it.

## Layout

```
.claude/agents/      three subagents
.claude/commands/    three slash commands
.claude/skills/      the stage 0 skill
src/bfgm/clients/    KEGG, UniProt, literature
src/bfgm/stages/     s0 to s4
src/bfgm/lexicon.py  the term-agnostic mechanism
tests/               offline + live-marked
runs/                run outputs, gitignored
```

## Licensing

MIT for code. KO columns make an output KEGG-derived and need a Pathway Solutions
licence commercially; UniProt sequences and taxonomy are CC BY 4.0. See `NOTICE`.
Dropping KO columns leaves a freely usable asset.

## Style

Match what is here: type hints, docstrings that explain why rather than what, module
docstrings that record empirically established constraints. Tests encode failure modes
we actually hit, not happy paths. No em-dashes in generated documents.
