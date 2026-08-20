# bfgm — bacterial function-to-gene mapper

Give it a function term. Get back the bacterial genes that perform it, their KEGG
Orthology assignments, their UniProt sequences, and a taxonomic map — with every
step audited and every rejection logged.

```bash
bfgm init-term "iron sequestration"
# run the function-gene-seed-agent skill to build seed_genes.csv
bfgm all --run runs/iron_sequestration/
```

The term is the only thing that changes. `"bile acid 7-alpha-dehydroxylation"`,
`"mucin glycan degradation"`, `"tryptophan catabolism to indole"` and
`"beta-glucuronidase activity"` all run through the same code.

## How it works

```
  "iron sequestration"
          |
   [stage 0] seed agent .............. literature -> genes
          |                            audited, tiered, critic-reviewed
          v  seed_genes.csv
   [stage 1] ko-map .................. gene symbols -> KEGG KOs
          |                            exact token match, collision curation
          v  gene_ko_map.csv
   [stage 2] uniprot-meta ............ metadata only, two axes
          |                            gene symbol + discovered Pfam
          v  uniprot_curated.tsv
   [stage 3] sequences ............... sequences + taxonomic lineage
          |
          v  proteins.fasta, taxonomy_mapping.tsv
   [stage 4] kegg-anchor ............. accession -> KEGG gene -> KO
          |                            validates the chain, finds missed KOs
          v  proteins_anchored.tsv, discovered_kos.csv --> feeds back to stage 0
   [stage 5] ko-sequences ............ KO -> KEGG genes -> aaseq
          |                            closes gaps stage 3 could not reach
          v  ko_sequences.fasta
   [stage 6] build-db ................ merge, dedupe, BLAST/DIAMOND database
          |
          v  <term>_reference.fasta + _reference_map.tsv + HOW_TO_SEARCH.md
```

## What makes it term-agnostic

Earlier versions of this pipeline had the target function baked in as regular
expressions. Two mechanisms replaced that:

**The term lexicon.** Stage 0 produces `term_lexicon.json` holding synonyms,
mechanism vocabulary, and — importantly — negative terms. Every downstream
classification decision reads from it. Stage 1 then grows the lexicon automatically
from the KO definitions of confirmed hits, so the classifier learns vocabulary the
user never supplied.

**Pfam discovery.** Stage 2's domain axis is not a hardcoded family list. It ranks
Pfam families by enrichment among on-term gene-symbol hits and queries the top ones.
Point the pipeline at a different function and it finds that function's domains.

Negative terms matter more than they look. Without excluding iron-sulfur cluster
assembly from an iron run, `iscR`, `sufB` and `iscS` swamp the extraction and the set
drifts off target.

## Two things learned the hard way

**UniProt has no KEGG Orthology cross-reference.** Verified against P06971 (*E. coli*
FhuA): the entry carries KEGG gene IDs, eggNOG, InterPro and Pfam, but no KO line, and
`xref:ko-` queries return zero. Stage 4 therefore builds the KO link in two hops,
through `xref_kegg` and then KEGG's per-organism `/link/ko/<org>` endpoint. The global
`/link/genes/ko` endpoint is unreliable at scale and is not used anywhere in this
package.

**Gene symbols collide constantly and the collisions are not obvious.** In a single
iron run: `desA` matched acyl-lipid desaturases, `acsA` returned 145 acetyl-CoA
synthetases, `hasA` returned six hyaluronan synthases against one real hemophore,
`tbpA` returned thiamine-binding proteins, `shr` matched a plant transcription factor,
and `entA` through `entH` matched staphylococcal enterotoxins. Substring matching makes
this dramatically worse, so matching is exact-token throughout, and every collision is
resolved against the lexicon and written to a rejection log.

## Run it in Claude Code

Open the repo in Claude Code and type:

```
/map-function "iron sequestration"
```

That scaffolds the run, delegates stage 0 to the **function-seed-researcher** subagent,
which hands off to the **microbiologist-critic** subagent for independent review, pauses
for your confirmation, then runs stages 1 to 4 via the **pipeline-runner** subagent.

`/seed-only "<term>"` stops after the gene list. `/verify-ids` checks whether PMIDs,
gene symbols, or organism names actually resolve.

The critic runs as a separate subagent deliberately: an agent that produced a gene set
is a poor judge of whether it over-claimed, and the isolation keeps the noisy review
work out of the main transcript.

## Install

```bash
git clone <your-repo-url> && cd bfgm
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
make test
```

No API keys required. Set `NCBI_API_KEY` to raise E-utilities rate limits if running
many terms.

## Commands

| command | does |
|---|---|
| `bfgm init-term "<term>"` | scaffold a run directory and lexicon stub |
| `bfgm seed --run <dir>` | stage 0 machine pass (**unreviewed** — hand to the skill) |
| `bfgm ko-map --run <dir>` | stage 1: symbols to KOs |
| `bfgm uniprot-meta --run <dir>` | stage 2: metadata harvest, two axes |
| `bfgm sequences --run <dir>` | stage 3: sequences and taxonomy |
| `bfgm kegg-anchor --run <dir>` | stage 4: KO anchoring and validation |
| `bfgm ko-sequences --run <dir>` | stage 5: sequences per KO, closes coverage gaps |
| `bfgm build-db --run <dir>` | stage 6: merged, deduplicated BLAST/DIAMOND reference |
| `bfgm report --run <dir>` | build the workbook |
| `bfgm all --run <dir>` | stages 1 to 4 plus report |

Useful flags: `--taxonomy 2157` for Archaea, `--include-unreviewed` to add TrEMBL
(expect very large result sets), `--keep-symbols keep.json` to retain genes whose
protein names legitimately lack the term vocabulary.

## The seed agent skill

`skills/function-gene-seed-agent/` is the stage 0 skill. Install it to
`~/.claude/skills/` or point Claude Code at the repo.

It does the part that needs judgement: term expansion, careful extraction that does
not upgrade hedged language, a four-check hallucination audit, evidence tiering, and
a skeptical microbiologist critic review with a hard three-iteration cap.

`bfgm seed` gives you a machine first pass, but it is high-recall and low-precision by
design and has no way to tell a demonstrated function from a hedged one. **Do not ship
it unreviewed.**

The audit's four checks: every PMID must resolve, every gene symbol must resolve in
KEGG or UniProt or NCBI Gene, colliding symbols get flagged (not dropped — stage 1
handles them), and organism names must be real taxa. Failures go to `quarantine.csv`
with a reason. Nothing is silently dropped; the quarantine file is a deliverable.

## The searchable output

Stage 6 is the point of the whole thing if you are screening genomes. It merges the
UniProt (gene-first) and KEGG (KO-first) sequences, collapses exact duplicates, and
writes one FASTA with headers you can parse straight out of `-outfmt 6`:

    >bfgm|<n>|<KO>|<gene>|<source_acc> <description> [<organism>]

So a BLAST or DIAMOND hit already carries the KO and the gene symbol in `sseqid`; no
join needed to get a KO profile out of a hit table. `HOW_TO_SEARCH.md` is written into
the run directory with the exact commands and sensible thresholds.

The two sequence directions matter. Stage 3 is gene-first and misses any KO no seed gene
happened to retrieve; stage 5 goes KO-first and fills those in. Running only stage 3
left 49 of 176 iron KOs with no sequence at all.

## The feedback loop

Stage 4 classifies every anchored KO against the lexicon. `ANCHOR_NEW_ON_TERM` means
the KO is on-term but the seed missed it — these land in `discovered_kos.csv` and
should go back into the lexicon. In the iron run this surfaced `K25283` (iron
siderophore permease) and `K28698` (putrebactin synthase).

`ANCHOR_ADJACENT_OFF_TERM` is the other useful signal: proteins that share domains
with the target function but do something else. In the iron run this caught 103
vitamin B12 transport components that entered through the shared FecCD and
periplasmic-binding folds. They are labelled rather than deleted, because they are the
correct negative controls for anything you train downstream.

## Licensing

MIT for the code. The data is another matter — KO columns make an output KEGG-derived
and commercial use needs a Pathway Solutions licence, while UniProt sequences and
taxonomy are CC BY 4.0 and unrestricted. See `NOTICE`.

## Push to GitHub

```bash
./scripts/init_repo.sh https://github.com/<user>/bfgm.git
# or, with the gh CLI installed:
./scripts/init_repo.sh
```
