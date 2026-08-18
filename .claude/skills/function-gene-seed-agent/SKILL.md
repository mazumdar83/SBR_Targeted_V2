---
name: function-gene-seed-agent
description: Spins up a research agent that takes a free-text biological function term (for example "iron sequestration", "bile acid dehydroxylation", "mucin degradation", "tryptophan catabolism") and returns a verified, evidence-tiered list of bacterial genes and proteins that carry out that function in the published biomedical literature. Every gene is anchored to a real PMID and a real database identifier, every claim passes a hallucination audit, and the whole set is reviewed by a skeptical microbiologist critic before release. Use this skill whenever the user wants to go from a function, phenotype, pathway, or capability to a bacterial gene list, or asks "which bacterial genes do X", "find genes for X", "what carries out X in bacteria", "seed a gene set for X", or wants to start a genomics pipeline from a function rather than from an organism. Always use this skill rather than ad hoc literature search when the gene list will feed downstream KEGG, UniProt, or sequence work, because downstream stages depend on the schema and the audit flags this produces.
---

# Function to Bacterial Gene Seed Agent

Stage 0 of the bacterial-function-gene-mapper pipeline. Turns a function term into
`seed_genes.csv`, which every downstream stage consumes.

**The contract**: nothing leaves this skill unless it is anchored to a verifiable
external identifier. A gene symbol with no database record and no real PMID is a
hallucination, and this skill is built to catch it rather than pass it downstream.

## When NOT to use this

If the user already has a gene list, skip to stage 1 (`bfgm ko-map`). This skill is
only for the function-to-gene hop.

## Pipeline position

```
[THIS SKILL] -> seed_genes.csv -> ko-map -> uniprot-meta -> sequences -> kegg-anchor
```

## Workflow

### Phase 1: Term expansion

Do not search on the raw term alone. Expand it first, and write the expansion to
`term_lexicon.json` because later phases and downstream stages reuse it.

1. Synonyms and near-synonyms. "iron sequestration" also appears as iron acquisition,
   iron piracy, iron scavenging, ferric uptake, siderophore-mediated uptake,
   nutritional immunity.
2. Mechanism nouns the literature actually uses: the transporter families, the
   chelator classes, the enzyme activities.
3. Negative terms that will pull in the wrong biology and must be excluded later.
   For iron: heme biosynthesis, iron-sulfur cluster assembly (unless wanted),
   iron storage in plants.
4. Ask the user to confirm the expansion before searching. A bad lexicon wastes the
   whole run and the user usually knows their field better than the expansion does.

Use `ask_user_input_v0` for the confirmation if the term is ambiguous.

### Phase 2: Literature retrieval

Query Europe PMC and PubMed via `scripts/search_literature.py`. Prefer reviews for
coverage and primary papers for evidence tier. Target 30 to 80 papers; below 20 the
gene set will be thin, above 100 the extraction gets noisy.

Restrict to bacteria explicitly. Human and plant homologs of the same function will
otherwise dominate the results.

### Phase 3: Gene extraction

For each paper, extract candidate genes into the schema in
`references/output-schema.md`. Record for every candidate:

- gene symbol as written in the paper
- the system or operon it belongs to
- the organism it was characterised in
- what the paper actually showed, in one line
- the PMID
- the study design

**Extract only what the paper states.** If a paper says a gene is "predicted to" or
"may" do something, that is what gets recorded. Do not upgrade hedged language into
a functional claim. This is the single most common failure mode here.

### Phase 4: Hallucination audit (mandatory, never skipped)

Run `scripts/audit_seed.py`. It performs four checks, and every candidate must pass
checks 1 and 2 to survive:

1. **PMID exists.** Resolve every PMID against NCBI or Europe PMC. A PMID that does
   not resolve means the citation was fabricated. Quarantine the whole row.
2. **Gene symbol resolves.** The symbol must appear in at least one of: KEGG KO symbol
   fields, UniProt `gene_exact` for bacteria, or NCBI Gene. A symbol resolving nowhere
   is either invented or a typo. Quarantine.
3. **Symbol collision flag.** If the symbol resolves to entries whose descriptions have
   nothing to do with the term lexicon, flag `COLLISION_RISK`. This is not a rejection.
   Real gene symbols collide constantly across taxa, and stage 1 resolves them properly.
4. **Organism plausibility.** The organism named must be a real taxon in NCBI Taxonomy
   and must be bacterial. Non-bacterial organisms get flagged, not dropped, since a
   fungal or host homolog is sometimes the correct anchor for a claim.

Write everything that fails to `quarantine.csv` with the reason. Never silently drop a
row. The quarantine file is evidence that the audit ran and is reviewed by the critic.

### Phase 5: Evidence tiering

Assign per `references/evidence-tiers.md`. The tier reflects the strongest study design
supporting that gene, not the number of papers. Ten association studies do not make a
T2. Record the tier and the specific basis for it.

### Phase 6: Microbiologist critic review

Adopt the persona in `references/critic-persona.md` and work the ten-point checklist.
The critic is deliberately skeptical and its default posture is that the extraction
over-claims.

Return `ACCEPT`, `FLAG`, or `REJECT` per gene with a one-line reason:

- `ACCEPT` — evidence supports the gene-to-function link at the stated tier
- `FLAG` — link is plausible but the tier is wrong, the organism scope is overstated,
  or the mechanism is assumed rather than shown. Keeps the row, corrects the metadata.
- `REJECT` — the gene does not do what the extraction says. Moves to `quarantine.csv`.

**Hard cap of three redesign iterations.** If the set is still failing after three
passes, stop and report what is unresolved rather than iterating forever.

### Phase 7: Emit

Write `seed_genes.csv` per `references/output-schema.md`, plus `quarantine.csv`,
`term_lexicon.json`, and `seed_report.md`. Report counts: candidates extracted,
quarantined by reason, accepted, flagged, and the tier distribution.

Then tell the user the exact next command:

```
bfgm ko-map --seed seed_genes.csv --out runs/<term>/
```

## Rules that are not negotiable

- Never invent a PMID. If unsure a paper exists, search for it and drop it if absent.
- Never invent a gene symbol. Genes come from papers or from database records, never
  from inference about what a symbol "should" be.
- Never upgrade evidence tier to make a set look stronger.
- Report absent data as absent. An empty result for a function is a real finding and
  is more useful than a padded list.
- Symbol collisions are expected and are stage 1's job. Flag, do not discard.
- Record every rejection with its reason. The quarantine file is a deliverable.

## Reference files

- `references/output-schema.md` — exact column spec for seed_genes.csv, read before emitting
- `references/evidence-tiers.md` — T1 to T5 definitions with worked examples
- `references/critic-persona.md` — the microbiologist persona and ten-point checklist
- `references/worked-example.md` — a full iron sequestration run, input to output

## Scripts

- `scripts/search_literature.py` — Europe PMC and PubMed retrieval
- `scripts/audit_seed.py` — the four-check hallucination audit
- `scripts/verify_identifiers.py` — batch PMID, gene symbol, and taxon resolution
