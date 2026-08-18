# seed_genes.csv schema

Stage 1 (`bfgm ko-map`) reads this file. Column names are load-bearing; do not rename.

| column | required | description |
|---|---|---|
| `gene` | yes | Gene symbol exactly as written in the source. No case normalisation. |
| `system` | yes | Operon, cluster, or system the gene belongs to (e.g. `enterobactin`, `Isd haem relay`). Free text. |
| `function_role` | yes | One line on what the gene product does within the system. |
| `organism_characterised` | yes | Organism the function was demonstrated in. Binomial where possible. |
| `organism_scope` | yes | `strain`, `species`, `genus`, or `broad`. How far the claim generalises. |
| `evidence_tier` | yes | `T1` to `T5`. See evidence-tiers.md. |
| `tier_basis` | yes | The specific study design justifying the tier. Not a citation. |
| `pmids` | yes | Semicolon-separated, verified to resolve. |
| `audit_status` | yes | `PASS`, `COLLISION_RISK`, `ORGANISM_FLAG`. Set by audit_seed.py. |
| `critic_verdict` | yes | `ACCEPT` or `FLAG`. `REJECT` rows go to quarantine.csv instead. |
| `critic_note` | no | One line, required when verdict is `FLAG`. |
| `term_match` | no | Which lexicon terms this gene matched. Populated by the extractor. |

## quarantine.csv

Same columns plus `quarantine_reason`, one of:

- `PMID_UNRESOLVED` - cited PMID does not exist
- `SYMBOL_UNRESOLVED` - symbol found in no reference database
- `CRITIC_REJECT` - critic determined the gene does not perform the stated function
- `ORGANISM_UNRESOLVED` - named organism is not a real taxon

Never delete quarantined rows. They are the audit trail.

## term_lexicon.json

```json
{
  "term": "iron sequestration",
  "synonyms": ["iron acquisition", "iron piracy", "ferric uptake"],
  "mechanism_terms": ["siderophore", "TonB-dependent", "ferritin"],
  "negative_terms": ["iron-sulfur cluster assembly", "heme biosynthesis"],
  "confirmed_at": "2026-08-18"
}
```

Downstream stages use `synonyms` and `mechanism_terms` to build the protein-name
classifier, and `negative_terms` to reject domain-adjacent false positives.
