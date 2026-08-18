# Worked example: "iron sequestration"

## Phase 1 output (term_lexicon.json)

```json
{
  "term": "iron sequestration",
  "synonyms": ["iron acquisition", "iron piracy", "iron scavenging", "ferric uptake",
               "nutritional immunity"],
  "mechanism_terms": ["siderophore", "TonB-dependent", "hemophore", "ferritin",
                      "transferrin", "lactoferrin", "ferrous transport", "NEAT domain"],
  "negative_terms": ["iron-sulfur cluster assembly", "heme biosynthesis",
                     "phytosiderophore"]
}
```

Note `negative_terms`. Without excluding Fe-S assembly, `iscR`, `sufB` and `iscS`
dominate the extraction and the set drifts away from acquisition.

## Phase 4 audit, real catches

| candidate | check | outcome |
|---|---|---|
| `iucA` | symbol resolves in KEGG (K03894) | PASS |
| `desA` | resolves, but to acyl-lipid desaturases in most taxa | `COLLISION_RISK` |
| `hasA` | resolves to both hemophore and hyaluronan synthase | `COLLISION_RISK` |
| `shr` | resolves to a plant transcription factor as well | `COLLISION_RISK` |

None of these are rejections. Stage 1 disambiguates them against KO definitions.

## Phase 6 critic, real corrections

- `hupB` in M. tuberculosis: FLAG. Moonlights as a carboxymycobactin receptor, but its
  primary annotation is a nucleoid-associated protein. Scope narrowed to genus, note
  added that KO annotation will not capture the iron role.
- `zupT`: FLAG. ZIP-family transporter with broad divalent specificity. Real, but not
  iron-specific. `function_role` corrected.
- `bfrA` as beta-fructosidase: REJECT. Extraction hit a symbol collision with the
  R. inulinivorans fructan gene. Quarantined as `CRITIC_REJECT`.

## Final counts

165 candidates extracted, 12 quarantined, 153 accepted or flagged.
Tier distribution T1:1, T2:41, T3:38, T4:52, T5:21.

## Downstream result

Stage 1 mapped 153 symbols to 122 distinct KOs, of which 34 needed manual collision
resolution. Stage 4 later surfaced two iron KOs the seed missed (`K25283`, `K28698`),
which were fed back into the lexicon. **That feedback loop is the intended use.**
