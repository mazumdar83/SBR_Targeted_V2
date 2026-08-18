# Evidence tiers

The tier reflects the **strongest single study design** supporting the gene-to-function
link. It is not a count. Ten association studies do not add up to a T2.

| tier | definition | worked example |
|---|---|---|
| **T1** | Human in vivo evidence that the gene is required for the function in a human host | Neisseria gonorrhoeae `tbpA`: human volunteer urethral challenge showed the receptor is required for efficient colonisation |
| **T2** | Animal model with an isogenic mutant showing a causal defect in the function | S. aureus `sfa`/`sbn` double mutant loses growth on holo-transferrin and shows reduced organ burden in mice |
| **T3** | Ex vivo human material, or in vitro growth on the physiological substrate as sole source | C. jejuni growth on ferri-lactoferrin as sole iron source |
| **T4** | Biochemical, structural, or genetic characterisation without a functional in vivo test | HasA to HasR haem transfer solved structurally; no in vivo requirement shown |
| **T5** | Genomic or bioinformatic prediction only: homology, operon context, or domain architecture | A `feoB` homolog called from an assembly with no experimental follow-up |

## Rules

- **Never promote a tier.** An animal model that measures colonisation but not the
  function itself is not T2 for that function.
- **T5 is legitimate output.** Predicted genes are useful seeds. They must be labelled
  so downstream stages do not treat them as confirmed.
- **Record `tier_basis` as the design**, not the citation. "RCT, n=32 mice, isogenic
  deletion, grip strength endpoint" is a basis. "Smith 2024" is not.
- **Null results are recorded, not dropped.** A gene shown not to perform a function
  belongs in the seed with a critic note, because it stops someone re-deriving it later.
