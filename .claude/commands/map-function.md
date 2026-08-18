---
description: Run the full bfgm pipeline for a biological function term
argument-hint: "<function term>, e.g. iron sequestration"
allowed-tools: Bash, Read, Write, Edit, Glob, Grep, Task, WebSearch, WebFetch
---

Run the complete bfgm pipeline for: **$ARGUMENTS**

Environment check first:

- repo state: !`ls pyproject.toml src/bfgm/cli.py 2>&1 | head -2`
- bfgm importable: !`python -c "import bfgm; print(bfgm.__version__)" 2>&1 | head -1`

If bfgm is not importable, run `pip install -e ".[dev]"` before continuing.

## Steps

1. Scaffold the run:
   ```bash
   python -m bfgm.cli init-term "$ARGUMENTS"
   ```

2. Delegate stage 0 to the **function-seed-researcher** subagent. Pass it the term and
   the run directory. It expands the lexicon, retrieves literature, extracts genes, runs
   the hallucination audit, and hands off to the **microbiologist-critic** subagent
   before returning.

   Do not do stage 0 yourself. The separation between extractor and critic is what makes
   the output trustworthy, and it collapses if one context does both.

3. Show the user the seed summary and **stop for confirmation** before spending API
   budget on stages 1 to 4. Report the gene count, tier distribution, quarantine count,
   and any `COLLISION_RISK` symbols.

4. On confirmation, delegate stages 1 to 4 to the **pipeline-runner** subagent.

5. Report: the workbook path, and anything in `discovered_kos.csv` — those are on-term
   KOs the seed missed and they should go back into the lexicon.

## Notes

- If the term is too broad to build a specific lexicon from ("metabolism", "virulence"),
  say so and ask the user to narrow it before spending a retrieval budget.
- The seed file is the contract between stage 0 and everything downstream. If the user
  already has one, skip to step 4.
