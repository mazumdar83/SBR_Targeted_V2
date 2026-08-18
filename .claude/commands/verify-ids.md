---
description: Check whether PMIDs, gene symbols, or organism names actually resolve
argument-hint: "--pmids / --genes / --organisms, then values"
allowed-tools: Bash
---

Verify identifiers against live databases. Arguments: **$ARGUMENTS**

Run the resolver with the appropriate flag:

```bash
python .claude/skills/function-gene-seed-agent/scripts/verify_identifiers.py \
  --pmids <ids> --genes <symbols> --organisms <names>
```

Report plainly which resolve and which do not. A PMID that does not resolve is a
fabricated citation. A gene symbol that resolves nowhere is invented or a typo. A
symbol that resolves to several unrelated KOs is a collision and needs the lexicon to
disambiguate it, not a guess.
