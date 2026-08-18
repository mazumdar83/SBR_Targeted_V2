---
description: Run only bfgm stage 0 (literature to audited gene list), no downstream stages
argument-hint: "<function term>"
allowed-tools: Bash, Read, Write, Edit, Glob, Grep, Task, WebSearch, WebFetch
---

Run **only stage 0** of bfgm for: **$ARGUMENTS**

Scaffold with `python -m bfgm.cli init-term "$ARGUMENTS"`, then delegate to the
**function-seed-researcher** subagent, which will hand off to **microbiologist-critic**
before returning.

Stop after the seed. Do not run ko-map or anything downstream.

Report the gene count, tier distribution, quarantine count with reasons, and the
`COLLISION_RISK` symbols that stage 1 will need to resolve. Then tell the user the
command to continue:

```
python -m bfgm.cli all --run runs/<term>/
```
