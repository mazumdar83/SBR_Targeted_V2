# Contributing

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
make test
```

## Before opening a PR

```bash
make test    # offline suite, must pass
make lint
make smoke   # live API tests, run if you touched a client
```

## Rules specific to this codebase

**Do not hardcode a biology.** Anything term-specific belongs in `term_lexicon.json`,
not in a regex in the source. The pipeline's value is that swapping the lexicon retargets
it; a hardcoded pattern breaks that.

**Symbol matching stays exact-token.** Never substring. If you are tempted, read
`tests/test_kegg_matching.py` first, which encodes why.

**Rejections are logged, never dropped.** Any filtering step writes what it removed and
why, to a file. If you add a filter, add its audit output.

**Record empirical constraints in the module docstring.** When you discover that an API
behaves unexpectedly, write it down where the next person will hit it. The existing
docstrings on `clients/uniprot.py` and `clients/kegg.py` are the pattern.

**New tests encode real failure modes.** The suite is not there for coverage; it exists
so that bugs we already paid for do not come back. A test named after the bug beats a
test named after the function.

## Adding a stage

Stages live in `src/bfgm/stages/` and take (input path or DataFrame, lexicon, out_dir).
They write their outputs to the run directory and return a DataFrame. Wire the CLI
subcommand in `cli.py` and add the output to `SHEETS` in `report.py`.
