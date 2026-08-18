"""Stage 1 must disambiguate collisions using the lexicon, never by position."""
import pandas as pd
import pytest
from bfgm.lexicon import TermLexicon
from bfgm.stages import s1_ko
from test_kegg_matching import FakeKegg


@pytest.fixture
def lex():
    return TermLexicon(term="iron sequestration",
                       synonyms=["iron acquisition"],
                       mechanism_terms=["siderophore", "dihydroxybenzoate",
                                        "hydroxylysine", "enterobactin"],
                       negative_terms=["enterotoxin", "desaturase"])


def test_collision_resolved_to_on_term_ko(tmp_path, lex):
    seed = tmp_path / "seed.csv"
    pd.DataFrame([{"gene": "entA", "system": "enterobactin", "pmids": "1"}]).to_csv(seed, index=False)
    out = s1_ko.run(seed, lex, tmp_path, FakeKegg())
    row = out.iloc[0]
    assert row.kegg_ko == "K00216", "must pick the DHB dehydrogenase, not the enterotoxin"
    assert row.status_code == "CURATED"


def test_rejected_collisions_are_logged(tmp_path, lex):
    seed = tmp_path / "seed.csv"
    pd.DataFrame([{"gene": "entA", "system": "", "pmids": ""}]).to_csv(seed, index=False)
    s1_ko.run(seed, lex, tmp_path, FakeKegg())
    rej = pd.read_csv(tmp_path / "rejected_ko_collisions.csv")
    assert "K11059" in set(rej.rejected_ko), "the rejected KO must be auditable"


def test_unique_hit_marked_unique(tmp_path, lex):
    seed = tmp_path / "seed.csv"
    pd.DataFrame([{"gene": "iucA", "system": "aerobactin", "pmids": ""}]).to_csv(seed, index=False)
    out = s1_ko.run(seed, lex, tmp_path, FakeKegg())
    assert out.iloc[0].status_code == "UNIQUE"


def test_symbol_with_no_ko_is_flagged_not_dropped(tmp_path, lex):
    seed = tmp_path / "seed.csv"
    pd.DataFrame([{"gene": "xusB", "system": "", "pmids": ""}]).to_csv(seed, index=False)
    out = s1_ko.run(seed, lex, tmp_path, FakeKegg())
    assert len(out) == 1 and out.iloc[0].status_code == "NOT_IN_KO"
