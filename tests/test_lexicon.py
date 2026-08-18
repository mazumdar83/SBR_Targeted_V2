from bfgm.lexicon import TermLexicon


def lex():
    return TermLexicon(
        term="iron sequestration",
        synonyms=["iron acquisition", "ferric uptake"],
        mechanism_terms=["siderophore", "TonB-dependent", "ferritin"],
        negative_terms=["iron-sulfur cluster assembly", "heme biosynthesis"],
    )


def test_on_term():
    assert lex().classify_text("ferric enterobactin receptor, TonB-dependent") == "ON_TERM"


def test_negative_terms_win_over_positive():
    # contains "iron" but is the wrong biology; must not be accepted
    assert lex().classify_text("iron-sulfur cluster assembly protein IscR") == "OFF_TERM"


def test_no_match():
    assert lex().classify_text("acetyl-coenzyme A synthetase") == "NO_MATCH"


def test_empty_lexicon_matches_nothing():
    l = TermLexicon(term="x")
    assert l.classify_text("anything at all") == "NO_MATCH"


def test_expansion_adds_vocabulary():
    l = lex()
    before = len(l.mechanism_terms)
    added = l.expand_from_ko_definitions([
        "ferric enterobactin transport system permease",
        "ferric enterobactin transport system ATP-binding",
        "hemophore receptor HasR",
        "hemophore secretion protein",
    ])
    assert added and len(l.mechanism_terms) > before


def test_roundtrip(tmp_path):
    p = tmp_path / "lex.json"
    lex().save(p)
    assert TermLexicon.load(p).negative_terms == lex().negative_terms
