"""Stage 5 must fetch by KO, deduplicate by organism, and never lose an uncovered KO."""
import pandas as pd
from bfgm.stages import s5_ko_sequences as s5


def test_representative_picking_is_one_per_organism():
    genes = ["eco:b1", "eco:b2", "eco:b3", "sty:c1", "mtu:Rv1", "mtu:Rv2"]
    reps = s5.pick_representatives(genes, per_ko=5)
    orgs = [g.split(":")[0] for g in reps]
    assert len(orgs) == len(set(orgs)), "must not return two genes from one organism"
    assert set(orgs) == {"eco", "sty", "mtu"}


def test_per_ko_cap_respected():
    genes = [f"org{i}:g" for i in range(50)]
    assert len(s5.pick_representatives(genes, per_ko=5)) == 5


def test_preferred_organisms_come_first():
    genes = ["aaa:g1", "eco:g2", "bbb:g3"]
    reps = s5.pick_representatives(genes, per_ko=2, prefer_orgs=["eco"])
    assert reps[0].startswith("eco"), "preferred organism should be taken first"


def test_empty_ko_list_returns_empty_frame(tmp_path):
    df = s5.run(tmp_path, kos=[], gaps_only=False)
    assert df.empty


def test_uncovered_kos_are_recorded_not_dropped(tmp_path, monkeypatch):
    monkeypatch.setattr(s5, "genes_for_ko", lambda ko: [] if ko == "K99999" else ["eco:b1"])
    monkeypatch.setattr(s5, "fetch_aaseq", lambda ids, batch=10: ">eco:b1\nMKV\n")
    monkeypatch.setattr(s5, "is_bacterial", lambda org, cache: True)
    s5.run(tmp_path, kos=["K00216", "K99999"], gaps_only=False)
    unc = pd.read_csv(tmp_path / "ko_still_uncovered.csv")
    assert "K99999" in set(unc.KO), "a KO with no members must be logged, not silently dropped"


_FAKE_BRITE_TREE = """+E	KEGG Organism
!
AEukaryotes (2)
B  Fungi (1)
C    Sordariomycetes (1)
D      nte  Neurospora tetrasperma FGSC 2508
B  Plants (1)
C    Poaceae (1)
D      osa  Oryza sativa
AProkaryotes (3)
B  Bacteria (2)
C    Enterobacteria (2)
D      eco  Escherichia coli K-12 MG1655
D      sty  Salmonella enterica
B  Archaea (1)
C    Halobacteria (1)
D      hma  Haloarcula marismortui
"""


def test_domain_map_classifies_bacteria_only(monkeypatch):
    monkeypatch.setattr(s5, "_get", lambda path, throttle=0.12: _FAKE_BRITE_TREE)
    cache = {}
    assert s5.is_bacterial("eco", cache) is True
    assert s5.is_bacterial("sty", cache) is True
    assert s5.is_bacterial("nte", cache) is False, "fungi must not pass as bacterial"
    assert s5.is_bacterial("osa", cache) is False, "plants must not pass as bacterial"
    assert s5.is_bacterial("hma", cache) is False, "archaea must not pass as bacterial"


def test_domain_map_loads_once_and_is_cached(monkeypatch):
    calls = []
    def fake_get(path, throttle=0.12):
        calls.append(path)
        return _FAKE_BRITE_TREE
    monkeypatch.setattr(s5, "_get", fake_get)
    cache = {}
    s5.is_bacterial("eco", cache)
    s5.is_bacterial("sty", cache)
    s5.is_bacterial("nte", cache)
    assert calls == ["get/br:br08601"], "the bulk tree must be fetched exactly once"


def test_unknown_organism_defaults_to_non_bacterial(monkeypatch):
    monkeypatch.setattr(s5, "_get", lambda path, throttle=0.12: _FAKE_BRITE_TREE)
    assert s5.is_bacterial("zzz_not_a_real_code", {}) is False


def test_non_bacterial_kos_end_up_uncovered_not_silently_swapped(tmp_path, monkeypatch):
    """A KO whose only KEGG members are non-bacterial (e.g. a fungal-only KO like
    K10531/SID1) must be reported as uncovered, never quietly filled with a
    eukaryotic representative."""
    monkeypatch.setattr(s5, "genes_for_ko", lambda ko: ["nte:g1", "pan:g2"])
    monkeypatch.setattr(s5, "fetch_aaseq", lambda ids, batch=10: ">should-not-be-called\nMKV\n")
    monkeypatch.setattr(s5, "is_bacterial", lambda org, cache: False)
    df = s5.run(tmp_path, kos=["K10531"], gaps_only=False)
    assert df.empty, "no bacterial representative exists, so nothing should be fetched"
    unc = pd.read_csv(tmp_path / "ko_still_uncovered.csv")
    assert "K10531" in set(unc.KO)
    assert "non-bacterial" in unc.set_index("KO").loc["K10531", "reason"]
    rej = pd.read_csv(tmp_path / "rejected_non_bacterial_kegg.csv")
    assert set(rej.organism_code) == {"nte", "pan"}


def test_bacterial_organisms_pass_the_filter_unmolested(tmp_path, monkeypatch):
    monkeypatch.setattr(s5, "genes_for_ko", lambda ko: ["eco:b1"])
    monkeypatch.setattr(s5, "fetch_aaseq", lambda ids, batch=10: ">eco:b1\nMKV\n")
    monkeypatch.setattr(s5, "is_bacterial", lambda org, cache: True)
    df = s5.run(tmp_path, kos=["K00216"], gaps_only=False)
    assert set(df.organism_code) == {"eco"}
    rej = pd.read_csv(tmp_path / "rejected_non_bacterial_kegg.csv")
    assert rej.empty
