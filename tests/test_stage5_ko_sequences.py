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
    s5.run(tmp_path, kos=["K00216", "K99999"], gaps_only=False)
    unc = pd.read_csv(tmp_path / "ko_still_uncovered.csv")
    assert "K99999" in set(unc.KO), "a KO with no members must be logged, not silently dropped"
