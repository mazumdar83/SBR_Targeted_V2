"""Symbol matching must be exact-token. These guard the main historical failure mode."""
import pytest
from bfgm.clients.kegg import KeggClient


class FakeKegg(KeggClient):
    def __init__(self):
        super().__init__(cache_dir="/tmp/bfgm-test-cache")
        self._ko_sym = {
            "K00216": ["entA"],
            "K11059": ["sea", "entA"],
            "K03894": ["iucA", "iucD"],
            "K10255": ["desA", "FAD2"],
            "K28682": ["desA", "bibA"],
        }
        self._ko_def = {
            "K00216": "2,3-dihydro-2,3-dihydroxybenzoate dehydrogenase [EC:1.3.1.28]",
            "K11059": "enterotoxin type A",
            "K03894": "N(2)-citryl-N(6)-acetyl-N(6)-hydroxylysine synthase",
            "K10255": "acyl-lipid omega-6 desaturase [EC:1.14.19.-]",
            "K28682": "lysine decarboxylase [EC:4.1.1.18]",
        }
        self._sym_index = None


def test_exact_match_not_substring():
    k = FakeKegg()
    # "ent" must not match "entA"; substring matching was the historical bug
    assert k.match_symbol("ent") == []
    assert set(k.match_symbol("entA")) == {"K00216", "K11059"}


def test_case_insensitive():
    k = FakeKegg()
    assert set(k.match_symbol("ENTA")) == set(k.match_symbol("entA"))


def test_collision_is_reported_not_silently_resolved():
    k = FakeKegg()
    assert len(k.match_symbol("desA")) == 2, "desA collides; both hits must surface"


def test_unknown_symbol_returns_empty():
    assert FakeKegg().match_symbol("notarealgene") == []
