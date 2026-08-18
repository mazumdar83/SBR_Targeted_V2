"""Live API smoke tests. Run with: pytest -m live

These verify the empirical constraints the pipeline is built on. If one fails, an
upstream API changed and the pipeline logic needs revisiting.
"""
import pytest
from bfgm.clients.kegg import KeggClient
from bfgm.clients.uniprot import UniProtClient

pytestmark = pytest.mark.live


def test_kegg_ko_catalogue_loads():
    sym, dfn = KeggClient().load_ko_catalogue()
    assert len(sym) > 20000
    assert "K02014" in dfn


def test_uniprot_has_no_ko_xref():
    """The constraint that forces the two-hop anchor in stage 4."""
    assert UniProtClient().count("xref:ko-K02014") == 0, \
        "UniProt now exposes KO xrefs; stage 4 can be simplified"


def test_gene_exact_works():
    assert UniProtClient().count("gene_exact:feoB") > 0


def test_per_organism_ko_route():
    kos = KeggClient().organism_kos("eco")
    assert len(kos) > 1000
