"""Stage 6 must merge both sources, dedupe by sequence, and keep headers BLAST-safe."""
import pandas as pd
import pytest
from bfgm.stages import s6_reference_db as s6

SEQ_A = "MKRTLIAAALLSGLSAWSANAAE" * 3
SEQ_B = "MSDLNQAVIAGVEHLLGKDAVYA" * 3


@pytest.fixture
def run_dir(tmp_path):
    pd.DataFrame([
        {"Entry": "P00001", "Sequence": SEQ_A, "Gene Names": "entA",
         "Organism": "Escherichia coli", "Organism (ID)": 562,
         "Protein names": "DHB dehydrogenase", "kegg_KOs": "K00216",
         "matched_query_terms": "entA"},
        {"Entry": "P00002", "Sequence": SEQ_A, "Gene Names": "entA",   # exact duplicate
         "Organism": "Shigella flexneri", "Organism (ID)": 623,
         "Protein names": "DHB dehydrogenase", "kegg_KOs": "K00216",
         "matched_query_terms": "entA"},
    ]).to_csv(tmp_path / "proteins_anchored.tsv", sep="\t", index=False)
    (tmp_path / "ko_sequences.fasta").write_text(
        f">eco:b0150 K02014 iron complex receptor | (RefSeq) fhuA; ferrichrome\n{SEQ_B}\n")
    return tmp_path


def test_merges_both_sources(run_dir):
    st = s6.run(run_dir, make_blast_db=False, make_diamond_db=False)
    assert st["from_uniprot_only"] == 1
    assert st["from_kegg_only"] == 1


def test_identical_sequences_collapse(run_dir):
    st = s6.run(run_dir, make_blast_db=False, make_diamond_db=False)
    assert st["duplicates_collapsed"] == 1, "identical sequences must collapse to one entry"
    assert st["sequences"] == 2


def test_header_ids_have_no_whitespace(run_dir):
    s6.run(run_dir, make_blast_db=False, make_diamond_db=False)
    fa = (run_dir / f"{run_dir.name}_reference.fasta").read_text()
    for line in fa.splitlines():
        if line.startswith(">"):
            sid = line[1:].split(" ")[0]
            assert "|" in sid and len(sid.split("|")) == 5
            assert not any(c.isspace() for c in sid)


def test_ko_recoverable_from_seqid(run_dir):
    s6.run(run_dir, make_blast_db=False, make_diamond_db=False)
    m = pd.read_csv(run_dir / f"{run_dir.name}_reference_map.tsv", sep="\t")
    assert set(m.ko) == {"K00216", "K02014"}
    for _, r in m.iterrows():
        assert r.seq_id.split("|")[2] == r.ko, "KO in the ID must match the map"


def test_short_sequences_filtered(run_dir):
    st = s6.run(run_dir, make_blast_db=False, make_diamond_db=False, min_length=10_000)
    assert st["sequences"] == 0 or st["sequences"] < 2


def test_usage_doc_written(run_dir):
    s6.run(run_dir, make_blast_db=False, make_diamond_db=False)
    doc = (run_dir / "HOW_TO_SEARCH.md").read_text()
    assert "diamond blastp" in doc and "reciprocal best hits" in doc


def test_no_sequences_raises(tmp_path):
    with pytest.raises(RuntimeError):
        s6.run(tmp_path)
