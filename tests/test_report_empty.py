"""Regression: report must survive header-less empty stage outputs."""
import pandas as pd
from bfgm.report import build


def test_build_skips_empty_files(tmp_path):
    (tmp_path / "ko_coverage_gaps.csv").write_text("")          # zero bytes
    (tmp_path / "discovered_kos.csv").write_text("KO,ko_definition\n")  # header only
    pd.DataFrame([{"gene": "nanH", "kegg_ko": "K23550", "status_code": "UNIQUE"}]) \
        .to_csv(tmp_path / "gene_ko_map.csv", index=False)
    out = build(tmp_path)
    assert out.exists()
