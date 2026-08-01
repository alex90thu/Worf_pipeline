import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from plots.plot_by_n_s import build_control_outlier_report


def test_build_control_outlier_report_sorts_and_exposes_fields():
    df = pd.DataFrame(
        {
            "target_id": [1, 2, 3],
            "gene_name": ["A", "B", "C"],
            "on_target": [True, False, False],
            "read_count_ctrl": [157, 12, 8],
            "read_count_sample": [10, 10, 10],
        }
    )

    report = build_control_outlier_report(df, topk=2, min_ctrl=10)

    assert len(report) == 2
    assert report.iloc[0]["gene_name"] == "A"
    assert report.iloc[0]["read_count_ctrl"] == 157
    assert report.iloc[1]["gene_name"] == "B"
    assert bool(report.iloc[1]["on_target"]) is False
