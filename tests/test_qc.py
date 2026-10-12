import json

import pytest

from plateqc.qc import load_plate, run_qc, qc_report, load_layout

GRID = "examples/control_plate_grid.csv"
TIDY = "examples/control_plate_tidy.csv"
LAYOUT = {"high": ["A1:H1"], "low": ["A12:H12"]}


def test_load_grid():
    wells, info = load_plate(GRID)
    assert info["format"] == "grid"
    assert info["n_wells"] == 96
    assert wells["A1"] is not None
    assert wells["H12"] is not None


def test_load_tidy():
    wells, info = load_plate(TIDY)
    assert info["format"] == "tidy"
    assert info["value_column"] == "channel_A"
    assert info["n_wells"] == 96
    # missing cell and censored token parse as None
    assert wells["A10"] is None
    assert wells["G6"] is None


def test_load_tidy_explicit_column():
    wells, info = load_plate(TIDY, value_column="channel_A")
    assert info["n_wells"] == 96


def test_load_rejects_empty(tmp_path):
    p = tmp_path / "empty.csv"
    p.write_text("\n")
    with pytest.raises(ValueError):
        load_plate(str(p))


def test_load_rejects_garbage(tmp_path):
    p = tmp_path / "junk.csv"
    p.write_text("hello,world\nfoo,bar\n")
    with pytest.raises(ValueError):
        load_plate(str(p))


def test_run_qc_assay_metrics():
    wells, _ = load_plate(GRID)
    report = run_qc(wells, LAYOUT)
    assay = report["assay"]
    assert assay["zprime"] > 0.5
    assert assay["zprime_grade"] == "excellent"
    assert assay["ssmd"] > 3
    assert assay["signal_window"] > 0
    assert report["groups"]["high"]["n"] == 8
    assert report["groups"]["low"]["n"] == 8


def test_run_qc_missing_controls_warn():
    wells, _ = load_plate(GRID)
    report = run_qc(wells, {"high": ["A1:H1"], "low": ["A12:H12", "Z99"]},
                    min_controls=8)
    assert any("missing" in w for w in report["warnings"])


def test_run_qc_too_few_controls_warn():
    wells, _ = load_plate(GRID)
    report = run_qc(wells, {"high": ["A1"], "low": ["H12"]})
    assert any("only 1 control values" in w for w in report["warnings"])
    # metrics degrade gracefully
    assert report["assay"]["zprime"] is None


def test_run_qc_overlap_warn():
    wells, _ = load_plate(GRID)
    report = run_qc(wells, {"high": ["A1:H1"], "low": ["H1:H2"]})
    assert any("more than one group" in w for w in report["warnings"])


def test_run_qc_spatial():
    wells, _ = load_plate(GRID)
    report = run_qc(wells, LAYOUT)
    spatial = report["spatial"]
    # control columns (1 and 12) are excluded from spatial diagnostics,
    # so the assay region spans 8x11 by well coordinates
    assert spatial["rows"] == 8 and spatial["cols"] == 11
    assert spatial["edge_ratio"] is not None
    assert abs(spatial["column_drift"]) < 0.5


def test_qc_report_text():
    wells, _ = load_plate(GRID)
    text = qc_report(run_qc(wells, LAYOUT))
    assert "[high]" in text
    assert "zprime" in text
    assert "excellent" in text


def test_load_layout(tmp_path):
    p = tmp_path / "layout.json"
    p.write_text(json.dumps({"high": ["A1:A12"]}))
    assert load_layout(str(p)) == {"high": ["A1:A12"]}
    p.write_text("[1,2]")
    with pytest.raises(ValueError):
        load_layout(str(p))
