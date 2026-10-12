import json

import pytest

from plateqc.cli import main

GRID = "examples/control_plate_grid.csv"
TIDY = "examples/control_plate_tidy.csv"
LAYOUT = "examples/layout.json"


def test_cli_text(capsys):
    assert main(["qc", GRID, "--layout", LAYOUT]) == 0
    out = capsys.readouterr().out
    assert "plate:" in out
    assert "zprime" in out


def test_cli_tidy(capsys):
    assert main(["qc", TIDY, "--layout", LAYOUT]) == 0
    out = capsys.readouterr().out
    assert "tidy" in out


def test_cli_json(capsys):
    assert main(["qc", GRID, "--layout", LAYOUT, "--format", "json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["assay"]["zprime"] > 0.5
    assert report["groups"]["high"]["n"] == 8


def test_cli_fail_below():
    # grid z' is ~0.9, so failing below 0.99 exits 1
    assert main(["qc", GRID, "--layout", LAYOUT, "--fail-below", "0.99"]) == 1
    assert main(["qc", GRID, "--layout", LAYOUT, "--fail-below", "0.1"]) == 0


def test_cli_bad_plate():
    assert main(["qc", "nope.csv", "--layout", LAYOUT]) == 2


def test_cli_bad_layout(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{oops")
    assert main(["qc", GRID, "--layout", str(p)]) == 2


def test_cli_missing_args():
    with pytest.raises(SystemExit):
        main(["qc", GRID])
