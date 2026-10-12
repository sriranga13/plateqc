# plateqc

Microplate assay quality control from plate reader exports. Feed it a CSV straight out of the reader, point it at your control wells, and get the metrics that decide whether a plate is trustworthy: **Z'-factor**, **SSMD**, signal window, and spatial diagnostics (edge effects, column drift) that catch plate problems before they become bad screening data.

Zero dependencies. Stdlib only.

## Quick start

```bash
pip install plateqc   # (once published; until then: pip install -e .)
plateqc qc examples/control_plate_grid.csv --layout examples/layout.json
```

Example output:

```
plate: examples/control_plate_grid.csv (grid, 96 wells)
[high] n=8
  mean: 999.87
  stdev: 38.12
  cv: 0.0381
  min: 944.2
  max: 1056.3
[low] n=8
  ...
[assay metrics]
  zprime: 0.94
  ssmd: 28.31
  signal_window: 752.41
  signal_to_background: 10.02
  zprime_grade: excellent
[spatial] 8x12 plate
  edge_ratio: -0.0213
  column_drift: 0.0832
  row_cv: 0.0176
  col_cv: 0.4121
```

JSON output (`--format json`) is machine-readable for LIMS pipelines, and `--fail-below 0.5` turns it into a CI gate: exit 1 when the plate's Z'-factor is below your cutoff.

## Input formats

plateqc auto-detects the export shape:

* **Grid CSV** - the classic reader matrix: a header row of column numbers (`1..12`, possibly under an instrument preamble) and row-labeled data rows (`A..H`). Trailing footer/note rows are ignored.
* **Tidy CSV** - one row per well with a `well` (or `well_id`) column and a numeric value column (auto-detected; override with `--value-column`).

Missing cells (`""`, `NA`, `N/A`) and censored tokens (`ND`, `BQL`, `<`, `>`, `TNTC`) parse as missing and are excluded from statistics rather than silently coerced to zero.

## Control layouts

A layout is a small JSON file mapping group names to well or range specs:

```json
{
  "high": ["A1:H1"],
  "low": ["A12:H12"]
}
```

Ranges are rectangles between two corners, so `"A1:H12"` is the whole plate and `"H12:A1"` works identically. Any group names work, but metrics need `high` and `low` groups:

| Metric | Meaning | Rule of thumb |
|---|---|---|
| `zprime` | Z'-factor, Zhang et al. 1999 | > 0.5 excellent, 0-0.5 marginal, < 0 poor |
| `ssmd` | Strictly standardized mean difference | larger is better |
| `signal_window` | (mh - 3sh) - (ml + 3sl) | positive means separation |
| `signal_to_background` | mean(high)/mean(low) | assay-dependent |
| `edge_ratio` | mean(edge)/mean(interior) - 1 | near 0 is healthy |
| `column_drift` | fractional signal change across columns | near 0 is healthy |

Warnings are raised for missing control wells, too few controls (`--min-controls`, default 4), wells listed in more than one group, and non-numeric cells.

## Why this exists

Z' is the first number reviewers ask for on any plate-based screen, and it is usually computed by hand in a spreadsheet after copy-pasting reader exports. plateqc makes it a one-liner that reads the export directly, flags spatial artifacts (evaporation edge effects, dispense drift) that a single Z' hides, and gates pipelines on plate quality. It pairs well with [lims-etl](https://github.com/sriranga13/lims-etl) for normalizing the export first.

## License

MIT.
