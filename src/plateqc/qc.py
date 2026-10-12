"""Plate loading and QC orchestration.

Two input shapes are supported, auto-detected:

* Grid CSV: a header row of column numbers (1..N) and row-labeled data
  rows (A..P), the classic plate-reader matrix export.
* Tidy CSV: one row per well with a well-id column (``well``, ``well_id``,
  ``Well``) and a numeric value column.

Control layouts map group names to lists of well/range specs, e.g.
``{"high": ["A1:A12"], "low": ["H1:H12"]}``. ``run_qc`` returns a report
dict with per-group stats, assay metrics, spatial diagnostics, and any
warnings (missing wells, non-numeric cells, too few controls).
"""

import csv
import io
import json
import re

from .metrics import (
    mean, stdev, cv, zprime, ssmd, signal_window, s2b,
    edge_ratio, drift, row_cv, col_cv,
)
from .wells import parse_well, expand_layout, plate_dims, edge_interior_split

_WELL_HEADER_RE = re.compile(r"^(well|well_id)$", re.IGNORECASE)
_MISSING_TOKENS = {"", "na", "n/a", "null", "none", "-", "--"}


def _parse_number(text):
    text = text.strip()
    if text.lower() in _MISSING_TOKENS:
        return None
    try:
        return float(text.replace(",", ""))
    except ValueError:
        return None


def _looks_like_grid_header(row):
    """A grid header row is mostly integer column numbers."""
    cells = [c.strip() for c in row if c.strip()]
    if len(cells) < 2:
        return False
    ints = 0
    for c in cells:
        try:
            int(c)
            ints += 1
        except ValueError:
            pass
    return ints >= len(cells) / 2


def _load_grid(rows, path):
    header_idx = None
    for i, row in enumerate(rows):
        if _looks_like_grid_header(row):
            header_idx = i
            break
    if header_idx is None:
        raise ValueError(f"{path}: no grid header row of column numbers found")
    header = [c.strip() for c in rows[header_idx]]
    # find first column whose header is an integer -> column-number offset
    col_offset = 0
    for j, c in enumerate(header):
        try:
            int(c)
            col_offset = j
            break
        except ValueError:
            continue
    n_cols = sum(1 for c in header[col_offset:] if c.strip())
    wells = {}
    warnings = []
    row_label = None
    expected_cols = None
    for row in rows[header_idx + 1:]:
        if not any(c.strip() for c in row):
            continue
        label = row[col_offset - 1].strip() if col_offset > 0 else row[0].strip()
        if not re.fullmatch(r"[A-Za-z]+", label):
            # footer / note row; stop reading data
            break
        row_label = label
        values = row[col_offset:col_offset + n_cols]
        for j, cell in enumerate(values):
            v = _parse_number(cell)
            if v is None and cell.strip():
                warnings.append(f"{path}: non-numeric cell at {label}{j + 1}: {cell.strip()!r}")
            wells[f"{label.upper()}{j + 1}"] = v
    if not wells:
        raise ValueError(f"{path}: no data rows found after the header row")
    return wells, warnings


def _load_tidy(rows, path, value_column=None):
    header = [c.strip() for c in rows[0]]
    well_idx = None
    for j, c in enumerate(header):
        if _WELL_HEADER_RE.match(c):
            well_idx = j
            break
    if well_idx is None:
        raise ValueError(f"{path}: no well-id column found (expected 'well' or 'well_id')")
    if value_column is None:
        # first numeric-ish column that is not the well column
        for j, c in enumerate(header):
            if j == well_idx:
                continue
            sample = [r[j] for r in rows[1:6] if len(r) > j]
            if sample and all(_parse_number(s) is not None for s in sample):
                value_column = c
                break
    if value_column is None:
        raise ValueError(f"{path}: could not find a numeric value column")
    if value_column not in header:
        raise ValueError(f"{path}: value column {value_column!r} not in header")
    val_idx = header.index(value_column)
    wells = {}
    warnings = []
    for row in rows[1:]:
        if len(row) <= max(well_idx, val_idx) or not row[well_idx].strip():
            continue
        wid = row[well_idx].strip().upper()
        try:
            parse_well(wid)
        except ValueError:
            warnings.append(f"{path}: skipping bad well id {row[well_idx].strip()!r}")
            continue
        v = _parse_number(row[val_idx])
        if v is None and row[val_idx].strip():
            warnings.append(f"{path}: non-numeric value for {wid}: {row[val_idx].strip()!r}")
        wells[wid] = v
    if not wells:
        raise ValueError(f"{path}: no well data parsed")
    return wells, warnings, value_column


def load_plate(path, value_column=None):
    """Load a plate CSV. Returns (wells, info) where info describes format."""
    with open(path, newline="") as f:
        rows = [row for row in csv.reader(f)]
    rows = [r for r in rows if any(c.strip() for c in r)]
    if not rows:
        raise ValueError(f"{path}: file is empty")
    header = [c.strip() for c in rows[0]]
    has_well_col = any(_WELL_HEADER_RE.match(c) for c in header)
    if has_well_col:
        wells, warnings, col = _load_tidy(rows, path, value_column)
        info = {"format": "tidy", "value_column": col, "n_wells": len(wells)}
    elif _looks_like_grid_header(header) or any(
        _looks_like_grid_header(r) for r in rows[:4]
    ):
        wells, warnings = _load_grid(rows, path)
        info = {"format": "grid", "n_wells": len(wells)}
    else:
        raise ValueError(
            f"{path}: unrecognized format (need a grid matrix or a well column)"
        )
    info["warnings"] = warnings
    return wells, info


def _group_stats(values):
    return {
        "n": len(values),
        "mean": mean(values),
        "stdev": stdev(values),
        "cv": cv(values),
        "min": min(values) if values else None,
        "max": max(values) if values else None,
    }


def run_qc(wells, layout, min_controls=4):
    """Compute QC report for well->value mapping against a control layout.

    ``layout``: ``{group: [specs...]}`` as in expand_layout. Values may be
    None (missing cells); they are excluded from stats and reported.
    """
    warnings = []
    groups = expand_layout(layout)
    values = {g: [] for g in groups}
    missing = {g: [] for g in groups}
    for group, wids in groups.items():
        for w in wids:
            v = wells.get(w)
            if v is None:
                missing[group].append(w)
            else:
                values[group].append(v)
        if missing[group]:
            warnings.append(
                f"group {group!r}: {len(missing[group])} control wells missing from plate"
            )
        if len(values[group]) < min_controls:
            warnings.append(
                f"group {group!r}: only {len(values[group])} control values (min {min_controls})"
            )
    overlap = set()
    seen = {}
    for group, wids in groups.items():
        for w in wids:
            if w in seen:
                overlap.add(w)
            seen[w] = group
    if overlap:
        warnings.append(
            f"wells in more than one group: {', '.join(sorted(overlap))}"
        )

    report = {"groups": {g: _group_stats(v) for g, v in values.items()},
              "missing": {g: m for g, m in missing.items() if m},
              "warnings": warnings}
    if "high" in values and "low" in values:
        h, l = values["high"], values["low"]
        report["assay"] = {
            "zprime": zprime(h, l),
            "ssmd": ssmd(h, l),
            "signal_window": signal_window(h, l),
            "signal_to_background": s2b(h, l),
            "zprime_grade": _zgrade(zprime(h, l)),
        }

    # Spatial diagnostics run on assay wells only: control wells are excluded
    # so the control signal range does not masquerade as plate drift.
    control_wells = set()
    for wids in groups.values():
        control_wells.update(wids)
    assay_wells = {w: v for w, v in wells.items() if w not in control_wells}
    spatial_wells = assay_wells if assay_wells else wells
    dims = plate_dims(spatial_wells)
    if dims:
        rows, cols = dims
        all_wells = sorted(spatial_wells, key=parse_well)
        edge, interior = edge_interior_split(all_wells, rows, cols)
        edge_v = [spatial_wells[w] for w in edge if spatial_wells[w] is not None]
        interior_v = [spatial_wells[w] for w in interior if spatial_wells[w] is not None]
        by_row = {}
        by_col = {}
        for w in all_wells:
            v = spatial_wells[w]
            if v is None:
                continue
            r, c = parse_well(w)
            by_row.setdefault(r, []).append(v)
            by_col.setdefault(c, []).append(v)
        col_means = [mean(by_col[c]) for c in sorted(by_col)]
        report["spatial"] = {
            "rows": rows,
            "cols": cols,
            "edge_ratio": edge_ratio(edge_v, interior_v),
            "column_drift": drift(sorted(by_col), [mean(by_col[c]) for c in sorted(by_col)]),
            "row_cv": row_cv([mean(v) for v in by_row.values()]),
            "col_cv": col_cv([mean(v) for v in by_col.values()]),
        }
    return report


def _zgrade(z):
    if z is None:
        return "unknown"
    if z >= 0.5:
        return "excellent"
    if z >= 0:
        return "marginal"
    return "poor"


def qc_report(report):
    """Render a human-readable text summary of a run_qc report."""
    lines = []
    for group, stats in report["groups"].items():
        lines.append(f"[{group}] n={stats['n']}")
        for key in ("mean", "stdev", "cv", "min", "max"):
            v = stats[key]
            lines.append(f"  {key}: {'n/a' if v is None else round(v, 4)}")
    assay = report.get("assay")
    if assay:
        lines.append("[assay metrics]")
        for key in ("zprime", "ssmd", "signal_window", "signal_to_background"):
            v = assay[key]
            lines.append(f"  {key}: {'n/a' if v is None else round(v, 4)}")
        lines.append(f"  zprime_grade: {assay['zprime_grade']}")
    spatial = report.get("spatial")
    if spatial:
        lines.append(f"[spatial] {spatial['rows']}x{spatial['cols']} plate")
        for key in ("edge_ratio", "column_drift", "row_cv", "col_cv"):
            v = spatial[key]
            lines.append(f"  {key}: {'n/a' if v is None else round(v, 4)}")
    for w in report["warnings"]:
        lines.append(f"WARNING: {w}")
    return "\n".join(lines) + "\n"


def load_layout(path):
    with open(path) as f:
        layout = json.load(f)
    if not isinstance(layout, dict):
        raise ValueError(f"{path}: layout must be a JSON object of group -> [well specs]")
    return layout
