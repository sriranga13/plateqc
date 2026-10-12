"""Well-id parsing, range expansion, and plate geometry helpers.

Well ids are the standard lab notation: row letter(s) + column number,
e.g. ``A1``, ``H12``, ``P24``. Ranges like ``A1:H1`` or ``A1:H12`` expand
to the full rectangle of wells between the two corners (inclusive).
"""

import re

_WELL_RE = re.compile(r"^([A-Za-z]+)(\d+)$")


def parse_well(well):
    """Return (row_index, col_index) for a well id like 'B07' (0-based)."""
    m = _WELL_RE.match(well.strip())
    if not m:
        raise ValueError(f"not a well id: {well!r}")
    row = 0
    for ch in m.group(1).upper():
        if not ("A" <= ch <= "Z"):
            raise ValueError(f"not a well id: {well!r}")
        row = row * 26 + (ord(ch) - ord("A") + 1)
    row -= 1
    col = int(m.group(2)) - 1
    if row < 0 or col < 0:
        raise ValueError(f"not a well id: {well!r}")
    return row, col


def well_id(row, col):
    """Return the canonical well id for 0-based (row, col)."""
    if row < 0 or col < 0:
        raise ValueError("row/col must be >= 0")
    letters = ""
    r = row + 1
    while r:
        r, rem = divmod(r - 1, 26)
        letters = chr(ord("A") + rem) + letters
    return f"{letters}{col + 1}"


def expand_range(spec):
    """Expand a well or range spec into a sorted list of canonical well ids.

    Accepts single wells (``"C4"``) and rectangular ranges (``"A1:H1"``,
    ``"A1:H12"``). The two corners of a range define the rectangle,
    so ``"H12:A1"`` works the same as ``"A1:H12"``.
    """
    spec = spec.strip()
    if ":" not in spec:
        row, col = parse_well(spec)
        return [well_id(row, col)]
    start, _, end = spec.partition(":")
    r1, c1 = parse_well(start)
    r2, c2 = parse_well(end)
    wells = [
        well_id(r, c)
        for r in range(min(r1, r2), max(r1, r2) + 1)
        for c in range(min(c1, c2), max(c1, c2) + 1)
    ]
    return sorted(wells, key=lambda w: parse_well(w))


def expand_layout(layout):
    """Expand a ``{group_name: [specs...]}`` layout mapping to well-id sets.

    Duplicate wells across specs are deduplicated; across groups they are
    kept (a well can legally be both a control and a sample in QC math,
    but run_qc flags overlaps as warnings).
    """
    out = {}
    for group, specs in layout.items():
        if isinstance(specs, str):
            specs = [specs]
        wells = []
        for spec in specs:
            wells.extend(expand_range(spec))
        out[group] = sorted(set(wells), key=lambda w: parse_well(w))
    return out


def plate_dims(wells):
    """Infer (rows, cols) from a collection of well ids; None if unknown."""
    rows = cols = 0
    for w in wells:
        r, c = parse_well(w)
        rows = max(rows, r + 1)
        cols = max(cols, c + 1)
    return (rows, cols) if rows and cols else None


def is_edge_well(well, rows, cols):
    """True when the well sits on the plate's outer rim."""
    r, c = parse_well(well)
    return r == 0 or c == 0 or r == rows - 1 or c == cols - 1


def edge_interior_split(wells, rows, cols):
    """Split well ids into (edge_wells, interior_wells)."""
    edge, interior = [], []
    for w in wells:
        (edge if is_edge_well(w, rows, cols) else interior).append(w)
    return edge, interior
