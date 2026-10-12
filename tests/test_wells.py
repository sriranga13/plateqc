import pytest

from plateqc.wells import (
    parse_well, well_id, expand_range, expand_layout, plate_dims,
    is_edge_well, edge_interior_split,
)


def test_parse_well_basic():
    assert parse_well("A1") == (0, 0)
    assert parse_well("B07") == (1, 6)
    assert parse_well("H12") == (7, 11)
    assert parse_well("P24") == (15, 23)


def test_parse_well_rejects_junk():
    for bad in ("", "1A", "A0", "A", "12", "A-1", "A1B"):
        with pytest.raises(ValueError):
            parse_well(bad)


def test_well_id_roundtrip():
    assert well_id(0, 0) == "A1"
    assert well_id(7, 11) == "H12"
    assert well_id(*parse_well("P24")) == "P24"


def test_expand_range_single():
    assert expand_range("C4") == ["C4"]


def test_expand_range_row():
    assert expand_range("A1:A3") == ["A1", "A2", "A3"]


def test_expand_range_rectangle():
    got = expand_range("A1:B2")
    assert got == ["A1", "A2", "B1", "B2"]


def test_expand_range_reversed_corners():
    assert expand_range("H12:A1") == expand_range("A1:H12")


def test_expand_layout_groups():
    layout = expand_layout({"high": ["A1:H1"], "low": ["A12:H12"]})
    assert len(layout["high"]) == 8
    assert layout["high"][0] == "A1"
    assert layout["low"][-1] == "H12"


def test_expand_layout_dedupes():
    layout = expand_layout({"g": ["A1:A2", "A1"]})
    assert layout["g"] == ["A1", "A2"]


def test_plate_dims():
    assert plate_dims(["A1", "H12"]) == (8, 12)
    assert plate_dims([]) is None


def test_edge_wells():
    assert is_edge_well("A1", 8, 12)
    assert is_edge_well("H12", 8, 12)
    assert not is_edge_well("D6", 8, 12)


def test_edge_interior_split():
    edge, interior = edge_interior_split(["A1", "D6", "H12"], 8, 12)
    assert edge == ["A1", "H12"]
    assert interior == ["D6"]
