"""plateqc: microplate assay quality-control metrics from plate reader exports."""

from .metrics import mean, stdev, cv, zprime, ssmd, signal_window, s2b
from .wells import parse_well, expand_range, expand_layout, plate_dims, is_edge_well
from .qc import load_plate, run_qc, qc_report

__all__ = [
    "mean", "stdev", "cv", "zprime", "ssmd", "signal_window", "s2b",
    "parse_well", "expand_range", "expand_layout", "plate_dims", "is_edge_well",
    "load_plate", "run_qc", "qc_report",
]

__version__ = "0.1.0"
