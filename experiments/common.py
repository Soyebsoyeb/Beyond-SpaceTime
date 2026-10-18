"""Shared experiment utilities: output directories and JSON result recording."""
import json
import os

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIG_DIR = os.path.join(ROOT, "results", "figures")
DATA_DIR = os.path.join(ROOT, "results", "data")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return float(o)
    if isinstance(o, (np.integer, int)):
        return int(o)
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    if isinstance(o, complex):
        return {"re": o.real, "im": o.imag}
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    return o


def record(name, data):
    path = os.path.join(DATA_DIR, f"{name}.json")
    with open(path, "w") as fh:
        json.dump(_clean(data), fh, indent=2, sort_keys=True)
    return path
