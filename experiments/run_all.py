"""Run every experiment, regenerate all figures and result JSON files, and write docs/results_summary.md."""
import importlib
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))
from common import DATA_DIR, ROOT, record  # noqa: E402

MODULES = [
    ("exp01_geometry", "exp01_geometry"),
    ("exp02_infall_tides", "exp02_infall_tides"),
    ("exp03_bridge_causality", "exp03_bridge_causality"),
    ("exp04_regular_bh", "exp04_regular_bh"),
    ("exp05_holography", "exp05_holography"),
    ("exp06_branching", "exp06_branching"),
    ("exp07_time_reversal", "exp07_time_reversal"),
    ("exp08_wormholes_chronology", "exp08_wormholes_chronology"),
    ("exp09_qnm_multiverse", "exp09_qnm_multiverse"),
    ("exp10_dependency_graph", "exp10_dependency_graph"),
]


def main():
    summary = {}
    for mod, name in MODULES:
        t0 = time.time()
        res = importlib.import_module(mod).run()
        record(name, res)
        summary[name] = round(time.time() - t0, 1)
        print(f"{name}: {summary[name]} s")
    with open(os.path.join(DATA_DIR, "run_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
