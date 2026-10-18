"""Structured comparison of four different meanings of "multiple worlds".

Each entry records (i) the mathematical object that supplies the multiplicity, (ii) whether the
multiplicity is a feature of the equations, an interpretive proposal, or a prediction that rests
on extra assumptions, and (iii) whether regions can exchange signals.  No claim is made that
Einstein-Rosen bridges connect many-worlds branches: no concrete model establishing this is known
to this repository, and none is implemented.
"""
from __future__ import annotations

EQUATIONS = "feature of the equations"
INTERPRETIVE = "interpretive proposal"
EXTRA = "prediction dependent on extra assumptions"

CONCEPTS = [
    {
        "key": "mwi",
        "name": "Quantum branching (many-worlds)",
        "object": "Orthogonal components of a state in a tensor-product Hilbert space",
        "status": INTERPRETIVE,
        "status_detail": "Unitary evolution (equations) yields entangled components; calling them worlds is interpretive.",
        "setting": "Non-relativistic or QFT Hilbert space; no extra spacetime regions",
        "signalling": "None between decohered components (no-signalling for unitary evolution on subsystems)",
        "implemented_in": "quantum.py, tfd.py",
    },
    {
        "key": "inflation",
        "name": "Eternal inflation / pocket universes",
        "object": "Disconnected regions of a spacetime with a metastable de Sitter vacuum and bubble nucleation",
        "status": EXTRA,
        "status_detail": "Requires a specific inflaton potential or landscape and a choice of measure for probabilities.",
        "setting": "Semiclassical gravity coupled to a scalar field with false-vacuum decay",
        "signalling": "Bubbles are causally disconnected by cosmological horizons (generically)",
        "implemented_in": "documentation only (no model implemented)",
    },
    {
        "key": "kruskal",
        "name": "Second exterior of maximally extended Schwarzschild",
        "object": "Region III of the analytic Kruskal extension",
        "status": EQUATIONS,
        "status_detail": "Exact vacuum solution; absent in collapse spacetimes, so not generic for astrophysical black holes.",
        "setting": "Eternal, exact vacuum, analytically extended spacetime",
        "signalling": "Not causally connected to region I (no future-directed causal curve I to III)",
        "implemented_in": "coordinates.py, causal.py",
    },
    {
        "key": "bounce",
        "name": "Regions beyond black-hole bounce / regular cores",
        "object": "Extension through a regular core or black-to-white-hole transition in an effective model",
        "status": EXTRA,
        "status_detail": "Depends on the chosen effective dynamics and on inner-horizon stability, which is contested.",
        "setting": "Effective metrics such as Hayward; tunnelling proposals",
        "signalling": "Model dependent; inner-horizon instability can invalidate the extension",
        "implemented_in": "regular.py (diagnostics only)",
    },
]


def markdown_table() -> str:
    """Return the comparison of the four multiverse concepts as a Markdown table string."""
    head = "| Concept | Mathematical object | Status of the multiplicity | Signalling between regions |\n"
    sep = "| --- | --- | --- | --- |\n"
    rows = "".join(f"| {c['name']} | {c['object']} | {c['status']} | {c['signalling']} |\n" for c in CONCEPTS)
    return head + sep + rows
