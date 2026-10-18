"""Shared matplotlib style and helpers (Agg backend, colour-blind-safe palette)."""
from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PALETTE = {
    "blue": "#1f77b4",
    "orange": "#d95f02",
    "green": "#1b9e77",
    "purple": "#7570b3",
    "red": "#c0392b",
    "grey": "#555555",
    "gold": "#b8860b",
}
CYCLE = [PALETTE[k] for k in ("blue", "orange", "green", "purple", "red", "gold", "grey")]


def apply_style() -> None:
    """Apply the project-wide Matplotlib style (resolution, grid, colour cycle, constrained layout)."""
    plt.rcParams.update(
        {
            "figure.dpi": 110,
            "savefig.dpi": 170,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "axes.prop_cycle": matplotlib.cycler(color=CYCLE),
            "legend.frameon": False,
            "figure.constrained_layout.use": True,
            "mathtext.default": "regular",
        }
    )


apply_style()

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "figures")


def save(fig, name: str, fig_dir: str | None = None) -> str:
    """Save a figure as PNG in the figures directory (default results/figures) and close it; returns the path."""
    d = os.path.abspath(fig_dir or FIG_DIR)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, name if name.endswith(".png") else name + ".png")
    fig.savefig(path)
    plt.close(fig)
    return path
