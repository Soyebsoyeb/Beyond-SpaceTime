"""Kruskal and Penrose diagram drawing helpers (matplotlib)."""
from __future__ import annotations

import numpy as np

from .coordinates import kruskal_to_r, schwarzschild_to_kruskal, to_penrose
from .plotting import PALETTE

PI = np.pi
HEX = np.array([[PI, 0], [PI / 2, PI / 2], [-PI / 2, PI / 2], [-PI, 0], [-PI / 2, -PI / 2], [PI / 2, -PI / 2], [PI, 0]])


def draw_penrose_frame(ax, labels=True, regions=True):
    """Draw the outline of the Penrose diagram: boundaries, both horizons, both singularities, region labels."""
    ax.plot(HEX[:, 0], HEX[:, 1], color="black", lw=1.2)
    # horizons (diagonals through the bifurcation sphere)
    for s in (+1, -1):
        ax.plot([-PI / 2, PI / 2], [-s * PI / 2, s * PI / 2], color=PALETTE["blue"], lw=1.4)
    # singularities
    x = np.linspace(-PI / 2, PI / 2, 200)
    for s in (+1, -1):
        ax.plot(x, s * (PI / 2 + 0.04 * np.sin(24 * x)), color=PALETTE["red"], lw=2.0)
    if labels and regions:
        for txt, xy in {"I": (PI * 0.62, 0), "III": (-PI * 0.62, 0), "II": (0, PI * 0.28), "IV": (0, -PI * 0.28)}.items():
            ax.text(*xy, txt, ha="center", va="center", fontsize=13)
        ax.text(0, PI / 2 + 0.17, "future singularity r = 0", color=PALETTE["red"], ha="center", fontsize=9)
        ax.text(0, -PI / 2 - 0.27, "past singularity r = 0", color=PALETTE["red"], ha="center", fontsize=9)
        ax.text(PI + 0.08, 0.0, "i$^0$", fontsize=9, va="center")
        ax.text(-PI - 0.28, 0.0, "i$^0$", fontsize=9, va="center")
        ax.text(PI * 0.55, PI * 0.35, "$\\mathscr{I}^+$", fontsize=9)
        ax.text(PI * 0.55, -PI * 0.38, "$\\mathscr{I}^-$", fontsize=9)
    ax.set_aspect("equal")
    ax.set_xlim(-PI - 0.3, PI + 0.3)
    ax.set_ylim(-PI / 2 - 0.45, PI / 2 + 0.45)
    ax.set_xticks([]), ax.set_yticks([])
    ax.grid(False)


def penrose_curve(ax, T, X, **kw):
    """Plot a curve given in Kruskal (T, X) on a Penrose-diagram axis."""
    Tp, Xp = to_penrose(T, X)
    ax.plot(Xp, Tp, **kw)


def draw_constant_r_and_t(ax, M=1.0):
    """Constant-r curves (exterior and interior) and constant-t rays in Penrose coordinates."""
    for r in (2.6, 4.0, 8.0, 20.0):
        for region in ("I", "III"):
            tt = np.linspace(-6, 6, 400)
            T, X = schwarzschild_to_kruskal(tt, np.full_like(tt, r), M, region)
            penrose_curve(ax, T, X, color=PALETTE["grey"], lw=0.7, alpha=0.8)
    for r in (1.7, 1.2, 0.6):
        for region in ("II", "IV"):
            tt = np.linspace(-6, 6, 400)
            T, X = schwarzschild_to_kruskal(tt, np.full_like(tt, r), M, region)
            penrose_curve(ax, T, X, color=PALETTE["orange"], lw=0.7, alpha=0.8)
    for tval in np.linspace(-3.0, 3.0, 7) * 4 * M * 0.6:
        rr = 2.0 + np.logspace(-8, np.log10(400.0), 800)
        for region in ("I", "III"):
            T, X = schwarzschild_to_kruskal(np.full_like(rr, tval), rr, M, region)
            penrose_curve(ax, T, X, color=PALETTE["green"], lw=0.6, alpha=0.7)


def draw_kruskal(ax, M=1.0, lim=3.0):
    """Draw the Kruskal-Szekeres diagram with singularities, horizons, constant-r and constant-t curves."""
    X = np.linspace(-lim, lim, 400)
    for s in (+1, -1):
        ax.plot(X, s * np.sqrt(1 + X**2), color=PALETTE["red"], lw=2.0)
    ax.plot([-lim, lim], [-lim, lim], color=PALETTE["blue"], lw=1.2)
    ax.plot([-lim, lim], [lim, -lim], color=PALETTE["blue"], lw=1.2)
    for r in (2.4, 3.0, 5.0, 10.0):
        tt = np.linspace(-6, 6, 300)
        for region in ("I", "III"):
            T, Xk = schwarzschild_to_kruskal(tt, np.full_like(tt, r), M, region)
            ax.plot(Xk, T, color=PALETTE["grey"], lw=0.7)
    for r in (1.6, 1.0, 0.4):
        tt = np.linspace(-6, 6, 300)
        for region in ("II", "IV"):
            T, Xk = schwarzschild_to_kruskal(tt, np.full_like(tt, r), M, region)
            ax.plot(Xk, T, color=PALETTE["orange"], lw=0.7)
    for k in np.linspace(-4, 4, 9):
        ax.plot([-lim, lim], [np.tanh(k / 4) * -lim, np.tanh(k / 4) * lim], color=PALETTE["green"], lw=0.5, alpha=0.7)
    ax.set_xlim(-lim, lim), ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Kruskal X / (2M)")
    ax.set_ylabel("Kruskal T / (2M)")
    for txt, xy in {"I": (2.0, 0.0), "III": (-2.0, 0.0), "II": (0, 1.6), "IV": (0, -1.6)}.items():
        ax.text(*xy, txt, ha="center", va="center", fontsize=13)
