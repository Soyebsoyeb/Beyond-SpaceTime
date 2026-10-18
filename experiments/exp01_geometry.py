"""Experiment 1: Schwarzschild geometry, curvature, Kruskal and Penrose diagrams (Part III-IV)."""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

from beyond_spacetime import coordinates as C
from beyond_spacetime import diagrams as D
from beyond_spacetime import symbolic as S
from beyond_spacetime.plotting import PALETTE, save

from common import FIG_DIR, record


def run():
    out = {}
    sch = S.schwarzschild()
    K = sch.kretschmann()
    out["kretschmann_symbolic"] = str(K)
    out["ricci_vanishes"] = bool(sch.ricci() == sp.zeros(4, 4))

    # Figure 1: metric function and curvature invariant
    r = np.linspace(0.05, 10, 2000)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(r, C.f(r), color=PALETTE["blue"])
    ax[0].axhline(0, color="black", lw=0.8)
    ax[0].axvline(2, color=PALETTE["red"], ls="--", lw=1)
    ax[0].set_xlabel("r / M"), ax[0].set_ylabel("f(r) = 1 - 2M/r")
    ax[0].set_ylim(-12, 1.2)
    ax[0].set_title("Metric function: sign change at r = 2M")
    ax[0].annotate("horizon", (2, 0), xytext=(3.3, -3), arrowprops=dict(arrowstyle="->"))
    rr = np.logspace(-1, 1.5, 300)
    Kn = sp.lambdify((S.r, S.M), K, "numpy")(rr, 1.0)
    ax[1].loglog(rr, Kn, color=PALETTE["orange"], label="K = 48 M$^2$/r$^6$ (symbolic)")
    ax[1].loglog(rr, 48 / rr**6, "k:", lw=1.2, label="closed form")
    ax[1].axvline(2, color=PALETTE["red"], ls="--", lw=1)
    ax[1].set_xlabel("r / M"), ax[1].set_ylabel("Kretschmann scalar  K M$^4$")
    ax[1].set_title("Curvature is finite at r = 2M, divergent at r = 0")
    ax[1].legend()
    save(fig, "fig01_metric_and_curvature", FIG_DIR)
    out["K_at_horizon_M1"] = float(48 / 2.0**6)

    # Figure 2: Kruskal and Penrose
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    D.draw_kruskal(ax[0])
    ax[0].set_title("Kruskal-Szekeres diagram (grey: constant r outside, orange: inside, green: constant t)")
    D.draw_penrose_frame(ax[1])
    D.draw_constant_r_and_t(ax[1])
    ax[1].set_title("Penrose diagram of the maximally extended solution")
    save(fig, "fig02_kruskal_penrose", FIG_DIR)

    # Removability of the coordinate singularity: EF v is finite while t diverges for an infalling ray
    out["tortoise_divergence_check"] = float(C.tortoise(2.0 * (1 + 1e-12)))
    return out


if __name__ == "__main__":
    print(record("exp01_geometry", run()))
