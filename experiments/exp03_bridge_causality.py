"""Experiment 3: Einstein-Rosen bridge, slicings, and causal structure (Part IV)."""
import matplotlib.pyplot as plt
import numpy as np

from beyond_spacetime import causal as CA
from beyond_spacetime import coordinates as C
from beyond_spacetime import diagrams as D
from beyond_spacetime import embedding as E
from beyond_spacetime.plotting import PALETTE, save

from common import FIG_DIR, record

PI = np.pi


def run():
    out = {}
    # ---- Figure 8: Flamm paraboloid (3D) ----
    r = np.linspace(2.0, 8.0, 120)
    phi = np.linspace(0, 2 * PI, 120)
    R, P = np.meshgrid(r, phi)
    Z = E.flamm_z(R)
    X, Y = R * np.cos(P), R * np.sin(P)
    fig = plt.figure(figsize=(11, 4.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.plot_surface(X, Y, Z, color=PALETTE["blue"], alpha=0.6, linewidth=0)
    ax.plot_surface(X, Y, -Z, color=PALETTE["orange"], alpha=0.6, linewidth=0)
    th = np.linspace(0, 2 * PI, 200)
    ax.plot(2 * np.cos(th), 2 * np.sin(th), 0 * th, color="black", lw=2)
    ax.set_title("Flamm paraboloid: two exteriors joined at r = 2M")
    ax.set_xlabel("x / M"), ax.set_ylabel("y / M"), ax.set_zlabel("z / M")
    ax.view_init(elev=14, azim=-60)
    ax2 = fig.add_subplot(1, 2, 2)
    rr = np.linspace(2.0, 10, 400)
    ax2.plot(rr, E.flamm_z(rr), color=PALETTE["blue"], label="z(r) = 2 sqrt(2M(r - 2M))")
    ax2.plot(rr, -E.flamm_z(rr), color=PALETTE["orange"])
    ax2.set_xlim(1.0, 10.5)
    ax2.set_xlabel("r / M"), ax2.set_ylabel("z / M")
    ax2.set_title("Meridian section; throat at r = 2M (area 16 $\\pi$ M$^2$)")
    ax2.axvline(2.0, color=PALETTE["red"], ls="--", lw=1)
    ax2.legend()
    save(fig, "fig08_flamm_paraboloid", FIG_DIR)

    # ---- Figure 9: slicing dependence ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
    for T0, c in zip((0.0, 0.5, 0.8, 0.95), (PALETTE["blue"], PALETTE["green"], PALETTE["orange"], PALETTE["purple"])):
        Xk, rk, lk = E.kruskal_slice_profile(T0, Xmax=3.5)
        ax[0].plot(lk, rk, color=c, label=f"Kruskal T = {T0}")
    rpg, lpg = E.painleve_gullstrand_profile(rmax=8)
    ax[0].plot(lpg, rpg, "k--", label="Painleve-Gullstrand (flat)")
    ax[0].plot(-lpg, rpg, "k--", alpha=0.0)
    ax[0].axhline(2, color=PALETTE["red"], ls=":", lw=1)
    ax[0].set_xlim(-6, 6), ax[0].set_ylim(0, 9)
    ax[0].set_xlabel("proper distance along slice / M"), ax[0].set_ylabel("areal radius r / M")
    ax[0].set_title("Areal radius on different spacelike slices")
    ax[0].legend(fontsize=8)
    T0s = np.linspace(0, 0.999, 200)
    rth = np.array([E.throat_radius_kruskal(t) for t in T0s])
    ax[1].plot(T0s, rth, color=PALETTE["blue"])
    ax[1].axhline(2, color=PALETTE["red"], ls=":", lw=1)
    ax[1].set_xlabel("Kruskal slice time T / (2M)"), ax[1].set_ylabel("minimal areal radius / M")
    ax[1].set_title("The throat radius is slicing dependent; it pinches to 0 at T = 1")
    save(fig, "fig09_slicings", FIG_DIR)
    out["throat_radius"] = {str(t): E.throat_radius_kruskal(t) for t in (0.0, 0.5, 0.8, 0.95)}

    # ---- Figure 10: causal structure ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    # (a) future of a point in region I
    D.draw_penrose_frame(ax[0], labels=False, regions=False)
    pu, pv = -0.9, 0.9          # point p in region I: (u~, v~) with u~ < 0 < v~
    g = np.linspace(-PI / 2, PI / 2, 700)
    UU, VV = np.meshgrid(g, g)
    mask = (UU >= pu) & (VV >= pv) & (UU + VV <= PI / 2)
    Tp, Xp = VV + UU, VV - UU
    ax[0].scatter(Xp[mask], Tp[mask], s=1.0, color=PALETTE["orange"], alpha=0.35, linewidths=0)
    ax[0].plot([(pv - pu)], [(pv + pu)], "ko")
    ax[0].text(pv - pu + 0.1, pv + pu - 0.25, "p", fontsize=11)
    for txt, xy in {"I": (PI * 0.62, -0.5), "III": (-PI * 0.62, 0), "II": (0, PI * 0.28), "IV": (0, -PI * 0.28)}.items():
        ax[0].text(*xy, txt, ha="center", va="center", fontsize=13)
    ax[0].set_title("(a) Causal future of p in region I never reaches III")
    # (b) collapse of a null shell onto a flat interior: only regions I and II survive
    v0 = 0.6
    A = (v0 + PI / 2, v0 - PI / 2)          # shell meets scri^- (Xp, Tp)
    Bp = (PI, 0.0)                          # i^0
    Cp = (PI / 2, PI / 2)                   # i^+
    Ep = (2 * v0 - PI / 2, PI / 2)          # shell reaches r = 0 (singularity forms)
    ax[1].plot(D.HEX[:, 0], D.HEX[:, 1], color="black", lw=0.6, alpha=0.25)
    ax[1].fill([A[0], Bp[0], Cp[0], Ep[0]], [A[1], Bp[1], Cp[1], Ep[1]], color=PALETTE["blue"], alpha=0.15)
    ax[1].plot([A[0], Ep[0]], [A[1], Ep[1]], color=PALETTE["purple"], lw=2.2, label="null shell v = v0")
    ax[1].plot([A[0], Bp[0]], [A[1], Bp[1]], color="black", lw=1.2)
    ax[1].plot([Bp[0], Cp[0]], [Bp[1], Cp[1]], color="black", lw=1.2)
    xs = np.linspace(Ep[0], Cp[0], 120)
    ax[1].plot(xs, PI / 2 + 0.04 * np.sin(24 * xs), color=PALETTE["red"], lw=2.0, label="singularity")
    ax[1].plot([v0, Cp[0]], [v0, Cp[1]], color=PALETTE["blue"], lw=1.4, label="event horizon")
    P_ = (Ep[0], A[1])
    ax[1].fill([A[0], Ep[0], P_[0]], [A[1], Ep[1], P_[1]], facecolor="none", edgecolor=PALETTE["grey"], hatch="//", lw=1.0, alpha=0.6)
    ax[1].plot([Ep[0], P_[0]], [Ep[1], P_[1]], color=PALETTE["grey"], ls="--", lw=1.0)
    ax[1].text(P_[0] + 0.55, A[1] + 0.9, "flat interior\n(schematic,\nnot to scale)", ha="center", fontsize=8, color=PALETTE["grey"])
    ax[1].text(2.15, -0.1, "I", fontsize=13), ax[1].text(1.0, 1.05, "II", fontsize=13)
    ax[1].text(-2.5, 0.9, "no region III or IV", fontsize=10, color=PALETTE["red"])
    ax[1].set_xlim(-PI - 0.3, PI + 0.3), ax[1].set_ylim(-PI / 2 - 0.45, PI / 2 + 0.45), ax[1].set_aspect("equal")
    ax[1].set_xticks([]), ax[1].set_yticks([]), ax[1].grid(False)
    ax[1].set_title("(b) Collapse (Vaidya null shell): only I and II survive")
    ax[1].legend(loc="lower left", fontsize=8)
    # (c) Ellis wormhole in (l, t): the future light cone of an event on the left side reaches the right side
    l = np.linspace(-8, 8, 200)
    l0 = -4.0
    ax[2].fill_between(l, np.abs(l - l0), 14, color=PALETTE["orange"], alpha=0.12, label="causal future of p")
    ax[2].plot(l, np.abs(l - l0), color=PALETTE["orange"], lw=1.6)
    tw = (l[l >= l0] - l0) / 0.6
    ax[2].plot(l[l >= l0], tw, color=PALETTE["blue"], lw=2, label="timelike traveller, v = 0.6")
    ax[2].plot([l0], [0], "ko")
    ax[2].axvline(0, color=PALETTE["grey"], ls="--", lw=1)
    ax[2].text(0.2, 12.3, "throat l = 0", color=PALETTE["grey"])
    ax[2].set_xlim(-8, 8), ax[2].set_ylim(0, 14)
    ax[2].set_xlabel("proper length l / b0"), ax[2].set_ylabel("t / b0")
    ax[2].set_title("(c) Ellis wormhole: the future of p crosses l = 0")
    ax[2].legend(loc="lower right", fontsize=8)
    save(fig, "fig10_causal_structure", FIG_DIR)

    bad, total = CA.count_signalling_pairs(150, seed=11)
    out["I_to_III_possible_pairs"] = {"bad": bad, "total": total}
    return out


if __name__ == "__main__":
    print(record("exp03_bridge_causality", run()))
