"""Experiment 8: Morris-Thorne wormholes, closed timelike curves, Goedel (Part VIII, Question G)."""
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

from beyond_spacetime import chronology as CH
from beyond_spacetime import symbolic as S
from beyond_spacetime import wormholes as W
from beyond_spacetime.plotting import CYCLE, PALETTE, save

from common import FIG_DIR, record


def run():
    out = {}
    b0 = 1.0
    # ---- Figure 23: Ellis wormhole geometry and exotic matter ----
    l = np.linspace(-4, 4, 400)
    r = W.ellis_areal_radius(l, b0)
    z = W.ellis_embedding_z(l, b0)
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.9))
    ax[0].plot(r, z, color=PALETTE["blue"]), ax[0].plot(r, -z * 0 + z, alpha=0)
    ax[0].plot(-r, z, color=PALETTE["blue"])
    ax[0].set_xlabel("$\\pm$ r / b0"), ax[0].set_ylabel("embedding height z / b0")
    ax[0].set_title("Ellis wormhole embedding: throat at r = b0")
    ax[0].set_aspect("equal")
    rho_sym, tau_sym, p_sym = W.mt_matter(sp.Symbol("b0", positive=True) ** 2 / S.r, sp.Integer(0))
    rr = np.linspace(b0, 5, 300)
    f_rho = sp.lambdify((S.r, sp.Symbol("b0", positive=True)), rho_sym, "numpy")
    f_tau = sp.lambdify((S.r, sp.Symbol("b0", positive=True)), tau_sym, "numpy")
    f_p = sp.lambdify((S.r, sp.Symbol("b0", positive=True)), p_sym, "numpy")
    ax[1].plot(rr, f_rho(rr, b0), color=PALETTE["blue"], lw=5, alpha=0.45, label="$\\rho$ (equals $p_r$ here)")
    ax[1].plot(rr, -f_tau(rr, b0), color=PALETTE["orange"], label="$p_r = -\\tau$")
    ax[1].plot(rr, f_p(rr, b0) * np.ones_like(rr), color=PALETTE["green"], label="$p_t$")
    ax[1].plot(rr, f_rho(rr, b0) - f_tau(rr, b0), "k--", label="$\\rho + p_r$ (radial NEC)")
    ax[1].axhline(0, color="black", lw=0.7)
    ax[1].set_xlabel("r / b0"), ax[1].set_ylabel("in units of 1/b0$^2$")
    ax[1].set_title("Matter needed: radial NEC violated everywhere")
    ax[1].legend(fontsize=8)
    vs = (0.1, 0.5, 0.9, 0.99)
    lt = np.linspace(-3, 3, 300)
    for v, c in zip(vs, CYCLE):
        ax[2].plot(lt, W.ellis_tidal_lateral(lt, v, b0), color=c, label=f"v = {v}")
    ax[2].set_xlabel("proper length l / b0"), ax[2].set_ylabel("lateral tidal component  (1/b0$^2$)")
    ax[2].set_title("Tides on a radial traveller (static observers feel none)")
    ax[2].legend(fontsize=8)
    save(fig, "fig23_ellis_wormhole", FIG_DIR)
    out["throat_radial_NEC_over_inv_b0sq"] = W.nec_radial_at_throat(b0)
    out["throat_NEC_formula_minus_1_over_4pi"] = -1.0 / (4 * np.pi)

    # ---- Figure 24: chronology horizon and CTC region ----
    betas = np.linspace(0.05, 0.99, 200)
    Ds = np.linspace(0.1, 6.0, 200)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.0))
    for D, c in zip((0.5, 1.0, 2.0), CYCLE):
        Th = CH.chronology_horizon_time(D, betas)
        ax[0].semilogy(betas, Th, color=c, label=f"separation D = {D}")
    ax[0].set_xlabel("mouth speed $\\beta$ (c = 1)"), ax[0].set_ylabel("trip duration before first CTC")
    ax[0].set_title("Chronology horizon: earliest closed timelike curve"), ax[0].legend(fontsize=8)
    Tt = np.linspace(0, 12, 300)
    for beta, c in zip((0.3, 0.6, 0.9), CYCLE):
        ax[1].plot(Tt, CH.clock_lag(Tt, beta), color=c, label=f"$\\beta$ = {beta}")
    ax[1].axhline(1.0, color=PALETTE["red"], ls="--", lw=1)
    ax[1].text(0.2, 1.1, "D = 1: CTC when lag > D", color=PALETTE["red"], fontsize=9)
    ax[1].set_xlabel("lab time of trip"), ax[1].set_ylabel("clock lag $\\Delta$ = T (1 - sqrt(1 - $\\beta^2$))")
    ax[1].set_title("Time offset accumulated between the mouths"), ax[1].legend(fontsize=8)
    m = CH.ctc_map(Ds, betas, 5.0)
    ax[2].pcolormesh(Ds, betas, m.astype(float), cmap="RdYlGn_r", shading="auto", alpha=0.85)
    ax[2].set_xlabel("mouth separation D"), ax[2].set_ylabel("speed $\\beta$")
    ax[2].set_title("CTC exist (red) iff $\\Delta$ > D (trip time 5)")
    ax[2].grid(False)
    save(fig, "fig24_chronology", FIG_DIR)
    out["chronology_horizon_D1_beta0.8"] = CH.chronology_horizon_time(1.0, 0.8)

    # ---- Figure 25: Goedel closed timelike circles ----
    rg = np.linspace(0, 2.0, 400)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax[0].plot(rg, CH.godel_gphiphi(rg), color=PALETTE["blue"])
    ax[0].axhline(0, color="black", lw=0.8)
    ax[0].axvline(CH.R_CRIT_GODEL, color=PALETTE["red"], ls="--", lw=1)
    ax[0].text(CH.R_CRIT_GODEL + 0.03, 3, "$r_c$ = ln(1 + $\\sqrt{2}$) = %.4f" % CH.R_CRIT_GODEL, color=PALETTE["red"])
    ax[0].set_xlabel("r"), ax[0].set_ylabel("g$_{\\phi\\phi}$ / a$^2$")
    ax[0].set_title("Goedel: g$_{\\phi\\phi}$ < 0 beyond r$_c$ (closed timelike circles)")
    th = np.linspace(0, 2 * np.pi, 400)
    ax[1].fill(CH.R_CRIT_GODEL * np.cos(th), CH.R_CRIT_GODEL * np.sin(th), color=PALETTE["green"], alpha=0.25)
    for rr_ in (0.5, CH.R_CRIT_GODEL, 1.4):
        col = PALETTE["red"] if rr_ > CH.R_CRIT_GODEL * 0.999 else PALETTE["grey"]
        ax[1].plot(rr_ * np.cos(th), rr_ * np.sin(th), color=col, lw=1.4)
    ax[1].set_aspect("equal"), ax[1].set_title("Circles outside r$_c$ are closed timelike curves")
    ax[1].text(0, 0, "causal", ha="center", color=PALETTE["green"])
    save(fig, "fig25_godel", FIG_DIR)
    out["godel_r_crit"] = CH.R_CRIT_GODEL
    return out


if __name__ == "__main__":
    print(record("exp08_wormholes_chronology", run()))
