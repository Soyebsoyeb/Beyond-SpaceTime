"""Experiment 5: AdS3/BTZ holography, RT, ER-bridge growth, thermofield-double toy (Part VI)."""
import matplotlib.pyplot as plt
import numpy as np

from beyond_spacetime import holography as Hh
from beyond_spacetime import tfd
from beyond_spacetime.plotting import CYCLE, PALETTE, save

from common import FIG_DIR, record


def run():
    out = {}
    rc = 1e5
    # ---- Figure 15: RT entropy: numerics vs closed form vs CFT; homology phase transition ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.9))
    dphis = np.linspace(0.2, np.pi, 22)
    maxdiff = 0.0
    for rp, c in zip((1.0, 1.5, 3.0), CYCLE):
        Snum, Sref = [], []
        for d in dphis:
            J = Hh.geodesic_J_from_angle(d, rp)
            Snum.append(Hh.length_numeric(J, rp, rc) / 4.0)
            Sref.append(Hh.entanglement_entropy_CFT(d, rp, rc))
            maxdiff = max(maxdiff, abs(Snum[-1] - Sref[-1]))
        dd = np.linspace(0.2, np.pi, 300)
        ax[0].plot(dd, Hh.entanglement_entropy_CFT(dd, rp, rc), color=c, label=f"$r_+$ = {rp} (T = {rp/(2*np.pi):.3f})")
        ax[0].plot(dphis, Snum, "o", color=c, ms=3.5)
    ax[0].set_xlabel("interval angle $\\Delta\\phi$ ($\\leq \\pi$)"), ax[0].set_ylabel(f"S$_A$, cutoff r$_c$ = {rc:.0e}, G$_N$ = 1")
    ax[0].set_title("Connected RT geodesic: numerical (dots) vs CFT formula (lines)", fontsize=10)
    ax[0].legend(fontsize=8)
    out["rt_vs_cft_max_abs_diff"] = maxdiff
    rp = 1.5
    ds = Hh.transition_angle(rp)
    dd = np.linspace(0.3, 2 * np.pi - 0.3, 600)
    conn = Hh.geodesic_length_closed(dd, rp, rc) / 4.0
    disc = (Hh.geodesic_length_closed(2 * np.pi - dd, rp, rc) + 2 * np.pi * rp) / 4.0
    ax[1].plot(dd, conn, color=PALETTE["blue"], ls="--", lw=1.2, label="connected geodesic only")
    ax[1].plot(dd, disc, color=PALETTE["orange"], ls="--", lw=1.2, label="complement geodesic + horizon")
    ax[1].plot(dd, Hh.rt_entropy(dd, rp, rc), "k", lw=2.2, label="RT: minimum of the two")
    ax[1].axvline(ds, color=PALETTE["red"], ls=":", lw=1)
    ax[1].text(ds + 0.08, conn.min() + 0.2, f"$\\Delta\\phi^*$ = {ds:.3f}", color=PALETTE["red"], fontsize=9)
    ax[1].set_xlabel("interval angle $\\Delta\\phi$ in (0, 2$\\pi$)"), ax[1].set_ylabel("S$_A$")
    ax[1].set_title("Homology constraint: the minimal surface changes topology", fontsize=10)
    ax[1].legend(fontsize=8)
    out["transition_angle_rplus_1.5"] = ds
    Js = np.linspace(1.01, 8, 200) * 1.5
    ang = np.array([Hh.angle_from_J_numeric(J, 1.5) for J in Js])
    ax[2].plot(Js, ang, color=PALETTE["purple"], label="numerical integral")
    ax[2].plot(Js, 2 * np.arctanh(1.5 / Js) / 1.5, "k--", lw=1.0, label="closed form  (2/r$_+$) arctanh(r$_+$/J)")
    ax[2].set_xlabel("turning point J = r*"), ax[2].set_ylabel("$\\Delta\\phi$")
    ax[2].set_title("Geodesic subtended angle (r$_+$ = 1.5)", fontsize=10)
    ax[2].legend(fontsize=8)
    save(fig, "fig15_rt_btz", FIG_DIR)
    out["angle_closed_form_max_diff"] = float(np.max(np.abs(ang - 2 * np.arctanh(1.5 / Js) / 1.5)))
    out["araki_lieb_saturation_error"] = float(abs(
        (Hh.rt_entropy(0.5 * (ds + 2 * np.pi), rp, rc) - Hh.rt_entropy(2 * np.pi - 0.5 * (ds + 2 * np.pi), rp, rc))
        - Hh.bh_entropy(rp)))

    # ---- Figure 16: ER bridge length growth ----
    rp = 1.0
    t = np.linspace(0, 8, 200)
    L = np.array([Hh.two_sided_length(tt, rp, rc) for tt in t])
    Lc = Hh.two_sided_length_closed(t, rp, rc)
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax[0].plot(t, L - L[0], color=PALETTE["blue"], lw=3, alpha=0.6, label="exact embedding distance")
    ax[0].plot(t, Lc - Lc[0], "k--", lw=1.2, label="2 ln cosh($r_+$ t)")
    ax[0].plot(t, 2 * rp * t, color=PALETTE["red"], ls=":", label="late-time slope 2 $r_+$")
    ax[0].set_xlabel("boundary time t (t$_L$ = t$_R$ = t)"), ax[0].set_ylabel("L(t) - L(0)")
    ax[0].set_title("Einstein-Rosen bridge keeps growing after thermalisation")
    ax[0].legend(fontsize=8)
    ax[1].semilogy(t, np.abs(L - Lc) + 1e-16, color=PALETTE["purple"])
    ax[1].set_xlabel("t"), ax[1].set_ylabel("|exact - closed form|")
    ax[1].set_title("Agreement at the level of the cutoff correction")
    save(fig, "fig16_er_bridge_growth", FIG_DIR)
    out["er_bridge_max_abs_diff"] = float(np.max(np.abs(L - Lc)))
    rps = np.linspace(0.3, 4, 20)
    out["bh_vs_cardy_entropy_max_diff"] = float(np.max(np.abs(Hh.bh_entropy(rps) - Hh.cardy_entropy(rps))))

    # ---- Figure 17: TFD toy model ----
    N = 5
    H, op = tfd.ising_hamiltonian(N, h=1.0)
    d = 2**N
    betas = np.linspace(0.05, 4.0, 25)
    S_ee, S_th = [], []
    for b in betas:
        psi, E, V, Z = tfd.tfd_state(H, b)
        S_ee.append(tfd.entropy_of_reduced_L(psi, d)[0])
        S_th.append(tfd.thermal_entropy(E, b))
    gs = np.linspace(0, 0.6, 31)
    psi, *_ = tfd.tfd_state(H, 1.0)
    sig = [tfd.signal_with_coupling(psi, N, g, 0.7, op) for g in gs]
    rng = np.random.default_rng(0)
    dev = [tfd.no_signalling_deviation(tfd.tfd_state(H, b)[0], d, rng) for b in betas]
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8))
    ax[0].plot(betas, S_th, color=PALETTE["blue"], lw=3, alpha=0.6, label="thermal entropy of H")
    ax[0].plot(betas, S_ee, "k.", label="S(L) of TFD (QuTiP partial trace)")
    ax[0].axhline(N * np.log(2), color=PALETTE["red"], ls=":", label="N ln 2")
    ax[0].set_xlabel("$\\beta$"), ax[0].set_ylabel("entropy"), ax[0].legend(fontsize=8)
    ax[0].set_title("TFD: entanglement entropy = thermal entropy")
    ax[1].semilogy(betas, np.array(dev) + 1e-18, color=PALETTE["green"])
    ax[1].set_xlabel("$\\beta$"), ax[1].set_ylabel("max |$\\delta\\rho_R$| after random unitary on L")
    ax[1].set_title("No signalling: $\\rho_R$ is unchanged to round-off")
    ax[2].plot(gs, sig, color=PALETTE["orange"])
    ax[2].set_xlabel("L-R coupling strength g"), ax[2].set_ylabel("change in receiver $\\langle X \\rangle$")
    ax[2].set_title("Signal appears only with an explicit coupling")
    save(fig, "fig17_tfd_toy", FIG_DIR)
    out["tfd_entropy_max_diff"] = float(np.max(np.abs(np.array(S_ee) - np.array(S_th))))
    out["tfd_no_signalling_max_deviation"] = float(max(dev))
    out["tfd_signal_at_g0"] = float(sig[0])
    out["tfd_signal_at_g0.3"] = float(sig[15])
    return out


if __name__ == "__main__":
    print(record("exp05_holography", run()))
