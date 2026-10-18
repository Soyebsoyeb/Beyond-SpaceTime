"""Experiment 9: Schwarzschild QNMs and perturbative stability (Part VIII, item 17); multiverse comparison (Part VII)."""
import matplotlib.pyplot as plt
import numpy as np

from beyond_spacetime import multiverse as MV
from beyond_spacetime import perturbations as P
from beyond_spacetime.plotting import CYCLE, PALETTE, save

from common import FIG_DIR, record


def run():
    out = {}
    ref = {"s2_l2_n0": P.leaver_qnm(2, 2, 0), "s2_l2_n1": P.leaver_qnm(2, 2, 1), "s0_l2_n0": P.leaver_qnm(2, 0, 0),
           "s2_l3_n0": P.leaver_qnm(3, 2, 0)}
    out["leaver"] = {k: {"re": v.real, "im": v.imag} for k, v in ref.items()}
    published = {"s2_l2_n0": (0.37367, -0.08896), "s2_l2_n1": (0.34671, -0.27391), "s0_l2_n0": (0.48364, -0.09677)}
    out["published_comparison_max_abs_diff"] = max(
        max(abs(ref[k].real - p[0]), abs(ref[k].imag - p[1])) for k, p in published.items())

    # ---- Figure 26: time-domain ringdown vs Leaver, for scalar and gravitational (RW) potentials ----
    f = lambda r: 1 - 2 / r
    fp = lambda r: 2 / r**2
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.9))
    V = P.potential_function_rs(f, fp, 2.0, l=2)
    t, psi = P.evolve(V, rs_min=-300, rs_max=800, dx=0.05, t_end=250, r0s=20, sigma=3, r_obs_s=10)
    w_td = P.ringdown_frequency(t, psi, 60, 110, order=8)
    w_ex = ref["s0_l2_n0"]
    ax[0].semilogy(t, np.abs(psi) + 1e-20, color=PALETTE["blue"])
    tt = np.linspace(60, 110, 100)
    amp = np.abs(psi[np.argmin(np.abs(t - 60))])
    ax[0].semilogy(tt, amp * np.exp(w_ex.imag * (tt - 60)) * 1.0, "k--", lw=1.2,
                   label=f"exp($\\omega_I$ t), $\\omega_I$ = {w_ex.imag:.5f} (Leaver)")
    ax[0].axvspan(60, 110, color=PALETTE["orange"], alpha=0.12, label="Prony fit window")
    ax[0].set_ylim(1e-10, 1), ax[0].set_xlabel("t / M"), ax[0].set_ylabel("|$\\psi$| at r* = 10 M")
    ax[0].set_title("Scalar l = 2 ringdown and late-time power-law tail"), ax[0].legend(fontsize=8)
    # gravitational (Regge-Wheeler) potential
    from scipy.interpolate import interp1d
    from beyond_spacetime.perturbations import tortoise_grid, rw_potential_schwarzschild
    r, rs = tortoise_grid(f, 2.0)
    Vi = interp1d(rs, rw_potential_schwarzschild(r, 2, 2), bounds_error=False, fill_value=(0, 0))
    t2, psi2 = P.evolve(lambda x: Vi(x), rs_min=-300, rs_max=800, dx=0.05, t_end=250, r0s=20, sigma=3, r_obs_s=10)
    w_td2 = P.ringdown_frequency(t2, psi2, 60, 110, order=8)
    ax[1].semilogy(t2, np.abs(psi2) + 1e-20, color=PALETTE["green"])
    ax[1].axvspan(60, 110, color=PALETTE["orange"], alpha=0.12)
    ax[1].set_ylim(1e-10, 1), ax[1].set_xlabel("t / M")
    ax[1].set_title("Regge-Wheeler (gravitational, l = 2) ringdown")
    labels = ["scalar l=2", "gravitational l=2"]
    tds = [w_td, w_td2]
    exs = [ref["s0_l2_n0"], ref["s2_l2_n0"]]
    ax[2].scatter([e.real for e in exs], [e.imag for e in exs], marker="x", s=80, color="black", label="Leaver")
    ax[2].scatter([w.real for w in tds], [w.imag for w in tds], marker="o", s=40, facecolors="none",
                  edgecolors=PALETTE["red"], label="time domain + Prony")
    for e, lab in zip(exs, labels):
        ax[2].annotate(lab, (e.real, e.imag), xytext=(e.real + 0.01, e.imag - 0.01), fontsize=8)
    ax[2].set_xlabel("M Re $\\omega$"), ax[2].set_ylabel("M Im $\\omega$")
    ax[2].set_title("Quasinormal frequencies: independent methods agree"), ax[2].legend(fontsize=8)
    save(fig, "fig26_qnm_ringdown", FIG_DIR)
    out["time_domain"] = {"scalar_l2": {"re": w_td.real, "im": w_td.imag, "ref_re": w_ex.real, "ref_im": w_ex.imag},
                          "grav_l2": {"re": w_td2.real, "im": w_td2.imag, "ref_re": ref["s2_l2_n0"].real,
                                      "ref_im": ref["s2_l2_n0"].imag}}

    # ---- Figure 27: classification matrix of "multiple worlds" concepts ----
    fig, ax = plt.subplots(figsize=(11, 3.6))
    cats = [MV.EQUATIONS, MV.INTERPRETIVE, MV.EXTRA]
    colors = {MV.EQUATIONS: PALETTE["green"], MV.INTERPRETIVE: PALETTE["purple"], MV.EXTRA: PALETTE["orange"]}
    for i, c in enumerate(MV.CONCEPTS):
        ax.barh(i, 1.0, color=colors[c["status"]], alpha=0.85)
        ax.text(0.02, i, c["name"], va="center", color="white", fontsize=10, fontweight="bold")
        ax.text(1.02, i, c["status"], va="center", fontsize=9)
    ax.set_yticks([]), ax.set_xticks([]), ax.set_xlim(0, 1.9), ax.invert_yaxis(), ax.grid(False)
    ax.set_title("Status of the multiplicity in each proposal (no ER-bridge/many-worlds link is established)")
    save(fig, "fig27_multiverse_classification", FIG_DIR)
    out["multiverse_status"] = {c["key"]: c["status"] for c in MV.CONCEPTS}
    return out


if __name__ == "__main__":
    print(record("exp09_qnm_multiverse", run()))
