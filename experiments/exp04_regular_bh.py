"""Experiment 4: Hayward effective regular black hole versus Schwarzschild (Part V, Question H)."""
import matplotlib.pyplot as plt
import numpy as np

from beyond_spacetime import perturbations as P
from beyond_spacetime import regular as H
from beyond_spacetime.plotting import CYCLE, PALETTE, save

from common import FIG_DIR, record


def run():
    out = {}
    M = 1.0
    ls = [0.3, 0.5, 0.7, H.l_critical_over_M() * 0.999]
    r = np.linspace(0.0005, 6, 4000)

    # ---- Figure 11: metric function and curvature invariants ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
    ax[0].plot(r, 1 - 2 * M / r, "k--", lw=1.2, label="Schwarzschild")
    for l, c in zip(ls, CYCLE):
        ax[0].plot(r, H.f(r, M, l), color=c, label=f"l/M = {l:.3f}")
    ax[0].axhline(0, color="black", lw=0.7)
    ax[0].set_ylim(-3, 1.3), ax[0].set_xlim(0, 6)
    ax[0].set_xlabel("r / M"), ax[0].set_ylabel("f(r)")
    ax[0].set_title("Metric function: de Sitter core (f -> 1) and two horizons")
    ax[0].legend(fontsize=8)
    rr = np.logspace(-2, 1.2, 400)
    ax[1].loglog(rr, 48 * M**2 / rr**6, "k--", lw=1.2, label="Schwarzschild")
    for l, c in zip(ls, CYCLE):
        ax[1].loglog(rr, H.kretschmann(rr, M, l), color=c, label=f"l/M = {l:.3f}")
        ax[1].axhline(H.kretschmann_core_limit(l), color=c, ls=":", lw=0.6)
    ax[1].set_xlabel("r / M"), ax[1].set_ylabel("Kretschmann K M$^4$")
    ax[1].set_title("Curvature saturates at 24/l$^4$ (dotted)")
    ax[1].legend(fontsize=8)
    save(fig, "fig11_hayward_metric_curvature", FIG_DIR)
    out["core_curvature"] = {str(l): H.kretschmann_core_limit(l) for l in ls}

    # ---- Figure 12: effective stress-energy and energy conditions ----
    l = 0.5
    rr = np.linspace(0.01, 8, 3000)
    rho, pr, pt = H.stress_energy(rr, M, l)
    ec = H.energy_conditions(rr, M, l)
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
    ax[0].plot(rr, rho, color=PALETTE["blue"], label="$\\rho$")
    ax[0].plot(rr, pr, color=PALETTE["orange"], label="$p_r = -\\rho$")
    ax[0].plot(rr, pt, color=PALETTE["green"], label="$p_t$")
    ax[0].set_xlabel("r / M"), ax[0].set_ylabel("effective T$^\\mu{}_\\nu$ components (G = 1)")
    ax[0].set_title(f"Effective stress-energy, l/M = {l}")
    ax[0].legend()
    names = ["NEC", "WEC", "SEC", "DEC"]
    for i, n in enumerate(names):
        ok = ec[n]
        ax[1].fill_between(rr, i, i + 0.7, where=ok, color=PALETTE["green"], alpha=0.8, step="mid")
        ax[1].fill_between(rr, i, i + 0.7, where=~ok, color=PALETTE["red"], alpha=0.8, step="mid")
    ax[1].set_yticks([i + 0.35 for i in range(4)]), ax[1].set_yticklabels(names)
    ax[1].set_xlabel("r / M"), ax[1].set_title("Energy conditions (green: satisfied, red: violated)")
    ax[1].grid(False)
    save(fig, "fig12_hayward_energy_conditions", FIG_DIR)
    out["energy_conditions_violation_radii"] = {
        n: ([float(rr[~ec[n]].min()), float(rr[~ec[n]].max())] if (~ec[n]).any() else None) for n in names}

    # ---- Figure 13: horizons and temperature ----
    ratio = np.linspace(0.02, H.l_critical_over_M() * 0.9999, 300)
    rp = np.array([H.horizons(M, x)[-1] for x in ratio])
    rm = np.array([H.horizons(M, x)[0] for x in ratio])
    T = np.array([H.hawking_temperature(M, x) for x in ratio])
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax[0].plot(ratio, rp, color=PALETTE["blue"], label="outer $r_+$")
    ax[0].plot(ratio, rm, color=PALETTE["orange"], label="inner $r_-$")
    ax[0].axvline(H.l_critical_over_M(), color=PALETTE["red"], ls="--", lw=1)
    ax[0].text(H.l_critical_over_M() - 0.28, 0.2, "extremal\nl = 4M/(3$\\sqrt{3}$)", color=PALETTE["red"], fontsize=9)
    ax[0].set_xlabel("l / M"), ax[0].set_ylabel("horizon radius / M"), ax[0].legend()
    ax[0].set_title("Horizon structure of the effective model")
    ax[1].plot(ratio, T * 8 * np.pi * M, color=PALETTE["purple"])
    ax[1].set_xlabel("l / M"), ax[1].set_ylabel("T / T$_{Schw}$")
    ax[1].set_title("Temperature vanishes at extremality")
    save(fig, "fig13_hayward_horizons_temperature", FIG_DIR)
    out["l_critical_over_M"] = H.l_critical_over_M()

    # ---- Figure 14: radial geodesics and ringdown comparison ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.9))
    lh = 0.5
    from scipy.integrate import solve_ivp
    for E, c in zip((0.5, 0.8, 1.2), CYCLE):
        rp_ = H.horizons(M, lh)[-1]
        # integrate r'' = -f'(r)/2 inward from the outer horizon with rdot = -sqrt(E^2 - f) (= -|E|)
        rhs = lambda tau, y: [y[1], -0.5 * float(H.fprime(max(y[0], 1e-8), M, lh))]
        ev = lambda tau, y: y[1]
        ev.direction = 1
        ev.terminal = False
        r_start = rp_ + 0.0
        centre = lambda tau, y: y[0] - 1e-4
        centre.terminal = True
        centre.direction = -1
        sol = solve_ivp(rhs, [0, 12], [r_start, -abs(E)], rtol=1e-10, atol=1e-12, dense_output=True,
                        events=[centre])
        tt = np.linspace(0, sol.t[-1], 600)
        ax[0].plot(tt, sol.sol(tt)[0], color=c, label=f"E = {E}")
        if len(sol.t_events[0]):
            ax[0].plot([sol.t[-1]], [1e-4], "o", color=c, ms=6)
            ax[0].annotate("reaches regular centre r = 0", (sol.t[-1], 0.0), xytext=(sol.t[-1] + 1.0, 1.1),
                           arrowprops=dict(arrowstyle="->"), fontsize=8)
    for rh in H.horizons(M, lh):
        ax[0].axhline(rh, color=PALETTE["grey"], ls=":", lw=1)
    ax[0].set_xlabel("proper time from $r_+$ / M"), ax[0].set_ylabel("r / M")
    ax[0].set_title("E < 1: bounce inside the core; E >= 1: reaches the centre", fontsize=10)
    ax[0].legend()
    # ringdown comparison for the l = 2 scalar mode
    results = {}
    for lval, c in [(1e-4, "k"), (0.3, PALETTE["blue"]), (0.5, PALETTE["orange"]), (0.7, PALETTE["green"])]:
        rpl = H.horizons(M, lval)[-1]
        V = P.potential_function_rs(lambda x, lv=lval: H.f(x, M, lv), lambda x, lv=lval: H.fprime(x, M, lv), rpl, l=2)
        t, psi = P.evolve(V, rs_min=-300, rs_max=800, dx=0.05, t_end=220, r0s=20, sigma=3, r_obs_s=10)
        ax[1].semilogy(t, np.abs(psi) + 1e-20, color=c, label=f"l/M = {lval:g}")
        w = P.ringdown_frequency(t, psi, 60, 110, order=8)
        results[str(lval)] = {"re": float(w.real), "im": float(w.imag)}
    ax[1].set_ylim(1e-9, 1), ax[1].set_xlim(0, 220)
    ax[1].set_xlabel("t / M"), ax[1].set_ylabel("|$\\psi$| at r* = 10 M")
    ax[1].set_title("Scalar l = 2 ringdown outside the outer horizon (decays)")
    ax[1].legend(fontsize=8)
    save(fig, "fig14_hayward_geodesics_ringdown", FIG_DIR)
    out["ringdown_scalar_l2"] = results
    out["tau_outer_horizon_to_turning_point_E0.5"] = H.tau_outer_horizon_to_turning_point(0.5, M, 0.5)[0]
    out["tau_to_center_E1.2"] = H.tau_to_center(1.2, M, 0.5)
    return out


if __name__ == "__main__":
    print(record("exp04_regular_bh", run()))
