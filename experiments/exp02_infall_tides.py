"""Experiment 2: radial infall, proper times, tidal dynamics (Part III, Questions A)."""
import matplotlib.pyplot as plt
import numpy as np

from beyond_spacetime import geodesics as G
from beyond_spacetime import tides as T
from beyond_spacetime.conventions import C_SI, G_EARTH, geometric_mass_m, M_SUN_KG
from beyond_spacetime.plotting import CYCLE, PALETTE, save

from common import FIG_DIR, record


def run():
    out = {}
    # ---- Figure 3: infall trajectories for several energies ----
    r0 = 12.0
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9))
    rows = []
    for E, c in zip([G.energy_from_rest(r0), 1.0, 1.5], CYCLE):
        tr = G.integrate_infall(r0, E, r_end=1e-3, n=3000)
        lab = f"E = {E:.3f}"
        ax[0].plot(tr["tau"], tr["r"], color=c, label=lab)
        out_side = tr["r"] > 2.0 * 1.0002
        ax[1].plot(G.schwarzschild_time(tr["r"][out_side], tr["v"][out_side]), tr["r"][out_side], color=c, label=lab)
        ax[2].plot(tr["tau"], tr["v"], color=c, label=lab)
        rows.append({"E": E, "tau_horizon": tr["tau_horizon"], "tau_end": tr["tau_end"],
                     "tau_inside": tr["tau_end"] - tr["tau_horizon"]})
    ax[0].axhline(2, color=PALETTE["red"], ls="--", lw=1)
    ax[0].set_xlabel("proper time $\\tau$ / M"), ax[0].set_ylabel("r / M"), ax[0].set_title("Areal radius vs proper time")
    ax[0].legend()
    ax[1].axhline(2, color=PALETTE["red"], ls="--", lw=1)
    ax[1].set_xlabel("Schwarzschild t / M"), ax[1].set_ylabel("r / M")
    ax[1].set_title("Coordinate time diverges at the horizon")
    ax[2].set_xlabel("proper time $\\tau$ / M"), ax[2].set_ylabel("ingoing EF v / M")
    ax[2].set_title("EF v is finite across the horizon")
    save(fig, "fig03_infall_trajectories", FIG_DIR)
    out["infall"] = rows

    # ---- Figure 4: remaining proper time vs energy ----
    Es = np.linspace(0, 4, 200)
    taus = np.array([G.tau_inside(E) for E in Es])
    fig, ax = plt.subplots(figsize=(6, 3.8))
    ax.plot(Es, taus / np.pi, color=PALETTE["blue"])
    ax.axhline(1, color=PALETTE["red"], ls="--", lw=1)
    ax.text(2.0, 1.01, "maximum $\\pi M$ (E = 0)", color=PALETTE["red"])
    ax.scatter([1.0], [4 / (3 * np.pi)], color=PALETTE["orange"], zorder=3)
    ax.annotate("E = 1: 4M/3", (1.0, 4 / (3 * np.pi)), xytext=(1.5, 0.6), arrowprops=dict(arrowstyle="->"))
    ax.set_xlabel("conserved energy E"), ax.set_ylabel("proper time horizon to singularity / ($\\pi$ M)")
    ax.set_title("Time left after crossing the horizon")
    save(fig, "fig04_proper_time_inside", FIG_DIR)
    out["tau_inside_E0_over_pi"] = G.tau_inside(0.0) / np.pi
    out["tau_inside_E1"] = G.tau_inside(1.0)

    # ---- Figure 5: curvature and tidal eigenvalues along the infall ----
    tr = G.integrate_infall(r0, 1.0, r_end=0.05, n=3000)
    r = tr["r"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    ax[0].semilogy(tr["tau"], G.kretschmann(r), color=PALETTE["orange"], label="K M$^4$")
    ax[0].semilogy(tr["tau"], 2 / r**3, color=PALETTE["blue"], label="radial tide 2M/r$^3$")
    ax[0].semilogy(tr["tau"], 1 / r**3, color=PALETTE["green"], label="transverse tide M/r$^3$")
    ax[0].axvline(tr["tau_horizon"], color=PALETTE["red"], ls="--", lw=1)
    ax[0].set_xlabel("proper time $\\tau$ / M"), ax[0].set_title("Invariant growth along an E = 1 infall")
    ax[0].legend()
    dev = T.integrate_deviation(r0, 1.0, xi0=(1.0, 1.0), r_end=0.05)
    ax[1].semilogy(dev["tau"], np.abs(dev["xi_r"]), color=PALETTE["blue"], label="radial |$\\xi_1$|")
    ax[1].semilogy(dev["tau"], np.abs(dev["xi_perp"]), color=PALETTE["green"], label="transverse |$\\xi_\\perp$|")
    ax[1].set_xlabel("proper time $\\tau$ / M"), ax[1].set_ylabel("separation (initial value 1)")
    ax[1].set_title("Geodesic deviation: stretching and compression")
    ax[1].legend()
    save(fig, "fig05_tides_and_deviation", FIG_DIR)
    out["deviation_final"] = {"xi_r": float(dev["xi_r"][-1]), "xi_perp": float(dev["xi_perp"][-1]),
                              "tau_end": float(dev["tau_end"])}

    # ---- Figure 6: validation of the deviation equation against neighbouring geodesics ----
    r0v, Ev, dE = 10.0, 0.9, 1e-6
    tau, sep = T.neighbouring_geodesic_separation(r0v, Ev, dE, n=400, r_end=0.2)
    xidot0 = sep[1] / tau[1]
    devv = T.integrate_deviation(r0v, Ev, xi0=(0.0, 0.0), xidot0=(xidot0, 0.0), r_end=0.2, n=400)
    xi_ref = np.interp(tau, devv["tau"], np.abs(devv["xi_r"]))
    m = tau < devv["tau_end"] * 0.95
    rel = np.abs(xi_ref[m][1:] - sep[m][1:]) / xi_ref[m][1:]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.7))
    ax[0].plot(tau, sep / sep[-1], color=PALETTE["orange"], lw=3, alpha=0.6, label="exact neighbouring geodesic")
    ax[0].plot(tau[m], xi_ref[m] / sep[-1], "k--", lw=1.2, label="geodesic deviation ODE")
    ax[0].set_xlabel("proper time $\\tau$ / M"), ax[0].set_ylabel("proper separation (arb.)")
    ax[0].legend(), ax[0].set_title("Independent validation of $\\xi_1$ evolution")
    ax[1].semilogy(tau[m][1:], rel + 1e-16, color=PALETTE["purple"])
    ax[1].set_xlabel("proper time $\\tau$ / M"), ax[1].set_ylabel("relative difference")
    ax[1].set_title("Residual (finite-difference limited)")
    save(fig, "fig06_deviation_validation", FIG_DIR)
    out["deviation_validation_max_rel_err"] = float(rel.max())

    # ---- Figure 7: spaghettification versus mass ----
    masses = np.logspace(0, 10, 200)
    ratio = np.array([T.spaghettification_radius(m)[1] for m in masses])
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.loglog(masses, ratio, color=PALETTE["blue"])
    ax.axhline(1.0, color=PALETTE["red"], ls="--", lw=1)
    ax.set_xlabel("black-hole mass / M$_\\odot$")
    ax.set_ylabel("radius where tide = 10 g across 2 m / $r_s$")
    ax.set_title("Tides become fatal outside the horizon only for small masses")
    msun_cross = np.interp(1.0, ratio[::-1], masses[::-1])
    ax.axvline(msun_cross, color=PALETTE["grey"], ls=":")
    ax.text(msun_cross * 1.3, 30, f"crossover ~ {msun_cross:.1e} M$_\\odot$", fontsize=9)
    save(fig, "fig07_spaghettification", FIG_DIR)
    out["spaghettification_crossover_solar_masses"] = float(msun_cross)
    out["tidal_acceleration_2m_at_horizon_10Msun_m_s2"] = None
    from beyond_spacetime.conventions import tidal_accel_si
    Mg = geometric_mass_m(10 * M_SUN_KG)
    out["tidal_acceleration_2m_at_horizon_10Msun_m_s2"] = tidal_accel_si(Mg, 2 * Mg, 2.0)
    Mg = geometric_mass_m(4e6 * M_SUN_KG)
    out["tidal_acceleration_2m_at_horizon_4e6Msun_m_s2"] = tidal_accel_si(Mg, 2 * Mg, 2.0)
    return out


if __name__ == "__main__":
    print(record("exp02_infall_tides", run()))
