"""Experiment 6: quantum branching, decoherence and coherence measures (Part VII)."""
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

from beyond_spacetime import quantum as Q
from beyond_spacetime.plotting import CYCLE, PALETTE, save

from common import FIG_DIR, record


def run():
    out = {}
    rng = np.random.default_rng(7)
    a, b = np.sqrt(0.3), np.sqrt(0.7)
    ts = np.linspace(0, 3.0, 61)

    # ---- Figure 18: decoherence vs environment size (exact dephasing model) ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.9))
    maxerr = 0.0
    for N, c in zip((1, 2, 4, 8), CYCLE):
        g = rng.uniform(0.5, 1.5, N)
        H = Q.dephasing_hamiltonian(g)
        states = Q.evolve_states(H, Q.initial_state(a, b, N), ts)
        coh = np.array([abs(Q.reduced_system(s).full()[0, 1]) for s in states])
        ana = a * b * np.abs([Q.decoherence_factor(g, t) for t in ts])
        maxerr = max(maxerr, float(np.max(np.abs(coh - ana))))
        ax[0].plot(ts, coh / (a * b), color=c, label=f"N = {N}")
        ax[0].plot(ts, ana / (a * b), "k:", lw=0.8)
        ent = [float(qt.entropy_vn(Q.reduced_system(s))) for s in states]
        ax[1].plot(ts, ent, color=c, label=f"N = {N}")
    ax[0].set_xlabel("time"), ax[0].set_ylabel("|$\\rho_{01}$| / |$\\rho_{01}(0)$|")
    ax[0].set_title("Coherence decay (dotted: analytic prod cos(2 g t))")
    ax[0].legend(loc="upper right", fontsize=8, ncol=2)
    H_b = -(0.3 * np.log(0.3) + 0.7 * np.log(0.7))
    ax[1].axhline(H_b, color=PALETTE["red"], ls="--", lw=1)
    ax[1].text(0.1, H_b + 0.01, "Shannon entropy of Born weights", color=PALETTE["red"], fontsize=8)
    ax[1].set_xlabel("time"), ax[1].set_ylabel("S(rho$_S$) / nats")
    ax[1].set_title("System entropy rises to the branch-weight entropy")
    ax[1].set_ylim(0, 0.75)
    ax[1].legend(loc="lower center", fontsize=8, ncol=4)
    out["dephasing_max_abs_error_vs_analytic"] = maxerr
    out["branch_weight_entropy"] = H_b

    # ---- recurrences for small environments (many samples, fixed g) ----
    tl = np.linspace(0, 40, 800)
    for N, c in zip((2, 6, 12), CYCLE):
        g = np.linspace(0.6, 1.4, N)
        D = np.array([abs(Q.decoherence_factor(g, t)) for t in tl])
        ax[2].plot(tl, D, color=c, label=f"N = {N}", lw=1.0)
    ax[2].set_xlabel("time"), ax[2].set_ylabel("|D(t)| = |$\\langle E_0|E_1\\rangle$|")
    ax[2].set_title("Finite environments recohere; recurrences vanish as N grows")
    ax[2].set_ylim(0, 1.25)
    ax[2].legend(loc="upper center", fontsize=8, ncol=3)
    save(fig, "fig18_decoherence", FIG_DIR)

    # ---- Figure 19: non-commuting interaction; positivity checks; coherence measures ----
    N = 6
    g = rng.uniform(0.4, 1.2, N)
    H2 = Q.noncommuting_hamiltonian(g, 0.15)
    psi0 = Q.initial_state(a, b, N)
    states = Q.evolve_states(H2, psi0, ts)
    rhos = [Q.reduced_system(s) for s in states]
    ent = [float(qt.entropy_vn(r)) for r in rhos]
    pur = [float((r * r).tr().real) for r in rhos]
    l1 = [Q.l1_coherence(r) for r in rhos]
    rec = [Q.relative_entropy_coherence(r) for r in rhos]
    checks = [Q.check_density_matrix(r) for r in rhos]
    mineig = np.array([c["min_eig"] for c in checks])
    herm = np.array([c["hermiticity_err"] for c in checks])
    trerr = np.array([c["trace_err"] for c in checks])
    dense_err = float(np.abs(states[-1].full()[:, 0] - Q.dense_evolve(H2, psi0, ts[-1])[:, 0]).max())
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8))
    ax[0].plot(ts, ent, color=PALETTE["blue"], label="S(rho$_S$)")
    ax[0].plot(ts, pur, color=PALETTE["orange"], label="purity Tr rho$^2$")
    ax[0].set_xlabel("time"), ax[0].set_title("Non-commuting coupling: mixedness of the system")
    ax[0].legend()
    ax[1].plot(ts, l1, color=PALETTE["green"], label="l1 coherence")
    ax[1].plot(ts, rec, color=PALETTE["purple"], label="relative-entropy coherence")
    ax[1].set_xlabel("time"), ax[1].set_title("Coherence measures in the energy-agnostic Z basis")
    ax[1].legend()
    ax[2].semilogy(ts, np.maximum(herm, 1e-18), color=PALETTE["red"], label="Hermiticity error")
    ax[2].semilogy(ts, np.maximum(trerr, 1e-18), color=PALETTE["blue"], label="|Tr rho - 1|")
    ax[2].plot(ts, np.maximum(mineig, 1e-18), color=PALETTE["green"], label="min eigenvalue (positive)")
    ax[2].set_xlabel("time"), ax[2].set_title("Density-matrix validity checks")
    ax[2].legend(fontsize=8)
    save(fig, "fig19_noncommuting_checks", FIG_DIR)
    out["noncommuting"] = {"min_eigenvalue_over_time": float(mineig.min()), "max_hermiticity_error": float(herm.max()),
                           "max_trace_error": float(trerr.max()), "qutip_vs_dense_expm_max_diff": dense_err,
                           "max_entropy": float(max(ent)), "ln2": float(np.log(2))}
    return out


if __name__ == "__main__":
    print(record("exp06_branching", run()))
