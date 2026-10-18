"""Experiment 7: time reversal, irreversibility and black/white holes (Part VIII)."""
import matplotlib.pyplot as plt
import numpy as np

from beyond_spacetime import coordinates as C
from beyond_spacetime import diagrams as D
from beyond_spacetime import geodesics as G
from beyond_spacetime import timereversal as TR
from beyond_spacetime.plotting import PALETTE, save

from common import FIG_DIR, record


def run():
    out = {}
    # ---- Figure 20: cat-map entropy production and Loschmidt reversal ----
    res = TR.cat_experiment(n=20000, t_max=50, seed=2)
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.9))
    ax[0].plot(res["t"], res["S"], "o-", ms=3, color=PALETTE["blue"], label="coarse-grained S(t)")
    ax[0].axhline(res["S_max"], color=PALETTE["red"], ls="--", lw=1, label="ln(cells)")
    ax[0].set_xlabel("map iterations t"), ax[0].set_ylabel("entropy / nats")
    ax[0].set_title("Entropy grows from a special initial macrostate"), ax[0].legend()
    ax[1].semilogy(res["t"], np.maximum(res["err_float"], 1e-17), "o-", ms=3, color=PALETTE["orange"], label="float64 reversal")
    ax[1].semilogy(res["t"], np.maximum(res["err_int"], 1e-17), "s-", ms=3, color=PALETTE["green"], label="exact integer lattice reversal")
    ref = 1e-16 * np.exp(TR.LYAPUNOV_CAT * res["t"])
    ax[1].semilogy(res["t"], np.minimum(ref, 1), "k:", label="1e-16 exp($\\lambda$ t)")
    ax[1].set_xlabel("iterations before reversal"), ax[1].set_ylabel("max position error after round trip")
    ax[1].set_title("Reversibility is exact in principle, fragile in practice"), ax[1].legend(fontsize=8)
    # reversal at t_r: entropy decreases back
    t_r = 20
    rng = np.random.default_rng(2)
    L = 2**40
    x0, y0 = rng.integers(0, L // 16, 20000, dtype=np.int64), rng.integers(0, L // 16, 20000, dtype=np.int64)
    S_path = []
    x, y = x0.copy(), y0.copy()
    for k in range(2 * t_r + 1):
        S_path.append(TR.coarse_grained_entropy(x / L, y / L, 8))
        x, y = TR.cat_map_int(x, y, 1, L, inverse=(k >= t_r))
    ax[2].plot(range(2 * t_r + 1), S_path, "o-", ms=3, color=PALETTE["purple"])
    ax[2].axvline(t_r, color=PALETTE["red"], ls="--", lw=1)
    ax[2].text(t_r + 0.5, 0.3, "reversal (Loschmidt)", color=PALETTE["red"])
    ax[2].set_xlabel("step"), ax[2].set_ylabel("coarse-grained S")
    ax[2].set_title("Exact reversal makes the entropy fall back to its initial value")
    save(fig, "fig20_irreversibility", FIG_DIR)
    out["cat_map"] = {"lyapunov": TR.LYAPUNOV_CAT, "S_final": float(res["S"][-1]), "S_max": res["S_max"],
                      "float_err_t10": float(res["err_float"][10]), "float_err_t40": float(res["err_float"][40]),
                      "integer_err_max": int(res["err_int"].max()), "S_after_roundtrip": float(S_path[-1])}

    # ---- Figure 21: wave vs heat equation, quantum reversal ----
    steps = np.arange(5, 130, 5)
    werr = [TR.wave_reversal_error(steps=int(s)) for s in steps]
    herr = [TR.heat_reversal_error(steps=int(s)) for s in steps]
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    ax[0].semilogy(steps, np.maximum(werr, 1e-17), "o-", ms=3, color=PALETTE["blue"], label="wave equation (second order in t)")
    ax[0].semilogy(steps, np.maximum(herr, 1e-17), "s-", ms=3, color=PALETTE["red"], label="heat equation (first order, dissipative)")
    ax[0].set_xlabel("steps before reversal"), ax[0].set_ylabel("error after round trip")
    ax[0].set_title("Reversible dynamics versus diffusion"), ax[0].legend(fontsize=8)
    Ts = np.linspace(0.1, 6, 40)
    err_th = [TR.quantum_reversal_check([0.3, 0.5, 0.7], t, np.array([0.6, 0.8j]))[0] for t in Ts]
    ax[1].semilogy(Ts, np.maximum(err_th, 1e-18), "o-", ms=3, color=PALETTE["green"])
    ax[1].set_xlabel("t"), ax[1].set_ylabel("|$\\Theta U(t)\\Theta^{-1}$ - exp(+i H' t)|")
    ax[1].set_title("Spin-1/2 in a field: $\\Theta$ maps H to H' = -H (field reversed)")
    save(fig, "fig21_wave_heat_quantum", FIG_DIR)
    out["wave_error_max"] = float(max(werr))
    out["heat_error_at_120_steps"] = float(herr[-1])

    # ---- Figure 22: black hole and white hole as time reverses ----
    M = 1.0
    tr = G.integrate_infall(8.0, 1.0, r_end=1e-3, n=400)
    # build the Kruskal trajectory: r(tau) known, v = t + r*; U = -exp(-u/4M), V = exp(v/4M) for region I; extend
    r, v = tr["r"], tr["v"]
    V = np.exp(v / (4 * M))
    X2mT2 = (r / (2 * M) - 1) * np.exp(r / (2 * M))            # = UV' sign convention: X^2 - T^2
    U = -X2mT2 / V                                             # U V = T^2 - X^2 = -(r/2M-1) e^{r/2M}
    Tk, Xk = (V + U) / 2, (V - U) / 2
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    for k, axi in enumerate(ax):
        D.draw_penrose_frame(axi)
        sgn = 1 if k == 0 else -1
        Tp, Xp = C.to_penrose(sgn * Tk, Xk)
        axi.plot(Xp, Tp, color=PALETTE["orange"], lw=2.2)
        axi.annotate("", xy=(Xp[len(Xp) // 2 + 4], Tp[len(Tp) // 2 + 4]), xytext=(Xp[len(Xp) // 2], Tp[len(Tp) // 2]),
                     arrowprops=dict(arrowstyle="->", color=PALETTE["orange"]))
    ax[0].set_title("Infalling observer: black-hole region II")
    ax[1].set_title("Time-reversed trajectory (T -> -T): emerges from white-hole region IV")
    save(fig, "fig22_black_white_time_reversal", FIG_DIR)
    out["kruskal_time_reversal_invariance_error"] = TR.kruskal_time_reversal_invariance(
        np.array([0.2, 0.5, -0.3]), np.array([0.4, 0.1, 0.8]))
    return out


if __name__ == "__main__":
    print(record("exp07_time_reversal", run()))
