"""Experiment 10: dependency diagram linking equations, code modules and figures (docs/verification_matrix.md)."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from beyond_spacetime.plotting import PALETTE, save

from common import FIG_DIR

MATH = [("E1-E2", "Vacuum and Kretschmann"), ("E3-E5", "Kruskal and Penrose"), ("E6-E9", "Geodesics and tides"),
        ("E10-E12", "Bridge and causality"), ("E13-E18", "Hayward and QNMs"), ("E19-E25", "BTZ, RT and TFD"),
        ("E26-E27", "Decoherence"), ("E28-E30", "Time reversal"), ("E31-E35", "Wormholes and CTCs")]
# (math index, module label, figure index)
MODULES = ["symbolic.py", "coordinates.py", "geodesics.py + tides.py (*)", "embedding.py + causal.py",
           "regular.py (*) + perturbations.py", "holography.py + tfd.py", "quantum.py", "timereversal.py",
           "wormholes.py (*) + chronology.py"]
FIGS = ["fig01-02", "fig01-02", "fig03-07", "fig08-10", "fig11-14, 26", "fig15-17", "fig18-19", "fig20-22", "fig23-25"]


def run():
    fig, ax = plt.subplots(figsize=(13, 6.4))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 10.6)
    ax.axis("off")
    ys = [9.2 - 1.0 * i for i in range(9)]
    cols = {"math": (0.2, 3.6), "mod": (4.7, 3.8), "fig": (10.4, 2.4)}
    colour = {"math": PALETTE["blue"], "mod": PALETTE["orange"], "fig": PALETTE["green"]}
    titles = {"math": "Mathematical results (verification matrix IDs)", "mod": "Code modules (src/beyond_spacetime)",
              "fig": "Figures (results/figures)"}
    for k, (x, w) in cols.items():
        ax.text(x + w / 2, 10.2, titles[k], ha="center", va="center", fontsize=10, fontweight="bold")
    for i, y in enumerate(ys):
        items = {"math": f"{MATH[i][0]}  {MATH[i][1]}", "mod": MODULES[i], "fig": FIGS[i]}
        for k, (x, w) in cols.items():
            ax.add_patch(FancyBboxPatch((x, y - 0.32), w, 0.64, boxstyle="round,pad=0.02,rounding_size=0.1",
                                        fc=colour[k], ec="none", alpha=0.18))
            ax.add_patch(FancyBboxPatch((x, y - 0.32), w, 0.64, boxstyle="round,pad=0.02,rounding_size=0.1",
                                        fc="none", ec=colour[k], lw=1.2))
            ax.text(x + w / 2, y, items[k], ha="center", va="center", fontsize=9)
        ax.annotate("", xy=(cols["mod"][0], y), xytext=(cols["math"][0] + cols["math"][1], y),
                    arrowprops=dict(arrowstyle="->", color="0.3"))
        ax.annotate("", xy=(cols["fig"][0], y), xytext=(cols["mod"][0] + cols["mod"][1], y),
                    arrowprops=dict(arrowstyle="->", color="0.3"))
    ax.text(6.6, 0.35, "(*) validated against the independent symbolic curvature pipeline in symbolic.py",
            ha="center", fontsize=9, color=PALETTE["red"])
    save(fig, "fig28_dependency_graph", FIG_DIR)
    return {"nodes": len(MATH) + len(MODULES) + len(set(FIGS))}

