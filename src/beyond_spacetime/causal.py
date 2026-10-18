"""Causal structure of the maximally extended Schwarzschild solution and of
a traversable (Ellis-Bronnikov) wormhole.

In Kruskal null coordinates u = T - X, v = T + X the radial (and, by conformal
flatness of the (T, X) sector, the projection of every) future-directed causal
curve satisfies du >= 0 and dv >= 0.  Hence q is in the causal future of p only
if u_q >= u_p and v_q >= v_p -- a necessary condition used here to exclude
signalling.  Regions I -> III would require v to decrease from v>0 to v<0.
"""
from __future__ import annotations

import numpy as np

from .coordinates import region_of


def null_coords(T, X):
    """Kruskal null coordinates (u, v) = (T - X, T + X)."""
    return np.asarray(T) - np.asarray(X), np.asarray(T) + np.asarray(X)


def may_be_in_causal_future(p, q):
    """Necessary condition for q in J+(p); p, q are (T, X) pairs."""
    up, vp = null_coords(*p)
    uq, vq = null_coords(*q)
    return bool((uq >= up) and (vq >= vp))


def sample_region(region, n, rng, M=1.0, tmax=3.0, xmax=3.0):
    """Random Kruskal points with T^2 - X^2 > -1 (r > 0) inside a region."""
    pts = []
    while len(pts) < n:
        T = rng.uniform(-tmax, tmax)
        X = rng.uniform(-xmax, xmax)
        if T**2 - X**2 >= 1.0:
            continue  # beyond the singularity
        if str(region_of(T, X)) == region:
            pts.append((T, X))
    return pts


def count_signalling_pairs(n, seed=0):
    """Count pairs (p in I, q in III) with q possibly in J+(p): must be zero."""
    rng = np.random.default_rng(seed)
    P = sample_region("I", n, rng)
    Q = sample_region("III", n, rng)
    return sum(may_be_in_causal_future(p, q) for p in P for q in Q), len(P) * len(Q)


def ellis_crossing_time(l_start, l_end):
    """Radial null crossing of the Ellis-Bronnikov wormhole ds^2 = -dt^2 + dl^2 + (b0^2+l^2) dOmega^2:
    dl/dt = +1, so the coordinate time is simply l_end - l_start (> 0 for l_end > l_start)."""
    return l_end - l_start


def causal_relation_table():
    """Structured summary used in docs and the README table."""
    return [
        {
            "geometry": "Eternal Schwarzschild (Kruskal extension)",
            "I_to_III_causal": False,
            "reason": "Bridge pinches off: v must decrease from >0 to <0 along a future-directed causal curve",
            "assumptions": "Exact vacuum, eternal (past of region IV exists), analytic extension",
            "matter": "None",
        },
        {
            "geometry": "Collapse (Vaidya thin shell onto flat interior)",
            "I_to_III_causal": None,
            "reason": "Regions III and IV are absent: the past of the shell is flat space",
            "assumptions": "Spherical collapse, null shell at v = v0",
            "matter": "Null dust shell",
        },
        {
            "geometry": "Ellis-Bronnikov wormhole",
            "I_to_III_causal": True,
            "reason": "Static throat, radial null rays have dl/dt = 1",
            "assumptions": "Radial NEC violated (ghost scalar); linearly unstable",
            "matter": "Exotic (phantom) scalar field",
        },
        {
            "geometry": "AdS2/BTZ eternal bridge (no coupling)",
            "I_to_III_causal": False,
            "reason": "Same Kruskal causal structure as Schwarzschild",
            "assumptions": "Two decoupled CFTs in the thermofield-double state",
            "matter": "Negative cosmological constant",
        },
        {
            "geometry": "Traversable holographic wormhole (double-trace coupling)",
            "I_to_III_causal": True,
            "reason": "Boundary coupling supplies negative averaged null energy; effect is quantum and small",
            "assumptions": "Gao-Jafferis-Wall setup; semiclassical bulk; transit time exceeds the free-signal time outside",
            "matter": "Bulk quantum fields with induced ANEC violation",
        },
    ]
