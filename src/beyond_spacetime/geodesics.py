"""Radial timelike geodesics of Schwarzschild.

First integral (E = -u_t = f dt/dtau):  (dr/dtau)^2 = E^2 - f(r).
Differentiating gives d^2 r/dtau^2 = -M/r^2 (valid for every E and through the
horizon).  Ingoing Eddington-Finkelstein v = t + r* obeys dv/dtau = 1/(E - rdot),
regular at r = 2M.  Analytic references:
  * E = 1:            tau(r1 -> r2) = 2 (r1^{3/2} - r2^{3/2}) / (3 sqrt(2M))
  * rest at r0 (E^2 = f(r0)): cycloid r = (r0/2)(1 + cos eta),
        tau = (r0/2) sqrt(r0/2M) (eta + sin eta)
  * horizon -> singularity proper time = integral_0^{2M} dr / sqrt(E^2 - f),
        maximal value pi M at E = 0.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad, solve_ivp

from .coordinates import f as fS, tortoise


def tau_inside(E, M=1.0):
    """Proper time from r = 2M to r = 0 for conserved energy E (E^2 >= 0)."""
    E2 = float(E) ** 2
    g = lambda r: 1.0 / np.sqrt(E2 - (1 - 2 * M / r))
    val, _ = quad(g, 0.0, 2 * M, limit=400, epsabs=1e-13, epsrel=1e-13)
    return val


def tau_E1(r1, r2, M=1.0):
    """Analytic proper time between radii r1 > r2 for the E = 1 radial geodesic: 2 (r1^{3/2} - r2^{3/2}) / (3 sqrt(2M))."""
    return 2.0 * (r1**1.5 - r2**1.5) / (3.0 * np.sqrt(2 * M))


def tau_rest_at(r0, r, M=1.0):
    """Proper time from rest at r0 to radius r (cycloid solution), r <= r0."""
    eta = np.arccos(2 * r / r0 - 1.0)
    return 0.5 * r0 * np.sqrt(r0 / (2 * M)) * (eta + np.sin(eta))


def tau_rest_to_singularity(r0, M=1.0):
    """Proper time from rest at r0 to r = 0 (cycloid): (pi/2) sqrt(r0^3 / 2M)."""
    return 0.5 * np.pi * np.sqrt(r0**3 / (2 * M))


def energy_from_rest(r0, M=1.0):
    """Conserved energy E = sqrt(f(r0)) of a geodesic released from rest at r0."""
    return float(np.sqrt(fS(r0, M)))


def integrate_infall(r0, E, M=1.0, r_end=1e-3, rtol=1e-12, atol=1e-14, n=2000):
    """Integrate r'' = -M/r^2 from r0 with rdot = -sqrt(E^2 - f(r0)); also v and tau.

    Returns dict with arrays tau, r, rdot, v and the horizon-crossing proper time.
    Requires E^2 >= f(r0).
    """
    rdot0 = -np.sqrt(max(E**2 - fS(r0, M), 0.0))

    def rhs(tau, y):
        r, rd, v = y
        return [rd, -M / r**2, 1.0 / (E - rd)]

    horizon = lambda tau, y: y[0] - 2 * M
    horizon.direction = -1
    stop = lambda tau, y: y[0] - r_end
    stop.terminal = True
    v0 = 0.0
    sol = solve_ivp(rhs, [0, 1e4], [r0, rdot0, v0], method="DOP853", events=[horizon, stop],
                    rtol=rtol, atol=atol, dense_output=True)
    tau_h = sol.t_events[0][0] if len(sol.t_events[0]) else np.nan
    tau_end = sol.t[-1]
    tau = np.linspace(0, tau_end, n)
    y = sol.sol(tau)
    return {"tau": tau, "r": y[0], "rdot": y[1], "v": y[2], "tau_horizon": tau_h, "tau_end": tau_end,
            "r_end": r_end, "E": E, "M": M}


def schwarzschild_time(r, v, M=1.0):
    """t = v - r* (valid outside the horizon only)."""
    return v - tortoise(r, M)


def kretschmann(r, M=1.0):
    """Schwarzschild Kretschmann scalar 48 M^2 / r^6."""
    return 48 * M**2 / np.asarray(r) ** 6
