"""Tidal tensor and geodesic deviation for radial Schwarzschild infall.

In the orthonormal frame (u, e1 = radial, e2, e3 = transverse) comoving with a
radial geodesic, the electric Weyl components are independent of E:
    E_11 = R_{1010} = -2M/r^3,   E_22 = E_33 = R_{2020} = +M/r^3.
Geodesic deviation d^2 xi^i/dtau^2 = -E_ij xi^j gives radial stretching
(+2M/r^3 xi_1) and transverse compression (-M/r^3 xi_perp).
"""
from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from . import symbolic as S
from .conventions import C_SI, G_SI, geometric_mass_m, tidal_accel_si
from .geodesics import integrate_infall

_cache = {}


def _lower_riemann_func():
    """Lambdified lower-index Riemann tensor of Schwarzschild, cached after the first call."""
    if "R" not in _cache:
        sch = S.schwarzschild()
        Rl = sch.riemann_lower()
        args = (S.r, S.th, S.M)
        _cache["R"] = [[[[sp.lambdify(args, Rl[a][b][c][d], "numpy") for d in range(4)] for c in range(4)]
                        for b in range(4)] for a in range(4)]
    return _cache["R"]


def frame_tidal_tensor(r, E, M=1.0, theta=np.pi / 2, sign=-1.0):
    """Tidal tensor E_ij = R_{i0j0} in the comoving frame, built from the full
    Riemann tensor (symbolic, lambdified), independent of the closed form."""
    R = _lower_riemann_func()
    f = 1 - 2 * M / r
    rd = sign * np.sqrt(E**2 - f)
    u = np.array([E / f, rd, 0.0, 0.0])
    e1 = np.array([rd / f, E, 0.0, 0.0])
    e2 = np.array([0.0, 0.0, 1.0 / r, 0.0])
    e3 = np.array([0.0, 0.0, 0.0, 1.0 / (r * np.sin(theta))])
    es = [e1, e2, e3]
    Rn = np.empty((4, 4, 4, 4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    Rn[a, b, c, d] = R[a][b][c][d](r, theta, M)
    T = np.empty((3, 3))
    for i in range(3):
        for j in range(3):
            T[i, j] = np.einsum("abcd,a,b,c,d->", Rn, es[i], u, es[j], u)
    return T


def tidal_closed_form(r, M=1.0):
    """Closed-form tidal tensor diag(-2M/r^3, M/r^3, M/r^3) in the comoving orthonormal frame."""
    return np.diag([-2 * M / r**3, M / r**3, M / r**3])


def integrate_deviation(r0, E, M=1.0, xi0=(1.0, 1.0), xidot0=(0.0, 0.0), r_end=1e-2, n=2000):
    """Integrate geodesic deviation (radial, transverse) along a radial infall."""
    traj = integrate_infall(r0, E, M, r_end=r_end)
    tau_end = traj["tau_end"]

    def rhs(tau, y):
        r = y[0]
        return [y[1], -M / r**2, y[3], 2 * M / r**3 * y[2], y[5], -M / r**3 * y[4]]

    y0 = [r0, traj["rdot"][0], xi0[0], xidot0[0], xi0[1], xidot0[1]]
    sol = solve_ivp(rhs, [0, tau_end], y0, method="DOP853", rtol=1e-11, atol=1e-13, dense_output=True)
    tau = np.linspace(0, tau_end, n)
    y = sol.sol(tau)
    return {"tau": tau, "r": y[0], "xi_r": y[2], "xi_perp": y[4], "tau_end": tau_end}


def neighbouring_geodesic_separation(r0, E, dE, M=1.0, n=400, r_end=0.05):
    """Independent validation: proper radial separation (orthogonal to u) between
    two exact radial geodesics with energies E and E+dE released from the same
    event; compared in tests with the solution of the deviation equation."""
    def run(En):
        return integrate_infall(r0, En, M, r_end=r_end, n=n)

    a, b = run(E), run(E + dE)
    tau_end = min(a["tau_end"], b["tau_end"]) * 0.999
    tau = np.linspace(0, tau_end, n)
    from scipy.interpolate import interp1d
    ra, va = interp1d(a["tau"], a["r"]), interp1d(a["tau"], a["v"])
    rb, vb = interp1d(b["tau"], b["r"]), interp1d(b["tau"], b["v"])
    rda = interp1d(a["tau"], a["rdot"])
    out = []
    for ta in tau:
        r_, v_, rd_ = float(ra(ta)), float(va(ta)), float(rda(ta))
        f_ = 1 - 2 * M / r_
        vd = 1.0 / (E - rd_)
        dv, dr = (float(vb(ta)) - v_) / dE, (float(rb(ta)) - r_) / dE
        # EF metric: -f dv^2 + 2 dv dr.  Remove component along u = (vd, rd_).
        g = lambda p, q: -f_ * p[0] * q[0] + p[0] * q[1] + p[1] * q[0]
        u = (vd, rd_)
        xi = (dv, dr)
        xi_perp = (xi[0] + g(xi, u) * u[0], xi[1] + g(xi, u) * u[1])  # u.u = -1
        out.append(np.sqrt(max(g(xi_perp, xi_perp), 0.0)))
    return tau, np.array(out) * dE


def spaghettification_radius(M_solar, body_length_m=2.0, limit_g=10.0):
    """Areal radius (m) where tidal acceleration across body_length_m equals
    limit_g * g_earth, and the ratio to the horizon radius."""
    from .conventions import G_EARTH, M_SUN_KG
    Mg = geometric_mass_m(M_solar * M_SUN_KG)
    target = limit_g * G_EARTH
    r = (C_SI**2 * 2 * Mg * body_length_m / target) ** (1 / 3)
    return r, r / (2 * Mg)
