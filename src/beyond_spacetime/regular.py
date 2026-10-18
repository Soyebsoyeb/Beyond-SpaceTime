"""Hayward regular black hole as a documented *effective* (non-fundamental) model.

    f(r) = 1 - 2 m(r)/r,   m(r) = M r^3 / (r^3 + 2 M l^2),   ds^2 = -f dt^2 + dr^2/f + r^2 dOmega^2.

This is a phenomenological metric (Hayward 2006).  It is NOT derived from a
fundamental quantum-gravity theory here; its effective stress-energy is read off
from the Einstein tensor, G^mu_nu = 8 pi T^mu_nu, and is diagnosed against the
classical energy conditions.  A finite curvature invariant is not claimed to be a
complete resolution of the singularity problem.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def m_of_r(r, M=1.0, l=0.5):
    """Hayward mass function m(r) = M r^3 / (r^3 + 2 M l^2)."""
    r = np.asarray(r, dtype=float)
    return M * r**3 / (r**3 + 2 * M * l**2)


def f(r, M=1.0, l=0.5):
    """Hayward metric function f(r) = 1 - 2 M r^2 / (r^3 + 2 M l^2)."""
    r = np.asarray(r, dtype=float)
    return 1.0 - 2.0 * M * r**2 / (r**3 + 2.0 * M * l**2)


def fprime(r, M=1.0, l=0.5):
    """Analytic derivative f'(r) of the Hayward metric function."""
    r = np.asarray(r, dtype=float)
    a = 2.0 * M * l**2
    return -2.0 * M * (2 * r * (r**3 + a) - 3 * r**4) / (r**3 + a) ** 2


def fsecond(r, M=1.0, l=0.5, h=1e-5):
    """Second derivative f''(r) by central finite difference of fprime (step h)."""
    return (fprime(r + h, M, l) - fprime(r - h, M, l)) / (2 * h)


_K_FUNC = None


def kretschmann(r, M=1.0, l=0.5):
    """K = f''^2 + 4 f'^2/r^2 + 4 (1-f)^2/r^4, simplified symbolically then lambdified so
    that r -> 0 is evaluated without cancellation (formula validated against the full
    Riemann computation in tests)."""
    global _K_FUNC
    if _K_FUNC is None:
        import sympy as sp
        from . import symbolic as S
        K = sp.simplify(S.kretschmann_static_formula(S.hayward_f()))
        _K_FUNC = sp.lambdify((S.r, S.M, S.l), K, "numpy")
    return _K_FUNC(np.asarray(r, dtype=float), M, l)


def kretschmann_core_limit(l):
    """de Sitter core: K(r -> 0) = 24 / l^4."""
    return 24.0 / l**4


def l_critical_over_M():
    """Extremal case (double horizon): l/M = 4/(3 sqrt 3)."""
    return 4.0 / (3.0 * np.sqrt(3.0))


def horizons(M=1.0, l=0.5):
    """Positive real roots of r^3 - 2M r^2 + 2 M l^2 = 0, ascending."""
    roots = np.roots([1.0, -2.0 * M, 0.0, 2.0 * M * l**2])
    pos = sorted(float(x.real) for x in roots if abs(x.imag) < 1e-9 and x.real > 0)
    return pos


def density(r, M=1.0, l=0.5):
    """rho = m'(r)/(4 pi r^2) = 3 M a / (4 pi (r^3 + a)^2), a = 2 M l^2."""
    r = np.asarray(r, dtype=float)
    a = 2.0 * M * l**2
    return 3.0 * M * a / (4 * np.pi * (r**3 + a) ** 2)


def stress_energy(r, M=1.0, l=0.5):
    """Effective rho, p_r, p_t from G^mu_nu = 8 pi T^mu_nu (anisotropic fluid)."""
    r = np.asarray(r, dtype=float)
    a = 2.0 * M * l**2
    rho = density(r, M, l)
    drho = -6.0 * 3.0 * M * a * r**2 / (4 * np.pi * (r**3 + a) ** 3)  # d rho/dr
    pr = -rho
    pt = -rho - 0.5 * r * drho
    return rho, pr, pt


def energy_conditions(r, M=1.0, l=0.5):
    """Pointwise booleans for NEC, WEC, SEC, DEC of the effective fluid."""
    rho, pr, pt = stress_energy(r, M, l)
    nec = (rho + pr >= -1e-14) & (rho + pt >= -1e-14)
    wec = nec & (rho >= 0)
    sec = nec & (rho + pr + 2 * pt >= -1e-14)
    dec = (rho >= 0) & (np.abs(pr) <= rho + 1e-14) & (np.abs(pt) <= rho + 1e-14)
    return {"NEC": nec, "WEC": wec, "SEC": sec, "DEC": dec}


def hawking_temperature(M=1.0, l=0.5):
    """T = f'(r_+)/(4 pi) at the outer horizon (0 at extremality, 1/(8 pi M) as l -> 0)."""
    rh = horizons(M, l)
    if not rh:
        return np.nan
    return float(fprime(rh[-1], M, l)) / (4 * np.pi)


def tau_outer_horizon_to_turning_point(E, M=1.0, l=0.5):
    """Proper time from r_+ to the radial turning point f(r_t) = E^2 inside (E^2 < 1),
    and the radius r_t.  Uses dtau = dr / sqrt(E^2 - f)."""
    rh = horizons(M, l)
    r_plus = rh[-1]
    g = lambda r: f(r, M, l) - E**2
    # turning point lies between r = 0 (g = 1 - E^2 > 0) and the inner horizon (g < 0)
    r_in = rh[0]
    rt = brentq(g, 1e-9, r_in)
    integrand = lambda r: 1.0 / np.sqrt(max(E**2 - f(r, M, l), 1e-300))
    val, _ = quad(integrand, rt, r_plus, limit=400, points=[r_in])
    return val, rt


def tau_to_center(E, M=1.0, l=0.5):
    """For E >= 1 the radial geodesic reaches the regular centre r = 0; proper time from r_+."""
    rh = horizons(M, l)
    integrand = lambda r: 1.0 / np.sqrt(E**2 - f(r, M, l))
    val, _ = quad(integrand, 0.0, rh[-1], limit=400)
    return val
