"""Morris-Thorne traversable wormholes and the Ellis-Bronnikov example.

ds^2 = -e^{2 Phi(r)} dt^2 + dr^2/(1 - b(r)/r) + r^2 dOmega^2.
Throat at r0: b(r0) = r0;  flare-out: b'(r0) < 1.  Einstein equations G = 8 pi T give
    rho = b'/(8 pi r^2),
    tau = [b/r - 2 (r - b) Phi'] / (8 pi r^2)   (radial tension, p_r = -tau),
    p   = (r/2) [(rho - tau) Phi' - tau'] - tau   (tangential pressure).
At the throat  rho - tau = (b'(r0) - 1)/(8 pi r0^2) < 0: radial NEC violation.
Ellis-Bronnikov: Phi = 0, b = b0^2/r, equivalently ds^2 = -dt^2 + dl^2 + (b0^2 + l^2) dOmega^2.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from . import symbolic as S

rr = sp.Symbol("r", positive=True)


def mt_metric(b, Phi):
    """Spacetime object for the Morris-Thorne metric with shape function b(r) and redshift function Phi(r)."""
    f = 1 - b / S.r
    return S.Spacetime([S.t, S.r, S.th, S.ph],
                       sp.diag(-sp.exp(2 * Phi), 1 / f, S.r**2, S.r**2 * sp.sin(S.th) ** 2))


def mt_orthonormal_einstein(b, Phi):
    """Orthonormal-frame Einstein tensor components (G_tt, G_rr, G_thth) computed from the metric."""
    sp_ = mt_metric(b, Phi)
    G = sp_.einstein()
    f = 1 - b / S.r
    gtt = sp.exp(2 * Phi)
    Gtt = sp.simplify(G[0, 0] / gtt)
    Grr = sp.simplify(G[1, 1] * f)
    Gthth = sp.simplify(G[2, 2] / S.r**2)
    return Gtt, Grr, Gthth


def mt_matter(b, Phi):
    """(rho, tau, p) from the closed-form Morris-Thorne expressions."""
    r = S.r
    rho = sp.diff(b, r) / (8 * sp.pi * r**2)
    tau = (b / r - 2 * (r - b) * sp.diff(Phi, r)) / (8 * sp.pi * r**2)
    p = r / 2 * ((rho - tau) * sp.diff(Phi, r) - sp.diff(tau, r)) - tau
    return rho, tau, p


def ellis_b(b0):
    """Ellis-Bronnikov shape function b(r) = b0^2 / r."""
    return b0**2 / S.r


def nec_radial_at_throat(b0):
    """rho - tau at the throat for the Ellis wormhole: -1/(4 pi b0^2)."""
    r0 = b0
    return (-(b0**2) / r0**2 - 1.0) / (8 * np.pi * r0**2)


def ellis_areal_radius(l, b0=1.0):
    """Areal radius sqrt(b0^2 + l^2) as a function of proper length l."""
    return np.sqrt(b0**2 + np.asarray(l, dtype=float) ** 2)


def ellis_embedding_z(l, b0=1.0):
    """Embedding height z(l) = sign(l) b0 arccosh(r/b0) of the Ellis wormhole (dz/dr = 1/sqrt(r^2/b0^2 - 1)); the sign selects the side."""
    r = ellis_areal_radius(l, b0)
    return np.sign(l) * b0 * np.arccosh(r / b0)


def ellis_tidal_lateral(l, v, b0=1.0):
    """Lateral tidal component gamma^2 v^2 R_{l theta l theta} felt by a radial traveller of speed v at proper length l.

    The orthonormal-frame curvature is R_{l theta l theta} = -b0^2 / (b0^2 + l^2)^2 and g_tt = -1 gives no
    time-time contribution, so a static observer (v = 0) feels no lateral tide. Checked against the symbolic
    Riemann tensor by ellis_boosted_tidal_symbolic.
    """
    g2 = 1.0 / (1.0 - v**2)
    Rltlt = -b0**2 / (b0**2 + np.asarray(l, dtype=float) ** 2) ** 2
    return g2 * v**2 * Rltlt


def ellis_boosted_tidal_symbolic():
    """Independent symbolic check: orthonormal-frame lateral tidal component of a boosted radial observer in the Ellis metric, returned with its symbols (l, b0, v)."""
    l, b0, v = sp.symbols("l b0 v", positive=True)
    th = S.th
    st = S.Spacetime([S.t, l, th, S.ph], sp.diag(-1, 1, b0**2 + l**2, (b0**2 + l**2) * sp.sin(th) ** 2))
    Rl = st.riemann_lower()
    g = 1 / sp.sqrt(1 - v**2)
    u = [g, g * v, 0, 0]
    el = [g * v, g, 0, 0]
    eth = [0, 0, 1 / sp.sqrt(b0**2 + l**2), 0]
    comp = sum(Rl[a][b][c][d] * eth[a] * u[b] * eth[c] * u[d]
               for a in range(4) for b in range(4) for c in range(4) for d in range(4))
    return sp.simplify(comp), (l, b0, v)
