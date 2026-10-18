"""AdS3/CFT2 holography with the non-rotating BTZ black hole (AdS radius ell = 1).

Bulk:   ds^2 = -(r^2 - r_+^2) dt^2 + dr^2/(r^2 - r_+^2) + r^2 dphi^2,  phi ~ phi + 2 pi.
Boundary theory: a 1+1 dimensional CFT on a circle with central charge c = 3 ell/(2 G_N)
(Brown-Henneaux).  State: thermal at T = r_+/(2 pi) on one side; the thermofield-double
on two copies for the two-sided geometry.  Regime: large c, classical gravity, intervals
whose RT geodesic is the connected one (no phase transition considered).

Ryu-Takayanagi in this static setting:  S_A = Length(gamma_A) / (4 G_N).
Cross-checks implemented here:
  1. numerical geodesic length vs closed form  L = 2 ln[(2 r_c/r_+) sinh(r_+ dphi/2)],
  2. agreement with the Calabrese-Cardy thermal formula  S = (c/3) ln[(beta/(pi eps)) sinh(pi l/beta)],
  3. exact AdS3 embedding-space geodesic distances,
  4. two-sided wormhole length growth L(t) = 2 ln[(2 r_c/r_+) cosh(r_+ t)] (Hartman-Maldacena),
  5. horizon entropy = Cardy entropy = thermal entropy of the boundary CFT.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def central_charge(G_N=1.0, ell=1.0):
    """Brown-Henneaux central charge c = 3 ell / (2 G_N)."""
    return 3.0 * ell / (2.0 * G_N)


def geodesic_J_from_angle(dphi, r_plus):
    """Turning-point radius r* = J = r_+ coth(r_+ dphi / 2)."""
    return r_plus / np.tanh(r_plus * dphi / 2.0)


def angle_from_J_numeric(J, r_plus):
    """dphi = 2 int_0^inf dx / (cosh x sqrt(J^2 cosh^2 x - r_+^2))."""
    g = lambda x: 1.0 / (np.cosh(x) * np.sqrt(J**2 * np.cosh(x) ** 2 - r_plus**2))
    val, _ = quad(g, 0.0, 60.0, limit=400, epsabs=1e-13, epsrel=1e-13)
    return 2.0 * val


def length_numeric(J, r_plus, r_cut):
    """Regulated geodesic length 2 int_0^{x_c} J cosh x / sqrt(J^2 cosh^2 x - r_+^2) dx."""
    import warnings
    xc = np.arccosh(r_cut / J)
    g = lambda x: J * np.cosh(x) / np.sqrt(J**2 * np.cosh(x) ** 2 - r_plus**2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")   # round-off warnings at the 1e-13 request; accuracy is checked in tests
        val, _ = quad(g, 0.0, xc, limit=400, epsabs=1e-13, epsrel=1e-13)
    return 2.0 * val


def geodesic_length_closed(dphi, r_plus, r_cut):
    """Closed-form regulated BTZ geodesic length 2 ln[(2 r_c / r_+) sinh(r_+ dphi / 2)] for an interval of angle dphi."""
    return 2.0 * np.log(2.0 * r_cut / r_plus * np.sinh(r_plus * dphi / 2.0))


def entanglement_entropy_RT(dphi, r_plus, r_cut, G_N=1.0):
    """Ryu-Takayanagi entropy of the connected geodesic, length / (4 G_N)."""
    return geodesic_length_closed(dphi, r_plus, r_cut) / (4.0 * G_N)


def entanglement_entropy_CFT(dphi, r_plus, r_cut, G_N=1.0):
    """Calabrese-Cardy: S = (c/3) ln[(beta/(pi eps)) sinh(pi l / beta)], beta = 2 pi/r_+, eps = 1/r_cut."""
    c = central_charge(G_N)
    beta = 2.0 * np.pi / r_plus
    eps = 1.0 / r_cut
    return c / 3.0 * np.log(beta / (np.pi * eps) * np.sinh(np.pi * dphi / beta))


def rt_entropy(dphi, r_plus, r_cut, G_N=1.0):
    """Full RT prescription for an interval A of angular size dphi in (0, 2 pi) on the thermal circle:
    S_A = min[ L_conn(dphi), L_conn(2 pi - dphi) + 2 pi r_+ ] / (4 G_N).
    The second candidate (the geodesic homologous to the complement plus the horizon) is required by the
    homology constraint and gives S_A - S_{A^c} = S_BH once it dominates (Araki-Lieb saturation).
    Valid for the BTZ saddle (r_+ > 1 is the Hawking-Page-dominant phase), large c, static slice."""
    conn = geodesic_length_closed(dphi, r_plus, r_cut)
    disc = geodesic_length_closed(2.0 * np.pi - dphi, r_plus, r_cut) + 2.0 * np.pi * r_plus
    return np.minimum(conn, disc) / (4.0 * G_N)


def transition_angle(r_plus):
    """Angle dphi* where the connected and horizon-wrapping candidates exchange dominance
    (cutoff independent):  sinh(r_+ dphi/2) = e^{pi r_+} sinh(pi r_+ - r_+ dphi/2)."""
    g = lambda d: np.log(np.sinh(r_plus * d / 2.0)) - np.log(np.sinh(r_plus * (2 * np.pi - d) / 2.0)) - np.pi * r_plus
    return brentq(g, np.pi, 2 * np.pi - 1e-9)


def vacuum_limit_length(dphi, r_cut):
    """r_+ -> 0 (planar vacuum) limit: L = 2 ln(r_cut dphi)."""
    return 2.0 * np.log(r_cut * dphi)


# ---- embedding-space geodesic distances (exact AdS3) ---------------------------------

def btz_embedding(r, t, phi, r_plus, side=+1):
    """Embedding coordinates (T1, T2, X1, X2) of BTZ exterior points, -T1^2 - T2^2 + X1^2 + X2^2 = -1.
    side=+1: right exterior, side=-1: left exterior with Kruskal time running upward."""
    a = r / r_plus
    b = np.sqrt(a**2 - 1.0)
    T1 = a * np.cosh(r_plus * phi)
    X1 = a * np.sinh(r_plus * phi)
    T2 = b * np.sinh(r_plus * t)
    X2 = side * b * np.cosh(r_plus * t)
    return np.array([T1, T2, X1, X2])


def embedding_distance(P, Q):
    """Geodesic distance between spacelike-separated points: cosh d = -P.Q."""
    dot = -P[0] * Q[0] - P[1] * Q[1] + P[2] * Q[2] + P[3] * Q[3]
    return np.arccosh(-dot)


def two_sided_length(t, r_plus, r_cut):
    """Exact regulated ER-bridge length for t_L = t_R = t (Kruskal-forward), same phi."""
    P = btz_embedding(r_cut, t, 0.0, r_plus, +1)
    Q = btz_embedding(r_cut, t, 0.0, r_plus, -1)
    return embedding_distance(P, Q)


def two_sided_length_closed(t, r_plus, r_cut):
    """Closed-form two-sided bridge length 2 ln[(2 r_c / r_+) cosh(r_+ t)]."""
    return 2.0 * np.log(2.0 * r_cut / r_plus * np.cosh(r_plus * t))


# ---- thermodynamics -------------------------------------------------------------------

def bh_entropy(r_plus, G_N=1.0):
    """S_BH = Area/(4 G) = 2 pi r_+ / (4 G)."""
    return 2.0 * np.pi * r_plus / (4.0 * G_N)


def cardy_entropy(r_plus, G_N=1.0):
    """High-temperature CFT entropy S = (pi c / 3) L T with L = 2 pi, T = r_+/(2 pi)."""
    c = central_charge(G_N)
    return np.pi * c / 3.0 * 2.0 * np.pi * (r_plus / (2.0 * np.pi))
