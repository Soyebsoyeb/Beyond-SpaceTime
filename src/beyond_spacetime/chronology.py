"""Chronology violation: a wormhole time machine (Morris-Thorne-Yurtsever type
identification model) and the Goedel universe.

Wormhole model (c = 1, throat length ~ 0): mouths A, B at lab separation D.  Mouth B makes a
round trip at speed beta for lab time T_trip; its internal clock then lags the lab by
    Delta = T_trip (1 - sqrt(1 - beta^2)).
Entering the lagging mouth B at lab time t emerges from A at lab time t - Delta (time advance);
the reverse direction is a time delay.  A closed timelike curve through the wormhole and the
exterior exists iff the exterior path A -> B can be traversed in lab time strictly less than ...
precisely: iff Delta > D (exterior transit needs time >= D, through-wormhole advance is Delta).
The first CTC appears at the chronology horizon.  Existence of a CTC in an *assumed* identification
model is mathematical; physical realizability requires exotic matter and is expected to fail
because of semiclassical stress-energy divergences at the horizon (chronology protection
conjecture; Hawking 1992; Kim and Thorne 1991).

Goedel (1949): in cylindrical coordinates ds^2 = 4a^2 [ -dt^2 + dr^2 - (sinh^4 r - sinh^2 r) dphi^2
- 2 sqrt2 sinh^2 r dt dphi + dz^2 ]; circles r = const have g_phiphi < 0 for sinh r > 1, i.e.
r > r_c = ln(1 + sqrt 2): closed timelike curves.
"""
from __future__ import annotations

import numpy as np


def clock_lag(T_trip, beta):
    """Clock offset Delta = T (1 - sqrt(1 - beta^2)) accumulated by a wormhole mouth that travels at speed beta (c = 1) for lab time T."""
    return T_trip * (1.0 - np.sqrt(1.0 - beta**2))


def ctc_exists(delta, D):
    """True iff a closed timelike curve exists through the wormhole, i.e. the clock lag Delta exceeds the mouth separation D."""
    return delta > D


def chronology_horizon_time(D, beta):
    """Trip duration T at which the clock lag first equals the mouth separation D (solves Delta = D), the onset of closed timelike curves.

    The mouth separation at the end of the trip is taken to equal D.
    """
    return D / (1.0 - np.sqrt(1.0 - beta**2))


def ctc_map(D_values, beta_values, T_trip):
    """Boolean array over (beta, D) that is True where Delta(T_trip, beta) > D (closed timelike curves exist)."""
    out = np.zeros((len(beta_values), len(D_values)), dtype=bool)
    for i, b in enumerate(beta_values):
        d = clock_lag(T_trip, b)
        out[i] = d > np.asarray(D_values)
    return out


def minimal_loop_time(delta, D):
    """For Delta > D the closed loop consists of a wormhole transit (advance Delta) and an exterior walk
    of lab duration Delta (> D, so timelike). Returns the exterior walking time needed and the mean speed."""
    if delta <= D:
        return None
    return delta, D / delta


# ---- Goedel -----

R_CRIT_GODEL = float(np.log(1.0 + np.sqrt(2.0)))


def godel_gphiphi(r, a=1.0):
    """g_{phi phi} = 4 a^2 (sinh^2 r - sinh^4 r)."""
    s2 = np.sinh(r) ** 2
    return 4.0 * a**2 * (s2 - s2**2)


def godel_circle_is_closed_timelike(r, a=1.0):
    """True where the circle of radius r in the Goedel universe is a closed timelike curve (g_phiphi < 0)."""
    return godel_gphiphi(r, a) < 0.0
