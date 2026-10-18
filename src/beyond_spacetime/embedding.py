"""Time-symmetric spatial slice, Flamm paraboloid and comparison of slicings.

Equatorial t = const slice: dl^2 = dr^2/(1 - 2M/r) + r^2 dphi^2.
Flamm embedding in R^3:  z(r) = 2 sqrt(2M (r - 2M)),  (dz/dr)^2 = 2M/(r - 2M).
Two sheets z = +/- z(r) glue at the throat r = 2M (area 16 pi M^2).
Isotropic coordinate rho: r = rho (1 + M/(2 rho))^2, throat at rho = M/2 and
inversion isometry rho -> M^2/(4 rho).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad

from .coordinates import kruskal_conformal_factor, kruskal_to_r


def flamm_z(r, M=1.0):
    """Flamm embedding height z(r) = 2 sqrt(2M (r - 2M)) of the t = const equatorial slice."""
    return 2.0 * np.sqrt(2.0 * M * (np.asarray(r, dtype=float) - 2.0 * M))


def flamm_dzdr(r, M=1.0):
    """Slope dz/dr = sqrt(2M/(r - 2M)) of the Flamm embedding."""
    return np.sqrt(2.0 * M / (np.asarray(r, dtype=float) - 2.0 * M))


def proper_radial_distance(r, M=1.0):
    """Proper distance from the throat along the t = const slice (closed form)."""
    r = np.asarray(r, dtype=float)
    return np.sqrt(r * (r - 2 * M)) + 2 * M * np.log((np.sqrt(r) + np.sqrt(r - 2 * M)) / np.sqrt(2 * M))


def isotropic_r(rho, M=1.0):
    """Areal radius r = rho (1 + M/(2 rho))^2 in terms of the isotropic coordinate rho."""
    return rho * (1.0 + M / (2.0 * rho)) ** 2


def kruskal_slice_profile(T0, M=1.0, Xmax=4.0, n=4001):
    """Areal radius r and proper distance l along the slice T = T0 (|T0| < 1 has a throat)."""
    X = np.linspace(-Xmax, Xmax, n)
    r = kruskal_to_r(np.full_like(X, T0), X, M)
    # induced metric on T = const: dl^2 = F dX^2
    dl = np.sqrt(kruskal_conformal_factor(r, M)) * (X[1] - X[0])
    l = np.concatenate([[0.0], np.cumsum(0.5 * (dl[1:] + dl[:-1]))])
    l -= l[n // 2]
    return X, r, l


def throat_radius_kruskal(T0, M=1.0):
    """Minimal areal radius on the slice T = T0 (at X = 0): (r/2M - 1) e^{r/2M} = -T0^2."""
    return float(kruskal_to_r(np.array(T0), np.array(0.0), M))


def painleve_gullstrand_profile(M=1.0, rmax=8.0, n=800):
    """PG slice: induced spatial metric is flat, r = l."""
    r = np.linspace(1e-3, rmax, n)
    return r, r.copy()
