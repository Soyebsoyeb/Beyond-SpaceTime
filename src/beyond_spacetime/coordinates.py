"""Schwarzschild, Eddington-Finkelstein, Kruskal-Szekeres and Penrose coordinates.

Regions of the maximally extended solution (Kruskal T, X):
    I   : X >  |T|   right exterior          II : T >  |X|  future interior (black hole)
    III : X < -|T|   left exterior           IV : T < -|X|  past interior (white hole)
Singularities: T^2 - X^2 = 1 (II: T>0, IV: T<0). Kruskal null coordinates
U = T - X, V = T + X.  Penrose coordinates: u~ = arctan U, v~ = arctan V,
T_p = v~ + u~, X_p = v~ - u~.
"""
from __future__ import annotations

import numpy as np
from scipy.special import lambertw


def f(r, M=1.0):
    """Schwarzschild metric function f(r) = 1 - 2M/r."""
    return 1.0 - 2.0 * M / np.asarray(r, dtype=float)


def tortoise(r, M=1.0):
    """r* = r + 2M ln|r/2M - 1| (valid on both sides of the horizon)."""
    r = np.asarray(r, dtype=float)
    return r + 2 * M * np.log(np.abs(r / (2 * M) - 1.0))


def ingoing_ef_v(t, r, M=1.0):
    """Ingoing Eddington-Finkelstein advanced time v = t + r*."""
    return t + tortoise(r, M)


def outgoing_ef_u(t, r, M=1.0):
    """Outgoing Eddington-Finkelstein retarded time u = t - r*."""
    return t - tortoise(r, M)


def schwarzschild_to_kruskal(t, r, M=1.0, region="I"):
    """Map (t, r) in the stated region to Kruskal (T, X)."""
    t = np.asarray(t, dtype=float)
    r = np.asarray(r, dtype=float)
    ex = np.exp(r / (4 * M))
    if region in ("I", "III"):
        amp = np.sqrt(r / (2 * M) - 1.0) * ex
        X, T = amp * np.cosh(t / (4 * M)), amp * np.sinh(t / (4 * M))
        if region == "III":
            X = -X
    elif region in ("II", "IV"):
        amp = np.sqrt(1.0 - r / (2 * M)) * ex
        T, X = amp * np.cosh(t / (4 * M)), amp * np.sinh(t / (4 * M))
        if region == "IV":
            T = -T
    else:
        raise ValueError(region)
    return T, X


def kruskal_to_r(T, X, M=1.0):
    """r from (r/2M - 1) e^{r/2M} = X^2 - T^2 via the Lambert W function."""
    s = np.asarray(X, dtype=float) ** 2 - np.asarray(T, dtype=float) ** 2
    return 2 * M * (1.0 + np.real(lambertw(s / np.e)))


def kruskal_conformal_factor(r, M=1.0):
    """F(r) in ds^2 = -F (dT^2 - dX^2) + r^2 dOmega^2, F = (32 M^3 / r) e^{-r/2M}."""
    r = np.asarray(r, dtype=float)
    return 32 * M**3 / r * np.exp(-r / (2 * M))


def region_of(T, X):
    """Label ('I', 'II', 'III' or 'IV') of the Kruskal point (T, X) in the maximally extended spacetime."""
    T, X = np.asarray(T, dtype=float), np.asarray(X, dtype=float)
    out = np.full(np.broadcast(T, X).shape, "", dtype=object)
    out[np.abs(X) > np.abs(T)] = "III"
    out[(np.abs(X) > np.abs(T)) & (X > 0)] = "I"
    out[(np.abs(T) > np.abs(X)) & (T > 0)] = "II"
    out[(np.abs(T) > np.abs(X)) & (T < 0)] = "IV"
    return out


def to_penrose(T, X):
    """Map Kruskal (T, X) to compactified Penrose coordinates (T_p, X_p) via arctan of the null coordinates."""
    U = np.asarray(T, dtype=float) - np.asarray(X, dtype=float)
    V = np.asarray(T, dtype=float) + np.asarray(X, dtype=float)
    ut, vt = np.arctan(U), np.arctan(V)
    return vt + ut, vt - ut  # (T_p, X_p)


def penrose_to_kruskal(Tp, Xp):
    """Inverse of to_penrose: Penrose (T_p, X_p) to Kruskal (T, X)."""
    vt, ut = (Tp + Xp) / 2, (Tp - Xp) / 2
    V, U = np.tan(vt), np.tan(ut)
    return (V + U) / 2, (V - U) / 2  # (T, X)
