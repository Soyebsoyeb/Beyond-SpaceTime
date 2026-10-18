"""Symbolic differential geometry (SymPy) used as the independent reference
for every numerical curvature computation in the repository.

Conventions: signature (-,+,+,+); R^a_{bcd} = d_c G^a_{db} - d_d G^a_{cb}
+ G^a_{ce} G^e_{db} - G^a_{de} G^e_{cb}; R_{bd} = R^a_{bad}.
"""
from __future__ import annotations

from itertools import product

import sympy as sp


class Spacetime:
    """Curvature calculator for an arbitrary metric given as a SymPy matrix."""

    def __init__(self, coords, metric, simplify=True):
        """Store coordinates and metric (SymPy matrix), compute the inverse metric."""
        self.x = list(coords)
        self.g = sp.Matrix(metric)
        self.n = len(self.x)
        self._s = (lambda e: sp.simplify(e)) if simplify else (lambda e: e)
        self.ginv = self.g.inv().applyfunc(self._s)
        self._Gamma = None
        self._Riem = None

    # -- connection -------------------------------------------------------
    def christoffel(self):
        """Christoffel symbols Gamma^a_{bc} as a nested list [a][b][c] (cached)."""
        if self._Gamma is None:
            n, x, g, gi = self.n, self.x, self.g, self.ginv
            G = [[[0] * n for _ in range(n)] for _ in range(n)]
            for a, b, c in product(range(n), repeat=3):
                G[a][b][c] = self._s(
                    sum(
                        gi[a, d]
                        * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b]) - sp.diff(g[b, c], x[d]))
                        for d in range(n)
                    )
                    / 2
                )
            self._Gamma = G
        return self._Gamma

    # -- curvature ----------------------------------------------------------
    def riemann_mixed(self):
        """R^a_{bcd}."""
        if self._Riem is None:
            n, x, G = self.n, self.x, self.christoffel()
            R = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
            for a, b, c, d in product(range(n), repeat=4):
                val = sp.diff(G[a][d][b], x[c]) - sp.diff(G[a][c][b], x[d])
                val += sum(G[a][c][e] * G[e][d][b] - G[a][d][e] * G[e][c][b] for e in range(n))
                R[a][b][c][d] = self._s(val)
            self._Riem = R
        return self._Riem

    def riemann_lower(self):
        """R_{abcd} = g_{ae} R^e_{bcd}."""
        n, R = self.n, self.riemann_mixed()
        return [
            [
                [[self._s(sum(self.g[a, e] * R[e][b][c][d] for e in range(n))) for d in range(n)] for c in range(n)]
                for b in range(n)
            ]
            for a in range(n)
        ]

    def ricci(self):
        """Ricci tensor R_{bd} = R^a_{bad}."""
        n, R = self.n, self.riemann_mixed()
        return sp.Matrix(n, n, lambda b, d: self._s(sum(R[a][b][a][d] for a in range(n))))

    def ricci_scalar(self):
        """Ricci scalar R = g^{bd} R_{bd}."""
        Ric = self.ricci()
        return self._s(sum(self.ginv[a, b] * Ric[a, b] for a in range(self.n) for b in range(self.n)))

    def einstein(self):
        """Einstein tensor G_{ab} = R_{ab} - R g_{ab} / 2."""
        return (self.ricci() - self.g * self.ricci_scalar() / 2).applyfunc(self._s)

    def kretschmann(self):
        """K = R_{abcd} R^{abcd}, computed via R^{ab}_{cd} R^{cd}_{ab}."""
        n, R, gi = self.n, self.riemann_mixed(), self.ginv
        # R^{ab}_{cd} = g^{be} R^a_{ecd}
        M = [[[[sum(gi[b, e] * R[a][e][c][d] for e in range(n)) for d in range(n)] for c in range(n)]
              for b in range(n)] for a in range(n)]
        K = sum(M[a][b][c][d] * M[c][d][a][b] for a, b, c, d in product(range(n), repeat=4))
        return self._s(K)

    def vacuum_residual(self):
        """Ricci tensor (zero iff vacuum with Lambda = 0)."""
        return self.ricci()


# ---------------------------------------------------------------------------
# Standard metrics
# ---------------------------------------------------------------------------
t, r, th, ph = sp.symbols("t r theta phi", real=True)
M, l = sp.symbols("M l", positive=True)


def static_spherical(f):
    """ds^2 = -f dt^2 + dr^2/f + r^2 dOmega^2 for a (possibly symbolic) f(r)."""
    return Spacetime([t, r, th, ph], sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th) ** 2))


def schwarzschild():
    """Spacetime object for the Schwarzschild metric with symbols (t, r, theta, phi) and mass M."""
    return static_spherical(1 - 2 * M / r)


def hayward_f(Mv=M, lv=l):
    """Symbolic Hayward metric function with mass Mv and length lv (defaults: the module symbols M, l)."""
    return 1 - 2 * Mv * r**2 / (r**3 + 2 * Mv * lv**2)


def kretschmann_static_formula(f):
    """Closed form K = f''^2 + 4 f'^2/r^2 + 4 (1-f)^2/r^4 for f(r) metrics."""
    return sp.diff(f, r, 2) ** 2 + 4 * sp.diff(f, r) ** 2 / r**2 + 4 * (1 - f) ** 2 / r**4


def einstein_mixed_static(f):
    """G^t_t, G^r_r, G^theta_theta for the f(r) metric."""
    st = static_spherical(f)
    G = st.einstein()
    gi = st.ginv
    mixed = (gi * G).applyfunc(sp.simplify)
    return mixed[0, 0], mixed[1, 1], mixed[2, 2]


def godel_metric(a=sp.Symbol("a", positive=True)):
    """Goedel (1949) in coordinates (t, x, y, z):
    ds^2 = a^2 [ -(dt + e^x dz)^2 + dx^2 + dy^2 + (1/2) e^{2x} dz^2 ]."""
    T, X, Y, Z = sp.symbols("T X Y Z", real=True)
    e = sp.exp(X)
    one_form = sp.Matrix([1, 0, 0, e])
    g = -(one_form * one_form.T)
    g[1, 1] += 1
    g[2, 2] += 1
    g[3, 3] += e**2 / 2
    return Spacetime([T, X, Y, Z], a**2 * g), (T, X, Y, Z, a)
