"""Time reversal across four settings, and microscopic reversibility versus
thermodynamic irreversibility.

Classical:   (q, p) -> (q, -p); orbital angular momentum flips sign.
Wave eq.:    u_tt = c^2 u_xx is invariant under t -> -t (u_t -> -u_t); the heat equation is not.
Quantum:     Theta antiunitary; Theta e^{-iHt} Theta^{-1} = e^{+i (Theta H Theta^{-1}) t};
             a symmetry iff Theta H Theta^{-1} = H.  Spin-1/2: Theta = i sigma_y K flips B.
Gravity:     t -> -t is an isometry of Schwarzschild (static); in Kruskal it is T -> -T,
             exchanging black-hole region II with white-hole region IV.
Irreversibility: arises from special (low-entropy) initial conditions plus coarse graining,
             demonstrated with Arnold's cat map (exact reversal on an integer lattice, destroyed
             in floating point by exponential error growth with Lyapunov exponent ln((3+sqrt5)/2)).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from .coordinates import kruskal_conformal_factor, kruskal_to_r

# ---------------- classical ----------------

def leapfrog_reversal_error(q0, p0, force, dt, nsteps):
    """Integrate with velocity Verlet, flip p, integrate again; return |(q, -p) - (q0, p0)|."""
    q, p = np.array(q0, float), np.array(p0, float)
    for _ in range(nsteps):
        p = p + 0.5 * dt * force(q)
        q = q + dt * p
        p = p + 0.5 * dt * force(q)
    p = -p
    for _ in range(nsteps):
        p = p + 0.5 * dt * force(q)
        q = q + dt * p
        p = p + 0.5 * dt * force(q)
    return float(np.linalg.norm(np.concatenate([q - q0, -p - p0])))


# ---------------- waves vs diffusion ----------------

def wave_reversal_error(n=256, steps=400, c=0.5, seed=0):
    """Leapfrog wave equation on a periodic grid: forward, swap time levels, forward again."""
    x = np.linspace(0, 1, n, endpoint=False)
    u0 = np.exp(-((x - 0.5) ** 2) / 0.002)
    u_prev, u = u0.copy(), u0.copy()  # zero initial velocity
    def step(uo, uc):
        return 2 * uc - uo + c**2 * (np.roll(uc, -1) - 2 * uc + np.roll(uc, 1))
    a, b = u_prev, u
    for _ in range(steps):
        a, b = b, step(a, b)
    a, b = b, a            # time reversal: swap the two time levels
    for _ in range(steps):
        a, b = b, step(a, b)
    return float(np.abs(a - u0).max())


def heat_reversal_error(n=128, steps=40, nu=0.2):
    """Explicit diffusion forward, then the (formally reversed) anti-diffusion: errors explode."""
    x = np.linspace(0, 1, n, endpoint=False)
    u0 = np.exp(-((x - 0.5) ** 2) / 0.002)
    u = u0.copy()
    for _ in range(steps):
        u = u + nu * (np.roll(u, -1) - 2 * u + np.roll(u, 1))
    for _ in range(steps):
        u = u - nu * (np.roll(u, -1) - 2 * u + np.roll(u, 1))
    return float(np.abs(u - u0).max())


# ---------------- quantum ----------------

SIGMA = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]


def theta_spin_half(psi):
    """Theta = i sigma_y K."""
    return 1j * SIGMA[1] @ psi.conj()


def quantum_reversal_check(B, t, psi0):
    """Return (|Theta U(t) Theta^-1 psi - e^{+iH't} psi|, |U(-t)U(t) psi - psi|) for H = B.sigma.
    H' = Theta H Theta^-1 = -B.sigma, i.e. reversing the field direction."""
    H = sum(b * s for b, s in zip(B, SIGMA))
    Hp = -H
    U = expm(-1j * H * t)
    # Theta U Theta^-1 psi = Theta (U (Theta^-1 psi)); Theta^-1 = -Theta on spin-1/2
    psi = psi0
    lhs = theta_spin_half(U @ (-theta_spin_half(psi)))
    rhs = expm(+1j * Hp * t) @ psi
    err1 = float(np.linalg.norm(lhs - rhs))
    err2 = float(np.linalg.norm(expm(1j * H * t) @ (U @ psi) - psi))
    return err1, err2


def real_hamiltonian_reversal_error(H, psi0, t):
    """For real symmetric H, Theta = K is a symmetry: psi -> K U(t) psi evolved again returns K psi0."""
    U = expm(-1j * H * t)
    psi_t = U @ psi0
    back = U @ psi_t.conj()
    return float(np.linalg.norm(back - psi0.conj()))


# ---------------- gravity ----------------

def kruskal_time_reversal_invariance(T, X, M=1.0):
    """Metric function F(r(T, X)) and r are invariant under T -> -T (returns max difference)."""
    r1, r2 = kruskal_to_r(T, X, M), kruskal_to_r(-np.asarray(T), X, M)
    F1, F2 = kruskal_conformal_factor(r1, M), kruskal_conformal_factor(r2, M)
    return float(np.max(np.abs(r1 - r2)) + np.max(np.abs(F1 - F2)))


# ---------------- irreversibility ----------------

CAT = np.array([[2, 1], [1, 1]], dtype=np.int64)
CAT_INV = np.array([[1, -1], [-1, 2]], dtype=np.int64)
LYAPUNOV_CAT = float(np.log((3 + np.sqrt(5)) / 2))


def coarse_grained_entropy(x, y, cells=8):
    """Boltzmann coarse-grained entropy -sum p ln p over a cells x cells partition of the unit square."""
    H, _, _ = np.histogram2d(np.mod(x, 1.0), np.mod(y, 1.0), bins=cells, range=[[0, 1], [0, 1]])
    p = H.ravel() / H.sum()
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())


def cat_map_float(x, y, steps, inverse=False):
    """Iterate Arnold's cat map (or its inverse) in float64 on the unit torus."""
    A = CAT_INV if inverse else CAT
    for _ in range(steps):
        x, y = np.mod(A[0, 0] * x + A[0, 1] * y, 1.0), np.mod(A[1, 0] * x + A[1, 1] * y, 1.0)
    return x, y


def cat_map_int(x, y, steps, L, inverse=False):
    """Iterate Arnold's cat map (or its inverse) exactly in integer arithmetic on an L x L lattice."""
    A = CAT_INV if inverse else CAT
    for _ in range(steps):
        x, y = (A[0, 0] * x + A[0, 1] * y) % L, (A[1, 0] * x + A[1, 1] * y) % L
    return x, y


def cat_experiment(n=20000, t_max=30, L=2**40, seed=0, cells=8):
    """Entropy growth S(t) and Loschmidt-reversal errors (float vs exact integer lattice)."""
    rng = np.random.default_rng(seed)
    xi = rng.integers(0, L // 16, n, dtype=np.int64)  # low-entropy start: small corner cell
    yi = rng.integers(0, L // 16, n, dtype=np.int64)
    # Floating-point branch: full-mantissa random start in the same corner cell (dyadic lattice
    # points would be mapped without rounding error and hide the effect).
    x0, y0 = rng.random(n) / 16.0, rng.random(n) / 16.0
    S, err_float, err_int = [], [], []
    for t in range(t_max + 1):
        xt, yt = cat_map_float(x0.copy(), y0.copy(), t)
        S.append(coarse_grained_entropy(xt, yt, cells))
        xb, yb = cat_map_float(xt, yt, t, inverse=True)
        d = lambda a, b: np.minimum(np.abs(a - b), 1 - np.abs(a - b))
        err_float.append(float(max(d(xb, x0).max(), d(yb, y0).max())))
        xe, ye = cat_map_int(xi.copy(), yi.copy(), t, L)
        xe, ye = cat_map_int(xe, ye, t, L, inverse=True)
        err_int.append(int(max(np.abs(xe - xi).max(), np.abs(ye - yi).max())))
    return {"t": np.arange(t_max + 1), "S": np.array(S), "err_float": np.array(err_float),
            "err_int": np.array(err_int), "S_max": float(np.log(cells * cells))}
