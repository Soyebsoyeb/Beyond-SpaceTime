"""Quantum branching: system qubit S + environment of N qubits with unitary interaction.

Model A (exactly solvable pure dephasing):  H = sum_k g_k Z_S Z_k, environment qubits start in |+>.
    |psi(t)> = a|0>|E_0(t)> + b|1>|E_1(t)>,  <E_0|E_1> = prod_k cos(2 g_k t),
    rho_S(t)_{01} = a b* prod_k cos(2 g_k t)   (decoherence factor D(t)).
Model B (non-commuting):  H = sum_k g_k (Z_S Z_k + X_S X_k)/... evolves with QuTiP and is checked
    against dense matrix exponentials; used to test positivity/Hermiticity/purity bookkeeping.

Interpretive note: the reduced-state calculation shows suppression of interference between
environment-orthogonal components.  It does not by itself select an interpretation, and it
does not make the branches separate spacetimes.
"""
from __future__ import annotations

import numpy as np
import qutip as qt
from scipy.linalg import expm


def _ops(N):
    """Single-site operators embedded in the (1 + N)-qubit space; index 0 is the system."""
    I2 = qt.qeye(2)

    def emb(op, k):
        lst = [I2] * (N + 1)
        lst[k] = op
        return qt.tensor(lst)

    return {"X": lambda k: emb(qt.sigmax(), k), "Y": lambda k: emb(qt.sigmay(), k),
            "Z": lambda k: emb(qt.sigmaz(), k)}


def initial_state(alpha, beta, N):
    """Product state (alpha|0> + beta|1>) for the system qubit and |+> for each of N environment qubits."""
    sys = alpha * qt.basis(2, 0) + beta * qt.basis(2, 1)
    plus = (qt.basis(2, 0) + qt.basis(2, 1)).unit()
    return qt.tensor([sys.unit()] + [plus] * N)


def dephasing_hamiltonian(g):
    """H = sum_k g_k Z_S Z_k (commuting, exactly solvable pure dephasing)."""
    N = len(g)
    o = _ops(N)
    return sum(gk * o["Z"](0) * o["Z"](k + 1) for k, gk in enumerate(g))


def noncommuting_hamiltonian(g, h_s=0.0):
    """H = sum_k g_k (Z_S Z_k + X_S X_k / 2) + h_s X_S, a generic non-commuting coupling."""
    N = len(g)
    o = _ops(N)
    H = sum(gk * (o["Z"](0) * o["Z"](k + 1) + 0.5 * o["X"](0) * o["X"](k + 1)) for k, gk in enumerate(g))
    return H + h_s * o["X"](0)


def evolve_states(H, psi0, times):
    """Evolve psi0 under H with QuTiP sesolve and return the list of states at the requested times."""
    res = qt.sesolve(H, psi0, times, options={"atol": 1e-12, "rtol": 1e-11})
    return res.states


def reduced_system(state):
    """Reduced density matrix of the system qubit (partial trace over the environment)."""
    return state.ptrace(0)


def decoherence_factor(g, t):
    """Analytic decoherence factor D(t) = prod_k cos(2 g_k t) of the dephasing model."""
    return np.prod(np.cos(2.0 * np.asarray(g) * t))


def l1_coherence(rho):
    """l1-norm of coherence: sum of absolute off-diagonal elements."""
    a = np.abs(rho.full())
    return float(a.sum() - np.trace(a))


def relative_entropy_coherence(rho):
    """Relative-entropy coherence S(diag(rho)) - S(rho)."""
    diag = qt.Qobj(np.diag(np.diag(rho.full())))
    return float(qt.entropy_vn(diag) - qt.entropy_vn(rho))


def check_density_matrix(rho, tol=1e-10):
    """Validate a density matrix: returns trace error, Hermiticity error, minimum eigenvalue and an overall ok flag."""
    M = rho.full()
    herm = np.abs(M - M.conj().T).max()
    ev = np.linalg.eigvalsh((M + M.conj().T) / 2)
    return {"trace_err": abs(np.trace(M) - 1.0), "hermiticity_err": float(herm),
            "min_eig": float(ev.min()), "ok": bool(abs(np.trace(M) - 1) < tol and herm < tol and ev.min() > -tol)}


def dense_evolve(H, psi0, t):
    """Reference evolution exp(-iHt) psi0 by dense matrix exponential (independent of QuTiP)."""
    U = expm(-1j * H.full() * t)
    return U @ psi0.full()


def branch_weights(alpha, beta):
    """Born weights (|alpha|^2, |beta|^2) normalised to sum to one."""
    n = abs(alpha) ** 2 + abs(beta) ** 2
    return abs(alpha) ** 2 / n, abs(beta) ** 2 / n


def environment_overlap_numeric(g, t):
    """|<E_0(t)|E_1(t)>| from explicit conditional environment states (independent of the analytic product)."""
    N = len(g)
    plus = (qt.basis(2, 0) + qt.basis(2, 1)).unit()
    ov = 1.0 + 0j
    for gk in g:
        e0 = (-1j * gk * t * qt.sigmaz()).expm() * plus   # system |0>: Z_S = +1
        e1 = (+1j * gk * t * qt.sigmaz()).expm() * plus   # system |1>: Z_S = -1
        ov *= e0.dag() * e1
    return complex(ov)
