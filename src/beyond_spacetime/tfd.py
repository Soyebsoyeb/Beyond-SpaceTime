"""Finite-dimensional thermofield-double (TFD) toy model.

|TFD(beta)> = Z^{-1/2} sum_n e^{-beta E_n/2} |n>_L |n*>_R  for a real-symmetric
Hamiltonian H (transverse-field Ising chain).  Tracing out R gives the thermal state
rho_L = e^{-beta H}/Z, so S(L) equals the thermal entropy.  This is a quantum-information
toy, not a holographic computation: it illustrates (i) purification of a thermal state,
(ii) the no-signalling theorem for an entangled pair, and (iii) that an explicit L-R
interaction is needed for a signal (the qualitative role of the Gao-Jafferis-Wall coupling).
"""
from __future__ import annotations

import numpy as np
import qutip as qt
from scipy.linalg import expm


def ising_hamiltonian(N=5, h=1.0, J=1.0):
    """H = -J sum Z_i Z_{i+1} - h sum X_i (open chain), as a dense real matrix."""
    X = np.array([[0, 1], [1, 0]], dtype=float)
    Z = np.array([[1, 0], [0, -1]], dtype=float)
    I = np.eye(2)

    def op(o, i):
        mats = [I] * N
        mats[i] = o
        out = mats[0]
        for m in mats[1:]:
            out = np.kron(out, m)
        return out

    H = np.zeros((2**N, 2**N))
    for i in range(N - 1):
        H -= J * op(Z, i) @ op(Z, i + 1)
    for i in range(N):
        H -= h * op(X, i)
    return H, op


def tfd_state(H, beta):
    """Thermofield-double state vector for Hamiltonian H at inverse temperature beta; returns (psi, eigenvalues, eigenvectors, Z)."""
    E, V = np.linalg.eigh(H)
    w = np.exp(-beta * (E - E.min()) / 2.0)
    Z = np.sum(w**2)
    d = len(E)
    psi = np.zeros(d * d)
    for n in range(d):
        psi += w[n] * np.kron(V[:, n], V[:, n])
    return psi / np.sqrt(Z), E, V, Z


def thermal_entropy(E, beta):
    """Shannon (von Neumann) entropy of the Gibbs distribution exp(-beta E) / Z."""
    p = np.exp(-beta * (E - E.min()))
    p /= p.sum()
    p = p[p > 1e-300]
    return float(-(p * np.log(p)).sum())


def entropy_of_reduced_L(psi, d):
    """Von Neumann entropy of the left factor (QuTiP partial trace)."""
    ket = qt.Qobj(psi.reshape(-1, 1), dims=[[d, d], [1, 1]])
    rho_L = ket.ptrace(0)
    return float(qt.entropy_vn(rho_L)), rho_L


def reduced_R_numpy(psi, d):
    """Reduced state of the right factor of a bipartite pure state, by NumPy reshape (independent of QuTiP)."""
    M = psi.reshape(d, d)
    return M.T @ M.conj()


def no_signalling_deviation(psi, d, rng):
    """Max |rho_R(after random unitary on L) - rho_R(before)| : must vanish."""
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    Q, _ = np.linalg.qr(A)
    U = np.kron(Q, np.eye(d))
    psi2 = U @ psi.astype(complex)
    return float(np.abs(reduced_R_numpy(psi2, d) - reduced_R_numpy(psi.astype(complex), d)).max())


def signal_with_coupling(psi, N, g, theta, op):
    """Receiver signal  <X_1^R>  difference between 'sender acts' and 'sender idle', after a
    coupling exp(-i g Z_1^L Z_1^R).  Sender acts with exp(-i theta X_1^L)."""
    d = 2**N
    I = np.eye(d)
    Xl = np.kron(op(np.array([[0, 1], [1, 0]], dtype=float), 0), I)
    Zl = np.kron(op(np.array([[1, 0], [0, -1]], dtype=float), 0), I)
    Zr = np.kron(I, op(np.array([[1, 0], [0, -1]], dtype=float), 0))
    Xr = np.kron(I, op(np.array([[0, 1], [1, 0]], dtype=float), 0))
    Ysend = expm(-1j * theta * Xl)
    Ucoup = expm(-1j * g * Zl @ Zr)
    psi = psi.astype(complex)

    def readout(state):
        s = Ucoup @ state
        return float(np.real(np.vdot(s, Xr @ s)))

    return readout(Ysend @ psi) - readout(psi)
