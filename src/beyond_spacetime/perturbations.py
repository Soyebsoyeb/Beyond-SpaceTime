"""Black-hole perturbations: quasinormal modes (Leaver continued fraction),
time-domain evolution of  psi_tt - psi_{r*r*} + V psi = 0 on a static
background ds^2 = -f dt^2 + dr^2/f + ..., and Prony extraction of ringdown.

Leaver (1985) conventions: 2M = 1, time dependence e^{-i omega t}.
Reported frequencies are returned as M*omega.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import root

# --------------------------------------------------------------------------
# Leaver continued fraction (Schwarzschild)
# --------------------------------------------------------------------------

def _leaver_residual(w, l, s, depth=400, inversion=0):
    """n-th inversion of Leaver's continued fraction (2M = 1):
    beta_n - alpha_{n-1} gamma_n/(beta_{n-1} - ...) = alpha_n gamma_{n+1}/(beta_{n+1} - ...)."""
    rho = -1j * w
    eps = s * s - 1.0

    def alpha(n):
        return n * n + (2 * rho + 2) * n + 2 * rho + 1

    def beta(n):
        return -(2 * n * n + (8 * rho + 2) * n + 8 * rho**2 + 4 * rho + l * (l + 1) - eps)

    def gamma(n):
        return n * n + 4 * rho * n + 4 * rho**2 - eps - 1

    up = 0.0 + 0.0j  # infinite tail: alpha_n gamma_{n+1}/(beta_{n+1} - ...)
    for k in range(depth, inversion, -1):
        up = alpha(k - 1) * gamma(k) / (beta(k) - up)
    lhs = beta(inversion)
    down = 0.0 + 0.0j
    for k in range(1, inversion + 1):
        down = alpha(k - 1) * gamma(k) / (beta(k - 1) - down)
    # lhs = beta_n - alpha_{n-1} gamma_n / (beta_{n-1} - ...) ; 'down' holds that fraction
    return (lhs - down) - up


def leaver_qnm(l=2, s=2, n=0, guess=None, depth=400):
    """Return M*omega for the n-th overtone using the n-th inversion."""
    defaults = {(2, 2, 0): 0.7473 - 0.1779j, (2, 0, 0): 0.9672 - 0.1935j,
                (3, 2, 0): 1.1989 - 0.1854j, (2, 2, 1): 0.6934 - 0.5478j}
    w0 = guess if guess is not None else defaults.get((l, s, n), 0.75 - 0.18j)

    def fun(x):
        z = _leaver_residual(complex(x[0], x[1]), l, s, depth, inversion=n)
        return [z.real, z.imag]

    sol = root(fun, [w0.real, w0.imag], tol=1e-13)
    w = complex(sol.x[0], sol.x[1])
    return w / 2.0  # convert 2M = 1 -> M = 1


# --------------------------------------------------------------------------
# Time-domain evolution
# --------------------------------------------------------------------------

def rw_potential_schwarzschild(r, l=2, s=2, M=1.0):
    """Regge-Wheeler/scalar potential V = f [l(l+1)/r^2 + (1 - s^2) 2M/r^3] for spin weight s."""
    f = 1 - 2 * M / r
    return f * (l * (l + 1) / r**2 + (1 - s * s) * 2 * M / r**3)


def scalar_potential_static(r, f, fp, l=2):
    """Massless scalar field on -f dt^2 + dr^2/f + r^2 dOmega^2: V = f (l(l+1)/r^2 + f'/r)."""
    return f * (l * (l + 1) / r**2 + fp / r)


def tortoise_grid(f_func, r_plus, r_max=1500.0, n=400001):
    """Tabulate r*(r) = int dr/f for r > r_plus using x = ln(r - r_plus) sampling."""
    x = np.linspace(np.log(1e-12), np.log(r_max - r_plus), n)
    r = r_plus + np.exp(x)
    integrand = np.exp(x) / f_func(r)
    rs = cumulative_trapezoid(integrand, x, initial=0.0)
    # fix the additive constant: r*(r_ref) = r_ref + ... is convention-free for evolution
    return r, rs - rs[np.searchsorted(r, 4 * r_plus)]


def evolve(Vfunc_rs, rs_min=-400.0, rs_max=1000.0, dx=0.1, t_end=300.0, r0s=20.0, sigma=3.0,
           r_obs_s=40.0, cfl=0.5):
    """Leapfrog evolution of a Gaussian pulse; returns (t, psi at r*=r_obs_s)."""
    x = np.arange(rs_min, rs_max + dx, dx)
    dt = cfl * dx
    V = Vfunc_rs(x)
    psi_prev = np.exp(-((x - r0s) ** 2) / (2 * sigma**2))
    # right-moving initial data: psi(t=-dt) from d'Alembert translation
    psi_prev_t = np.exp(-((x + dt - r0s) ** 2) / (2 * sigma**2))
    psi = psi_prev.copy()
    psi_old = psi_prev_t
    nsteps = int(t_end / dt)
    iobs = int(round((r_obs_s - rs_min) / dx))
    rec_t, rec = [], []
    c2 = (dt / dx) ** 2
    for k in range(nsteps):
        lap = np.zeros_like(psi)
        lap[1:-1] = psi[2:] - 2 * psi[1:-1] + psi[:-2]
        new = 2 * psi - psi_old + c2 * lap - dt * dt * V * psi
        new[0] = new[-1] = 0.0
        psi_old, psi = psi, new
        if k % 5 == 0:
            rec_t.append((k + 1) * dt)
            rec.append(psi[iobs])
    return np.array(rec_t), np.array(rec)


def prony(t, y, order=4):
    """Matrix-pencil/linear-prediction fit y ~ sum a_k exp(-i w_k t); returns complex w_k."""
    dt = t[1] - t[0]
    N = len(y)
    L = N // 2
    H = np.array([y[i:i + L] for i in range(N - L)])
    U, s, Vh = np.linalg.svd(H, full_matrices=False)
    Ur = U[:, :order]
    z = np.linalg.eigvals(np.linalg.pinv(Ur[:-1]) @ Ur[1:])
    return 1j * np.log(z) / dt


def ringdown_frequency(t, psi, t_start, t_stop, order=6):
    """Dominant slowly-decaying positive-frequency mode in the window."""
    m = (t >= t_start) & (t <= t_stop)
    w = prony(t[m], psi[m], order=order)
    w = w[(w.real > 0) & (w.imag < 0)]
    return w[np.argmin(np.abs(w.imag))] if len(w) else np.nan


def potential_function_rs(f_func, fp_func, r_plus, l=2, r_max=3000.0):
    """Return V(r*) interpolated from the tabulated tortoise map of a static metric."""
    r, rs = tortoise_grid(f_func, r_plus, r_max=r_max)
    V = scalar_potential_static(r, f_func(r), fp_func(r), l)
    from scipy.interpolate import interp1d
    Vi = interp1d(rs, V, bounds_error=False, fill_value=(0.0, 0.0))
    return lambda x: Vi(x)
