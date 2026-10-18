import numpy as np
import pytest
import sympy as sp

from beyond_spacetime import chronology as CH
from beyond_spacetime import multiverse as MV
from beyond_spacetime import symbolic as S
from beyond_spacetime import timereversal as TR
from beyond_spacetime import wormholes as W


# ---------------- time reversal ----------------

def test_classical_reversal_returns_to_initial_state():
    assert TR.leapfrog_reversal_error([1.0, 0.3], [0.0, 0.5], lambda q: -q, 0.01, 800) < 1e-12


def test_wave_equation_reversible_heat_equation_not():
    assert TR.wave_reversal_error() < 1e-12
    assert TR.heat_reversal_error(steps=60) > 1e-3


def test_quantum_time_reversal_spin_half():
    e1, e2 = TR.quantum_reversal_check([0.3, 0.5, 0.7], 1.3, np.array([0.6, 0.8j]))
    assert e1 < 1e-12 and e2 < 1e-12     # Theta U Theta^-1 = exp(+i H' t) with H' = -H


def test_real_hamiltonian_has_time_reversal_symmetry():
    rng = np.random.default_rng(0)
    H = rng.normal(size=(8, 8))
    H = H + H.T
    psi0 = rng.normal(size=8) + 0j
    assert TR.real_hamiltonian_reversal_error(H, psi0, 0.9) < 1e-10


def test_kruskal_time_reversal_is_isometry():
    T = np.array([0.2, 0.5, -0.3])
    X = np.array([0.4, 0.1, 0.8])
    assert TR.kruskal_time_reversal_invariance(T, X) < 1e-12


def test_cat_map_entropy_growth_and_reversal():
    res = TR.cat_experiment(n=8000, t_max=45, seed=1)
    assert res["S"][0] < 0.1 and res["S"][-1] > 0.95 * res["S_max"]
    assert res["err_int"].max() == 0                        # exact reversal on the integer lattice
    assert res["err_float"][10] < 1e-8 and res["err_float"][-1] > 0.1   # float reversal fails
    # error growth rate approximates the Lyapunov exponent in the exponential regime
    t = np.arange(12, 28)
    slope = np.polyfit(t, np.log(res["err_float"][12:28]), 1)[0]
    assert slope == pytest.approx(TR.LYAPUNOV_CAT, rel=0.15)


# ---------------- wormholes ----------------

@pytest.fixture(scope="module")
def b0():
    return sp.Symbol("b0", positive=True)


def test_morris_thorne_equations_ellis(b0):
    b = b0**2 / S.r
    Gtt, Grr, Gthth = W.mt_orthonormal_einstein(b, sp.Integer(0))
    rho, tau, p = W.mt_matter(b, sp.Integer(0))
    assert sp.simplify(Gtt - 8 * sp.pi * rho) == 0
    assert sp.simplify(Grr + 8 * sp.pi * tau) == 0
    assert sp.simplify(Gthth - 8 * sp.pi * p) == 0


def test_morris_thorne_equations_with_redshift_function(b0):
    b = b0**2 / S.r
    Phi = -sp.Rational(1, 3) * b0 / S.r
    Gtt, Grr, Gthth = W.mt_orthonormal_einstein(b, Phi)
    rho, tau, p = W.mt_matter(b, Phi)
    pts = {S.r: 2.3, b0: 1.3}
    assert abs(float((Gtt - 8 * sp.pi * rho).subs(pts))) < 1e-12
    assert abs(float((Grr + 8 * sp.pi * tau).subs(pts))) < 1e-12
    assert abs(float((Gthth - 8 * sp.pi * p).subs(pts))) < 1e-12


def test_throat_nec_violation_and_flare_out(b0):
    b = b0**2 / S.r
    rho, tau, p = W.mt_matter(b, sp.Integer(0))
    val = float((rho - tau).subs({S.r: 1.3, b0: 1.3}))
    assert val < 0
    assert val == pytest.approx(W.nec_radial_at_throat(1.3), rel=1e-12)
    assert float(sp.diff(b, S.r).subs({S.r: 1.3, b0: 1.3})) < 1.0       # flare-out


def test_ellis_tidal_formula_matches_symbolic_riemann():
    expr, (l, b0s, v) = W.ellis_boosted_tidal_symbolic()
    for lv, vv in [(0.0, 0.3), (0.4, 0.6), (1.5, 0.9)]:
        assert float(expr.subs({l: lv, b0s: 1.0, v: vv})) == pytest.approx(W.ellis_tidal_lateral(lv, vv, 1.0), rel=1e-10)


def test_ellis_embedding_throat_and_area():
    assert W.ellis_areal_radius(0.0, 1.3) == pytest.approx(1.3)
    assert W.ellis_embedding_z(0.0, 1.3) == pytest.approx(0.0)


# ---------------- chronology ----------------

def test_ctc_condition_and_horizon():
    D, beta = 1.0, 0.8
    Th = CH.chronology_horizon_time(D, beta)
    assert CH.clock_lag(Th, beta) == pytest.approx(D)
    assert not CH.ctc_exists(CH.clock_lag(0.99 * Th, beta), D)
    assert CH.ctc_exists(CH.clock_lag(1.01 * Th, beta), D)


def test_ctc_map_monotone_in_speed():
    betas = np.linspace(0.1, 0.99, 10)
    m = CH.ctc_map([0.5], betas, 5.0)[:, 0]
    assert (np.diff(m.astype(int)) >= 0).all() and m[-1] and not m[0]


def test_godel_ctc_threshold():
    rc = CH.R_CRIT_GODEL
    assert rc == pytest.approx(np.log(1 + np.sqrt(2)))
    assert not CH.godel_circle_is_closed_timelike(0.9 * rc)
    assert CH.godel_circle_is_closed_timelike(1.1 * rc)


def test_godel_ricci_is_dust_form():
    g, (T, X, Y, Z, a) = S.godel_metric()
    Ric = g.ricci()
    ex = sp.exp(X)
    expected = sp.Matrix([[1, 0, 0, ex], [0, 0, 0, 0], [0, 0, 0, 0], [ex, 0, 0, ex**2]])
    assert sp.simplify(Ric - expected) == sp.zeros(4, 4)


# ---------------- multiverse table ----------------

def test_multiverse_statuses_are_classified():
    allowed = {MV.EQUATIONS, MV.INTERPRETIVE, MV.EXTRA}
    assert len(MV.CONCEPTS) == 4 and all(c["status"] in allowed for c in MV.CONCEPTS)
    status = {c["key"]: c["status"] for c in MV.CONCEPTS}
    assert status["kruskal"] == MV.EQUATIONS and status["mwi"] == MV.INTERPRETIVE
