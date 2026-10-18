import numpy as np
import pytest
import sympy as sp

from beyond_spacetime import causal as CA
from beyond_spacetime import regular as H
from beyond_spacetime import symbolic as S


def test_no_causal_signal_from_region_I_to_III():
    bad, total = CA.count_signalling_pairs(80, seed=4)
    assert total == 6400 and bad == 0


def test_causal_future_necessary_condition_sanity():
    # a point inside I does have the interior (II) in its future
    p = (0.0, 1.5)
    q = (2.0, 0.5)
    assert CA.may_be_in_causal_future(p, q)
    # a point in the left exterior is not in the future
    assert not CA.may_be_in_causal_future(p, (0.2, -1.5))


def test_ellis_wormhole_radial_null_crossing():
    assert CA.ellis_crossing_time(-5.0, 5.0) == pytest.approx(10.0)


def test_causal_table_has_expected_entries():
    t = {row["geometry"]: row["I_to_III_causal"] for row in CA.causal_relation_table()}
    assert t["Eternal Schwarzschild (Kruskal extension)"] is False
    assert t["Ellis-Bronnikov wormhole"] is True


# ---------------- Hayward ----------------

def test_hayward_kretschmann_matches_full_riemann_computation():
    st = S.static_spherical(S.hayward_f())
    K_full = sp.lambdify((S.r, S.M, S.l), st.kretschmann(), "numpy")
    for r in (0.3, 1.0, 3.0):
        assert H.kretschmann(r, 1.0, 0.5) == pytest.approx(K_full(r, 1.0, 0.5), rel=1e-10)


def test_hayward_core_is_de_sitter_and_finite():
    assert H.kretschmann(1e-4, 1.0, 0.5) == pytest.approx(H.kretschmann_core_limit(0.5), rel=1e-6)
    r = np.logspace(-4, 2, 200)
    assert np.all(np.isfinite(H.kretschmann(r, 1.0, 0.5)))


def test_hayward_classical_limit():
    assert H.kretschmann(5.0, 1.0, 1e-3) == pytest.approx(48 / 5.0**6, rel=1e-6)
    assert H.f(5.0, 1.0, 1e-6) == pytest.approx(1 - 2 / 5.0, rel=1e-9)
    assert H.hawking_temperature(1.0, 1e-6) == pytest.approx(1 / (8 * np.pi), rel=1e-6)


def test_hayward_horizon_structure_and_extremality():
    assert len(H.horizons(1.0, 0.5)) == 2
    assert len(H.horizons(1.0, 1.0)) == 0
    lc = H.l_critical_over_M()
    hz = H.horizons(1.0, lc * (1 - 1e-9))
    assert len(hz) == 2 and hz[0] == pytest.approx(hz[1], rel=1e-3)
    assert H.hawking_temperature(1.0, lc * (1 - 1e-9)) < 1e-3


def test_effective_stress_energy_matches_einstein_tensor():
    f = S.hayward_f()
    Gtt, Grr, Gthth = S.einstein_mixed_static(f)
    rho = sp.lambdify((S.r, S.M, S.l), -Gtt / (8 * sp.pi), "numpy")
    pr = sp.lambdify((S.r, S.M, S.l), Grr / (8 * sp.pi), "numpy")
    pt = sp.lambdify((S.r, S.M, S.l), Gthth / (8 * sp.pi), "numpy")
    r = np.array([0.3, 1.0, 3.0])
    a, b, c = H.stress_energy(r, 1.0, 0.5)
    assert np.allclose(a, rho(r, 1.0, 0.5)) and np.allclose(b, pr(r, 1.0, 0.5)) and np.allclose(c, pt(r, 1.0, 0.5))


def test_hayward_energy_condition_pattern():
    r = np.linspace(0.01, 20, 4000)
    ec = H.energy_conditions(r, 1.0, 0.5)
    assert ec["NEC"].all() and ec["WEC"].all()
    assert not ec["SEC"][0]                      # strong energy condition violated in the core
    assert not ec["DEC"][-1]                     # dominant energy condition violated at large r


def test_hayward_radial_geodesic_regular_through_core_for_E_above_1():
    assert np.isfinite(H.tau_to_center(1.2, 1.0, 0.5))
    tau, rt = H.tau_outer_horizon_to_turning_point(0.5, 1.0, 0.5)
    assert tau > 0 and 0 < rt < H.horizons(1.0, 0.5)[0]


def test_hayward_energy_condition_thresholds_are_analytic():
    """SEC fails for r < (M l^2)^{1/3}; DEC fails for r > (4 M l^2)^{1/3} (derived from rho = 3Ma/(4 pi (r^3+a)^2))."""
    M, l = 1.0, 0.5
    r_sec, r_dec = (M * l**2) ** (1 / 3), (4 * M * l**2) ** (1 / 3)
    r = np.linspace(0.01, 8, 200001)
    ec = H.energy_conditions(r, M, l)
    assert r[~ec["SEC"]].max() == pytest.approx(r_sec, abs=1e-3)
    assert r[~ec["DEC"]].min() == pytest.approx(r_dec, abs=1e-3)
