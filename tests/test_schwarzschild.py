import numpy as np
import pytest
import sympy as sp

from beyond_spacetime import coordinates as C
from beyond_spacetime import embedding as E
from beyond_spacetime import geodesics as G
from beyond_spacetime import symbolic as S
from beyond_spacetime import tides as T


@pytest.fixture(scope="module")
def sch():
    return S.schwarzschild()


def test_vacuum_einstein_equations(sch):
    assert sch.ricci() == sp.zeros(4, 4)
    assert sch.ricci_scalar() == 0


def test_kretschmann_scalar(sch):
    assert sp.simplify(sch.kretschmann() - 48 * S.M**2 / S.r**6) == 0


def test_generic_static_kretschmann_formula():
    f = sp.Function("f")(S.r)
    st = S.static_spherical(f)
    assert sp.simplify(st.kretschmann() - S.kretschmann_static_formula(f)) == 0


def test_christoffel_known_components(sch):
    G_ = sch.christoffel()
    f = 1 - 2 * S.M / S.r
    assert sp.simplify(G_[0][0][1] - sp.diff(f, S.r) / (2 * f)) == 0       # Gamma^t_{tr}
    assert sp.simplify(G_[1][1][1] + sp.diff(f, S.r) / (2 * f)) == 0       # Gamma^r_{rr}
    assert sp.simplify(G_[1][2][2] + S.r * f) == 0                         # Gamma^r_{theta theta}


def test_kruskal_pullback_metric():
    """Pull back -F(dT^2 - dX^2) to (t, r) in region I and recover the Schwarzschild metric."""
    t, r, M = sp.symbols("t r M", positive=True)
    amp = sp.sqrt(r / (2 * M) - 1) * sp.exp(r / (4 * M))
    X, Tk = amp * sp.cosh(t / (4 * M)), amp * sp.sinh(t / (4 * M))
    F = 32 * M**3 / r * sp.exp(-r / (2 * M))
    gtt = -F * (sp.diff(Tk, t) ** 2 - sp.diff(X, t) ** 2)
    grr = -F * (sp.diff(Tk, r) ** 2 - sp.diff(X, r) ** 2)
    gtr = -F * (sp.diff(Tk, t) * sp.diff(Tk, r) - sp.diff(X, t) * sp.diff(X, r))
    pts = {t: 1.7, r: 3.9, M: 1.0}
    f = 1 - 2 * M / r
    assert abs(float(gtt.subs(pts)) - float((-f).subs(pts))) < 1e-12
    assert abs(float(grr.subs(pts)) - float((1 / f).subs(pts))) < 1e-12
    assert abs(float(gtr.subs(pts))) < 1e-12


def test_ef_metric_regular_at_horizon():
    """Ingoing EF: ds^2 = -f dv^2 + 2 dv dr has nonvanishing determinant at r = 2M."""
    f = sp.Symbol("f")
    g = sp.Matrix([[-f, 1], [1, 0]])
    assert g.det() == -1                                   # independent of f, regular at f = 0
    # and the transformation t = v - r*: dt = dv - dr/f reproduces Schwarzschild
    v, r = sp.symbols("v r")
    dv, dr = sp.symbols("dv dr")
    dt = dv - dr / f
    expr = sp.expand(-f * dt**2 + dr**2 / f)
    assert sp.simplify(expr - (-f * dv**2 + 2 * dv * dr)) == 0


@pytest.mark.parametrize("region", ["I", "II", "III", "IV"])
def test_kruskal_roundtrip_r(region):
    r = 3.3 if region in ("I", "III") else 1.2
    Tk, Xk = C.schwarzschild_to_kruskal(0.9, r, 1.0, region)
    assert abs(C.kruskal_to_r(Tk, Xk, 1.0) - r) < 1e-12
    assert C.region_of(Tk, Xk) == region


def test_penrose_roundtrip_and_singularity_line():
    Tk, Xk = np.array([0.3, -0.7, 1.1]), np.array([0.2, 0.5, 0.4])
    Tp, Xp = C.to_penrose(Tk, Xk)
    T2, X2 = C.penrose_to_kruskal(Tp, Xp)
    assert np.allclose(T2, Tk) and np.allclose(X2, Xk)
    # singularity T^2 - X^2 = 1 maps to the straight line T_p = pi/2
    X = np.linspace(-0.9, 0.9, 9)
    Tsing = np.sqrt(1 + X**2)
    Tp, _ = C.to_penrose(Tsing, X)
    assert np.allclose(Tp, np.pi / 2)


def test_flamm_embedding_induces_spatial_metric():
    r, M = sp.symbols("r M", positive=True)
    z = 2 * sp.sqrt(2 * M * (r - 2 * M))
    assert sp.simplify(1 + sp.diff(z, r) ** 2 - 1 / (1 - 2 * M / r)) == 0
    rr = np.linspace(2.5, 9, 7)
    assert np.allclose(E.flamm_dzdr(rr) ** 2, 2 / (rr - 2))


def test_proper_distance_closed_form_matches_quadrature():
    from scipy.integrate import quad
    r1 = 7.0
    val, _ = quad(lambda r: 1 / np.sqrt(1 - 2 / r), 2.0, r1)
    assert abs(val - E.proper_radial_distance(r1)) < 1e-9


def test_throat_area_and_slicing_dependence():
    assert E.throat_radius_kruskal(0.0) == pytest.approx(2.0)
    radii = [E.throat_radius_kruskal(T0) for T0 in (0.0, 0.4, 0.8, 0.95)]
    assert all(a > b for a, b in zip(radii, radii[1:]))      # throat shrinks on later slices
    assert radii[-1] < 2.0
    r, l = E.painleve_gullstrand_profile()
    assert r.min() < 2.0 and np.all(np.diff(r) > 0)          # PG slice has no throat


def test_horizon_proper_time_analytic_E1():
    tr = G.integrate_infall(8.0, 1.0, r_end=1e-4)
    assert tr["tau_horizon"] == pytest.approx(G.tau_E1(8.0, 2.0), rel=1e-9)
    assert tr["tau_end"] == pytest.approx(G.tau_E1(8.0, 1e-4), rel=1e-8)


def test_rest_infall_cycloid():
    r0 = 10.0
    tr = G.integrate_infall(r0, G.energy_from_rest(r0), r_end=1e-4)
    assert tr["tau_horizon"] == pytest.approx(G.tau_rest_at(r0, 2.0), rel=1e-9)
    assert tr["tau_end"] == pytest.approx(G.tau_rest_to_singularity(r0), rel=1e-5)


def test_max_proper_time_inside_horizon_is_pi_M():
    assert G.tau_inside(0.0) == pytest.approx(np.pi, rel=1e-9)
    vals = [G.tau_inside(E) for E in (0.0, 0.3, 0.7, 1.0, 2.0, 5.0)]
    assert all(a > b for a, b in zip(vals, vals[1:]))
    assert G.tau_inside(1.0) == pytest.approx(4.0 / 3.0, rel=1e-9)


def test_coordinate_time_diverges_but_ef_v_is_finite():
    tr = G.integrate_infall(8.0, 1.0, r_end=1.0)
    outside = tr["r"] > 2.0001
    tcoord = G.schwarzschild_time(tr["r"][outside], tr["v"][outside])
    assert tcoord[-1] > tcoord[len(tcoord) // 2]
    assert np.all(np.isfinite(tr["v"]))


def test_tidal_tensor_from_full_riemann_matches_closed_form_independent_of_E():
    for r, Ev in [(3.0, 1.3), (5.0, 0.8), (0.7, 0.4), (1.5, 2.0)]:
        Tt = T.frame_tidal_tensor(r, Ev, 1.0)
        assert np.allclose(Tt, T.tidal_closed_form(r, 1.0), atol=1e-10)


def test_tidal_trace_free_in_vacuum():
    assert np.trace(T.tidal_closed_form(2.7, 1.0)) == pytest.approx(0.0, abs=1e-14)


def test_deviation_matches_neighbouring_geodesics():
    r0, E, dE = 10.0, 0.9, 1e-6
    tau, sep = T.neighbouring_geodesic_separation(r0, E, dE, n=300, r_end=0.2)
    # initial relative velocity from the first nonzero separation
    xidot0 = sep[1] / tau[1]
    dev = T.integrate_deviation(r0, E, xi0=(0.0, 0.0), xidot0=(xidot0, 0.0), r_end=0.2, n=300)
    from scipy.interpolate import interp1d
    xi = interp1d(dev["tau"], np.abs(dev["xi_r"]))
    k = np.arange(1, len(tau), 7)
    k = k[tau[k] < dev["tau_end"] * 0.95]
    ref = np.array([float(xi(tau[i])) for i in k])
    rel = np.abs(ref - sep[k]) / ref
    assert rel.max() < 5e-3


def test_transverse_deviation_equals_r_times_dtheta():
    d = T.integrate_deviation(10.0, 0.9, xi0=(0.0, 10.0 * 1e-3), xidot0=(0.0, 0.0), r_end=0.5)
    # exact: xi_perp = r(tau) dtheta (nonzero initial rdot component: xi_dot = rdot dtheta)
    traj = G.integrate_infall(10.0, 0.9, r_end=0.5)
    # reintegrate with the correct initial velocity
    d = T.integrate_deviation(10.0, 0.9, xi0=(0.0, 10.0 * 1e-3), xidot0=(0.0, traj["rdot"][0] * 1e-3), r_end=0.5)
    from scipy.interpolate import interp1d
    r_of_tau = interp1d(traj["tau"], traj["r"])
    m = d["tau"] < min(traj["tau_end"], d["tau_end"]) * 0.99
    assert np.allclose(d["xi_perp"][m], 1e-3 * r_of_tau(d["tau"][m]), rtol=1e-6)


def test_spaghettification_radius_inside_horizon_for_large_mass():
    r_big, ratio_big = T.spaghettification_radius(1e9)
    r_small, ratio_small = T.spaghettification_radius(10.0)
    assert ratio_big < 1.0 < ratio_small


def test_si_tidal_acceleration_dimensional_check():
    """a(r_s) = c^6 xi / (4 G^2 M^2) restored from geometrized units."""
    from beyond_spacetime.conventions import C_SI, G_SI, M_SUN_KG, geometric_mass_m, tidal_accel_si
    M_kg, xi = 10 * M_SUN_KG, 2.0
    Mg = geometric_mass_m(M_kg)
    expected = C_SI**6 * xi / (4 * G_SI**2 * M_kg**2)
    assert tidal_accel_si(Mg, 2 * Mg, xi) == pytest.approx(expected, rel=1e-12)
