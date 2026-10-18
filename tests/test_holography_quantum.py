import numpy as np
import pytest
import qutip as qt

from beyond_spacetime import holography as Hh
from beyond_spacetime import perturbations as P
from beyond_spacetime import quantum as Q
from beyond_spacetime import tfd


# ---------------- holography ----------------

@pytest.mark.parametrize("rp,dphi", [(0.8, 0.7), (1.7, 1.5), (3.0, 3.0)])
def test_btz_geodesic_numeric_equals_closed_form(rp, dphi):
    rc = 1e4
    J = Hh.geodesic_J_from_angle(dphi, rp)
    assert Hh.angle_from_J_numeric(J, rp) == pytest.approx(dphi, rel=1e-10)
    assert Hh.length_numeric(J, rp, rc) == pytest.approx(Hh.geodesic_length_closed(dphi, rp, rc), abs=5e-7)


def test_rt_equals_calabrese_cardy():
    for rp in (0.5, 2.0):
        for dphi in (0.3, 2.0):
            assert Hh.entanglement_entropy_RT(dphi, rp, 1e5) == pytest.approx(
                Hh.entanglement_entropy_CFT(dphi, rp, 1e5), rel=1e-12)


def test_vacuum_limit():
    assert Hh.geodesic_length_closed(1.0, 1e-6, 1e4) == pytest.approx(Hh.vacuum_limit_length(1.0, 1e4), abs=1e-9)


def test_rt_homology_phase_transition_and_araki_lieb_saturation():
    rp, rc = 1.5, 1e5
    ds = Hh.transition_angle(rp)
    assert np.pi < ds < 2 * np.pi
    # below the transition the connected geodesic dominates; above it S_A - S_{A^c} = S_BH
    d_hi = 0.5 * (ds + 2 * np.pi)
    diff = Hh.rt_entropy(d_hi, rp, rc) - Hh.rt_entropy(2 * np.pi - d_hi, rp, rc)
    assert diff == pytest.approx(Hh.bh_entropy(rp), rel=1e-10)
    d_lo = 1.0
    assert Hh.rt_entropy(d_lo, rp, rc) == pytest.approx(Hh.geodesic_length_closed(d_lo, rp, rc) / 4.0)
    # continuity at the transition
    e = 1e-7
    assert Hh.rt_entropy(ds - e, rp, rc) == pytest.approx(Hh.rt_entropy(ds + e, rp, rc), abs=1e-5)


def test_embedding_points_lie_on_AdS3_hyperboloid():
    P_ = Hh.btz_embedding(5.0, 0.3, 0.4, 1.7, +1)
    n = -P_[0] ** 2 - P_[1] ** 2 + P_[2] ** 2 + P_[3] ** 2
    assert n == pytest.approx(-1.0, abs=1e-12)


def test_embedding_same_side_distance_equals_geodesic():
    rc, rp, dphi = 1e5, 1.3, 0.9
    P_ = Hh.btz_embedding(rc, 0.0, 0.0, rp, +1)
    Q_ = Hh.btz_embedding(rc, 0.0, dphi, rp, +1)
    assert Hh.embedding_distance(P_, Q_) == pytest.approx(Hh.geodesic_length_closed(dphi, rp, rc), abs=1e-7)


def test_two_sided_wormhole_length_growth():
    rc, rp = 1e4, 1.7
    for t in (0.0, 0.5, 2.0):
        assert Hh.two_sided_length(t, rp, rc) == pytest.approx(Hh.two_sided_length_closed(t, rp, rc), abs=1e-6)
    slope = (Hh.two_sided_length_closed(6.0, rp, rc) - Hh.two_sided_length_closed(5.0, rp, rc))
    assert slope == pytest.approx(2 * rp, rel=1e-3)


def test_bh_entropy_equals_cardy_entropy():
    assert Hh.bh_entropy(1.7) == pytest.approx(Hh.cardy_entropy(1.7), rel=1e-12)


# ---------------- TFD toy ----------------

def test_tfd_reduced_state_is_thermal_and_entropy_matches():
    N = 4
    Hm, op = tfd.ising_hamiltonian(N)
    d = 2**N
    for beta in (0.3, 1.5):
        psi, E, V, Z = tfd.tfd_state(Hm, beta)
        assert np.linalg.norm(psi) == pytest.approx(1.0, abs=1e-12)
        S, rho = tfd.entropy_of_reduced_L(psi, d)
        assert S == pytest.approx(tfd.thermal_entropy(E, beta), abs=1e-10)
        th = np.linalg.matrix_power(np.eye(d), 1)  # placeholder to keep numpy import used
        rho_th = qt.Qobj(__import__("scipy.linalg", fromlist=["expm"]).expm(-beta * Hm))
        rho_th = rho_th / rho_th.tr()
        assert (rho - rho_th).norm() < 1e-9


def test_tfd_annihilated_by_HL_minus_HR():
    N = 4
    Hm, _ = tfd.ising_hamiltonian(N)
    d = 2**N
    psi, *_ = tfd.tfd_state(Hm, 0.8)
    diff = np.kron(Hm, np.eye(d)) @ psi - np.kron(np.eye(d), Hm) @ psi
    assert np.linalg.norm(diff) < 1e-10


def test_no_signalling_without_coupling_and_signal_with_coupling():
    N = 3
    Hm, op = tfd.ising_hamiltonian(N)
    d = 2**N
    psi, *_ = tfd.tfd_state(Hm, 0.6)
    rng = np.random.default_rng(5)
    assert tfd.no_signalling_deviation(psi, d, rng) < 1e-12
    assert abs(tfd.signal_with_coupling(psi, N, 0.0, 0.7, op)) < 1e-12
    assert abs(tfd.signal_with_coupling(psi, N, 0.3, 0.7, op)) > 1e-3


# ---------------- quantum branching ----------------

def test_dephasing_matches_analytic_decoherence_factor():
    rng = np.random.default_rng(3)
    N = 5
    g = rng.uniform(0.5, 1.5, N)
    a, b = np.sqrt(0.3), np.sqrt(0.7)
    psi0 = Q.initial_state(a, b, N)
    ts = np.linspace(0, 1.5, 4)
    states = Q.evolve_states(Q.dephasing_hamiltonian(g), psi0, ts)
    for t, s in zip(ts, states):
        rho = Q.reduced_system(s)
        assert abs(rho.full()[0, 1]) == pytest.approx(a * b * abs(Q.decoherence_factor(g, t)), abs=1e-8)
        assert abs(Q.environment_overlap_numeric(g, t)) == pytest.approx(abs(Q.decoherence_factor(g, t)), abs=1e-10)
        # Born weights unchanged by pure dephasing
        assert rho.full()[0, 0].real == pytest.approx(0.3, abs=1e-9)
        assert Q.check_density_matrix(rho)["ok"]


def test_noncommuting_model_qutip_equals_dense_expm_and_density_matrix_valid():
    rng = np.random.default_rng(8)
    g = rng.uniform(0.4, 1.2, 4)
    psi0 = Q.initial_state(0.6, 0.8, 4)
    H = Q.noncommuting_hamiltonian(g, 0.2)
    t = 1.1
    s = Q.evolve_states(H, psi0, [0, t])[-1]
    assert np.abs(s.full()[:, 0] - Q.dense_evolve(H, psi0, t)[:, 0]).max() < 1e-7
    rho = Q.reduced_system(s)
    assert Q.check_density_matrix(rho)["ok"]
    assert qt.entropy_vn(rho) >= -1e-12


def test_entropy_bounded_by_ln2_and_total_state_stays_pure():
    g = np.full(4, 0.9)
    psi0 = Q.initial_state(np.sqrt(0.5), np.sqrt(0.5), 4)
    s = Q.evolve_states(Q.noncommuting_hamiltonian(g), psi0, [0, 1.3])[-1]
    assert qt.entropy_vn(Q.reduced_system(s)) <= np.log(2) + 1e-10
    assert abs((s * s.dag()).tr() - 1) < 1e-10 and abs(s.norm() - 1) < 1e-10


def test_coherence_measures_consistent():
    rho = qt.Qobj(np.array([[0.5, 0.2], [0.2, 0.5]]))
    assert Q.l1_coherence(rho) == pytest.approx(0.4)
    assert Q.relative_entropy_coherence(rho) > 0
    assert Q.l1_coherence(qt.Qobj(np.diag([0.5, 0.5]))) == pytest.approx(0.0)


# ---------------- perturbations ----------------

def test_leaver_qnm_known_values():
    w = P.leaver_qnm(2, 2, 0)
    assert w.real == pytest.approx(0.37367, abs=2e-5) and w.imag == pytest.approx(-0.08896, abs=2e-5)
    w = P.leaver_qnm(2, 0, 0)
    assert w.real == pytest.approx(0.48364, abs=3e-5) and w.imag == pytest.approx(-0.09676, abs=3e-5)


def test_time_domain_ringdown_matches_leaver():
    f = lambda r: 1 - 2 / r
    fp = lambda r: 2 / r**2
    V = P.potential_function_rs(f, fp, 2.0, l=2)
    t, psi = P.evolve(V, rs_min=-300, rs_max=800, dx=0.05, t_end=250, r0s=20, sigma=3, r_obs_s=10)
    w = P.ringdown_frequency(t, psi, 60, 110, order=8)
    ref = P.leaver_qnm(2, 0, 0)
    assert abs(w.real - ref.real) < 2e-3 and abs(w.imag - ref.imag) < 3e-3
