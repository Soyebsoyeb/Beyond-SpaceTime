# Verification Matrix and Dependency Graph

Each entry links an important equation to its physical meaning, assumptions, derivation or reference, dimensional check,
independent verification, implementation and automated test. File names are relative to `src/beyond_spacetime/` for
code and `tests/` for tests.

## 1. Equation-level matrix

### Classical geometry (Parts III and IV): E1 to E12

#### E1. $R_{\mu\nu} = 0$ for $f = 1 - 2M/r$

- **Meaning and assumptions:** Vacuum solution; spherical, static
- **Derivation or reference:** Direct curvature computation (MTW ch. 31)
- **Dimensional check:** $f$ dimensionless, $M/r$ dimensionless
- **Independent verification:** Full symbolic Ricci tensor simplifies to zero
- **Code:** `symbolic.py`
- **Test:** `test_schwarzschild.py::test_vacuum_einstein_equations`

#### E2. $K = 48M^2/r^6$

- **Meaning and assumptions:** Curvature invariant, coordinate independent
- **Derivation or reference:** Symbolic contraction $R^{ab}{}_{cd}R^{cd}{}_{ab}$
- **Dimensional check:** $[K] = L^{-4}$, $M^2/r^6 = L^{-4}$
- **Independent verification:** Generic formula $f''^2 + 4f'^2/r^2 + 4(1-f)^2/r^4$ for all $f(r)$
- **Code:** `symbolic.py`, `geodesics.py`
- **Test:** `test_kretschmann_scalar`, `test_generic_static_kretschmann_formula`

#### E3. Kruskal metric $-F(dT^2 - dX^2) + r^2d\Omega^2$, $F = 32M^3e^{-r/2M}/r$

- **Meaning and assumptions:** Regular at $r = 2M$; region I chart
- **Derivation or reference:** Coordinate transformation from Schwarzschild
- **Dimensional check:** $F$ has dimension $L^2$ with $T, X$ dimensionless in units of $2M$
- **Independent verification:** Pull-back to $(t, r)$ reproduces $g_{tt}$, $g_{rr}$, $g_{tr} = 0$
- **Code:** `coordinates.py`
- **Test:** `test_kruskal_pullback_metric`

#### E4. $(r/2M - 1)e^{r/2M} = X^2 - T^2$

- **Meaning and assumptions:** Implicit radius
- **Derivation or reference:** Kruskal definition; inverted with Lambert $W$
- **Dimensional check:** dimensionless
- **Independent verification:** Round trip $(t, r) \to (T, X) \to r$ in all four regions
- **Code:** `coordinates.py`
- **Test:** `test_kruskal_roundtrip_r`

#### E5. $\tan\tilde u\tan\tilde v = 1$ on the singularity

- **Meaning and assumptions:** Penrose compactification; singularity is a straight line
- **Derivation or reference:** Algebra
- **Dimensional check:** dimensionless
- **Independent verification:** Numerical map of $T^2 - X^2 = 1$ gives $T_p = \pi/2$
- **Code:** `coordinates.py`
- **Test:** `test_penrose_roundtrip_and_singularity_line`

#### E6. $(dr/d\tau)^2 = E^2 - f$

- **Meaning and assumptions:** First integral of radial geodesics
- **Derivation or reference:** Killing energy conservation
- **Dimensional check:** $E$ dimensionless
- **Independent verification:** Analytic proper times: $E = 1$, cycloid, $\pi M$ bound
- **Code:** `geodesics.py`
- **Test:** `test_horizon_proper_time_analytic_E1`, `test_rest_infall_cycloid`, `test_max_proper_time_inside_horizon_is_pi_M`

#### E7. $\mathcal E = \operatorname{diag}(-2M/r^3, M/r^3, M/r^3)$

- **Meaning and assumptions:** Tidal tensor, independent of $E$
- **Derivation or reference:** Frame computation from full Riemann
- **Dimensional check:** $[\mathcal E] = L^{-2}$
- **Independent verification:** Full symbolic Riemann evaluated in the boosted frame at several $(r, E)$; trace-free
- **Code:** `tides.py`
- **Test:** `test_tidal_tensor_from_full_riemann_matches_closed_form_independent_of_E`, `test_tidal_trace_free_in_vacuum`

#### E8. $\ddot\xi^i = -\mathcal E_{ij}\xi^j$

- **Meaning and assumptions:** Linear geodesic deviation
- **Derivation or reference:** Differentiating the geodesic family
- **Dimensional check:** consistent with E7
- **Independent verification:** Exact neighbouring-geodesic separation; transverse exact solution $\xi_\perp = r\,\delta\theta$
- **Code:** `tides.py`
- **Test:** `test_deviation_matches_neighbouring_geodesics`, `test_transverse_deviation_equals_r_times_dtheta`

#### E9. SI tide $a = c^6\xi/(4G^2M^2)$ at the horizon

- **Meaning and assumptions:** Restores $G$ and $c$
- **Derivation or reference:** Dimensional analysis
- **Dimensional check:** $\mathrm{m\,s^{-2}}$
- **Independent verification:** Independent closed form
- **Code:** `conventions.py`
- **Test:** `test_si_tidal_acceleration_dimensional_check`

#### E10. $(dz/dr)^2 = 2M/(r - 2M)$

- **Meaning and assumptions:** Flamm embedding of the $t$ = const slice
- **Derivation or reference:** Induced metric of a surface of revolution
- **Dimensional check:** $z$ has dimension $L$
- **Independent verification:** Symbolic check $1 + z'^2 = 1/f$
- **Code:** `embedding.py`
- **Test:** `test_flamm_embedding_induces_spatial_metric`

#### E11. $(r_{\min}/2M - 1)e^{r_{\min}/2M} = -T_0^2$

- **Meaning and assumptions:** Throat radius on the Kruskal slice $T = T_0$
- **Derivation or reference:** Evaluate E4 at $X = 0$
- **Dimensional check:** dimensionless
- **Independent verification:** Monotonic shrinking; PG slice has no throat
- **Code:** `embedding.py`
- **Test:** `test_throat_area_and_slicing_dependence`

#### E12. $du \ge 0$, $dv \ge 0$ for causal curves

- **Meaning and assumptions:** No causal curve from region I to III
- **Derivation or reference:** Conformal flatness of the $(T, X)$ sector
- **Dimensional check:** dimensionless
- **Independent verification:** Monte Carlo falsification test
- **Code:** `causal.py`
- **Test:** `test_no_causal_signal_from_region_I_to_III`


### Singularity models and perturbations (Parts V and VIII): E13 to E18

#### E13. $f_H = 1 - 2Mr^2/(r^3 + 2Ml^2)$

- **Meaning and assumptions:** Effective regular black hole
- **Derivation or reference:** Hayward (2006), phenomenological
- **Dimensional check:** $l$ has dimension $L$
- **Independent verification:** Kretschmann against the full Riemann computation
- **Code:** `regular.py`
- **Test:** `test_hayward_kretschmann_matches_full_riemann_computation`

#### E14. $K(r\to0) = 24/l^4$

- **Meaning and assumptions:** de Sitter core
- **Derivation or reference:** Series of E13
- **Dimensional check:** $L^{-4}$
- **Independent verification:** Numerical evaluation at $r = 10^{-4}$
- **Code:** `regular.py`
- **Test:** `test_hayward_core_is_de_sitter_and_finite`

#### E15. $\rho = m'/4\pi r^2$, $p_r = -\rho$, $p_t = -\rho - r\rho'/2$

- **Meaning and assumptions:** Effective stress-energy of $G^\mu{}_\nu = 8\pi T^\mu{}_\nu$
- **Derivation or reference:** Einstein tensor of the $f(r)$ metric
- **Dimensional check:** $L^{-2}$
- **Independent verification:** Symbolic Einstein tensor evaluated at sample radii
- **Code:** `regular.py`
- **Test:** `test_effective_stress_energy_matches_einstein_tensor`

#### E16. SEC fails for $r < (Ml^2)^{1/3}$; DEC fails for $r > (4Ml^2)^{1/3}$

- **Meaning and assumptions:** Energy-condition thresholds
- **Derivation or reference:** Algebra on E15
- **Dimensional check:** $L$
- **Independent verification:** Grid scan against the analytic thresholds
- **Code:** `regular.py`
- **Test:** `test_hayward_energy_condition_thresholds_are_analytic`

#### E17. $l_c/M = 4/(3\sqrt3)$

- **Meaning and assumptions:** Extremal degenerate horizon
- **Derivation or reference:** $f = f' = 0$
- **Dimensional check:** dimensionless
- **Independent verification:** Double root found numerically; $T_H \to 0$
- **Code:** `regular.py`
- **Test:** `test_hayward_horizon_structure_and_extremality`

#### E18. $\ddot\psi - \psi'' + V\psi = 0$, $V_s = f[\ell(\ell+1)/r^2 + (1-s^2)2M/r^3]$

- **Meaning and assumptions:** Linear perturbations on Schwarzschild
- **Derivation or reference:** Regge-Wheeler, Zerilli
- **Dimensional check:** $V$ has $L^{-2}$
- **Independent verification:** Leaver continued fraction vs literature; time-domain vs Leaver
- **Code:** `perturbations.py`
- **Test:** `test_leaver_qnm_known_values`, `test_time_domain_ringdown_matches_leaver`


### Holography and quantum information (Parts VI and VII): E19 to E27

#### E19. $L = 2\ln[(2r_c/r_+)\sinh(r_+\Delta\phi/2)]$

- **Meaning and assumptions:** Connected BTZ geodesic; static slice
- **Derivation or reference:** Geodesic equation with conserved $J = r^2\dot\phi$
- **Dimensional check:** $\ell = 1$
- **Independent verification:** Numerical quadrature of the geodesic integral; exact embedding distance
- **Code:** `holography.py`
- **Test:** `test_btz_geodesic_numeric_equals_closed_form`, `test_embedding_same_side_distance_equals_geodesic`

#### E20. $S_A = L/4G_N$ equals Calabrese-Cardy

- **Meaning and assumptions:** Ryu-Takayanagi; large $c$, BTZ saddle
- **Derivation or reference:** Ryu and Takayanagi (2006)
- **Dimensional check:** nats
- **Independent verification:** Direct comparison at several $r_+$, $\Delta\phi$
- **Code:** `holography.py`
- **Test:** `test_rt_equals_calabrese_cardy`, `test_vacuum_limit`

#### E21. $S_A = \min[L(\Delta\phi), L(2\pi-\Delta\phi) + 2\pi r_+]/4G_N$

- **Meaning and assumptions:** Homology constraint
- **Derivation or reference:** RT with homology
- **Dimensional check:** nats
- **Independent verification:** Continuity at $\Delta\phi^*$; $S_A - S_{A^c} = S_{BH}$
- **Code:** `holography.py`
- **Test:** `test_rt_homology_phase_transition_and_araki_lieb_saturation`

#### E22. $L(t) = 2\ln[(2r_c/r_+)\cosh(r_+t)]$

- **Meaning and assumptions:** ER-bridge growth for $t_L = t_R = t$
- **Derivation or reference:** Hartman and Maldacena (2013)
- **Dimensional check:** $\ell = 1$
- **Independent verification:** Exact embedding-space distance
- **Code:** `holography.py`
- **Test:** `test_two_sided_wormhole_length_growth`

#### E23. $S_{BH} = 2\pi r_+/4G_N = \pi c LT/3$

- **Meaning and assumptions:** Bekenstein-Hawking equals Cardy
- **Derivation or reference:** Brown-Henneaux central charge
- **Dimensional check:** nats
- **Independent verification:** Algebraic identity checked numerically
- **Code:** `holography.py`
- **Test:** `test_bh_entropy_equals_cardy_entropy`

#### E24. $\rho_L = e^{-\beta H}/Z$ for the TFD

- **Meaning and assumptions:** Purification of a thermal state
- **Derivation or reference:** Israel (1976)
- **Dimensional check:** dimensionless
- **Independent verification:** QuTiP partial trace vs thermal state; $(H_L - H_R)\lvert\mathrm{TFD}\rangle = 0$
- **Code:** `tfd.py`
- **Test:** `test_tfd_reduced_state_is_thermal_and_entropy_matches`, `test_tfd_annihilated_by_HL_minus_HR`

#### E25. $\delta\rho_R = 0$ after any operation on $L$

- **Meaning and assumptions:** No signalling
- **Derivation or reference:** Unitary invariance of partial trace
- **Dimensional check:** none
- **Independent verification:** Random unitary on $L$; coupling produces a signal
- **Code:** `tfd.py`
- **Test:** `test_no_signalling_without_coupling_and_signal_with_coupling`

#### E26. $\rho_{01}(t) = \alpha\beta^*\prod_k\cos(2g_kt)$

- **Meaning and assumptions:** Exact pure dephasing
- **Derivation or reference:** Direct solution
- **Dimensional check:** none
- **Independent verification:** QuTiP `sesolve` vs analytic; explicit environment overlap
- **Code:** `quantum.py`
- **Test:** `test_dephasing_matches_analytic_decoherence_factor`

#### E27. Density-matrix validity

- **Meaning and assumptions:** Trace 1, Hermitian, $\rho \ge 0$
- **Derivation or reference:** Axioms of quantum mechanics
- **Dimensional check:** none
- **Independent verification:** QuTiP vs dense `expm`; eigenvalue check
- **Code:** `quantum.py`
- **Test:** `test_noncommuting_model_qutip_equals_dense_expm_and_density_matrix_valid`


### Time reversal, wormholes and chronology (Part VIII): E28 to E35

#### E28. $\Theta e^{-iHt}\Theta^{-1} = e^{+i\Theta H\Theta^{-1}t}$

- **Meaning and assumptions:** Antiunitary time reversal
- **Derivation or reference:** Wigner
- **Dimensional check:** none
- **Independent verification:** Spin-1/2 with $\Theta = i\sigma_yK$
- **Code:** `timereversal.py`
- **Test:** `test_quantum_time_reversal_spin_half`, `test_real_hamiltonian_has_time_reversal_symmetry`

#### E29. $T \to -T$ is an isometry of Kruskal

- **Meaning and assumptions:** Black/white hole exchange
- **Derivation or reference:** Metric depends on $X^2 - T^2$
- **Dimensional check:** none
- **Independent verification:** Numerical invariance of $r$ and $F$
- **Code:** `timereversal.py`
- **Test:** `test_kruskal_time_reversal_is_isometry`

#### E30. Cat-map entropy growth and Lyapunov error growth $\lambda = \ln\frac{3+\sqrt5}{2}$

- **Meaning and assumptions:** Coarse-grained irreversibility; finite-precision reversal
- **Derivation or reference:** Arnold
- **Dimensional check:** none
- **Independent verification:** Exact integer reversal; float error growth rate fitted
- **Code:** `timereversal.py`
- **Test:** `test_cat_map_entropy_growth_and_reversal`

#### E31. $\rho = b'/8\pi r^2$, $\tau = [b/r - 2(r-b)\Phi']/8\pi r^2$, $p = \tfrac r2[(\rho-\tau)\Phi' - \tau'] - \tau$

- **Meaning and assumptions:** Morris-Thorne sources
- **Derivation or reference:** Morris and Thorne (1988)
- **Dimensional check:** $L^{-2}$
- **Independent verification:** Orthonormal Einstein tensor computed from the metric
- **Code:** `wormholes.py`
- **Test:** `test_morris_thorne_equations_ellis`, `test_morris_thorne_equations_with_redshift_function`

#### E32. $\rho - \tau\lvert_{\rm throat} = (b'-1)/8\pi r_0^2 < 0$

- **Meaning and assumptions:** Radial NEC violation
- **Derivation or reference:** Flare-out condition
- **Dimensional check:** $L^{-2}$
- **Independent verification:** Evaluated for the Ellis wormhole
- **Code:** `wormholes.py`
- **Test:** `test_throat_nec_violation_and_flare_out`

#### E33. Lateral tide of a boosted traveller

- **Meaning and assumptions:** Ellis wormhole
- **Derivation or reference:** Riemann in the boosted frame
- **Dimensional check:** $L^{-2}$
- **Independent verification:** Symbolic Riemann
- **Code:** `wormholes.py`
- **Test:** `test_ellis_tidal_formula_matches_symbolic_riemann`

#### E34. CTC iff $\Delta > D$ with $\Delta = T(1 - \sqrt{1-\beta^2})$

- **Meaning and assumptions:** Wormhole time-shift model
- **Derivation or reference:** Morris, Thorne, Yurtsever (1988)
- **Dimensional check:** time $= L$ ($c = 1$)
- **Independent verification:** Chronology-horizon time solves $\Delta = D$ exactly
- **Code:** `chronology.py`
- **Test:** `test_ctc_condition_and_horizon`

#### E35. $g_{\phi\phi} < 0$ for $r > \ln(1+\sqrt2)$

- **Meaning and assumptions:** Goedel closed timelike circles
- **Derivation or reference:** Goedel (1949)
- **Dimensional check:** dimensionless
- **Independent verification:** Ricci tensor has dust form
- **Code:** `chronology.py`, `symbolic.py`
- **Test:** `test_godel_ctc_threshold`, `test_godel_ricci_is_dust_form`


## 2. Dependency graph

The graph connects mathematical results to the code modules, tests and figures that use them. Figure numbers refer to
files in `results/figures/`.

![Dependency graph linking equations, modules and figures](../results/figures/fig28_dependency_graph.png)

**Figure 28.** Each group of equations (left) is implemented in the listed modules (centre) and visualised in the listed figures (right).

## 3. Cross-validation pairs

Whenever a library or a numerical method produces a number, a second, independent path produces the same number:

| Quantity | Path 1 | Path 2 |
| --- | --- | --- |
| Tidal tensor | Closed form $\operatorname{diag}(-2M/r^3, M/r^3, M/r^3)$ | Symbolic Riemann in the comoving frame |
| Deviation growth | ODE for $\xi$ | Separation of exact neighbouring geodesics |
| QNM frequency | Leaver continued fraction | Time-domain evolution with Prony fit |
| BTZ geodesic length | Quadrature | Closed form and embedding-space distance |
| Reduced density matrix | QuTiP partial trace | NumPy reshape and contraction |
| State evolution | QuTiP `sesolve` | SciPy `expm` |
