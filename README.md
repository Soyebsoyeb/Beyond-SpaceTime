# Beyond Spacetime

**A computational investigation of Einstein-Rosen bridges, black-hole singularities, holography, quantum branching, temporal symmetry, and quantum-gravity models.**

[![Python](https://img.shields.io/badge/python-3.13%20tested-3776AB.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Tests](https://img.shields.io/badge/tests-71%20passing-brightgreen.svg)](tests/)
[![Figures](https://img.shields.io/badge/figures-28-blue.svg)](docs/results.md)
[![Units](https://img.shields.io/badge/units-G%20%3D%20c%20%3D%201-informational.svg)](docs/mathematical_conventions.md)
[![Signature](https://img.shields.io/badge/signature-mostly%20plus-informational.svg)](docs/mathematical_conventions.md)
[![Symbolic](https://img.shields.io/badge/symbolic-SymPy-3B5526.svg)](https://www.sympy.org/)
[![Quantum](https://img.shields.io/badge/quantum-QuTiP-orange.svg)](https://qutip.org/)
[![Status](https://img.shields.io/badge/status-research%20code-lightgrey.svg)](docs/limitations_and_open_problems.md)
[![Reproducible](https://img.shields.io/badge/reproducible-one%20command-success.svg)](#reproducing-the-results)
<!-- After pushing to GitHub, add a live CI badge by replacing OWNER:
[![CI](https://github.com/OWNER/beyond-spacetime/actions/workflows/ci.yml/badge.svg)](https://github.com/OWNER/beyond-spacetime/actions/workflows/ci.yml)
-->

---

## Contents

1. [Overview](#overview)
2. [Summary of results](#summary-of-results)
3. [The research questions](#the-research-questions)
4. [Repository layout](#repository-layout)
5. [Installation](#installation)
6. [Reproducing the results](#reproducing-the-results)
7. [How the code fits together](#how-the-code-fits-together)
8. [Module-by-module guide](#module-by-module-guide)
9. [Experiments, figures and data](#experiments-figures-and-data)
10. [Tests](#tests)
11. [Using the library](#using-the-library)
12. [Extending the code](#extending-the-code)
13. [Methods and models](#methods-and-models)
14. [Selected figures](#selected-figures)
15. [Verification strategy](#verification-strategy)
16. [Software and dependencies](#software-and-dependencies)
17. [Scientific scope and honest limits](#scientific-scope-and-honest-limits)
18. [Documentation index](#documentation-index)
19. [Citing this work](#citing-this-work)
20. [License](#license)

---

## Overview

This repository asks what the mathematical structure of black-hole spacetimes, combined with quantum information theory
and carefully specified quantum-gravity models, can tell us about the relationship between spacetime connectivity,
singularities, quantum branching, and the direction of time.

The question is split into eight separately testable problems (A to H). For each one the repository states the theory,
the equations, the assumptions and the regime of validity, implements a computation, and validates it against an
independent method. The project does not try to prove a unified theory. Its purpose is to provide a **verified
computational base** and a precise account of where established results end and speculation begins.

Design principles:

- **Derive and verify, do not assert.** Curvature identities are computed symbolically from the metric. Numerical
  results are compared with closed forms or a second numerical method.
- **One convention set.** Signature, curvature signs, units, and Fourier conventions are fixed in
  [docs/mathematical_conventions.md](docs/mathematical_conventions.md).
- **Separate mathematics from interpretation.** Each claim is labelled as a feature of the equations, an interpretive
  proposal, or a prediction that depends on extra assumptions.
- **Report failures.** Experiments whose expected outcome is failure (floating-point time reversal, for example) are
  reported as such.

## Summary of results

| Topic | Result | Evidence |
| --- | --- | --- |
| Schwarzschild vacuum | $R_{\mu\nu} = 0$ and $K = 48M^2/r^6$ derived symbolically from the metric | `tests/test_schwarzschild.py` |
| Coordinate singularity | Ingoing Eddington-Finkelstein and Kruskal charts are regular at $r = 2M$; $K(2M) = 3/4M^4$ is finite | Figures 1 to 3 |
| Time left inside | Proper time from the horizon to $r = 0$ is $4M/3$ for $E = 1$ and at most $\pi M$ for any timelike worldline | Figure 4, exact to $10^{-12}$ |
| Tides | $\mathcal E = \operatorname{diag}(-2M/r^3, M/r^3, M/r^3)$, independent of $E$; ODE matches exact neighbouring geodesics to $9\times10^{-4}$ | Figures 5 and 6 |
| Einstein-Rosen bridge | Throat of area $16\pi M^2$ on the $T = 0$ slice; radius shrinks on later slices and vanishes at $T = 1$ | Figures 8 and 9 |
| Traversability | No future-directed causal curve from region I to III (22,500 sampled pairs, 0 violations of the necessary condition) | Figure 10 |
| Regular black hole | Hayward core: $K \to 24/l^4$; two horizons for $l < 0.770M$; NEC and WEC hold, SEC fails inside $r < (Ml^2)^{1/3}$, DEC fails outside $r > (4Ml^2)^{1/3}$ | Figures 11 to 14 |
| Quasinormal modes | Leaver: $M\omega = 0.37367 - 0.08896i$ ($\ell = 2$ gravitational), matches the literature to about $10^{-5}$ and a time-domain run to $10^{-5}$ | Figure 26 |
| Holographic entropy | Ryu-Takayanagi in BTZ equals the Calabrese-Cardy formula to $1.4\times10^{-9}$; homology transition at $\Delta\phi^* = 5.821$ for $r_+ = 1.5$ | Figure 15 |
| ER bridge growth | Length follows $2\ln\cosh(r_+t)$, slope $2r_+$ at late times | Figure 16 |
| Entanglement and signalling | Thermofield-double reduced state is thermal; no signalling without an explicit coupling (deviation $3\times10^{-16}$) | Figure 17 |
| Quantum branching | Decoherence matches the exact product formula to $7\times10^{-10}$; density matrices stay valid | Figures 18 and 19 |
| Time reversal | Wave equation reversible, heat equation not; cat map reversal exact on a lattice, fails in floating point with rate $\ln\frac{3+\sqrt5}{2}$ | Figures 20 and 21 |
| Black and white holes | $T \to -T$ is an isometry of the Kruskal extension and exchanges regions II and IV | Figure 22 |
| Wormholes and CTCs | Ellis throat violates the radial NEC ($\rho - \tau = -1/4\pi b_0^2$); CTC iff clock lag exceeds separation; Goedel CTCs for $r > 0.8814$ | Figures 23 to 25 |

All numbers are produced by `experiments/run_all.py` and stored in `results/data/*.json`.

## The research questions

| Question | Topic | Where answered |
| --- | --- | --- |
| A | What an infalling observer experiences after crossing the horizon | [Question A](docs/research_questions.md#question-a-black-hole-interiors) |
| B | Why the Einstein-Rosen bridge is not traversable | [Question B](docs/research_questions.md#question-b-einstein-rosen-bridges) |
| C | Geodesic incompleteness versus curvature divergence | [Question C](docs/research_questions.md#question-c-singularities) |
| D | When entanglement can be related to geometry | [Question D](docs/research_questions.md#question-d-quantum-information-and-holography) |
| E | Many-worlds, inflationary and geometric multiverses compared | [Question E](docs/research_questions.md#question-e-quantum-branching-and-multiverse-proposals) |
| F | Time-reversal symmetry versus irreversibility | [Question F](docs/research_questions.md#question-f-time-reversal) |
| G | Spacetimes with closed timelike curves | [Question G](docs/research_questions.md#question-g-time-machine-proposals) |
| H | A regular-black-hole effective model | [Question H](docs/research_questions.md#question-h-quantum-gravity-extensions) |

## Repository layout

```text
beyond-spacetime/
|-- README.md
|-- LICENSE
|-- CITATION.cff
|-- pyproject.toml
|-- requirements.txt
|-- environment.yml
|-- Makefile
|-- docs/
|   |-- mathematical_conventions.md
|   |-- software_and_dependencies.md
|   |-- research_questions.md
|   |-- verification_matrix.md
|   |-- api_reference.md       generated from docstrings
|   |-- results.md
|   |-- multiverse_comparison.md
|   |-- limitations_and_open_problems.md
|   |-- references.md
|   `-- ai_assistance.md
|-- src/beyond_spacetime/
|   |-- symbolic.py          symbolic curvature pipeline (SymPy)
|   |-- coordinates.py       Schwarzschild, EF, Kruskal, Penrose
|   |-- geodesics.py         radial timelike geodesics
|   |-- tides.py             tidal tensor and geodesic deviation
|   |-- embedding.py         Flamm paraboloid and slicings
|   |-- causal.py            causal-structure tests
|   |-- regular.py           Hayward effective regular black hole
|   |-- perturbations.py     Leaver QNMs and time-domain evolution
|   |-- holography.py        BTZ, Ryu-Takayanagi, ER-bridge growth
|   |-- tfd.py               thermofield-double toy model
|   |-- quantum.py           system-environment branching (QuTiP)
|   |-- timereversal.py      time reversal and irreversibility
|   |-- wormholes.py         Morris-Thorne and Ellis wormholes
|   |-- chronology.py        closed timelike curves
|   |-- multiverse.py        classification of multiverse concepts
|   |-- diagrams.py          Kruskal and Penrose drawing helpers
|   |-- conventions.py       constants and tolerances
|   `-- plotting.py          shared plot style
|-- scripts/
|   `-- build_api_docs.py    regenerates docs/api_reference.md
|-- tests/                   71 automated tests (4 files)
|-- experiments/             exp01 to exp10 and run_all.py
|-- notebooks/               quick-start notebook
`-- results/
    |-- figures/             28 PNG figures
    `-- data/                JSON result files
```

## Installation

Developed and tested on Python 3.13. `pyproject.toml` declares Python 3.10 or newer and the CI workflow covers 3.10 to 3.13, but only 3.13 has been run by the author of this build.

```bash
git clone https://github.com/OWNER/beyond-spacetime.git
cd beyond-spacetime
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

Conda users can instead run `conda env create -f environment.yml`.

## Reproducing the results

```bash
make test        # run the 71 automated tests
make experiments # regenerate all figures and result files (about one minute)
make all         # both
```

Without `make`:

```bash
pytest
python experiments/run_all.py
```

The pinned versions in `requirements.txt` are those used to produce the committed figures and numbers.

## How the code fits together

The package `beyond_spacetime` (in `src/`) is a flat library of pure functions plus one class. Nothing holds global
state, so every function can be called on its own from a notebook. Experiments and tests are thin clients of the library.

```text
                      +-----------------------+
                      |  conventions.py       |  constants, SI conversions
                      |  plotting.py          |  style and save()
                      +-----------+-----------+
                                  |
   +------------------------------+--------------------------------+
   |                              |                                |
   v                              v                                v
 symbolic.py  (SymPy Spacetime)   coordinates.py / diagrams.py      multiverse.py
   |  independent curvature          |  charts: EF, Kruskal, Penrose   (pure data)
   |  reference for everything       v
   |                              geodesics.py -> tides.py
   |                              embedding.py,  causal.py
   +--> regular.py -> perturbations.py          (Hayward model, QNMs)
   +--> wormholes.py -> chronology.py           (Morris-Thorne, CTCs)
        holography.py -> tfd.py                 (BTZ, RT, toy TFD)
        quantum.py, timereversal.py             (QuTiP branching, reversal)
                                  |
                                  v
            experiments/expNN_*.py  --> results/figures/*.png, results/data/*.json
                                  |
                                  v
            tests/test_*.py  (check every module against an independent path)
```

Data flow for one result, using the tidal tensor as an example:

1. `symbolic.Spacetime` builds the Christoffel symbols, Riemann tensor and Kretschmann scalar of Schwarzschild.
2. `tides.frame_tidal_tensor` lambdifies the lower-index Riemann tensor, builds the comoving orthonormal frame
   $(u, e_1, e_2, e_3)$ of a radial geodesic and contracts $\mathcal E_{ij} = R_{\mu\nu\alpha\beta}\,e_i^\mu u^\nu e_j^\alpha u^\beta$.
3. `tides.tidal_closed_form` returns $\operatorname{diag}(-2M/r^3, M/r^3, M/r^3)$. The test suite requires the two to agree
   for several $(r, E)$, which also proves the tensor is independent of $E$.
4. `experiments/exp02_infall_tides.py` uses `geodesics.integrate_infall` and `tides.integrate_deviation` to produce
   Figures 3 to 7 and `results/data/exp02_infall_tides.json`.

Conventions used everywhere: $G = c = 1$, signature $(-,+,+,+)$, $R^a{}_{bcd} = \partial_c\Gamma^a_{db} - \dots$,
$R_{bd} = R^a{}_{bad}$. Details are in [docs/mathematical_conventions.md](docs/mathematical_conventions.md).

## Module-by-module guide

Each subsection states what the module computes, the mathematics behind its main functions, and where it is tested and
used. Signatures and full docstrings of every public function are in [docs/api_reference.md](docs/api_reference.md).

### `symbolic.py`: symbolic curvature

Class `Spacetime(coords, metric)` takes any SymPy metric matrix and provides `christoffel()`, `riemann_mixed()`,
`riemann_lower()`, `ricci()`, `ricci_scalar()`, `einstein()`, `kretschmann()` and `vacuum_residual()`. It implements

$$
\Gamma^a_{bc} = \tfrac12 g^{ad}(\partial_c g_{db} + \partial_b g_{dc} - \partial_d g_{bc}), \qquad
K = R^{ab}{}_{cd}R^{cd}{}_{ab}.
$$

Factories: `schwarzschild()`, `static_spherical(f)` for $-f\,dt^2 + dr^2/f + r^2d\Omega^2$, `hayward_f()`, `godel_metric()`.
Closed forms used as cross-checks: `kretschmann_static_formula(f)` $= f''^2 + 4f'^2/r^2 + 4(1-f)^2/r^4$ and
`einstein_mixed_static(f)`. This module is the independent reference that validates `tides`, `regular` and `wormholes`.
Tests: `test_schwarzschild.py`, `test_causal_and_regular.py`, `test_timereversal_wormholes.py`. Experiments: 1, 4, 8.

### `coordinates.py`: charts of the maximal extension

Maps between Schwarzschild $(t, r)$, ingoing and outgoing Eddington-Finkelstein, Kruskal-Szekeres $(T, X)$ and Penrose
coordinates. Key functions: `tortoise`, `ingoing_ef_v`, `outgoing_ef_u`, `schwarzschild_to_kruskal(t, r, M, region)`,
`kruskal_to_r` (solves $(r/2M - 1)e^{r/2M} = X^2 - T^2$ with the Lambert $W$ function), `kruskal_conformal_factor`
$F = (32M^3/r)e^{-r/2M}$, `region_of(T, X)` (labels I to IV), `to_penrose` and `penrose_to_kruskal`
(compactification by $\arctan$ of the null coordinates). Tests: Kruskal pull-back reproduces the Schwarzschild metric, the
$(t, r) \to (T, X) \to r$ round trip holds in all four regions, and the singularity maps to a straight line.

### `diagrams.py`: drawing helpers

`draw_kruskal`, `draw_penrose_frame`, `penrose_curve` and `draw_constant_r_and_t` draw spacetime diagrams on Matplotlib
axes. They contain no physics beyond calls into `coordinates.py`. Used by experiments 1, 3, 7.

### `geodesics.py`: radial timelike geodesics

Implements the first integral and the equation of motion

$$
\left(\frac{dr}{d\tau}\right)^2 = E^2 - \left(1 - \frac{2M}{r}\right), \qquad \ddot r = -\frac{M}{r^2}.
$$

`integrate_infall(r0, E, M, r_end)` integrates with DOP853 and returns the arrays `tau`, `r`, `rdot`, `v` (advanced time),
plus `tau_end` and `tau_horizon`. Closed forms for validation: `tau_E1` $= 2(r_1^{3/2} - r_2^{3/2})/(3\sqrt{2M})$,
`tau_rest_at` and `tau_rest_to_singularity` (cycloid, $\tfrac\pi2\sqrt{r_0^3/2M}$), `tau_inside` (horizon to singularity,
maximum $\pi M$ at $E = 0$), `energy_from_rest`, `schwarzschild_time` and `kretschmann`. Tests: analytic proper times to
$10^{-9}$ or better; the $\pi M$ bound.

### `tides.py`: tidal tensor and geodesic deviation

`frame_tidal_tensor(r, E, M)` evaluates $\mathcal E_{ij}$ from the full symbolic Riemann tensor in the comoving frame,
`tidal_closed_form` gives $\operatorname{diag}(-2M/r^3, M/r^3, M/r^3)$, and `integrate_deviation` solves

$$
\ddot\xi_\parallel = \frac{2M}{r^3}\,\xi_\parallel, \qquad \ddot\xi_\perp = -\frac{M}{r^3}\,\xi_\perp
$$

along the infall. `neighbouring_geodesic_separation` is the independent validation: the separation of two exact geodesics
with energies $E$ and $E + dE$. `spaghettification_radius(M_solar, body_length_m, limit_g)` converts to SI through
`conventions.tidal_accel_si`.

### `embedding.py`: the Einstein-Rosen bridge geometry

`flamm_z` and `flamm_dzdr` give the Flamm paraboloid $z = 2\sqrt{2M(r - 2M)}$; `proper_radial_distance` and `isotropic_r`
describe the time-symmetric slice; `kruskal_slice_profile` and `throat_radius_kruskal` compute the throat on slices
$T = T_0$, where $(r_{\min}/2M - 1)e^{r_{\min}/2M} = -T_0^2$; `painleve_gullstrand_profile` shows that a flat slicing has no
throat at all, which demonstrates the slicing dependence.

### `causal.py`: why the bridge is not traversable

`null_coords` returns $(u, v) = (T - X, T + X)$; `may_be_in_causal_future` implements the necessary condition
$\Delta u \ge 0$ and $\Delta v \ge 0$; `sample_region` and `count_signalling_pairs` test random pairs from region I to
III (the count must be zero); `ellis_crossing_time` shows that the traversable Ellis wormhole does connect both sides.

### `regular.py`: the Hayward regular black hole (phenomenological)

With mass function $m(r) = Mr^3/(r^3 + 2Ml^2)$ and $f = 1 - 2m/r$, the module provides `f`, `fprime`, `fsecond`,
`kretschmann` and `kretschmann_core_limit` ($24/l^4$), `horizons` (roots of $r^3 - 2Mr^2 + 2Ml^2 = 0$),
`l_critical_over_M` ($4/3\sqrt3 \approx 0.770$), `hawking_temperature` $f'(r_+)/4\pi$, and the effective fluid
$\rho = m'/4\pi r^2$, $p_r = -\rho$, $p_t = -\rho - r\rho'/2$ in `density`, `stress_energy` and `energy_conditions`
(NEC, WEC, SEC, DEC as pointwise booleans). `tau_outer_horizon_to_turning_point` and `tau_to_center` give proper times
of radial geodesics. Test: the Kretschmann scalar equals the full symbolic Riemann computation, and the effective stress
tensor equals the symbolic Einstein tensor.

### `perturbations.py`: quasinormal modes

`leaver_qnm(l, s, n)` solves Leaver's continued fraction with the $n$-th inversion. The time-domain pipeline is
`rw_potential_schwarzschild` (or `scalar_potential_static`), `tortoise_grid`, `potential_function_rs`, `evolve`
(leapfrog for $\partial_t^2\psi - \partial_{r_*}^2\psi + V\psi = 0$), `prony` and `ringdown_frequency`. The test compares both
paths, and the continued fraction with published values.

### `holography.py`: BTZ, Ryu-Takayanagi and the ER bridge

`central_charge` ($c = 3\ell/2G_N$); the geodesic data `geodesic_J_from_angle`, `angle_from_J_numeric`, `length_numeric`
(quadrature) and `geodesic_length_closed` $= 2\ln[(2r_c/r_+)\sinh(r_+\Delta\phi/2)]$; entropies `entanglement_entropy_RT`,
`entanglement_entropy_CFT` (Calabrese-Cardy) and `rt_entropy` (with the homology constraint
$\min[L(\Delta\phi), L(2\pi - \Delta\phi) + 2\pi r_+]/4G_N$); `transition_angle`; the vacuum limit; embedding-space tools
`btz_embedding` and `embedding_distance`; two-sided growth `two_sided_length` and `two_sided_length_closed`
$= 2\ln[(2r_c/r_+)\cosh(r_+t)]$; and `bh_entropy` versus `cardy_entropy`.

### `tfd.py`: thermofield-double toy model

`ising_hamiltonian` builds a transverse-field Ising chain; `tfd_state` returns
$|\mathrm{TFD}\rangle = Z^{-1/2}\sum_n e^{-\beta E_n/2}|n\rangle_L|n^*\rangle_R$; `entropy_of_reduced_L` uses the QuTiP partial trace
and `reduced_R_numpy` an independent reshape; `no_signalling_deviation` applies a random unitary on $L$ and measures the
change of $\rho_R$; `signal_with_coupling` shows that an explicit $L$-$R$ coupling does transmit a signal. It is a quantum
information toy, not a holographic computation.

### `quantum.py`: system-environment branching

`initial_state`, `dephasing_hamiltonian` ($H = \sum_k g_k Z_SZ_k$) and `noncommuting_hamiltonian` build the models;
`evolve_states` runs QuTiP `sesolve`; `reduced_system` traces out the environment; `decoherence_factor` gives the exact
$\prod_k\cos(2g_kt)$; `environment_overlap_numeric` rebuilds it from explicit environment states; `l1_coherence` and
`relative_entropy_coherence` quantify coherence; `check_density_matrix` verifies trace, Hermiticity and positivity;
`dense_evolve` is the SciPy `expm` reference; `branch_weights` returns Born weights.

### `timereversal.py`: reversal and irreversibility

`leapfrog_reversal_error` (classical), `wave_reversal_error` versus `heat_reversal_error`, `theta_spin_half` and
`quantum_reversal_check` (antiunitary $\Theta = i\sigma_yK$), `real_hamiltonian_reversal_error`,
`kruskal_time_reversal_invariance` ($T \to -T$ isometry), and the Arnold cat map: `coarse_grained_entropy`, `cat_map_float`,
`cat_map_int` and `cat_experiment`, which contrasts exact integer reversal with the floating-point error growing at
$\lambda = \ln\frac{3+\sqrt5}{2}$.

### `wormholes.py`: Morris-Thorne and Ellis-Bronnikov

`mt_metric(b, Phi)` builds the metric; `mt_orthonormal_einstein` computes the frame Einstein tensor from it and
`mt_matter` the closed forms

$$
\rho = \frac{b'}{8\pi r^2}, \qquad \tau = \frac{b/r - 2(r - b)\Phi'}{8\pi r^2}, \qquad
p = \frac r2\left[(\rho - \tau)\Phi' - \tau'\right] - \tau .
$$

`nec_radial_at_throat` gives $\rho - \tau = -1/4\pi b_0^2$ for Ellis; `ellis_b`, `ellis_areal_radius`,
`ellis_embedding_z`, `ellis_tidal_lateral` and `ellis_boosted_tidal_symbolic` (independent symbolic check) complete it.

### `chronology.py`: closed timelike curves

`clock_lag` $\Delta = T(1 - \sqrt{1 - \beta^2})$, `ctc_exists` ($\Delta > D$), `chronology_horizon_time(D, beta)`
(solves $\Delta = D$), `ctc_map`, `minimal_loop_time`, and for the Goedel universe `godel_gphiphi` and
`godel_circle_is_closed_timelike` (true for $r > \ln(1 + \sqrt2)$).

### `multiverse.py`, `conventions.py`, `plotting.py`

`multiverse.CONCEPTS` and `markdown_table()` hold the classification of four notions of multiple worlds.
`conventions` holds SI constants, `geometric_mass_m` and `tidal_accel_si`. `plotting` provides `apply_style`, `save` and the
colour-blind-safe `PALETTE`.

## Experiments, figures and data

`experiments/run_all.py` imports each `expNN_*.py`, calls its `run()` function, writes the returned dictionary to
`results/data/<name>.json` (through `experiments/common.record`) and logs timings in `results/data/run_summary.json`.
Each experiment writes its figures itself with `plotting.save`.

| Experiment | Topic | Figures | Library modules |
| --- | --- | --- | --- |
| `exp01_geometry.py` | Metric, curvature, Kruskal and Penrose | 1, 2 | symbolic, coordinates, diagrams |
| `exp02_infall_tides.py` | Infall, proper time, tides | 3 to 7 | geodesics, tides, conventions |
| `exp03_bridge_causality.py` | Flamm paraboloid, slicings, causal structure | 8 to 10 | embedding, causal, diagrams |
| `exp04_regular_bh.py` | Hayward model | 11 to 14 | regular, perturbations |
| `exp05_holography.py` | RT, ER-bridge growth, TFD | 15 to 17 | holography, tfd |
| `exp06_branching.py` | Decoherence and coherence | 18, 19 | quantum |
| `exp07_time_reversal.py` | Reversal, cat map, Kruskal reversal | 20 to 22 | timereversal, diagrams |
| `exp08_wormholes_chronology.py` | Ellis wormhole, CTCs, Goedel | 23 to 25 | wormholes, chronology, symbolic |
| `exp09_qnm_multiverse.py` | QNMs and multiverse table | 26, 27 | perturbations, multiverse |
| `exp10_dependency_graph.py` | Equation, module and figure map | 28 | none (drawing only) |

## Tests

Run `pytest -v`. The 71 tests are grouped by topic:

| File | Tests | Covers |
| --- | --- | --- |
| `tests/test_schwarzschild.py` | 24 | Vacuum equations, Kretschmann, Kruskal and Penrose maps, proper times, tides, deviation, Flamm embedding, SI units |
| `tests/test_causal_and_regular.py` | 12 | Causal structure, Hayward curvature, stress-energy, energy conditions, horizons |
| `tests/test_holography_quantum.py` | 19 | Leaver and time-domain QNMs, BTZ geodesics, RT and Cardy, thermofield double, decoherence, density matrices |
| `tests/test_timereversal_wormholes.py` | 16 | Time reversal, cat map, Morris-Thorne, Ellis, CTCs, Goedel |

The test-to-equation mapping is in [docs/verification_matrix.md](docs/verification_matrix.md).

## Using the library

All snippets below were executed against this repository.

```python
import numpy as np
from beyond_spacetime import symbolic as S, coordinates as C, geodesics as G
from beyond_spacetime import regular as R, holography as H, perturbations as P, quantum as Q, chronology as CH

# 1. Curvature from the metric: Ricci tensor vanishes, Kretschmann is 48 M^2 / r^6
sch = S.schwarzschild()
sch.ricci()                  # zero 4x4 matrix
sch.kretschmann()            # 48*M**2/r**6

# 2. Kruskal coordinates and region label for a point at r = 3M, t = 0
T, X = C.schwarzschild_to_kruskal(t=0.0, r=3.0, M=1.0, region="I")
C.region_of(T, X)            # 'I'

# 3. Radial infall from r = 10M with E = 1; proper time from the horizon to the singularity
traj = G.integrate_infall(r0=10.0, E=1.0, M=1.0, r_end=1e-3)
G.tau_inside(1.0, M=1.0)     # 4/3 for E = 1

# 4. Regular black hole: horizons for l = 0.5 M and the extremal ratio
R.horizons(M=1.0, l=0.5)     # [0.597, 1.855]
R.l_critical_over_M()        # 0.7698

# 5. Holographic entropy of an interval, including the homology constraint
H.rt_entropy(dphi=1.0, r_plus=1.5, r_cut=1e5)

# 6. Fundamental gravitational quasinormal mode
P.leaver_qnm(l=2, s=2, n=0)  # (0.37367 - 0.08896j)

# 7. Exact dephasing factor and chronology horizon
Q.decoherence_factor(np.array([0.3, 0.5, 0.7]), 0.4)
CH.chronology_horizon_time(D=1.0, beta=0.8)   # 2.5
```

A longer walk-through is in `notebooks/quickstart.ipynb`.

## Extending the code

1. **New metric.** Create it with `symbolic.static_spherical(f)` (or `Spacetime(coords, metric)`), check
   `ricci()` and `kretschmann()`, and compare with `kretschmann_static_formula(f)`.
2. **New model module.** Add `src/beyond_spacetime/<name>.py` with a module docstring that states the equations, give each
   public function a docstring, and write at least one test that compares the result with an independent path.
3. **New experiment.** Add `experiments/expNN_<name>.py` with `run()` returning a JSON-serialisable dictionary, save figures
   with `plotting.save(fig, name, FIG_DIR)`, and append the module to `MODULES` in `experiments/run_all.py`.
4. **Documentation.** Add a row to [docs/verification_matrix.md](docs/verification_matrix.md) for every new equation and
   regenerate the API reference with `python scripts/build_api_docs.py`.

## Methods and models

### Classical relativity

The Schwarzschild metric

$$
ds^2 = -\left(1 - \frac{2M}{r}\right)dt^2 + \left(1 - \frac{2M}{r}\right)^{-1}dr^2 + r^2 d\Omega^2
$$

is analysed with a custom SymPy curvature pipeline (`symbolic.py`) that fixes the sign conventions of the project. The
Kretschmann scalar $K = R_{abcd}R^{abcd} = 48M^2/r^6$ is derived, not assumed. The maximal extension uses Kruskal-Szekeres
coordinates with the radius recovered through the Lambert $W$ function, and Penrose diagrams use the compactification
$\tilde u = \arctan(T - X)$, $\tilde v = \arctan(T + X)$.

For radial geodesics with conserved energy $E$,

$$
\left(\frac{dr}{d\tau}\right)^2 = E^2 - \left(1 - \frac{2M}{r}\right), \qquad
\frac{d^2 r}{d\tau^2} = -\frac{M}{r^2},
$$

integrated with DOP853 and validated against $\tau = 2(r_1^{3/2} - r_2^{3/2})/(3\sqrt{2M})$ for $E = 1$ and the cycloid
solution for infall from rest. Tidal effects follow from

$$
\frac{D^2\xi^\mu}{d\tau^2} = -R^\mu{}_{\nu\alpha\beta}\,u^\nu\xi^\alpha u^\beta .
$$

### The Einstein-Rosen bridge

The equatorial slice of the time-symmetric geometry embeds in $\mathbb R^3$ as

$$
z(r) = 2\sqrt{2M(r - 2M)}, \qquad \left(\frac{dz}{dr}\right)^2 = \frac{2M}{r - 2M}.
$$

Causal non-traversability is argued from the null coordinates of the Kruskal chart, not from the embedding, and is tested
numerically. Slicing dependence of the throat is computed explicitly.

### Singularity resolution model

The Hayward metric

$$
f(r) = 1 - \frac{2Mr^2}{r^3 + 2Ml^2}
$$

replaces the singularity with a de Sitter core. Its effective stress-energy is read from the Einstein tensor and tested
against the classical energy conditions. The model is treated as **phenomenological**; it is not presented as a prediction of
any fundamental theory, and finite curvature is not taken as proof of a complete resolution.

### Quasinormal modes

Leaver's continued fraction (with the $n$-th inversion for overtones) is compared with a time-domain leapfrog evolution of
$\partial_t^2\psi - \partial_{r_*}^2\psi + V\psi = 0$ analysed with a Prony fit.

### Holography

BTZ black hole in AdS$_3$ with Brown-Henneaux central charge $c = 3\ell/2G_N$. The Ryu-Takayanagi formula
$S_A = \operatorname{Length}(\gamma_A)/4G_N$ is applied only to static slices in the large-$c$ regime, including the
homology constraint

$$
S_A = \frac{1}{4G_N}\min\Big[L(\Delta\phi),\; L(2\pi - \Delta\phi) + 2\pi r_+\Big].
$$

The two-sided bridge length is checked with exact embedding-space distances. A finite-dimensional thermofield-double
model illustrates purification of a thermal state and the no-signalling theorem. It is explicitly labelled a toy.

### Quantum branching

A system qubit interacts unitarily with an environment of up to 13 qubits. For pure dephasing the exact off-diagonal element
is $\rho_{01}(t) = \alpha\beta^*\prod_k\cos(2g_kt)$; QuTiP is validated against it and against dense matrix exponentials.
Density-matrix normalization, Hermiticity and positivity are verified at every time step.

### Time reversal and chronology

Time reversal is defined separately for trajectories, wave equations, quantum evolution (antiunitary $\Theta$) and
the Kruskal geometry. Irreversibility is demonstrated with a coarse-grained entropy on a chaotic map. Closed timelike curves
are analysed in a wormhole clock-lag model (exists iff $\Delta = T(1 - \sqrt{1 - \beta^2}) > D$) and the Goedel solution.

## Selected figures

![Kruskal and Penrose diagrams](results/figures/fig02_kruskal_penrose.png)

*Kruskal-Szekeres and Penrose diagrams of the maximally extended Schwarzschild solution.*

![Tides and geodesic deviation](results/figures/fig05_tides_and_deviation.png)

*Curvature and tidal growth along a radial infall, and the geodesic-deviation solution.*

![Causal structure](results/figures/fig10_causal_structure.png)

*The causal future of an event in region I never reaches region III; collapse spacetimes lack the second exterior;
a traversable Ellis wormhole does connect both sides.*

![Hayward energy conditions](results/figures/fig12_hayward_energy_conditions.png)

*Effective stress-energy and energy conditions for the regular black hole.*

![Ryu-Takayanagi in BTZ](results/figures/fig15_rt_btz.png)

*Holographic entanglement entropy with the homology phase transition.*

![Quasinormal modes](results/figures/fig26_qnm_ringdown.png)

*Time-domain ringdown against the Leaver continued fraction.*

![Irreversibility](results/figures/fig20_irreversibility.png)

*Entropy growth and fragile reversibility in a chaotic map.*

All 28 figures with captions are in [docs/results.md](docs/results.md).

## Verification strategy

Every important equation has a row in [docs/verification_matrix.md](docs/verification_matrix.md) stating its meaning,
assumptions, source, dimensional check, independent verification, implementing function and test. Highlights:

| Quantity | Path 1 | Path 2 |
| --- | --- | --- |
| Curvature invariants | Symbolic contraction of the full Riemann tensor | Generic closed form for $f(r)$ metrics |
| Tidal tensor | Closed form | Symbolic Riemann projected onto the comoving frame |
| Geodesic deviation | ODE | Separation of exact neighbouring geodesics |
| Quasinormal modes | Continued fraction | Time-domain evolution |
| Geodesic length in BTZ | Numerical quadrature | Closed form and embedding-space distance |
| Quantum evolution | QuTiP `sesolve` | SciPy `expm` and analytic decoherence factor |

Run `pytest -v` to see the 71 tests by name.

## Software and dependencies

| Package | Role | Version used |
| --- | --- | --- |
| NumPy | arrays and linear algebra | 2.5.3 |
| SciPy | ODEs, quadrature, root finding, Lambert $W$, matrix exponentials | 1.18.1 |
| SymPy | symbolic tensor calculus | 1.14.0 |
| Matplotlib | figures | 3.11.2 |
| QuTiP | open quantum systems | 5.3.1 |
| pytest | testing | 9.1.1 |

Tools considered but not used (SageMath, xAct, Cadabra, EinsteinPy, Julia DifferentialEquations.jl, the Black Hole
Perturbation Toolkit, Qiskit) are discussed in [docs/software_and_dependencies.md](docs/software_and_dependencies.md),
together with licenses and the status of their verification.

## Scientific scope and honest limits

This repository does **not**:

- claim that Einstein-Rosen bridges connect many-worlds branches. No model establishing such a relation is implemented;
- present ER=EPR as established for all entangled systems. It is examined in the BTZ bridge and a toy model only;
- present the Hayward metric as a prediction of quantum gravity, or finite curvature as proof that the singularity
  problem is solved;
- demonstrate a working time machine or any reversal of an observer's experienced time. Closed timelike curves are shown to
  exist in assumed idealized models that require exotic matter;
- apply the Ryu-Takayanagi formula outside static, large-$c$ settings.

The full list of limitations and open problems is in
[docs/limitations_and_open_problems.md](docs/limitations_and_open_problems.md). The reference list in
[docs/references.md](docs/references.md) was written without live database access and must be checked before citation.

## Documentation index

| Document | Content |
| --- | --- |
| [Mathematical conventions](docs/mathematical_conventions.md) | Signature, curvature signs, units, coordinates, quantum and holographic conventions, tolerances |
| [Software and dependencies](docs/software_and_dependencies.md) | Selected and rejected tools with versions and licenses |
| [Research questions A to H](docs/research_questions.md) | Theory, equations, assumptions, methods, findings and limits for each question |
| [API reference](docs/api_reference.md) | Signature and docstring of every public function, with the tests and experiments that use each module |
| [Verification matrix](docs/verification_matrix.md) | Equation-to-code-to-test entries and dependency graph |
| [Results](docs/results.md) | All figures with captions and an uncertainty table |
| [Multiverse comparison](docs/multiverse_comparison.md) | Four notions of multiple worlds compared |
| [Limitations and open problems](docs/limitations_and_open_problems.md) | Scope limits and suggested extensions |
| [References](docs/references.md) | Bibliography |
| [AI assistance and provenance](docs/ai_assistance.md) | How AI-generated material was verified |

## Citing this work

A `CITATION.cff` file is provided, so GitHub shows a "Cite this repository" button. Replace the placeholder author entry with
the project owner's details before publishing.

## License

Released under the MIT License. See [LICENSE](LICENSE).
