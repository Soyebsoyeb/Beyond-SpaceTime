# Software and Dependencies

This document records what was selected, what was considered and not used, and why. Versions and licenses of the
packages that were actually installed were read from the package metadata in the build environment
(Python 3.13, Linux). Statements about packages that were **not** installed are marked as unverified: this build
environment had no access to verify their current API, maintenance status or licensing, and those entries must be
checked before anyone relies on them.

## 1. Selection principles

1. Use a mature library when the mathematics matches, and write own code only where the goal is to *verify* a
   result independently (for example a continued fraction for quasinormal modes, checked against a time-domain solver).
2. Every library result that enters a figure has an independent check in `tests/` (closed form, second method, or
   symbolic derivation).
3. Keep the dependency set small. Five scientific packages plus `pytest` are enough for everything in this repository.

## 2. Required dependencies (installed and used)

### NumPy 2.5.3

- **License:** BSD-3-Clause (plus bundled permissive components)
- **Purpose:** Arrays, linear algebra
- **Input/output conventions:** Row-major `ndarray`, float64
- **Known limitations:** None relevant
- **Validation in this repository:** Used by all tests
- **Alternative:** None needed

### SciPy 1.18.1

- **License:** BSD
- **Purpose:** ODE solvers (`solve_ivp`, DOP853), quadrature (`quad`), root finding, matrix exponentials, Lambert W
- **Input/output conventions:** Callables with state vector `y`, event functions
- **Known limitations:** Adaptive quadrature emits round-off warnings at $10^{-13}$ requests; handled and checked
- **Validation in this repository:** Geodesics vs analytic proper times; BTZ lengths vs closed forms; Lambert-W round trips
- **Alternative:** Julia DifferentialEquations.jl (not used)

### SymPy 1.14.0

- **License:** BSD
- **Purpose:** Symbolic tensor calculus: Christoffel symbols, Riemann, Ricci, Einstein, Kretschmann
- **Input/output conventions:** Custom `Spacetime` class in `symbolic.py` (SymPy has no built-in tensor-curvature pipeline that fixes our sign conventions)
- **Known limitations:** Slow for large metrics; simplification can be incomplete, so numerical spot checks supplement symbolic zeros
- **Validation in this repository:** $K = 48M^2/r^6$, vacuum equations, generic $f(r)$ formula, Goedel Ricci tensor, Morris-Thorne equations
- **Alternative:** SageMath, Cadabra, Mathematica with xAct (not used)

### Matplotlib 3.11.2

- **License:** Python Software Foundation License
- **Purpose:** All figures
- **Input/output conventions:** Agg backend, PNG output
- **Known limitations:** 3D surface plots are rasterized
- **Validation in this repository:** Figures reviewed by inspection; underlying numbers are tested separately
- **Alternative:** Plotly (not used)

### QuTiP 5.3.1

- **License:** BSD-3-Clause
- **Purpose:** State evolution (`sesolve`), partial traces, von Neumann entropy for the decoherence and thermofield-double models
- **Input/output conventions:** `Qobj` with explicit `dims`
- **Known limitations:** Dense memory growth with environment size; here at most $2^{13}$ dimensions
- **Validation in this repository:** Every QuTiP result cross-checked against NumPy/SciPy `expm` and against the analytic decoherence factor
- **Alternative:** NumPy/SciPy only (implemented as the cross-check)


## 3. Development dependencies

| Package | Version used | License | Purpose |
| --- | --- | --- | --- |
| pytest | 9.1.1 | MIT | Automated tests (`tests/`) |

## 4. Considered and not used

The following tools were named in the project brief. They were **not** installed or exercised in this build. The
license and capability notes are general knowledge and have not been verified here.

| Tool | Reason not used | Status of verification |
| --- | --- | --- |
| SageMath | Heavy dependency; SymPy suffices for the metrics studied | Unverified |
| Mathematica with xAct | Proprietary license required; open-source SymPy pipeline provides the cross-check | Unverified |
| Cadabra | Field-theory tensor manipulation not needed at this scope | Unverified |
| EinsteinPy | Own symbolic pipeline was needed to control conventions; a cross-check against EinsteinPy's Schwarzschild module is a worthwhile extension | Not exercised |
| Julia DifferentialEquations.jl | SciPy DOP853 meets the tolerances; no stiff systems were encountered | Not exercised |
| Black Hole Perturbation Toolkit | Own Leaver implementation reproduces published Schwarzschild QNMs to about $10^{-5}$; the toolkit is the right choice for Kerr | Not exercised |
| Qiskit | No circuit-level experiment adds scientific value here | Not exercised |
| Plotly, Jupyter Book, Sphinx, Quarto | The primary documentation is plain GitHub Markdown and must stay readable without a build step | Not exercised |

## 5. Environment specification

| File | Content |
| --- | --- |
| `pyproject.toml` | Lower bounds for runtime dependencies and the `dev` extra |
| `requirements.txt` | Exact versions used to produce the committed results |
| `environment.yml` | Conda specification with the same pinned versions |

Reproduce the exact environment with `pip install -r requirements.txt && pip install -e .`.
Newer versions of the libraries are expected to work but have not been tested; if a tolerance in `tests/` fails after
an upgrade, inspect the residual before loosening it.

## 6. Optional environments

No optional heavyweight environment is needed. If a future extension adds Kerr perturbations, create a separate
environment so that specialist packages do not enter the core `requirements.txt`.

## 7. Computational assistance record

See [AI assistance and provenance](ai_assistance.md).
