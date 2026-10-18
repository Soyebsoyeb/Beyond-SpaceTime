# Comparison of Multiverse Concepts

The phrase "multiple worlds" is used for four mathematically different things. This document states, for each, the object that supplies the multiplicity, whether
that multiplicity is a feature of the equations, an interpretive proposal, or a prediction that depends on extra
assumptions, and whether the regions can exchange signals. The same table is encoded in
`src/beyond_spacetime/multiverse.py` and checked in `tests/test_timereversal_wormholes.py`.

| Concept | Mathematical object | Status of the multiplicity | Signalling between regions | Implemented in |
| --- | --- | --- | --- | --- |
| Quantum branching (many-worlds) | Orthogonal components of a state in a tensor-product Hilbert space | Interpretive proposal | None between decohered components | `quantum.py`, `tfd.py` |
| Eternal inflation and pocket universes | Disconnected regions of a spacetime with a metastable de Sitter vacuum and bubble nucleation | Prediction dependent on extra assumptions | Generically none (cosmological horizons) | Documentation only |
| Second exterior of maximally extended Schwarzschild | Region III of the analytic Kruskal extension | Feature of the equations | None; no future-directed causal curve from I to III | `coordinates.py`, `causal.py` |
| Regions beyond a black-hole bounce or regular core | Extension through a regular core or a black-to-white-hole transition in an effective model | Prediction dependent on extra assumptions | Model dependent | `regular.py` (diagnostics only) |

## 1. Quantum branching

The state vector evolves unitarily. If a system interacts with an environment, the state takes the form
$\sum_i c_i\lvert s_i\rangle\lvert E_i(t)\rangle$ and the environment states become nearly orthogonal,
$\langle E_i\vert E_j\rangle \to 0$. The reduced density matrix then has vanishing off-diagonal elements in the $s_i$ basis.
This is a statement about the equations (see `docs/results.md`, Figure 18). Describing each term as an independent "world" is
an interpretive step that the equations do not force; collapse and hidden-variable interpretations accept the same
mathematics. Nothing in the model assigns a spacetime to a branch.

## 2. Eternal inflation and pocket universes

In inflationary cosmology with a potential containing a metastable vacuum, quantum fluctuations can keep inflation going in
some regions while others nucleate bubbles. The multiplicity is a consequence of a specific field content and potential, and
the assignment of probabilities (the measure problem) is unsettled. The spacetime regions are parts of one global
geometry. No model of this type is implemented in the repository.

## 3. Second exterior of the Kruskal extension

The analytic extension of the Schwarzschild solution contains a second asymptotic region. Its existence follows from the
vacuum equations and analyticity. It is absent from a black hole formed by collapse (Figure 10b), so it is not generic for
astrophysical black holes. Regions I and III are not causally connected.

## 4. Bounce and regular-core models

Effective models in which a black hole interior is replaced by a regular core or a transition to a white hole predict
additional regions only if one adopts the dynamics of that model and extends it beyond the inner horizon. The Hayward
metric implemented here has such a chain of regions in its maximal analytic extension, but the inner horizon is expected to
be unstable, so the extension should not be taken as physical without further analysis.

## 5. What is not claimed

- That Einstein-Rosen bridges connect many-worlds branches. No concrete model in this repository establishes any relation
  between a spacetime bridge and quantum branching.
- That the thermofield-double toy model or the BTZ computation implies a many-worlds structure.
- That any of these multiplicities is observable.
