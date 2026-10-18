# Limitations and Open Problems

This document lists what the repository does not do, where the models stop being valid, and which questions remain open.

## 1. Scope limits

| Area | Limitation |
| --- | --- |
| Classical geometry | Only non-rotating, uncharged, eternal or analytically extended Schwarzschild is treated. Rotating (Kerr) and charged (Reissner-Nordstrom) interiors, where inner-horizon instabilities matter, are not computed. |
| Test-body approximation | Geodesics ignore back-reaction, Hawking radiation and the finite size of the observer beyond linear tidal effects. |
| Singularity resolution | Only one effective model (Hayward) is implemented. It is phenomenological. No derivation from a fundamental theory, no dynamical formation from collapse and no inner-horizon stability analysis are included. |
| Energy conditions | Computed pointwise for an effective anisotropic fluid. Averaged energy conditions and quantum energy inequalities are not examined. |
| Holography | AdS$_3$/CFT$_2$ in the BTZ saddle, large $c$, static slices. No dynamical (HRT) surfaces, no quantum corrections, no higher dimensions, no rotating BTZ, no Hawking-Page transition computation. |
| ER=EPR | Demonstrated only through the BTZ bridge length and a finite-dimensional thermofield-double toy. The toy has no bulk geometry. The conjecture is not tested for generic entangled states and is not claimed to hold for them. |
| Traversable wormholes | Ellis-Bronnikov wormhole and a clock-lag model are analysed classically. The Gao-Jafferis-Wall mechanism is described qualitatively through a coupling in the toy model; its bulk calculation is not performed. |
| Quantum branching | Qubit system and qubit environment of at most 13 spins. Large-environment decoherence is extrapolated from the analytic product formula. No measure of branch weights beyond the Born weights is considered. |
| Time reversal | Cat map and finite-difference wave/heat equations are used to illustrate principles. They do not model a realistic gas or the cosmological arrow of time. |
| Chronology | Existence of closed timelike curves in assumed models is shown. Chronology protection is not tested; it remains a conjecture. |

## 2. Numerical limitations

- Quasinormal-mode extraction from time series with a Prony fit is sensitive to the window and the order; the stated
  uncertainty of about $10^{-3}$ is empirical.
- Adaptive quadrature raises round-off warnings at the requested $10^{-13}$ tolerance. They are suppressed in the geodesic
  length routine, and accuracy is verified by test.
- Dense state-vector methods limit the quantum models to about $2^{13}$ dimensions.
- Cutoff-dependent holographic entropies are compared only at equal cutoff.
- The collapse diagram in Figure 10(b) is schematic: the flat interior is drawn to scale only qualitatively.

## 3. Interpretational limits

- Finite curvature is not equivalent to a complete resolution of the singularity problem.
- The existence of closed timelike curves in a model is a mathematical fact, not a statement about technology.
- Unitary branching does not by itself decide between interpretations of quantum mechanics.
- A geometric interpretation of entanglement does not imply that arbitrary entangled systems can exchange signals.

## 4. Open problems relevant to the investigation

1. Which, if any, effective description of a black-hole interior is a controlled limit of a fundamental theory of quantum
   gravity, and does it have a stable inner horizon?
2. How is the interior of a black hole that formed from collapse encoded in the boundary theory, beyond the eternal two-sided case?
3. Can the relation between entanglement and connectivity be stated as a theorem for a class of states beyond the large-$c$
   regime?
4. What is the status of the white-hole region in a realistic, evaporating black hole?
5. Is there any concrete model in which quantum branching and spacetime connectivity are related, as opposed to
   analogies between them?
6. Is chronology protection a theorem of a complete quantum theory of gravity?

## 5. Suggested extensions

- Cross-check the Schwarzschild computations against an independent relativity package and a Black Hole Perturbation Toolkit
  calculation of Kerr quasinormal modes.
- Add a dynamical collapse model (Vaidya) with numerical geodesics through the shell.
- Implement a Jackiw-Teitelboim gravity model with the Gao-Jafferis-Wall coupling for a controlled traversable wormhole.
- Extend the regular black hole to inner-horizon stability and mass-inflation analysis.
- Replace the cat map with a many-particle classical gas to quantify the dependence of the relaxation time on system size.
