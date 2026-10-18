# Research Questions A to H

For each question this document states the theory, the equations, the assumptions, the numerical method, the outputs
produced, the limitations, and what the computations do and do not establish. Nothing here claims that a single model
or repository can answer all eight questions together.

Notation: $M$ is the black-hole mass in geometrized units, $f(r) = 1 - 2M/r$.

---

## Question A. Black-hole interiors

**What does general relativity predict for an infalling observer after crossing the event horizon?**

**Theory.** Vacuum general relativity, Schwarzschild solution, radial timelike geodesics.

**Equations.** With conserved energy $E = f\,dt/d\tau$,

$$
\left(\frac{dr}{d\tau}\right)^2 = E^2 - f(r), \qquad \frac{d^2 r}{d\tau^2} = -\frac{M}{r^2}, \qquad
\frac{dv}{d\tau} = \frac{1}{E - \dot r}.
$$

The last relation (ingoing Eddington-Finkelstein time $v$) is regular at $r = 2M$, while the Schwarzschild time diverges.
Inside $r < 2M$ the coordinate $r$ is timelike, so decreasing $r$ is unavoidable for every future-directed worldline.

**Assumptions.** Spherical symmetry, vacuum, test-body motion, no back-reaction, no Hawking radiation.

**Method.** DOP853 integration of $\ddot r = -M/r^2$ with events at $r = 2M$ and a small $r$; comparison with analytic
proper times; quadrature for the remaining proper time.

**Outputs.** Figures 3 to 7 in [results.md](results.md).

**Findings (all verified against closed forms in `tests/test_schwarzschild.py`).**

- The curvature invariant at the horizon is finite: $K(2M) = 3/(4M^4)$. The horizon is not a curvature singularity.
- Proper time from the horizon to $r = 0$ is $\tau = \int_0^{2M} dr/\sqrt{E^2 - f}$. It equals $4M/3$ for $E = 1$ and is
  maximal, $\pi M$, for $E = 0$. Every timelike worldline that crosses the horizon reaches $r = 0$ in at most $\pi M$.
  For a $4\times10^6\,M_\odot$ black hole, $\pi M \approx 62\ \mathrm{s}$; for $10\,M_\odot$ it is about $0.15\ \mathrm{ms}$.
- The tidal tensor along a radial geodesic is independent of $E$: $\mathcal E = \operatorname{diag}(-2M/r^3, M/r^3, M/r^3)$,
  computed here from the full Riemann tensor and not assumed. The deviation equation was validated against exact
  neighbouring geodesics (maximum relative error about $9\times10^{-4}$, limited by finite differencing).
- The tide across 2 m at the horizon is about $2\times10^{8}\ \mathrm{m/s^2}$ for $10\,M_\odot$ and about
  $1.3\times10^{-3}\ \mathrm{m/s^2}$ for $4\times10^{6}\,M_\odot$. A tide of 10 g across 2 m is reached outside the horizon only
  for masses below about $1.5\times10^{4}\,M_\odot$.

**Limitations.** The classical, eternal, test-body description fails near $r \sim \ell_P$. The interior of a realistic
rotating black hole differs qualitatively (inner horizon, mass inflation).

---

## Question B. Einstein-Rosen bridges

**Why does the original Einstein-Rosen bridge fail to provide a traversable passage, and how does its causal structure differ from a traversable wormhole?**

**Theory.** Maximal analytic extension of Schwarzschild (Kruskal-Szekeres), Flamm embedding, Morris-Thorne wormholes.

**Equations.** Time-symmetric slice $dl^2 = dr^2/(1 - 2M/r) + r^2 d\phi^2$ embeds as

$$
z(r) = 2\sqrt{2M(r - 2M)}, \qquad \left(\frac{dz}{dr}\right)^2 = \frac{2M}{r - 2M}.
$$

In Kruskal null coordinates $u = T - X$, $v = T + X$, a future-directed causal curve has $du \ge 0$ and $dv \ge 0$.

**Argument.** Region I has $v > 0$ and region III has $v < 0$. A future-directed causal curve cannot decrease $v$, so
no such curve connects them. The embedding picture is a statement about one spacelike slice, not about causal curves.

**Method.** Symbolic verification of the induced metric on the Flamm surface; numerical Kruskal slicings;
a Monte Carlo falsification test of the $u, v$ monotonicity statement (22,500 sampled pairs, zero violations);
Penrose diagrams for the eternal solution, for collapse, and an Ellis wormhole light-cone diagram.

**Findings.**

- The throat on the $T = 0$ slice has area radius $2M$. On the slice $T = 0.5$ (in units of $2M$) it is
  $1.796M$, on $T = 0.95$ it is $0.775M$, and it vanishes at $T = 1$. The bridge **pinches off**; the throat radius is
  slicing dependent while causal structure is not.
- In a collapse spacetime (a Vaidya null shell onto a flat interior) regions III and IV do not exist.
- The Ellis-Bronnikov wormhole has $dl/dt = 1$ for radial light rays, so the future of an event on one side reaches the
  other side. It requires radial null energy condition violation (see Question G).

**Limitations.** The collapse diagram is schematic for the flat interior. A fully traversable wormhole needs matter
that is not classical.

---

## Question C. Singularities

**What does geodesic incompleteness establish, and what happens to predictions when curvature becomes extremely large?**

**Theory.** Penrose (1965) and Hawking-Penrose (1970) singularity theorems; Schwarzschild at $r = 0$.

**Statement.** A spacetime is geodesically incomplete if some inextendible causal geodesic has finite affine length.
The Penrose theorem states, under the null energy condition, a non-compact Cauchy surface and a closed trapped surface,
that the spacetime is null geodesically incomplete. It does **not** by itself say that curvature diverges or what happens
at the endpoint.

**Distinction.** In Schwarzschild both hold at $r = 0$: radial timelike geodesics terminate at finite proper time and
$K = 48M^2/r^6 \to \infty$. In general the two are independent properties; examples with incomplete geodesics but bounded
curvature exist in the literature (not implemented here).

**Non-uniqueness of continuation.** Classical general relativity does not single out a continuation through $r = 0$: the
vacuum solution ends there, and the literature reports that the Schwarzschild metric admits no continuous ($C^0$) extension
through the singularity (Sbierski; to be verified before citing). Any regular replacement therefore requires additional
physics, and different replacements (such as the Hayward model of Question H) are different models, not predictions of
classical relativity.

**Outputs.** Symbolic $K$ (`tests/test_schwarzschild.py`), proper-time integrals (Question A).

---

## Question D. Quantum information and holography

**Under what assumptions can entanglement structure be related to spacetime geometry?**

**Concrete setting.** AdS$_3$/CFT$_2$ with the non-rotating BTZ black hole, $\ell = 1$.

- Bulk geometry: BTZ, $ds^2 = -(r^2 - r_+^2)dt^2 + dr^2/(r^2 - r_+^2) + r^2 d\phi^2$.
- Boundary theory: a large-$c$ CFT on a circle, $c = 3\ell/2G_N$.
- State: thermal at $T = r_+/2\pi$ (one-sided), or the thermofield double (two-sided).
- Observables: entanglement entropy of an interval; wormhole length between the two boundaries.
- Regime: classical gravity (large $c$), static slice, BTZ saddle dominant ($r_+ > 1$).
- Correspondence tested: Ryu-Takayanagi against the Calabrese-Cardy formula.

**Equations.**

$$
S_A = \frac{\operatorname{Length}(\gamma_A)}{4G_N}, \qquad
L(\Delta\phi) = 2\ln\!\left[\frac{2r_c}{r_+}\sinh\frac{r_+\Delta\phi}{2}\right], \qquad
S = \frac{c}{3}\ln\!\left[\frac{\beta}{\pi\epsilon}\sinh\frac{\pi \Delta\phi}{\beta}\right].
$$

For a thermal state the minimal surface is the minimum of two candidates, the connected geodesic and the geodesic of
the complement together with the horizon, because of the homology constraint:

$$
S_A = \frac{1}{4G_N}\min\!\big[L(\Delta\phi),\; L(2\pi - \Delta\phi) + 2\pi r_+\big].
$$

**Method.** Direct numerical integration of the geodesic equation, closed-form comparison, and an exact AdS$_3$ embedding
distance. The geodesic turning point is $r_* = r_+\coth(r_+\Delta\phi/2)$.

**Findings.** The numerical geodesic length equals the closed form to the cutoff correction. RT equals the CFT formula
(maximum difference about $1.4\times10^{-9}$ at cutoff $r_c = 10^5$). The phase transition between the two candidates occurs at
$\Delta\phi^* = 5.821$ for $r_+ = 1.5$. Beyond it $S_A - S_{A^c} = S_{BH}$ exactly (Araki-Lieb saturation). The horizon entropy
$2\pi r_+/4G_N$ equals the Cardy entropy $\pi c L T/3$.

**Limitations.** Valid only in the stated regime. The formula is not applied to arbitrary states or
dynamical spacetimes, where the covariant HRT prescription is needed. This repository does not derive AdS/CFT.

---

## Question E. Quantum branching and multiverse proposals

**How do many-worlds branching, cosmological multiverse models and additional spacetime regions differ mathematically?**

**Model.** A system qubit coupled to $N$ environment qubits by $H = \sum_k g_k Z_S Z_k$, environment qubits initially in
$\lvert +\rangle$, system in $\alpha\lvert 0\rangle + \beta\lvert 1\rangle$. The state is

$$
\lvert\psi(t)\rangle = \alpha\lvert 0\rangle\lvert E_0(t)\rangle + \beta\lvert 1\rangle\lvert E_1(t)\rangle, \qquad
\langle E_0\vert E_1\rangle = \prod_k \cos(2 g_k t),
$$

and $\rho_S$ has off-diagonal element $\alpha\beta^*\prod_k\cos(2g_kt)$. A second model with non-commuting coupling tests
generic behaviour.

**Findings.**

- The numerical QuTiP evolution matches the analytic decoherence factor to $7\times10^{-10}$.
- Born weights are unchanged by the dephasing coupling; the reduced entropy rises to the Shannon entropy of the weights
  ($0.611$ nats for weights $0.3, 0.7$).
- A finite environment recoheres; recurrences become rarer as $N$ grows.
- The non-commuting model keeps all density matrices valid (trace error $10^{-15}$, Hermiticity error $10^{-17}$,
  smallest eigenvalue non-negative to round-off) and agrees with dense matrix exponentials to $4\times10^{-10}$.

**What this does and does not show.** Unitary evolution produces components that no longer interfere. Whether these
components are "worlds" is an interpretive proposal. Nothing here makes a branch a separate spacetime.

**Classification** (see [multiverse comparison](multiverse_comparison.md)): many-worlds branching is an interpretive
proposal; the second exterior of Kruskal is a feature of the equations; eternal-inflation pocket universes and bounce
regions are predictions that depend on additional assumptions. No concrete model linking Einstein-Rosen bridges to
many-worlds branches is implemented or claimed.

---

## Question F. Time reversal

**How do time-reversal symmetry, white-hole solutions, thermodynamic irreversibility, and the reversal of an observer's experienced time differ?**

| Setting | Time-reversal operation | What is tested |
| --- | --- | --- |
| Classical trajectories | $(q, p) \to (q, -p)$, angular momentum flips sign | Velocity-Verlet round trip returns to the initial state to $10^{-13}$ |
| Wave equation | $u_t \to -u_t$ (swap time levels) | Round-trip error about $10^{-15}$; the heat equation fails with errors of order $10^{14}$ |
| Quantum evolution | $\Theta$ antiunitary, $\Theta e^{-iHt}\Theta^{-1} = e^{+i\Theta H\Theta^{-1}t}$ | Spin-1/2 in a field: $\Theta H \Theta^{-1} = -H$, so the field must be reversed; real $H$ has $\Theta = K$ as a symmetry |
| Gravity | $t \to -t$ (Schwarzschild), $T \to -T$ (Kruskal) | Exact isometry; black-hole region II maps to white-hole region IV |

**Microscopic reversibility versus irreversibility.** An area-preserving chaotic map (Arnold's cat map) with a low-entropy
initial macrostate shows coarse-grained entropy rising from about 0 to its maximum $\ln 64 = 4.159$ (4.157 reached). Exact
reversal on an integer lattice returns every particle to its initial cell; the entropy then falls back to 0. In
floating point, the reversal error grows as $e^{\lambda t}$ with $\lambda = \ln\frac{3+\sqrt5}{2} = 0.962$ and exceeds $0.4$
by $t = 40$. Irreversibility comes from special initial conditions, coarse graining, and finite precision, not from the
dynamics.

**Experienced time.** A time-reversed solution is a different solution with reversed boundary data. It does not describe an
observer who experiences time in reverse; the arrow of experienced time is tied to the thermodynamic and psychological
arrow, which the dynamics alone do not fix.

---

## Question G. Time-machine proposals

**Which explicitly defined spacetime models admit closed timelike curves, and what assumptions or limitations accompany them?**

**Model 1: wormhole with a time shift (Morris-Thorne-Yurtsever type).** Two mouths at lab separation $D$ (units $c=1$). Mouth B
travels at speed $\beta$ for lab time $T$ and returns. Its clock lags by

$$
\Delta = T\left(1 - \sqrt{1 - \beta^2}\right).
$$

Entering the lagging mouth and leaving the other advances arrival by $\Delta$ in lab time. A closed timelike curve through
the wormhole and the exterior exists iff $\Delta > D$. For $D = 1$, $\beta = 0.8$ the chronology horizon appears at lab
time $T = 2.5$.

**Model 2: Goedel universe.** Closed timelike circles exist for $r > r_c = \ln(1 + \sqrt2) = 0.8814$ because
$g_{\phi\phi} = 4a^2(\sinh^2 r - \sinh^4 r)$ becomes negative. The Ricci tensor was verified symbolically to equal the
dust form $(dt + e^x dz)^2$.

**Morris-Thorne conditions (Ellis case, $\Phi = 0$, $b = b_0^2/r$).** Verified against the symbolically computed Einstein tensor:
$\rho = b'/8\pi r^2$, $\tau = [b/r - 2(r-b)\Phi']/8\pi r^2$. At the throat $\rho - \tau = (b' - 1)/8\pi r_0^2 = -1/(4\pi b_0^2) < 0$,
so the radial null energy condition is violated.

**Limitations and physical status.** The existence of a CTC in an assumed identification model is a mathematical statement.
It requires exotic matter, and semiclassical calculations suggest that vacuum stress-energy diverges at the chronology horizon
(Hawking's chronology protection conjecture; Kim-Thorne). The Ellis wormhole is linearly unstable in the literature.
**These calculations do not produce or demonstrate a working time machine, and they do not show that any observer can
reverse their experienced time.**

---

## Question H. Quantum-gravity extensions

**Can a specified semiclassical, effective or holographic model replace a classical singularity or change the causal structure of the interior?**

**Chosen model.** The Hayward effective regular black hole,

$$
f(r) = 1 - \frac{2Mr^2}{r^3 + 2Ml^2}.
$$

It is a phenomenological metric, **not** a derived prediction of any fundamental quantum-gravity theory.

**Comparison with the Schwarzschild baseline (all verified).**

| Diagnostic | Schwarzschild | Hayward |
| --- | --- | --- |
| Curvature at $r \to 0$ | $K = 48M^2/r^6 \to \infty$ | finite, $K \to 24/l^4$ (de Sitter core) |
| Horizons | one | two for $l < 4M/(3\sqrt3) \approx 0.770M$; degenerate at that value |
| Radial geodesics ($E \ge 1$) | incomplete at $r = 0$ | reach the regular centre in finite proper time |
| Radial geodesics ($E < 1$) | incomplete at $r = 0$ | bounce inside the core and re-emerge through the inner horizon of the extended spacetime |
| Effective stress-energy | none (vacuum) | $\rho = 3Ma/[4\pi(r^3 + a)^2]$, $a = 2Ml^2$, $p_r = -\rho$ |
| NEC, WEC | satisfied | satisfied |
| SEC | satisfied | violated for $r < (Ml^2)^{1/3}$ |
| DEC | satisfied | violated for $r > (4Ml^2)^{1/3}$ |
| Classical limit | exact | $l \to 0$ reproduces Schwarzschild ($K$, $f$, $T_H \to 1/8\pi M$) |
| Hawking temperature | $1/8\pi M$ | goes to zero at extremality |
| Exterior stability diagnostic | $\ell=2$ scalar mode decays | decays; frequency shifts (for $l = 0.5M$: $M\omega = 0.4936 - 0.0923i$ versus $0.4836 - 0.0968i$) |

**Conclusions the computation supports.** In this effective model the singularity is replaced by a de Sitter core and the
maximal radial geodesics are complete in the sense of reaching a regular centre or a turning point.
The causal structure of the extended spacetime changes: there is an inner horizon, hence a Reissner-Nordstrom-like chain of regions.

**Conclusions the computation does not support.** Finite curvature does not prove that the model is a complete or correct
physical resolution. The model violates the strong and dominant energy conditions in different regions, which is
precisely what its geometry requires. The inner horizon of regular black holes is believed to be unstable
(mass-inflation arguments); this repository does not analyse inner-horizon stability, so the extended causal structure
should not be taken as physically realized. The value $l = 0.5M$ is chosen for visibility; a Planck-scale $l$ for an astrophysical black hole gives corrections far below any observation.

**Other candidate directions** (semiclassical dilaton gravity, loop-quantum-gravity-inspired interiors, holographic models)
were considered but not implemented.

---

## Overall assessment

The classical results (Questions A to C and the baseline of H) are reproduced with symbolic and numerical cross-checks.
The holographic and quantum-information results (D, E) are verified inside explicitly stated toy and large-$c$ regimes.
Questions F and G are answered at the level of mathematical definitions and explicit models. The repository does not
establish a new physical result about spacetime and quantum mechanics; it provides a verified computational base on which
such investigations can be built. The limitations of each part are collected in
[limitations_and_open_problems.md](limitations_and_open_problems.md).
