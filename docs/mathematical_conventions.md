# Mathematical Conventions

Every numerical result in this repository depends only on the conventions stated here. The same conventions are
encoded in `src/beyond_spacetime/conventions.py` and `src/beyond_spacetime/symbolic.py`.

## 1. Spacetime conventions

| Item | Convention |
| --- | --- |
| Metric signature | $(-,+,+,+)$ |
| Units | Geometrized, $G = c = 1$, unless a function states otherwise. Mass $M$ is the length $GM/c^2$. |
| Christoffel symbols | $\Gamma^a{}_{bc} = \tfrac12 g^{ad}(\partial_c g_{db} + \partial_b g_{dc} - \partial_d g_{bc})$ |
| Riemann tensor | $R^a{}_{bcd} = \partial_c \Gamma^a{}_{db} - \partial_d \Gamma^a{}_{cb} + \Gamma^a{}_{ce}\Gamma^e{}_{db} - \Gamma^a{}_{de}\Gamma^e{}_{cb}$ |
| Ricci tensor and scalar | $R_{bd} = R^a{}_{bad}$, $R = g^{bd}R_{bd}$ |
| Einstein equations | $G_{\mu\nu} + \Lambda g_{\mu\nu} = 8\pi T_{\mu\nu}$ |
| Kretschmann scalar | $K = R_{abcd}R^{abcd}$, evaluated as $R^{ab}{}_{cd}R^{cd}{}_{ab}$ |
| Geodesic deviation | $\dfrac{D^2\xi^\mu}{d\tau^2} = -R^\mu{}_{\nu\alpha\beta}\,u^\nu \xi^\alpha u^\beta$ |
| Tidal tensor | $\mathcal{E}_{ij} = R_{i\hat 0 j \hat 0}$ in an orthonormal frame with $e_{\hat 0} = u$; then $\ddot\xi^i = -\mathcal{E}_{ij}\xi^j$ |

These are the Misner-Thorne-Wheeler conventions. A different Riemann sign flips the sign of $\mathcal{E}_{ij}$ and of the
geodesic-deviation equation together; the physical statement (radial stretching, transverse compression) does not change.

## 2. Coordinates

| Name | Definition |
| --- | --- |
| Schwarzschild | $(t, r, \theta, \phi)$, $f(r) = 1 - 2M/r$ |
| Tortoise coordinate | $r_* = r + 2M\ln\lvert r/2M - 1\rvert$ |
| Ingoing Eddington-Finkelstein | $v = t + r_*$, $ds^2 = -f\,dv^2 + 2\,dv\,dr + r^2 d\Omega^2$ |
| Outgoing Eddington-Finkelstein | $u = t - r_*$, $ds^2 = -f\,du^2 - 2\,du\,dr + r^2 d\Omega^2$ |
| Kruskal-Szekeres (region I) | $X = \sqrt{r/2M - 1}\,e^{r/4M}\cosh(t/4M)$, $T = \sqrt{r/2M - 1}\,e^{r/4M}\sinh(t/4M)$ |
| Kruskal-Szekeres (region II) | $T = \sqrt{1 - r/2M}\,e^{r/4M}\cosh(t/4M)$, $X = \sqrt{1 - r/2M}\,e^{r/4M}\sinh(t/4M)$ |
| Kruskal metric | $ds^2 = -\dfrac{32M^3}{r}e^{-r/2M}(dT^2 - dX^2) + r^2 d\Omega^2$ |
| Implicit radius | $\left(\dfrac{r}{2M} - 1\right)e^{r/2M} = X^2 - T^2$, inverted with the Lambert $W$ function |
| Kruskal null coordinates | $U = T - X$, $V = T + X$ |
| Penrose coordinates | $\tilde u = \arctan U$, $\tilde v = \arctan V$, $T_p = \tilde v + \tilde u$, $X_p = \tilde v - \tilde u$ |

Regions of the maximally extended solution: I ($X > \lvert T\rvert$), II ($T > \lvert X\rvert$), III ($X < -\lvert T\rvert$),
IV ($T < -\lvert X\rvert$). Singularities are $T^2 - X^2 = 1$. In Penrose coordinates both singularities are the
straight lines $T_p = \pm\pi/2$ because $\tan\tilde u\tan\tilde v = 1$ implies $\tilde u + \tilde v = \pi/2$.

Coordinate components (for example $g_{tt}$, $\Gamma^r{}_{tt}$) are never plotted as if they were physical. Invariants
($K$, proper time, tidal eigenvalues, areal radius) are used for physical statements.

## 3. Index and summation rules

Greek indices run over spacetime, Latin hatted indices over an orthonormal frame, repeated indices are summed,
and indices are raised and lowered with $g_{\mu\nu}$ of the stated signature. Frame vectors are stored with upper
(contravariant) components in coordinate bases.

## 4. Complex frequencies and Fourier conventions

Time dependence is $e^{-i\omega t}$. Quasinormal modes have $\operatorname{Im}\omega < 0$ (decay). Leaver's continued
fraction is implemented in units $2M = 1$ and the returned values are converted to $M\omega$. The perturbation equation is

$$
\frac{\partial^2\psi}{\partial t^2} - \frac{\partial^2\psi}{\partial r_*^2} + V(r)\,\psi = 0,
\qquad
V_s(r) = f\left[\frac{\ell(\ell+1)}{r^2} + (1 - s^2)\frac{2M}{r^3}\right],
$$

with spin weight $s = 0$ (scalar) or $s = 2$ (Regge-Wheeler, gravitational axial).

## 5. Quantum conventions

| Item | Convention |
| --- | --- |
| State normalization | $\langle\psi\vert\psi\rangle = 1$; density matrices satisfy $\operatorname{Tr}\rho = 1$, $\rho = \rho^\dagger$, $\rho \ge 0$ |
| Pauli matrices | $X = \sigma_x$, $Y = \sigma_y$, $Z = \sigma_z$ with eigenvalues $\pm 1$ |
| Von Neumann entropy | $S(\rho) = -\operatorname{Tr}\rho\ln\rho$, natural logarithm (nats) |
| Reduced state | $\rho_A = \operatorname{Tr}_B\lvert\psi\rangle\langle\psi\rvert$ |
| Coherence measures | $C_{\ell_1}(\rho) = \sum_{i\ne j}\lvert\rho_{ij}\rvert$; $C_{\mathrm{rel}}(\rho) = S(\Delta\rho) - S(\rho)$ with $\Delta$ the dephasing map in the computational basis |
| Time reversal | $\Theta$ antiunitary; for spin-1/2, $\Theta = i\sigma_y K$ with $K$ complex conjugation |
| Thermofield double | $\lvert\mathrm{TFD}(\beta)\rangle = Z^{-1/2}\sum_n e^{-\beta E_n/2}\lvert n\rangle_L\lvert n^*\rangle_R$ |

## 6. Holographic conventions

AdS radius $\ell = 1$. The non-rotating BTZ metric is

$$
ds^2 = -(r^2 - r_+^2)\,dt^2 + \frac{dr^2}{r^2 - r_+^2} + r^2 d\phi^2, \qquad \phi \sim \phi + 2\pi,
$$

with Hawking temperature $T = r_+/2\pi$, Bekenstein-Hawking entropy $S_{BH} = 2\pi r_+/4G_N$, and Brown-Henneaux central
charge $c = 3\ell/2G_N$. The Ryu-Takayanagi relation $S_A = \operatorname{Length}(\gamma_A)/4G_N$ is used only for
static slices of a classical bulk in the large-$c$ regime. The BTZ saddle dominates over thermal AdS$_3$ for $r_+ > 1$
(Hawking-Page); statements about a physical thermal state therefore use $r_+ \ge 1$. The cutoff surface is
$r = r_c$ and corresponds to the CFT lattice cutoff $\epsilon = 1/r_c$.

## 7. Entropy and information conventions

All entropies are in nats. Cutoff-dependent quantities are compared only between configurations with the same
cutoff, or after subtracting a reference (vacuum) value.

## 8. Numerical precision requirements

| Quantity | Requirement used in tests |
| --- | --- |
| Symbolic identities (vacuum equations, $K = 48M^2/r^6$) | Exactly zero after simplification |
| Closed-form vs numerical (proper times, geodesic lengths) | Relative error below $10^{-8}$ unless a cutoff correction is stated |
| ODE integration | DOP853, `rtol = 1e-12`, `atol = 1e-14` for geodesics |
| Density matrices | Trace error, Hermiticity error and negative eigenvalue below $10^{-10}$ |
| QNM extraction | Leaver continued fraction to $10^{-13}$; time-domain agreement to about $10^{-3}$ |
| Floating-point reversal experiments | Failure is expected and reported, not hidden (see the cat-map experiment) |

## 9. Dimensional-consistency rules

In geometrized units a mass, a length and a time all have dimension of length. Results that are converted to SI
(tidal accelerations, spaghettification radii) pass through `conventions.py`, which restores $G$ and $c$. The
tidal acceleration across a body of length $\xi$ at areal radius $r$ is
$a = c^2\,(2 G M/c^2)\,\xi / r^3$ in SI units. At the horizon $r = 2GM/c^2$ this equals
$c^6\xi/(4G^2M^2)$, which is used as an independent test.
