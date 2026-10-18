# AI Assistance and Provenance

This repository was produced in a single working session with an AI coding assistant (Claude, Anthropic) acting on
a research brief written by the project owner. The brief required that AI-generated material never be accepted as a
derivation or a result without independent verification. This file records how that requirement was applied.

## What the assistant did

- Wrote the Python package, the experiment scripts, the tests and the Markdown documentation.
- Ran the tests and experiments and inspected the resulting figures.

## How results were verified

| Claim type | Independent check |
| --- | --- |
| Curvature identities | Symbolic computation from the metric (SymPy), no closed form assumed |
| Proper times, tidal tensors | Closed-form analytic expressions, plus an independent full-Riemann evaluation and an exact neighbouring-geodesic experiment |
| Quasinormal modes | Continued-fraction solver compared with published literature values and with a separate time-domain evolution |
| Holographic lengths | Numerical geodesic integral, closed form, and exact AdS$_3$ embedding-space distance |
| Quantum models | QuTiP results compared with dense matrix exponentials and with exact analytic decoherence factors |
| Morris-Thorne equations | Einstein tensor computed from the metric and compared with the closed-form density, tension and pressure |

## Limits of this verification

- The published QNM values used for comparison were recalled from the literature and were not re-fetched during the
  build.
- The reference list in [references.md](references.md) was written from memory without live lookup. Bibliographic
  details (volume, page, year) must be checked before the list is cited in a paper.
- No peer review of the physics has taken place. The repository is a computational study, not a published result.

No credentials, private data or confidential material were involved.
