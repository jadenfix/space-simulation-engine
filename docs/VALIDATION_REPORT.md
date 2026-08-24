# Spacewind Native Validation Report

**Engine:** 0.2.0
**Validation date (America/Los_Angeles):** 2026-08-16
**Scope:** bounded software and numerical checks; not propulsion validation

## Build and analysis matrix

| Layer | Result |
|---|---:|
| Deterministic checks | 3,912 pass |
| GCC warning-as-error C11 | pass |
| Clang warning-as-error C11 | pass |
| CMake + CTest | pass |
| AddressSanitizer | pass |
| UndefinedBehaviorSanitizer | pass |
| Clang static analyzer | 13 translation units pass |
| GCC `-fanalyzer` | 13 translation units pass |
| Source-line coverage | 86.26% of 2467 executable lines |
| Branch sites executed | 93.95% of 1255 |
| Branch outcomes taken | 67.57% of 1255 |

## Numerical receipts

- Kerr circular geodesic: relative radius drift `1.971e-16`, interval drift `-1.443e-14` over one orbit.
- Schwarzschild circular geodesic: relative radius drift `1.794e-14`.
- Sun-Earth-Jupiter N-body: 12-year relative energy error `4.391e-14` at the selected six-hour step.
- 41^3 Poisson sphere: residual `7.928e-09` after 180 iterations; maximum sampled potential error `1.516%`.
- Boris magnetic rotation: relative speed drift `4.366e-16`.
- Periodic FDTD pulse: relative field-energy change `9.824e-04` in the finite run.
- Electrostatic PIC: relative total-energy change `1.720e-09` in the finite run.

## Mission scenarios

| Scenario | Integrator | Accepted | Rejected trials | Final radius (AU) | Final speed (km/s) |
|---|---|---:|---:|---:|---:|
| `baseline` | rk4 | 43,200 | 0 | 1.007777959 | 29.691201221 |
| `dynamic_soaring_shear` | rk4 | 120,960 | 0 | 1.000551361 | 2.256884934 |
| `cme_encounter` | dopri54 | 2,049 | 905 | 1.000748282 | 29.804228817 |
| `magnetic_sail` | rk4 | 43,200 | 0 | 1.011519301 | 29.640160039 |
| `hybrid_glider` | dopri54 | 1,385 | 397 | 1.000343707 | 29.812073971 |

The dynamic-soaring and hybrid rows are conditional model outputs. Their force coefficients and effective interaction radii are explicit inputs, not measured spacecraft performance.

## Reproducibility boundary

Every mission receipt records the source-config FNV-1a hash, seed, compiler identity, C standard, double precision parameters, floating-point evaluation method, fast-math status, accepted steps, and internally rejected adaptive trials. Adaptive trajectories may differ slightly across compiler and optimization profiles even when physical inputs are identical; the receipt exposes that context instead of asserting bitwise portability.

## Nonclaims

- This is not a solver for the full nonlinear Einstein field equations.
- Mission-scale electric- and magnetic-sail coefficients are phenomenological until calibrated against multidimensional kinetic calculations or experiment.
- The planar shear and CME shell are synthetic controlled environments, not live space-weather products.
- The local PIC, FDTD, and ideal-MHD solvers are one-dimensional and cannot establish three-dimensional tether or plasma-wing performance.
- A successful synthetic dynamic-soaring trajectory does not prove a realizable vehicle or net energy extraction in the heliosphere.
