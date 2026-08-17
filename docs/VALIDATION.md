# Validation, Claims, and Nonclaims

## Validation principle

The test suite verifies bounded numerical statements. It does not certify a
physical propulsion concept.

Each result belongs to one of four classes:

- **Equation implementation:** a declared equation is represented in code.
- **Finite check:** one bounded numerical experiment passes a tolerance.
- **Model result:** a trajectory follows from selected parameters.
- **Open problem:** measurement or higher-fidelity work is still required.

## Build matrix exercised

| Check | Result |
|---|---|
| GCC, strict C11 warnings | pass |
| Clang, strict C11 warnings | pass |
| Make build and test | pass |
| CMake build and CTest | pass |
| AddressSanitizer | pass |
| UndefinedBehaviorSanitizer | pass |
| Clang static analyzer, 13 translation units | pass |
| Combined test/example source-line coverage | about 86.3% |
| Deterministic unit/integration checks | 3,912 pass |

The machine-readable and human-readable recorded run are in
[`VALIDATION_REPORT.json`](VALIDATION_REPORT.json) and
[`VALIDATION_REPORT.md`](VALIDATION_REPORT.md).

The exact count can change when tests are added. CI should treat any failure as
a regression but should not interpret a passing workflow as new physical
evidence.

## Current finite checks

### Linear algebra and geometry

- Cross products are orthogonal to both inputs.
- Quaternion normalization is stable in the mission loop.
- A representative 4x4 matrix inverse multiplies back to identity.
- Minkowski Christoffel symbols are numerically zero.
- A Minkowski geodesic remains linear.
- The metric interval remains at the expected value in that finite case.
- Kerr with zero spin agrees with Schwarzschild at a sampled point.
- A prograde equatorial circular Kerr geodesic at `r = 10M`, `a/M = 0.8`
  closes one finite orbit while preserving radius and the selected constants to
  floating-point tolerance.
- FLRW with unit scale factor and zero Hubble rate agrees with Minkowski.

### Gravity

- Velocity Verlet maintains a one-AU circular test orbit within its declared
  energy/radius bounds.
- The barycentric two-body leapfrog test maintains energy, momentum, center of
  mass, and separation within declared bounds over one nominal year.
- The 12-year Sun-Earth-Jupiter example reports relative energy drift of about
  `4.4e-14` for its six-hour step. This unusually small value is conditional on
  the smooth idealized initial state and should not be generalized to close
  encounters.
- A manufactured sine-mode Poisson problem converges below `1e-9` relative
  residual on a 17^3 grid and has below `0.4%` discrete L2 error.
- A 41^3 uniform-sphere demonstration converges in roughly 180 iterations at
  SOR relaxation 1.85. Worst sampled potential error is about 1.6%, dominated by
  the coarse stair-step representation of the sphere.

### Solar wind

- The Parker branch accelerates from 0.1 AU through 1 AU to 5 AU for the default
  isothermal parameters.
- Density falls while conserving the declared radial mass flux.
- The one-AU reference density is reproduced.
- Debye length and magnetic/electric field magnitudes are finite and positive.

These checks do not validate the Parker model against a particular solar cycle
or event.

### Particle and field solvers

- The nonrelativistic Boris pusher conserves speed in a uniform magnetic field
  with zero electric field over 10,000 finite steps.
- The relativistic pusher remains subluminal and conserves the selected speed in
  the same magnetic-only case.
- PIC kinetic and field energies remain finite in the finite periodic run.
- The FDTD pulse's total field energy changes by less than the declared 3%
  bound in the selected periodic grid/time step.

### MHD

- The Brio-Wu run keeps density and pressure positive.
- Total mass is conserved within `1e-10` relative error in the selected run.
- The test does not compare every wave location against a high-resolution
  reference solution yet.

### Integrators and mission coupling

- RK4 integrates a harmonic oscillator to a strict analytic tolerance.
- Dormand-Prince rejects an intentionally oversized trial step and returns an
  accepted smaller step.
- Mission state and quaternion remain finite/normalized over the finite coupled
  run.
- Radiation and electric-sail force paths are exercised.
- Adaptive receipts report internal rejected trial steps.
- The all-force `hybrid_glider.cfg` integration path completes with every
  mission-level force switch enabled. Its output is an integration check only.

## Scenario results are conditional

The `hybrid_glider.cfg` configuration is a software stress case, not a
hardware design. Its charge, magnetic moment, tether, area, and force
coefficients are explicit inputs.

The included dynamic-soaring configuration is deliberately synthetic. Its
planar shear, force coefficients, maximum sheath radius, and controller are
explicit inputs. A speed increase or closed crossing pattern in that file is a
**model result** only.

It does not imply:

- that the specified lift coefficient exists in solar-wind plasma;
- that a compact device can generate the requested virtual wing;
- that voltage/current/power/mass/deployment constraints close;
- that the same shear remains spatially coherent;
- that the controller is stable under real turbulence and sensor delay;
- that the vehicle can reach the modeled boundary.

## Highest-value missing validation

1. Measure transverse and axial force on a candidate electric/magnetic/plasma
   wing in a calibrated plasma chamber.
2. Reproduce the setup in a multidimensional kinetic code.
3. Build a response surface with uncertainty, not one best-fit coefficient.
4. Propagate that uncertainty through the mission simulator.
5. Demonstrate deployment and high-voltage operation in orbit.
6. Detect and recross a real plasma boundary while measuring net mechanical
   energy gain over a complete cycle.

## Convergence tests still needed

- Mission trajectory vs time-step sequence.
- Geodesic invariant error vs affine step and metric derivative step.
- Poisson solution vs grid spacing and domain size.
- PIC dispersion/noise vs particles per cell and grid size.
- FDTD phase velocity vs Courant number.
- MHD shock locations vs grid resolution and higher-order reference.
- Field-sail response uncertainty vs local plasma resolution.

## Nonclaims

Version 0.2.0 does not claim:

- a solution of the full Einstein equations;
- a complete heliosphere forecast;
- an experimentally validated electric- or magnetic-sail force law;
- a proven dynamic-soaring vehicle;
- supersonic, relativistic, or interstellar mission feasibility;
- flight-qualified deployment, charging, thermal, radiation, or structural
  hardware;
- safety or mission certification.
