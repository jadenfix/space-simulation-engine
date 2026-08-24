# Space Simulation Engine

[![native-spacewind](https://github.com/jadenfix/space-simulation-engine/actions/workflows/test.yml/badge.svg)](https://github.com/jadenfix/space-simulation-engine/actions/workflows/test.yml)

**Spacewind** is the native physics engine in this repository: a dependency-light
C11 simulation environment for testing whether a deployable
spacecraft can exchange momentum with photons, solar-wind plasma, magnetic
fields, velocity shear, and gravity.

The project is intentionally built from low-level numerical pieces rather than
wrapping a large astrodynamics framework. Its purpose is to make every
assumption inspectable, every force separable, and every result falsifiable.
It does **not** assume that spacetime is a material. It models spacetime geometry,
matter-generated weak gravitational fields, the heliospheric environment, and
spacecraft coupling as distinct physical layers.

## What is implemented

### Geometry and gravity

- Minkowski spacetime in Cartesian coordinates.
- Schwarzschild spacetime in spherical coordinates.
- Kerr spacetime in Boyer-Lindquist coordinates.
- Flat FLRW/de Sitter expansion.
- Numerical metric derivatives and Christoffel symbols.
- Timelike/null geodesic integration with RK4.
- Weak-field proper-time accumulation.
- Schwarzschild first post-Newtonian orbital correction.
- Point-mass gravity and optional J2 acceleration.
- Direct-summation barycentric N-body gravity with symplectic leapfrog.
- A three-dimensional finite-difference Poisson solver that maps mass density
  to gravitational potential, acceleration, and a weak-field metric.

### Solar and plasma environment

- Isothermal transonic Parker solar wind.
- Mass-flux-conserving proton density.
- Parker-spiral magnetic field.
- Ideal-MHD motional electric field, `E = -v x B`.
- Photon flux, plasma thermal pressure, dynamic pressure, magnetic pressure,
  Debye length, sound speed, Alfvén speed, and fast magnetosonic speed.
- Rotating fast/slow stream modulation.
- Deterministic Fourier-mode magnetic and Alfvénic velocity turbulence.
- Configurable planar velocity shear for controlled dynamic-soaring tests.
- Configurable propagating CME-like density/field/speed shell.

### Kinetic and field solvers

- Nonrelativistic Boris charged-particle pusher.
- Relativistic Boris charged-particle pusher.
- One-dimensional electrostatic particle-in-cell solver:
  cloud-in-cell deposition, neutralizing background, spectral Poisson solve,
  field interpolation, and kick-drift-kick stepping.
- One-dimensional Yee-style FDTD Maxwell solver.
- One-dimensional finite-volume ideal-MHD solver with HLL fluxes, CFL time
  stepping, positivity floors, and the Brio-Wu shock-tube initial condition.

### Spacecraft and control

- Six-degree-of-freedom position, velocity, quaternion attitude, and rigid-body
  angular dynamics.
- Solar-radiation pressure with absorption and specular reflection.
- Parameterized electric-sail interaction area derived from Debye length,
  electron temperature, tether voltage, and tether length.
- Parameterized magnetic-sail pressure-balance radius.
- Charged-spacecraft Lorentz force.
- Force and torque decomposition in every mission trace.
- Fixed, Sun-cone, and shear-crossing dynamic-soaring controllers.
- Fixed-step RK4 and adaptive Dormand-Prince 5(4) integration.
- Reproducible CSV traces and JSON receipts with the source-config hash,
  accepted steps, rejected adaptive trials, seed, enabled models, and nonclaim.

## Build

Requirements:

- A C11 compiler (`gcc`, `clang`, or compatible)
- `make` or CMake 3.18+
- The standard C library and `libm`

No runtime Python dependency is required.

```sh
make
make test
make examples
```

The strict default flags are:

```text
-std=c11 -O2 -g -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wstrict-prototypes
```

CMake works independently:

```sh
cmake -S . -B cmake-build -DCMAKE_BUILD_TYPE=Release
cmake --build cmake-build
ctest --test-dir cmake-build --output-on-failure
```

Sanitizer validation:

```sh
make asan
make ubsan
```

Run the complete locally available validation matrix or create source archives:

```sh
make validate-all
make package
```

## Command-line interface

```text
spacewind simulate <config.cfg> <trajectory.csv> [receipt.json]
spacewind wind <profile.csv>
spacewind geodesic <orbit.csv>
spacewind kerr <orbit.csv>
spacewind pic <plasma.csv>
spacewind fdtd <wave.csv>
spacewind mhd <shock.csv>
spacewind boris <particle.csv>
spacewind nbody <orbits.csv>
spacewind poisson <gravity.csv>
```

Examples:

```sh
./build/spacewind simulate \
  configs/baseline_1au.cfg \
  output/baseline.csv \
  output/baseline.receipt.json

./build/spacewind simulate \
  configs/dynamic_soaring_shear.cfg \
  output/dynamic_soaring.csv \
  output/dynamic_soaring.receipt.json

./build/spacewind simulate \
  configs/hybrid_glider.cfg \
  output/hybrid_glider.csv \
  output/hybrid_glider.receipt.json

./build/spacewind poisson output/gravity_grid.csv
./build/spacewind nbody output/nbody.csv
```

`make examples` runs all solver demonstrations from one command.

## Fidelity hierarchy

Spacewind does not pretend that one numerical method can resolve the entire
problem from particle scales to interplanetary distances. It separates the
problem into three levels:

1. **Mission scale:** analytic heliosphere fields and spacecraft ODEs over days,
   months, and astronomical units.
2. **Mesoscale:** ideal-MHD discontinuities and stream interfaces.
3. **Kinetic/local scale:** charged particles, electrostatic PIC, and Maxwell
   field propagation.

The intended calibration loop is:

```text
laboratory or kinetic simulation
        -> force/torque response surface
        -> mission-scale field-sail coefficients
        -> trajectory and control simulation
        -> new high-value kinetic or chamber experiment
```

This prevents an impossible full-PIC simulation of an entire interplanetary
mission while keeping the uncertain coupling coefficients explicit.

## Configuration

Configuration files use strict `key = value` syntax. All physical values are SI
unless the key explicitly says otherwise. Unknown keys are errors rather than
being silently ignored.

A minimal mission configuration:

```ini
integrator = rk4
controller = sun_cone
duration_s = 2592000
step_s = 60
output_every_steps = 60

position_m = 149597870700, 0, 0
velocity_m_s = 0, 29784.6918, 0
mass_kg = 12
sail_area_m2 = 100
sun_cone_angle_deg = 25

tether_length_m = 20000
tether_voltage_v = 20000
electric_sail_cd = 1.0
electric_sail_cl = 0.15
electric_sail_max_radius_m = 200

enable_gravity = true
enable_1pn = true
enable_radiation_pressure = true
enable_electric_sail = true
enable_magnetic_sail = false
```

Included configurations:

- `baseline_1au.cfg`: heliocentric photon/electric-sail baseline.
- `dynamic_soaring_shear.cfg`: explicit synthetic shear-layer falsification
  experiment.
- `cme_encounter.cfg`: adaptive integration through a synthetic CME shell.
- `magnetic_sail.cfg`: parameterized magnetic-sail mission case.
- `hybrid_glider.cfg`: all mission-level force paths enabled together as a
  software-integration stress case, not a realizability claim.

## Reproducibility

Every loaded config is hashed with FNV-1a 64-bit and included in the receipt.
The deterministic turbulence generator is seeded. Output values use 17 decimal
digits so double-precision state can be independently replayed within ordinary
floating-point/compiler limitations. Adaptive trajectories can still differ by
a few accepted/rejected steps across compiler and optimization profiles, so the
receipt records the build-level floating-point context rather than pretending
bitwise portability.

A receipt records:

- engine version, compiler identity, C standard, and floating-point model;
- config hash and random seed;
- selected integrator and controller;
- accepted steps and rejected adaptive trial steps;
- final position, velocity, proper time, radius, and speed;
- every enabled force/environment model;
- an explicit warning that field-sail response parameters are not flight
  validation.

## Current validation snapshot

At version 0.2.0:

- 3,912 deterministic checks pass.
- GCC strict build passes.
- Clang strict build passes.
- CMake/Ninja/CTest passes.
- Clang static analysis passes for all 13 translation units.
- The combined test and example matrix executes about 86.3% of native source
  lines; this is coverage evidence, not a correctness proof.
- AddressSanitizer and UndefinedBehaviorSanitizer pass.
- A two-body Sun-Earth leapfrog run conserves energy to the test threshold over
  one year.
- The 12-year Sun-Earth-Jupiter demonstration has relative energy drift around
  `4.4e-14` for the selected six-hour step and initial condition.
- The manufactured three-dimensional Poisson problem converges below `1e-9`
  relative residual and stays below `0.4%` discrete L2 error on a 17^3 grid.
- The 41^3 uniform-sphere demonstration converges in roughly 180 SOR iterations
  to a residual near `8e-9`; its coarse-grid potential differs from the analytic
  sphere by about 1.6% at the worst sampled point.
- Boris magnetic rotation conserves particle speed in the no-electric-field
  test.
- The FDTD pulse conserves field energy within the declared finite-grid bound.
- The MHD shock tube maintains positive density/pressure and conserves total
  mass within the declared bound.
- A prograde equatorial Kerr circular-geodesic demonstration closes one orbit
  at `r = 10M` for spin `a/M = 0.8` while preserving the selected invariants to
  floating-point tolerance.
- The hybrid mission scenario exercises radiation, electric-sail, magnetic-sail,
  Lorentz, Newtonian-gravity, and 1PN paths in one adaptive run.

See [`docs/VALIDATION.md`](docs/VALIDATION.md) for the validation policy and
[`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md) for the recorded build,
coverage, solver, and mission receipts.

## What this does not establish

Spacewind is a research environment, not evidence that the proposed spacecraft
can yet be built.

In particular:

- It does not treat spacetime as a fluid or material surface.
- It does not solve the full nonlinear Einstein field equations.
- The Poisson lattice is a weak-field approximation.
- The mission environment is not a live WSA-ENLIL or OMNI data product yet.
- The CME and shear layers are controlled synthetic models.
- Electric- and magnetic-sail lift coefficients are phenomenological inputs.
- A positive dynamic-soaring trajectory in the synthetic layer is a model
  result, not a demonstrated plasma-wing lift-to-drag ratio.
- One-dimensional PIC/FDTD/MHD solvers cannot resolve three-dimensional tether,
  sheath, reconnection, turbulence, or deployment physics.
- Positivity floors in the first-order MHD solver trade exact total-energy
  conservation for robustness near nonphysical states.

The decisive physical experiment remains a repeatable, commandable transverse
force measurement on a candidate field-sail device in a calibrated plasma flow.

## Project layout

```text
include/spacewind/   public C API
src/                 implementation and CLI
tests/               deterministic finite checks
configs/             mission scenarios
docs/                equations, architecture, validation, roadmap
claims/              machine-readable claims and nonclaims
scripts/             release, analyzer, sanitizer, and packaging helpers
output/               generated traces; ignored except .gitkeep
```

## Research discipline

Every major statement should be treated as one of:

- **Equation implementation:** code implementing a stated physical/numerical
  model.
- **Finite check:** a bounded numerical test with declared grid, step, and
  tolerance.
- **Model result:** a result conditional on a configuration and calibration.
- **Open problem:** a quantity requiring higher-fidelity simulation or physical
  measurement.

A passing test demonstrates that code behaves as declared on that finite case.
It does not upgrade a field-sail model into an experimentally established
propulsion system.
