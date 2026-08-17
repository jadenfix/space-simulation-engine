# Roadmap

## Phase 1: data-grounded heliosphere

- Add a versioned NASA OMNI importer for near-Earth magnetic field, plasma,
  energetic-particle, and geomagnetic-index data.
- Add WSA-ENLIL/CCMC snapshot ingestion and interpolation.
- Add SPICE kernels for ephemerides, frames, occultations, and spacecraft
  geometry.
- Preserve the analytic Parker model as the controlled baseline.
- Add units/schema validation and input-file hashes to every imported dataset.

## Phase 2: calibrated plasma-wing response

- Add two-dimensional electrostatic/electromagnetic PIC.
- Add conducting tether boundary conditions and spacecraft charging.
- Add electron emission/current closure and finite tether geometry.
- Add magnetic dipole/plasma-magnet source models.
- Sweep density, temperature, Mach number, magnetic-field angle, voltage,
  current, geometry, and incidence.
- Fit a physics-constrained force/torque surrogate with uncertainty.
- Replace free mission-scale lift/drag coefficients with that calibrated map.

## Phase 3: higher-fidelity heliospheric structures

- Add multidimensional finite-volume MHD with constrained transport.
- Add second-order reconstruction and HLLD/Roe-type solvers.
- Add adaptive mesh refinement around shocks and shear layers.
- Couple kinetic patches to MHD where sheath/reconnection physics matters.
- Validate against published shock tubes and NASA event data.

## Phase 4: spacecraft realism

- Flexible tether and membrane finite elements.
- Deployment transients, spin stabilization, and entanglement/failure modes.
- High-voltage power, electron/ion gun current, thermal, radiation, and
  micrometeoroid budgets.
- Sensor models, latency, estimator uncertainty, actuator saturation, and
  fault-tolerant control.
- Hardware-in-the-loop interface and plasma-chamber calibration files.

## Phase 5: trajectory search and falsification

- Multiple shooting and direct collocation.
- Robust/model-predictive control over uncertain plasma boundaries.
- Bayesian uncertainty propagation and adversarial parameter sweeps.
- Compare delivered mass, time, power, and risk against photon sails, Hall
  thrusters, gravity assists, and conventional chemical stages.
- Define automatic go/no-go gates for minimum measured transverse force,
  lift-to-drag ratio, full-system acceleration, power-to-mass, and deployment
  reliability.

## Phase 6: genuine dynamical spacetime, only where justified

The current metric and Poisson paths are enough for solar-system mission design.
If the research expands to compact objects or dynamical strong gravity, add a
separate numerical-relativity subsystem rather than mislabeling the weak-field
lattice:

- 3+1 ADM/BSSN variables;
- constraint-preserving evolution;
- gauge conditions;
- finite-difference or discontinuous-Galerkin evolution;
- black-hole and gravitational-wave benchmarks;
- relativistic MHD/PIC coupling.

This is a distinct major program and should not be presented as necessary for a
solar-wind glider near ordinary planets and stars.

## Performance roadmap

- Structure-of-arrays particle storage.
- SIMD Boris/PIC kernels.
- Multigrid Poisson.
- Barnes-Hut/FMM gravity.
- OpenMP deterministic mode.
- MPI domain decomposition.
- CUDA/HIP backends.
- Binary chunked traces with checksums.
