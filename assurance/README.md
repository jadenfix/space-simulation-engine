# Spacewind assurance kernel

This directory is an independent C11 verification layer for the Spacewind
simulation engine. It is intentionally small, dependency-light, and separable
from the mission and plasma solvers it audits.

Its purpose is not to produce a more impressive trajectory. Its purpose is to
make incorrect trajectories, unclosed ledgers, unsupported extrapolations,
biased force experiments, and over-promoted scientific claims fail closed.

## Implemented assurance surfaces

### Numerical foundations

- Seven-base-dimension runtime SI algebra.
- Rejection of invalid mixed-dimension addition.
- Outward-rounded interval arithmetic.
- Kahan compensated summation.
- Welford online statistics.
- Deterministic Halton design sequences and transcript hashes.
- Three-level Richardson extrapolation and grid-convergence index checks.

### Conservation and invariants

- Mass, charge, momentum, and energy control-volume ledgers.
- Minkowski four-velocity normalization.
- Lorentz-force decomposition and the zero-work magnetic-force identity.
- Electromagnetic field energy, Poynting flux, and the two field invariants.
- Parker-wind mass-flux, radial magnetic-flux, and motional-field identities.
- Photon momentum and plasma momentum-flux upper bounds.
- Kepler energy, angular momentum, eccentricity, and apsis invariants.
- N-body center-of-mass, momentum, angular momentum, and total-energy audits.

### Cross-environment validity

The plasma similarity contract compares twelve dimensionless or scale-sensitive
coordinates:

1. normalized electrode voltage;
2. electrode-radius to Debye-length ratio;
3. ion-gyroradius to interaction-length ratio;
4. ion-inertial-length to interaction-length ratio;
5. sonic Mach number;
6. Alfvénic Mach number;
7. fast-magnetosonic Mach number;
8. plasma beta;
9. Knudsen number;
10. circuit to flow timescale ratio;
11. available-power margin;
12. Debye length.

A response surface is outside contract when any declared coordinate exceeds its
log-space tolerance or when its evidence level is below the consumer's required
level.

### Closed-cycle accounting

A cycle audit requires:

- state closure in position, velocity, stored energy, and controller phase;
- the moving-medium identity

  `source work = craft work + relative-flow dissipation`;

- actuator, thermal, and numerical energy accounting;
- a net-gain interval whose lower bound remains positive;
- an evidence level high enough for the requested promotion.

A numerical closed cycle can pass while physical promotion remains false.

### Controlled force reversal

The force protocol evaluates interleaved positive-command, negative-command,
and zero-command measurements. It reports:

- reversible force estimate;
- common-mode bias;
- command-reversal asymmetry;
- lag-one autocorrelation and effective sample size;
- zero-control drift;
- an empirical-Bernstein effect interval;
- an exact sign-flip test for at most 20 pairs, or a deterministic 65,536-sample
  sign-flip audit for larger series;
- independent-replication consistency.

This is designed to reject thermal drift, facility bias, nonreversing drag,
serial correlation, and opposite-sign replications before a local transverse
force is promoted.

### Claim tiers

The claim gate is monotone and fail closed:

1. `software_finite_check`
2. `converged_model_result`
3. `chamber_local_force`
4. `closed_cycle_candidate`
5. `flight_propulsion`

Missing evidence never defaults to true. A local chamber force cannot become a
propulsion claim without robust closed-cycle, deployment, and flight evidence.

## Build and test

```sh
cmake -S assurance -B build-assurance \
  -DCMAKE_BUILD_TYPE=Release \
  -DSPACEWIND_ASSURANCE_WERROR=ON
cmake --build build-assurance
ctest --test-dir build-assurance --output-on-failure
```

Sanitizers:

```sh
cmake -S assurance -B build-assurance-sanitized \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_COMPILER=clang \
  -DSPACEWIND_ASSURANCE_SANITIZERS=ON
cmake --build build-assurance-sanitized
ctest --test-dir build-assurance-sanitized --output-on-failure
```

The tests generate deterministic JSON receipts under `assurance/output/`.

## Scientific boundary

Passing this suite establishes only that the implemented finite checks behaved
as declared for the tested inputs and toolchains. It does not establish:

- that the electric- or magnetic-sail force law is experimentally correct;
- that a chamber reproduces the solar wind outside its similarity envelope;
- that a measured force survives spacecraft power, thermal, deployment, and
  structural constraints;
- that a complete soaring cycle extracts positive net energy;
- that a flight-qualified propulsion system exists.
