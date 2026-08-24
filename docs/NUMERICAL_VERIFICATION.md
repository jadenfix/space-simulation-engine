# Numerical verification and validation

## Verification versus validation

**Verification** asks whether the code solves the declared equations correctly.

**Validation** asks whether those equations represent the physical system well
enough for the intended use.

A converged solution to the wrong force law is verified but not validated. A
chamber measurement without a numerically resolved control volume is physical
data but may not validate the intended model.

## Required solver evidence

Every solver promoted into a mission-level response model should provide the
following evidence bundle.

### 1. Equation and boundary declaration

Record:

- governing equations;
- normalization;
- dependent variables;
- initial conditions;
- boundary conditions;
- source terms;
- constitutive closures;
- floors, limiters, filters, and caps;
- discretization in space and time;
- solver tolerances;
- stopping conditions.

A numerical stabilizer is part of the model and must not be omitted from the
receipt.

### 2. Method of manufactured solutions

Where practical, choose a smooth manufactured state \(u^*(x,t)\), substitute it
into the governing equation, and derive the forcing required to make it exact.
The solver should recover the designed order on this forced problem.

Manufactured solutions test implementation and boundary handling. They do not
validate the physical closure.

### 3. Analytic and limiting cases

Use independently known boundaries such as:

- vacuum Maxwell waves;
- uniform-field Boris rotation;
- cold-plasma oscillation;
- Debye shielding in a declared limit;
- two-body circular motion;
- zero-voltage and zero-field force limits;
- infinite-mass or zero-charge limits;
- symmetry planes;
- conservation under periodic boundaries.

The test suite should also verify that intentionally broken cases fail.

### 4. Mesh, timestep, and particle convergence

Refine independently:

- spatial grid;
- timestep;
- particles per cell;
- domain size;
- boundary distance;
- nonlinear iteration tolerance;
- field-solver tolerance;
- output cadence.

A result is not converged if only one combined refinement path is used. Error
cancellation can make a diagonal refinement appear stable.

For three levels, report observed order, Richardson extrapolation, fine-grid
error estimate, grid-convergence index, and whether the solutions are in the
asymptotic range.

### 5. Statistical convergence

Particle methods require both discretization and sampling evidence. Report:

- independent seeds or deterministic low-discrepancy designs;
- mean and variance of the observable;
- autocorrelation time;
- effective sample size;
- confidence or concentration interval;
- bias-versus-variance behavior under particle refinement;
- rare-event diagnostics for collection, reflection, and loss channels.

Increasing particle count on one fixed realization is not equivalent to
independent replication.

### 6. Conservation

Audit local and global:

- mass;
- charge;
- momentum;
- field plus particle energy;
- angular momentum when appropriate;
- divergence constraints;
- circuit charge and energy;
- control-system work;
- thermal loss.

A bounded global residual can hide large compensating local errors. Publish
spatial and temporal residual distributions for high-fidelity runs.

### 7. Invariance and metamorphic tests

When an exact reference solution is unavailable, transform the input in a way
that implies a known output transformation.

Examples:

- rotate all vectors and require scalar invariants to remain unchanged;
- translate an isolated N-body system and require energy and momentum to remain
  unchanged;
- permute particles or bodies and require aggregate observables to remain
  unchanged;
- reverse a steering command and require the controlled transverse component
  to reverse;
- scale heliocentric radius and require photon pressure and radial gravity to
  follow inverse-square behavior;
- apply a Galilean boost where the model is nonrelativistic and require the
  corresponding momentum and kinetic-energy transformation;
- set magnetic field to zero and recover the electric-only limit.

### 8. Independent implementation

At least one critical quantity should be computed by an implementation that
shares neither discretization nor code path with the primary solver.

Examples:

- particle momentum flux versus Maxwell stress;
- finite-volume versus particle control-volume force;
- direct pairwise gravity versus grid potential;
- analytic orbit invariants versus integrated trajectory;
- circuit current integration versus capacitor voltage change.

Agreement between two wrappers around the same hidden routine is not
independent evidence.

### 9. Precision study

Compare at least:

- strict GCC and Clang;
- optimization levels;
- fused-multiply-add settings where relevant;
- float64 against a higher-precision reference for reduced problems;
- deterministic replay;
- sanitizer builds;
- static analyzers.

Bitwise equality is not required for adaptive floating-point trajectories, but
scientific conclusions and declared tolerances must remain stable.

### 10. Boundary-condition sensitivity

Open plasma boundaries can dominate a local force result. Vary:

- domain extent;
- injection distribution;
- absorbing versus reflecting boundaries;
- potential boundary conditions;
- magnetic boundary treatment;
- neutralizing background;
- electron closure;
- reinjection method.

A response surface should include a boundary-sensitivity uncertainty, not only
a nominal solver error.

## Coupled-system verification

When plasma, circuit, thermal, structural, and flight dynamics are coupled,
verify each interface ledger.

For one coupling interval, require:

```text
plasma charge transfer
    == circuit collected charge

plasma momentum loss
    == spacecraft impulse + field momentum change + boundary flux

bus energy
    == field/storage change + collection work + conversion loss + heat

force/torque work
    == mechanical energy change + damping + structural storage change
```

Coupling cadence must be refined independently. Stable subsystem solvers can
still produce an unstable or energy-generating partitioned coupling.

## Uncertainty decomposition

Do not report one undifferentiated error bar. Separate:

- numerical discretization;
- particle sampling;
- solver tolerance;
- boundary modeling;
- closure/model-form uncertainty;
- input measurement uncertainty;
- calibration uncertainty;
- surrogate interpolation;
- extrapolation outside similarity support;
- control and state-estimation error;
- hardware variability.

Correlations must be declared. Adding all errors in quadrature silently assumes
independence.

## Promotion rule

A local response value can be consumed by the mission simulator only when its
receipt contains:

- solver version and source hash;
- complete normalized configuration;
- convergence result;
- conservation result;
- uncertainty decomposition;
- similarity support;
- evidence level;
- explicit nonclaims.

Outside that support, the mission layer must fail closed, widen the uncertainty
interval, or label the result an extrapolation. It must not silently reuse the
nearest attractive coefficient.
