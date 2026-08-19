# One-dimensional, three-velocity-component electromagnetic PIC reference

## Purpose

This module is an independent C11 reference for the first fully dynamic
particle-and-electromagnetic-field coupling layer in Spacewind. It is designed
to expose charge, field, energy, momentum, staggering, and stability mistakes
before a higher-dimensional plasma-wing solver consumes the same assumptions.

The implementation is intentionally bounded. It is a periodic **1D3V** model:
all spatial variation is along `x`, while particles and electromagnetic fields
retain three vector components.

The code lives in:

```text
assurance/include/spacewind/em_pic1d.h
assurance/src/em_pic1d.c
assurance/tests/test_em_pic1d.c
```

## Grid topology

For `N` cells of width `dx = L/N`:

- `Ex`, `Ey`, and `Ez` live on periodic faces `x_i = i dx`;
- `By` and `Bz` live at cell centers `x_(i+1/2)`;
- `Bx` is represented as one uniform scalar, because in one dimension
  `div(B) = dBx/dx` and a spatially constant `Bx` satisfies the constraint;
- charge density is deposited at cell centers;
- all three current components are deposited on faces.

This arrangement makes the longitudinal continuity and Gauss equations use the
same discrete divergence.

## Maxwell update

The transverse field step is a symmetric Yee-style split. First, magnetic
fields advance by half a step:

```text
By^(n+1/2)_i = By^n_i
             + dt/(2 dx) [Ez^n_(i+1) - Ez^n_i]

Bz^(n+1/2)_i = Bz^n_i
             - dt/(2 dx) [Ey^n_(i+1) - Ey^n_i]
```

Electric fields then advance for a full step:

```text
Ex^(n+1)_i = Ex^n_i - dt Jx_i / epsilon0

Ey^(n+1)_i = Ey^n_i
             - c^2 dt/dx [Bz^(n+1/2)_i - Bz^(n+1/2)_(i-1)]
             - dt Jy_i / epsilon0

Ez^(n+1)_i = Ez^n_i
             + c^2 dt/dx [By^(n+1/2)_i - By^(n+1/2)_(i-1)]
             - dt Jz_i / epsilon0
```

The magnetic fields then receive the second half-step using the new electric
field. Periodic topology gives zero net boundary Poynting flux, but the module
still reports particle work, field-energy change, and their total residual
separately.

The field Courant number is:

```text
C = c dt / dx
```

and is required to remain below both one and the caller's stricter declared
limit.

## Relativistic particle update

Particles are advanced with a relativistic Boris rotation. The solver evolves
proper velocity `u = gamma v`, applies half of the electric impulse, performs a
magnetic rotation with the relativistic gamma factor, and applies the second
electric half-impulse. Every accepted particle velocity must satisfy
`|v| < c`.

The accepted unwrapped displacement is based on the average of initial and
final longitudinal velocity:

```text
x^(n+1) = x^n + dt [vx^n + vx^(n+1)] / 2
```

Unwrapped positions remain authoritative. Wrapped positions are used only for
periodic interpolation and deposition.

## Charge-compatible longitudinal current

Initial and final cell-centered CIC charge are deposited independently. The
longitudinal face current is reconstructed from the exact finite-volume
continuity equation:

```text
[rho^(n+1)_i - rho^n_i]/dt
+ [Jx_(i+1) - Jx_i]/dx = 0.
```

The divergence determines current up to one uniform harmonic mode. That mode is
fixed by the full unwrapped particle transport:

```text
mean(Jx) = sum_p q_p w_p [x_p^(n+1)-x_p^n] / (L A dt).
```

The mean longitudinal electric field is audited against the corresponding
harmonic Ampere identity:

```text
mean(Ex^(n+1)) - mean(Ex^n)
+ dt mean(Jx)/epsilon0 = 0.
```

## Transverse path current

`Jy` and `Jz` use the midpoint transverse particle velocity and an exact
integral of the periodic linear face shape over the accepted unwrapped path.
The path is split only where it crosses grid faces. A complete periodic turn is
represented by its uniform shape integral rather than being aliased into a
shortest-image displacement.

The module separately checks that the sum of deposited transverse face current
matches the charge-weighted transverse particle current. This catches missing
path weight, double deposition, and periodic-seam errors even when the field
update remains finite.

## Per-step ledgers

Every step reports:

- field and particle Courant numbers;
- total charge and charge change;
- local discrete continuity residual;
- initial and final Gauss residuals;
- magnetic-divergence residual;
- harmonic longitudinal Ampere residual;
- transverse-current partition error;
- initial and final relativistic particle energy;
- initial and final electromagnetic field energy;
- midpoint `J dot E` work;
- particle-work and field-work residuals;
- total particle-plus-field energy residual;
- initial and final particle momentum;
- initial and final field momentum `epsilon0 E cross B`;
- total particle-plus-field momentum residual;
- maximum speed as a fraction of `c`;
- maximum particle displacement in cells;
- maximum number of grid segments traversed by one particle path.

The step has configurable fail-closed tolerances. A finite trajectory is not a
passing trajectory when any declared conservation, Gauss, stability,
subluminal, energy, or momentum gate fails.

## Deterministic tests

The bounded test suite includes:

1. A periodic right-going vacuum wave for 500 field steps.
2. Timestep refinement requiring the one-step field-energy defect to decrease
   under two successive halvings.
3. A neutral mixed-sign relativistic particle ensemble evolved for 120 coupled
   steps with longitudinal, transverse, and magnetic motion.
4. A node-crossing case that forces path segmentation and carries a nonzero
   harmonic longitudinal current.
5. Periodic nonneutral-charge rejection.
6. Electromagnetic CFL rejection.
7. Luminal-particle input rejection.
8. Deliberate initial Gauss-law corruption, which must remain a failed final
   Gauss audit rather than being hidden by the field update.
9. GCC, Clang, address-sanitizer, undefined-behavior-sanitizer, Clang static
   analyzer, GCC analyzer, deterministic-replay, and macOS ARM CI paths.

The standalone test performs more than three thousand deterministic checks and
emits `assurance/output/em_pic1d_vacuum_receipt.json`.

## Scientific boundary

This is a reference implementation, not a propulsion result. In particular:

- the geometry is one-dimensional and periodic;
- transverse current uses path-integrated linear weighting, but the complete
  particle-field scheme is explicit rather than exactly energy-conserving;
- energy and momentum are therefore audited against declared tolerances and
  refinement behavior, not asserted to close to roundoff for arbitrary states;
- there are no collisions, ionization, recombination, material surfaces,
  electron emitters, secondary emission, arcing, radiation, or sheath physics;
- there is no multidimensional magnetic topology or plasma-wing body;
- passing this solver does not establish a transverse plasma force, chamber
  similarity, dynamic-soaring net gain, or flight propulsion.

The next fidelity gate is a multidimensional electromagnetic kinetic solver
with a local charge-conserving current, divergence-preserving magnetic update,
open-boundary Poynting flux, particle-plus-field momentum closure, grid and
particle-count refinement, and chamber-calibrated force and torque surfaces.
