# Root PIC transport integration audit

## Purpose

Independent mathematical checks are most valuable when they are exercised
against the actual engine rather than only against manufactured arrays. The
cross-layer test

```text
assurance/tests/test_pic_transport_integration.c
```

links the root electrostatic particle-in-cell implementation to the independent
periodic charge/current reference in the assurance kernel.

This is intentionally a one-way test dependency. The mission and plasma engine
do not depend on the assurance library at runtime.

## 1. Grid-topology reconciliation

The root `sw_pic1d` implementation stores charge on nodes with coordinate

\[
\xi = \frac{x}{\Delta x}.
\]

The assurance reference stores charge at cell centers with coordinate

\[
\xi_{\mathrm{ref}}
=
\frac{x_{\mathrm{ref}}}{\Delta x}-\frac12.
\]

The two shapes are identical under the explicit topology transform

\[
x_{\mathrm{ref}} = x + \frac{\Delta x}{2}.
\]

The integration test first constructs a deliberately nonuniform three-particle
state. It requires the half-cell transform to reproduce the root density and
then verifies that omitting the transform produces a detectable mismatch.
This prevents an apparent verification success caused by comparing different
staggerings as though they were the same discretization.

## 2. Unwrapped trajectory recovery

The root PIC state stores wrapped periodic positions. Current, however, depends
on transported charge and therefore on the unwrapped path.

For the bounded integration scenario, every particle displacement is required
to remain much smaller than half the domain. The test therefore recovers the
unique nearest periodic image and accumulates an unwrapped position history.
The receipt records the maximum displacement in cell widths, and the test fails
if the nearest-image assumption ceases to be unambiguous.

A production solver should store unwrapped displacement or boundary-crossing
counts directly rather than reconstructing them after the step.

## 3. Effective macro-particle charge

The root density deposition uses

\[
q_{p,\mathrm{effective}}
=
q_p w_p,
\]

where `w_p` is the macro-particle weight. The independent reference receives
this effective charge explicitly. The stationary neutralizing background is
removed before comparing particle density and has zero transport current.

## 4. Per-step audit

For one hundred accepted root PIC steps, the integration test:

1. saves the root initial deposited charge;
2. advances the actual electrostatic PIC step;
3. reconstructs the bounded unwrapped particle displacement;
4. deposits independent initial and final CIC charge;
5. reconstructs a periodic face current with the correct harmonic component;
6. audits local and global continuity;
7. compares both independent charge states against the root charge states;
8. accumulates the worst density, continuity, and mean-current errors;
9. records, but does not over-interpret, total electrostatic energy change.

The generated receipt is:

```text
assurance/output/pic_transport_integration_receipt.json
```

## 5. What this closes

The audit provides direct evidence that, for the bounded scenario:

- root macro-particle charge and the independent reference agree;
- root nodal CIC and the reference cell-centered CIC agree after an explicit
  half-cell topology transform;
- the actual accepted particle motion admits a face current that closes the
  declared finite-volume continuity equation;
- periodic transported charge is retained through an unwrapped trajectory;
- the result is deterministic across replayed CI runs.

## 6. What remains open

The current root PIC solver is electrostatic. It solves Poisson's equation from
charge at each step and does not advance a dynamic magnetic field or consume an
explicit deposited current. Therefore this integration does not yet establish:

- a local charge-conserving electromagnetic current deposition;
- discrete Ampere-Maxwell consistency;
- particle-plus-field Poynting closure;
- magnetic-divergence preservation;
- agreement between Maxwell stress and particle momentum transfer;
- three-dimensional sheath or plasma-wing force fidelity.

The next solver-level milestone is a multidimensional electromagnetic PIC step
whose accepted particle trajectory, current deposition, field update, and force
receipt all share the same discrete control volume and conservation ledger.

## Scientific boundary

Passing this integration test establishes compatibility between the bounded
root electrostatic PIC trajectory and the independent charge/current reference.
It is not chamber evidence, a validated field-sail force law, a closed soaring
cycle, or a propulsion demonstration.
