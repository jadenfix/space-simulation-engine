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

## 2. Explicit periodic transport state

The root PIC state now keeps three separate position concepts:

- wrapped particle position used for field interpolation and deposition;
- unwrapped position after the accepted drift;
- previous unwrapped position before that drift.

For particle `p` and accepted step `n`,

\[
\Delta x_p^n
=
x_{p,\mathrm{unwrapped}}^{n+1}
-
x_{p,\mathrm{unwrapped}}^n.
\]

The wrapped state is always reconstructed from the unwrapped state:

\[
x_{p,\mathrm{wrapped}}
=
x_{p,\mathrm{unwrapped}}\bmod L.
\]

This removes the nearest-image assumption from the audit and preserves
multi-wrap transport even when a particle moves through several periodic
copies in one accepted step.

The engine exposes:

```c
bool sw_pic1d_sync_unwrapped_positions(sw_pic1d *pic);
bool sw_pic1d_unwrapped_positions_consistent(
    const sw_pic1d *pic,
    double absolute_tolerance_m
);
```

Manual particle initialization must be followed by explicit synchronization.
If external code later changes wrapped positions without updating the transport
state, the consistency check fails rather than silently inventing a path.

A dedicated test advances a zero-charge particle by three and a half domain
lengths in one step. It requires the full unwrapped displacement to survive
while the wrapped state remains inside the periodic domain.

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
3. verifies wrapped/unwrapped consistency;
4. reads the engine's previous and final unwrapped positions directly;
5. deposits independent initial and final CIC charge;
6. reconstructs a periodic face current with the correct harmonic component;
7. audits local and global continuity;
8. compares both independent charge states against the root charge states;
9. accumulates the worst density, continuity, and mean-current errors;
10. records, but does not over-interpret, total electrostatic energy change.

The generated receipt is:

```text
assurance/output/pic_transport_integration_receipt.json
```

The receipt schema is version 2 and explicitly states that accepted transport
came from stored unwrapped positions rather than inferred nearest images.

## 5. What this closes

The audit provides direct evidence that, for the bounded scenario:

- root macro-particle charge and the independent reference agree;
- root nodal CIC and the reference cell-centered CIC agree after an explicit
  half-cell topology transform;
- the actual accepted particle motion admits a face current that closes the
  declared finite-volume continuity equation;
- periodic transported charge is retained through an explicit unwrapped
  trajectory;
- wrapped state remains consistent with its unwrapped source;
- multi-wrap motion is not aliased into a shorter displacement;
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

The next solver-level milestone is a local trajectory current deposited from
the same explicit unwrapped path, followed by a multidimensional
electromagnetic PIC step whose particle push, field update, and force receipt
share one discrete control volume and conservation ledger.

## Scientific boundary

Passing this integration test establishes compatibility between the bounded
root electrostatic PIC trajectory and the independent charge/current reference.
It is not chamber evidence, a validated field-sail force law, a closed soaring
cycle, or a propulsion demonstration.
