# Periodic CIC charge and current reference

## Purpose

A particle-in-cell implementation must not deposit charge at two timesteps and
then invent a current that is only approximately compatible with those charge
states. Inconsistent current deposition violates the discrete continuity
equation, causes Gauss-law drift, and can create false electromagnetic forces or
energy transfer.

Spacewind now includes an independent one-dimensional periodic reference under:

```text
assurance/include/spacewind/charge_deposition.h
assurance/src/charge_deposition.c
assurance/tests/test_charge_deposition.c
```

It is designed as a transparent verification oracle for bounded PIC tests. It
is not presented as the final local current-deposition kernel for a production
multidimensional electromagnetic PIC solver.

## 1. Charge assignment

For a periodic domain of length `L`, cross-sectional area `A`, `N` cells, and
cell width

\[
\Delta x = \frac{L}{N},
\]

each particle charge is assigned to the two nearest cell centers with the
linear cloud-in-cell shape. The density scale is

\[
\frac{q_p}{A\Delta x}.
\]

The same shape and the same compensated accumulation are used for the initial
and final particle positions.

Unwrapped particle positions are retained. Wrapping is used only to evaluate
the periodic charge shape. This distinction is essential: two particles can
have identical wrapped endpoints while one has transported charge around the
domain several times.

## 2. Discrete continuity determines current divergence

For cell-centered charge and face current,

\[
\frac{\rho_i^{n+1}-\rho_i^n}{\Delta t}
+
\frac{J_{i+1/2}-J_{i-1/2}}{\Delta x}
=0.
\]

Starting with an arbitrary reference face, the nonuniform part of the current
can be reconstructed by a prefix integration:

\[
J_{i+1/2}
=
J_{i-1/2}
-
\frac{\Delta x}{\Delta t}
\left(\rho_i^{n+1}-\rho_i^n\right).
\]

This makes the declared finite-volume continuity operator and the deposited
charge algebraically compatible up to floating-point accumulation error.

## 3. The harmonic current is not determined by continuity

On a periodic domain, adding the same constant to every face current leaves its
divergence unchanged. Charge continuity therefore cannot recover the uniform,
or harmonic, component of current.

Spacewind fixes this null-space component using the unwrapped particle
transport:

\[
\bar J_{\mathrm{target}}
=
\frac{1}{A L\Delta t}
\sum_p q_p
\left(x_p^{n+1}-x_p^n\right).
\]

Equivalently, define transported charge-turns

\[
Q_{\mathrm{turns}}
=
\sum_p q_p\frac{x_p^{n+1}-x_p^n}{L}.
\]

Then

\[
\bar J_{\mathrm{target}}
=
\frac{Q_{\mathrm{turns}}}{A\Delta t}.
\]

After reconstructing the divergent current, the implementation shifts every
periodic face by one constant so that its discrete mean equals this target.
This preserves local continuity while retaining multi-wrap transport that
would be lost from wrapped endpoints alone.

## 4. Audited outputs

Each run reports:

- particle charge sum;
- charge recovered from the initial and final CIC densities;
- both charge-deposition errors;
- transported charge-turns;
- target and deposited mean current;
- mean-current error;
- maximum unwrapped displacement in cell widths;
- the complete local and global charge-continuity audit;
- a versioned deterministic JSON receipt.

The implementation rejects invalid dimensions, nonfinite particles, zero
physical scales, undersized output buffers, and negative tolerances.

## 5. Manufactured and adversarial checks

The deterministic suite includes:

1. **Zero motion.** Initial and final densities agree and every face current is
   zero.
2. **Multi-wrap motion.** A particle returns to the same wrapped location after
   two full domain traversals. Charge density is unchanged, but the correct
   uniform current remains nonzero.
3. **Long random motion.** More than one thousand mixed-sign particles move by
   more than one hundred cell widths while local and global continuity close.
4. **Time reversal.** Swapping initial and final unwrapped positions swaps the
   charge fields and negates every current face.
5. **Permutation stability.** Reversing particle order changes the accumulated
   fields only within the predeclared floating-point tolerance.
6. **Invalid input rejection.** Zero timesteps, nonfinite positions, and
   undersized buffers fail before producing evidence.

## 6. Why this remains a reference algorithm

The prefix reconstruction is intentionally global. It proves that the two CIC
charge states and a periodic face current can satisfy one declared discrete
continuity contract while retaining the current null-space component.

A production electromagnetic PIC solver should instead use a local
charge-conserving trajectory deposition, such as a rigorously implemented
Villasenor-Buneman or Esirkepov-family scheme, sharing:

- the particle shape function;
- the field staggering;
- the boundary topology;
- the particle substepping and cell-crossing decomposition;
- the same accepted trajectory used by the pusher.

The reference reconstruction is valuable because that future local algorithm
must agree with it on charge states, continuity residuals, total transported
charge, and periodic mean current. Disagreement exposes either a topology,
shape, trajectory, or sign-convention error.

## 7. Required next promotion gates

Before current deposition can support a plasma-wing force claim, the local
kinetic solver must demonstrate:

- charge continuity at every accepted timestep;
- Gauss-law preservation without untracked charge cleaning;
- convergence under timestep, cell, particle-count, and substep refinement;
- agreement between local deposition and this independent periodic reference;
- particle-plus-field energy closure;
- field-plus-particle-plus-spacecraft momentum closure;
- Maxwell-stress and direct particle-momentum agreement;
- chamber-calibrated boundary conditions inside a declared similarity envelope.

## Scientific boundary

Passing the reference tests establishes a bounded discrete bookkeeping result.
It does not establish that the current distribution is a high-fidelity model of
solar-wind plasma, that a multidimensional kinetic solver is converged, or that
a field-sail produces useful thrust or lift.
