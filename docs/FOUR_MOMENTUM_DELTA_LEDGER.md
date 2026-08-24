# Incremental four-momentum ledger

## Why a delta ledger is necessary

A conventional conservation audit forms total initial and final four-momentum and subtracts them. That is physically valid, but it can be numerically weak when a tiny same-timestep exchange is represented as the difference between very large totals.

For particles, the constant rest-energy contribution can dominate the time component:

```text
P_matter^0 = [sum_p w_p m_p c^2 + K_particle] / c.
```

For fields, a large nearly constant background energy can create the same problem. Floating-point subtraction can lose the small dynamic residual that the audit is supposed to measure.

The incremental ledger stores changes directly:

```text
Delta P_field
Delta P_matter
I_external
P_in
P_out.
```

Its residual is

```text
R = Delta P_field + Delta P_matter
    - I_external - P_in + P_out.
```

This avoids reconstructing a small change by subtracting large initial and final totals.

## API

```text
assurance/include/spacewind/four_momentum_delta.h
assurance/src/four_momentum_delta.c
```

The audit reports:

- the residual four-vector;
- Euclidean residual norm used for the bounded numerical gate;
- a declared physical dynamic scale;
- relative residual;
- finite and pass flags.

The physical scale is built from the field change, matter change, supplied boundary/external transport, and an explicit caller floor. Constant particle rest energy is not included in the acceptance scale.

## Lorentz transformation

Every term in the delta ledger is a four-vector. The module transforms each term independently with the existing Lorentz matrix. The 3D3V integration test requires the residual reconstructed in a boosted frame to agree with the Lorentz transform of the laboratory-frame residual.

This checks covariance of the numerical bookkeeping. It does not make a noncovariant discretization covariant, and it does not prove Lorentz invariance of the complete grid algorithm.

## Relationship to total ledgers

The total and incremental ledgers are both retained:

- the total ledger checks the conventional initial/final accounting contract;
- the incremental ledger protects the same-timestep acceptance gate from catastrophic cancellation;
- the solver's separately accumulated energy and three-momentum residuals must map exactly into the incremental four-residual.

A deliberate untracked four-impulse must fail both the intended physical-scale gate and the boosted residual comparison.

## Scientific boundary

A closed incremental ledger proves only that the declared numerical terms reconcile. It does not establish that all physically relevant source, boundary, material, radiation, or uncertainty terms have been modeled.
