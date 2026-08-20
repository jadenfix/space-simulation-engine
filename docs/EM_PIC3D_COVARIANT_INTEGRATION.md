# 3D3V PIC, stress-energy, and covariant four-momentum integration

## Purpose

This integration layer prevents the 3D3V solver from promoting energy and momentum diagnostics that are internally inconsistent with the independent relativistic assurance modules.

The test lives in:

```text
assurance/tests/test_em_pic3d_covariant_integration.c
```

and emits:

```text
assurance/output/em_pic3d_covariant_receipt.json
```

## Same-timestep mapping

For every accepted 3D3V timestep, the field four-momentum is reconstructed as

```text
P_field = (U_field/c, p_field).
```

Matter four-momentum includes rest energy in the physical total:

```text
P_matter = ([E_rest + K_particle]/c, p_particle).
```

A conventional total ledger is evaluated, but its acceptance scale is tied to the dynamic exchange rather than the large constant rest-energy term.

The independent incremental ledger uses

```text
Delta P_field  = (Delta U_field/c, Delta p_field)
Delta P_matter = (Delta K_particle/c, Delta p_particle).
```

It must reproduce the solver's separately accumulated residual exactly:

```text
R_solver = (R_energy/c, R_momentum_x, R_momentum_y, R_momentum_z).
```

This directly catches unit, sign, component-order, and rest-energy-scaling mistakes.

## Boosted residual audit

The complete incremental ledger is Lorentz-transformed into an oblique inertial frame. The test then forms the residual again from the transformed terms and compares it with the direct Lorentz transform of the laboratory residual.

```text
R_prime(reconstructed) = Lambda R(laboratory).
```

The check is performed on the actual finite numerical defect, not only on an analytic test vector.

## Independent stress-energy integration

The staggered 3D fields are interpolated to cell centers. Each cell is passed through the separately implemented electromagnetic stress-energy tensor. Volume integration produces a second field four-momentum estimate.

The test requires:

- integrated stress-energy spatial momentum to agree with the 3D Yee `epsilon0 E cross B` diagnostic;
- the dominant-energy ratio not to exceed its relativistic bound beyond roundoff;
- the integrated electromagnetic stress-energy trace to vanish within the declared floating-point scale.

The stress-energy implementation does not share the 3D solver's field-momentum accumulation code.

## Current acceptance before covariant promotion

A timestep is not admitted to the covariant layer unless:

- the coordinate-split charge-conserving current passes local continuity;
- the independent direct-DFT spectral oracle requires only a near-roundoff correction;
- Gauss law and magnetic divergence pass;
- particle and field work gates pass;
- the nonlinear midpoint iteration converges;
- particle speeds remain subluminal.

Thus covariance does not mask a failed charge or field update.

## Adversarial test

After a valid step, the test inserts an untracked spatial four-impulse into the incremental field change. The physical-scale ledger must reject it.

## Scientific boundary

This layer verifies a bounded mapping among implemented diagnostics. It does not establish that the periodic 3D3V discretization is a complete relativistic Vlasov-Maxwell model, that every boundary or material source has been represented, or that a plasma wing produces useful force.

The next covariant control-volume milestone is to connect a genuinely open 3D kinetic solver to area-integrated Poynting flux, Maxwell-stress impulse, particle inward/outward four-momentum, and material force/torque on every accepted timestep.
