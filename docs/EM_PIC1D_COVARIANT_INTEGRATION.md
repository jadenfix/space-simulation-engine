# EM PIC to covariant four-momentum integration

## Purpose

This cross-layer test connects the bounded periodic 1D3V electromagnetic PIC
reference directly to the independent special-relativistic four-momentum
ledger. The goal is to make one accepted kinetic timestep accountable in the
same energy-momentum representation used by the covariant field tests.

The test lives in:

```text
assurance/tests/test_em_pic1d_covariant_integration.c
```

and consumes:

```text
assurance/include/spacewind/em_pic1d.h
assurance/include/spacewind/covariant_em.h
```

## Same-timestep construction

For every accepted EM-PIC step, the field four-momentum is constructed as

```text
P_field = (U_field / c, p_field_x, p_field_y, p_field_z)
```

and the matter four-momentum is

```text
P_matter = ((E_rest + K_particle) / c,
            p_particle_x,
            p_particle_y,
            p_particle_z).
```

The weighted rest energy is constant because particle mass and macro-particle
weight do not change in this reference. It is included in the physical
four-vector but deliberately excluded from the numerical acceptance scale.
Otherwise a very large constant rest-energy component could make a poor kinetic
energy or momentum residual appear artificially small.

For a periodic closed domain with no declared external source, inflow, or
outflow, the ledger evaluates

```text
R^mu = P_final^mu - P_initial^mu.
```

Its time component must reproduce the EM-PIC total energy residual:

```text
c R^0 = Delta E_particle + Delta E_field,
```

while its spatial components must reproduce the EM-PIC total momentum residual.
This is checked independently rather than inferred from a shared pass flag.

## Dynamic acceptance scale

The pass scale is built from dynamic rather than constant rest-energy terms:

```text
max(
    |K_initial + U_field_initial| / c,
    |K_final + U_field_final| / c,
    |p_total_initial|,
    |p_total_final|
).
```

The covariant ledger is therefore unable to pass merely because rest energy is
large compared with the numerical exchange being audited.

## Lorentz-covariance check

Every field, matter, external, inflow, and outflow four-vector in the same-step
ledger is independently transformed into a moving inertial frame. The
numerical residual from the boosted ledger must equal the Lorentz transform of
the laboratory-frame residual. This tests covariance of the bookkeeping itself,
not just covariance of isolated electromagnetic tensors.

## Bounded deterministic experiment

The current test uses a neutral mixed-sign relativistic particle ensemble,
transverse electromagnetic wave content, a uniform longitudinal magnetic
field, and 64 coupled timesteps. At each step it requires:

- finite EM-PIC diagnostics;
- local charge-continuity closure;
- initial and final Gauss-law closure;
- subluminal particle motion;
- the existing EM-PIC energy and momentum gates;
- agreement between three-plus-one ledger residuals and EM-PIC residuals;
- a dynamic-scale four-momentum gate;
- covariance of the numerical residual under an oblique Lorentz boost.

A deliberate untracked spatial four-impulse is added after the sequence and
must be rejected.

## Scientific boundary

This is a bookkeeping and covariance integration test for the bounded periodic
1D3V solver. It does not establish exact energy-momentum conservation for all
states. The underlying explicit particle-field algorithm retains finite
numerical defects and is accepted only inside declared tolerances.

It also does not establish:

- a two- or three-dimensional kinetic plasma interaction;
- open-boundary Poynting or Maxwell-stress transport;
- a physical plasma-wing lift coefficient;
- chamber-to-solar-wind similarity;
- useful dynamic-soaring net gain;
- flight propulsion.

The next fidelity step is a multidimensional electromagnetic kinetic control
volume with local charge-conserving current, open field and particle fluxes,
Maxwell-stress momentum transport, refinement in grid/timestep/particle count,
and force/torque response surfaces calibrated against experiment.
