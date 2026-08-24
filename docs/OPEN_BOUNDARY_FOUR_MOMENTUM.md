# Open-boundary electromagnetic four-momentum assurance

## Purpose

Periodic kinetic tests cannot expose mistakes in energy or momentum that leave a
finite simulation domain. This assurance layer therefore converts local
electromagnetic stress-energy at an arbitrary oriented surface into an outward
four-momentum flux, and combines it with discrete particle crossings in the
existing covariant control-volume ledger.

The implementation lives in:

```text
assurance/include/spacewind/open_boundary.h
assurance/src/open_boundary.c
assurance/tests/test_open_boundary.c
```

## Electromagnetic boundary flux

For a unit outward spatial normal `n_i`, the outward electromagnetic
four-momentum flux density is

```text
F^nu = n_i T^{i nu},
```

where `T^{mu nu}` is the independently implemented electromagnetic
stress-energy tensor using the repository's `(-,+,+,+)` metric and `x^0 = c t`
convention.

The time component satisfies

```text
c F^0 = S dot n,
```

so it is the outward Poynting power density. The three spatial components are
the outward field momentum flux. Under this stress-energy convention they are
the negative of the conventional Maxwell traction returned by the independent
`field_ledger` implementation:

```text
F_spatial = - sigma_Maxwell dot n.
```

The test evaluates both constructions independently over thousands of bounded
field/normal states.

For constant flux over area `A` and time interval `dt`, the outward boundary
four-impulse is

```text
Delta P_boundary^nu = A dt F^nu.
```

## Particle transport

A particle that crosses an oriented control surface carries its relativistic
four-momentum

```text
p^mu = (gamma m c, gamma m v).
```

Macro-particle weight is applied before the crossing is classified. The sign of
`v dot n` determines whether the crossing contributes to outward or inward
four-momentum transport. Luminal/superluminal states, nonpositive mass or
weight, and degenerate normals fail closed.

## Open control-volume identity

The audit uses the existing covariant ledger with explicit boundary terms:

```text
P_final
=
P_initial
+ I_external
+ P_inward
- P_outward.
```

`P_outward` contains both electromagnetic stress-energy transport and particle
crossings. No boundary term defaults to an implicit numerical loss.

## Manufactured tests

The deterministic suite includes:

- 3,000 arbitrary field and surface-normal states checking `S dot n` and
  Maxwell-traction consistency;
- a right-moving plane-wave slab escaping a right boundary;
- a left-moving plane-wave slab escaping a left boundary;
- relativistic outward and inward particle crossings;
- a combined field-plus-particle open control volume;
- rejection of a deliberately corrupted Maxwell-stress boundary impulse;
- zero-normal and luminal-particle rejection.

For the plane-wave tests, the initial four-momentum is computed independently
from stress-energy integrated over a slab of length `c dt`; the boundary
four-impulse is computed from flux integrated over area and time. Their closure
is therefore not obtained by copying the same value into both sides of the
ledger.

## Scientific boundary

This module validates the sign, units, and bookkeeping contract for local
stress-energy and discrete particle transport through a control surface. It is
not yet a multidimensional kinetic solver and does not establish a plasma-wing
force.

The next fidelity gate is to make a 2D3V/3D3V electromagnetic PIC domain emit
these field and particle boundary terms on every accepted timestep, then demand
joint volume-plus-boundary four-momentum convergence under grid, timestep,
particle-count, shape-order, and domain refinement before any force/torque
response surface can be promoted.
