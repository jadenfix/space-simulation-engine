# Field and particle conservation ledger

## Purpose

A plasma-wing simulation is not credible merely because particles deflect or a
surface-force integral is nonzero. Every reported force must be compatible with
the discrete transport of charge, field energy, field momentum, particle
momentum, boundary fluxes, and explicitly modeled sources.

The independent assurance module in:

```text
assurance/include/spacewind/field_ledger.h
assurance/src/field_ledger.c
assurance/tests/test_field_ledger.c
```

implements this bookkeeping without depending on the mission-scale force model.

The first implementation is a **two-dimensional, per-unit-depth audit** for
cell-centered charge and staggered face fields. It is intentionally an audit
surface, not a replacement particle-in-cell solver.

## Discrete storage contract

For a grid with `nx` by `ny` cells:

- cell-centered arrays contain `nx * ny` values;
- x-face arrays contain `(nx + 1) * ny` values;
- y-face arrays contain `nx * (ny + 1)` values;
- `dx_m`, `dy_m`, and `dt_s` are positive SI values;
- all charge densities are in coulombs per cubic metre;
- all current densities are in amperes per square metre;
- the resulting integrated two-dimensional charge and current are reported per
  metre of unresolved depth.

Invalid dimensions, non-finite values, zero spacing, zero timesteps, and
negative tolerances are rejected rather than converted into a failed metric.

## 1. Local charge continuity

For every cell,

\[
R_Q^{i,j}
=
\frac{\rho_{i,j}^{n+1}-\rho_{i,j}^{n}}{\Delta t}
+
\frac{J_{x,i+1/2,j}-J_{x,i-1/2,j}}{\Delta x}
+
\frac{J_{y,i,j+1/2}-J_{y,i,j-1/2}}{\Delta y}.
\]

The audit reports:

- mean absolute residual;
- root-mean-square residual;
- maximum local residual;
- RMS and maximum residual normalized by the declared transport terms;
- signed residual integrated over the domain.

A global charge check is evaluated independently:

\[
R_{Q,\mathrm{global}}
=
Q^{n+1}-Q^n
+
\Delta t\,I_{\mathrm{out}}.
\]

The outward current uses the actual staggered boundary faces and their outward
normals. The audit also checks the discrete divergence theorem by comparing the
global balance with the cell-integrated local residual.

A simulation can therefore fail in three distinct ways:

1. local charge is moved between cells incorrectly;
2. the domain gains or loses unexplained total charge;
3. the local and boundary ledgers disagree.

## 2. Gauss's law

For cell-centered deposited charge and face-centered electric field,

\[
R_G
=
\epsilon_0
\left(
\frac{E_{x,i+1/2,j}-E_{x,i-1/2,j}}{\Delta x}
+
\frac{E_{y,i,j+1/2}-E_{y,i,j-1/2}}{\Delta y}
\right)
-
\rho_{i,j}.
\]

The residual is reported in charge-density units. This avoids presenting an
arbitrary normalized Poisson residual as if it were a physical Gauss-law
error.

A small iterative-solver residual does not imply a small Gauss residual unless
the deposited charge, boundary conditions, and field reconstruction use the
same discrete operators.

## 3. Magnetic divergence

For staggered magnetic field,

\[
R_B
=
\frac{B_{x,i+1/2,j}-B_{x,i-1/2,j}}{\Delta x}
+
\frac{B_{y,i,j+1/2}-B_{y,i,j-1/2}}{\Delta y}.
\]

The manufactured test constructs a nontrivial field from a discrete stream
function:

\[
B_x = \frac{\partial_h \psi}{\partial y},
\qquad
B_y = -\frac{\partial_h \psi}{\partial x},
\]

so the mixed finite differences cancel algebraically. A deliberately corrupted
face then verifies that the gate detects an injected numerical magnetic
monopole.

## 4. Poynting-energy ledger

The sign convention is:

- `particle_work_J` is positive when fields transfer energy to particles or
  matter;
- `outward_poynting_energy_J` is positive when field energy leaves the control
  volume;
- `impressed_source_energy_J` is positive when an explicit source injects
  energy.

The residual is

\[
R_E
=
U_f-U_i
+
W_{J\cdot E}
+
E_{S,\mathrm{out}}
-
E_{\mathrm{source}}.
\]

A numerical field-heating term must not be hidden as generic dissipation. It
must appear in the residual or be represented as a separately named modeled
channel.

## 5. Field-plus-mechanical momentum

The field momentum density is

\[
\mathbf g
=
\epsilon_0\mathbf E\times\mathbf B
=
\frac{\mathbf S}{c^2}.
\]

For an integrated control volume the audit uses

\[
\mathbf R_P
=
\left(
\mathbf P_{f,\mathrm{field}}
-
\mathbf P_{i,\mathrm{field}}
\right)
+
\mathbf J_{\mathrm{mechanical}}
+
\mathbf J_{\mathrm{Maxwell,out}}
-
\mathbf J_{\mathrm{external}}.
\]

A spacecraft impulse is not accepted merely because a tether or sail force
routine returns it. The corresponding momentum must be removed from plasma,
fields, an impressed source, or a declared boundary flux.

## 6. Maxwell stress and local electromagnetic identities

The module evaluates

\[
u_{\mathrm{EM}}
=
\frac{1}{2}
\left(
\epsilon_0 E^2+\frac{B^2}{\mu_0}
\right),
\]

\[
\mathbf S
=
\frac{\mathbf E\times\mathbf B}{\mu_0},
\qquad
\mathbf g=\frac{\mathbf S}{c^2},
\]

and the Maxwell traction

\[
\mathbf t
=
\left[
\epsilon_0
\left(
\mathbf E\mathbf E-\frac{1}{2}E^2\mathbf I
\right)
+
\frac{1}{\mu_0}
\left(
\mathbf B\mathbf B-\frac{1}{2}B^2\mathbf I
\right)
\right]\hat{\mathbf n}.
\]

It also records the two local scalar combinations

\[
B^2-\frac{E^2}{c^2},
\qquad
\frac{\mathbf E\cdot\mathbf B}{c},
\]

and checks the vacuum dominant-energy inequality

\[
\frac{|\mathbf S|}{c\,u_{\mathrm{EM}}}\le 1.
\]

The tests rotate the electric field, magnetic field, and surface normal
together and require the calculated traction to rotate covariantly.

## 7. Pass semantics

Every residual gate accepts both:

- an absolute tolerance in the physical units of the residual; and
- a relative tolerance against the magnitude of the terms being balanced.

A gate passes when the residual is below its absolute tolerance or its relative
tolerance. Charge continuity additionally requires both local and global
closure.

Tolerances are part of the model contract. They must be selected before looking
at the desired force result and should be tightened through timestep, grid,
particle-count, and domain-size refinement.

## 8. Manufactured and adversarial checks

The deterministic test suite includes:

- charge at the next timestep manufactured from an arbitrary staggered current;
- a one-cell charge corruption that must fail local and global continuity;
- electric field with charge manufactured from its discrete divergence;
- a one-cell Gauss-law corruption;
- nontrivial divergence-free magnetic field from a discrete stream function;
- an injected magnetic-face error;
- exact Poynting-energy closure and an untracked-energy corruption;
- exact field-plus-mechanical momentum closure and an untracked impulse;
- thousands of low-discrepancy electromagnetic states;
- Poynting/momentum consistency;
- dominant-energy checks;
- Maxwell-traction rotational covariance;
- rejection of a zero surface normal.

The receipt is written to:

```text
assurance/output/field_ledger_receipt.json
```

## 9. Required integration with a kinetic solver

The assurance API can audit arrays supplied by any solver, but a future
electromagnetic PIC implementation should make the following non-optional at
every accepted timestep:

1. deposit charge and current with one charge-conserving scheme;
2. audit local continuity before field advancement;
3. preserve or explicitly clean Gauss's law without hiding removed charge;
4. audit magnetic divergence after the Maxwell update;
5. account for particle kinetic-energy change using the same interpolated
   electric field used by the pusher;
6. integrate boundary Poynting and Maxwell-stress flux on the same control
   surface;
7. compare field stress, particle momentum change, and spacecraft impulse;
8. emit a versioned receipt containing raw scales, residuals, tolerances, grid,
   timestep, boundary conditions, and source terms.

A force response surface may only be promoted into the mission simulator after
these ledgers converge under refinement and remain compatible with the chamber
similarity contract.

## Scientific boundary

Passing these checks means that the declared discrete ledgers close within
their configured tolerances for a bounded run. It does not prove that:

- the underlying plasma model is complete;
- the boundary conditions represent the solar wind;
- a field-sail produces a useful transverse force;
- the power, thermal, material, deployment, or control systems are feasible;
- a closed dynamic-soaring cycle has positive physical net energy;
- flight propulsion has been demonstrated.
