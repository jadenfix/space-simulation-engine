# Assurance and falsification contract

## Purpose

The simulator has two different jobs:

1. propagate a declared mathematical model;
2. prevent the output of that model from being described as stronger evidence
   than it is.

The second job is handled by the independent C11 assurance kernel. A trajectory
that integrates successfully is not automatically a physically valid force
model, a closed energy cycle, an experimentally measured effect, or a
flight-qualified propulsion system.

## Dependency chain

The promotion chain is deliberately one-way:

```text
well-formed equations and units
        -> finite numerical execution
        -> conservation and invariant closure
        -> mesh/time/particle convergence
        -> independent implementation agreement
        -> similarity-qualified local physics
        -> finite power/thermal/structural feasibility
        -> reversible chamber force and replication
        -> robust closed-cycle net gain
        -> deployment and flight evidence
```

Failure at an earlier layer invalidates every stronger downstream claim.
Success at an earlier layer does not imply success at a later one.

## 1. Dimensional contract

Every physical equation must be dimensionally homogeneous. The assurance API
represents a quantity using a numerical value and seven SI base exponents:

\[
[M,L,T,I,\Theta,N,J].
\]

Multiplication adds exponent vectors, division subtracts them, and addition is
permitted only when vectors are equal. This catches errors such as adding a
force to an energy or treating a voltage as an electric field.

Dimension checks are necessary but not sufficient: a dimensionally valid
formula can still be physically wrong.

## 2. Interval arithmetic

Uncertain scalar inputs are represented by outward-rounded intervals. For
example,

\[
[a,b]+[c,d]=[\operatorname{down}(a+c),
             \operatorname{up}(b+d)].
\]

Multiplication evaluates all four endpoint products. Division is refused when
the denominator interval contains zero. Square root is refused for negative
lower bounds.

A robust positive-gain claim requires the entire resulting interval to be
positive. A point estimate above zero is insufficient when the interval still
crosses zero.

## 3. Conservation ledgers

For an open control volume, the kernel audits separate ledgers.

### Mass

\[
R_m=m_f-(m_i+m_{in}-m_{out}).
\]

### Charge

\[
R_Q=Q_f-(Q_i+Q_{in}-Q_{out}).
\]

### Momentum

\[
\mathbf R_p=
\mathbf p_f-
\left(
\mathbf p_i+\mathbf J_{external}
+\mathbf p_{in}-\mathbf p_{out}
\right).
\]

### Energy

\[
R_E=E_f-
\left(
E_i+W_{source}+Q_{heat}+E_{in}-E_{out}
\right).
\]

Absolute and relative tolerances are both recorded. A small relative residual
cannot hide a large absolute leak, and a tiny absolute residual is not divided
by an unstable zero scale.

## 4. Relativistic invariant

For a three-velocity \(\mathbf v\), the four-velocity is

\[
u^\mu=\gamma(c,\mathbf v),
\qquad
\gamma=\frac{1}{\sqrt{1-v^2/c^2}}.
\]

The Minkowski norm must remain

\[
\eta_{\mu\nu}u^\mu u^\nu=-c^2.
\]

Luminal and superluminal three-velocities are rejected before construction.
This check validates a declared special-relativistic identity; it is not a
solution of the Einstein field equations.

## 5. Electromagnetic invariants

The Lorentz force is decomposed as

\[
\mathbf F=q\mathbf E+q\mathbf v\times\mathbf B.
\]

The magnetic component must satisfy

\[
\mathbf v\cdot(q\mathbf v\times\mathbf B)=0,
\]

so a magnetic field alone cannot directly change particle kinetic energy.

The kernel also computes

\[
u_{EM}=\frac12\left(\epsilon_0E^2+\frac{B^2}{\mu_0}\right),
\qquad
\mathbf S=\frac{\mathbf E\times\mathbf B}{\mu_0},
\]

and the field invariants

\[
I_1=B^2-\frac{E^2}{c^2},
\qquad
I_2=\frac{\mathbf E\cdot\mathbf B}{c}.
\]

Rotating the full physical state must preserve scalar invariants and rotate
vector quantities covariantly.

## 6. Heliospheric consistency

For a steady radial wind, two samples at different radii must obey mass-flux
continuity:

\[
nv_rr^2=\text{constant}.
\]

A divergence-free radial magnetic field requires

\[
B_rr^2=\text{constant}.
\]

The ideal-MHD electric field must satisfy

\[
\mathbf E=-\mathbf v\times\mathbf B,
\quad
\mathbf E\cdot\mathbf B=0,
\quad
\mathbf E\cdot\mathbf v=0.
\]

These checks catch incompatible density, velocity, magnetic-field, and electric
field profiles before they are consumed by a spacecraft model.

## 7. Momentum-source upper bounds

### Photons

For irradiance \(I\), area \(A\), and reflectivity \(R\in[0,1]\), the normal
force cannot exceed

\[
F_{\gamma,max}=(1+R)\frac{IA}{c}.
\]

### Plasma

A declared interaction area and relative speed imply a momentum-flux scale

\[
F_{plasma,max}=\kappa\rho v_{rel}^2A,
\]

where \(\kappa\) is an explicit momentum multiplier. The simulator may explore
a configured \(\kappa\), but it may not silently exceed the resulting bound.
The effective interaction area itself remains a model requiring kinetic or
experimental validation.

## 8. Orbit and N-body invariants

For a two-body state,

\[
\epsilon=\frac{v^2}{2}-\frac{\mu}{r},
\qquad
\mathbf h=\mathbf r\times\mathbf v,
\]

and

\[
\mathbf e=
\frac{1}{\mu}
\left[
\left(v^2-\frac{\mu}{r}\right)\mathbf r
-(\mathbf r\cdot\mathbf v)\mathbf v
\right].
\]

For N-body snapshots the kernel independently computes total mass, center of
mass, momentum, angular momentum, kinetic energy, pairwise gravitational
potential energy, and total energy. Translation, rotation, permutation, and
Galilean-boost metamorphisms are checked where appropriate.

## 9. Convergence

Three resolutions \(f_h,f_{h/r},f_{h/r^2}\) produce an observed order

\[
p=\frac{\log\left|(f_h-f_{h/r})/(f_{h/r}-f_{h/r^2})\right|}
        {\log r}.
\]

The extrapolated value is

\[
f_{ext}=f_{h/r^2}+
\frac{f_{h/r^2}-f_{h/r}}{r^p-1}.
\]

The kernel records monotonicity, observed order, extrapolation, a fine-grid
error estimate, grid-convergence index, and asymptotic-range ratio. Merely
obtaining three finite answers is not convergence.

## 10. Similarity-qualified extrapolation

A local kinetic simulation or chamber result is consumable only inside its
validated similarity neighborhood. The signature includes voltage, Debye,
gyroradius, inertial-length, Mach, beta, Knudsen, circuit-timescale, and power
coordinates. Distance is evaluated in log space because many plasma scales
span orders of magnitude.

Evidence level is part of the contract. Numerically matching dimensionless
coordinates does not turn a simulation into an experiment.

## 11. Closed dynamic-soaring cycle

A valid cycle must return not just to a nearby position, but to the same
relevant state: position, velocity, stored energy, and controller phase.

For plasma force \(\mathbf F\), plasma velocity \(\mathbf u\), and spacecraft
velocity \(\mathbf v\), the moving-medium identity is

\[
\mathbf F\cdot\mathbf u=
\mathbf F\cdot\mathbf v+
\mathbf F\cdot(\mathbf u-\mathbf v).
\]

Integrated over a cycle,

\[
W_{source}=W_{craft}+D_{relative}.
\]

The net interval must then subtract actuator and thermal costs:

\[
G=W_{craft}-E_{actuator}-Q_{thermal}.
\]

Promotion requires the lower uncertainty bound of \(G\) to remain positive.

## 12. Experimental force reversal

A physical transverse-force claim must survive:

- positive command;
- negative command;
- zero command;
- common-mode subtraction;
- reversal-symmetry limits;
- drift limits;
- serial-correlation limits;
- a sign-flip randomization test;
- independent replication.

The reversible estimator is

\[
\hat F=\frac{\bar F_+-\bar F_-}{2},
\]

while common-mode bias is

\[
\hat b=\frac{\bar F_++\bar F_-}{2}-\bar F_0.
\]

A large positive reading that does not reverse with command is not accepted as
controlled lift.

## 13. Fail-closed claim promotion

The machine-readable tiers are:

| Tier | Minimum interpretation |
|---|---|
| none | no promotable result |
| software finite check | bounded code behavior only |
| converged model result | numerically qualified conditional model |
| chamber local force | reversible, replicated local physical force |
| closed-cycle candidate | robust positive closed cycle under qualified force law |
| flight propulsion | deployment and flight-qualified propulsion evidence |

The current project must remain below the propulsion tier until real evidence
satisfies every upstream gate.
