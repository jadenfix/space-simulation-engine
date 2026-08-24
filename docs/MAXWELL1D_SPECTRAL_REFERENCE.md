# Periodic 1D3V Maxwell spectral reference

This module is an independent, direct-DFT reference for the transverse
one-dimensional Maxwell system. It is intentionally slow and small so that a
finite-difference or particle-in-cell implementation can be checked against a
second numerical method rather than against itself.

The represented fields are

\[
E_y(x),\quad E_z(x),\quad B_y(x),\quad B_z(x).
\]

For variation only along \(x\), vacuum Maxwell evolution is

\[
\partial_t E_y=-c^2\partial_x B_z,\qquad
\partial_t B_z=-\partial_x E_y,
\]

\[
\partial_t E_z=+c^2\partial_x B_y,\qquad
\partial_t B_y=+\partial_x E_z.
\]

The code evolves the characteristic fields

\[
R_y^\pm=E_y\pm cB_z,
\]

\[
R_z^+=E_z-cB_y,\qquad R_z^-=E_z+cB_y.
\]

The plus characteristics translate by \(+c\Delta t\), and the minus
characteristics translate by \(-c\Delta t\). Each translation is performed by
a direct discrete Fourier transform and an exact modal phase rotation. No CFL
approximation is introduced by the vacuum propagator.

## Explicit null-space contract

Only odd sample counts are accepted. A collocated even grid contains a real
Nyquist mode whose continuously shifted sine partner is absent from the sampled
basis. Rejecting even grids is more auditable than silently discarding the
imaginary component of a translated Nyquist mode.

## Energy and momentum

The field energy is

\[
U=A\Delta x\sum_i\left[
\frac{\epsilon_0}{2}(E_y^2+E_z^2)+
\frac{1}{2\mu_0}(B_y^2+B_z^2)
\right].
\]

The longitudinal field momentum is

\[
P_x=\frac{A\Delta x}{c^2}\sum_i
\frac{E_yB_z-E_zB_y}{\mu_0}.
\]

Vacuum steps require both quantities to close. Right- and left-moving
characteristic energies are also reported separately.

A prescribed transverse current is applied with symmetric source splitting.
Each half-kick uses the exact midpoint work identity

\[
\Delta U_E+\Delta t\,A\Delta x\sum_i
\mathbf J_i\cdot\frac{\mathbf E_i^{\rm before}+
\mathbf E_i^{\rm after}}{2}=0.
\]

The vacuum translation between the kicks preserves field energy, so the full
step closes field-energy change plus current work to roundoff.

## Tests

The deterministic suite checks:

- exact one-cell translation of finite right- and left-moving modes;
- \(P_x=+U/c\) and \(P_x=-U/c\) for pure traveling waves;
- arbitrary-state positive/negative-time reversal;
- bitwise deterministic replay;
- transverse-polarization covariance;
- prescribed-current midpoint work closure;
- rejection of corrupted source work, nonfinite inputs, mutable-field aliasing,
  undersized grids, and even-grid Nyquist ambiguity.

## Boundary of the result

This is not yet a self-consistent electromagnetic particle-in-cell solver. The
current arrays are prescribed, longitudinal charge and current are absent, and
no particle pusher consumes these fields in this module. Its purpose is to act
as an independent electromagnetic oracle for the staggered 1D3V PIC and Yee
implementations. Passing it does not establish plasma-wing force, net dynamic
soaring, chamber similarity, or flight propulsion.
