# Maxwell reference to covariant stress-energy integration

This cross-layer test ensures that the periodic transverse Maxwell oracle and
the independently implemented covariant electromagnetic stress-energy layer
assign the same global energy and momentum to the same fields.

For each sample, the test constructs

\[
\mathbf E=(0,E_y,E_z),
\qquad
\mathbf B=(0,B_y,B_z),
\]

computes the local stress-energy tensor, and integrates

\[
P^0_{\rm field}=\frac{1}{c}\int T^{00}\,dV,
\]

\[
P^i_{\rm field}=\frac{1}{c}\int T^{0i}\,dV.
\]

The result must equal the Maxwell module's independently accumulated field
energy and Poynting momentum. The right/left characteristic decomposition must
also satisfy

\[
U=U_++U_-,
\qquad
P_x=\frac{U_+-U_-}{c}.
\]

A long vacuum trajectory is then treated as a field-only four-momentum control
volume. Initial and final integrated four-momentum must close with zero
external, inflow, outflow, and matter terms. Injecting an invented external
four-impulse must fail.

For a prescribed transverse current, the time component is checked against the
independently recorded midpoint current work:

\[
P^0_f-P^0_i+W_{J\cdot E}/c=0.
\]

The spatial source impulse is not promoted here because the current is
prescribed and the exact time-centered magnetic force density is not yet
returned by the Maxwell reference API.

This is a consistency test between independent numerical layers. It is not a
self-consistent particle-field momentum proof and is not plasma-wing force or
flight evidence.
