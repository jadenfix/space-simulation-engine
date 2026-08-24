# Covariant electromagnetic stress-energy assurance

This layer verifies that the low-level electromagnetic quantities used by the
space and plasma solvers obey special-relativistic transformation and
conservation contracts.

The module uses coordinates \(x^0=ct\) and metric

\[
g_{\mu\nu}=\operatorname{diag}(-1,+1,+1,+1).
\]

The contravariant field tensor is defined by

\[
F^{0i}=E_i/c,
\qquad
F^{i0}=-E_i/c,
\qquad
F^{ij}=\epsilon_{ijk}B_k.
\]

The electromagnetic stress-energy tensor is evaluated from

\[
T^{\mu\nu}=\frac{1}{\mu_0}
\left(
F^{\mu\alpha}F^{\nu}{}_{\alpha}
-\frac14 g^{\mu\nu}F^{\alpha\beta}F_{\alpha\beta}
\right).
\]

This construction yields

\[
T^{00}=u_{\rm EM},
\qquad
T^{0i}=S_i/c,
\]

where \(u_{\rm EM}\) is electromagnetic energy density and \(\mathbf S\) is
Poynting flux.

## Lorentz covariance

For a frame velocity \(\mathbf v\), the implementation constructs the full
Lorentz matrix and transforms four-vectors and rank-two contravariant tensors.
It independently transforms \(\mathbf E\) and \(\mathbf B\) using the standard
SI field formulas, then requires the reconstructed field tensor and
stress-energy tensor to match the tensor-transformed versions.

The tests check invariance of

\[
B^2-E^2/c^2
\]

and

\[
\mathbf E\cdot\mathbf B/c.
\]

They also require the stress-energy trace to vanish and the electromagnetic
dominant-energy bound

\[
|\mathbf S|\le c u_{\rm EM}
\]

to hold for every deterministic sample.

## Four-current and force density

The four-current is

\[
J^\mu=(c\rho,\mathbf J),
\]

and the four-force density is

\[
f^\mu=F^\mu{}_{\nu}J^\nu.
\]

Its time component is \(\mathbf J\cdot\mathbf E/c\), while its spatial
components are

\[
\rho\mathbf E+\mathbf J\times\mathbf B.
\]

The suite transforms the field, current, and force independently and requires
covariance under thousands of subluminal boosts.

## Matter mass shell and total four-momentum

A particle with rest mass \(m\) and three-velocity \(\mathbf v\) is represented
by

\[
p^\mu=(\gamma mc,\gamma m\mathbf v),
\]

with the invariant

\[
p^\mu p_\mu=-m^2c^2.
\]

A generic field-plus-matter four-momentum ledger separates initial, final,
external, inflow, and outflow terms. It uses a positive Euclidean residual norm
for numerical gating, while retaining the Minkowski metric for physical
invariants.

## Falsification boundary

These checks establish covariance and bounded conservation identities for the
implemented finite cases. They do not solve the Einstein field equations, do
not establish that the plasma model is complete, and do not demonstrate a
transverse plasma-wing force, net dynamic soaring, chamber equivalence, or
flight propulsion.
