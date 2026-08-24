# Yee FDTD versus exact spectral Maxwell cross-check

Spacewind now tests the root periodic Yee update against the independent
characteristic-field Maxwell reference rather than validating each method only
with its own implementation.

The root solver stores transverse electric field at integer time and magnetic
field at the interleaved half time. For a right-moving wave, the initial data
therefore uses

\[
E_i^0=f(x_i),
\]

\[
B_{i+1/2}^{-1/2}=\frac{1}{c}
 f\left(x_{i+1/2}+\frac{c\Delta t}{2}\right).
\]

The test advances one-half light-crossing time with Courant number 0.5 on 33,
65, and 129 cells. It compares the root electric field with both the analytic
traveling wave and the direct-DFT Maxwell oracle at the exact final time of
each run.

Required checks include:

- spectral-reference error near floating-point roundoff;
- monotonically decreasing Yee error;
- observed coarse/medium and medium/fine order between 1.75 and 2.25;
- improving staggered field-energy diagnostic under refinement;
- rejection of a deliberately reversed magnetic-wave sign.

This isolates spatial/time staggering and numerical dispersion from plasma
particle noise. It is a bounded one-polarization periodic wave test, not a
proof of the complete electromagnetic PIC solver, chamber force, or
propulsion.
