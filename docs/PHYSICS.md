# Physics and Numerical Model

## 1. Conventions

Spacewind uses SI units throughout. Four-coordinate APIs use

\[
x^\mu=(ct,x^1,x^2,x^3)
\]

and metric signature \((-+++ )\). Mission-scale spacecraft propagation uses
heliocentric Cartesian coordinates. Analytic Schwarzschild and Kerr geodesic
examples use their conventional spherical/Boyer-Lindquist coordinates.

The engine intentionally separates:

- spacetime geometry;
- matter-generated gravity;
- electromagnetic fields and plasma;
- spacecraft constitutive/interaction models;
- numerical integration.

No code path assumes that spacetime itself supplies aerodynamic lift.

## 2. Analytic spacetime metrics

### Minkowski

\[
ds^2=-d(ct)^2+dx^2+dy^2+dz^2.
\]

### Schwarzschild

With \(r_s=2GM/c^2\),

\[
ds^2=-\left(1-\frac{r_s}{r}\right)d(ct)^2
+\left(1-\frac{r_s}{r}\right)^{-1}dr^2
+r^2d\theta^2+r^2\sin^2\theta\,d\phi^2.
\]

### Kerr

The Kerr implementation uses Boyer-Lindquist coordinates with geometric mass
\(M_g=GM/c^2\), dimensionless spin clipped to \([-1,1]\), and
\(a=\chi M_g\). It includes the frame-dragging \(g_{t\phi}\) term. Setting
\(\chi=0\) is tested against Schwarzschild in the same coordinates.

### Flat FLRW/de Sitter

\[
ds^2=-d(ct)^2+a(t)^2(dx^2+dy^2+dz^2),
\qquad a(t)=a_0e^{Ht}.
\]

Setting \(H=0\) and \(a_0=1\) is tested against Minkowski.

## 3. Christoffel symbols and geodesics

The metric derivatives are evaluated with centered finite differences. The
connection is then assembled from

\[
\Gamma^\mu_{\alpha\beta}
=\frac{1}{2}g^{\mu\nu}
\left(
\partial_\alpha g_{\nu\beta}
+\partial_\beta g_{\nu\alpha}
-\partial_\nu g_{\alpha\beta}
\right).
\]

Geodesics satisfy

\[
\frac{d^2x^\mu}{d\lambda^2}
+\Gamma^\mu_{\alpha\beta}
\frac{dx^\alpha}{d\lambda}
\frac{dx^\beta}{d\lambda}=0
\]

and are integrated as an eight-dimensional first-order system with RK4.
Numerical metric differentiation is general but expensive and
coordinate-sensitive. It is not a substitute for analytic derivatives in a
high-precision black-hole code.

## 4. Weak-field gravity from matter

### Direct N-body system

For body \(i\),

\[
\mathbf a_i
=G\sum_{j\ne i}m_j
\frac{\mathbf r_j-\mathbf r_i}
{\left(|\mathbf r_j-\mathbf r_i|^2+\epsilon^2\right)^{3/2}}.
\]

Pairs are accumulated symmetrically. Time integration uses kick-drift-kick
leapfrog, a second-order symplectic method. Optional overlap merging preserves
mass and linear momentum and uses volume-additive radii.

### Three-dimensional Poisson lattice

In the Newtonian/weak-field regime,

\[
\nabla^2\Phi=4\pi G\rho.
\]

Spacewind discretizes the Laplacian on a uniform Cartesian grid with the
seven-point stencil and solves the resulting system with red-black successive
over-relaxation. Cloud-in-cell mass deposition and trilinear potential sampling
are provided.

Acceleration is

\[
\mathbf a=-\nabla\Phi.
\]

The corresponding static first-order weak-field metric is exposed as

\[
g_{00}=-\left(1+\frac{2\Phi}{c^2}\right),
\qquad
 g_{ij}=\left(1-\frac{2\Phi}{c^2}\right)\delta_{ij}.
\]

This is not a nonlinear numerical-relativity solver. It is valid only when
\(|\Phi|/c^2\ll1\) and source velocities/stresses do not require gravitomagnetic
or radiative terms.

### First post-Newtonian correction

For a test particle around one central mass, the implemented Schwarzschild 1PN
correction is

\[
\mathbf a_{1PN}
=\frac{\mu}{c^2r^3}
\left[
\left(\frac{4\mu}{r}-v^2\right)\mathbf r
+4(\mathbf r\cdot\mathbf v)\mathbf v
\right].
\]

It is added to Newtonian gravity when enabled.

## 5. Parker solar wind

The baseline environment solves the isothermal Parker wind. With isothermal
sound speed \(a\) and critical radius

\[
r_c=\frac{GM_\odot}{2a^2},
\]

the transonic solution satisfies

\[
y-\ln y
=4\ln\left(\frac{r}{r_c}\right)
+\frac{4r_c}{r}-3,
\qquad y=\left(\frac{u}{a}\right)^2.
\]

The subsonic and supersonic branches are selected on opposite sides of the
critical point and solved by safeguarded bisection. Number density follows
steady radial mass-flux conservation:

\[
n(r)u(r)r^2=n(r_0)u(r_0)r_0^2.
\]

The implementation is an isothermal baseline. It omits heat conduction,
wave-pressure acceleration, multi-fluid temperatures, pickup ions, shocks, and
kinetic departures unless represented separately.

## 6. Parker spiral and plasma diagnostics

The radial magnetic field scales as

\[
B_r(r)=B_{r0}\left(\frac{r_0}{r}\right)^2.
\]

Outside a source surface, solar rotation generates an azimuthal component of
the form

\[
B_\phi
=-\frac{\Omega_\odot(r-r_s)\sin\theta}{u_r}B_r.
\]

The ideal-MHD electric field is

\[
\mathbf E=-\mathbf u\times\mathbf B.
\]

The environment also computes:

\[
P_{dyn}=\rho |\mathbf u|^2,
\qquad
P_B=\frac{B^2}{2\mu_0},
\]

\[
\lambda_D
=\sqrt{\frac{\epsilon_0 k_BT_e}{n_e e^2}},
\qquad
v_A=\frac{B}{\sqrt{\mu_0\rho}},
\]

plus thermal, sound, and fast-magnetosonic quantities.

## 7. Streams, turbulence, shear, and CME shells

These additions are explicit controlled models:

- Fast/slow streams use rotating longitudinal and latitude modulation.
- Turbulence is a deterministic sum of transverse Fourier modes with seeded
  phases and Alfvénic velocity perturbations.
- The shear layer uses a smooth hyperbolic-tangent transition across a plane.
- The CME uses a Gaussian radial shell that scales density, field, and radial
  speed.

They are useful for controller and sensitivity experiments. They are not live
space-weather reconstructions. Future adapters should ingest NASA OMNI and
WSA-ENLIL/CCMC data while preserving the analytic model as a controlled
baseline.

## 8. Solar-radiation pressure

For photon flux \(I\), sail normal \(\hat n\), radial photon direction
\(\hat s\), and \(\cos\alpha=|\hat n\cdot\hat s|\), the implementation uses

\[
\mathbf F_{abs}
=\frac{I}{c}A\alpha_{abs}\cos\alpha\,\hat s,
\]

\[
\mathbf F_{refl}
=\frac{2I}{c}A\rho_{spec}\cos^2\alpha\,\hat n.
\]

No thermal re-radiation or membrane deformation model is included yet.

## 9. Electric-sail surrogate

The mission-scale electric-sail model estimates a sheath radius from Debye
length and voltage-to-electron-thermal-energy ratio:

\[
r_{eff}
=\lambda_D
\sqrt{\ln\left(1+\frac{e|V|}{k_BT_e}\right)},
\]

then clamps it to an explicit user maximum. The interaction area is

\[
A_{eff}=2L_{tether}r_{eff}.
\]

Force is a pressure-law surrogate:

\[
\mathbf F
=P_{dyn}A_{eff}
\left(C_D\hat d+C_L\hat l\right).
\]

This is deliberately parameterized. It is not a self-consistent charged-tether
PIC solution and it does not establish a real transverse lift coefficient.
Local PIC or chamber data should replace the coefficients with a response
surface over density, velocity distribution, field orientation, voltage,
geometry, and spacecraft charge.

## 10. Magnetic-sail surrogate

A pressure-balance radius is estimated from dipole moment \(m\):

\[
r_m
=\left(
\frac{\mu_0m^2}{8\pi^2P_{dyn}}
\right)^{1/6}.
\]

The effective area is \(\pi r_m^2\), followed by the same explicit drag/lift
coefficient model. The radius is clamped to prevent an unconstrained dipole
formula from silently generating an arbitrarily large virtual sail.

## 11. Charged-particle dynamics

The Lorentz equation is

\[
\frac{d\mathbf p}{dt}
=q(\mathbf E+\mathbf v\times\mathbf B).
\]

The Boris pusher splits electric acceleration around a magnetic rotation. In a
static magnetic field with zero electric field, the method preserves particle
speed to floating-point/numerical-rotation error. A relativistic variant evolves
momentum and reconstructs velocity with the Lorentz factor.

## 12. Electrostatic particle-in-cell

The one-dimensional PIC path performs:

1. cloud-in-cell charge deposition;
2. subtraction of the domain-average background charge;
3. periodic spectral solution of Poisson's equation;
4. electric-field reconstruction;
5. field interpolation to particles;
6. kick-drift-kick particle stepping.

The direct spectral transform is intentionally dependency-free and
\(O(N_g^2)\). It is transparent but not the production path for large grids.

## 13. Maxwell FDTD

The one-dimensional field solver advances a transverse electromagnetic wave on
a staggered Yee-like grid. Its time step must satisfy the one-dimensional CFL
condition

\[
c\Delta t/\Delta x\le1.
\]

Periodic boundaries make it a clean finite check for wave speed and field
energy, not a complete antenna or plasma-wave solver.

## 14. Ideal MHD

The MHD state contains density, three momenta, transverse magnetic components,
and total energy, with constant longitudinal magnetic field. Interface fluxes
use the Harten-Lax-van Leer approximation bounded by local fast-magnetosonic
signal speeds. The solver is first-order in space and time.

The Brio-Wu shock tube is included because it exercises fast waves, slow waves,
contacts, and compound structures. Positivity floors are explicit numerical
safeguards and can perturb energy conservation when activated.

## 15. Spacecraft dynamics and control

The mission state has 14 scalars:

- 3 position;
- 3 velocity;
- 4 attitude quaternion;
- 3 body angular velocity;
- 1 accumulated proper time.

Rigid-body dynamics use Euler's equation with diagonal inertia. The attitude
controller is a bounded-complexity quaternion/vector PD controller. The
`dynamic_soar` mode reverses desired crossing direction at configurable shear
coordinates and requests a lift direction transverse to relative flow.

The control law is a research baseline. It is not a globally optimal policy and
has no guarantee of extracting net energy from an arbitrary measured plasma
boundary.

## 16. Multiscale coupling rule

A full interplanetary PIC calculation is computationally intractable. Spacewind
therefore requires uncertain microphysics to be surfaced as calibrated
parameters rather than hidden:

\[
(\mathbf F,\boldsymbol\tau)
=f(n,T_i,T_e,\mathbf v,\mathbf B,\mathbf E,V,m,
\text{geometry},\text{attitude}).
\]

High-fidelity local calculations or experiments should populate this response
map. The mission integrator then evaluates it along a trajectory. Uncertainty
bands should be propagated over the response map before any performance claim
is made.
