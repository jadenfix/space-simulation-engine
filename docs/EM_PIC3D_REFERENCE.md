# Periodic 3D3V electromagnetic PIC reference

## Purpose

This module is Spacewind's first bounded electromagnetic particle-in-cell reference with three spatial dimensions and three particle-velocity components. Its purpose is not to produce an attractive plasma-wing trajectory. Its purpose is to make the most consequential numerical assumptions explicit and independently falsifiable before those assumptions are reused in an open-boundary spacecraft interaction model.

The implementation lives in:

```text
assurance/include/spacewind/em_pic3d.h
assurance/src/em_pic3d.c
assurance/src/em_pic3d_support.inc
assurance/src/em_pic3d_deposition.inc
assurance/src/em_pic3d_spectral.inc
assurance/src/em_pic3d_diagnostics.inc
assurance/src/em_pic3d_step.inc
assurance/tests/test_em_pic3d.c
```

The current model is periodic, first-order CIC, deterministic, dependency-light C11, and deliberately small enough that an independent direct-DFT oracle can run in the test path.

## State and Yee topology

For a Cartesian periodic domain with `Nx * Ny * Nz` cells, cell-centered charge density is stored at

```text
rho(i,j,k) : ((i+1/2) dx, (j+1/2) dy, (k+1/2) dz).
```

Electric field and current share face locations:

```text
Ex, Jx : (i dx,       (j+1/2)dy, (k+1/2)dz)
Ey, Jy : ((i+1/2)dx, j dy,       (k+1/2)dz)
Ez, Jz : ((i+1/2)dx, (j+1/2)dy, k dz)
```

Magnetic components use the dual edge staggering:

```text
Bx : ((i+1/2)dx, j dy,       k dz)
By : (i dx,       (j+1/2)dy, k dz)
Bz : (i dx,       j dy,       (k+1/2)dz)
```

This arrangement makes the finite-volume electric divergence use the same face currents as the discrete continuity equation. The paired forward/backward differences also make the discrete identities

```text
div(curl E) = 0
div(curl B) = 0
```

hold to floating-point roundoff on the periodic grid.

## Maxwell update

A same-time state is advanced by a symmetric split:

```text
B(n+1/2) = B(n)   - (dt/2) curl E(n)
E(n+1)   = E(n)   + c^2 dt curl B(n+1/2) - dt J(n+1/2)/epsilon0
B(n+1)   = B(n+1/2) - (dt/2) curl E(n+1)
```

The electromagnetic Courant gate is

```text
C = c dt sqrt(1/dx^2 + 1/dy^2 + 1/dz^2) <= 1.
```

The caller may impose a stricter limit. A finite result is not accepted when this gate fails.

The tests independently exercise:

- random-field `div(curl)` identities;
- propagation of Gauss law from continuity;
- preservation of magnetic divergence;
- right-going plane waves in each Cartesian direction and both relevant staggerings;
- grid refinement with approximately second-order field error;
- cyclic permutation of the three axes;
- positive-time followed by negative-time vacuum reversal;
- rejection of a deliberately incorrect magnetic-wave sign.

## Relativistic particle advance

Particles retain authoritative unwrapped positions and ordinary velocity. The velocity update uses the Higuera-Cary relativistic pusher. Internally, it advances proper velocity `u = gamma v`, performs a time-centered electric impulse and relativistic magnetic rotation, and converts back to `v` only after the update.

For static fields, the gather is corrected to the predicted path midpoint before the final push. The accepted displacement is trapezoidal:

```text
x(n+1) = x(n) + dt [v(n) + v(n+1)] / 2.
```

Tests cover:

- zero-field identity;
- positive/negative-time reversal;
- cyclic rotation covariance;
- exact proper-momentum change in a uniform electric field;
- long pure-magnetic speed conservation;
- the relativistic `E cross B` drift;
- comparison with an independently implemented high-resolution RK4 trajectory;
- second-order timestep refinement;
- rejection of luminal and nonfinite states.

## Charge deposition

Endpoint charge uses trilinear cloud-in-cell weighting at cell centers. The shape weights form a partition of unity, so total macro-particle charge is preserved independently at the initial and final endpoints.

A periodic Gauss initializer solves

```text
div E = rho / epsilon0
```

with the same discrete divergence and a direct 3D DFT Poisson solve. A nonzero periodic charge zero mode is rejected rather than silently neutralized. Requested uniform electric-field modes are preserved because they lie outside the periodic scalar potential.

## Three independent current constructions

The solver does not rely on one current implementation and then test that implementation against itself. It constructs and compares three currents.

### 1. Path-integrated face current

A raw physical-current diagnostic integrates each component's trilinear face shape along the full unwrapped particle path. The path is split at every shape-function breakpoint, including periodic crossings and complete wraps. Within each smooth segment, the product of the three linear shape factors is cubic in path parameter, so two-point Gauss-Legendre quadrature is exact for that segment.

The resulting mean current must reproduce complete unwrapped transport:

```text
mean(J) = sum_p q_p w_p [x_p(n+1)-x_p(n)] / (V_domain dt).
```

This current is used for the particle-work diagnostic because it follows the gathered path. It is not assumed to satisfy the finite-volume continuity equation exactly.

### 2. Symmetric coordinate-split charge-conserving current

The field update uses a separate current that closes discrete continuity directly.

For each particle trajectory, the algorithm evaluates all six coordinate orders:

```text
xyz, xzy, yxz, yzx, zxy, zyx.
```

Each order decomposes the 3D move into three axis-aligned moves. For an axis move, endpoint CIC charge is deposited before and after that move. The corresponding face current is reconstructed from the one-dimensional finite-volume continuity recurrence along every transverse grid line. Its undetermined periodic harmonic mode is fixed by the complete unwrapped displacement for that axis.

Every coordinate order therefore telescopes from the same initial 3D charge deposit to the same final 3D charge deposit. Averaging all six orders preserves exact discrete continuity while removing a preferred coordinate ordering.

The tests require:

- continuity to roundoff before any global correction;
- complete multi-wrap harmonic current;
- trajectory-reversal antisymmetry;
- particle-order invariance;
- integer-cell periodic translation covariance;
- cyclic-axis covariance of the full coupled step.

This construction is a bounded coordinate-split or zigzag-style reference. It is not claimed to be algebraically identical to the general Esirkepov arbitrary-shape deposition.

### 3. Independent spectral continuity oracle

A direct long-double-complex 3D DFT computes the continuity defect of the local current, solves the matching discrete periodic Poisson problem, and forms the minimum-norm longitudinal correction.

Because the coordinate-split current already satisfies continuity, the oracle correction should be near roundoff. The oracle additionally checks that its own correction is curl-free under the Yee curl operator.

The spectral current is not used to advance the fields. It is an independent, intentionally expensive verification oracle. This separation prevents a nonlocal projection from hiding an error in the current that actually drives Ampere's law.

## Midpoint particle-field coupling

The magnetic half-step is fixed first. Electric midpoint fields are then iterated:

```text
E_mid guess
  -> particle push
  -> endpoint charge
  -> path and charge-conserving currents
  -> Ampere update
  -> new E_mid = [E(n) + E(n+1)] / 2.
```

The iteration stops only when the declared relative gate closes or the bounded iteration count is exhausted. Nonconvergence blocks step promotion.

## Energy diagnostics

The module intentionally keeps two current-work identities separate:

```text
Delta K_particle  versus  integral J_path dot E_mid dt dV
Delta U_field     versus -integral J_local dot E_mid dt dV.
```

The first checks particle gathering and path integration. The second checks the Maxwell update driven by the charge-conserving current. The difference between the two work integrals exposes the divergence-free difference between the two current constructions instead of hiding it.

The total energy residual is also reported:

```text
R_E = Delta K_particle + Delta U_field.
```

The explicit scheme is not asserted to conserve energy to roundoff for arbitrary states. The tests require bounded defects and strong reduction under timestep refinement.

## Momentum and stress-energy diagnostics

Particle momentum is relativistic:

```text
p_particle = sum_p w_p gamma_p m_p v_p.
```

Field momentum is reconstructed from cell-centered interpolation of the staggered fields:

```text
p_field = integral epsilon0 E cross B dV.
```

The periodic step reports the change in total particle-plus-field momentum. A separate covariant integration test independently reconstructs cellwise electromagnetic stress-energy and requires its integrated spatial momentum to agree with the Yee diagnostic.

## Fail-closed gates

A step can be finite and still fail. Promotion requires all declared gates:

- electromagnetic CFL;
- particle displacement in cells;
- global charge closure;
- local finite-volume continuity;
- initial and final Gauss law;
- initial and final magnetic divergence;
- bounded path-current versus local-current discrepancy;
- near-roundoff agreement with the independent spectral oracle;
- exact harmonic mean current;
- midpoint nonlinear convergence;
- subluminal particle motion;
- particle-work, field-work, total-energy, and total-momentum tolerances.

Adversarial tests deliberately trigger CFL rejection, over-strict current-promotion rejection, nonneutral periodic Gauss rejection, magnetic-divergence corruption, Gauss corruption, current zero-mode corruption, wrong wave polarization, luminal input, and nonfinite input.

## Validation matrix

The registered test is designed for:

- Linux GCC warning-as-error;
- Linux Clang warning-as-error;
- macOS ARM64 Clang;
- AddressSanitizer and UndefinedBehaviorSanitizer;
- Clang static analyzer;
- GCC `-fanalyzer`;
- deterministic receipt replay;
- standalone Make and CMake/CTest paths.

The deterministic core suite currently performs more than seventy-eight thousand bounded checks and emits:

```text
assurance/output/em_pic3d_reference_receipt.json
assurance/output/em_pic3d_step_receipt.json
```

## Scientific boundary

This module materially improves the numerical evidence substrate, but it does not establish a plasma-wing force or propulsion system.

Remaining fidelity blockers include:

- a production local Esirkepov or independently equivalent higher-order deposition;
- higher-order particle shapes and matching gather/scatter operators;
- open electromagnetic boundaries and PML or characteristic boundary treatment;
- particle injection, absorption, reflection, and secondary emission;
- collisions, ionization, recombination, and radiation;
- material surfaces, sheath formation, field emission, arcing, and spacecraft charging;
- adaptive or nonuniform meshes;
- scalable FFT or domain-decomposed solvers;
- particle-count, grid, timestep, shape-order, domain-size, and boundary-condition refinement;
- independent implementation agreement;
- chamber-calibrated force and torque closure.

## Method references

- K. S. Yee, *Numerical solution of initial boundary value problems involving Maxwell's equations in isotropic media*, IEEE Transactions on Antennas and Propagation 14 (1966), DOI `10.1109/TAP.1966.1138693`.
- J. Villasenor and O. Buneman, *Rigorous charge conservation for local electromagnetic field solvers*, Computer Physics Communications 69 (1992), DOI `10.1016/0010-4655(92)90169-Y`.
- T. Umeda et al., *A new charge conservation method in electromagnetic particle-in-cell simulations*, Computer Physics Communications 156 (2003), DOI `10.1016/S0010-4655(03)00437-5`.
- T. Zh. Esirkepov, *Exact charge conservation scheme for Particle-in-Cell simulation with an arbitrary form-factor*, Computer Physics Communications 135 (2001), DOI `10.1016/S0010-4655(00)00228-9`.
- A. V. Higuera and J. R. Cary, *Structure-preserving second-order integration of relativistic charged particle trajectories in electromagnetic fields*, Physics of Plasmas 24 (2017), DOI `10.1063/1.4979989`.
