# Plasma-sail validation campaign

## Status and scope

This change implements **preflight and numerical screening**, not a validated
plasma-sail force law. The native engine still uses phenomenological electric-
and magnetic-sail coefficients. No chamber measurements, flight data, GPU PIC
benchmarks, or calibrated response surfaces are supplied by this change.

The first scientific milestone is a force-vector, torque, current and power
prediction for **one explicitly specified device**, checked against independent
numerics and withheld physical measurements. More trajectory resolution alone
cannot establish that a real device produces the assumed lift.

## Run the implemented work

Requirements for these optional tools: the native C11 build and Python 3.9+.
The native executable still has no runtime Python dependency.

```sh
make all test validation-tools

# Arithmetic only: no GPUs are provisioned, purchased or used.
python3 scripts/sail_validation.py compute configs/validation/compute_3d_example.json

# Bounded sensitivity experiment, not a closed soaring cycle or mission proof.
# The output directory must not already exist.
python3 scripts/sail_validation.py sweep --output output/validation/smoke

# A longer version of the supplied seven-day synthetic experiment.
# 8 cases * 120,960 steps = 967,680 native ODE steps, not PIC steps.
python3 scripts/sail_validation.py sweep \
  --duration 604800 --step 5 --max-total-steps 1000000 \
  --timeout-seconds 300 --output output/validation/seven-day

# Separate sensitivity run restoring the existing gravity and 1PN paths.
# This does not repair the supplied non-orbital initial state or validate the device.
python3 scripts/sail_validation.py sweep --enable-gravity \
  --output output/validation/with-gravity

# Manufactured scalar data; explicitly NOT force measurements.
python3 scripts/sail_validation.py convergence \
  configs/validation/convergence_manufactured.json
```

`compute` and `convergence` return 0 for a passed **screen**, 2 for invalid input
or a file error, and 3 for an exceeded budget/failed numerical screen. `sweep`
returns 1 when any case fails or times out, records partial results, and never
marks a failed execution as a successful experiment. Reports cannot overwrite
existing files; sweeps require a new output directory. This tooling does not
launch external solvers or paid infrastructure.

The Make target exercises real CLI runs. CMake can enable the same optional
suite with `-DSPACEWIND_VALIDATION_TOOLS=ON`; the default CMake build remains
C-only. Both are connected to the existing native GitHub Actions matrix.

## 1. Conditional force-budget screen

The engine's convention is `F = rho*u^2*A*(CD*d_hat + CL*l_hat)`, with no factor
of one-half. Let `A` equal the *actual intercepted mass-flux area*, so that
`mass_flow = rho*u*A`. In a steady device frame, a cold uniform incoming stream
has speed `u`. Momentum conservation fixes the mean outgoing velocity components
as `u*(1-CD)` and `-u*CL`. Nonnegative velocity variance then gives

```text
outgoing kinetic power / incoming kinetic power >= (1-CD)^2 + CL^2
```

For a passive steady device with no additional particle/field momentum or
power flux, and no stored-energy change, a necessary condition is

```text
(1-CD)^2 + CL^2 <= 1
```

The unchanged soaring hypothesis `CD=0.4, CL=1.2` gives a lower bound of **1.8**.
Under these assumptions, it needs at least **0.8 times the incoming cold-stream
kinetic power** from additional accounting. This is **not** automatically the
spacecraft electrical-power requirement. The area, thermal distributions,
externally powered fields, emitted particles, stored energy and control-volume
boundaries can invalidate the simple assumptions and must be treated explicitly.
A high lift-to-drag ratio by itself is not disproved by this calculation.

`sw_sail_screen_cold_passive_stream()` implements this necessary condition with
only a floating-point roundoff allowance. It supports signed lift and rejects
nonfinite/overflowing results. It neither clips the forces nor changes the
coefficients to make an experiment pass.

Every native v3 mission receipt now has an additive `field_sail_validation`
section, versioned independently as `spacewind.field-sail-screen.v1`, containing:

- the enabled electric/magnetic coefficient screens, or `disabled`;
- explicit `physical_validation_established: false` and
  `flight_validation_established: false`;
- `coupled_device_power_budget: false` and the conditional screen assumptions.

A completed trajectory does not promote those fields. The existing assurance
claim tiers remain independent and unchanged. Receipts report supplied
coefficients, not a full trajectory-wide particle/field energy ledger. Even an
enabled sail with zero voltage is screened as a coefficient hypothesis; the
screen is not a measurement of its instantaneous force or power.

The native regression also checks actual electric-sail forces in a prescribed
uniform environment against `(v-U) dot F = -D*|v-U| <= 0`. It verifies the
wind-frame work identity at the force-function level, not a complete integrated
closed-cycle experiment. Lift cancels from that identity. It is possible for
this identity to pass while the separate cold-stream momentum budget fails.

## 2. Mission requirements and controls

The sweep automatically includes zero-lift, synthetic-shear-on/off and
field-sail-disabled controls. It holds the other declared conditions fixed,
uses fixed-step RK4 with a total-step cap, and disables photon, magnetic-sail
and Lorentz forces to isolate the assumed electric-sail response. Gravity is
inherited from the source configuration unless explicitly enabled.

Each case retains the generated config, CSV, native receipt and process logs.
The summary binds the original config, generated configs, executable, CSVs and
receipts with SHA-256 hashes, and records actual elapsed CPU wall time. This
wall time is not a GPU PIC throughput benchmark. Source configurations are
never modified by a sweep.

A **no-synthetic-shear** control is NOT an exactly uniform flow: the existing
Parker profile and seeded turbulence remain. The default 600-second run is a
smoke/sensitivity test and may not traverse the shear layer. Even the full
seven-day run is not a certified repeatable soaring cycle. No mission-success
threshold or power/thermal model is supplied. The reported velocity change is
not delta-v delivered by flight hardware.

Before mission promotion, predeclare a requirement envelope for useful force,
steering bandwidth, allowable power, mass, orbital initial conditions and
flow-structure lifetime. Restore gravity and uncertainty in those quantities.
Check absolute force and acceleration, not only `L/D`. Compare matched runs
against no-shear, no-lift, no-device and simpler control strategies.

## 3. Local device model: next implementation stage

Specify tether count, length, radius, layout, material/surface behavior, bias,
current closure, electron emission/collection, spin, structural support and
power electronics. Record the actual geometry and configuration content hashes.
A 50-km total tether length is not an adequate geometry specification.

Start with an appropriate electrostatic local sheath benchmark, then selected
3D device geometries. Derive `(Fx,Fy,Fz, torque, current, power, response_time)`
from particle/field/device interaction rather than prescribing `CL`. Include
open boundaries, physically specified particle injection and surface behavior.
Compute force both on the device and from a control-volume particle-and-field
momentum balance, including electromagnetic stress and changes in field
momentum. Include current-source/emitter recoil rather than mislabeling it as
solar-wind coupling.

WarpX's ion-extraction example supplies relevant building blocks (biased
embedded boundaries and injection), not a validated electric-sail input deck.
Its documented CI settings explicitly prioritize speed over production accuracy.
Do not claim an independent comparison merely by copying Spacewind's force law
into another program. A hybrid ion model is an option only where omitted
electron kinetics do not control the observable. Hoshi et al. report a case
where electron shielding changes the inferred electric-sail thrust.

Keep the current direct-DFT 3D oracle on tiny verification meshes. Its forward
transform does `N_cells^2` cell-mode contributions; the inverse is similarly
expensive. At `128^3`, one direction alone has about `4.4e12` contributions.
A scalable production solver must be independently checked against that oracle,
using FFTs only for compatible boundaries, or an appropriate open-boundary
solver. This change does **not** replace the oracle with an untested FFT or
claim to implement production open-boundary PIC.

## 4. Numerical and experimental acceptance

The `convergence` tool consumes three coarse-to-fine samples of **one scalar
observable** with a constant refinement ratio, declared absolute/relative
tolerances, estimated standard errors and distinct run identifiers. Its input
schema is demonstrated by `convergence_manufactured.json` (test data only).
It requires monotone differences resolvable above the supplied sampling error,
an observed order within the declared range, and a Richardson/GCI estimate
with safety factor at least 1.25 plus three supplied fine-level standard errors
inside the tolerance. Flat or noise-dominated outputs do not fabricate an
observed convergence order. Exact solutions require an independent analytic
oracle instead. This deliberately conservative scalar screen is not a universal
criterion for every oscillatory or stochastic discretization.

Hashes identify caller-supplied runs; this CLI does not authenticate the solver
or inspect a chamber's calibration. A digest is not evidence by itself. Sampling
errors must account for autocorrelation/effective sample size and independent
seeds. The `3*SE` allowance is a screening choice, not a distribution-free
confidence guarantee. Apply meaningful absolute tolerances near zero.

A passing scalar screen cannot set physical or flight validation to true.
Still required before promotion:

1. Separate mesh, timestep, particle-count, seed, domain-size, shape-order and
   boundary studies, including all useful force/torque/current outputs.
2. Agreement with an independent implementation and closed particle, field,
   device momentum and energy budgets. Error tolerances must be below the
   useful effect, not merely small relative to a huge background energy.
3. A preregistered chamber experiment with plasma-on/off, matched electrical and
   thermal controls, reversals where physically meaningful, dummy-device
   checks and independent replications. Freeze predictions for held-out
   operating points; do not tune on the validation set.
4. An explicit chamber-to-space similarity envelope (Debye and gyroradius
   ratios, Mach numbers, collisionality, geometry and circuit/flow timescales).
5. Coupled spacecraft mass, power, thermal, attitude, structural, deployment and
   environment uncertainty, followed by appropriate integrated/flight evidence.

A control/charging cycle should close the states that must be periodic and
account for stored-energy changes; it should NOT require the spacecraft's
velocity to return to its initial value when acceleration is the desired output.

## 5. Computation and budget contract

`compute` uses explicit, reviewable inputs, with macroparticles per cell defined
as the **total across species**. Raw particle memory and declared working
copies/field arrays are separate. Users supply effective particle updates per
GPU-second, GPU count, parallel efficiency and optional capacity/hour budgets.
A rate labeled `user_supplied_measurement` must identify its benchmark receipt
and solver revision, but is not independently authenticated by this script.

```text
particle_updates = cells * total_particles_per_cell * steps
ideal_GPU_hours = particle_updates / effective_rate / 3600
allocated_GPU_hours = ideal_GPU_hours / parallel_efficiency
elapsed_hours = allocated_GPU_hours / GPUs
```

A declared effective rate must include deposition, field solves, synchronization
and diagnostics. Aggregate memory fit does not guarantee load balance or a
working decomposition. Dollar pricing is intentionally excluded: no cloud
quote or charge is implied.

For explicit EM, the tool includes the multidimensional light-speed CFL bound.
Both modes also screen electron-plasma and maximum-particle-speed sampling.
Electrostatic mode removes the light-wave bound, not electron/sheath accuracy
requirements. Supply a maximum speed covering bias-accelerated electrons.
These ceilings do not prove numerical convergence, geometry resolution or an
adequate physical settling time.

The supplied illustrative 1-km, `256^3`, 64-total-particles/cell, 0.1-second EM
case gives about 14.77 million steps and 44,053 ideal GPU-hours at an **assumed**
`1e8` effective updates/GPU-second. It has 64 GiB raw particle storage and 132 GiB
of declared working storage (two copies plus 256 field bytes/cell). These are
arithmetic results, not measured Spacewind or WarpX performance. The example
does not establish wire/sheath resolution or a sufficient computational domain.

Benchmark a small representative local problem first. Use reduced/local models
for screening and mission response surfaces; reserve expensive 3D runs for
uncertainty and cross-code checks that could change the design decision.

## Primary references

- NASA-STD-7009, *Standard for Models and Simulations*:
  https://standards.nasa.gov/standard/NASA/NASA-STD-7009
- NASA HERTS/E-sail charged-wire/PIC investigation:
  https://ntrs.nasa.gov/citations/20180003494
- Hoshi et al. (2016), electric-sail full-3D PIC and electron shielding:
  https://angeo.copernicus.org/articles/34/845/2016/
- Larrouturou, Higgins and Greason (2022), dynamic-soaring trajectory concept:
  https://www.frontiersin.org/journals/space-technologies/articles/10.3389/frspt.2022.1017442/full
- WarpX ion-beam extraction example (versioned by the chosen external solver
  revision when used, not treated as a sail validation):
  https://warpx.readthedocs.io/en/latest/usage/examples/ion_beam_extraction/README.html
- Existing independent ledgers and claim gates:
  [assurance kernel](../assurance/README.md),
  [3D reference limitations](EM_PIC3D_REFERENCE.md),
  [chamber protocol](PLASMA_CHAMBER_PROTOCOL.md).
