# Architecture

## Design goals

1. Keep the native core independent of Python and external physics frameworks.
2. Make units and state ownership explicit.
3. Separate exact equations from constitutive approximations.
4. Permit independent tests of geometry, gravity, fields, particles, fluids,
   forces, controls, and integrators.
5. Produce replayable traces and receipts.
6. Leave clean adapters for SPICE, OMNI, WSA-ENLIL, chamber data, and higher
   fidelity solvers.

## Module map

| Module | Responsibility |
|---|---|
| `constants.h` | SI constants and solar/planetary reference values |
| `math3.c` | vectors, quaternions, 4x4 matrices |
| `rng.c` | deterministic small-state random generator |
| `metric.c` | analytic metrics, Christoffels, geodesics, 1PN |
| `gravity_grid.c` | 3-D Poisson lattice and weak-field metric |
| `nbody.c` | pairwise gravity, leapfrog, conservation diagnostics |
| `fields.c` | Parker wind/spiral, synthetic structures, gravity helpers |
| `plasma.c` | Boris, electrostatic PIC, Maxwell FDTD |
| `mhd1d.c` | finite-volume ideal MHD |
| `spacecraft.c` | photon/electric/magnetic/Lorentz forces and attitude |
| `integrator.c` | RK4, Dormand-Prince 5(4), velocity Verlet |
| `simulation.c` | config, coupled mission RHS, diagnostics, receipts |
| `cli.c` | isolated examples and mission runner |

## Ownership and allocation

- Small mathematical values are passed by value.
- Dynamically allocated solver arrays are owned by their solver structs.
- Every `*_init` function has a matching `*_destroy` function.
- Allocation failure leaves objects destroyable.
- The mission loop allocates its ODE workspace once and reuses it.
- CSV writing is streaming; mission histories are not retained in RAM.

## Determinism

- The turbulence generator takes a 64-bit seed.
- Config files are hashed byte-for-byte.
- The same compiler, flags, architecture, config bytes, and seed are expected to
  reproduce the same trace within ordinary IEEE-754 behavior.
- Cross-compiler results should be physically equivalent but are not promised to
  be bit-identical because transcendental functions and optimization choices can
  vary.

## Error behavior

The public API generally returns `bool` for operations that can fail. The CLI
uses nonzero exits for invalid config, I/O failure, solver nonconvergence, or
nonfinite state. Unknown config keys are rejected.

The engine does not use exceptions, hidden global state, or process-wide mutable
singletons.

## Coupling boundaries

### Environment to vehicle

`sw_environment_sample` is the stable interface between heliosphere models and
spacecraft forces. A future data-backed model can fill the same structure.

### Local plasma to mission model

Local PIC/MHD simulations should not be called at every mission RHS evaluation.
They should generate calibrated force/torque tables or a constrained surrogate.
The current `C_D`, `C_L`, and radius caps are explicit placeholders for that
interface.

### Gravity sources

Three gravity paths coexist because they answer different questions:

- analytic metric/geodesic for relativistic geometry;
- direct N-body for moving compact sources over orbital times;
- Poisson grid for extended weak-field mass distributions.

They are not silently summed. A caller must choose or deliberately couple them.

## Output contracts

Mission CSV files contain state, environment, force decomposition, effective
interaction radii, and energy diagnostics. JSON receipts contain run identity
and completion metadata.

CSV is used for portability and inspection. A future high-throughput binary
format should include:

- a fixed endian/version header;
- schema hash;
- unit metadata;
- chunk checksums;
- optional zstd compression;
- deterministic conversion to CSV/Arrow/Parquet outside the core.

## Threading and acceleration

Version 0.2.0 is single-threaded by design. This keeps finite checks simple.
Natural acceleration targets are:

- SIMD particle pushes;
- OpenMP/MPI domain decomposition for PIC/MHD/Poisson;
- CUDA/HIP particle and field kernels;
- Barnes-Hut or fast-multipole N-body gravity;
- multigrid Poisson solves.

Any parallel implementation should preserve a deterministic validation mode
with stable reduction ordering.
