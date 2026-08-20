# Physical assurance and falsification

Spacewind's independent assurance kernel lives under [`assurance/`](assurance/).
It checks dimensions, intervals, convergence, control-volume conservation,
discrete field-and-particle ledgers, charge-compatible particle transport,
periodic spectral Poisson and Gauss-law identities, exact periodic Maxwell
references, bounded 1D3V and 3D3V electromagnetic PIC references, covariant
stress-energy and total/incremental four-momentum identities, open-boundary
Poynting/Maxwell-stress transport, relativistic and electromagnetic invariants,
heliospheric flux identities, force upper bounds, closed-cycle energy accounting,
chamber force reversals, independent replication, and fail-closed scientific
claim promotion.

The 3D3V reference now separates three current constructions: an unwrapped
path-integrated current for particle work, a six-permutation coordinate-split
charge-conserving current that drives Ampere's law, and an independent direct-DFT
minimum-norm continuity oracle that must leave the field-driving current unchanged
to near roundoff. Same-timestep field and matter changes are additionally mapped
into an incremental four-momentum ledger so constant rest energy cannot hide a
small numerical residual.

The detailed contract is documented in:

- [`assurance/README.md`](assurance/README.md)
- [`docs/ASSURANCE_AND_FALSIFICATION.md`](docs/ASSURANCE_AND_FALSIFICATION.md)
- [`docs/FIELD_AND_PARTICLE_LEDGER.md`](docs/FIELD_AND_PARTICLE_LEDGER.md)
- [`docs/CHARGE_CONSERVING_DEPOSITION.md`](docs/CHARGE_CONSERVING_DEPOSITION.md)
- [`docs/PIC_TRANSPORT_INTEGRATION.md`](docs/PIC_TRANSPORT_INTEGRATION.md)
- [`docs/SPECTRAL_POISSON_AUDIT.md`](docs/SPECTRAL_POISSON_AUDIT.md)
- [`docs/MAXWELL1D_SPECTRAL_REFERENCE.md`](docs/MAXWELL1D_SPECTRAL_REFERENCE.md)
- [`docs/EM_PIC1D_REFERENCE.md`](docs/EM_PIC1D_REFERENCE.md)
- [`docs/EM_PIC1D_COVARIANT_INTEGRATION.md`](docs/EM_PIC1D_COVARIANT_INTEGRATION.md)
- [`docs/EM_PIC3D_REFERENCE.md`](docs/EM_PIC3D_REFERENCE.md)
- [`docs/EM_PIC3D_COVARIANT_INTEGRATION.md`](docs/EM_PIC3D_COVARIANT_INTEGRATION.md)
- [`docs/FOUR_MOMENTUM_DELTA_LEDGER.md`](docs/FOUR_MOMENTUM_DELTA_LEDGER.md)
- [`docs/OPEN_BOUNDARY_FOUR_MOMENTUM.md`](docs/OPEN_BOUNDARY_FOUR_MOMENTUM.md)
- [`docs/FDTD_SPECTRAL_CROSSCHECK.md`](docs/FDTD_SPECTRAL_CROSSCHECK.md)
- [`docs/COVARIANT_STRESS_ENERGY.md`](docs/COVARIANT_STRESS_ENERGY.md)
- [`docs/MAXWELL_COVARIANT_INTEGRATION.md`](docs/MAXWELL_COVARIANT_INTEGRATION.md)
- [`docs/PLASMA_CHAMBER_PROTOCOL.md`](docs/PLASMA_CHAMBER_PROTOCOL.md)
- [`docs/NUMERICAL_VERIFICATION.md`](docs/NUMERICAL_VERIFICATION.md)

The current repository is a simulation and verification environment. It does
not claim that a flight-qualified plasma-wing propulsion system exists.
