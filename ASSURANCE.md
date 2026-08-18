# Physical assurance and falsification

Spacewind's independent assurance kernel lives under [`assurance/`](assurance/).
It checks dimensions, intervals, convergence, control-volume conservation,
discrete field-and-particle ledgers, charge-compatible particle transport,
relativistic and electromagnetic invariants, heliospheric flux identities,
force upper bounds, closed-cycle energy accounting, chamber force reversals,
independent replication, and fail-closed scientific claim promotion.

The detailed contract is documented in:

- [`assurance/README.md`](assurance/README.md)
- [`docs/ASSURANCE_AND_FALSIFICATION.md`](docs/ASSURANCE_AND_FALSIFICATION.md)
- [`docs/FIELD_AND_PARTICLE_LEDGER.md`](docs/FIELD_AND_PARTICLE_LEDGER.md)
- [`docs/CHARGE_CONSERVING_DEPOSITION.md`](docs/CHARGE_CONSERVING_DEPOSITION.md)
- [`docs/PLASMA_CHAMBER_PROTOCOL.md`](docs/PLASMA_CHAMBER_PROTOCOL.md)
- [`docs/NUMERICAL_VERIFICATION.md`](docs/NUMERICAL_VERIFICATION.md)

The current repository is a simulation and verification environment. It does
not claim that a flight-qualified plasma-wing propulsion system exists.
