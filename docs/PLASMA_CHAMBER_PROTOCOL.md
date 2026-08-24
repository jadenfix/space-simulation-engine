# Plasma-wing chamber protocol

## Objective

The first decisive physical milestone is not high speed. It is a repeatable,
commandable force component perpendicular to a calibrated plasma flow whose
sign reverses with the commanded field orientation and survives null,
common-mode, drift, and momentum-budget checks.

The protocol below is designed to distinguish controlled transverse momentum
exchange from ordinary drag, thermal drift, electrostatic attraction to the
facility, magnetic coupling to chamber hardware, sensor bias, cable force, and
analysis flexibility.

## 1. Declare the estimand before data collection

Primary estimand:

\[
F_{rev}=\frac{E[F\mid +u]-E[F\mid -u]}{2},
\]

where `+u` and `-u` are equal-magnitude opposite steering commands.

Primary null:

\[
H_0:F_{rev}=0.
\]

Predeclare:

- minimum scientifically relevant force;
- maximum common-mode bias;
- maximum reversal asymmetry;
- maximum acceptable zero-control drift;
- maximum lag-one autocorrelation;
- significance level;
- stopping rule;
- excluded-run rules;
- all calibration corrections.

Do not select these thresholds after observing the result.

## 2. Facility characterization

Measure and archive before each run block:

- ion and electron density;
- ion and electron temperature;
- ion species and charge state;
- full velocity distribution or its declared approximation;
- magnetic-field vector throughout the test volume;
- electric potential throughout the test volume;
- neutral pressure and composition;
- beam nonuniformity and divergence;
- chamber-wall potential;
- vibration spectrum;
- thermal gradients;
- force-sensor transfer function and noise spectrum.

Every measurement needs units, calibration date, calibration uncertainty,
sampling cadence, location, and raw-data provenance.

## 3. Similarity vector

Publish the chamber and target-space values for at least:

\[
\frac{eV}{k_BT_e},
\quad
\frac{a}{\lambda_D},
\quad
\frac{\rho_i}{L},
\quad
\frac{d_i}{L},
\]

\[
M_s,
\quad
M_A,
\quad
M_f,
\quad
\beta,
\quad
Kn,
\]

\[
\frac{\tau_{circuit}}{\tau_{flow}},
\quad
\frac{P_{available}}{P_{collection}},
\quad
\lambda_D.
\]

A chamber result outside the target neighborhood must be labeled an
extrapolation. Matching one attractive ratio does not establish global
similarity.

## 4. Mechanical measurement

Use at least two independent force-estimation paths when practical:

1. a calibrated torsion balance or flexure;
2. momentum-flux reconstruction from downstream plasma diagnostics.

The mechanical path should record all axes, not only the expected lift axis.
Resolve:

- axial drag;
- commanded transverse force;
- torque;
- structural resonances;
- cable and feedthrough loads;
- thermal expansion;
- electrostatic interaction with grounded surfaces.

Sensor bandwidth must exceed the command schedule bandwidth while remaining
below dominant mechanical resonance unless a validated dynamic model is used.

## 5. Randomized command schedule

Use balanced, concealed blocks containing:

- positive steering command;
- negative steering command;
- zero steering command;
- high voltage with steering disabled;
- steering electronics active with high voltage disabled;
- plasma off;
- dummy electrode or electrically equivalent noninteracting load.

Randomize order within constraints that prevent unsafe switching. The analysis
should not know command labels until calibration and exclusion decisions are
frozen.

Interleave zero controls so thermal and mechanical drift are observable rather
than inferred from a separate day.

## 6. Reversal requirements

A valid controlled transverse-force result must show:

- positive-command mean on the declared positive axis;
- negative-command mean on the declared negative axis;
- a zero-command mean consistent with the bias limit;
- similar absolute response under opposite commands;
- a paired reversed-force interval above the minimum relevant force;
- a sign-flip randomization result below the predeclared alpha;
- acceptable autocorrelation and zero-control drift.

A force that remains positive under both steering polarities is common-mode
force, not demonstrated lift.

## 7. Momentum accounting

Construct a control volume around the device. Compare mechanical force with:

\[
\mathbf F_{device}=-\oint
\left(ho\mathbf vv+\mathbf T_{EM}+p\mathbf I\right)\cdot d\mathbf A
-\frac{d}{dt}\int\mathbf g_{field}\,dV.
\]

At minimum, record:

- incoming ion momentum flux;
- outgoing ion momentum flux;
- absorbed-particle momentum;
- electron contribution or explicit bound;
- Maxwell stress;
- field momentum storage;
- force on nearby electrodes and chamber walls.

Do not add particle-impact force and Maxwell-stress force if they describe the
same control-volume transfer.

## 8. Electrical accounting

Record synchronized:

- tether/electrode voltage;
- supply current;
- collected ion current;
- collected electron current;
- emitter or plasma-contactor current;
- leakage current;
- arc events;
- stored capacitive energy;
- bus power and conversion loss.

Audit charge:

\[
C(V_f-V_i)=\int(I_s+I_i-I_e-I_{leak})dt.
\]

Audit electrical energy:

\[
\Delta\left(\frac12CV^2\right)=
E_{bus}+E_i-E_e-E_{leak}-E_{conversion}.
\]

A force obtained only while the voltage collapses is a transient discharge,
not sustained propulsion.

## 9. Thermal and structural controls

Measure temperature at the electrode, support, sensor, feedthrough, and nearby
reference structure. Run a thermal dummy with equivalent electrical
heating but no intended plasma interaction.

Track:

- creep;
- outgassing;
- tether elongation;
- electrostatic deflection;
- magnetic torque on conductors;
- Lorentz force on supply leads;
- vibration and acoustic coupling;
- plasma-induced erosion.

## 10. Replication ladder

Promotion requires progressively stronger evidence:

1. repeated blocks on one apparatus;
2. rebuild with new hardware on the same facility;
3. independent operator and analysis replay;
4. second facility with independently calibrated diagnostics;
5. free-flying or externally isolated demonstration;
6. flight experiment in the intended plasma regime.

Two analyses of the same raw run are not independent physical replications.

## 11. Required release package

A credible result should release or escrow:

- preregistration;
- command schedule generation code;
- raw force and diagnostic streams;
- calibration data;
- exclusion log;
- exact analysis code and environment;
- similarity vector;
- conservation ledgers;
- uncertainty budget;
- negative and null results;
- machine-readable claim receipt.

## 12. Stop conditions

The plasma-wing thesis should be narrowed or rejected for a tested design when:

- transverse force does not reverse;
- common-mode bias dominates the reversible component;
- force disappears after thermal or cable controls;
- momentum accounting cannot support the measured force;
- power collection collapses the field;
- the effect fails rebuild or independent-facility replication;
- the complete hardware mass and power remove mission-level advantage.

A failed design is useful evidence. It should update the model response surface
rather than be hidden by changing the estimand after the fact.
