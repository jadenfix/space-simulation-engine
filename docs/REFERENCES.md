# Primary References and Data Sources

The code is independently implemented. These references define or motivate the
models and validation problems.

## Heliophysics

- E. N. Parker, “Dynamics of the Interplanetary Gas and Magnetic Fields,”
  *The Astrophysical Journal* 128, 664 (1958).
  DOI: https://doi.org/10.1086/146579
- NASA CCMC, WSA-ENLIL model documentation:
  https://ccmc.gsfc.nasa.gov/models/WSA-Enlil-at-SWPC~3/
- NASA Space Physics Data Facility, OMNI datasets:
  https://data.nasa.gov/dataset/omni-daily-data-set
- NASA NAIF, SPICE Toolkit documentation:
  https://naif.jpl.nasa.gov/naif/

## Numerical electromagnetics and plasma

- J. P. Boris, “Relativistic Plasma Simulation—Optimization of a Hybrid Code,”
  Proceedings of the Fourth Conference on Numerical Simulation of Plasmas
  (1970), pp. 3-67.
- K. S. Yee, “Numerical Solution of Initial Boundary Value Problems Involving
  Maxwell's Equations in Isotropic Media,” *IEEE Transactions on Antennas and
  Propagation* 14(3), 302-307 (1966).
  DOI: https://doi.org/10.1109/TAP.1966.1138693
- M. Brio and C. C. Wu, “An Upwind Differencing Scheme for the Equations of
  Ideal Magnetohydrodynamics,” *Journal of Computational Physics* 75(2),
  400-422 (1988).
  DOI: https://doi.org/10.1016/0021-9991(88)90120-9

## Field sails and dynamic soaring

- P. Janhunen, “Electric Sail for Spacecraft Propulsion,” *Journal of
  Propulsion and Power* 20(4), 763-764 (2004).
  DOI: https://doi.org/10.2514/1.8580
- P. Janhunen and A. Sandroos, “Simulation Study of Solar Wind Push on a
  Charged Wire,” *Annales Geophysicae* 25, 755-767 (2007).
  DOI: https://doi.org/10.5194/angeo-25-755-2007
- D. G. Andrews and R. M. Zubrin, “Use of Magnetic Sails for Advanced
  Exploration Missions,” NASA NTRS 19910012840 (1990).
  https://ntrs.nasa.gov/citations/19910012840
- M. N. Larrouturou, A. J. Higgins, and J. K. Greason, “Dynamic Soaring as a
  Means to Exceed the Solar Wind Speed,” *Frontiers in Space Technologies* 3
  (2022).
  DOI: https://doi.org/10.3389/frspt.2022.1017442

## Interpretation

Citing a concept paper does not validate the surrogate force law in this code.
The references motivate models and test problems; Spacewind's coefficients and
caps remain explicit calibration inputs until tied to experiment or a verified
higher-fidelity calculation.
