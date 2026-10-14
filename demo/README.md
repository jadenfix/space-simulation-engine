# Periodic plasma benchmark demo

The GIF and MP4 render actual numerical convergence and conservation outputs from
source commit `d0d546ee5b519a739a45989c5a28524bf3a91c55`.
These are explanatory animations of executed C results, not screen recordings,
particle trajectories, laboratory measurements or device validation.

## Reproduce

On Linux with Python 3/Pillow, DejaVu fonts, a C11 compiler, make and ffmpeg:

```sh
python3 demo/render.py
```

The renderer verifies 92 source/dependency Git blob hashes, executes the full
20-binary assurance suite (`make -C assurance -j2 test`), then renders the freshly
emitted `em_pic3d_plasma_modes_receipt.json`. All six cases must pass.
`plasma-results.json` contains raw numeric output; `core-run.txt` records executed
command output. Rendering refuses source drift under the recorded revision.
Change the two font paths for non-Linux rendering. Media are 960×540, 16 seconds
at 12fps; MP4 is H.264/yuv420p with fast-start. Encoding is single-threaded and
uses a two-pass palette to limit memory use in shared environments.

## Scenario and plotted quantities

The benchmark uses a periodic 8×1×1 grid with two initially stationary charged
macro-particles per cell, a uniform 1000 V/m electric field and opposite charges.
Base density is 1.2468181534416748e15 m⁻³. Six cases vary timestep, density or
particle masses. This bounded, effectively one-dimensional mode is embedded in
the 3D3V reference solver; it is not a multidimensional device simulation.

Frequency error bars are `abs(measured_omega/continuum_omega - 1)` from coarse,
medium and fine runs. Second-order refinement values come unchanged from the
actual receipt. Conservation values are maximum relative residuals of the
pair-plasma case, not measured physical accuracy or speed.

Independent verification: full assurance and separately compiled UBSan suites
passed. Current-main hosted GCC/Clang/macOS strict builds, ASan+UBSan,
deterministic replay, analyzers and native regression passed in run 38328636603.
No source behavior, tolerance, iteration budget or CI gate changes in this PR.

## Boundaries

Synthetic initial conditions; no live space-weather integration, calibrated
propulsion, laboratory measurements, flight validation or performance claim.
See root README and assurance documentation for all scientific limits.

The October 14 refresh preserves the same numerical solver and media, and re-executed
all 20 assurance binaries after adding the 1D included-fragment build dependency.
All four scientific hosted workflows and workflow lint passed at the recorded source.
