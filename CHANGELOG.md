# Changelog

Changes of each release of rest-for-physics-toolkit (rest-framework). The pull request that sets the
`VERSION` of a release moves the entries of "Unreleased" under its version (CONTRIBUTING.md, Versions and
releases); the notes of the GitHub release list its pull requests.

## Unreleased

- Continuous integration on the environment pinned by each commit: `.github/rest-env` holds the tag of the
  `rest-docker` repository (image `rest-env` and `rest/Dockerfile`), so a release is validated again with
  its own environment.
- restGarfield builds with the Garfield++ versions after 4.0 (e.g. that of the LCG views).
- `scripts/install.sh`: Garfield++ 4.0 builds on the LCG views (C++20) and is the one used even when the
  view has another (`--lcg-garfield` to use that of the view); the packages of the system checked by the
  version available to install (Debian 13); the script is read whole before running. Source mode checked on
  Ubuntu 24.04, Debian 13, AlmaLinux 9 and macOS 15: ROOT 6.36 with the recent CMake versions and without
  Davix, Geant4 11.0 with Xcode 16 (zlib of the system, g4tools typo patched), packages of dnf by what
  they provide, bash 3.2 of macOS.
- restG4 example test of the full chain decay: the ground state of Pa234 is a rare branch.
- restG4: the histogram files of the generator (TH1D, TH2D and cosmic distributions) are found when restG4 is
  launched from another directory than its configuration (the generator of each thread read them again
  relative to the launch directory), and reading them no longer leaves the output file without the geometry
  (the current ROOT directory is kept). The tracks kept of each event (removeUnwantedTracks) are printed only
  with verbose debug.
- restGarfield drift responses in CST field maps: the nodes with the field 0 (out of the solved volume, in
  the conductors) are out of the gas; the grid of start points by its spacing (`gridStep`, instead of
  `gridNodes`) and, by default, over the readout channels plus `margin`; `--jobs` splits the start points;
  the time per electron printed. The response file has a summary per start point (survival, mean
  displacement, drift time, spread) and the macro `REST_Garfield_DrawDriftResponse` draws it; counts stored
  as 16-bit integers with LZMA compression (`electrons` at most 65535). Response files computed before must
  be computed again.
- Gaussian sigma of tracks and hits (`gaussFit`): the drift coordinate is fitted as continuous, with a
  resolution per hit (`gaussResolution`, default [0, 0, 0.1 mm]: the strips along X and Y), instead of
  strips as wide as the smallest distance between hits (slow on long tracks); three positions with more
  than 1% of the energy of the largest one (not of the total: long tracks fit); the strip fit evaluates
  only the strips near the centre (faster, same results).
- TRestRawSignalPeaksFinderProcess: all the peaks of raw signals over a threshold in sigmas of the
  baseline of each signal, and their coincidence windows, so channels with different gains (e.g. vetoes)
  need no threshold in ADC nor a calibration.

## 3.0.0 (2026-11-24)

First release of rest-for-physics-toolkit, the new version of REST-for-Physics (whose last release is
2.4.3). It does not read the files of REST-for-Physics; MIGRATION_GUIDE.md describes the move.

- YAML configurations with registered fields (unknown keys are errors), `parameters` and `--set`,
  `include`; event classes with flat branches and a format version per class; metadata stored as YAML.
- Libraries: tools and core (runs, processes, pipelines of restManager, analysis tree read with
  RDataFrame, cuts and plots), analysis, raw, detector (readouts, media and their charge, light and
  wavelength shifter responses, drift, photosensors), track, geant4, connectors.
- Datasets (`TRestDataSet`) and dataset processes: calibration, gain map, log odds, TMVA, signal
  selection, combine, Tardis.
- Executables: restManager, restRoot, restG4 (with CRY), restGarfield (transport and gain tables, CST
  field maps, drift responses), restDAQ (FEMDAQ).
- Version of the code (commit and clean tree) stored in the output files and given by `--version`.
- Tests: unit tests per component, restG4 examples and pipelines on data and simulations of
  rest-common-data, with results independent of the platform (Linux and macOS).
- Scripts: `scripts/install.sh` (dependencies of the system, an LCG view or built from sources) and
  `restSubmit` (restG4 and restManager jobs on HTCondor or Slurm).
- Docker image `ghcr.io/rest-for-physics-toolkit/rest:v3.0.0` and documentation of each version in
  GitHub Pages.
