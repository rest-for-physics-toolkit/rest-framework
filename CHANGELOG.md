# Changelog

Changes of each release of rest-for-physics-toolkit (rest-framework). The pull request that sets the
`VERSION` of a release moves the entries of "Unreleased" under its version (CONTRIBUTING.md, Versions and
releases); the notes of the GitHub release list its pull requests.

## Unreleased

- Continuous integration on the environment pinned by each commit: `.github/rest-env` holds the tag of the
  `rest-docker` repository (image `rest-env` and `rest/Dockerfile`), so a release is validated again with
  its own environment.

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
