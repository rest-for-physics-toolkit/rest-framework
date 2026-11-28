# Changelog

Changes of each release of rest-for-physics-toolkit (rest-framework). The pull request that sets the
`VERSION` of a release moves the entries of "Unreleased" under its version (CONTRIBUTING.md, Versions and
releases); the notes of the GitHub release list its pull requests.

## Unreleased

- `examples/DataSetLogOdds.yaml`: ranges of `MaxTrackSigmaZ2` ([0, 0.3] mm2, 60 bins) and
  `trackAnalysis_MaxTrackEnergyRatio` ([0, 1], 20 bins); with the legacy ranges ([0, 20] and [0, 0.1])
  their probability densities did not discriminate (all in the first bin, or out of range). References of
  the IAXO dataset chain test regenerated.
- restG4: the units of the content of the histograms of the sources from the title of its axis (z of a
  TH2, y of a TH1): counts / s / cm2 per bin or per eV, keV, MeV or GeV (each bin weighted by its width, also
  to sample it); the flux and the exposure of the `cosmic` generator are not computed without them.
  `activity` of the generator (e.g. `10 mBq`, units `uBq` to `MBq`): the exposure of the other generators
  (events / activity, e.g. a contamination). `flux` of the `surface` generator with a shape (counts / s /
  cm2; for a sphere the omnidirectional flux, times its cross section): the exposure of e.g. the
  environmental gammas from a sphere, also with volume biasing.
- `exposureTime` of the runs (`TRestRunInfo`): the real time a simulation represents (not its computing
  time), written by restG4 (cosmic and cry generators), kept by the runs processed from it and summed by
  `REST_MergeRuns`; datasets (`GetTotalExposureTime`) and `TRestAnalysisPlot` (`exposureTime`, hours, as
  `runLength`) to give the rates of the simulations. The computing time of restG4 is stored as
  `computingTime` of the `TRestGeant4Metadata` (`GetComputingTime`; it was not stored).
- restG4, formula `CosmicMuons` (Guan): the azimuth factor 2 pi was counted twice (as legacy before
  geant4lib#151): its flux is 6.28 times lower (0.0172 counts / s / cm2); the sampled distribution is the
  same.
- restG4, `cosmic` generator: the flux of its source (the integral of its distribution of energy and angle
  in `range`) is stored (`flux`), and `TRestGeant4Metadata::GetEquivalentExposureTime` gives the time a
  simulation represents (events / (flux * surface); the exposure of CRY for the `cry` generator).
- restG4, `cry` generator: as the `cosmic` one, only the particles crossing the world are simulated, from
  where they enter it (their time delayed by the flight); the showers without any are skipped (they count in
  the exposure of CRY, not in the events). The CRY square (`subBoxLength`) can be larger than the world.
  Each thread has its own CRY (seed of the simulation plus the thread index), not one shared under a lock:
  the simulation scales with the threads; the exposure is the sum of all.
- TRestAnalysisPlot: `allEntries` in the labels of the panels (the entries with the global cuts only),
  e.g. `{entries/allEntries}` for the fraction of the events that pass the cuts of the panel.
- TRestAnalysisPlot: the integer values of the labels (e.g. entries) are divided as real numbers
  (`{entries/allEntries}` was an integer division).
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
- TRestDataSetSignalSelection: `applyCuts` (false by default): the cuts select the events of the reference
  only, unless `applyCuts` is true or the reference is the processed dataset itself (e.g. a fiducial cut
  of the calibration does not restrict the background).
- TRestRawToDetectorSignalProcess: `keepEmptyEvents` keeps the events without signals (e.g. vetoes without
  any pulse), so the output has an entry per input event.
- restGarfield: amplification responses in the CST field maps of a unit cell of the amplification
  structure (`amplification` section, `TRestDetectorAmplificationResponse`): microscopic avalanches from the
  plane where the drift response ends; transparency, avalanche size and time histograms, and where the lost
  electrons end (conductor, dielectric, out of the cell). Several field maps combined (a fine map of the
  mesh and a coarse one above it), dielectrics from an export of the displacement field D, periodic cells.
- `TRestDetectorHitsDriftResponseProcess` removed: the drift response is the `responseFile` of
  `TRestDetectorHitsDriftProcess` (instead of the drift field of the medium); the amplification response the
  `responseFile` of `TRestDetectorAvalancheProcess` (instead of `amplificationVoltage` and `gap`).
- Amplification response depending on the position in the cell: restGarfield maps the transparency and the
  mean avalanche size by the start point (`cellBins`, TH2D maps in the response file) and counts the
  electrons at the anode on the pad of another cell; `TRestDetectorAvalancheProcess` applies the bin of each
  hit in its pixel of the readout (`readout`; `TRestDetectorReadout::GetPixelPosition`: the position in the
  frame of the volume containing it, whatever the rotations). Responses computed before have no map (the
  cell average). restGarfield prints the rms of the gain factors of the bins and its statistical value
  (~500 primaries per bin for ~3 %).

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
