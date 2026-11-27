# The restGarfield package

This package is used to compute with Garfield++ the transport and gain properties of the detector media,
stored in tables that the detector library reads without depending on Garfield++.

Optional application that computes with [Garfield++](https://garfieldpp.web.cern.ch) the transport and gain
tables of a medium defined in a `TRestDetectorTransportManager` YAML configuration, and stores them in the
medium file (`mediumFile` of the `TRestDetectorMedium`, by default `<name>.root`). The REST libraries do not
depend on Garfield++: the tables are read by `TRestDetectorMedium` (detector library).

- **Transport tables** (`TRestDetectorTransportTable`): drift velocity, longitudinal and transversal
  diffusion, Townsend and attachment coefficients versus the electric field, for electrons, holes and ions
  (the carriers available in the Garfield medium). Gases use Magboltz (`MediumMagboltz`, up to 6
  components); solids of a single component `Si`, `GaAs`, `GaN`, `CdTe` or `diamond` use the Garfield solid
  media.
- **Gain tables** (`TRestDetectorGainTable`, gases only): mean gain, Polya theta and the avalanche size
  distributions, from microscopic avalanches of one electron in a uniform field across a gap (parallel
  plates).
- **Drift responses in field maps** (`TRestDetectorDriftResponse`, gases only, special simulations with a
  non-uniform drift field): electrons drifted with `AvalancheMC` from the nodes of a grid, in a field map
  exported by CST Studio Suite, to a plane above the readout (the start of the amplification), with the
  gas of a transport table of the medium file; per node the arrived fraction and the histograms of the
  displacement on the readout plane and of the drift time. Each response goes to its own `responseFile`,
  read by `TRestDetectorHitsDriftProcess` (detector library, `responseFile`) instead of the uniform drift
  field of the medium.
- **Amplification responses in field maps** (`TRestDetectorAmplificationResponse`, gases only, special
  simulations): microscopic avalanches (`AvalancheMicroscopic`, Magboltz with the composition of the medium)
  in the field maps of a unit cell of the amplification structure (e.g. a microbulk Micromegas: a fine map
  of the mesh, the holes and the kapton, a coarse map of the drift field above it, the cell repeated along
  x and y), started at the plane where the drift response ends; the transparency (fraction of the electrons
  amplified), the histogram of the avalanche size and of the time to the anode. Each response goes to its
  own `responseFile`, read by `TRestDetectorAvalancheProcess` (`responseFile`) instead of the gain table.

Tables already in the medium file, and response files that exist, are skipped unless `overwrite: true`.

## Build

CMake option `REST_GARFIELD` (default `ON`). It needs Garfield++ (found with `GARFIELD_INSTALL` or
`GARFIELD_HOME`) and the detector library; otherwise the package is skipped with a message.

```
cmake -DREST_GARFIELD=ON ..
```

## Usage

```
restGarfield <config.yaml> [--set name=value ...] [--jobs N]
restGarfield --version
```

`--version` prints the version of the code, also printed at the start; each table keeps in its `generator`
(entry of the table in the metadata of the medium file) the configuration used, the `program`, the `date`
and the `code` (version, commit and clean state of restGarfield).

`--set` overrides a parameter of the `parameters` section for this execution (repeatable). `--jobs N`
splits the field nodes of each gas transport table in up to N processes, each computing a contiguous part
of the grid of at least 2 nodes (Magboltz is not thread safe, so processes instead of threads); the table is
the same. Several media can also run at once as separate `restGarfield` processes, even writing to the same
medium file (e.g. `Ar99_iC4H10_1`, `Ar98_iC4H10_2`... in `ArIso.root`, each one in `tables/<medium>/`): the
file is locked while it is read or written. Example based on
`examples/ArIso.yaml` (installed in `examples/restGarfield`; its `gain` and `drift`
sections are commented out there):

```yaml
transportManager:
  class: "TRestDetectorTransportManager"
  medium:
    class: "TRestDetectorMedium"
    name: "ArIso"
    state: "gas"
    components: {Ar: 99, iC4H10: 1}
  transport:
    temperatures: [293.15]
    fieldRange: [30 V/cm, 500 V/cm]
    fieldNodes: 30
  gain:
    gap: 50 um
    reducedFieldRange: [50 kV/cm/bar, 70 kV/cm/bar]
```

## Configuration keys

`TRestDetectorTransportManager`: `medium` (a `TRestDetectorMedium`), `transport`, `gain`, `drift` and
`amplification` (each a map or a list of maps, one per set of tables or response). The drift responses are
computed after the tables, so the transport table they use can be computed in the same run; the
amplification responses last.

`transport`:

| Key | Default | Description |
|---|---|---|
| `pressure` | 1 bar | Pressure of the tables (gas tables are scaled to other pressures with E/p). |
| `temperatures` | medium temperature | One table per temperature (K). |
| `fieldRange` | [10 V/cm, 1000 V/cm] | Electric field range. |
| `fieldNodes` | 20 | Number of field nodes. |
| `logField` | true | Logarithmic field spacing. |
| `nCollisions` | 10 | Magboltz collisions (x 10^7). |
| `maxElectronEnergy` | 400 eV | Maximum electron energy in Magboltz. |
| `usePenning` | false | Penning transfer in the Townsend coefficient (gases). |
| `penningRate` | 0 | Penning transfer probability, (0, 1]. |
| `overwrite` | false | Compute again the tables already in the file. |

`gain`:

| Key | Default | Description |
|---|---|---|
| `gap` | (required) | Amplification gap. |
| `pressures` | medium pressure | One set of tables per pressure. |
| `temperatures` | medium temperature | One table per temperature (K). |
| `reducedFieldRange` | (required) | Reduced field range E/p; the fields are E = (E/p) * p. |
| `fieldNodes` | 5 | Number of equally spaced field nodes. |
| `avalanches` | 1000 | Avalanches (one primary electron each) per field node. |
| `usePenning` | false | Penning transfer. |
| `penningRate` | 0 | Penning transfer probability. |
| `overwrite` | false | Compute again the tables already in the file. |

`drift` (lengths in mm unless a unit is given):

| Key | Default | Description |
|---|---|---|
| `fieldMap` | (required) | ASCII export of the electric field by CST (header `x [mm] y [mm] z [mm] Ex [V/m] Ey [V/m] Ez [V/m]`, the units read from it; x fastest). Nodes with the field exactly 0 are out of the gas (CST writes 0 out of the solved volume, e.g. the corners of the box around a cylinder, and in the conductors): an electron reaching a cell with one of them is lost, as at a wall; start points there are not drifted (survival 0). |
| `fieldMapOffset` | [0, 0, 0] | Position of the origin of the CST frame in the detector frame. |
| `readouts` | (required) | YAML file with the `TRestDetectorReadoutManager` of the readout (distance and normal; the normal must be along x, y or z). |
| `readout` | data readout | Name of the readout in the manager. |
| `boundary` | 0 | Distance from the readout to the plane where the drift ends (above the mesh); the field map must reach above it. |
| `transportTable` | the only one | Key of the transport table of the medium file (the gas at the pressure of the medium). |
| `xRange`, `yRange`, `zRange` | the readout | Grid of start points (detector frame). A range not given along an axis of the readout plane is the extent of the readout channels plus `margin`; the range along the readout normal is required. The response file keeps the ranges used. |
| `margin` | 5 mm | Added on each side of the extent of the readout channels. |
| `gridStep` | (required) | Largest spacing of the nodes along x, y and z: each range is divided into the fewest equal steps not larger than it (at least 2 nodes per axis), e.g. a multiple of the readout pitch. |
| `electrons` | 1000 | Electrons drifted from each node (at most 65535: the counts are stored as 16-bit integers, the file compressed with LZMA). |
| `stepSize` | 1 mm | Distance step of `AvalancheMC`. |
| `positionBins` | 20 | Bins per axis of the displacement histogram of each node. |
| `timeBins` | 20 | Bins of the drift time histogram of each node. |
| `responseFile` | (required) | Output file of the response. |
| `overwrite` | false | Compute the response again if the file exists. |

`--jobs N` also splits the start points of the drift response in N processes (every N-th point each,
so that each process gets points at all the drift lengths; each process with its own random seed).

`amplification` (lengths in mm unless a unit is given; z along the drift, the electrons moving to -z):

| Key | Default | Description |
|---|---|---|
| `fieldMaps` | (required) | ASCII exports of the electric field by CST, as `fieldMap` of `drift` (nodes with the field 0, the conductors, stop the electrons). Where they overlap, the first one wins: the fine map of the amplification structure first, then the coarse map of the drift field above it (different steps, not one regular grid). |
| `dielectricMaps` | none | Exports of the displacement field D (header `... Dx [C/m^2] ...`) on nodes of each field map, in the same order (`""` for a map without): a node with D > 1.5 eps0 E (a dielectric, e.g. kapton; 1 in the gas) stops the electrons. It may cover part of its field map (e.g. only the layers with kapton). |
| `fieldMapOffset` | [0, 0, 0] | Position of the origin of the CST frame in the frame of the response. |
| `periodicity` | "xy" | Axes along which the cell repeats; the period is the extent of the maps. |
| `mirror` | false | Mirror periodicity (the cell reflected) instead of a translation. |
| `start` | (required) | z of the plane where the primary electrons start, at random over the cell: the `boundary` of the drift response, measured from the mesh. |
| `anode` | (required) | z below which an electron of the avalanche is collected by the anode (a little above the anode surface: the electrons stop in the last cell before the conductor). |
| `electrons` | 1000 | Avalanches, one primary electron each. |
| `usePenning`, `penningRate` | false, 0 | Penning transfer of the gas. |
| `sizeBins` | 100 | Bins of the avalanche size histogram (over [0, the largest size]). |
| `timeBins` | 50 | Bins of the histogram of the time from the start plane to the anode (mean time of the electrons of an avalanche). |
| `responseFile` | (required) | Output file of the response. |
| `overwrite` | false | Compute the response again if the file exists. |

A primary electron whose avalanche brings no electron to the anode (stopped on the mesh, on a dielectric
or attached) is lost: the transparency is the fraction of primaries with an avalanche at the anode. Where
the lost ones ended is counted and printed (out of the cell, on a conductor, on a dielectric, attached,
other: `losses` in the metadata of the response file). One
response per operating point (the voltages of the CST model). `--jobs N` splits the avalanches in N
processes. E.g. the pixel cell of a microbulk (mesh at z = 0, anode below -55 um):

```yaml
  amplification:
    fieldMaps: ["pixel_mesh320_E_fine.txt", "pixel_mesh320_E_coarse.txt"]
    dielectricMaps: ["pixel_mesh320_D_fine.txt"]
    periodicity: "xy"
    start: 1 mm
    anode: -54 um
    electrons: 2000
    usePenning: true
    penningRate: 0.3
    responseFile: "amplification_320V.root"
```

See the TRestDetectorTransportManager, TRestDetectorMedium, TRestDetectorDriftResponse and
TRestDetectorAmplificationResponse documentation for the validation of these keys.

## Macros

- `REST_Garfield_DrawMedium(fileName, medium = "", pressure = 0)`: draws the transport tables (versus E,
  one color per table, one line style per carrier; scaled to `pressure` in bar if > 0) and the gain tables
  (mean gain and Polya theta) of a medium file.

- `REST_Garfield_DrawDriftResponse(fileName, z = -1)`: draws a drift response, for the layer of start
  points at `z` (by default the farthest from the readout): maps on the readout plane of the survival, the
  mean displacement of the arrival (the distortion of the field, with arrows), the mean drift time and the
  transverse spread; the radial displacement versus the radius; the drift time versus z at the centre and
  the edges. The tree `DriftResponse` of the file has the same summary per start point (`x`, `y`, `z`,
  `survival`, `meanU`, `meanV`, `meanT`, `sigmaU`, `sigmaV`, `sigmaT`) for other plots.

```
restRoot
root [0] REST_Garfield_DrawMedium("ArIso.root")
root [1] REST_Garfield_DrawMedium("ArIso.root", "", 1.4)
root [2] REST_Garfield_DrawDriftResponse("AlphaCAMM_drift.root")
root [3] DriftResponse->Draw("meanU:x", "z == 90 && abs(y) < 5")   // after TFile f("AlphaCAMM_drift.root")
```

## Tests

With `BUILD_TESTING` (`test/`, ctest label `restGarfield`): the Garfield++ gas of a medium and the transport
and gain tables computed with small grids (a silicon table, a gas table with one Magboltz collision in two
processes, a gain table of a few avalanches at a low field), written to medium files of the build directory
and read back; a small CST field map (with a node out of the gas), the Garfield++ gas of a transport table
and the drift response in a uniform field map (straight drift to the boundary plane at the drift velocity of
the table); two field maps combined, with a dielectric map and periodicity; the amplification response of
a few avalanches in a periodic cell at a low field, in two processes. The suites `*Long` (the transport
table in two processes) are skipped with `ctest -L restGarfield -E Long`.
