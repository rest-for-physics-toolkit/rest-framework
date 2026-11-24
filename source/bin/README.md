# REST executables (source/bin)

Each `.cxx` file in this directory is built as an executable of the same name, linked against RestCore and
every REST library built, and installed in `$REST_PATH/bin`.

## restManager

Runs a YAML configuration with a `TRestManager` section: event processing, dataset generation, dataset
processing and analysis plots.

```
restManager --c|--config <config.yaml> [--i|--input <input.root>] [--o|--output <output.root>]
            [--v|--verbose <level>] [--e|--events <n>] [--j|--threads <n>]
            [--r|--reference <process>=<dataset>] [--set <name>=<value>]
            [--maxRepeated <n>] [--repeatedInterval <s>]
restManager --version
```

| Option | Description |
|---|---|
| `--c`, `--config` | Configuration file (required). |
| `--i`, `--input` | Input: `inputFileName` of the TRestRun; with a TRestDataSet its `filePattern`; with a TRestDataSetProcessManager and no TRestDataSet its `inputFileName`; with no run nor dataset, the `inputFileName` of every TRestAnalysisPlot. |
| `--o`, `--output` | Output: `outputFileName` of the TRestRun, and of the TRestDataSetProcessManager if any, otherwise of the TRestDataSet. |
| `--v`, `--verbose` | Global verbose level. |
| `--e`, `--events` | `eventsToProcess` of the TRestEventProcessManager. |
| `--j`, `--threads` | `nThreads` of the TRestEventProcessManager. |
| `--r`, `--reference` | `referenceFileName` of a process of the TRestDataSetProcessManager (repeatable), e.g. `--r calib=D_00120_00122_Calibration.root`. |
| `--set` | Sets a parameter of the `parameters` section for this execution (repeatable; the file is not changed), e.g. `--set meshVoltage=300V`. A parameter not used is an error. |
| `--maxRepeated` | Similar messages (warning and info) from the same code line within the time window after which the next ones are suppressed (default 10, 0: no suppression; never with verbose debug). |
| `--repeatedInterval` | Time window of `--maxRepeated` in seconds (default 10, 0: a single window for the whole execution). |
| `--version` | Prints the version of the code (`VERSION` of the project, commit of the build and whether it had uncommitted changes) and exits. It is also printed at the start and stored in the output files (`codeVersion`, `codeCommit`, `codeClean`). |

SIGINT, SIGTERM, SIGHUP and SIGQUIT stop the processing cleanly. The exit code is 0 on success and 1 on error.

## restRoot

ROOT prompt with REST loaded; the version of the code is printed at its start (`REST_Version::GetSummary()`).

```
restRoot [ROOT options] [file.root] [macro.C ...]
```

- Requires `REST_PATH`. Loads `$REST_PATH/lib/libRest*.so` (except `libRestG4.so`), adds
  `$REST_PATH/include` to the include path and loads the compiled macros of `$REST_PATH/macros`
  (see [source/macros](../macros/README.md)).
- A REST run file or dataset among the arguments is opened with `REST_OpenInputFile` (globals `_run`,
  `_ev` or `_dataSet`); if several are given only the last one is opened (with a warning). The other
  arguments are passed to ROOT (TRint).

## restSubmit

Runs restG4 simulations or restManager processings as jobs of a batch system, HTCondor or Slurm (a Python 3
script, `restSubmit.py`, installed as `$REST_PATH/bin/restSubmit`).

```
restSubmit g4 <simulation.yaml> --jobs <n> [--seed <s>] [--process <analysis.yaml>] [options]
           [-- <restG4 options>]
restSubmit manager <config.yaml> --input <pattern> [options] [-- <restManager options>]
```

- `g4`: n restG4 jobs of the same configuration, job i with the seed s + i (s random by default), each
  one optionally processed by restManager (`--process`).
- `manager`: a restManager job per file of the pattern (output `<input name>_<name>.root`).
- `--merge`: a last job, run when all have finished, merges the outputs with `REST_MergeRuns` (in `g4`
  into a single simulation, also the processings; in `manager` the subruns of a run).

| Option | Description |
|---|---|
| `--batch` | `condor` (a job cluster; with `--merge` a DAG, `condor_submit_dag`) or `slurm` (a job array; the merge job depends on it). By default the one installed (`condor_submit` or `sbatch`). |
| `--name`, `--dir` | Name of the submission (default: the configuration and the date) and its directory (default `restSubmit/<name>`, must be empty). |
| `--time` | Time per job, e.g. `8h`, `1h30m`: restG4 `--time` (it stops and writes) and the limit of the batch system, one hour more (condor `+RequestRuntime` of DESY NAF and `+MaxRuntime` of CERN, Slurm `--time`). |
| `--memory`, `--threads` | Memory (MB) and threads of each job (restG4 `--threads`, restManager `--j`). |
| `--setup` | Script loading the REST environment in the jobs (default `$REST_PATH/thisREST.sh`). |
| `--env NAME=VALUE` | Environment variable of the jobs (repeatable). |
| `--batch-option` | A line added to the condor submit files or an option of sbatch (repeatable), e.g. a queue or an account of the site. |
| `--dry-run` | Writes the scripts and submit files and prints the submission commands without running them. |

The jobs run in the directory where restSubmit was run (as restG4 or restManager there) and write in the
submission directory: `output/` (the finished outputs; the processings of `g4` in `output/analysis/`),
`logs/` (stdout, stderr and logs of the batch system), `job.sh` and `merge.sh` (the scripts of the jobs, the
argument of job.sh is the job index) and the submit files. The merged outputs are `<dir>/<name>.root` (and
`<name>_analysis.root`).

```bash
restSubmit g4 Muons.yaml --jobs 100 --time 8h --merge --process Analysis.yaml -- -n 1E5
restSubmit manager FemDAQToTracks.yaml --input "data/R0122*.aqs" --memory 4000
```

Other executables are built by the packages (`restG4`, `restGarfield`, `femclient` of restDAQ).
