# Scripts

Scripts of rest-for-physics-toolkit, not built: `install.sh` (installer) and `restSubmit.py` (jobs of a
batch system, installed in `$REST_PATH/bin` as `restSubmit`).

## install.sh

Clones a release of rest-framework, builds and installs it with its dependencies, and writes
`<prefix>/setup.sh`, the environment to source in each shell (dependencies and REST). It can be run without
cloning the repository:

```bash
curl -fsSLO \
  https://raw.githubusercontent.com/rest-for-physics-toolkit/rest-framework/master/scripts/install.sh
bash install.sh --mode cvmfs --prefix /data/$USER/rest
source /data/$USER/rest/setup.sh
```

| Mode | Dependencies |
|---|---|
| `system` (default) | ROOT (required), Geant4 and Garfield++ (optional) already in the environment, e.g. a server; `--env <script>` sources a script loading them first (also in setup.sh). |
| `cvmfs` | The LCG view of `/cvmfs/sft.cern.ch` of `--lcg` (default `LCG_109/x86_64-el9-gcc13-opt`), e.g. lxplus or NAF. ROOT and Geant4 are those of the view; Garfield++ 4.0, the version of the CI, is built on top of it (about 10 minutes): the view has another one (LCG_109: a later commit, not tested with restGarfield), used instead with `--lcg-garfield`. |
| `source` | The versions of the CI (rest-env image), any Linux or macOS: the packages of the system checked with apt (Debian, Ubuntu), dnf (RHEL, AlmaLinux, Fedora) or brew (macOS), those the distribution has: the missing ones are installed with sudo only with `--install-packages`, otherwise the command to install them is printed (for the administrator) and the installation goes on (brew needs no privileges and installs them); ROOT (on Ubuntu 24.04 x86_64 the binary of root.cern used by the CI, elsewhere the same version built from its sources, about an hour), Geant4, Garfield++ and CRY. |

In every mode CRY (cry generator of restG4) and yaml-cpp are built if missing. Options: `--prefix` (default
`~/rest`), `--version` (a release `vX.Y.Z` or `master`; default the last release), `--jobs` (default the
CPUs, at most one per 2 GB of memory), `--daq` (restDAQ), `--no-g4`, `--no-garfield`; `--lcg-garfield`
(cvmfs) uses the Garfield++ of the view.

`<prefix>` holds `src/rest-framework` (the clone), `build/`, `install/` (`REST_PATH`), `deps/` (the
dependencies built), `setup.sh` and `install.log`. Run again with the same prefix it continues (the
dependencies built are kept) or moves to another `--version` (the clone must have no local changes). The
build does not treat warnings as errors (`REST_WARNINGS_AS_ERRORS=OFF`): the CI checks them on its compilers,
another compiler may warn of something new.

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
| `--setup` | Script loading the REST environment in the jobs (default `$REST_SETUP`, the `setup.sh` of `install.sh`, otherwise `$REST_PATH/thisREST.sh`). |
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
