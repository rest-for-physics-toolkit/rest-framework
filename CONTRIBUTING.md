# Contributing to rest-framework

## Workflow

`master` is protected: it only changes through pull requests whose validation passes.

1. Create a branch from `master` (`git switch -c <user>-<topic> master`).
2. Commit the changes (one-line commit messages) and push the branch.
3. Open a pull request to `master` (as draft while the work is in progress), filling the template.
4. The validation (`.github/workflows/validation.yml`) runs on the pull request, and again at each push to
   it: the checks **Build**, **Unit tests**, **restG4 examples** and **Format** must pass, and the branch must
   be up to date with `master`, before merging. The **macOS** check (build and unit tests with the ROOT
   version of the Docker image from conda-forge, without Geant4, Garfield++ or restDAQ) also runs, but it
   is not required yet.
5. Merge the pull request; the validation runs again on `master` and the documentation is published.

## Before pushing

- Format: `pre-commit run --all-files` (clang-format with the line limit of 110, cmake-format, black and
  basic checks); `pre-commit install` runs it at each commit.
- Tests: build with `-DBUILD_TESTING=ON` and run `ctest` (all), `ctest -L <component>` (one component:
  tools, core, analysis, raw, detector, track, geant4, connectors, restG4, restGarfield, restDAQ) or
  `make check` (build and run all). New features need unit tests in the `test/` directory of their
  component (see below).
- Documentation: Doxygen in the `.cxx` of each class (keys, results, example), the README of the
  component, the user guide (doc/UserGuide.md) for the formats and keys, and MIGRATION_GUIDE.md for what
  comes from legacy. Each process is used in an example configuration (`examples/` of its component or of
  the repository; an optional process may be commented), listed in the table "Processes in the examples"
  of the README of its component.

## Tests

Each component has a `test/` directory with a GoogleTest executable (`<component>Tests`, sources in
`test/src/*.cxx`) added with `rest_add_tests(<component> <libraries>)` (cmake/RestTests.cmake): each test
has the label of its component; the packages restGarfield and restDAQ compile their sources into their test
executable (with Garfield++, readline and prometheus-cpp); restGarfield lists its tests from the sources
(`SCAN_SOURCES`, gtest_add_tests: Garfield++ prints at load time, which the discovery of the tests in the
executable would take as a test). The unit tests build their inputs (events,
metadata, datasets) in memory or in temporary files, without external data.

The pipeline tests (`examples/test`, label `pipelines`) run restManager with the YAML of the examples, as a
user would, on real data cut to a few hundred events, datasets of whole runs (for the dataset processes)
and fixed restG4 simulations, published in
[rest-common-data](https://github.com/rest-for-physics-toolkit/rest-common-data) (a fixed release, downloaded
by URL and cached), each test in its own work directory. The random processes of the examples have fixed
seeds (`seed`, a different one per process), so their outputs are reproducible. As the legacy pipeline
validation, the `AnalysisTree` of each output is compared event by event with a reference of
`examples/test/reference/` (every numeric column, scalars and vectors, relative tolerance 1e-5; the results
of fits, columns with `Gauss` in their name, depend on the platform: tolerance 1e-2, up to 2% of their
values may differ, and a fit without result, -1, on either side is not compared; the scores of classifiers
trained by TMVA, columns ending in `_score`, change on arm64, where ROOT uses fused multiply-adds: absolute
tolerance 0.1, up to 5% of their values may differ, and a selection by a score is compared on its number of
entries, within 30%). Large outputs (the datasets of whole runs) are compared on their first entries and
their number of entries, so the references stay small. After an intended change of the results, the
references are written again and committed with the change:

```
REST_UPDATE_REFERENCE=1 ctest -L pipelines
```

The results must be the same on every platform (they may change with the version of ROOT, not with the
compiler or the standard library), so the code follows these rules:

- random numbers from a `TRandom3` of the process (seeded with `TRestEvent::EventSeed`), with its own
  distributions (`Gaus`, `Poisson`, `Rndm`...) or algorithms written on them; never the distributions of
  `<random>` (`std::binomial_distribution`, `std::gamma_distribution`...), whose algorithms differ between
  libstdc++ and libc++;
- one random number per statement: the operands of an expression (`Gaus() * u + Gaus() * v`) and the
  arguments of a call are evaluated in an order that depends on the compiler (braced lists `{a, b, c}` are
  evaluated in order);
- `std::stable_sort` when the comparison does not order all the elements (by a coordinate, a time, an
  energy): `std::sort` leaves the equal ones in an order that depends on the library;
- no output depending on the iteration order of `std::unordered_map` or `std::unordered_set` (use
  `std::map`, `std::set` or a vector);
- the build disables fused multiply-adds (`-ffp-contract=off`), which clang uses by default on arm64.

Tests of a minute or more (e.g. Magboltz in restGarfield) go in a test suite whose name ends in `Long`
(`RestGarfieldLong`), so they can be skipped while developing with `ctest -L <component> -E Long`; the CI
runs them all.

## Versions and releases

The version of the code is the `VERSION` of `project()` in CMakeLists.txt; the build adds the commit and
whether the tracked files had uncommitted changes (cmake/RestVersion.cmake, at each build). A release is
the tag `vX.Y.Z` of the commit whose `VERSION` is `X.Y.Z` (a different tag on the commit is a warning of
the build). Builds without git give the commit with `-DREST_GIT_COMMIT=<commit>`.

A release, once `VERSION` is `X.Y.Z` in master (changed by a pull request, which also moves the entries of
"Unreleased" in CHANGELOG.md under `X.Y.Z`; a pull request with a notable change adds its entry there):

```bash
git switch master && git pull
git tag vX.Y.Z && git push origin vX.Y.Z      # documentation vX.Y.Z (docs.yml)
gh release create vX.Y.Z --generate-notes     # validation and image ghcr.io/<owner>/rest:vX.Y.Z
```

The validation of the release fails if its tag is not `vX.Y.Z` of `VERSION` (or `vX.Y.Z-<suffix>`, a
pre-release: `gh release create vX.Y.Z-rc1 --prerelease`, image without documentation), and the image is
pushed only if `restManager --version` in it gives the commit of the tag, the tag and a clean tree.

## Continuous integration

The jobs run in the Docker image of the `rest-docker` repository (ROOT, Geant4, Garfield++, CRY,
prometheus-cpp) pinned by each commit: `.github/rest-env` holds a tag of `rest-docker` (`vN`), and the jobs
use `ghcr.io/<owner>/rest-env:vN`, with `<owner>` the owner of this repository, and the `rest/Dockerfile` of
that tag for the `rest` image. A release is thus validated, and validated again, with its own environment,
and a new environment reaches master only through a pull request that changes `.github/rest-env`, after its
tag is pushed in `rest-docker` and its image is built. The repository variable `RESTV2_IMAGE` (Settings,
Secrets and variables, Actions, Variables) sets another image. The package of the image must be readable by
this repository (public, or Package settings, Manage Actions access).

| Workflow | When | What |
|---|---|---|
| `validation.yml` | Pull requests to master, pushes to master, releases, by hand | Build (all the packages, compiler warnings as errors, install, `restRoot`, Doxygen), unit tests (a step per component), restG4 examples, format (pre-commit). |
| `docs.yml` | Pushes to master, release tags `vX.Y.Z`, by hand | Doxygen documentation in GitHub Pages, a folder per version (`master/`, `vX.Y.Z/`) with a version selector; the site opens the newest release. |

The documentation site is the `gh-pages` branch, written by `.github/scripts/publish-docs.sh` (a single
commit, replaced at each publication); GitHub Pages publishes it (Settings, Pages, Source: Deploy from a
branch, `gh-pages`, `/ (root)`). The documentation of an old release is published again by running
`docs.yml` by hand on its tag.

Other files of `.github`: the pull request template, CODEOWNERS (reviewers requested automatically),
dependabot (weekly updates of the actions) and pr-badge.yml (badges of author, size and validation of each
pull request, with the [PR Badge](https://pullrequestbadge.com) GitHub App installed in the repository).

## Protection of master (administrators)

The protection is a ruleset of the repository, not a file: import `.github/rulesets/master.json` in
Settings, Rules, Rulesets, New ruleset, Import a ruleset. It applies to the default branch (master):

- no deletion and no force push;
- changes only through pull requests (no approvals required: a single developer cannot approve their own
  pull requests; raise `required_approving_review_count` when there are reviewers);
- the checks Environment, Build, Unit tests, restG4 examples and Format must pass, with the branch up to
  date (Environment too: the jobs after a failed one are skipped, and a skipped check counts as passed).

With the GitHub CLI: `gh api repos/<owner>/<repository>/rulesets --method POST --input .github/rulesets/master.json`.
The checks can only be selected once they have run at least once (open a first pull request).

## Build internals

- **Libraries and registration.** Events, metadata and processes register themselves (EventRegistry,
  MetadataClassRegistry) from static initializers, so the executables, macros and tests link the REST
  libraries with `--no-as-needed` (`rest_link_registered_libraries` in CMakeLists.txt): with `--as-needed`
  (the default of Ubuntu and others) a library whose symbols are not used would be dropped and its classes
  would not exist in `restManager`. The build-tree path of the libraries is written as RPATH (searched before
  `LD_LIBRARY_PATH`), so the programs of the build load its libraries and not those of another REST.
- **ROOT dictionaries.** Each library generates its dictionary with a LinkDef written at build time
  (`<library>_LinkDef.h`, cmake/scripts/RestGenerateLinkDef.cmake) selecting the class or namespace named as
  each header. Without it rootcling scans every declaration (also std and ROOT templates) and crashed with
  ROOT 6.36. The libraries find their headers and sources with `file(GLOB)`: a new file needs re-running
  `cmake`.
- **Macros.** Each macro (`source/macros`, `libraries/*/macros`, `packages/restGarfield/macros`) is compiled
  into `macros/lib<macro>_C.so`, without dictionary, and CMake writes `macros/REST_Macros.h` with the
  declarations of its `REST_*` functions and global pointers (`_run`, `_ev`, `_dataSet`), taken from the
  definitions starting at column 0. `restRoot` loads the libraries and declares that header, so the macros
  are available at start-up and run compiled. The `.C` files are installed as reference but are not in the
  macro path (`.x`/`.L` would redefine the declared functions); after changing a macro, rebuild and install.

## Publishing a package in its own repository

`packages/restG4`, `packages/restGarfield` and `packages/restDAQ` are developed in this repository. To publish
one of them alone (e.g. with its own tags), `git subtree` extracts the commits that touch it, with its files
at the root, and pushes them to another repository; the work is always done here (commits in the other
repository would need `git subtree pull` and conflict resolution).

```bash
# First and next publications (incremental)
git subtree push --prefix packages/restG4 git@github.com:<owner>/restG4.git master
# A tag in the other repository: on the branch generated by split
git subtree split --prefix packages/restG4 -b restg4-split
git push git@github.com:<owner>/restG4.git restg4-split:refs/tags/v3.0.0
git branch -D restg4-split
```

The history of the split starts with the first commit of rest-framework (the legacy history of restG4 is in
its original repository), so pushing to a repository with that legacy history is rejected unless forced.
