#!/usr/bin/env python3
"""restSubmit: runs restG4 simulations or restManager processings as jobs of a batch system.

Two modes:

    restSubmit g4 <simulation.yaml> --jobs <n> [options] [-- <restG4 options>]
    restSubmit manager <config.yaml> --input <pattern> [options] [-- <restManager options>]

- g4: n jobs of restG4 with the same configuration and different seeds (--seed, the first one, job i
  uses seed + i; random by default), each one optionally processed by restManager (--process); --merge
  adds a job, run when all have finished, merging the simulations (and their processings) with
  REST_MergeRuns into a single simulation.
- manager: a job of restManager per input file (files of the pattern); --merge merges the outputs (the
  subruns of a run, see REST_MergeRuns).

Batch systems: HTCondor (a job cluster, and a DAG with the merge job: condor_submit_dag) and Slurm (a job
array, and the merge job with a dependency on it); --batch chooses one, by default the one installed
(condor_submit or sbatch). Each job loads the REST environment (--setup, by default
$REST_PATH/thisREST.sh) and runs in the directory where restSubmit was run, as restG4 or restManager would
run there. Everything goes to the submission directory (--dir, by default restSubmit/<name>): the
outputs (output/, the processings of the simulations in output/analysis/), the scripts and submit
files, and the stdout, stderr and logs of each job (logs/). --dry-run writes everything without submitting.

Examples:

    restSubmit g4 Muons.yaml --jobs 100 --time 8h --merge --process Analysis.yaml -- -n 1E5
    restSubmit manager FemDAQToTracks.yaml --input "data/R0122*.aqs" --memory 4000
"""

import argparse
import datetime
import glob
import os
import random
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

# Margin of the time limit of the batch system over the --time of restG4 (it stops at --time and writes)
TIME_MARGIN = 3600


def parse_time(expression):
    """Seconds of a time like 1h20m30s, 8h or 45m (the format of restG4 --time)."""
    seconds, number = 0, ""
    units = {"h": 3600, "m": 60, "s": 1}
    for character in expression:
        if character.isdigit():
            number += character
        elif character in units and number:
            seconds += int(number) * units[character]
            number = ""
        else:
            raise argparse.ArgumentTypeError(
                f"invalid time '{expression}' (e.g. 1h20m30s)"
            )
    if number or seconds <= 0:
        raise argparse.ArgumentTypeError(f"invalid time '{expression}' (e.g. 1h20m30s)")
    return seconds


def environment_variable(text):
    """An environment variable NAME=VALUE for the jobs."""
    name, equal, value = text.partition("=")
    if not equal or not name.isidentifier():
        raise argparse.ArgumentTypeError(f"invalid variable '{text}' (NAME=VALUE)")
    return name, value


def arguments(argv):
    """The options of restSubmit and those after --, passed to restG4 or restManager."""
    passed = []
    if "--" in argv:
        passed = argv[argv.index("--") + 1 :]
        argv = argv[: argv.index("--")]

    parser = argparse.ArgumentParser(
        prog="restSubmit",
        description="Runs restG4 or restManager as jobs of a batch system (HTCondor, Slurm).",
        epilog="Options after -- are passed to restG4 or restManager.",
    )
    modes = parser.add_subparsers(dest="mode", required=True)

    def common(mode):
        mode.add_argument("--batch", choices=["condor", "slurm"], help="batch system")
        mode.add_argument("--name", help="name of the submission and its files")
        mode.add_argument(
            "--dir", help="submission directory (default restSubmit/<name>)"
        )
        mode.add_argument(
            "--time", type=parse_time, help="time per job, e.g. 8h, 1h30m"
        )
        mode.add_argument("--memory", type=int, default=0, help="memory per job in MB")
        mode.add_argument("--threads", type=int, default=1, help="threads per job")
        mode.add_argument("--merge", action="store_true", help="merge the outputs")
        mode.add_argument("--setup", help="script loading the REST environment")
        mode.add_argument(
            "--env",
            type=environment_variable,
            action="append",
            default=[],
            help="environment variable NAME=VALUE of the jobs (repeatable)",
        )
        mode.add_argument(
            "--batch-option",
            action="append",
            default=[],
            help="line of the condor submit file or sbatch option (repeatable)",
        )
        mode.add_argument("--dry-run", action="store_true", help="write, do not submit")

    g4 = modes.add_parser("g4", help="restG4 jobs with different seeds")
    g4.add_argument("config", help="restG4 configuration (YAML)")
    g4.add_argument("--jobs", type=int, required=True, help="number of jobs")
    g4.add_argument(
        "--seed", type=int, help="seed of the first job (random by default)"
    )
    g4.add_argument(
        "--process", help="restManager configuration run on each simulation"
    )
    common(g4)

    manager = modes.add_parser("manager", help="a restManager job per input file")
    manager.add_argument("config", help="restManager configuration (YAML)")
    manager.add_argument("--input", required=True, help="input files (pattern)")
    common(manager)

    args = parser.parse_args(argv)
    args.passed = passed
    if args.threads < 1 or args.memory < 0:
        parser.error("--threads must be at least 1 and --memory positive")
    if args.mode == "g4":
        if args.jobs < 1:
            parser.error("--jobs must be at least 1")
        if args.merge and args.jobs < 2:
            parser.error("--merge needs at least 2 jobs")
        if args.seed is not None and args.seed < 1:
            parser.error("--seed must be positive")
    if not Path(args.config).is_file():
        parser.error(f"configuration {args.config} not found")
    if getattr(args, "process", None) and not Path(args.process).is_file():
        parser.error(f"configuration {args.process} not found")
    return args


def batch_system(requested):
    """The batch system requested, or the one installed."""
    if requested:
        return requested
    if shutil.which("condor_submit"):
        return "condor"
    if shutil.which("sbatch"):
        return "slurm"
    sys.exit(
        "restSubmit: neither condor_submit nor sbatch found, choose one with --batch"
    )


def write_script(path, lines):
    path.write_text("\n".join(lines) + "\n")
    path.chmod(0o755)


def job_header(args, setup):
    """Start of the script of a job: the REST environment and the working directory."""
    lines = ["#!/bin/bash", "set -euo pipefail", f"source {shlex.quote(str(setup))}"]
    lines += [f"export {name}={shlex.quote(value)}" for name, value in args.env]
    lines.append(f"cd {shlex.quote(os.getcwd())}")
    return lines


def g4_jobs(args, setup, directory):
    """Script of the restG4 jobs (argument: the job index) and the merge script."""
    output = directory / "output"
    seed = args.seed if args.seed is not None else random.randint(1, 2**30)
    config = Path(args.config).resolve()
    q = shlex.quote
    simulation = f"{q(str(output))}/{args.name}_${{i}}.root"
    command = f'restG4 {q(str(config))} --output "$tmp" --seed $(({seed} + i))'
    if args.threads > 1:
        command += f" --threads {args.threads}"
    if args.time:
        command += f" --time {args.time}s"
    command += "".join(" " + q(option) for option in args.passed)

    lines = job_header(args, setup) + [
        "i=$1",
        # Written in tmp/ and moved at the end: output/ has only finished files
        f'tmp={q(str(directory / "tmp"))}/{args.name}_${{i}}.root',
        command,
    ]
    if args.process:
        process = q(str(Path(args.process).resolve()))
        analysis = f'{q(str(directory / "tmp"))}/{args.name}_${{i}}_analysis.root'
        lines += [
            f'restManager --c {process} --i "$tmp" --o {analysis}',
            f'mv {analysis} {q(str(output / "analysis"))}/',
        ]
    lines.append(f'mv "$tmp" {simulation}')
    write_script(directory / "job.sh", lines)

    merge = []
    if args.merge:
        merge = [(f"{output}/{args.name}_*.root", f"{directory}/{args.name}.root")]
        if args.process:
            merge.append(
                (
                    f"{output}/analysis/{args.name}_*_analysis.root",
                    f"{directory}/{args.name}_analysis.root",
                )
            )
    print(f"restG4 jobs: seeds {seed} to {seed + args.jobs - 1}")
    return args.jobs, merge


def manager_jobs(args, setup, directory):
    """Script of the restManager jobs (argument: the job index, the line of inputs.txt)."""
    inputs = sorted(str(Path(name).resolve()) for name in glob.glob(args.input))
    if not inputs:
        sys.exit(f"restSubmit: no file matches {args.input}")
    (directory / "inputs.txt").write_text("\n".join(inputs) + "\n")

    output = directory / "output"
    q = shlex.quote
    command = (
        f'restManager --c {q(str(Path(args.config).resolve()))} --i "$input" --o "$tmp"'
    )
    if args.threads > 1:
        command += f" --j {args.threads}"
    command += "".join(" " + q(option) for option in args.passed)
    lines = job_header(args, setup) + [
        "i=$1",
        f'input=$(sed -n "$((i + 1))p" {q(str(directory / "inputs.txt"))})',
        'name=$(basename "${input%.*}")',
        f'tmp={q(str(directory / "tmp"))}/"$name"_{args.name}.root',
        command,
        f'mv "$tmp" {q(str(output))}/',
    ]
    write_script(directory / "job.sh", lines)

    merge = []
    if args.merge:
        merge = [(f"{output}/*_{args.name}.root", f"{directory}/{args.name}.root")]
    print(f"restManager jobs: {len(inputs)} input files")
    return len(inputs), merge


def merge_script(args, setup, directory, merges):
    """Script of the merge job: REST_MergeRuns of each pattern, failing if a merge fails."""
    lines = job_header(args, setup)
    for pattern, merged in merges:
        call = f'gSystem->Exit(REST_MergeRuns("{pattern}", "{merged}"))'
        lines.append(f"restRoot -b -q -e {shlex.quote(call)}")
    write_script(directory / "merge.sh", lines)


def submit_condor(args, directory, jobs, merging):
    logs = directory / "logs"

    def submit_file(name, executable, arguments, queue):
        lines = [
            f"executable = {executable}",
            f"arguments = {arguments}",
            f"output = {logs}/{name}.$(Process).out",
            f"error = {logs}/{name}.$(Process).err",
            f"log = {logs}/{name}.log",
            f"request_cpus = {args.threads}",
        ]
        if args.memory:
            lines.append(f"request_memory = {args.memory}")
        if args.time:
            # Attributes of the runtime of the sites: DESY NAF (RequestRuntime), CERN (MaxRuntime)
            runtime = args.time + TIME_MARGIN
            lines += [f"+RequestRuntime = {runtime}", f"+MaxRuntime = {runtime}"]
        lines += args.batch_option + [f"queue {queue}"]
        path = directory / f"{name}.sub"
        path.write_text("\n".join(lines) + "\n")
        return path

    jobs_file = submit_file("job", directory / "job.sh", "$(Process)", jobs)
    if not merging:
        command = ["condor_submit", str(jobs_file)]
    else:
        merge_file = submit_file("merge", directory / "merge.sh", "", 1)
        dag = directory / "jobs.dag"
        dag.write_text(
            f"JOB jobs {jobs_file}\nJOB merge {merge_file}\nPARENT jobs CHILD merge\n"
        )
        command = ["condor_submit_dag", "-outfile_dir", str(logs), str(dag)]
    run(args, command)


def submit_slurm(args, directory, jobs, merging):
    logs = directory / "logs"

    def options(name):
        result = [
            f"--job-name={args.name}_{name}",
            f"--cpus-per-task={args.threads}",
        ]
        if args.memory:
            result.append(f"--mem={args.memory}M")
        if args.time:
            seconds = args.time + TIME_MARGIN
            result.append(f"--time={seconds // 3600}:{seconds % 3600 // 60:02d}:00")
        return result + args.batch_option

    command = ["sbatch", "--parsable", f"--array=0-{jobs - 1}"]
    command += [f"--output={logs}/job.%a.out", f"--error={logs}/job.%a.err"]
    command += options("job") + ["--wrap", f"{directory}/job.sh $SLURM_ARRAY_TASK_ID"]
    array = run(args, command) or "<job array>"
    if merging:
        command = ["sbatch", "--parsable", f"--dependency=afterok:{array}"]
        command += [f"--output={logs}/merge.out", f"--error={logs}/merge.err"]
        command += options("merge") + [str(directory / "merge.sh")]
        run(args, command)


def run(args, command):
    """Runs a submission command (printed; only printed with --dry-run), returning its output."""
    print(" ".join(shlex.quote(part) for part in command))
    if args.dry_run:
        return ""
    result = subprocess.run(
        command, check=True, stdout=subprocess.PIPE, universal_newlines=True
    )
    print(result.stdout, end="")
    return result.stdout.strip().split(";")[0]


def main(argv):
    args = arguments(argv)
    batch = batch_system(args.batch)

    setup = args.setup
    if not setup:
        if "REST_PATH" not in os.environ:
            sys.exit(
                "restSubmit: REST_PATH not set (source thisREST.sh) and no --setup"
            )
        setup = Path(os.environ["REST_PATH"]) / "thisREST.sh"
    setup = Path(setup).resolve()
    if not setup.is_file():
        sys.exit(f"restSubmit: setup script {setup} not found")

    if not args.name:
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        args.name = f"{Path(args.config).stem}_{stamp}"
    directory = Path(args.dir or Path("restSubmit") / args.name).resolve()
    if directory.exists() and any(directory.iterdir()):
        sys.exit(f"restSubmit: {directory} is not empty")
    for sub in ("output", "tmp", "logs") + (
        ("output/analysis",) if args.mode == "g4" and args.process else ()
    ):
        (directory / sub).mkdir(parents=True, exist_ok=True)

    if args.mode == "g4":
        jobs, merges = g4_jobs(args, setup, directory)
    else:
        jobs, merges = manager_jobs(args, setup, directory)
    if merges:
        merge_script(args, setup, directory, merges)

    if batch == "condor":
        submit_condor(args, directory, jobs, bool(merges))
    else:
        submit_slurm(args, directory, jobs, bool(merges))
    print(f"Submission directory: {directory}")
    for _, merged in merges:
        print(f"Merged output: {merged}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
