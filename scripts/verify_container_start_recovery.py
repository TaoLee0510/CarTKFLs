#!/usr/bin/env python3
"""Check targeted container-start recovery before rerunning the fit audit."""

import argparse
import csv
import hashlib
import re
import subprocess
from pathlib import Path


RAW_FILES = (
    "bootstrap_res.Rds",
    "landscape.Rds",
    "landscape_posterior_samples.Rds",
    "xval.Rds",
)


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def command(*args):
    return subprocess.run(
        args, check=True, universal_newlines=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    ).stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve(strict=True)
    evidence = rows(run_dir / "manifests" / "container_start_recovery_evidence.tsv")
    if not evidence or len({row["task_id"] for row in evidence}) != len(evidence):
        raise RuntimeError("Recovery evidence is empty or has duplicate task IDs")
    submissions = [row for row in rows(run_dir / "submissions.tsv")
                   if row["stage"].startswith("container_start_recovery_MINOBS_")]
    if sum(int(row["tasks"]) for row in submissions) != len(evidence):
        raise RuntimeError("Recovery submissions do not match the evidence manifest")
    job_ids = [row["job_id"] for row in submissions]
    if command("squeue", "-r", "-h", "-j", ",".join(job_ids), "-o", "%A").strip():
        raise RuntimeError("Recovery tasks are still active")
    observed = {}
    for line in command(
        "sacct", "-X", "-j", ",".join(job_ids), "--noheader",
        "--parsable2", "-o", "JobID,State,ExitCode",
    ).splitlines():
        job_id, state, exit_code, *_ = line.split("|")
        match = re.fullmatch(r"(\d+)_(\d+)", job_id)
        if match and match.group(1) in job_ids:
            observed[job_id] = (state, exit_code)
    if len(observed) != len(evidence) or any(
        value != ("COMPLETED", "0:0") for value in observed.values()
    ):
        raise RuntimeError(f"Recovery Slurm states are incomplete or failed: {observed}")

    for row in evidence:
        stderr = Path(row["stderr_path"])
        if hashlib.sha256(stderr.read_bytes()).hexdigest() != row["stderr_sha256"]:
            raise RuntimeError(f"Original error log changed for task {row['task_id']}")
        status = rows(Path(row["status_path"]))
        if len(status) != 1 or status[0]["task_id"] != row["task_id"] or (
            status[0]["state"] != "COMPLETE"
        ):
            raise RuntimeError(f"Recovered task {row['task_id']} is not COMPLETE")
        outdir = Path(row["output_dir"])
        outputs = [outdir / name for name in RAW_FILES]
        outputs.append(outdir.parent / f"{row['patient']}.Rds")
        if any(not path.is_file() or path.stat().st_size <= 0 for path in outputs):
            raise RuntimeError(f"Recovered task {row['task_id']} lacks fit output")
        if (outdir / "xval_recovered.txt").exists():
            raise RuntimeError(f"Task {row['task_id']} used the old xval workaround")
    print(f"VERIFIED {len(evidence)} targeted fits with complete status and raw/flat output")


if __name__ == "__main__":
    main()
