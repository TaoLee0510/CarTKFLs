#!/usr/bin/env python3
"""Identify only ALFA-K tasks lost before container startup."""

import argparse
import csv
import hashlib
import re
import subprocess
import tempfile
from pathlib import Path


ERROR = re.compile(r"FATAL:\s+user: unknown userid 107865")
RAW_FILES = (
    "bootstrap_res.Rds",
    "landscape.Rds",
    "landscape_posterior_samples.Rds",
    "xval.Rds",
)


def rows(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def command(*args):
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout


def write_tsv(path, fields, records):
    if path.exists():
        raise RuntimeError(f"Refusing to overwrite {path}")
    with tempfile.NamedTemporaryFile(
        mode="w", newline="", dir=path.parent, prefix=f".{path.name}.",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    run_dir = args.run_dir.resolve(strict=True)
    if not (run_dir / "fit_audit.tsv").is_file():
        raise RuntimeError("Original fit audit must finish before recovery")

    _, submissions = rows(run_dir / "submissions.tsv")
    jobs = {}
    manifests = {}
    manifest_fields = None
    for min_obs in (5, 10, 20):
        matching = [row for row in submissions
                    if row["stage"] == f"fit_MINOBS_{min_obs}"]
        if len(matching) != 1:
            raise RuntimeError(f"Expected one original MINOBS {min_obs} array")
        entry = matching[0]
        if not entry["manifest"].endswith("#offset=0"):
            raise RuntimeError("Recovery requires zero-offset original arrays")
        jobs[entry["job_id"]] = min_obs
        fields, manifest = rows(Path(entry["manifest"].split("#", 1)[0]))
        if manifest_fields is None:
            manifest_fields = fields
        elif fields != manifest_fields:
            raise RuntimeError("Original task manifests have different columns")
        manifests[min_obs] = manifest

    audits = [row for row in submissions if row["stage"] == "audit"]
    if len(audits) != 1:
        raise RuntimeError("Expected exactly one original audit submission")
    parent_ids = list(jobs) + [audits[0]["job_id"]]
    if command("squeue", "-r", "-h", "-j", ",".join(parent_ids),
               "-o", "%A").strip():
        raise RuntimeError("Original fit arrays or audit are still active")

    parent_states = {}
    failed_start = {}
    for line in command(
        "sacct", "-X", "-j", ",".join(parent_ids), "--noheader",
        "--parsable2", "-o", "JobID,State,ExitCode,NodeList",
    ).splitlines():
        job_id, state, exit_code, node, *_ = line.split("|")
        if job_id in parent_ids:
            parent_states[job_id] = state
            continue
        match = re.fullmatch(r"(\d+)_(\d+)", job_id)
        if match and match.group(1) in jobs and state == "FAILED" and exit_code == "127:0":
            failed_start[(match.group(1), int(match.group(2)))] = node
    if set(parent_states) != set(parent_ids):
        raise RuntimeError(f"Missing Slurm parent accounting: {parent_states}")
    if any(parent_states[job] != "COMPLETED" for job in jobs):
        raise RuntimeError(f"Original fit arrays are not complete: {parent_states}")
    if parent_states[audits[0]["job_id"]] not in ("COMPLETED", "FAILED"):
        raise RuntimeError(f"Original audit is not terminal: {parent_states}")

    _, audit = rows(run_dir / "fit_audit.tsv")
    if len(audit) != 8892:
        raise RuntimeError(f"Expected 8892 audit rows; observed {len(audit)}")
    missing = {row["task_id"]: row for row in audit
               if row["recorded_state"] == "MISSING_STATUS"}
    _, inputs = rows(run_dir / "input_provenance.tsv")
    staged_inputs = {(row["patient"], row["high_cn"]): row["staged_path"]
                     for row in inputs}
    evidence = []
    recovery_tasks = {min_obs: [] for min_obs in (5, 10, 20)}
    for (job_id, array_index), node in sorted(failed_start.items()):
        min_obs = jobs[job_id]
        manifest = manifests[min_obs]
        if array_index < 1 or array_index > len(manifest):
            raise RuntimeError(f"Invalid array index {job_id}_{array_index}")
        task = manifest[array_index - 1]
        task_id = task["task_id"]
        if task_id not in missing:
            raise RuntimeError(f"Exit-127 task {task_id} is not MISSING_STATUS")
        audited = missing[task_id]
        if (task["patient"] != audited["patient"]
                or task["high_cn"] != audited["high_cn"]
                or task["min_obs"] != audited["min_obs"]
                or task["pm_label"] != audited["pm_label"]):
            raise RuntimeError(f"Manifest and audit disagree for task {task_id}")
        if audited["raw_valid"] != "FALSE" or audited["flat_valid"] != "FALSE":
            raise RuntimeError(f"Task {task_id} has valid fit output")
        status_path = Path(audited["status_path"])
        outdir = Path(audited["output_dir"])
        flat_path = outdir.parent / f'{task["patient"]}.Rds'
        if status_path.exists() or flat_path.exists() or any(
            (outdir / name).exists() for name in RAW_FILES
        ):
            raise RuntimeError(f"Task {task_id} already has a status or output")
        stderr = run_dir / "logs" / f"fit{min_obs}_{job_id}_{array_index}.err"
        if not stderr.is_file() or not ERROR.search(stderr.read_text(errors="replace")):
            raise RuntimeError(f"No exact container startup error in {stderr}")
        staged_input = staged_inputs.get((task["patient"], task["high_cn"]))
        if not staged_input or not Path(staged_input).is_file():
            raise RuntimeError(f"Staged input missing for task {task_id}")
        evidence.append({
            "task_id": task_id, "patient": task["patient"],
            "high_cn": task["high_cn"], "pm": task["pm"],
            "pm_label": task["pm_label"], "min_obs": task["min_obs"],
            "original_job_id": job_id, "original_array_index": array_index,
            "original_state": "FAILED", "original_exit_code": "127:0",
            "node": node, "stderr_path": str(stderr),
            "stderr_sha256": hashlib.sha256(stderr.read_bytes()).hexdigest(),
            "status_path": str(status_path), "staged_input_path": staged_input,
            "output_dir": str(outdir),
        })
        recovery_tasks[min_obs].append(task)
    if {row["task_id"] for row in evidence} != set(missing):
        extra = set(missing) - {row["task_id"] for row in evidence}
        raise RuntimeError(f"Missing statuses without exact startup evidence: {sorted(extra)}")
    if not evidence:
        raise RuntimeError("No container startup failures require recovery")

    print(f"Container startup failures: {len(evidence)}")
    for min_obs in (5, 10, 20):
        print(f"MINOBS {min_obs}: {len(recovery_tasks[min_obs])}")
    print("Original task IDs:", ",".join(row["task_id"] for row in evidence))
    if args.write:
        directory = run_dir / "manifests"
        files = [directory / "container_start_recovery_evidence.tsv"]
        files += [directory / f"container_start_recovery_MINOBS_{min_obs}.tsv"
                  for min_obs in (5, 10, 20) if recovery_tasks[min_obs]]
        if any(path.exists() for path in files):
            raise RuntimeError("Recovery manifest already exists; refusing duplicate submission")
        write_tsv(files[0], list(evidence[0]), evidence)
        for min_obs in (5, 10, 20):
            if recovery_tasks[min_obs]:
                write_tsv(directory / f"container_start_recovery_MINOBS_{min_obs}.tsv",
                          manifest_fields, recovery_tasks[min_obs])


if __name__ == "__main__":
    main()
