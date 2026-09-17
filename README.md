# Penn CarT ALFA-K analysis

This repository runs the four paired-resection GSE296419 patients (P1, P4,
P7, P13) through ALFA-K, outcome-blind PM/MINOBS selection, and the compatible
PANcanKFLs downstream analysis. `high_cn_6` is primary; `high_cn_5` and
`high_cn_8` are separate chromosome-state sensitivity analyses. These are
RNA-derived, diploid-relative HMM summaries, not validated absolute DNA copy
numbers. Time zero is the first resection and the interval includes time before
and after CAR-T infusion.

The HPC result root is exactly
`/share/lab_crd/taoli/Project/CarTKFLs/results`. Each new run writes its
inputs, ALFA-K fits, and downstream results under
`results/runs/<run_id>/analysis`, separate from the canceled run. Source RDS
files remain immutable in the CarTData handoff. The workflow copies and hashes
them into the run directory for the PANcanKFLs consumer contract.

## HPC run

After cloning/pulling this repository on RED, run:

```bash
bash scripts/submit.sh
```

The controller submits preflight as a Slurm job, waits for its result, and
submits one fit task per patient, high-CN mapping, PM, and MINOBS combination.
The fit arrays use `xxlarge`, 12 hours, one CPU, and separate memory requests
of 48/16/8 GB for MINOBS 5/10/20. No per-array concurrency limit is set;
Slurm schedules tasks according to available resources.
The controller then submits a dependent audit and, if it passes, the KFL and
downstream stages. `results/runs/<run_id>/submissions.tsv` records all job IDs.

PANcanKFLs code is a pinned external source dependency at the commit in
`config/analysis.yaml`; no files in that repository are changed. The specified
alfakR dev2 SIF is also pinned by SHA-256. Failed ALFA-K combinations remain
visible in task status and are not automatically retried. The final selection
uses KFL accuracy/stability with `selection_mode=kfl_only`; no clinical
outcome, ABM, or expression annotation is used to select parameters.

The pinned alfakR dev2 image contains package commit `f93b2d2`, which fixes
the `xval()` prediction-column naming defect. This run calls the installed
package directly and fits every combination from the staged input; it does not
recover or reuse output from the canceled run. The earlier compatibility and
repair scripts remain in the repository solely for historical provenance. The
dev2 `xval_data` fields are extracted directly into PANcanKFLs' patient-fit
format.

Outputs beneath `results/` are ignored by Git. Inspect `submissions.tsv`,
`fit_audit.tsv`, `selection_summary.tsv`, and Slurm accounting for progress.
`selection_coverage.tsv` lists every patient/mapping pair and its count of
significant KFL candidates. Pairs without a qualifying candidate are recorded
in `selection_exclusions.tsv` and receive no selected PM/MINOBS tuple.

## Recovery of container startup failures

For the uncapped `20260917_f93b2d2_uncapped_f4633d2` run only, after its three
fit arrays and original audit are terminal, run
`scripts/submit_container_start_recovery.sh` with that existing run directory.
The builder requires an exact `unknown userid 107865` container error,
Slurm exit 127, a missing task status, and no saved fit output. It writes an
evidence manifest and submits only matching original task rows with their
original MINOBS memory request and no array throttle. The original audit is
copied before any recovery work. Run `scripts/verify_container_start_recovery.py`
after the recovery arrays finish and before submitting a fresh fit audit in
the same run directory; that audit
continues into KFL selection and downstream analysis if it passes.

## Temporary Pearson-only downstream branch

For the existing `20260917_f93b2d2_uncapped_f4633d2` run, the alternative
selection requested after the strict KFL screen is documented in
[`SELECTION_METHOD.md`](SELECTION_METHOD.md). Submit its isolated downstream
branch with `bash scripts/submit_correlation_downstream.sh` on RED. It writes
under that same run's `correlation_only/` directory, links the completed
ALFA-K inputs and KFL sample results, and leaves the original strict-selection
artifacts intact. The branch uses `xxlarge` and a 12-hour limit for its Slurm
jobs. `selection_comparison.tsv` records both methods for all 12 pairs.
