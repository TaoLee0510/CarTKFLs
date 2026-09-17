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
