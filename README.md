# Penn CarT ALFA-K analysis

This repository runs the four paired-resection GSE296419 patients (P1, P4,
P7, P13) through ALFA-K, outcome-blind PM/MINOBS selection, and the compatible
PANcanKFLs downstream analysis. `high_cn_6` is primary; `high_cn_5` and
`high_cn_8` are separate chromosome-state sensitivity analyses. These are
RNA-derived, diploid-relative HMM summaries, not validated absolute DNA copy
numbers. Time zero is the first resection and the interval includes time before
and after CAR-T infusion.

The HPC result root is exactly
`/share/lab_crd/taoli/Project/CarTKFLs/results`. Source RDS files remain
immutable in the CarTData handoff. The workflow copies and hashes them into
the result root for the PANcanKFLs consumer contract.

## HPC run

After cloning/pulling this repository on RED, run:

```bash
bash scripts/submit.sh
```

The controller submits preflight as a Slurm job, waits for its result, and
submits one fit task per patient, high-CN mapping, PM, and MINOBS combination.
The fit arrays use `xxlarge`, 12 hours, one CPU, and separate memory requests
of 48/16/8 GB for MINOBS 5/10/20. Each array has a concurrency cap of 32.
The controller then submits a dependent audit and, if it passes, the KFL and
downstream stages. `results/runs/<run_id>/submissions.tsv` records all job IDs.

PANcanKFLs code is a pinned external source dependency at the commit in
`config/analysis.yaml`; no files in that repository are changed. The specified
alfakR dev2 SIF is also pinned by SHA-256. Failed ALFA-K combinations remain
visible in task status and are not automatically retried. The final selection
uses KFL accuracy/stability with `selection_mode=kfl_only`; no clinical
outcome, ABM, or expression annotation is used to select parameters.

The pinned alfakR dev2 `xval()` has a confirmed result-column naming defect:
some Krig predictions occupy the second column under a name other than
`est_f`, while the package reads `est_f` by name. `scripts/xval_compat.R`
retains the installed function body and changes only that column lookup to
position 2, failing if the expected function body changes. When ALFA-K has
already written a valid bootstrap, landscape, and posterior but failed at this
lookup, `fit_one.R` finishes cross-validation from the saved bootstrap. The
repair scripts do the same for earlier failed tasks without rerunning fitting,
and preserve each original status in `*.before_xval_repair.tsv`. The dev2
`xval_data` fields are extracted directly into PANcanKFLs' patient-fit format.

Outputs beneath `results/` are ignored by Git. Inspect `submissions.tsv`,
`fit_audit.tsv`, `selection_summary.tsv`, and Slurm accounting for progress.
