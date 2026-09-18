# GSE210079 multiple myeloma: temporary parameter selection

This note describes the P16, P32, and P33 ALFA-K analysis under
`multiple_myeloma/runs/20260917_175858_7d2f0257/`. The active downstream
outputs **do not use the original PANcanKFLs accuracy and stability parameter
selection**. They use a temporary, exploratory Pearson-correlation rule. The
2,223 completed ALFA-K fits and their KFL evaluation rows are reused; no
fitness model is refitted for this change. The selection, downstream results,
and audits occupy the existing run paths. `selection_comparison.tsv` records
the per-patient comparison in machine-readable form.

## What the original process would select

The pinned PANcanKFLs implementation first requires every candidate to have
a finite positive signed Huber-weighted R², Pearson `p_value < 0.05`, Huber
permutation `p < 0.05`, Huber slope in `[0.6, 1.6667]`, relative intercept
`< 0.3`, outlier ratio `<= 0.2`, and bias score `< 2`. It then combines
accuracy with stability and posterior scores to rank qualifying candidates.

All three patients have 741 completed PM/MINOBS combinations. Under the
original process, P16 has **zero** qualifying candidates, P32 has **six**,
and P33 has **zero**. The original selection therefore chooses only
P32 `pm_0.0001 / MINOBS_20`. It cannot provide selected-parameter downstream
results for P16 or P33. Their absence under that rule reflects selection
failure, not missing ALFA-K fits.

## Temporary rule used for the active results

For each patient, use the `correlation` field in the saved KFL evaluation
rows. PANcanKFLs computes it as the Pearson correlation between ALFA-K's
cross-validation fitness estimate `f_est` and held-out value `f_xv`. Among
rows with finite correlation, choose the **largest signed Pearson r**. Exact
ties use smaller PM, then smaller MINOBS, then lexical PM label. We do not
require the original significance, permutation, slope, intercept, outlier,
bias, stability, or posterior gates. Clinical outcomes are not used in the
choice. The selected ALFA-K fit file must exist before downstream execution.

| Patient | Pearson-selected PM | MINOBS | Pearson r | CV points | Predictive R² | Passes original rule? |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| P16 | `pm_0.00335` | 20 | 0.999871 | 3 | -24.025 | No |
| P32 | `pm_0.000005` | 20 | 0.916335 | 13 | 0.059 | No |
| P33 | `pm_0.00505` | 5 | 0.99999978 | 3 | -60.205 | No |

These values come from the completed KFL grid before resubmission and are
checked by the selection job. The previous temporary rule ranked predictive
R²; its tuples are superseded in the same output paths. The original strict
choice for P32 also differs from the Pearson choice.

## Interpretation limits

Pearson r measures a linear trend, not agreement in the **absolute fitness
values** used downstream. In particular, the Pearson winners for P16 and
P33 have near-perfect r from only three cross-validation points but strongly
negative predictive R². Even P32's predictive R² is only about 0.06. All
three Pearson winners fail at least one original accuracy gate. Selecting
the maximum from 741 combinations per patient also makes the observed r
optimistic; the reported nominal Pearson p-values are not adjusted for that
search. These downstream results are exploratory and must not be described
as having passed the original PANcanKFLs parameter-selection procedure.

P33's Numbat tumor assignment has separate cell-identity uncertainty.
Sampling times are day 0, day 28 where available, and nominal day 90 for
month 3, with ALFA-K `dt=1` day. Completion should be checked in the run's
`final_audit.tsv`, `final_cancer_audit.tsv`, `msr_sample_audit.tsv`, and
`msr_cancer_audit.tsv` rather than inferred from selection alone.
