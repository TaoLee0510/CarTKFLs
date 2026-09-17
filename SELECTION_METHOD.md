# Temporary Pearson-correlation parameter selection

This note applies to GSE296419 patients P1, P4, P7, and P13 in the existing
`20260917_f93b2d2_uncapped_f4633d2` run. The active downstream analysis under
`correlation_only/` **does not follow the original PANcanKFLs
KFL accuracy/stability PM/MINOBS selection procedure**. It is an explicitly
exploratory, outcome-blind alternative requested after the original selection
failed to cover all four patients. The original ALFA-K fits and KFL metrics are
reused; no fitness model is refitted for this alternative.

## What the original procedure produced

The pinned PANcanKFLs code applies all of these gates to each PM/MINOBS
candidate: finite positive `signed_r2_Huber_weighted`, Pearson
`p_value < 0.05`, permutation `huber_perm_p < 0.05`, Huber slope within
`[0.6, 1.6667]`, relative intercept `< 0.3`, outlier ratio `<= 0.2`,
and bias score `< 2`. Among candidates that pass, it uses stability and
posterior scores to choose one tuple per patient/mapping pair.

The original fit audit accepted 8,820 fitted grid combinations and recorded
72 model errors. The strict KFL rule found 10 qualifying rows across only
five of 12 patient/mapping pairs. It selected:

| Mapping | Patient | Original qualifying rows | Original selected PM | Original MINOBS |
| --- | --- | ---: | --- | ---: |
| high_cn_5 | P7 | 2 | pm_0.00325 | 20 |
| high_cn_6 | P4 | 1 | pm_0.0005 | 10 |
| high_cn_6 | P7 | 1 | pm_0.0122 | 5 |
| high_cn_8 | P1 | 1 | pm_0.01025 | 5 |
| high_cn_8 | P7 | 5 | pm_0.0055 | 10 |

There was no qualifying row for high_cn_5/P1, P4, or P13;
high_cn_6/P1 or P13; or high_cn_8/P4 or P13. Thus the strict procedure
cannot produce a primary high_cn_6 result for P1 or P13, and cannot produce
any mapping result for P13. The seven absent pairs have evaluable KFL rows;
they are not missing fit jobs. Most rows that passed the two significance
tests and positive signed metric were then rejected by the positive Huber
slope interval. high_cn_8/P4 had no row that simultaneously passed those
first three gates.

The source run retains the original `manifests/selected_samples.tsv`,
`final_parameters_opt.Rds`, and KFL metric files without modification.
Its select job failed because a CarTKFLs summary script incorrectly required
all 12 pairs to pass the strict rule; the PANcanKFLs selection itself
produced the five rows listed above.

## Current temporary rule

For each patient and each high-copy mapping independently, consider every
completed PM/MINOBS row with a finite `correlation` value in the existing
KFL output. This column is the Pearson correlation between ALFA-K's
cross-validation estimated fitness `f_est` and held-out fitness `f_xv`.
Choose the row with the **largest Pearson correlation**. An exact numerical
tie is resolved by smaller PM, then smaller MINOBS, then lexical PM label.
No `p_value`, permutation, slope, intercept, outlier, bias, stability, or
clinical outcome threshold is used to choose the tuple. The selected fit
file must exist before downstream submission.

This rule selects all 12 pairs:

| Mapping | Patient | Pearson-only PM | MINOBS | Pearson r | Passes original rule? |
| --- | --- | --- | ---: | ---: | --- |
| high_cn_5 | P1 | pm_0.0055 | 20 | 0.733 | No |
| high_cn_5 | P4 | pm_0.009 | 10 | 0.725 | No |
| high_cn_5 | P7 | pm_0.0093 | 20 | 0.866 | No |
| high_cn_5 | P13 | pm_0.007 | 5 | 0.829 | No |
| high_cn_6 | P1 | pm_0.00125 | 20 | 0.807 | No |
| high_cn_6 | P4 | pm_0.0005 | 10 | 0.788 | Yes |
| high_cn_6 | P7 | pm_0.0106 | 20 | 0.704 | No |
| high_cn_6 | P13 | pm_0.0052 | 20 | 0.830 | No |
| high_cn_8 | P1 | pm_0.0086 | 20 | 0.911 | No |
| high_cn_8 | P4 | pm_0.00885 | 20 | 0.549 | No |
| high_cn_8 | P7 | pm_0.01215 | 20 | 0.834 | No |
| high_cn_8 | P13 | pm_0.00555 | 10 | 0.893 | No |

These values were read from the frozen KFL output before submitting the
Pearson-only downstream work. The selection job independently recomputes
them from the saved RDS metrics and writes `selection_comparison.tsv` and
`selection_summary.tsv` for machine-readable provenance.

## Interpretation

Pearson correlation measures concordance in ordering and linear trend, not
agreement in the absolute fitness values consumed downstream. For example,
the highest-correlation high_cn_6/P13 row has `r = 0.830` from nine finite
cross-validation points, but weighted predictive R² is `-17.30` and the
relative intercept is `1.30`. The highest-correlation high_cn_6/P1 row has
`r = 0.807` but a Huber slope of `2.48`, beyond the original upper bound.
Eleven of the 12 Pearson maxima fail at least one original accuracy gate.
Taking the maximum over hundreds of grid rows can also make the observed
correlation optimistic. These downstream outputs must therefore be labeled
as exploratory Pearson-only analyses, not as results that passed the
original KFL accuracy/stability selection.

The input chromosomes are unfiltered RNA-derived HMM summaries, not
validated absolute DNA copy numbers. Time zero is the first tumor resection;
the inferred interval includes time before and after CAR-T infusion.
