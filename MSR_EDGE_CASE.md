# MSR threshold edge case in the Pearson-only branch

The Pearson-only downstream analysis completed all 12 patient-level MSR
steady-state calculations. The high_cn_5 and high_cn_6 cancer-level MSR
reductions completed with the pinned PANcanKFLs code. The high_cn_8 reduction
stopped while computing six error-threshold estimates for P7.

The cached high_cn_8/P7 original-karyotype steady-state curve has its minimum
at the **first** assayed error rate, `p = 0`. The PANcanKFLs threshold
function fits a line to the curve from the first rate through its minimum.
Here that interval contains only one point, so its regression slope is
undefined. The original function tests `if (b != 0)` when `b` is `NA`,
which causes an R error.

For this one verified edge case, `scripts/msr_reduce_first_minimum.R`
records all six pre-minimum threshold estimates as `NA`; no threshold is
imputed or inferred from a different interval. Every curve with at least two
points on the pre-minimum interval still uses the unmodified PANcanKFLs
function. The script reruns only the high_cn_8 cancer-level MSR reduction
from its completed patient-level caches, writes `msr_edge_case.tsv`, and
reruns the standard PANcanKFLs MSR audit afterward.

An `NA` threshold for this curve means the requested pre-minimum threshold
is not estimable on this error-rate grid. It should not be interpreted as a
threshold of zero or as evidence that a threshold exists.
