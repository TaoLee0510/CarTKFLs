# Wang GBM–ovarian T-cell harmonization

This workflow harmonizes T/NK annotations across the Wang GBM and ovarian Seurat objects and evaluates whether chromosome-fitness associations with total T-cell fraction are sensitive to T-cell subtype composition or state.

## Annotation contract

- Mutually exclusive subtype: `CD8`, `Conventional_CD4`, `Treg`, `Unresolved_T`; NK is retained as a distinct lineage/subtype (`NK`).
- Continuous programs: cytotoxicity, activation, naive/memory, exhaustion/dysfunction, interferon response, and stress.
- UCell scores are the primary state variables. Robust z-scores and primary/supported-state labels are secondary summaries.
- Ovarian labels map from the approved hierarchical annotation. GBM candidates are the 924 cells in the final coarse `Normal_celltype == "T cell"` category; they are re-clustered and assigned by robust ovarian subtype-reference centroids. Failed distance, lineage, margin, or size checks remain unresolved.
- Marker genes, citations, coverage, assignment source, confidence, and uncertainty reasons are written to the result directory.

## Execution contract

The workflow must run directly on `hpctpa3pc0028` inside the configured SIF. It does not invoke Slurm. `hpc/run_direct_hpctpa3pc0028.sh` starts exactly one detached `nohup` process after the run-specific UCell library has been prepared.

Completed harmonization outputs are checkpointed. A later attempt in the same run directory reuses valid GBM/ovarian objects and cell metadata instead of recomputing or rewriting them.

Pipeline order:

1. `01_audit_inputs.R` (pre-run audit, already separable from the formal pipeline)
2. `02_harmonize_annotations.R`
3. `03_qc_and_comparison.R`
4. `04_fitness_association.R`
5. `05_build_report.R`
6. `06_final_audit.R` (post-run contract, result inventory, and key-findings summary)
7. `07_validate_objects.R` (sequential readback validation of the two large annotated Seurat objects)

## Interpretation

Chromosome summaries are marginal effects derived from each selected ALFA-K patient landscape. Cross-patient correlations and one-feature adjustment models are exploratory, use patient—not cell—as the statistical unit, and are not causal mediation tests.
