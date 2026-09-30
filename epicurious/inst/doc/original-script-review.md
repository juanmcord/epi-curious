# Review of the original epi-curious script

## What the screenshot shows

The supplied screenshot contains a study-specific R workflow that regresses
methylation at individual CpGs on phenotype/covariate terms, optionally adding
predefined SNPs for each CpG. The CpG–SNP associations are read from a mapping
file; the script does not discover mQTLs in this step.

The image is 553 × 2048 pixels. Function structure and major operations are
readable, but some fine text, paths, and variable names are not unambiguous.
This rewrite therefore preserves the visible modelling approach, rather than
claiming a byte-for-byte transcription. The original R source and study inputs
were not available, so exact reproduction of the original study output was
not possible.

The visible workflow is:

1. Load packages including data.table, bigsnpr, parallel, and stats; define study
   folders, filenames, sample-ID column, formulas, and an EWAS output name.
2. Read the phenotype/covariate data and wide sample-by-CpG methylation data.
   Intersect sample sets, attach a bigsnpr genotype dataset, and align rows.
3. Read a CpG–SNP map of predefined, apparently independent SNPs. The visible
   default selection appears to intersect methylation column names with mapped
   CpGs before constructing the CpG list. It should not be assumed to analyze
   every array CpG by default.
4. Build the phenotype/covariate design matrix once.
5. For each CpG, obtain its methylation vector, look up its SNPs, extract genotype
   columns, append them to the shared design matrix, and attempt a linear fit.
6. Call R's internal QR least-squares routine through `stats:::` and `.Call`,
   then calculate coefficient estimates, residual variance, standard errors,
   t statistics, and two-sided p-values manually.
7. Run CpGs with a hard-coded multicore setting, reshape coefficient lists
   into a wide table, add an external CpG annotation table, and save the result.

For a CpG j, the intended model is:

```text
methylation_ij = intercept_j
               + exposure_i × exposure_effect_j
               + covariates_i × covariate_effects_j
               + sum(SNP_dosage_ik × SNP_effect_jk for k mapped to CpG j)
               + error_ij
```

This makes the exposure effect conditional on the selected SNPs and covariates.
It does not residualize methylation as a separate preprocessing step. The model
has the same scientific structure whether the response is an M value or beta
value, but the scale and inferential behavior differ. The rewrite requires an
explicitly documented scale and does no automatic transformation.

## Findings that matter for reuse

| Finding in the visible script | Why it matters | Rewrite behavior |
|---|---|---|
| Paths, study variables, output names, genotype backing files, and worker count are embedded in the script. | Reuse requires editing analysis internals and assumes a particular machine. | Function arguments plus an optional configuration file. No global working-directory change or genotype backing-file writes. |
| Sample order is established through several intersections and matches. | Duplicate IDs or an incomplete match can corrupt row alignment. | Unique-ID validation, explicit strict/intersection policies, deterministic reordering, and sample audit. |
| The shared model matrix is created before CpG-specific missing-value filtering. | Default model-frame NA handling can remove covariate rows independently of the methylation vector. | Preserve missing rows during design construction, then form one joint complete-case mask. |
| When CpG methylation has missing values, its response/design are shortened before full-length SNP columns are appended. | The genotype block can have a different number of rows, causing an error or incorrect alignment. | Append aligned SNP columns first, then subset response and every design column together. |
| SNP inclusion is guarded by a permissive availability check, with an unchanged-model branch when none are found. | A mapped CpG can receive an unadjusted fit without an explicit adjustment-status record. Partial matches also need careful handling. | Missing mapped SNPs are an error by default; explicit alternatives report skipped or reduced-adjustment models. |
| A per-CpG formula is updated, but the visible fit uses the numeric design matrix. | The formula text is not what actually controls the QR fit; they can diverge. | A single shared design matrix plus a directly appended, audited SNP block. |
| Internal R QR symbols and manual coefficient extraction are used. | Reliance on an internal interface complicates maintenance and error handling. | Public `stats::lm.fit()` with tested QR-based inference. |
| SNP names are rewritten with a prefix while extracting genotype columns. | Requested and actually available SNP sets must stay in the same order and have the same size. | Exact ID matching and explicit matrix indexing; original variant IDs remain in the audit. |
| Coefficients are widened over the union of variables/SNPs. | Different CpGs have different SNPs; a genome-wide table can become mostly empty and very large. | Compact CpG results plus a long coefficient table, exposure-only by default. |
| The final visible p-values are coefficient-level raw p-values. | A reusable EWAS needs a clearly defined testing family and exposure test. | One primary exposure-term test per CpG, BH or another chosen correction, and Bonferroni values. |
| Annotation comes from a fixed external file. | Annotation is platform/build dependent and can duplicate result rows on a many-to-many join. | Annotation is an explicit downstream join controlled by the user. |

## Intentional changes from the original workflow

The rewrite analyzes **all supplied methylation CpGs by default**, including
unmapped sites. To reproduce a mapping-restricted analysis, use:

```r
cpgs_to_test <- intersect(rownames(methylation), unique(cpg_snp_map$cpg_id))
fit <- ewas_mqtl(
  methylation, phenotypes, ~ exposure + age + sex, "exposure",
  genotypes, cpg_snp_map, cpgs = cpgs_to_test, methylation_scale = "M"
)
```

It also defaults to one worker, fails on unavailable mapped SNPs, records skipped
CpGs instead of silently discarding them, and uses a single primary exposure-term
test. Multiple independent SNPs remain separate additive covariates; no pruning,
variable selection, genotype imputation, or empirical-Bayes shrinkage is added.

The original all-coefficient output is available in long form with
`keep_coefficients = "all"`. The primary results table gives the exposure result
most EWAS users need. For a multilevel factor, the primary test is joint across
its design columns, with individual contrasts in the coefficient table.

To add a unique annotation table while preserving row order:

```r
stopifnot(!anyDuplicated(annotation$cpg_id))
idx <- match(fit$results$cpg_id, annotation$cpg_id)
extra_cols <- setdiff(names(annotation), "cpg_id")
stopifnot(!any(extra_cols %in% names(fit$results)))
annotated <- cbind(fit$results, annotation[idx, extra_cols, drop = FALSE])
```

No annotation resource is guessed or downloaded. Use the study's appropriate
array platform and genome build. Missing annotations remain missing rather than
removing the corresponding CpG result.

## Interpretation and next validation step

This is an implementation of conditional association analysis, not a guarantee
that genotype is the only explanation for a methylation association. The SNP
mapping, cell composition, ancestry adjustment, phenotype coding, and sample
design remain scientific inputs to be specified before a real analysis.

The synthetic regression suite establishes agreement with ordinary R models
under supported conditions. The remaining study-specific check is to supply
the original R text and a de-identified input subset, choose the same CpGs,
formula, scale, samples, and SNP map, and reconcile each model against the
original expected output. That comparison cannot be made from a screenshot
alone.
