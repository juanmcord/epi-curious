# Validation record — epi-curious 0.1.0

Validation performed on 2026-09-30 with R 4.3.3 on Ubuntu 24.04.3, x86_64 Linux.
This is validation of the new implementation against standard R models on
synthetic data. It is not a reproduction of a real epi-curious study analysis.

## Completed

| Check | Result |
|---|---|
| Source package build | Passed: `R CMD build epicurious`. |
| Source installation | Passed, including temporary and final namespace loading. |
| Standard package check | `R CMD check --no-manual epicurious_0.1.0.tar.gz`: **Status: OK**, zero errors, warnings, or notes. |
| Numerical and behavior regression suite | **34 of 34 checks passed**, including the suite run within the package check. |
| Installed documentation examples | Passed as part of the package check. |
| Source-tree loader | Loaded functions successfully; synthetic EWAS completed. |
| Single-file standalone script | Loaded independently and completed its synthetic EWAS. |
| CLI demo | Completed a 200-sample, 12-CpG analysis and exported results. |
| CLI configuration route | Read all four RDS inputs using paths relative to the configuration file; exported 12 successful CpG models. |
| Export checks | RDS object round-trip, expected output files, overwrite protection, and removal of an obsolete paired-comparison CSV passed. |
| Parallel execution | Two-worker results matched sequential results on Linux. |

The statistical reference checks compare estimates, standard errors, confidence
intervals, and primary p-values with `lm()`. Multilevel exposure F tests are
compared with nested-model `anova()`. The comparison tolerance is `1e-8` in
R's `all.equal()`; this is a numerical-equivalence test, not a speed benchmark.

## Cases covered by the regression suite

- Zero, one, and multiple SNPs per CpG; CpG-specific SNP assignment.
- Different input sample orders and explicit intersected sample sets.
- Simultaneous missing covariates, methylation, and genotypes, with correct
  per-CpG sample counts and same-sample SNP-unadjusted comparisons.
- Multilevel exposure factors and single-coefficient continuous exposures.
- Duplicate sample/CpG IDs, unknown requested CpGs, duplicate map pairs,
  unavailable SNPs under each policy, and irrelevant mapping entries.
- Fractional genotype dosage, invalid dosage values, and punctuated SNP IDs.
- Exact SNP collinearity, monomorphic SNPs, constant methylation, insufficient
  sample counts, insufficient residual degrees of freedom, and all-failed runs.
- M and beta response scales, explicitly unsupported formulas, and single-CpG
  matrix subsetting.
- Multiple-testing family size including failed models.
- Wide-table conversion, a matrix-backed bigSNP adapter fixture, reproducible
  simulation, coefficient-retention options, and different output chunk sizes.

## Limits of this validation

- No original `.R` source or real input dataset was supplied. Screenshot-level
  analysis cannot establish exact reproduction of the original study results.
- The optional bigSNP adapter was tested against a matrix-backed structural
  fixture. Its interface was checked against official bigsnpr documentation,
  but a real file-backed bigsnpr dataset was not available for integration testing.
- Windows/macOS execution, other R versions, and genome-wide runtime/memory
  performance were not tested. The package explicitly rejects parallel workers
  on Windows and defaults to sequential execution everywhere.
- PDF manual generation was excluded with `--no-manual`; R help and examples
  were checked. The delivery does not claim a CRAN release or submission review.
- The implementation provides conventional OLS inference. No real-cohort
  residual diagnostics, genomic calibration, sample QC, or model-choice
  validation can be performed without study data.

The bundle includes `validation/regression.log` and `validation/R-CMD-check.log`
as evidence of the executed checks. Runtime setup files and original image
inspection intermediates are not part of the deliverable.
