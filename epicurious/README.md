# epi-curious

**A reusable R workflow for EWAS with predefined, CpG-specific SNP adjustment.**

For each CpG, epi-curious fits:

```text
methylation ~ exposure + covariates + SNPs mapped to this CpG
```

The exposure is the association of interest. Mapped SNP dosages are adjustment
variables. An unmapped CpG receives the exposure-and-covariate model. The supplied
mapping determines adjustment; this software does not discover or select mQTLs.

This is version 0.1.0: an installable R package, a standalone R script, and a
configuration-driven command-line workflow. The core needs only packages shipped
with R. See `inst/doc/original-script-review.md` for the screenshot analysis and
the intentional changes from the original workflow.

## Start with the example

Requires R 4.1 or later. From the extracted bundle, choose either approach.

**Install the package:**

```r
install.packages("epicurious_0.1.0.tar.gz", repos = NULL, type = "source")
library(epicurious)
```

**Or load the standalone file:**

```r
source("epicurious.R")
```

Then run:

```r
d <- simulate_epicurious()

fit <- ewas_mqtl(
  methylation = d$methylation,
  phenotypes = d$phenotypes,
  formula = ~ exposure + age + sex,
  exposure = "exposure",
  genotypes = d$genotypes,
  cpg_snp_map = d$cpg_snp_map,
  methylation_scale = "M",
  compare_unadjusted = TRUE
)

print(fit)
fit$results[, c("cpg_id", "estimate", "p_value", "p_adjust", "status")]
write_ewas(fit, "results/example")
```

The simulation includes unmapped CpGs, one-SNP and two-SNP adjustments, known
exposure effects, and a CpG with a SNP effect but no true exposure effect.
`write_ewas()` refuses to overwrite an existing analysis unless explicitly asked.

## Supply your own data

| Input | Required structure | Meaning |
|---|---|---|
| `methylation` | Numeric matrix: **CpGs × samples** | Row names are exact CpG IDs; column names are unique sample IDs. |
| `phenotypes` | Data frame: **samples × variables** | Contains a unique `sample_id` column and every formula variable. |
| `genotypes` | Numeric matrix: **samples × SNPs** | Row names are sample IDs; column names are exact SNP IDs; values are dosage 0–2 or `NA`. |
| `cpg_snp_map` | Data frame with `cpg_id`, `snp_id` | One row per predefined CpG–SNP pair; multiple rows per CpG are allowed. |

For example, a map can contain:

| cpg_id | snp_id |
|---|---|
| cg00000001 | rs1001 |
| cg00000001 | rs1002 |
| cg00000002 | rs2001 |

The first CpG gets both rs1001 and rs1002. The second gets rs2001. Other CpGs
are fitted without SNP terms. A missing map entry is different from a mapped
SNP that is unavailable in genotype data.

```r
methylation <- readRDS("data/methylation.rds")
phenotypes <- readRDS("data/phenotypes.rds")
genotypes <- readRDS("data/genotypes.rds")
cpg_snp_map <- readRDS("data/cpg_snp_map.rds")

# Set categories and their reference levels intentionally.
phenotypes$sex <- factor(phenotypes$sex, levels = c("female", "male"))

fit <- ewas_mqtl(
  methylation = methylation,
  phenotypes = phenotypes,
  formula = ~ smoking + age + sex + cell_CD4 + cell_CD8 + PC1 + PC2,
  exposure = "smoking",
  genotypes = genotypes,
  cpg_snp_map = cpg_snp_map,
  methylation_scale = "M",
  sample_id = "sample_id"
)
```

Replace the example formula with your study's covariates. This package does not
choose confounders, cell estimates, ancestry PCs, technical variables, or a
population-specific mQTL reference. Including every cell proportion together
with an intercept can make the design singular when proportions sum to one.

Keep identifiers as character strings when importing files, particularly when
sample IDs have leading zeros. Do not rely on file order. SNP IDs can include
punctuation such as `1:12345:A:G`; matching is exact and does not require rsIDs.
Prepare consistent variant definitions, genome builds, and dosage coding upstream.

### Methylation tables and bigsnpr

To convert the original script's sample-by-CpG table:

```r
methylation <- methylation_from_table(
  wide_methylation,
  sample_id = "Sample_Name",
  cpg_cols = c("cg00000001", "cg00000002")
)
```

If `cpg_cols = NULL`, every column except the ID is treated as methylation.
Specify CpG columns explicitly when the table also contains numeric metadata.
The conversion transposes and copies the selected values, so it needs memory.

An optional adapter supports an already attached `bigsnpr` object:

```r
big <- bigsnpr::snp_attach("data/genotypes.rds")
snps <- unique(cpg_snp_map$snp_id[
  cpg_snp_map$cpg_id %in% rownames(methylation)
])
genotypes <- genotypes_from_bigsnp(
  big, snps = snps, sample_ids = phenotypes$sample_id,
  sample_col = "sample.ID", snp_col = "marker.ID"
)
```

This reads only the selected sample/SNP subset into memory and requires every
requested ID to exist. It does not impute missing values or perform genotype QC.
If no selected CpGs have mapped SNPs, omit the genotype and mapping inputs.

### M values and beta values

Set `methylation_scale = "M"` or `"beta"` to describe the supplied values.
Nothing is transformed automatically. Beta values must lie in [0,1]. The
coefficient is on the supplied methylation scale. M values are not beta-value
percentage changes. A coefficient of 0.02 on the beta scale means a difference
of 0.02, or two percentage points, per exposure unit or category contrast.

Preprocessing, normalization, probe filtering, methylation transformation,
genotype QC, and sample QC belong upstream of this regression workflow.

## Choose the analysis behavior

| Option | Default | Behavior |
|---|---|---|
| `sample_match` | `"strict"` | All supplied sample-ID sets must agree; inputs are reordered into phenotype order. `"intersect"` explicitly drops unmatched IDs and records them. |
| `missing_snps` | `"error"` | A selected CpG's unavailable mapped SNP stops the run. `"skip_cpg"` leaves an untested row. `"drop"` warns and reports reduced adjustment. |
| `cpgs` | `NULL` | Analyze every methylation row, including unmapped CpGs. Supply an explicit vector to restrict the analysis. |
| `min_samples` | `20` | Minimum complete samples at a CpG. This is a configurable safeguard, not a study-design or power rule. |
| `min_residual_df` | `3` | Minimum residual degrees of freedom. |
| `p_adjust_method` | `"BH"` | Adjust primary exposure p-values over all requested CpGs; Bonferroni values are also returned. |
| `keep_coefficients` | `"exposure"` | Retain exposure contrast rows. `"all"` additionally retains covariate/SNP/intercept rows; `"none"` keeps only CpG summaries. |
| `compare_unadjusted` | `FALSE` | Fit a model without SNP terms on exactly the same complete samples as each adjusted model. |
| `n_cores` | `1` | Sequential, portable execution. More workers use Unix forking; use a non-GUI Rscript session. Windows requires 1. |
| `chunk_size` | `1000` | Collect output tables in chunks; does not load input matrices from disk in chunks. |

Missing values are handled independently for each CpG. A sample must have
observed methylation, every phenotype design value, and all SNP dosages used
for that CpG. Missingness in a SNP does not remove a sample from unrelated CpGs.
There is no automatic imputation or global removal of samples with any missing
methylation value.

If a selected SNP is constant or exactly dependent on another predictor, the
affected model is marked `rank_deficient`. Review the predefined set or input
coding and rerun deliberately. The software does not automatically prune LD
or select a representative SNP. Near-collinearity can still yield imprecise
estimates even when the model is full rank.

## Read the results

`fit$results` has one row per requested CpG, in requested order:

| Columns | Interpretation |
|---|---|
| `cpg_id`, `status`, `message` | CpG identifier and whether inference succeeded. Start here before interpreting p-values. |
| `estimate`, `std_error`, `conf_low`, `conf_high` | Conditional exposure coefficient and its interval for a one-column exposure. `NA` for a multicolumn term. |
| `f_statistic`, `exposure_df`, `p_value` | Joint F test of all columns belonging to the selected exposure term. With one column, this is equivalent to a two-sided t test. |
| `p_adjust`, `p_bonferroni` | Corrected primary p-values across all requested CpGs. Failed models remain `NA` and count toward the testing-family size. |
| `n_aligned`, `n`, `n_excluded` | Aligned samples, complete samples for this model, and their difference. `n` is unknown when a model is skipped before its mask is formed. |
| `n_parameters`, `rank`, `df_residual` | Attempted design width, QR rank, and residual degrees of freedom. |
| `n_snps_requested`, `n_snps_available`, `n_snps_in_model` | Mapping size, available genotype columns, and columns included in the attempted design. |
| `snps_requested`, `snps_in_model`, `snps_missing` | Semicolon-separated SNP IDs. The structured CpG–SNP mapping remains in `snp_audit`. |
| `adjustment` | `full`, `partial`, `unavailable`, or `unmapped`; describes map coverage, independently of model success. |

Failure statuses include `missing_snps`, `invalid_methylation`,
`insufficient_samples`, `insufficient_df`, `rank_deficient`,
`constant_methylation`, and `degenerate_fit`. These rows are retained with
missing inferential results rather than silently disappearing.

`fit$coefficients` contains the requested coefficient rows, including the actual
factor contrast names, signed t statistics, confidence intervals, and raw
coefficient p-values. Its p-values are not corrected for multiplicity. For a
factor with more than two levels, the primary CpG table reports a single joint
exposure test; this coefficient table provides the individual contrasts.

`fit$sample_audit` records input membership, retention, and complete phenotype
design data. `fit$snp_audit` records the deduplicated map and selected/available
flags. Metadata stores the formula, model settings, factor levels and contrasts,
retained sample IDs, timestamps, and R session information. It does not store
raw input matrices. Per-CpG excluded sample IDs are not stored; sample counts are.

If requested, `fit$unadjusted` uses the same rows as the adjusted design, so an
effect change is not caused by comparing different sample subsets. Here
"unadjusted" means **without the mapped SNPs**; phenotype covariates remain.
Each table has its own corrected primary p-values. When a missing-SNP policy
skips a CpG, no paired comparison is fitted for that CpG. A paired baseline can
still succeed when the adjusted design is rank deficient.

## Statistical scope

The fitted exposure coefficient is conditional on the included covariates and
SNPs. SNP adjustment may change the scientific estimand and does not establish
causality, remove every genetic effect, or replace population-structure control.
The choice of which SNPs to condition on remains a study-design decision.

This implementation uses ordinary least squares with conventional standard
errors. Exact t/F inference assumes independent normal errors with common
variance. Related subjects, repeated measurements, or substantial residual
heteroskedasticity may require a different model or variance estimator. There
are no mixed effects, cluster-robust errors, or limma moderation in this version.

Use a one-sided formula and name an exact exposure term. For example,
`~ log(exposure) + age` requires `exposure = "log(exposure)"`. With interaction
terms, a selected main effect is conditional on those interactions, not an
omnibus test of all terms containing that variable. Intercepts are required;
offsets and `~ .` are rejected explicitly. Define transformations before calling
the function if they require custom missing-value handling.

BH correction is the default; `"BY"` is available for a more conservative
dependence allowance. Correction applies to the requested analysis set. If you
analyze separate batches in separate calls, combine the raw CpG p-values and
recompute correction over the complete intended family.

## Command-line use

From the extracted bundle:

```bash
Rscript epicurious/run_ewas.R --demo results/demo
Rscript epicurious/run_ewas.R epicurious/config.R
```

Copy and edit `config.R` for a real run. Input paths are RDS files and resolve
relative to the configuration file. The output path follows the same rule.
Genotype and map paths can be `NULL` for an analysis without SNP adjustment.
The CLI is separate from the modelling functions; it never installs packages
or modifies input data. Configuration files are R code, so use your own or a
trusted configuration.

The export contains `analysis.rds`, `results.csv`, `coefficients.csv`,
`sample_audit.csv`, `snp_audit.csv`, `metadata.txt`, `sessionInfo.txt`, and
optionally `unadjusted_same_samples.csv`.

The core uses in-memory dense matrices. A double methylation matrix alone
requires approximately `8 × n_CpGs × n_samples` bytes, before working memory
and output tables. For example, 850,000 CpGs × 1,000 samples is about 6.8 GB
(6.3 GiB). Start with a small CpG subset and one worker to assess memory and
runtime. No genome-wide performance claim is made for this first version.

## Development and validation

The source package separates validation, fitting, orchestration, import/export,
and simulation under `R/`. It includes installed R help pages and a base-R
regression suite comparing inference with `lm()` and nested-model `anova()`.

```bash
R CMD build epicurious
R CMD check --no-manual epicurious_0.1.0.tar.gz
```

See `VALIDATION.md` in the bundle for the checks actually performed and their
limits. Public-release maintainer metadata and the project license are left
for the project owner to set; the current package can be installed locally.

Reference documentation used for the implementation:
[R linear model fitting](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/lmfit.html),
[design matrices](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/model.matrix.html),
[multiple testing](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html),
and [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html).
