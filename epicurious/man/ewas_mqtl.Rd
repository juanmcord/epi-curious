\name{ewas_mqtl}
\alias{ewas_mqtl}
\title{EWAS with CpG-specific predefined SNP adjustment}
\description{Fits a separate ordinary least squares methylation model for each
CpG, with a common phenotype design and that CpG's mapped additive SNP dosages.}
\usage{
ewas_mqtl(methylation, phenotypes, formula, exposure,
  genotypes = NULL, cpg_snp_map = NULL, sample_id = "sample_id",
  cpgs = NULL, methylation_scale = c("M", "beta"),
  sample_match = c("strict", "intersect"),
  missing_snps = c("error", "skip_cpg", "drop"),
  min_samples = 20L, min_residual_df = 3L, conf_level = 0.95,
  p_adjust_method = "BH", keep_coefficients = c("exposure", "all", "none"),
  compare_unadjusted = FALSE, n_cores = 1L, chunk_size = 1000L,
  verbose = TRUE)
}
\arguments{
\item{methylation}{Numeric matrix: CpGs in rows, sample IDs in column names.
Unique row and column names are required. Missing values are allowed.}
\item{phenotypes}{Data frame with one row per sample, a unique ID column, and all formula variables.}
\item{formula}{One-sided formula, e.g. \code{~ exposure + age + sex}.
An intercept is required. Offsets and dot expansion are not supported.
All variables must be present in phenotypes.}
\item{exposure}{Exact term label from \code{terms(formula)}, e.g. \code{"exposure"},
\code{"group"}, or \code{"log(exposure)"}. All design columns assigned to this
one term are tested jointly, conditional on all other included columns.}
\item{genotypes}{Optional numeric matrix with samples in rows and SNPs in columns.
Unique row and column names are required. Values are allele dosages in [0,2]
or missing. Needed when selected CpGs have mapped SNPs.}
\item{cpg_snp_map}{Optional data frame with \code{cpg_id} and \code{snp_id}.
One row per CpG-SNP pair; multiple SNPs per CpG are allowed. Duplicate pairs
are removed and counted. Extra columns are ignored. No associations are selected
or discovered by this function.}
\item{sample_id}{Name of the phenotype sample ID column.}
\item{cpgs}{Unique character vector of CpGs to analyze, in the desired order.
NULL analyzes all methylation rows. Unknown CpGs are an error.}
\item{methylation_scale}{Declared response scale: M values (default) or beta
values in [0,1]. Data are never automatically transformed.}
\item{sample_match}{Strict requires equal sample ID sets and reorders by
phenotype order. Intersect explicitly retains common IDs and records exclusions.
If genotypes are supplied, their sample set participates in alignment.}
\item{missing_snps}{Error (default) stops for unavailable mapped SNPs.
Skip-cpg returns an untested row for affected CpGs. Drop warns and uses only
available mapped SNPs, marking partial or unavailable adjustment.}
\item{min_samples}{Minimum complete observations per CpG; default 20 is a
configurable computational safeguard, not a power recommendation.}
\item{min_residual_df}{Minimum residual degrees of freedom; default 3.}
\item{conf_level}{Coefficient confidence level strictly between 0 and 1.}
\item{p_adjust_method}{Method accepted by \code{stats::p.adjust}; default BH.}
\item{keep_coefficients}{Retain exposure coefficients (default), all model
coefficients, or none. Coefficient-level p-values are unadjusted.}
\item{compare_unadjusted}{Also fit the common phenotype model without SNPs
using exactly the adjusted model's complete samples for each CpG.}
\item{n_cores}{Workers; default 1 works on all platforms. Values greater than 1
use Unix forking, intended for non-GUI Rscript sessions, and are rejected on Windows.}
\item{chunk_size}{Number of CpGs collected per chunk. Input matrices stay in memory;
this is not an out-of-core implementation.}
\item{verbose}{Print alignment and chunk progress.}
}
\details{
For CpG j the response is methylation and the predictors are the supplied
phenotype design plus only its mapped SNPs. SNP effects are additive dosages.
A single complete-case mask covers methylation, the phenotype design, and the
available mapped SNPs. No imputation is performed. Infinite design or genotype
values are errors; invalid methylation rows are flagged individually.

The primary p-value tests all columns of the selected exposure term jointly
with an F test. For a one-column exposure this equals the two-sided t-test
p-value, and the main table also contains its effect estimate and interval.
For a multicolumn exposure those scalar effect fields are NA; see the coefficient
table for contrasts. Factor reference levels and contrasts are retained in metadata.
When interactions are included, the selected term is tested conditionally on
the interaction columns, not as an overall test of every occurrence of the variable.

Standard errors assume independent, homoskedastic errors; exact t/F inference
also assumes normal errors. There is no mixed-model, robust-standard-error,
limma moderation, or relatedness correction. Rank-deficient models are skipped,
including constant SNPs and exact linear dependence; no SNP is silently selected
to resolve linkage disequilibrium. Near-collinearity may still inflate uncertainty.

Multiplicity correction uses one primary p-value per requested CpG, with
family size equal to the number requested, including failed models.
BH and Bonferroni-adjusted values are returned (or the chosen method plus
Bonferroni). Failure p-values remain NA. A subset run corrects only that subset;
recompute globally after combining separately run batches.
}
\value{An \code{epicurious_ewas} list with \code{results} (one row per requested
CpG), \code{coefficients}, optional \code{unadjusted} (same-sample comparisons),
\code{sample_audit}, \code{snp_audit}, and \code{metadata}. Model failures retain
their row, status, message, and missing p-value. Results include sample counts,
requested/available/in-model SNP counts and IDs, rank, degrees of freedom,
conditional effect summaries, primary F-test p-values, and adjusted p-values.}
\seealso{\code{\link{simulate_epicurious}}, \code{\link{write_ewas}}}
\examples{
d <- simulate_epicurious(n = 80, n_cpgs = 4)
fit <- ewas_mqtl(d$methylation, d$phenotypes,
  ~ exposure + age + sex, exposure = "exposure",
  genotypes = d$genotypes, cpg_snp_map = d$cpg_snp_map,
  methylation_scale = "M", verbose = FALSE)
fit$results[, c("cpg_id", "estimate", "p_value", "p_adjust", "status")]
}
