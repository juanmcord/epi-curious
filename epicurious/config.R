# Copy this file and edit the paths, formula, and exposure term.
# Paths below are resolved relative to this configuration file.
config <- list(
  methylation_file = "data/methylation.rds",  # numeric matrix: CpGs x samples
  phenotypes_file = "data/phenotypes.rds",    # data.frame with sample_id
  genotypes_file = "data/genotypes.rds",      # numeric matrix: samples x SNPs
  cpg_snp_map_file = "data/cpg_snp_map.rds",   # data.frame: cpg_id, snp_id
  output_dir = "results/my_ewas",
  formula = ~ exposure + age + sex,
  exposure = "exposure",
  sample_id = "sample_id",
  methylation_scale = "M",                    # or "beta"; no conversion occurs
  sample_match = "strict",                   # explicitly use "intersect" if needed
  missing_snps = "error",
  cpgs = NULL,                                # NULL means all methylation rows
  min_samples = 20L,
  min_residual_df = 3L,
  p_adjust_method = "BH",
  keep_coefficients = "exposure",
  compare_unadjusted = TRUE,
  n_cores = 1L,
  chunk_size = 1000L,
  verbose = TRUE
)
