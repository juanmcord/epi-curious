# After installing epicurious, run this file from an R session.
library(epicurious)
d <- simulate_epicurious(n = 200, n_cpgs = 12)
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
print(fit$results[, c("cpg_id", "estimate", "p_value", "p_adjust", "status")])
# Use a new directory, or explicitly allow overwrite when appropriate.
# write_ewas(fit, "results/my_first_ewas")
