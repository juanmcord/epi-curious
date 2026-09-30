#' Generate a small synthetic dataset with known exposure and SNP effects
#' @export
simulate_epicurious <- function(n = 200L, n_cpgs = 12L, seed = 42L) {
  .check_scalar(n, "n", 30, TRUE)
  .check_scalar(n_cpgs, "n_cpgs", 4, TRUE)
  .check_scalar(seed, "seed", 0, TRUE)
  # Restore the caller's random-number state even when none existed before.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  ids <- sprintf("sample_%03d", seq_len(n))
  cpgs <- sprintf("cg_demo_%04d", seq_len(n_cpgs))
  g <- matrix(stats::rbinom(n * 4, 2, 0.3), nrow = n,
              dimnames = list(ids, paste0("rs_demo_", 1:4)))
  phenotype <- data.frame(sample_id = ids,
                          exposure = 0.7 * as.numeric(scale(g[, 1])) + stats::rnorm(n),
                          age = stats::runif(n, 30, 75),
                          sex = factor(rep(c("female", "male"), length.out = n)))
  m <- matrix(stats::rnorm(n * n_cpgs, sd = 0.6), nrow = n_cpgs,
              dimnames = list(cpgs, ids))
  m <- sweep(m, 2, 0.012 * (phenotype$age - 50), "+")
  m[1, ] <- m[1, ] + 0.5 * phenotype$exposure + 1.2 * g[, 1]
  m[2, ] <- m[2, ] + 0.3 * phenotype$exposure + 0.9 * g[, 2] - 0.6 * g[, 3]
  m[3, ] <- m[3, ] + 1.5 * g[, 1]  # SNP-associated, true exposure effect zero.
  m[4, ] <- m[4, ] - 0.4 * phenotype$exposure
  map <- data.frame(cpg_id = c(cpgs[1], cpgs[2], cpgs[2], cpgs[3]),
                     snp_id = c("rs_demo_1", "rs_demo_2", "rs_demo_3", "rs_demo_1"))
  truth <- data.frame(cpg_id = cpgs, exposure_effect = c(0.5, 0.3, 0, -0.4, rep(0, n_cpgs - 4)))
  list(methylation = m, phenotypes = phenotype, genotypes = g,
       cpg_snp_map = map, truth = truth, methylation_scale = "M")
}
