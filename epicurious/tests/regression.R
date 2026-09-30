# Dependency-free statistical regression tests, run by R CMD check.
library(epicurious)
checks <- 0L
check <- function(label, code) {
  force(code)
  checks <<- checks + 1L
  cat(sprintf("PASS %02d: %s\n", checks, label))
}
expect_error <- function(code, pattern) {
  e <- tryCatch({ force(code); NULL }, error = identity)
  if (!inherits(e, "error") || !grepl(pattern, conditionMessage(e), fixed = TRUE)) {
    stop("Expected an error containing: ", pattern)
  }
}
equal <- function(a, b, tolerance = 1e-8) {
  if (!isTRUE(all.equal(unname(a), unname(b), tolerance = tolerance, check.attributes = FALSE))) {
    stop("Values differ: ", paste(a, collapse = ", "), " versus ", paste(b, collapse = ", "))
  }
}
d <- simulate_epicurious(n = 160, n_cpgs = 8)
run <- function(data = d, ...) {
  do.call(ewas_mqtl, c(list(methylation = data$methylation,
                            phenotypes = data$phenotypes,
                            formula = ~ exposure + age + sex,
                            exposure = "exposure", genotypes = data$genotypes,
                            cpg_snp_map = data$cpg_snp_map,
                            methylation_scale = "M", verbose = FALSE), list(...)))
}
fit <- run(keep_coefficients = "all", compare_unadjusted = TRUE)
check("CpG-specific models match lm estimates, SEs, confidence intervals, and p-values", {
  for (i in seq_len(nrow(d$methylation))) {
    cpg <- rownames(d$methylation)[i]
    snps <- unique(d$cpg_snp_map$snp_id[d$cpg_snp_map$cpg_id == cpg])
    dat <- d$phenotypes
    dat$y <- as.numeric(d$methylation[i, ])
    for (snp in snps) dat[[snp]] <- d$genotypes[, snp]
    f <- reformulate(c("exposure", "age", "sex", snps), response = "y")
    oracle <- lm(f, data = dat)
    tab <- coef(summary(oracle))
    row <- fit$results[i, ]
    equal(row$estimate, tab["exposure", 1])
    equal(row$std_error, tab["exposure", 2])
    equal(row$f_statistic, tab["exposure", 3]^2)
    equal(row$p_value, tab["exposure", 4])
    equal(c(row$conf_low, row$conf_high), confint(oracle)["exposure", ])
    stored <- fit$coefficients[fit$coefficients$cpg_id == cpg, ]
    equal(stored$estimate, tab[, 1])
    equal(stored$std_error, tab[, 2])
  }
})
check("Different input sample orders produce identical results", {
  z <- d
  z$methylation <- z$methylation[, rev(seq_len(ncol(z$methylation)))]
  z$genotypes <- z$genotypes[c(2:nrow(z$genotypes), 1), ]
  z$phenotypes <- z$phenotypes[rev(seq_len(nrow(z$phenotypes))), ]
  equal(run(z)$results$estimate, fit$results$estimate)
  equal(run(z)$results$p_value, fit$results$p_value)
})
check("Joint phenotype, methylation, and genotype missingness matches lm complete cases", {
  z <- d
  z$phenotypes$age[1:3] <- NA
  z$methylation[2, 4:6] <- NA
  z$genotypes[7:10, "rs_demo_2"] <- NA
  got <- run(z, compare_unadjusted = TRUE)
  dat <- z$phenotypes
  dat$y <- z$methylation[2, ]
  dat$s1 <- z$genotypes[, "rs_demo_2"]
  dat$s2 <- z$genotypes[, "rs_demo_3"]
  ref <- lm(y ~ exposure + age + sex + s1 + s2, dat)
  equal(got$results$estimate[2], coef(ref)["exposure"])
  equal(got$results$p_value[2], coef(summary(ref))["exposure", 4])
  stopifnot(got$results$n[2] == nobs(ref), got$results$n[2] == 150)
  cc <- complete.cases(dat)
  base <- lm(y ~ exposure + age + sex, dat[cc, ])
  equal(got$unadjusted$estimate[2], coef(base)["exposure"])
  stopifnot(identical(got$results$n, got$unadjusted$n))
  # The SNP missing only for CpG 2 must not remove people from unmapped CpG 4.
  stopifnot(got$results$n[4] == 157)
})
check("Multilevel exposure gets a joint F test matching nested-model anova", {
  z <- d
  z$phenotypes$group <- factor(rep(c("control", "low", "high"), length.out = 160),
                               levels = c("control", "low", "high"))
  g <- ewas_mqtl(z$methylation, z$phenotypes, ~ group + age + sex, "group",
                 z$genotypes, z$cpg_snp_map, verbose = FALSE)
  dat <- z$phenotypes
  dat$y <- z$methylation[2, ]
  dat$s1 <- z$genotypes[, 2]
  dat$s2 <- z$genotypes[, 3]
  small <- lm(y ~ age + sex + s1 + s2, dat)
  full <- lm(y ~ group + age + sex + s1 + s2, dat)
  a <- anova(small, full)
  equal(g$results$p_value[2], a$`Pr(>F)`[2])
  equal(g$results$f_statistic[2], a$F[2])
  stopifnot(g$results$exposure_df[2] == 2, is.na(g$results$estimate[2]))
  equal(g$coefficients$estimate[g$coefficients$cpg_id == rownames(z$methylation)[2]],
        coef(full)[c("grouplow", "grouphigh")])
})
check("Strict matching rejects unequal sample sets", {
  z <- d
  z$genotypes <- z$genotypes[-1, ]
  expect_error(run(z), "Sample ID sets differ")
})
check("Intersection is explicit, reorders correctly, and audits exclusions", {
  z <- d
  z$genotypes <- z$genotypes[-1, ]
  g <- run(z, sample_match = "intersect")
  stopifnot(sum(g$sample_audit$retained) == 159, !g$sample_audit$retained[1])
  stopifnot(all(g$results$n == 159))
})
check("Duplicate sample identifiers are rejected", {
  z <- d
  z$phenotypes$sample_id[2] <- z$phenotypes$sample_id[1]
  expect_error(run(z), "unique")
})
check("Duplicate CpG identifiers are rejected", {
  z <- d
  rownames(z$methylation)[2] <- rownames(z$methylation)[1]
  expect_error(run(z), "unique")
})
check("Unknown CpGs are rejected instead of silently omitted", {
  expect_error(run(cpgs = "unknown"), "CpGs absent")
})
check("A requested SNP absent from genotype data is an error by default", {
  z <- d
  z$genotypes <- z$genotypes[, -2]
  expect_error(run(z), "Mapped SNPs absent")
})
check("Skip policy retains a failed result row and all valid other CpGs", {
  z <- d
  z$genotypes <- z$genotypes[, -2]
  g <- run(z, missing_snps = "skip_cpg", compare_unadjusted = TRUE)
  stopifnot(nrow(g$results) == 8, g$results$status[2] == "missing_snps")
  stopifnot(is.na(g$results$p_value[2]), all(g$results$status[-2] == "ok"))
  stopifnot(g$unadjusted$status[2] == "missing_snps")
})
check("Opt-in SNP dropping is flagged as partial adjustment", {
  z <- d
  z$genotypes <- z$genotypes[, -2]
  g <- suppressWarnings(run(z, missing_snps = "drop"))
  stopifnot(g$results$adjustment[2] == "partial", g$results$n_snps_in_model[2] == 1)
  stopifnot(g$results$snps_missing[2] == "rs_demo_2")
})
check("Unselected map entries do not block a requested CpG subset", {
  z <- d
  z$cpg_snp_map <- rbind(z$cpg_snp_map, data.frame(cpg_id = "unused", snp_id = "missing"))
  g <- run(z, cpgs = rownames(z$methylation)[4])
  stopifnot(nrow(g$results) == 1, g$results$status == "ok")
})
check("Duplicate map pairs are removed and counted", {
  z <- d
  z$cpg_snp_map <- rbind(z$cpg_snp_map, z$cpg_snp_map[1, ])
  g <- run(z)
  stopifnot(g$metadata$duplicate_pairs_removed == 1, g$results$n_snps_requested[1] == 1)
  equal(g$results$p_value, fit$results$p_value)
})
check("Dosages outside [0, 2] are rejected", {
  z <- d
  z$genotypes[1, 1] <- 9
  expect_error(run(z), "[0, 2]")
})
check("Fractional dosage is accepted and matches lm", {
  z <- d
  z$genotypes[, 1] <- 0.9 * z$genotypes[, 1] + 0.05
  g <- run(z)
  dat <- z$phenotypes
  dat$y <- z$methylation[1, ]
  dat$snp <- z$genotypes[, 1]
  equal(g$results$estimate[1], coef(lm(y ~ exposure + age + sex + snp, dat))["exposure"])
})
check("Dependent SNPs fail explicitly instead of yielding arbitrary coefficients", {
  z <- d
  z$genotypes[, 3] <- z$genotypes[, 2]
  g <- run(z, compare_unadjusted = TRUE)
  stopifnot(g$results$status[2] == "rank_deficient", is.na(g$results$p_value[2]))
  stopifnot(g$unadjusted$status[2] == "ok")
})
check("Monomorphic SNPs are diagnosed", {
  z <- d
  z$genotypes[, 1] <- 0
  g <- run(z)
  stopifnot(all(g$results$status[c(1, 3)] == "rank_deficient"))
})
check("Constant methylation is diagnosed", {
  z <- d
  z$methylation[4, ] <- 0.5
  stopifnot(run(z)$results$status[4] == "constant_methylation")
})
check("Too few complete observations return a status without aborting other CpGs", {
  z <- d
  z$methylation[1, 1:150] <- NA
  g <- run(z)
  stopifnot(g$results$status[1] == "insufficient_samples", g$results$n[1] == 10)
})
check("Residual degree-of-freedom threshold is enforced", {
  g <- run(min_residual_df = 200)
  stopifnot(all(g$results$status == "insufficient_df"))
})
check("SNP IDs containing punctuation are supported without formula parsing", {
  z <- d
  colnames(z$genotypes)[1] <- "1:12345:A:G"
  z$cpg_snp_map$snp_id[z$cpg_snp_map$snp_id == "rs_demo_1"] <- "1:12345:A:G"
  equal(run(z)$results$p_value, fit$results$p_value)
})
check("An EWAS without genotype input matches the unadjusted paired comparison", {
  z <- d
  z$genotypes <- NULL
  z$cpg_snp_map <- NULL
  equal(run(z)$results$p_value, fit$unadjusted$p_value)
})
check("A single selected CpG and SNP keep matrix dimensions", {
  g <- run(cpgs = rownames(d$methylation)[1])
  stopifnot(nrow(g$results) == 1, g$results$n_snps_in_model == 1)
  equal(g$results$estimate, fit$results$estimate[1])
})
check("Multiple testing includes requested CpGs with failed models", {
  z <- d
  z$methylation[1, ] <- NA
  g <- run(z)
  equal(g$results$p_adjust, p.adjust(g$results$p_value, "BH", n = 8))
  equal(g$results$p_bonferroni, p.adjust(g$results$p_value, "bonferroni", n = 8))
})
check("All CpGs failing still produce valid empty coefficient and result tables", {
  z <- d
  z$methylation[,] <- NA
  g <- run(z)
  stopifnot(nrow(g$coefficients) == 0, nrow(g$results) == 8, all(is.na(g$results$p_adjust)))
})
check("Beta scale is checked and never automatically transformed", {
  g <- ewas_mqtl(d$methylation, d$phenotypes, ~ exposure + age + sex, "exposure",
                 methylation_scale = "beta", verbose = FALSE)
  stopifnot(all(g$results$status == "invalid_methylation"))
  b <- plogis(d$methylation)
  g <- ewas_mqtl(b, d$phenotypes, ~ exposure + age + sex, "exposure",
                 methylation_scale = "beta", verbose = FALSE)
  dat <- d$phenotypes
  dat$y <- b[1, ]
  equal(g$results$estimate[1], coef(lm(y ~ exposure + age + sex, dat))["exposure"])
})
check("Unsupported formulas and missing phenotype columns fail clearly", {
  expect_error(ewas_mqtl(d$methylation, d$phenotypes, y ~ exposure, "exposure"), "one-sided")
  expect_error(ewas_mqtl(d$methylation, d$phenotypes, ~ exposure + missing, "exposure"), "explicit formula variables")
  expect_error(ewas_mqtl(d$methylation, d$phenotypes, ~ exposure + offset(age), "exposure"), "offset")
  expect_error(ewas_mqtl(d$methylation, d$phenotypes, ~ exposure - 1, "exposure"), "intercept")
  expect_error(ewas_mqtl(d$methylation, d$phenotypes, ~ ., "exposure"), "explicit formula variables")
})
check("Wide table conversion preserves values and identifiers", {
  tab <- data.frame(sample_id = colnames(d$methylation), t(d$methylation), check.names = FALSE)
  stopifnot(identical(methylation_from_table(tab), d$methylation))
})
check("bigSNP adapter preserves requested sample/SNP order on a matrix-backed fixture", {
  fixture <- structure(list(fam = data.frame(sample.ID = rownames(d$genotypes)),
                             map = data.frame(marker.ID = colnames(d$genotypes)),
                             genotypes = d$genotypes), class = "bigSNP")
  snps <- colnames(d$genotypes)[c(3, 1)]
  samples <- rownames(d$genotypes)[c(9, 3, 1)]
  stopifnot(identical(genotypes_from_bigsnp(fixture, snps, samples), d$genotypes[samples, snps]))
})
check("Simulation is deterministic and does not change caller RNG state", {
  set.seed(2026)
  before <- .Random.seed
  a <- simulate_epicurious()
  stopifnot(identical(before, .Random.seed), identical(a, simulate_epicurious()))
})
check("Output export round-trips the object and prevents silent overwrite", {
  directory <- tempfile("epicurious-test-")
  paths <- write_ewas(fit, directory)
  stopifnot(all(file.exists(paths)), identical(readRDS(paths[["analysis.rds"]]), fit))
  expect_error(write_ewas(fit, directory), "already exist")
  write_ewas(run(), directory, overwrite = TRUE)
  stopifnot(!file.exists(file.path(directory, "unadjusted_same_samples.csv")))
  unlink(directory, recursive = TRUE)
})
check("Chunking and coefficient storage do not affect inference", {
  g <- run(chunk_size = 1, keep_coefficients = "none")
  equal(g$results$p_value, fit$results$p_value)
  stopifnot(nrow(g$coefficients) == 0)
})
if (.Platform$OS.type != "windows") {
  check("Two-worker results match sequential results", {
    g <- run(n_cores = 2, chunk_size = 3)
    equal(g$results$p_value, fit$results$p_value)
    equal(g$results$estimate, fit$results$estimate)
  })
}
cat(sprintf("\nAll %d regression checks passed.\n", checks))
