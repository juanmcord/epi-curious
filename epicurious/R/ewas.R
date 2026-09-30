#' EWAS with predefined CpG-specific additive SNP adjustment
#'
#' See ?ewas_mqtl for input contracts, inference, and failure policies.
#' Methylation is always the response; formula describes phenotype predictors.
#' @export
ewas_mqtl <- function(methylation, phenotypes, formula, exposure,
                      genotypes = NULL, cpg_snp_map = NULL,
                      sample_id = "sample_id", cpgs = NULL,
                      methylation_scale = c("M", "beta"),
                      sample_match = c("strict", "intersect"),
                      missing_snps = c("error", "skip_cpg", "drop"),
                      min_samples = 20L, min_residual_df = 3L,
                      conf_level = 0.95, p_adjust_method = "BH",
                      keep_coefficients = c("exposure", "all", "none"),
                      compare_unadjusted = FALSE, n_cores = 1L,
                      chunk_size = 1000L, verbose = TRUE) {
  started <- Sys.time()
  methylation_scale <- match.arg(methylation_scale)
  sample_match <- match.arg(sample_match)
  missing_snps <- match.arg(missing_snps)
  keep_coefficients <- match.arg(keep_coefficients)
  p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)
  .check_scalar(min_samples, "min_samples", 2, TRUE)
  .check_scalar(min_residual_df, "min_residual_df", 1, TRUE)
  .check_scalar(n_cores, "n_cores", 1, TRUE)
  .check_scalar(chunk_size, "chunk_size", 1, TRUE)
  .check_scalar(conf_level, "conf_level", 0)
  if (conf_level <= 0 || conf_level >= 1) stop("conf_level must be between 0 and 1.", call. = FALSE)
  .check_flag(compare_unadjusted, "compare_unadjusted")
  .check_flag(verbose, "verbose")
  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    stop("Use n_cores = 1 on Windows; parallel execution uses Unix forking.", call. = FALSE)
  }
  prepared <- .prepare_inputs(methylation, phenotypes, genotypes, cpg_snp_map,
                              sample_id, cpgs, sample_match, missing_snps)
  design <- .prepare_design(formula, exposure, prepared$phenotypes)
  x_base <- design$matrix
  if (any(startsWith(colnames(x_base), "SNP:"))) {
    stop("The design-column prefix 'SNP:' is reserved for genotypes.", call. = FALSE)
  }
  n_total <- length(prepared$sample_ids)
  if (sum(design$complete) < min_samples) {
    stop("Too few samples have complete phenotype/covariate data.", call. = FALSE)
  }
  complete_lookup <- match(prepared$sample_audit$sample_id, prepared$sample_ids)
  prepared$sample_audit$complete_covariates <- design$complete[complete_lookup]
  n_cpgs <- length(prepared$cpgs)
  if (verbose) message("Aligned ", n_total, " samples; fitting ", n_cpgs, " CpGs.")

  fit_one <- function(i) {
    cpg <- prepared$cpgs[i]
    snps <- prepared$snps_by_cpg[[cpg]]
    available <- snps[snps %in% colnames(genotypes)]
    absent <- setdiff(snps, available)
    row <- .empty_result(cpg, n_total, snps, available, absent, length(design$tested))
    if (length(absent) && missing_snps == "skip_cpg") {
      row$status <- "missing_snps"
      row$message <- "At least one requested SNP is unavailable; CpG skipped."
      return(list(result = row, coefficients = NULL, unadjusted = NULL))
    }
    y <- as.numeric(methylation[prepared$cpg_index[i], prepared$methylation_index])
    if (any(is.infinite(y)) ||
        (methylation_scale == "beta" && any(y < 0 | y > 1, na.rm = TRUE))) {
      row$status <- "invalid_methylation"
      row$message <- "Infinite methylation, or beta-scale values outside [0, 1]."
      return(list(result = row, coefficients = NULL, unadjusted = NULL))
    }
    x <- x_base
    if (length(available)) {
      g <- genotypes[prepared$genotype_index, available, drop = FALSE]
      colnames(g) <- paste0("SNP:", available)
      x <- cbind(x_base, g)
    }
    # One shared mask is applied to response, predictors, and SNPs together.
    complete <- is.finite(y) & stats::complete.cases(x)
    y <- y[complete]
    x <- x[complete, , drop = FALSE]
    row$n <- length(y)
    row$n_excluded <- n_total - row$n
    row$n_parameters <- ncol(x)
    row$n_snps_in_model <- length(available)
    row$snps_in_model <- paste(available, collapse = ";")
    fitted <- .fit_ols(x, y, design$tested, conf_level, min_samples, min_residual_df)
    fields <- intersect(names(row), names(fitted))
    row[fields] <- fitted[fields]
    coef <- NULL
    if (fitted$status == "ok" && keep_coefficients != "none") {
      coef <- fitted$coefficients
      coef$role <- ifelse(seq_len(nrow(coef)) %in% design$tested, "exposure",
                          ifelse(seq_len(nrow(coef)) > ncol(x_base), "snp", "covariate"))
      coef$role[coef$coefficient == "(Intercept)"] <- "intercept"
      if (keep_coefficients == "exposure") coef <- coef[coef$role == "exposure", , drop = FALSE]
      coef$cpg_id <- cpg
      coef$n <- row$n
      coef$df_residual <- row$df_residual
      coef <- coef[, names(.empty_coefficients()), drop = FALSE]
    }
    baseline <- NULL
    if (compare_unadjusted) {
      baseline <- row
      baseline$n_parameters <- ncol(x_base)
      baseline$n_snps_in_model <- 0L
      baseline$snps_in_model <- ""
      baseline$adjustment <- "unadjusted_same_samples"
      base_fit <- .fit_ols(x_base[complete, , drop = FALSE], y, design$tested,
                           conf_level, min_samples, min_residual_df)
      base_fields <- intersect(names(baseline), names(base_fit))
      baseline[base_fields] <- base_fit[base_fields]
    }
    list(result = row, coefficients = coef, unadjusted = baseline)
  }

  # Collect only compact tables per chunk, rather than every model object.
  chunks <- split(seq_len(n_cpgs), ceiling(seq_len(n_cpgs) / chunk_size))
  result_chunks <- vector("list", length(chunks))
  coef_chunks <- vector("list", length(chunks))
  baseline_chunks <- vector("list", length(chunks))
  for (j in seq_along(chunks)) {
    fits <- if (n_cores == 1L) {
      lapply(chunks[[j]], fit_one)
    } else {
      parallel::mclapply(chunks[[j]], fit_one, mc.cores = min(n_cores, length(chunks[[j]])),
                         mc.preschedule = TRUE, mc.set.seed = FALSE)
    }
    if (any(vapply(fits, function(z) inherits(z, "try-error") || is.null(z), logical(1)))) {
      stop("A worker failed. Re-run with n_cores = 1 to diagnose.", call. = FALSE)
    }
    result_chunks[[j]] <- do.call(rbind, lapply(fits, `[[`, "result"))
    coef_chunks[j] <- list(do.call(rbind, lapply(fits, `[[`, "coefficients")))
    if (compare_unadjusted) {
      baseline_chunks[[j]] <- do.call(rbind, lapply(fits, function(z) {
        if (!is.null(z$unadjusted)) return(z$unadjusted)
        r <- z$result
        r$adjustment <- "unadjusted_same_samples"
        r$message <- paste("Paired comparison unavailable:", r$message)
        r
      }))
    }
    if (verbose) message("Completed ", max(chunks[[j]]), "/", n_cpgs, " CpGs.")
  }
  results <- .adjust_results(do.call(rbind, result_chunks), p_adjust_method)
  coefficients <- do.call(rbind, coef_chunks)
  if (is.null(coefficients)) coefficients <- .empty_coefficients()
  rownames(results) <- rownames(coefficients) <- NULL
  unadjusted <- if (compare_unadjusted) {
    .adjust_results(do.call(rbind, baseline_chunks), p_adjust_method)
  } else NULL
  if (!is.null(unadjusted)) rownames(unadjusted) <- NULL
  if (verbose && any(results$status != "ok")) {
    message(sum(results$status != "ok"), " CpGs have no valid adjusted test; inspect results$status and results$message.")
  }
  output <- list(
    results = results, coefficients = coefficients, unadjusted = unadjusted,
    sample_audit = prepared$sample_audit, snp_audit = prepared$snp_audit,
    metadata = list(version = "0.1.0", function_name = "ewas_mqtl",
                    formula = paste(deparse(formula), collapse = " "), exposure = exposure,
                    exposure_columns = colnames(x_base)[design$tested],
                    design_columns = colnames(x_base), contrasts = design$contrasts,
                    factor_levels = design$factor_levels, methylation_scale = methylation_scale,
                    sample_ids = prepared$sample_ids, sample_match = sample_match,
                    missing_snps = missing_snps, min_samples = min_samples,
                    min_residual_df = min_residual_df, conf_level = conf_level,
                    p_adjust_method = p_adjust_method, testing_family_size = n_cpgs,
                    duplicate_pairs_removed = prepared$duplicate_pairs_removed,
                    keep_coefficients = keep_coefficients,
                    compare_unadjusted = compare_unadjusted, n_cores = n_cores,
                    chunk_size = chunk_size, started = started, finished = Sys.time(),
                    session_info = utils::sessionInfo())
  )
  class(output) <- "epicurious_ewas"
  output
}

#' @export
print.epicurious_ewas <- function(x, ...) {
  cat("epi-curious EWAS\n")
  cat("  Formula: methylation", x$metadata$formula, "+ CpG-specific SNPs\n")
  cat("  Exposure:", x$metadata$exposure, "\n")
  cat("  Aligned samples:", length(x$metadata$sample_ids), "\n")
  cat("  CpGs:", nrow(x$results), "| Valid tests:", sum(x$results$status == "ok"), "\n")
  cat("  Adjustment:", x$metadata$p_adjust_method, "over", x$metadata$testing_family_size, "requested CpGs\n")
  print(table(x$results$status))
  invisible(x)
}
