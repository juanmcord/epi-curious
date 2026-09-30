# epi-curious 0.1.0 -- standalone entry point
# Source this file, then call ewas_mqtl(). Only standard R packages are used.
# Full input guide and source package: epicurious_bundle.zip.
# This file is assembled from the package R/ sources; edit those for development.
# Example:
# d <- simulate_epicurious()
# fit <- ewas_mqtl(d$methylation, d$phenotypes, ~ exposure + age + sex,
#                  "exposure", d$genotypes, d$cpg_snp_map,
#                  methylation_scale = "M")

# Module: validation.R
# Internal checks deliberately avoid fixing or guessing identifiers.
.check_ids <- function(x, label) {
  if (is.null(x) || !length(x) || anyNA(x) ||
      any(!nzchar(as.character(x))) || anyDuplicated(x)) {
    stop(label, " must be present, non-empty, non-missing, and unique.",
         call. = FALSE)
  }
  if (any(trimws(as.character(x)) != as.character(x))) {
    stop(label, " contain leading or trailing whitespace.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_scalar <- function(x, label, lower = 1, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower ||
      (integer && x != floor(x))) {
    stop(label, " must be a single ", if (integer) "integer " else "number ",
         ">= ", lower, ".", call. = FALSE)
  }
}

.check_flag <- function(x, label) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(label, " must be TRUE or FALSE.", call. = FALSE)
  }
}

.check_matrix <- function(x, label) {
  if (!is.matrix(x) || !is.numeric(x) || any(dim(x) < 1L)) {
    stop(label, " must be a non-empty numeric matrix.", call. = FALSE)
  }
  .check_ids(rownames(x), paste(label, "row names"))
  .check_ids(colnames(x), paste(label, "column names"))
}

.prepare_inputs <- function(methylation, phenotypes, genotypes, cpg_snp_map,
                            sample_id, cpgs, sample_match, missing_snps) {
  .check_matrix(methylation, "methylation (CpGs x samples)")
  if (!is.data.frame(phenotypes)) {
    stop("phenotypes must be a data.frame.", call. = FALSE)
  }
  phenotypes <- as.data.frame(phenotypes)
  .check_ids(names(phenotypes), "phenotypes column names")
  if (!is.character(sample_id) || length(sample_id) != 1L ||
      is.na(sample_id) || !sample_id %in% names(phenotypes)) {
    stop("sample_id must name a column in phenotypes.", call. = FALSE)
  }
  pids <- as.character(phenotypes[[sample_id]])
  .check_ids(pids, "phenotypes sample IDs")
  mids <- colnames(methylation)
  if (!is.null(genotypes)) .check_matrix(genotypes, "genotypes (samples x SNPs)")
  gids <- if (is.null(genotypes)) NULL else rownames(genotypes)

  if (is.null(cpgs)) cpgs <- rownames(methylation)
  if (!is.character(cpgs)) stop("cpgs must be a character vector.", call. = FALSE)
  .check_ids(cpgs, "requested CpG IDs")
  unknown <- setdiff(cpgs, rownames(methylation))
  if (length(unknown)) {
    stop("CpGs absent from methylation: ", paste(utils::head(unknown, 10), collapse = ", "),
         call. = FALSE)
  }
  if (is.null(cpg_snp_map)) {
    cpg_snp_map <- data.frame(cpg_id = character(), snp_id = character())
  }
  if (!is.data.frame(cpg_snp_map) ||
      !all(c("cpg_id", "snp_id") %in% names(cpg_snp_map))) {
    stop("cpg_snp_map must contain cpg_id and snp_id columns.", call. = FALSE)
  }
  map <- as.data.frame(cpg_snp_map)[, c("cpg_id", "snp_id"), drop = FALSE]
  map[] <- lapply(map, as.character)
  for (nm in names(map)) {
    v <- map[[nm]]
    if (anyNA(v) || any(!nzchar(v)) || any(trimws(v) != v)) {
      stop("cpg_snp_map IDs must be non-missing and non-blank, without surrounding whitespace.",
           call. = FALSE)
    }
  }
  duplicates <- sum(duplicated(map))
  map <- unique(map)
  map$selected <- map$cpg_id %in% cpgs
  map$available <- map$snp_id %in% colnames(genotypes)
  selected_map <- map[map$selected, , drop = FALSE]
  if (nrow(selected_map) && is.null(genotypes)) {
    stop("genotypes are required when selected CpGs have mapped SNPs.", call. = FALSE)
  }
  absent <- unique(selected_map$snp_id[!selected_map$available])
  if (length(absent) && missing_snps == "error") {
    stop("Mapped SNPs absent from genotypes: ", paste(utils::head(absent, 10), collapse = ", "),
         ". Supply them, or explicitly choose missing_snps = 'skip_cpg' or 'drop'.",
         call. = FALSE)
  }
  if (length(absent) && missing_snps == "drop") {
    warning(length(absent), " mapped SNP(s) unavailable; affected models will have reduced adjustment. See snp_audit.",
            call. = FALSE)
  }
  requested_snps <- unique(selected_map$snp_id[selected_map$available])
  # Check one SNP at a time, avoiding a full-size logical copy of genotype data.
  for (snp in requested_snps) {
    v <- genotypes[, snp]
    if (any(is.infinite(v)) || any(v < 0 | v > 2, na.rm = TRUE)) {
      stop("Genotype dosages must be in [0, 2] or NA; invalid SNP: ", snp,
           call. = FALSE)
    }
  }

  strict_ok <- setequal(pids, mids) &&
    (is.null(gids) || setequal(pids, gids))
  if (sample_match == "strict" && !strict_ok) {
    stop("Sample ID sets differ. Use sample_match = 'intersect' only if intentional; exclusions are then audited.",
         call. = FALSE)
  }
  retained <- pids[pids %in% mids & (is.null(gids) | pids %in% gids)]
  if (!length(retained)) stop("No common samples across the inputs.", call. = FALSE)
  all_ids <- unique(c(pids, mids, gids))
  audit <- data.frame(sample_id = all_ids,
                      in_phenotypes = all_ids %in% pids,
                      in_methylation = all_ids %in% mids,
                      in_genotypes = if (is.null(gids)) NA else all_ids %in% gids,
                      retained = all_ids %in% retained,
                      stringsAsFactors = FALSE)
  pheno <- phenotypes[match(retained, pids), , drop = FALSE]
  rownames(pheno) <- retained
  list(phenotypes = pheno, sample_ids = retained,
       methylation_index = match(retained, mids),
       genotype_index = if (is.null(gids)) NULL else match(retained, gids),
       cpgs = cpgs, cpg_index = match(cpgs, rownames(methylation)),
       snps_by_cpg = split(selected_map$snp_id, selected_map$cpg_id),
       sample_audit = audit, snp_audit = map, duplicate_pairs_removed = duplicates)
}

.prepare_design <- function(formula, exposure, phenotypes) {
  if (!inherits(formula, "formula") || length(formula) != 2L) {
    stop("formula must be one-sided, e.g. ~ exposure + age + sex.", call. = FALSE)
  }
  vars <- all.vars(formula)
  if ("." %in% vars || !all(vars %in% names(phenotypes))) {
    stop("Use explicit formula variables, all present in phenotypes; '.' is not supported.",
         call. = FALSE)
  }
  trm <- stats::terms(formula)
  if (length(attr(trm, "offset")) || attr(trm, "intercept") != 1L) {
    stop("The formula must include an intercept and cannot contain offsets.", call. = FALSE)
  }
  labels <- attr(trm, "term.labels")
  if (!is.character(exposure) || length(exposure) != 1L || is.na(exposure) ||
      !exposure %in% labels) {
    stop("exposure must exactly match one formula term: ", paste(labels, collapse = ", "),
         call. = FALSE)
  }
  frame <- tryCatch(
    stats::model.frame(trm, data = phenotypes, na.action = stats::na.pass,
                       drop.unused.levels = TRUE),
    error = function(e) stop("Cannot construct model frame: ", conditionMessage(e), call. = FALSE)
  )
  design <- tryCatch(stats::model.matrix(trm, data = frame),
                     error = function(e) stop("Cannot construct design matrix: ", conditionMessage(e), call. = FALSE))
  if (nrow(design) != nrow(phenotypes) ||
      !identical(rownames(design), rownames(phenotypes))) {
    stop("Design-matrix rows do not match aligned samples.", call. = FALSE)
  }
  if (any(is.infinite(design))) {
    stop("The design matrix contains infinite values; check covariates and transformations.", call. = FALSE)
  }
  .check_ids(colnames(design), "design matrix column names")
  tested <- which(attr(design, "assign") == match(exposure, labels))
  if (!length(tested)) stop("The exposure produced no design columns.", call. = FALSE)
  list(matrix = design, tested = tested, complete = stats::complete.cases(design),
       contrasts = attr(design, "contrasts"),
       factor_levels = lapply(frame[vapply(frame, is.factor, logical(1))], levels))
}


# Module: fit.R
# Fit a full-rank OLS model. Inference uses the public lm.fit API and its QR.
.fit_ols <- function(x, y, tested, conf_level, min_samples, min_residual_df) {
  n <- length(y)
  p <- ncol(x)
  out <- list(status = "ok", message = "", rank = NA_integer_,
              df_residual = NA_integer_, estimate = NA_real_, std_error = NA_real_,
              conf_low = NA_real_, conf_high = NA_real_, f_statistic = NA_real_,
              p_value = NA_real_, coefficients = NULL)
  fail <- function(status, message) {
    out$status <- status
    out$message <- message
    out
  }
  if (n < min_samples) return(fail("insufficient_samples", "Too few complete observations."))
  if (n < p) return(fail("insufficient_df", "Fewer observations than design columns."))
  fit <- stats::lm.fit(x, y, tol = 1e-7)
  out$rank <- fit$rank
  out$df_residual <- n - fit$rank
  if (fit$rank < p) {
    return(fail("rank_deficient", "Non-estimable design; check constant variables, absent factor levels, or dependent SNPs/covariates."))
  }
  if (out$df_residual < min_residual_df) {
    return(fail("insufficient_df", "Too few residual degrees of freedom."))
  }
  total_ss <- sum((y - mean(y))^2)
  if (total_ss == 0) return(fail("constant_methylation", "Methylation is constant among complete observations."))
  rss <- sum(fit$residuals^2)
  if (!is.finite(rss) || !is.finite(total_ss) ||
      rss <= total_ss * 100 * .Machine$double.eps^2) {
    return(fail("degenerate_fit", "Residual variation is zero, numerically negligible, or non-finite."))
  }
  variance <- rss / out$df_residual
  # QR is pivoted; put the covariance back into the original column order.
  inverse <- chol2inv(qr.R(fit$qr)[seq_len(p), seq_len(p), drop = FALSE])
  covariance <- matrix(0, p, p)
  covariance[fit$qr$pivot, fit$qr$pivot] <- inverse
  estimates <- unname(fit$coefficients)
  se <- sqrt(diag(covariance) * variance)
  t_statistic <- estimates / se
  critical <- stats::qt((1 + conf_level) / 2, out$df_residual)
  k <- length(tested)
  b <- estimates[tested]
  f_statistic <- as.numeric(crossprod(b, solve(covariance[tested, tested, drop = FALSE], b))) /
    (k * variance)
  out$f_statistic <- max(0, f_statistic)
  out$p_value <- stats::pf(out$f_statistic, k, out$df_residual, lower.tail = FALSE)
  out$coefficients <- data.frame(
    coefficient = colnames(x), estimate = estimates, std_error = se,
    t_statistic = t_statistic,
    p_value = 2 * stats::pt(abs(t_statistic), out$df_residual, lower.tail = FALSE),
    conf_low = estimates - critical * se, conf_high = estimates + critical * se,
    stringsAsFactors = FALSE
  )
  if (k == 1L) {
    out$estimate <- estimates[tested]
    out$std_error <- se[tested]
    out$conf_low <- out$coefficients$conf_low[tested]
    out$conf_high <- out$coefficients$conf_high[tested]
  }
  out
}

.empty_result <- function(cpg, n_total, snps, available, missing, k) {
  data.frame(cpg_id = cpg, n = NA_integer_, n_excluded = NA_integer_,
             n_aligned = n_total, n_parameters = NA_integer_, rank = NA_integer_,
             df_residual = NA_integer_, exposure_df = k,
             n_snps_requested = length(snps), n_snps_available = length(available),
             n_snps_in_model = 0L, snps_requested = paste(snps, collapse = ";"),
             snps_in_model = "", snps_missing = paste(missing, collapse = ";"),
             adjustment = if (!length(snps)) "unmapped" else if (!length(missing)) "full" else if (length(available)) "partial" else "unavailable",
             estimate = NA_real_, std_error = NA_real_, conf_low = NA_real_,
             conf_high = NA_real_, f_statistic = NA_real_, p_value = NA_real_,
             status = "ok", message = "", stringsAsFactors = FALSE)
}

.empty_coefficients <- function() {
  data.frame(cpg_id = character(), coefficient = character(), role = character(),
             estimate = numeric(), std_error = numeric(), t_statistic = numeric(),
             p_value = numeric(), conf_low = numeric(), conf_high = numeric(),
             n = integer(), df_residual = integer(), stringsAsFactors = FALSE)
}

.adjust_results <- function(results, method) {
  # The testing family includes every requested CpG, including failed models.
  results$p_adjust <- stats::p.adjust(results$p_value, method = method, n = nrow(results))
  results$p_bonferroni <- stats::p.adjust(results$p_value, method = "bonferroni", n = nrow(results))
  results
}


# Module: ewas.R
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


# Module: io.R
#' Convert a samples-by-CpGs table to the standard CpGs-by-samples matrix
#' @export
methylation_from_table <- function(data, sample_id = "sample_id", cpg_cols = NULL) {
  if (!is.data.frame(data)) stop("data must be a data.frame.", call. = FALSE)
  data <- as.data.frame(data)
  if (!is.character(sample_id) || length(sample_id) != 1L ||
      is.na(sample_id) || !sample_id %in% names(data)) {
    stop("sample_id must name an existing column.", call. = FALSE)
  }
  .check_ids(names(data), "table column names")
  .check_ids(as.character(data[[sample_id]]), "table sample IDs")
  if (is.null(cpg_cols)) cpg_cols <- setdiff(names(data), sample_id)
  .check_ids(cpg_cols, "cpg_cols")
  if (!all(cpg_cols %in% names(data)) || sample_id %in% cpg_cols) {
    stop("cpg_cols must select existing CpG columns, excluding sample_id.", call. = FALSE)
  }
  if (!all(vapply(data[cpg_cols], is.numeric, logical(1)))) {
    stop("All selected CpG columns must be numeric.", call. = FALSE)
  }
  out <- t(as.matrix(data[cpg_cols]))
  colnames(out) <- as.character(data[[sample_id]])
  out
}

#' Extract only requested SNPs from an attached bigsnpr bigSNP object
#' @export
genotypes_from_bigsnp <- function(x, snps, sample_ids = NULL,
                                 sample_col = "sample.ID", snp_col = "marker.ID") {
  if (!inherits(x, "bigSNP")) stop("x must be an attached bigsnpr bigSNP object.", call. = FALSE)
  if (!sample_col %in% names(x$fam) || !snp_col %in% names(x$map)) {
    stop("sample_col or snp_col is missing from the bigSNP metadata.", call. = FALSE)
  }
  ids <- as.character(x$fam[[sample_col]])
  variants <- as.character(x$map[[snp_col]])
  .check_ids(ids, "bigSNP sample IDs")
  .check_ids(variants, "bigSNP variant IDs")
  if (!is.character(snps)) stop("snps must be a character vector.", call. = FALSE)
  .check_ids(snps, "requested SNP IDs")
  if (is.null(sample_ids)) sample_ids <- ids
  sample_ids <- as.character(sample_ids)
  .check_ids(sample_ids, "requested sample IDs")
  if (!all(snps %in% variants) || !all(sample_ids %in% ids)) {
    stop("Some requested SNPs or samples are absent from the bigSNP object.", call. = FALSE)
  }
  if (!identical(as.integer(dim(x$genotypes)), as.integer(c(length(ids), length(variants))))) {
    stop("bigSNP genotype dimensions do not match its sample and variant metadata.", call. = FALSE)
  }
  out <- as.matrix(x$genotypes[match(sample_ids, ids), match(snps, variants), drop = FALSE])
  dimnames(out) <- list(sample_ids, snps)
  .check_matrix(out, "extracted genotypes")
  out
}

#' Export an analysis without silently overwriting previous results
#' @export
write_ewas <- function(x, directory, overwrite = FALSE) {
  if (!inherits(x, "epicurious_ewas")) stop("x must be an epicurious_ewas result.", call. = FALSE)
  .check_flag(overwrite, "overwrite")
  if (!is.character(directory) || length(directory) != 1L || is.na(directory) || !nzchar(directory)) {
    stop("directory must be a non-empty path.", call. = FALSE)
  }
  names <- c("analysis.rds", "results.csv", "coefficients.csv", "sample_audit.csv",
             "snp_audit.csv", "metadata.txt", "sessionInfo.txt")
  if (!is.null(x$unadjusted)) names <- c(names, "unadjusted_same_samples.csv")
  targets <- file.path(directory, names)
  if (any(file.exists(targets)) && !overwrite) {
    stop("Output files already exist; choose a new directory or overwrite = TRUE.", call. = FALSE)
  }
  # Prevent an old comparison table surviving a different overwritten analysis.
  stale_comparison <- file.path(directory, "unadjusted_same_samples.csv")
  if (is.null(x$unadjusted) && file.exists(stale_comparison)) {
    if (!overwrite) stop("A prior comparison file exists in directory.", call. = FALSE)
    if (!file.remove(stale_comparison)) stop("Could not remove prior comparison file.", call. = FALSE)
  }
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
    stop("Cannot create output directory.", call. = FALSE)
  }
  saveRDS(x, targets[1])
  for (nm in c("results", "coefficients", "sample_audit", "snp_audit")) {
    utils::write.csv(x[[nm]], file.path(directory, paste0(nm, ".csv")), row.names = FALSE, na = "")
  }
  if (!is.null(x$unadjusted)) {
    utils::write.csv(x$unadjusted, stale_comparison, row.names = FALSE, na = "")
  }
  metadata <- x$metadata
  metadata$session_info <- NULL
  writeLines(utils::capture.output(dput(metadata)), file.path(directory, "metadata.txt"))
  writeLines(utils::capture.output(print(x$metadata$session_info)), file.path(directory, "sessionInfo.txt"))
  invisible(stats::setNames(normalizePath(targets, mustWork = TRUE), names))
}


# Module: simulate.R
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
