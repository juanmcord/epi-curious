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
