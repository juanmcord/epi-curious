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
