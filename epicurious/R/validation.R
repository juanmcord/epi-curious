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
