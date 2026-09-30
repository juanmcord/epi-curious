#!/usr/bin/env Rscript
# Rscript epicurious/run_ewas.R /path/to/config.R
# Rscript epicurious/run_ewas.R --demo /path/to/new_output_directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || args[1] %in% c("-h", "--help")) {
  cat("Usage:\n  Rscript run_ewas.R config.R\n  Rscript run_ewas.R --demo output_directory\n")
  quit(status = if (length(args)) 0L else 1L)
}
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[startsWith(script_arg, "--file=")][1])
root <- dirname(normalizePath(script_path, mustWork = TRUE))
for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))) source(f)

if (identical(args[1], "--demo")) {
  if (length(args) != 2L) stop("--demo requires exactly one output directory.")
  d <- simulate_epicurious()
  fit <- ewas_mqtl(d$methylation, d$phenotypes, ~ exposure + age + sex, "exposure",
                   d$genotypes, d$cpg_snp_map, methylation_scale = "M",
                   compare_unadjusted = TRUE)
  write_ewas(fit, args[2])
} else {
  if (length(args) != 1L) stop("Provide exactly one configuration path.")
  config_file <- normalizePath(args[1], mustWork = TRUE)
  config_env <- new.env(parent = globalenv())
  sys.source(config_file, envir = config_env)
  if (!exists("config", config_env, inherits = FALSE) || !is.list(config_env$config)) {
    stop("Configuration must define a list named config.")
  }
  cfg <- config_env$config
  required <- c("methylation_file", "phenotypes_file", "output_dir", "formula", "exposure", "methylation_scale")
  if (!all(required %in% names(cfg))) stop("Missing config entries: ", paste(setdiff(required, names(cfg)), collapse = ", "))
  file_keys <- c("methylation_file", "phenotypes_file", "genotypes_file", "cpg_snp_map_file")
  allowed <- c(file_keys, "output_dir", setdiff(names(formals(ewas_mqtl)),
                                                c("methylation", "phenotypes", "genotypes", "cpg_snp_map")))
  if (any(!names(cfg) %in% allowed)) stop("Unknown config entries: ", paste(setdiff(names(cfg), allowed), collapse = ", "))
  resolve <- function(path) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) stop("Paths must be non-empty strings.")
    if (grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)) path.expand(path)
    else file.path(dirname(config_file), path)
  }
  read_input <- function(key) {
    if (is.null(cfg[[key]])) return(NULL)
    readRDS(resolve(cfg[[key]]))
  }
  inputs <- list(methylation = read_input("methylation_file"),
                 phenotypes = read_input("phenotypes_file"),
                 genotypes = read_input("genotypes_file"),
                 cpg_snp_map = read_input("cpg_snp_map_file"))
  fit <- do.call(ewas_mqtl, c(inputs, cfg[setdiff(names(cfg), c(file_keys, "output_dir"))]))
  write_ewas(fit, resolve(cfg$output_dir))
}
print(fit)
