# Source this from an extracted project without installing the package:
# source("epicurious/standalone.R")
# Functions are loaded into the environment in which this file is sourced.
local({
  source_path <- sys.frame(1)$ofile
  if (is.null(source_path)) stop("Load this file with source().")
  root <- dirname(normalizePath(source_path, mustWork = TRUE))
  for (f in sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))) {
    sys.source(f, envir = parent.env(environment()))
  }
})
