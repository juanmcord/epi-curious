\name{simulate_epicurious}
\alias{simulate_epicurious}
\title{Generate a synthetic epi-curious example}
\description{Creates reproducible methylation, phenotype, genotype, and mapping
inputs with known additive effects. Preserves the caller's random-number state.}
\usage{simulate_epicurious(n = 200L, n_cpgs = 12L, seed = 42L)}
\arguments{
\item{n}{Number of samples, at least 30.}
\item{n_cpgs}{Number of CpGs, at least 4.}
\item{seed}{Non-negative integer random seed.}
}
\value{A list of methylation, phenotypes, genotypes, cpg_snp_map, truth, and
methylation_scale. All identifiers and measurements are synthetic.}
\examples{d <- simulate_epicurious(); dim(d$methylation)}
