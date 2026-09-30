\name{genotypes_from_bigsnp}
\alias{genotypes_from_bigsnp}
\title{Extract selected SNPs from an attached bigSNP object}
\description{Reads only the specified sample and SNP subset into a numeric
matrix. The caller attaches the object with bigsnpr before calling this adapter.}
\usage{genotypes_from_bigsnp(x, snps, sample_ids = NULL,
  sample_col = "sample.ID", snp_col = "marker.ID")}
\arguments{
\item{x}{An attached object of class bigSNP, with fam, map, and genotypes components.}
\item{snps}{Unique character vector of exact variant IDs, in output order.}
\item{sample_ids}{Unique sample IDs, in output order. NULL uses all samples.}
\item{sample_col}{Sample identifier column in x$fam; identifiers must be unique.
For repeated individual IDs across families, prepare unique composite IDs first.}
\item{snp_col}{Variant identifier column in x$map; identifiers must be unique.}
}
\value{Numeric samples-by-SNPs matrix with matching dimension names.}
\details{Unavailable requested IDs are errors. This adapter does not impute,
recode alleles, resolve duplicate variant IDs, or perform genotype quality control.
Large selected subsets require corresponding memory.}
