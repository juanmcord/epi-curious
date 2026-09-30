\name{methylation_from_table}
\alias{methylation_from_table}
\title{Convert a samples-by-CpGs data frame}
\description{Converts a numeric wide methylation table to the CpGs-by-samples
matrix required by \code{ewas_mqtl}, preserving exact identifiers.}
\usage{methylation_from_table(data, sample_id = "sample_id", cpg_cols = NULL)}
\arguments{
\item{data}{A data frame with one row per sample.}
\item{sample_id}{Name of its unique sample ID column.}
\item{cpg_cols}{Character vector of CpG column names. NULL selects every column
except the ID column; provide explicit names when metadata are present.}
}
\value{A numeric CpGs-by-samples matrix with named dimensions.}
\examples{
tab <- data.frame(sample_id = c("s1", "s2"), cg1 = c(0.2, 0.5), cg2 = c(0.7, 0.3))
methylation_from_table(tab)
}
