\name{write_ewas}
\alias{write_ewas}
\title{Export EWAS results and audit information}
\description{Writes the complete analysis as RDS, result and audit tables as
CSV, and metadata and session information as text.}
\usage{write_ewas(x, directory, overwrite = FALSE)}
\arguments{
\item{x}{An epicurious_ewas object.}
\item{directory}{Output directory, created recursively if necessary.}
\item{overwrite}{Whether existing analysis files may be replaced. Default FALSE.}
}
\value{Invisibly, named absolute paths of written files.}
\details{CSV exports do not contain raw genotype or methylation matrices, but
sample audit and RDS metadata include sample IDs. A write is not transactional;
an I/O error can leave partial files. Unrelated files are preserved. An obsolete
comparison CSV is removed when overwriting with an analysis without comparison.}
