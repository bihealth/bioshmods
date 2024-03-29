% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{C19}
\alias{C19}
\alias{C19_gs}
\title{Example COVID-19 expression data set}
\format{
An object of class \code{list} of length 4.

An object of class \code{list} of length 4.
}
\source{
\url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103}
}
\usage{
C19

C19_gs
}
\description{
Example COVID-19 expression and analysis data set
}
\details{
This data set was derived from a standard analysis of GSE157103 data
set. The GSE157103 data set includes samples from COVID patients and
non-covid patients as well the
information whether patients had to be cared for at an intensive care
unit (ICU).

The C19 data set contains contrasts, covariates, expression matrix and annotation.

The C19_gs data set contains gene set collections (object \code{tmod_dbs}),
results of gene set enrichment analysis (object \code{tmod_res}).
}
\keyword{datasets}
