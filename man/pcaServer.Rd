% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny_pca.R
\name{pcaUI}
\alias{pcaUI}
\alias{pcaServer}
\title{Shiny Module – PCA plots}
\usage{
pcaUI(id, datasets = NULL)

pcaServer(
  id,
  pca,
  covar,
  idcol = "ID",
  threeD = FALSE,
  colorBy = NULL,
  symbolBy = NULL
)
}
\arguments{
\item{id}{identifier of the shiny module (character vector)}

\item{datasets}{if not NULL, a character vector specifying the data sets
(see Details)}

\item{pca}{pca matrix – columns correspond to principal components, rows
to observations. Rows must be named and must correspond to the ID column
of the covariate data frame.}

\item{covar}{data frame containing covariates. The identifiers of the
samples in the covariate data frame are taken from the ID column (by
default, "ID").}

\item{idcol}{name of the ID column in the covariate data frame.}

\item{threeD}{whether the plot should be three-dimensional by default}

\item{colorBy}{selected covariate to use for coloring the plot}

\item{symbolBy}{selected covariate to use for symbols on the plot}
}
\description{
Shiny Module – PCA plots
}
\section{Datasets}{


Rather than specifying pca as a matrix and covariates as a data frame,
these arguments can be, respectively, named lists of matrices and of data frames. This
allows to switch between different data sets (e.g. with different
samples, covariates) or between different representations (e.g.
different PCA variants or other transformations such as Umap). Note that
if the pca argument is a list, then the covar argument must be a list,
too, and all the names of the pca list must be also in the covar list.
}

\examples{
if(interactive()) {
  data(iris)
  covar <- iris
 
  pca <- prcomp(iris[,1:4], scale.=TRUE)
 
  ui <- fluidPage(pcaUI("pca"))
 
  server <- function(input, output, session) {
    pcaServer("pca", pca$x, covar, colorBy="Species")
  }
 
  shinyApp(ui, server)
}
}
