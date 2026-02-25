
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bioshmods

`bioshmods` provides reusable Shiny modules for different tasks in
exploration of different bioinformatics data, such as differential gene
expression results and other high-throughput datasets. The modules
include interactive tables, plots, and data export functionality to
facilitate creating of custom bioinformatic Shiny apps.

## Included modules

- `geneBrowserTableServer()` / `geneBrowserPlotServer()` for table-first
  gene browsing and expression profiling.
- `volcanoServer()` and `discoServer()` for interactive contrast
  visualization.
- `pcaServer()` for PCA exploration.
- `tmodBrowserTableServer()`, `tmodBrowserPlotServer()`, and
  `tmodPanelPlotServer()` for tmod gene-set results.
- `geneGroupSelectorServer()` to build gene lists by name, expression,
  or differential-expression criteria.
- `fileExportServer()` to export data frames and R objects (`xlsx`,
  `rds`, `zip`).

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pkg_install("bihealth/bioshmods")
```

Alternative:

``` r
# install.packages("remotes")
remotes::install_github("bihealth/bioshmods")
```

## Example data

The package ships with `C19`, a list containing differential-expression
contrasts, annotation, expression, and covariate data.

``` r
library(bioshmods)
data(C19)
names(C19)
#> [1] "contrasts"  "annotation" "expression" "covariates"
names(C19$contrasts)
#> [1] "COVID19_ID0" "ICU_ID1"
dim(C19$expression)
#> [1] 558  32
head(colnames(C19$annotation))
#> [1] "PrimaryID" "ENSEMBL"   "SYMBOL"    "ENTREZID"  "REFSEQ"    "GENENAME"
```

## Quick start

Launch the integrated example gene browser:

``` r
library(bioshmods)
data(C19)

annot_linkout <- list(
  SYMBOL = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s",
  ENTREZID = "https://www.ncbi.nlm.nih.gov/gene/?term=%s"
)

gene_browser(C19, annot_linkout = annot_linkout)
```

## Data expectations

Many of the modules use similar data structures, for example:

- `contrasts`: named list of data frames; each contrast should include
  `PrimaryID`, `log2FoldChange`, and a p-value column (`pvalue` or
  `padj`, depending on module).
- `annotation`: data frame containing at least `PrimaryID`.
- `expression`: expression matrix/data frame with gene IDs in row names.
- `covariates`: sample metadata aligned to expression columns.

## Communication between modules

Modules communicate through shared reactive objects. The general pattern
is:

- create one `reactiveValues()` object in the app server
- pass it to a “producer” module (writes selected IDs)
- pass the same object to a “consumer” module (reacts to changes)

For example, in the gene browser, `geneBrowserTableServer()` writes
`gene_id$id` and `gene_id$ds`, and `geneBrowserPlotServer()` observes
those values and updates the plot accordingly.

``` r
library(shiny)
library(bioshmods)
data(C19)

# creeate the UI
ui <- fluidPage(
  fluidRow(
    column(6, geneBrowserTableUI("tbl", names(C19$contrasts))),
    column(6, geneBrowserPlotUI("gplot", contrasts = TRUE))
  )
)

# create the server
server <- function(input, output, session) {

  # this is the reactive value that will be 
  # used for communication between the table and plot modules
  gene_id <- reactiveValues()

  geneBrowserTableServer(
    "tbl",
    cntr = C19$contrasts,
    annot = C19$annotation,
    gene_id = gene_id
  )

  geneBrowserPlotServer(
    "gplot",
    gene_id = gene_id,
    covar = C19$covariates,
    exprs = C19$expression,
    annot = C19$annotation,
    cntr = C19$contrasts
  )
}

shinyApp(ui, server)
```

The same approach works with `volcanoServer()` and `discoServer()`,
which can also publish selected genes into the same shared `gene_id`
object.

## Minimal module app

``` r
ui <- shiny::fluidPage(
  bioshmods::volcanoUI("volc", datasets = names(C19$contrasts))
)

server <- function(input, output, session) {
  bioshmods::volcanoServer(
    "volc",
    cntr = C19$contrasts,
    annot = C19$annotation
  )
}

shiny::shinyApp(ui, server)
```
