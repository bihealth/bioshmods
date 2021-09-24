.gethovertext <- function(covar, ids) {

  ids <- setdiff(colnames(covar), ids)
  if(length(ids) < 1) {
    ids <- colnames(covar)
  }

  tmp <- map_dfc(rlang::set_names(ids), ~ paste0(.x, ": ", covar[[.x]]))
  apply(tmp, 1, function(x) paste(x, collapse="\n"))
}

.getcol <- function(x, df) {
  if(x == "N/A" || ! x %in% colnames(df)) {
    return(NULL)
  }
  return(as.formula(paste("~", x)))
}


.interesting_covariates <- function(covar) {
  covar %>% summary_colorDF() %>%
    dplyr::filter(unique > 1 & (unique < nrow(covar) | .data[["Class"]] == '<dbl>')) %>%
    pull("Col")
}

.categorical_covs <- function(covar) {
  categorical_covs <- covar %>% summary_colorDF() %>%
    dplyr::filter(unique > 1 & unique < nrow(covar) & .data[["Class"]] != '<dbl>') %>%
    pull("Col")
}
 
#' @rdname pcaServer
#' @export
pcaUI <- function(id, datasets=NULL) {



  if(is.null(datasets)) {
    ds_selector <- hidden(selectInput(NS(id, "dataset"), "Dataset:", "default", selected="default"))
  } else {
    ds_selector <- selectInput(NS(id, "dataset"), "Dataset:", datasets, selected=datasets[1])
  }

    sidebarLayout(
      sidebarPanel(
        ds_selector,
        uiOutput(NS(id, "color_ui")),
        uiOutput(NS(id, "symbol_ui")),
          checkboxInput(NS(id, "threeD"), "3D "),
        fluidRow(
          column(width=4, uiOutput(NS(id, "pca_x_ui"))),
          column(width=4, uiOutput(NS(id, "pca_y_ui"))),
          column(width=4, uiOutput(NS(id, "pca_z_ui")))
        ),
      width=3),
      mainPanel(
       plotlyOutput(NS(id, "pca_plot"), width="100%", height=600),
       width=9, useShinyjs()
      )
    )
}

#' Shiny Module – PCA plots
#'
#' Shiny Module – PCA plots
#' 
#' @section Datasets:
#'
#' Rather than specifying pca as a matrix and covariates as a data frame, 
#' these arguments can be, respectively, named lists of matrices and of data frames. This
#' allows to switch between different data sets (e.g. with different
#' samples, covariates) or between different representations (e.g.
#' different PCA variants or other transformations such as Umap). Note that
#' if the pca argument is a list, then the covar argument must be a list,
#' too, and all the names of the pca list must be also in the covar list.
#' 
#' @param id identifier of the shiny module (character vector)
#' @param pca pca matrix – columns correspond to principal components, rows
#' to observations. Rows must be named and must correspond to the ID column 
#' of the covariate data frame. 
#' @param covar data frame containing covariates. The identifiers of the
#' samples in the covariate data frame are taken from the ID column (by
#' default, "ID").
#' @param datasets if not NULL, a character vector specifying the data sets
#' (see Details)
#' @param colorBy selected covariate to use for coloring the plot
#' @param symbolBy selected covariate to use for symbols on the plot
#' @param threeD whether the plot should be three-dimensional by default
#' @param idcol name of the ID column in the covariate data frame.
#' @importFrom plotly renderPlotly plotlyOutput
#' @importFrom shiny isTruthy
#' @importFrom stats as.formula prcomp
#' @importFrom shinyjs hidden
#' @importFrom shinyBS tipify
#' @examples
#' if(interactive()) {
#'   data(iris)
#'   covar <- iris
#'  
#'   pca <- prcomp(iris[,1:4], scale.=TRUE)
#'  
#'   ui <- fluidPage(pcaUI("pca"))
#'  
#'   server <- function(input, output, session) {
#'     pcaServer("pca", pca$x, covar, colorBy="Species")
#'   }
#'  
#'   shinyApp(ui, server)
#' }
#' @export
pcaServer <- function(id, pca, covar, idcol="ID", threeD=FALSE, colorBy=NULL, symbolBy=NULL) {

  stopifnot(is.character(id))
  stopifnot(is.logical(threeD))

  if(is(covar, "list") && !is.data.frame(covar) && is.data.frame(covar[[1]])) {
    message("pcaServer: running in multilevel mode")
  } else {
    pca   <- list(default=pca)
    covar <- list(default=covar)
  }

  data <- .pcaServer_check_args(pca, covar, idcol, colorBy, symbolBy)
  pca   <- data$pca
  covar <- data$covar
  df    <- data$df
  


  moduleServer(id, function(input, output, session) {
    disable("z")

    output$pca_x_ui <- renderUI({
      ds <- input$dataset
      if(is.null(ds)) { return("") }
      selectInput(NS(id, "x"), "X: ", choices=colnames(pca[[ds]]), 
                                                    selected=colnames(pca[[ds]])[1], width="100%")
    })
    output$pca_y_ui <- renderUI({
      ds <- input$dataset
      if(is.null(ds)) { return("") }
      selectInput(NS(id, "y"), "Y: ", choices=colnames(pca[[ds]]), 
                                                    selected=colnames(pca[[ds]])[2], width="100%")
    })

    output$pca_z_ui <- renderUI({
      ds <- input$dataset
      if(is.null(ds)) { return("") }
      selectInput(NS(id, "z"), "Z: ", choices=colnames(pca[[ds]]), 
                                                    selected=colnames(pca[[ds]])[3], width="100%")
    })

    output$color_ui <- renderUI({
      ds <- input$dataset
      if(is.null(ds)) { return("") }

      interesting_covariates <- .interesting_covariates(covar[[ds]])

       if(is.null(colorBy) || !colorBy %in% interesting_covariates) {
         colorBy <- interesting_covariates[1]
       }

      selectInput(NS(id, "color"), "Color by:",   
                    c("N/A", interesting_covariates), selected=colorBy)
    })

    output$symbol_ui <- renderUI({
      ds <- input$dataset
      if(is.null(ds)) { return("") }

      interesting_covariates <- .interesting_covariates(covar[[ds]])
      categorical_covs       <- .categorical_covs(covar[[ds]])

       if(is.null(symbolBy) || !symbolBy %in% interesting_covariates) {
         symbolBy <- "N/A"
       }

       selectInput(NS(id, "symbol"), "Symbol by:", 
                   c("N/A", categorical_covs), selected=symbolBy)
    })

    output$pca_plot <- renderPlotly({ 
      ds <- input$dataset

      x <- input$x
      y <- input$y
      z <- input$z
      if(!(isTruthy(x) && isTruthy(y) && isTruthy(z))) {
        return(NULL)
      }
 
      symbol <- .getcol(input$symbol, df[[ds]])
      color  <- .getcol(input$color, df[[ds]])
 
      if(isTruthy(input$threeD)) {
        enable("z")
        plot_ly(data=df[[ds]], type="scatter3d", x=df[[ds]][[x]], y=df[[ds]][[y]], z=df[[ds]][[z]], 
                mode="markers", color=color, symbol=symbol,
                hoverinfo="hoverinfo")
      } else {
        disable("z")
        plot_ly(data=df[[ds]], type="scatter", x=df[[ds]][[x]], y=df[[ds]][[y]], 
                mode="markers", color=color, symbol=symbol,
                hovertext=df[[ds]][["hoverinfo"]])
      }
    })

  })
}

## check arguments to the server for corectness
.pcaServer_check_args <- function(pca, covar, idcol, colorBy, symbolBy) {

  stopifnot(!is.null(names(pca)))
  stopifnot(!is.null(names(covar)))

  stopifnot(all(names(pca) %in% names(covar))) 

  df <- list()

  for(c in names(pca)) {
    stopifnot(all(c(colorBy, symbolBy) %in% colnames(covar[[c]])))

    if(!(is.matrix(pca[[c]]) || is.data.frame(pca[[c]]))) {
      stop(
           sprintf("PCA object must be either a matrix or a data frame [data set: %s]", c)
           )
    }

    if(is.null(colnames(pca[[c]]))) {
      colnames(pca[[c]]) <- paste0("PC", 1:ncol(pca[[c]]))
    }

    if(any(ops <- colnames(pca[[c]]) %in% colnames(covar[[c]]))) {
      colnames(covar[[c]])[ops] <- paste(colnames(covar[[c]]), "_covariate")
    }

    df[[c]] <- cbind(covar[[c]], pca[[c]])

    if(is.null(df[[c]][["hoverinfo"]])) {
      df[[c]][["hoverinfo"]] <- .gethovertext(covar[[c]], names(covar[[c]]))
    }

  }

  return(list(pca=pca, covar=covar, df=df))
}


#' Simple PCA plot in a browser
#'
#' Simple PCA plot in a browser using shiny
#'
#' @param covar covariate file
#' @param x matrix or data frame for PCA
#' @param colorBy name of the covariate column for coloring
#' @param symbolBy name of the covariate column for chosing symbols
#' @examples
#' if(interactive()) {
#'    plot_pca_shiny(iris[,1:4], iris[,5,drop=FALSE])
#' }
#' @export
plot_pca_shiny <- function(x, covar, colorBy=NULL, symbolBy=NULL) {

  pca <- prcomp(x, scale.=TRUE)

  ui <- fluidPage(pcaUI("pca")) 
                  

  server <- function(input, output, session) {
    pcaServer("pca", pca$x, covar, colorBy=colorBy, symbolBy=symbolBy)
  }

  shinyApp(ui, server)
}




