.gene_browser_plot_log <- function(...) {
  .bioshmods_log(..., .prefix="gene_browser_plot")
}

## prepare the additional gene info tab panel
.gene_browser_info_tab <- function(id, x, y, covar) {
     ret <- ""

     if(is.numeric(x)) {
       pearson.test  <- cor.test(x, y, use="p")
       spearman.test <- cor.test(x, y, use="p", method="s")
       ret <- paste0(ret,
         sprintf("Correlation: r=%.2f [p = %s], rho=%.2f [p = %s]",
                 pearson.test$estimate,
                 format.pval(pearson.test$p.value, digits=2),
                 spearman.test$estimate,
                 format.pval(spearman.test$p.value, digits=2)))
     }
     return(ret)
}

.default_covar <- function(covar, all_covars, default="group") {
  interesting_covars <- covar %>% 
      summary_colorDF() %>% 
      filter(unique < n()) %>% 
      pull(.data$Col)

  if(default %in% interesting_covars) {
    default_covar <- default
  } else {
    if(length(interesting_covars) > 0) {
      default_covar <- interesting_covars[1]
    } else {
      default_covar <- all_covars[1]
    }
  }

  return(default_covar)
}

.gene_browser_palette_scale <- function(palettes, dataset, color_by) {
  if(is.null(palettes) || !isTruthy(color_by) || identical(color_by, "N/A")) {
    return(NULL)
  }

  pal_all <- if(is.function(palettes)) palettes() else palettes
  if(!is.list(pal_all)) {
    return(NULL)
  }

  pal_ds <- pal_all
  if(isTruthy(dataset) && dataset %in% names(pal_all) && is.list(pal_all[[dataset]])) {
    pal_ds <- pal_all[[dataset]]
  } else if("default" %in% names(pal_all) && is.list(pal_all[["default"]])) {
    pal_ds <- pal_all[["default"]]
  }

  if(!is.list(pal_ds) || !color_by %in% names(pal_ds)) {
    return(NULL)
  }

  .cast_palette_to_ggplot(pal_ds[[color_by]])
}


## Wrapper around plot_gene, mainly to replace "N/A" with NA - the UI
## uses N/A because it is easier to understand than NA by some people, apparently
.gene_browser_plot <- function(covar, id, covarXName, covarYName, rld, annot, 
                               groupBy = "N/A", colorBy = "N/A", symbolBy = "N/A", trellisBy="N/A",
                               trellisScales = "fixed",
                               exprs_label = "Expression",
                               plotType = "box", colorScale=NULL,
                               transpose = FALSE) {
  if(covarYName == "Expression") { covarYName <- exprs_label }
  if(covarXName == "Expression") { covarXName <- exprs_label }

  .args <- list(id=id, xCovar=covarXName, yCovar=covarYName, covar=covar, exprs=rld, groupBy=groupBy, annot=annot,
                expressionLabel = exprs_label,
                colorBy=colorBy, symbolBy=symbolBy, trellisBy=trellisBy,
                trellisScales=trellisScales,
                categoricalPlot=plotType, colorScale=colorScale,
                transpose=transpose)

  ## weirdly, the line below is really, really slow
  #.args <- map(.args, ~ if(!is.na(.x) && .x == "N/A") { NA } else { .x })

  if(.args$groupBy   == "N/A") .args$groupBy   <- NA
  if(.args$colorBy   == "N/A") .args$colorBy   <- NA
  if(.args$symbolBy  == "N/A") .args$symbolBy  <- NA
  if(.args$trellisBy == "N/A") .args$trellisBy <- NA

  do.call(plot_gene, .args)
}

# Build the report-code chunk for the current gene-browser plot.
# Mirrors `.gene_browser_plot()` but returns code text instead of a ggplot.
.gene_browser_plot_chunk <- function(covar, id, rld, annot, input,
                                     exprs_label = "Expression",
                                     colorScale=NULL,
                                     env_map = NULL) {
  covarXName <- input$covarXName
  covarYName <- input$covarYName
  if(covarYName == "Expression") { covarYName <- exprs_label }
  if(covarXName == "Expression") { covarXName <- exprs_label }

  .args <- list(
    id=id,
    xCovar=covarXName,
    yCovar=covarYName,
    covar=covar,
    exprs=rld,
    groupBy=input$groupBy,
    annot=annot,
    expressionLabel=exprs_label,
    colorBy=input$colorBy,
    symbolBy=input$symbolBy,
    trellisBy=input$trellisBy,
    trellisScales=input$trellisScales,
    categoricalPlot=input$plotType,
    colorScale=colorScale,
    transpose=input$transposePlot,
    env_map=env_map
  )

  if(.args$groupBy   == "N/A") .args$groupBy   <- NA
  if(.args$colorBy   == "N/A") .args$colorBy   <- NA
  if(.args$symbolBy  == "N/A") .args$symbolBy  <- NA
  if(.args$trellisBy == "N/A") .args$trellisBy <- NA

  do.call(plot_gene_chunk, .args)$code
}



# create the interface dynamically, depending on the data - covariates
# and datasets
.dynamic_col_control <- function(id, covar, datasets, ds_selected) {

  covar_sum <- summary_colorDF(covar)
  all_covars         <- covar_sum %>% filter(unique > 1) %>% pull(.data$Col)
  non_unique         <- covar_sum %>% 
    filter(.data[["Class"]] %in% c("<dbl>", "<int>") | unique < nrow(covar)) %>% 
    pull(.data[["Col"]])
  non_unique         <- c(non_unique, "Expression")
  default_covar <- .default_covar(covar, all_covars, default="group")

  ds_selector <- selectInput(NS(id, "dataset"), "Dataset", choices=datasets, selected=ds_selected) 
  if(length(datasets) < 2L) {
    ds_selector <- hidden(ds_selector)
  }

  tagList(
      fluidRow(ds_selector),
      fluidRow(column(width=5, 
        fluidRow(
               tipify(selectInput(NS(id, "covarXName"), "X covariate", non_unique, selected=default_covar, width="100%"),
                      "Variable shown on the X axis", placement="right")),
        fluidRow(
               tipify(selectInput(NS(id, "covarYName"), "Y covariate", non_unique, selected="Expression", width="100%"),
                      "Variable shown on the Y axis", placement="right")),
        fluidRow(
               tipify(selectInput(NS(id, "colorBy"), "Color by", c("N/A", non_unique), selected="N/A", width="100%"),
                      "Variable coded as color", placement="right")),
        fluidRow(
               tipify(selectInput(NS(id, "symbolBy"), "Symbol by", c("N/A", non_unique), selected="N/A", width="100%"),
                      "Variable coded as symbol", placement="right")),
      ),
      column(width=5,
        fluidRow(
               tipify(selectInput(NS(id, "groupBy"), "Link data points by", c("N/A", non_unique), selected="N/A", width="100%"),
                      "Points with identical values will be linked by a line", placement="right")),
        fluidRow(
               tipify(selectInput(NS(id, "plotType"), "Categorical plot", c("box", "violin", "raincloud"), selected="box", width="100%"),
                      "Display style when X covariate is categorical and points are not linked", placement="right")),
        fluidRow(tipify(selectInput(NS(id, "trellisBy"), "Trellis by", c("N/A", non_unique), selected="N/A", width="100%"),
                      "Each unique value of the variable will be shown on a separate subplot", placement="right")),
        fluidRow(tipify(selectInput(NS(id, "trellisScales"), "Trellis scales", c("fixed", "free", "free_x", "free_y"), selected="fixed", width="100%"),
                      "Scale behavior for trellis facets; active only when Trellis by is set", placement="right")),
        fluidRow(tipify(numericInput(NS(id, "fontSize"),    label="Font size", min=6, value=14, step=1, width="50%"),
                      "Change the base font size of the figure", placement="right")),
        fluidRow(tipify(checkboxInput(NS(id, "transposePlot"), "Transpose plot", value=FALSE, width="100%"),
                      "Flip X/Y orientation of the whole plot", placement="right")),
       fluidRow(figsizeInput(NS(id, "figure_size"), width="100%", selected="800x600"),
                bsTooltip(NS(id, "figure_size"), 
                  "Change the figure size (in pixels), width x height. Press backspace to enter your own sizes.", placement="right")),
      offset=1),
      ),



      fluidRow(textOutput(NS(id, "addInfo"))),
      fluidRow(h3("Additional info:")),
      fluidRow(tableOutput(NS(id, "geneData")))
    )


}

#' @rdname geneBrowserPlotServer
#' @export
geneBrowserPlotUI <- function(id, contrasts=FALSE) {
  col_control <- 
    sidebarPanel(
                 uiOutput(NS(id, "col_control"))
                 )
  plot_ui <- 
      fluidRow(column(width=1, 
                      tipify(downloadButton(NS(id, "save"), "PDF", class="bg-success"),
                             "Save as PDF")),
               column(width=2, .report_add_button(id)),
               column(width=9,
                      withSpinner(plotOutput(NS(id, "countsplot"), height="100%", width="100%")))
      )

  if(contrasts) {
    return(sidebarLayout(col_control,
      mainPanel(
      column(9, style="padding:20px;", tabsetPanel(
      tabPanel("Plot", fluidRow(br(), plot_ui)),
      tabPanel("Contrast overview", fluidRow(br(), DTOutput(NS(id, "contr_sum")))),
      .report_code_tab(id)
      )))))
  } else {
    return(sidebarLayout(
                         col_control,
      mainPanel(
        tabsetPanel(
          tabPanel("Plot", fluidRow(br(), plot_ui)),
          .report_code_tab(id)
        )
      )))
  }

}

#' Shiny Module – gene browser expression profile plot
#'
#' Shiny Module – gene browser expression profile plot
#'
#' The `gene_id` parameter must be a reactive value, because that is the
#' whole point of the plotting module: observe changes to the gene ID and
#' update the plot accordingly.
#' 
#' In contrast, other parameters must not be reactive values. This may
#' change in future to allow for dynamic exchange of data sets.
#'
#' The parameter `annot_linkout` is a named list. Names must correspond to
#' columns from the annotation data frame. The elements of the list are
#' character strings containing URLs with the `%s` placeholder. For
#' example, if the column `ENSEMBL` contains ENSEMBL identifiers, you can
#' link out by specifying 
#'
#' ```
#' annot_linkout=list(ENSEMBL="https://www.ensembl.org/id/%s")
#' ```
#' @param gene_id primary identifier of the gene to show. This must be
#'        either a list containing at least the element `id` and possibly
#'        the element `ds` (if multiple datasets are used). Alternatively,
#'        it is a `reactiveValues` object with the same elements.
#' @param primary_id name of the column which holds the primary identifiers
#' @param exprs named list of dataset expression matrices/data frames; row names
#'        must correspond to the primary identifiers
#' @param contrasts (logical) whether or not create an additional panel
#'        next to the plot which can be used to show detailed contrast
#'        information for a gene
#' @param annot (optional) named list of dataset annotation data frames
#'        containing column 'PrimaryID'
#' @param annot_linkout a list; see Details. 
#' @param id module identifier (same as the one passed to geneBrowserTableUI)
#' @param covar named list of dataset covariate data frames
#' @param cntr (optional) named list of dataset contrast lists
#' @param exprs_label Label to be used for the expression values
#' @param palettes (optional) reactive value with current color palettes
#' @param report Optional report list. When supplied, clicking `Add to report`
#'        appends the current generated code chunk to `report$chunks`.
#' @return does not return anything useful
#' @importFrom shiny is.reactivevalues
#' @examples
#' mtx <- matrix(rnorm(40, mean=rep(c(0, 1), each=20)), nrow=1)
#' rownames(mtx) <- "MUZG"
#' covar <- data.frame(
#'                     em=rep(LETTERS[1:2], each=20),
#'                     pstrem=rep(letters[1:20], 2),
#'                     bzdrem=rnorm(40))
#'   
#' if(interactive()) {
#'    ui  <- fluidPage(geneBrowserPlotUI("gplot", FALSE))
#'    serv <- function(input, output, session) {
#'      geneBrowserPlotServer(
#'        "gplot",
#'        list(id="MUZG"),
#'        covar=list(default=covar),
#'        exprs=list(default=mtx)
#'      )
#'    }
#'    shinyApp(ui, serv)
#' }
#'
#' ## Example with the C19 dataset
#' data(C19)
#' if(interactive()) {
#'   ui <- fluidPage(
#'          fluidRow(selectizeInput("id", label="Search for a gene",
#'            choices=NULL),
#'          fluidRow(geneBrowserPlotUI("gplot", TRUE))
#'          ))
#'
#'   server <- function(input, output, session) {
#'     gene_id <- reactiveValues()
#'     updateSelectizeInput(session, "id", choices=C19$annotation$SYMBOL)
#'
#'     ## translate symbol to primary ID
#'     observeEvent(input$id, {
#'       nn <- match(input$id, C19$annotation$SYMBOL)
#'       gene_id$id <- C19$annotation$PrimaryID[ nn ]
#'     })
#'
#'     geneBrowserPlotServer("gplot", gene_id=gene_id, 
#'                           covar=list(default=C19$covariates),
#'                           exprs=list(default=C19$expression),
#'                           annot=list(default=C19$annotation),
#'                           cntr=list(default=C19$contrasts)
#'      )
#'   }
#'   shinyApp(ui, server)
#' }
#' @seealso [geneBrowserTableServer()], and [gene_browser()] for example
#' code.
#' @importFrom shiny updateSelectizeInput
#' @export
geneBrowserPlotServer <- function(id, gene_id, covar, exprs, annot=NULL, cntr=NULL, 
                                  primary_id="PrimaryID",
                                  annot_linkout=NULL,
                                  exprs_label = "Expression", palettes=NULL,
                                  report=NULL) {
  covar <- .check_dataset_data_frames(covar, "covar")
  datasets <- names(covar)
  exprs <- .check_dataset_matrices(exprs, "exprs", datasets=datasets)
  annot <- .check_dataset_data_frames(annot, "annot", datasets=datasets, allow_null=TRUE)
  cntr <- .check_dataset_contrasts(cntr, "cntr", datasets=datasets, allow_null=TRUE)

  # vector holding the names of all datasets
  names(datasets) <- datasets

  # Code chunks refer to the module-level dataset lists by selected dataset.
  env_maps <- lapply(datasets, function(dataset) {
    list(
      exprs=sprintf('exprs[["%s"]]', dataset),
      covar=sprintf('covar[["%s"]]', dataset)
    )
  })
  names(env_maps) <- datasets

  # start the module server
  moduleServer(id, function(input, output, session) {
    disable("save")

    # if gene_id is not a reactiveValues object, wrap it into one
    if(!is.reactivevalues(gene_id)) {
      tmp <- gene_id
      gene_id <- reactiveValues()
      gene_id$id <- tmp$id
      gene_id$ds <- tmp$ds %||% "default"
    }

    # ds holds the dataset; g_id holds the gene ID
    ds   <- reactiveVal()
    g_id <- reactiveVal()

    observe({
      # observe the "outside" gene_id object, and, if it changes, update
      # the internal ds and g_id objects
      if(!is.null(gene_id)) {
        if(isTruthy(gene_id$ds)) { ds(gene_id$ds) }
        if(isTruthy(gene_id$id)) { g_id(gene_id$id) }
      }
    })

    observe({
      # if the dataset changes, update the ds object
      if(isTruthy(input$dataset)) {
        ds(input$dataset)
      }})

    observeEvent(input$add_to_report, {
      .append_report_chunk(report, input$report_code)
    })

    observeEvent(input$add_to_report_main, {
      .append_report_chunk(report, input$report_code)
    })

    ## XXX for some reason this does not work
    observe({
      if(!isTruthy(input$trellisBy) || input$trellisBy == "N/A") {
        disable(NS(id, "trellisScales"))
      } else {
        enable(NS(id, "trellisScales"))
      }
    })

    fig_size <- reactiveValues(width=600, height=600)
    observeEvent(input$figure_size, {
      size <- sanitize_figsize(input$figure_size, default=c(800, 600))
      fig_size$width <- size$width
      fig_size$height <- size$height
    })

    ## Save figure as a PDF
    output$save <- downloadHandler(
      filename = function() {
        .id <- g_id()
        .ds <- ds()
        base <- sprintf(
          "expression_profile_ds_%s_%s_covarX_%s_covarY_%s_colorBy_%s_groupBy_%s_symbolBy_%s_trellisBy_%s_trellisScales_%s_plotType_%s_transpose_%s",
          .ds, .id,
          input$covarXName, input$covarYName, input$colorBy, input$groupBy, input$symbolBy, input$trellisBy,
          input$trellisScales, input$plotType, input$transposePlot
        )
        sprintf("%s.pdf", sanitize_filename(base, "expression_profile"))
      },
      content = function(file) {
        save_pdf(file=file, width=8, height=5, draw=function() {
          color_scale <- .gene_browser_palette_scale(
            palettes=palettes,
            dataset=ds(),
            color_by=input$colorBy
          )

          g <- .gene_browser_plot(covar[[ ds() ]], g_id(), 
                                  input$covarXName, 
                                  input$covarYName, 
                                  exprs[[ ds() ]], 
                                  annot[[ ds() ]], 
                                  input$groupBy, input$colorBy, input$symbolBy, input$trellisBy,
                                  trellisScales = input$trellisScales,
                                  exprs_label = exprs_label,
                                  plotType = input$plotType, colorScale=color_scale,
                                  transpose = input$transposePlot)
          print(g)
        })
      }
    )
 
    # Show a turbo card for a gene
    output$geneData <- renderTable({
      if(!isTruthy(ds()) || !isTruthy(g_id())) {
        return(NULL)
      }
 
      if(is.null(annot[[ ds() ]])) {
        ret <- data.frame(V1=primary_id, V2=g_id())
      } else {
        m <- match(g_id(), annot[[  ds()  ]][[ primary_id ]])
        ret <- annot[[ ds() ]][ m, , drop=FALSE ]

        if(!is.null(annot_linkout)) {
          ret <- .apply_annot_linkout(ret, annot_linkout[[ ds() ]])
        }

        ret <- data.frame(V1=colnames(ret), V2=t(ret))
      }

      colnames(ret) <- NULL
      return(ret)
    }, sanitize.text.function = function(x) x)

    ## summary contrasts table
    output$contr_sum <- renderDT({
      if(!isTruthy(ds()) || !isTruthy(g_id()) || is.null(cntr[[ ds() ]])) {
        return(NULL)
      }
      cn <- names(cntr[[ ds() ]])
      res <- imap_dfr(cntr[[ ds() ]], ~ .x %>% 
                      filter(.data[[ primary_id ]] == g_id()) %>% 
                  mutate(contrast=.y))
      res <- imap_dfr(cntr, ~ {
                        .ds <- .y
                        imap_dfr(.x, ~ {
                                   .x %>% filter(.data[[ primary_id ]] == g_id()) %>%
                                     mutate("Data set"=.ds, Contrast=.y)
                        })
                  })
 
      res <- res %>% relocate(all_of(c("Data set", "Contrast")))
      numcol <- map_lgl(res, is.numeric)
      res %>% datatable(escape=FALSE, selection='none', options=list(pageLength=5)) %>%
        formatSignif(columns=colnames(res)[numcol], digits=2)
    })
 
    ## Additional information - e.g. correlation coefficient if the
    ## covariate is numeric
    output$addInfo <- renderText({
      if(!isTruthy(ds()) || !isTruthy(g_id())) {
        return(NULL)
      }
      .gene_browser_info_tab(g_id(), covar[[g_id()]][[input$covarXName]], exprs[[ds()]][ g_id(), ])
    })
 
    ## reload the plot interface only if the data set (and covariates)
    ## changed
    output$col_control <- renderUI({
      .ds <- ds()
      if(!isTruthy(.ds)) { .ds <- 1 }
      .dynamic_col_control(id, covar[[.ds]], names(covar), datasets[.ds])
    })

    ## Keep the report-code tab in sync with the current plot inputs.
    ## This generates text only; the plot itself still uses `plot_gene()`.
    observe({
      if(!isTruthy(ds()) || !isTruthy(g_id())) { return(NULL) }
      if(!isTruthy(input$covarXName)) { return(NULL) }
      if(!isTruthy(input$covarYName)) { return(NULL) }
      if(is.na(g_id())) { return(NULL) }

      color_scale <- .gene_browser_palette_scale(
        palettes=palettes,
        dataset=ds(),
        color_by=input$colorBy
      )

      code <- .gene_browser_plot_chunk(
        covar[[ds()]],
        g_id(),
        exprs[[ds()]],
        annot[[ds()]],
        input,
        exprs_label=exprs_label,
        colorScale=color_scale,
        env_map=env_maps[[ds()]]
      )

      updateTextAreaInput(session, "report_code", value=code)
    })
 
    ## The actual plot. need to put inside "observe" to use the reactive
    ## figure size
    observe({ output$countsplot <- renderPlot({
      if(!isTruthy(ds()) || !isTruthy(g_id())) { return(NULL) }
      if(!isTruthy(input$covarXName)) { return(NULL) }
      if(!isTruthy(input$covarYName)) { return(NULL) }
      if(is.na(g_id())) { return(NULL) }
      enable("save")
      
      #message(sprintf("plotting started with dataset=%s and gene=%s", ds(), g_id()))
      color_scale <- .gene_browser_palette_scale(
        palettes=palettes,
        dataset=ds(),
        color_by=input$colorBy
      )

      .gene_browser_plot(covar[[ds()]], g_id(), input$covarXName, input$covarYName, exprs[[ ds() ]], annot[[ ds() ]], 
                         input$groupBy, input$colorBy, input$symbolBy, input$trellisBy,
                         trellisScales = input$trellisScales,
                         exprs_label = exprs_label, plotType = input$plotType, colorScale=color_scale,
                         transpose = input$transposePlot) +
                                    theme(text=element_text(size=input$fontSize))
    }, width=fig_size$width, height=fig_size$height) })
  })
}
