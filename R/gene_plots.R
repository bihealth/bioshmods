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


## Wrapper around plot_gene, mainly to replace "N/A" with NA
.gene_browser_plot <- function(covar, id, covarName, rld, annot, 
                               groupBy = "N/A", colorBy = "N/A", symbolBy = "N/A", trellisBy="N/A") {
  .args <- list(id=id, xCovar=covarName, covar=covar, exprs=rld, groupBy=groupBy, annot=annot,
                colorBy=colorBy, symbolBy=symbolBy, trellisBy=trellisBy)
  ## weirdly, the line below is really, really slow
  #.args <- map(.args, ~ if(!is.na(.x) && .x == "N/A") { NA } else { .x })
  if(.args$groupBy == "N/A") .args$groupBy <- NA
  if(.args$colorBy == "N/A") .args$colorBy <- NA
  if(.args$symbolBy == "N/A") .args$symbolBy <- NA
  if(.args$trellisBy == "N/A") .args$trellisBy <- NA
  do.call(plot_gene, .args)
}



.dynamic_col_control <- function(id, covar, datasets, ds_selected) {

  all_covars         <- covar %>% summary_colorDF() %>% filter(unique > 1) %>% pull(.data$Col)
  non_unique         <- covar %>% summary_colorDF() %>% filter(unique < nrow(covar)) %>% pull(.data$Col)
  default_covar <- .default_covar(covar, all_covars, default="group")

  ds_selector <- selectInput(NS(id, "dataset"), "Dataset", choices=datasets, selected=ds_selected) 
  if(length(datasets) < 2L) {
    ds_selector <- hidden(ds_selector)
  }

  tagList(
      fluidRow(ds_selector),
      column(width=5, 
      fluidRow(
               tipify(selectInput(NS(id, "covarName"), "X covariate", non_unique, selected=default_covar, width="100%"),
                      "Variable shown on the X axis", placement="right")),
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
      fluidRow(tipify(selectInput(NS(id, "trellisBy"), "Trellis by", c("N/A", non_unique), selected="N/A", width="100%"),
                      "Each unique value of the variable will be shown on a separate subplot", placement="right")),
      fluidRow(tipify(numericInput(NS(id, "fontSize"),    label="Font size", min=6, value=14, step=1, width="50%"),
                      "Change the base font size of the figure", placement="right")),
      offset=1),
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
               column(width=11,
                      withSpinner(plotOutput(NS(id, "countsplot")))))

  if(contrasts) {
    return(sidebarLayout(col_control,
      mainPanel(
      column(9, style="padding:20px;", tabsetPanel(
      tabPanel("Plot", fluidRow(br(), plot_ui)),
      tabPanel("Contrast overview", fluidRow(br(), DTOutput(NS(id, "contr_sum"))))
      )))))
  } else {
    return(sidebarLayout(
                         col_control,
      mainPanel(plot_ui)))
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
#' @param exprs expression matrix; row names must correspond to the primary identifiers
#' @param contrasts (logical) whether or not create an additional panel
#'        next to the plot which can be used to show detailed contrast
#'        information for a gene
#' @param annot (optional) annotation data frame containing column 'PrimaryID'
#'        corresponding to the rownames of the contrast data frames
#' @param annot_linkout a list; see Details. 
#' @param id module identifier (same as the one passed to geneBrowserTableUI)
#' @param covar data frame with all covariates
#' @param cntr (optional) list of contrasts
#' @param symbol_col name of the column in `annot` which contains the gene
#'        symbols; use NULL if no such column
#' @param description_col name of the column in `annot` which contains the gene
#'        title / description; use NULL if no such column
#' @return does not return anything useful
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
#'      geneBrowserPlotServer("gplot", list(id="MUZG"), covar, mtx)
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
#'                           covar=C19$covariates, 
#'                           exprs=C19$expression,
#'                           annot=C19$annotation, 
#'                           cntr=C19$contrasts
#'      )
#'   }
#'   shinyApp(ui, server)
#' }
#' @seealso [geneBrowserTableServer()], and [gene_browser()] for example
#' code.
#' @export
geneBrowserPlotServer <- function(id, gene_id, covar, exprs, annot=NULL, cntr=NULL, 
                                  primary_id="PrimaryID", symbol_col="SYMBOL", description_col="GENENAME", 
                                  annot_linkout=NULL) {
## XXX make checks
# stopifnot(is.reactive(gene_id))
# stopifnot(!is.reactive(covar))
# stopifnot(!is.reactive(exprs))
# stopifnot(!is.reactive(annot))
# stopifnot(is.null(annot) || is.data.frame(annot))
# if(!is.null(annot)) {
#   stopifnot(all(c(primary_id, symbol_col, description_col) %in% 
#                 colnames(annot)))
# }

  if(is.data.frame(covar)) {
    covar <- list(default=covar)
    exprs <- list(default=exprs)
    annot <- list(default=annot)
    cntr  <- list(default=cntr)
  } else {
    message("geneBrowserPlotServer in multi dataset mode")
  }

  ## vector holding the names of all datasets
  datasets        <- names(covar)
  names(datasets) <- datasets

  moduleServer(id, function(input, output, session) {
    disable("save")

    if(!is.reactivevalues(gene_id)) {
      tmp <- gene_id
      gene_id <- reactiveValues()
      gene_id$id <- tmp$id
      if(is.null(gene_id$ds <- tmp$ds)) {
        gene_id$ds <- "default"
      }
    }

    ds   <- reactiveVal()
    g_id <- reactiveVal()

    observe({
      if(!is.null(gene_id)) {
        if(isTruthy(gene_id$ds)) { ds(gene_id$ds) }
        if(isTruthy(gene_id$id)) { g_id(gene_id$id) }
      }
    })

    observe({
      if(isTruthy(input$dataset)) {
        ds(input$dataset)
      }})

    ## Save figure as a PDF
    output$save <- downloadHandler(
      filename = function() {
        .id <- g_id()
        .ds <- ds()
        ret <- sprintf("expression_profile_ds_%s_%s_covarX_%s_colorBy_%s_groupBy_%s_symbolBy_%s_trellisBy_%s.pdf",
                       .ds, .id,
                       input$covarName, input$colorBy, input$groupBy, input$symbolBy, input$trellisBy)
        ret <- gsub("[^0-9a-zA-Z_.-]", "", ret)
        return(ret)
      },
      content = function(file) {
        pdf(file=file, width=8, height=5)
        g <- .gene_browser_plot(covar[[ ds() ]], g_id(), 
                                input$covarName, 
                                exprs[[ ds() ]], 
                                annot[[ ds() ]], 
                           input$groupBy, input$colorBy, input$symbolBy, input$trellisBy)
        print(g)
        dev.off()
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
      .gene_browser_info_tab(g_id(), covar[[g_id()]][[input$covarName]], exprs[[ds()]][ g_id(), ])
    })
 
    ## reload the plot interface only if the data set (and covariates)
    ## changed
    output$col_control <- renderUI({
      .ds <- ds()
      if(!isTruthy(.ds)) { .ds <- 1 }
      .dynamic_col_control(id, covar[[.ds]], names(covar), datasets[.ds])
    })
 
    ## The actual plot
    output$countsplot <- renderPlot({
      if(!isTruthy(ds()) || !isTruthy(g_id())) { return(NULL) }
      if(!isTruthy(input$covarName)) { return(NULL) }
      if(is.na(g_id())) { return(NULL) }
      enable("save")
      
      .gene_browser_plot(covar[[ds()]], g_id(), input$covarName, exprs[[ ds() ]], annot[[ ds() ]], 
                         input$groupBy, input$colorBy, input$symbolBy, input$trellisBy) +
                                    theme(text=element_text(size=input$fontSize))
    })
  })
}

