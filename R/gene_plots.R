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

## figure out the default covariate to be pre-selected in the UI
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

.add_chunk <- function(id, code, rmd_var, gene_id, dataset, mode, fig_size) { 
      message(sprintf("generate rmd started with dataset=%s and gene=%s", dataset, gene_id))

      if(mode == "multi") {
        title <- sprintf("Expression profile of gene %s in dataset %s", gene_id, dataset)
      } else {
        title <- sprintf("Expression profile of gene %s", gene_id)
      }

      msg(".add_chunk: chunk title is", title)
      add_chunk(id=id, type="plot", code=code, label="gene_plot", rmd_var=rmd_var, title=title, 
                fig.cap=title, fig.width=fig_size$width, fig.height=fig_size$height)

}

.contrast_summary <- function(contrasts, dataset, gene_id, primary_id) {
  cn <- names(contrasts[[dataset]])
  res <- imap_dfr(contrasts, ~ {
                    .ds <- .y
                    imap_dfr(.x, ~ {
                               .x %>% filter(.data[[ primary_id ]] == gene_id) %>%
                                 mutate("Data set"=.ds, Contrast=.y)
                    })
              })

  res <- res %>% relocate(all_of(c("Data set", "Contrast")))
  numcol <- map_lgl(res, is.numeric)
  res %>% datatable(escape=FALSE, selection='none', options=list(pageLength=5)) %>%
    formatSignif(columns=colnames(res)[numcol], digits=2)
}

# returns the code to evaluate in order to generate the plot
# if substitute is TRUE, the variable names will be substituted
# by names of the variables in the environment which is parent to the shiny
# module
.gene_browser_generate_plot_code <- function(input, covar, exprs, id, dataset, mode, varnames, exprs_label=NULL, substitute=FALSE) {
  multi <- mode == "multi"

  ret <- ""
  covar_name <- "covar"
  exprs_name <- "exprs"

  if(is_na(input$covarXName)) { warning("X covariate is not defined") ; return(NULL) }
  if(is_na(input$covarYName)) { warning("Y covariate is not defined") ; return(NULL) }

  if(is.null(exprs_label)) { exprs_label <- "Expression" }

  if(substitute) {
    covar_name <- varnames$covar
    exprs_name <- varnames$exprs
  } else {
    # when running in own frame, multi is always true because we converted
    # the data
    multi <- TRUE
  }

  if(multi) {
    covar_name <- paste0(covar_name, '[["', dataset, '"]]')
    exprs_name <- paste0(exprs_name, '[["', dataset, '"]]')
  }

  # need to create this data frame to test whether variables are numerical
  # or categorical
  if(multi) {
    df <- data.frame(covar[[dataset]], Expression=exprs[[dataset]][id, , drop=TRUE])
    colnames(df)[ncol(df)] <- exprs_label
  } else {
    df <- data.frame(covar[[dataset]], Expression=exprs[[dataset]][id, , drop=TRUE])
    colnames(df)[ncol(df)] <- exprs_label
  }

  ret <- glue('df <- data.frame({covar_name}, Expression={exprs_name}["{id}", , drop=TRUE])')
  ret <- glue('{ret}\n\ncolnames(df)[ncol(df)] <- "{exprs_label}"')
  ret <- glue('{ret}\n\nggplot(df, aes(x={input$covarXName}, y={input$covarYName}')

  if(!is_na(input$colorBy))  { ret <- glue('{ret}, color={input$colorBy}') }
  if(!is_na(input$groupBy))  { ret <- glue('{ret}, group={input$groupBy}') }
  if(!is_na(input$symbolBy)) { ret <- glue('{ret}, shape={input$symbolBy}') }
  ret <- glue('{ret})) + ')

  if(!is.numeric(df[[input$covarXName]]) && is_na(input$groupBy)) {
    ret <- glue('
{ret}
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size=3, alpha=.5, width=.1)')
  } else {
    ret <- glue('
{ret}
  geom_point(size=3)')
  }

  if(!is_na(input$groupBy))   { ret <- glue('
{ret} + \n  geom_line()') }
  if(!is_na(input$trellisBy)) { ret <- glue('
{ret} + \n  facet_wrap( ~ {input$trellisBy})') }
  ret <- glue('
{ret} + \n  theme(text=element_text(size={input$fontSize}))')
  print(ret)
  return(ret)
}


## dynamically generate the UI controls for the gene plots from the
## actually selected data set
.dynamic_col_control <- function(id, covar, datasets, ds_selected) {

  covar_sum <- summary_colorDF(covar)
  all_covars         <- covar_sum %>% filter(unique > 1) %>% pull(.data$Col)
  non_unique         <- covar_sum %>% 
    filter(Class %in% c("<dbl>", "<int>") | unique < nrow(covar)) %>% 
    pull(.data$Col)
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
        fluidRow(tipify(selectInput(NS(id, "trellisBy"), "Trellis by", c("N/A", non_unique), selected="N/A", width="100%"),
                      "Each unique value of the variable will be shown on a separate subplot", placement="right")),
        fluidRow(tipify(numericInput(NS(id, "fontSize"),    label="Font size", min=6, value=14, step=1, width="50%"),
                      "Change the base font size of the figure", placement="right")),
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
                      tipify(downloadButton(NS(id, "save"), "PDF", class="bg-success"), "Save as PDF"),
                      tipify(actionButton(NS(id, "rmd"), "Add to report", class="bg-success"), "Rmarkdown")
                      ),
               column(width=11,
                      withSpinner(plotOutput(NS(id, "countsplot"), height="100%", width="100%")))
      )

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
#' @param exprs_label Label to be used for the expression values
#' @param rmd_var a reactive values object which will be used to store the
#'        generated markdown chunks
#' @return does not return anything useful
#' @importFrom shiny is.reactivevalues
#' @importFrom glue glue
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
                                  annot_linkout=NULL, rmd_var=NULL,
                                  exprs_label = "Expression") {
  parent_frame <- parent.frame()
  varnames <- list(
                   covar           = deparse(substitute(covar)),
                   exprs           = deparse(substitute(exprs)),
                   annot           = deparse(substitute(annot)),
                   exprs_label     = deparse(substitute(exprs_label)),
                   annot_linkout   = deparse(substitute(annot_linkout)),
                   description_col = deparse(substitute(description_col)),
                   symbol_col      = deparse(substitute(symbol_col)),
                   primary_id      = deparse(substitute(primary_id)),
                   cntr            = deparse(substitute(cntr))
                   )

  # if we have a single dataset, we need to wrap it into a list
  if(is.data.frame(covar)) {
    mode  <- "single"
    covar <- list(default=covar)
    exprs <- list(default=exprs)
    annot <- list(default=annot)
    cntr  <- list(default=cntr)

    if(!is.null(annot_linkout)) {
      annot_linkout <- list(default=annot_linkout)
    }
  } else {
    mode <- "multi"
    message("geneBrowserPlotServer in multi dataset mode")
  }

  # vector holding the names of all datasets
  datasets        <- names(covar)
  names(datasets) <- datasets

  # start the module server
  moduleServer(id, function(input, output, session) {
    disable("save")

    # if gene_id is not a reactiveValues object, wrap it into one
    if(!is.reactivevalues(gene_id)) {
      tmp <- gene_id
      gene_id <- reactiveValues()
      gene_id$id <- tmp$id
      if(is.null(gene_id$ds <- tmp$ds)) {
        gene_id$ds <- "default"
      }
    }

    # ds holds the dataset; g_id holds the gene ID
    ds        <- reactiveVal()
    g_id      <- reactiveVal()
    fig_size  <- reactiveValues(width=600, height=600)
    plot_code <- reactiveVal()

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


    observeEvent(input$figure_size, {
      fig_size$width <- 
        as.numeric(gsub(" *([0-9]+) *x *([0-9]+)", "\\1", input$figure_size))

      fig_size$height <- 
        as.numeric(gsub(" *([0-9]+) *x *([0-9]+)", "\\2", input$figure_size))
    })


    ## Save figure as a PDF
    output$save <- downloadHandler(
      filename = function() {
        .id <- g_id()
        .ds <- ds()
        ret <- sprintf("expression_profile_ds_%s_%s_covarX_%s_covarY_%s_colorBy_%s_groupBy_%s_symbolBy_%s_trellisBy_%s.pdf",
                       .ds, .id,
                       input$covarXName, input$covarYName, input$colorBy, input$groupBy, input$symbolBy, input$trellisBy)
        ret <- gsub("[^0-9a-zA-Z_.-]", "", ret)
        return(ret)
      },
      content = function(file) {
        pdf(file=file, width=8, height=5)
        code <- plot_code()
        print(eval(parse(text=code), envir=parent_frame))
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
      if(!isTruthy(ds()) || !isTruthy(g_id()) || is.null(cntr[[ ds() ]])) { return(NULL) }
      .contrast_summary(cntr, ds(), g_id(), primary_id)
    })
 
    ## Additional information - e.g. correlation coefficient if the
    ## covariate is numeric
    output$addInfo <- renderText({
      if(!isTruthy(ds()) || !isTruthy(g_id())) { return(NULL) }
      .gene_browser_info_tab(g_id(), covar[[g_id()]][[input$covarXName]], exprs[[ds()]][ g_id(), ])
    })
 
    ## reload the plot interface only if the data set (and covariates)
    ## changed
    output$col_control <- renderUI({
      .ds <- ds()
      if(!isTruthy(.ds)) { .ds <- 1 }
      .dynamic_col_control(id, covar[[.ds]], names(covar), datasets[.ds])
    })

    observe({
      if(!isTruthy(ds()) || !isTruthy(g_id())) { 
        msg("ds or g_id not truthy")
        return(NULL) 
      }
      if(!isTruthy(input$covarXName)) { 
        msg("covarXName not truthy")
        return(NULL) 
      }
      if(!isTruthy(input$covarYName)) { 
        msg("covarYName not truthy")
        return(NULL) 
      }
      if(is.na(g_id())) { 
        msg("g_id is NA")
        return(NULL) 
      }

      if(is.null(rmd_var)) { 
        msg("rmd_var is NULL")
        return(NULL) 
      }

      message(sprintf("generate code started with dataset=%s and gene=%s", ds(), g_id()))
      plot_code(.gene_browser_generate_plot_code(input=input, covar=covar, exprs=exprs,
                                      id=g_id(), dataset=ds(), 
                                      mode=mode, varnames=varnames, 
                                      exprs_label=exprs_label, substitute=TRUE))
    })

 
    # generate the markdown code
    observeEvent(input$rmd, {
      msg("rmd pressed")
      .add_chunk(id, plot_code(), rmd_var, g_id(), ds(), mode, fig_size)
    })
 
    ## The actual plot. need to put inside "observe" to use the reactive
    ## figure size
    observe({ output$countsplot <- renderPlot({
      code <- plot_code()
      if(!isTruthy(code)) { return(NULL) }
      enable("save")
      
      message(sprintf("plotting started with dataset=%s and gene=%s", ds(), g_id()))
      message("evaluating code:", code)
      eval(str2expression(code), envir=parent_frame)
      }, width=fig_size$width, height=fig_size$height) 
    })

  }) # end of moduleServer
}

