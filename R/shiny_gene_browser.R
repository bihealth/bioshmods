## prepare contrasts for the gene browser, adding action button, sorting,
## removing unnecessary columns etc.
.gene_browser_prep_res <- function(cntr, but, annot, annot_linkout, primary_id) {

	names(cntr_ids) <- cntr_ids <- names(cntr)

  cntr <- map(cntr_ids, ~
              .gene_browser_prep_res_single(
                                            .x,
                                            cntr[[.x]],
                                            but,
                                            annot[[.x]],
                                            annot_linkout[[.x]],
                                            primary_id))
  return(cntr) 
}

.apply_annot_linkout <- function(.x, annot_linkout) {

  for(n in intersect(colnames(.x), names(annot_linkout))) {
    fmt <- sprintf('<a href="%s" target="_blank">%%s</a>', annot_linkout[[n]])
    .x[[n]] <- ifelse(
                      is.na(.x[[n]]) | .x[[n]] == "", 
                      .x[[n]], 
                      sprintf(fmt, .x[[n]], .x[[n]]))
  }

  return(.x)
}

## this function runs the preparation for a single dataset
.gene_browser_prep_res_single <- function(ds_id, cntr, but, annot, annot_linkout, primary_id) {
  cntr   <- cntr %>% 
    map(~ { .x %>% mutate('>'= sprintf(but, ds_id, .data[[primary_id]])) }) %>%
    map(~ { .x %>% select(all_of(setdiff(colnames(.x), c("stat", "lfcSE", "symbol", "entrez")))) }) %>%
    map(~ { .x %>% { merge(annot, ., by=primary_id, all.x=TRUE) } %>% 
        relocate(all_of(">"), .before=1) %>% arrange(pvalue)})

  if(!is.null(annot_linkout)) {
    cntr <- map(cntr, ~ .apply_annot_linkout(.x, annot_linkout))
  }

  return(cntr)
}






## prepares IDs/titles of the contrasts for use with the UI
.prep_cntr_titles <- function(cntr_titles) {

  if(is.list(cntr_titles)) {
    cntr_titles <- imap(cntr_titles, ~ {
                        ret <- paste0(.y, '::', .x)
                        #names(ret) <- paste0(.y, ": ", names(.x))
                        names(ret) <- names(.x)
                        ret
         })
    if(length(cntr_titles) < 2) {
      cntr_titles <- cntr_titles[[1]]
    }
  } else {
    tmp <- cntr_titles
    cntr_titles <- paste0("default::", tmp)
    names(cntr_titles) <- names(tmp)
  }

  return(cntr_titles)
}

#' @rdname geneBrowserTableServer
#' @export
geneBrowserTableUI <- function(id, cntr_titles) {

  cntr_titles <- .prep_cntr_titles(cntr_titles)
  but <- actionButton("foo", label=" \U25B6 ", class = "btn-primary btn-sm")
  sidebarLayout(
    sidebarPanel(
        fluidRow(
                 tipify(selectInput(NS(id, "contrast"), label = "Contrast", choices = cntr_titles, width="100%"),
                        "Choose the contrast to show in the table", placement="right")
                        ),
        fluidRow(
                 column(6,
                 tipify(checkboxInput(NS(id, "filter"), "Filter", TRUE, width="100%"),
                        "Choose whether the output should be filtered"),
                 tipify(numericInput(NS(id, "f_lfc"),    label="Filter by abs(LFC)", min=0, value=0.5, step=.1, width="100%"),
                        "Show only genes which have an absolute log fold change greater than this value"),
                 tipify(numericInput(NS(id, "f_pval"),    label="Filter by FDR", min=0, max=1.0, value=0.05, step=.1, width="100%"),
                        "Show only genes which have a p-value smaller than this value")
                        ), 

                 column(6,
                 tipify(
                 selectInput(NS(id, "f_dir"), label="Direction", choices=c(Any="any", Up="up", "Down"="dw"), 
                             width="100%"),
                        "Show only genes with log fold change greater or smaller than 0", placement="right")
                 )),
      tagList(
        HTML(paste("Click on the", but, "buttons to view an expression profile<br/>"))
      ),
      width=3
    ),
    mainPanel(
      withSpinner(dataTableOutput(NS(id, "result_tbl"))),
      width=9
    )
  )
}

## check whether data is list of lists or just a list of dfs
.check_multilevel <- function(cntr) {
  !any(map_lgl(cntr, is.data.frame)) 
}

.check_params <- function(multilevel, annot=NULL, cntr=NULL, annot_linkout=NULL, primary_id="PrimaryID") {

  if(!is.null(cntr)) {
    stopifnot(!is.null(names(cntr)))

    if(multilevel) {
      for(c in names(cntr)) {
        stopifnot(all(map_lgl(cntr[[c]], is.data.frame)))
      }
    }
  }

  if(!is.null(annot)) {
    if(multilevel) {
      stopifnot(all(map_lgl(annot, is.data.frame)))
      stopifnot(all(map_lgl(annot, ~ primary_id %in% colnames(.x))))
    } else {
      stopifnot(is.data.frame(annot))
      stopifnot(primary_id %in% colnames(annot))
    }
  }

  if(!is.null(annot) && !is.null(cntr)) {
    if(multilevel) {
      stopifnot(!is.null(names(annot)))
      stopifnot(all(names(cntr) %in% names(annot)))
    }
  }
}


#' Shiny Module – gene browser table selection
#'
#' Shiny Module – gene browser table selection
#'
#' The basic data set structure that this module takes is a named list of data
#' frames. These data frames will be shown in the browser when the specific
#' contrast (corresponding to a name in the list) is selected from the
#' configuration sidebar. The data frames *must* contain a column called
#' "PrimaryID" (this identifier can be changed with the parameter
#' `primary_id`). This is necessary in order to link the table rows with
#' e.g. plotting genes with `geneBrowserPlotServer`.
#'
#' Log2 fold changes must be stored in a column called
#' "log2FoldChange", and p-values in a column called "padj". These are the
#' default column names returned by DESeq2.
#'
#' Alternatively, tmodBrowserTableServer takes a list of lists of data
#' frames; that is, it allows to group the results of differential gene
#' analysis.
#'
#' The linkout feature (parameter `annot_linkout`) allows to define how the
#' different columns from the annotation data frame are represented as
#' linkouts to external data bases. 
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
#' @param cntr a list of data frames containing the DE analysis results, or
#'             a list of lists of data frames
#' @param annot annotation data frame containing column 'PrimaryID' (
#'        or another specified by the parameter `primary_id`)
#'        corresponding to the rownames of the contrast data frames
#' @param id identifier for the namespace of the module
#' @param primary_id name of the column which holds the primary identifiers
#' @param cntr_titles named character vector for contrast choices
#' @param annot_linkout a list; see Details. 
#' @param gene_id must be a `reactiveValues` object. If not NULL, then
#' clicking on a gene identifier will modify this object (possibly
#' triggering an event in another module).
#' @return reactive value containing the gene ID
#' @importFrom rlang .data
#' @importFrom stats cor.test
#' @importFrom bslib bs_theme
#' @export
geneBrowserTableServer <- function(id, cntr, annot, annot_linkout=NULL,
                                   primary_id="PrimaryID", gene_id=NULL) {

  multilevel <- .check_multilevel(cntr)
  .check_params(multilevel, cntr=cntr, annot=annot, 
                annot_linkout=annot_linkout, primary_id=primary_id)

  if(!multilevel) {
    cntr <- list(default=cntr)
    annot <- list(default=annot)
    annot_linkout=list(default=annot_linkout)
  }

  moduleServer(id, function(input, output, session) {

    but <- actionButton("go~%s~%s", label=" \U25B6 ", 
                         onclick=sprintf('Shiny.onInputChange(\"%s-select_button\",  this.id)', id),  
                         class = "btn-primary btn-sm")

    cntr <- .gene_browser_prep_res(cntr, as.character(but), annot, 
                                   annot_linkout, primary_id=primary_id)

    observeEvent(input$select_button, {
      if(!is.null(gene_id)) {
        ids <- strsplit(input$select_button, '~')[[1]]
        gene_id$ds <- ids[2]
        gene_id$id <- ids[3]
      }
    })

    observeEvent(input$filter, {
                   if(length(input$filter) > 0 && input$filter) {
                     enable("f_dir")
                     enable("f_pval")
                     enable("f_lfc")
                   } else {
                     disable("f_dir")
                     disable("f_pval")
                     disable("f_lfc")
                   }
    })

    output$id_summary <- renderText({
      .cntr <- input$contrast
      sprintf("Contrast is %s", .cntr)
    })

    output$result_tbl <- DT::renderDataTable({

      .cntr <- gsub(".*::", "", input$contrast)
      .ds   <- gsub("::.*", "", input$contrast)

      res <- cntr[[ .ds ]][[ .cntr ]]
      if(input$filter) {
        if(input$f_dir == "up") {
          res <- res %>% filter(.data[["log2FoldChange"]] > 0)
        } else if(input$f_dir == "dw") {
          res <- res %>% filter(.data[["log2FoldChange"]] < 0)
        }

        res <- res %>% filter(.data[["padj"]] < input$f_pval & abs(.data[["log2FoldChange"]]) > input$f_lfc) 
      }

      res %>% datatable(escape=FALSE, selection='none', extensions="Buttons",
                options=list(pageLength=5, dom="Bfrtip", scrollX=TRUE)) %>%
        formatSignif(columns=intersect(colnames(res), 
                                       c("baseMean", "log2FoldChange", "pvalue", "padj")), digits=2)
    })
  })
}


##' Launch a browser of DE analysis results
##'
##' Launch a shiny-based browser of DE analysis results
##'
##' Launches a shiny app, web based, which allows to show gene expression
##' profiles in a pipeline. 
##'
##' To speed up launching, you can pre-load the contrasts with the
##' `get_contrasts` function.
##' @param pip pipeline object returned by `load_de_pipeline`
##' @param cntr (optional) pre-loaded contrasts (returned by `get_contrasts`)
##' @param annot (optional) pre-loaded annotation table (returned by `get_annot`)
##' @return does not return a value
##' @importFrom rlang .data
##' @importFrom stats cor.test
##' @importFrom bslib bs_theme
##' @examples
##' \dontrun{
##' pip <- load_de_pipeline(config_file="DE_config.yaml")
##' gene_browser(pip, tmod_dbs)
##' }
##' @export
# gene_browser <- function(pip, cntr=NULL, annot=NULL) {
#
#   message("preparing...")
#   if(is.null(annot)) {
#     message(" * Loading Annotation (consider using the annot option to speed this up)")
#     annot  <- get_annot(pip)
#   }
#
#   # prepare the contrast tables
#   if(is.null(cntr)) {
#     message(" * Loading contrasts... (consider using the cntr option to speed this up)")
#     cntr <- get_contrasts(pip)
#   }
#
#   config <- get_config(pip)
#   covar  <- get_covariates(pip)
#
#   rld    <- get_object(pip, step="DESeq2", extension="rld.blind.rds")
#   rld    <- rld@assays@data@listData[[1]]
#   
#   cntr_titles <- map_chr(config$contrasts$contrast_list, `[[`, "ID")
#   names(cntr_titles) <- map_chr(config$contrasts$contrast_list, `[[`, "title")
#
#   thematic_shiny(font="auto")
#
#   ## prepare the UI
#   ui <- fluidPage(
#     useShinyjs(),
#     theme = bs_theme(bootswatch = "united"),
#     fluidRow(titlePanel(h1("Gene browser")), class="bg-primary"),
#     fluidRow(HTML("<hr>")),
#     geneBrowserTableUI("geneTab", cntr_titles),
#     geneBrowserPlotUI("genePlot", contrasts=TRUE)
#   )
#
#   ## prepare the server
#   server <- function(input, output, session) {
#     gene_id <- reactiveValues()
#     geneBrowserTableServer("geneTab", cntr, annot, gene_id=gene_id)
#     geneBrowserPlotServer("genePlot", gene_id, covar, rld, annot)
#
#     message("launching!")
#   }
#
#   shinyApp(ui, server)
# }

