## prepare contrasts for the gene browser, adding action button, sorting,
## removing unnecessary columns etc.
.gene_browser_prep_res <- function(cntr, but, annot, annot_linkout, primary_id, cols_to_hide) {

	names(cntr_ids) <- cntr_ids <- names(cntr)

  cntr <- map(cntr_ids, ~
              .gene_browser_prep_res_single(
                                            .x,
                                            cntr[[.x]],
                                            but,
                                            annot[[.x]],
                                            annot_linkout[[.x]],
                                            primary_id,
                                            cols_to_hide))
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
.gene_browser_prep_res_single <- function(ds_id, cntr, but, annot, annot_linkout, primary_id,
                                          cols_to_hide) {

  cntr   <- cntr %>%
    map(~ { .x %>% select(all_of(setdiff(colnames(.x), cols_to_hide))) }) %>%
    map(~ { .x %>% { merge(annot, ., by=primary_id, all.x=TRUE) }} )

  if(!is.null(but) && length(but) == 1L) {
    cntr   <- cntr %>% 
      map(~ { .x %>% mutate('>'= sprintf(but, ds_id, .data[[primary_id]])) }) %>% 
      map(~ { .x %>% relocate(all_of(">"), .before=1) %>% arrange(pvalue)})
  }

  if(!is.null(annot_linkout)) {
    cntr <- map(cntr, ~ .apply_annot_linkout(.x, annot_linkout))
  }

  return(cntr)
}



.gene_browser_prep_res_single_contrast <- function(ds_id, cntr_id, cntr, but, annot, annot_linkout, primary_id, cols_to_hide) {

  res <- cntr[[ ds_id ]][[ cntr_id ]]

  res <- res %>% select(all_of(setdiff(colnames(res), cols_to_hide))) %>% 
    { merge(annot[[ ds_id ]], ., by=primary_id, all.x=TRUE) } %>% arrange(pvalue)

  if(!is.null(but) && length(but) == 1L) {
    res <- res %>%
      mutate('>'= sprintf(but, ds_id, .data[[primary_id]])) %>% 
      relocate(all_of(">"), .before=1)

  }
  

  if(!is.null(annot_linkout)) {
    res <-  .apply_annot_linkout(res, annot_linkout)
  }

  res
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
    if(is.null(names(cntr_titles))) {
      names(cntr_titles) <- cntr_titles
    }

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
      withSpinner(DTOutput(NS(id, "result_tbl"))),
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
        stopifnot(all(map_lgl(cntr[[c]], ~ is.data.frame(.x) || is(.x, "disk.frame"))))
      }
    }
  }

  if(!is.null(annot)) {
    if(multilevel) {
      stopifnot(all(map_lgl(annot, ~ is.data.frame(.x) || is(.x, "disk.frame"))))
      if(primary_id == 0) {
        stopifnot(all(map_lgl(annot, ~ !is.null(rownames(.x)))))
      } else {
        stopifnot(all(map_lgl(annot, ~ primary_id %in% names(.x))))
      }
    } else {
      stopifnot(is.data.frame(annot) || is(.x, "disk.frame"))
      if(primary_id == 0) {
        stopifnot(!is.null(rownames(annot)))
      } else {
        stopifnot(primary_id %in% names(annot))
      }
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
#' @param gene_id must be either NULL or a `reactiveValues` object. If not NULL, then
#' a button is displayed; clicking on it will modify the `gene_id` reactive
#' value (possibly triggering an event in another module).
#' @param cols_to_hide columns in the contrasts data frames to hide in the table. 
#'        Default is for DESeq2 derived contrasts.
#' @return `geneBrowserTableServer` returns NULL. `geneBrowserTableUI`
#'         returns the interface.
#' @importFrom rlang .data
#' @importFrom stats cor.test
#' @importFrom bslib bs_theme
#' @examples 
#' if(interactive()) {
#'   data(C19)
#'   ui <- fluidPage(
#'            geneBrowserTableUI("gb", names(C19$contrasts))
#'         )
#'
#'   server <- function(input, output) {
#'     geneBrowserTableServer("gb", cntr=C19$contrasts,
#'                                  annot=C19$annotation)
#'   }
#'
#'   shinyApp(ui, server)
#' }
#' @export
geneBrowserTableServer <- function(id, cntr, annot, annot_linkout=NULL,
                                   primary_id="PrimaryID", gene_id=NULL,
                                   cols_to_hide=c("stat", "lfcSE", "symbol", "entrez")
                                   ) {

  multilevel <- .check_multilevel(cntr)
  .check_params(multilevel, cntr=cntr, annot=annot, 
                annot_linkout=annot_linkout, primary_id=primary_id)

  if(!multilevel) {
    cntr <- list(default=cntr)
    annot <- list(default=annot)
    annot_linkout=list(default=annot_linkout)
  }

  ## in this case we take primary id from row names of annot / cntr
  if(primary_id == 0) {
    primary_id <- "__primary_id"
    annot <- map(annot, ~ {
                   stopifnot(!is.null(rownames(.x)))
                   .x[[primary_id]] <- rownames(.x) 
                   .x })
    cntr  <- map(cntr, ~ 
                 map(.x, ~ {
                       stopifnot(!is.null(rownames(.x)))
                       .x[[primary_id]] <- rownames(.x)
                       .x }))
  }


  moduleServer(id, function(input, output, session) {

    if(is.null(gene_id)) {
      but <- NULL
    } else {
      but <- actionButton("go~%s~%s", label=" \U25B6 ", 
                         onclick=sprintf('Shiny.onInputChange(\"%s-select_button\",  this.id)', id),  
                         class = "btn-primary btn-sm")
    }

   #cntr <- .gene_browser_prep_res(cntr, as.character(but), annot, 
   #                               annot_linkout, primary_id=primary_id,
   #                               cols_to_hide=cols_to_hide)

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

    output$result_tbl <- renderDT({

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

      res <- .gene_browser_prep_res_single_contrast(.ds, .cntr, cntr, 
                                                    as.character(but), 
                                                    annot, annot_linkout, primary_id, 
                                                    cols_to_hide)

      res %>% datatable(escape=FALSE, selection='none', extensions="Buttons",
                options=list(pageLength=5, dom="Bfrtip", scrollX=TRUE)) %>%
        formatSignif(columns=intersect(colnames(res), 
                                       c("baseMean", "log2FoldChange", "pvalue", "padj")), digits=2)
    })
    return(NULL)
  })
}


#' Launch a browser of DE analysis results
#'
#' Launch a simple shiny-based browser of DE analysis results
#'
#' Creates a shiny app, which allows to show gene expression
#' profiles.
#'
#' @param x object holding the DE analysis results. List containing the
#'          elements `contrasts`, `covariates`, `expression` and `annotation`.
#' @param primary_id The name of the column in `contrasts` and `annotation`
#'                   which holds the primary, unique identifier of genes.
#' @return does not return a value
#' @importFrom rlang .data
#' @importFrom stats cor.test
#' @importFrom bslib bs_theme
#' @inheritParams geneBrowserTableServer
#' @examples
#' if(interactive()) {
#'   data(C19)
#'   annot_linkout <- list(
#'     SYMBOL="https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s",
#'     ENTREZ="https://www.ncbi.nlm.nih.gov/gene/?term=%s"
#'   )
#'   gene_browser(C19, annot_linkout=annot_linkout)
#' }
#' @export
gene_browser <- function(x, annot_linkout=NULL, 
                            cols_to_hide=c("stat", "lfcSE", "symbol", "entrez"),
                            primary_id="PrimaryID") {

  thematic_shiny(font="auto")

  ## prepare the UI
  ui <- fluidPage(
    useShinyjs(),
    #theme = bs_theme(bootswatch = "united"),
    fluidRow(titlePanel(h1("Gene browser")), class="bg-primary"),
    fluidRow(HTML("<hr>")),
    geneBrowserTableUI("geneTab", names(x$contrasts)),
    geneBrowserPlotUI("genePlot", contrasts=TRUE)
  )

  ## prepare the server
  server <- function(input, output, session) {
    gene_id <- reactiveValues()
    geneBrowserTableServer("geneTab", 
                           cntr=x$contrasts, 
                           annot=x$annotation, 
                           gene_id=gene_id,
                           annot_linkout=annot_linkout)
    geneBrowserPlotServer("genePlot", 
                          gene_id, 
                          covar=x$covariates, 
                          exprs=x$expression, 
                          annot=x$annotation)
  }

  shinyApp(ui, server)
}

