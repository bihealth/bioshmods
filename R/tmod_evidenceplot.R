
.tmod_rev_db_map_ids <- function(ids, dbname, tmod_dbs_mapping_obj) {

  mp_id   <- tmod_dbs_mapping_obj$dbs[dbname]
  mapping <- tmod_dbs_mapping_obj$maps[[mp_id]]

  ret <- names(mapping)[ match(ids, mapping) ]
  
  if(all(is.na(ret))) {
    warning(glue("No IDs were found in the mapping... are you sure you are using the {dbname} IDs?\nCheck the annotation data frame (see `get_annot()`)"))
  }

  names(ret) <- ids

  return(ret)
}


## create a datatable with the genes from a gene set
.tmod_browser_gene_table <- function(but, ds, id, db_name, cntr_name, 
                                     sort_name, tmod_dbs, cntr, tmod_map, primary_id) {

  db        <- tmod_dbs[[db_name]]
  genes     <- db[["MODULES2GENES"]][[id]]
  genes_pid <- .tmod_rev_db_map_ids(ids=genes, dbname=db_name, tmod_dbs_mapping_obj=tmod_map)

  ret <- cntr[[cntr_name]] %>% filter(.data[[primary_id]] %in% genes_pid)
  ret <- ret %>% mutate('>' = sprintf(but, ds, .data[[primary_id]])) %>% relocate(all_of(">"), .before=1) %>%
    arrange(.data[["pvalue"]])

  num_cols <- colnames(ret)[ map_lgl(ret, is.numeric) ]

  datatable(ret, escape=FALSE, selection='none',
                options=list(pageLength=5, dom="Bfrtip", scrollX=TRUE, buttons=c("copy", "csv", "excel"))) %>%
        formatSignif(columns=num_cols, digits=2)
  
}

.plot_evidence <- function(mod_id, cntr_id, db_id, sort_id, cntr, tmod_dbs,
                           tmod_gl=NULL, tmod_map=NULL, annot=NULL, primary_id) {

  mset <- tmod_dbs[[db_id]]
  cntr <- cntr[[ cntr_id ]]

  if(!all(c("pvalue", "log2FoldChange", "padj") %in% colnames(cntr))) {
    stop("Parameter `cntr` must have columns log2FoldChange, pvalue and padj")
  }

  if(!primary_id %in% colnames(cntr)) {
    cntr[[primary_id]] <- rownames(cntr)
  }

  if(is.null(tmod_gl)) {
    ## create an ad hoc list
    cntr <- cntr[ order(cntr$pvalue), ]
    gl <- tmod_map$maps[[ tmod_map$dbs[[ db_id ]] ]][ cntr[[primary_id]] ]
  } else {
    gl <- tmod_gl[[cntr_id]][[db_id]][[sort_id]]
  }

  symbols <- names(gl)

  if("SYMBOL" %in% colnames(annot)) {
    symbols <- annot[["SYMBOL"]][ match(symbols, annot[[primary_id]]) ]
  }

  names(symbols) <- names(gl)

  ## We need to figure out what genes are in the gene set to select them as
  ## labels

  genes <- getModuleMembers(mod_id, mset)[[mod_id]]
  if(is.null(genes)) {
    stop(sprintf("Gene set %s not found in db %s", mod_id, db_id))
  }

  sel                <- gl %in% genes
  gene.labels        <- symbols[sel]
  names(gene.labels) <- gl[sel]
  
  ## now for some color
  # make sure that cntr is in the same order as gl
  cntr <- cntr[ match(names(gl), cntr[[primary_id]]), ]

  p <- cntr$padj
  l <- cntr$log2FoldChange

  colors <- ifelse(l < 0, 
      ifelse(p < .05, 'blue', '#000066'),
      ifelse(p < .05, 'red', '#660000'))
  names(colors) <- gl

  evidencePlot(gl, m=mod_id, mset=mset, gene.colors=colors, gene.labels=gene.labels)
}

## given a module, contrast, sorting: prepare a module info tab contents,
## including which genes are significant in the given contrast
.tmod_browser_mod_info <- function(id, ds, db_name, cntr_name, sort_name, tmod_dbs, cntr, tmod_map) {
  db <- tmod_dbs[[db_name]]
  ret <- sprintf("Module ID: %s\nDescription: %s\nData set: %s\nContrast: %s\nDatabase: %s\nSorting: %s",
          id,
          db[["MODULES"]][id, ][["Title"]],
          ds,
          cntr_name,
          db_name, 
          sort_name)

}


#' @rdname tmodBrowserPlotServer
#' @export
tmodBrowserPlotUI <- function(id) {
    sidebarLayout(
      sidebarPanel(
        fluidRow(downloadButton(NS(id, "save"), "Save plot", class="bg-success")),
        fluidRow(verbatimTextOutput(NS(id, "modinfo"))),
        width=5
      ),
      mainPanel(
        fluidRow(verbatimTextOutput(NS(id, "cmdline"))),
        fluidRow(
                 tabsetPanel(
                             tabPanel("Plot", withSpinner(plotOutput(NS(id, "evidencePlot"), height="100%"))),
                             tabPanel("Genes", withSpinner(DTOutput(NS(id, "moduleGenes"))))
                 )),
        width=7
      )
    )
}



## server module for viewing evidence plots

#' Shiny Module – tmod browser evidence plots
#'
#' Shiny Module – gene browser evidence plots
#'
#' This part is a bit complex, because a lot of different data go into
#' the evidence plots. To create a plot, following elements are necessary:
#'
#'  * ordered gene list: this is the same gene list that is used as an
#'    input to tmod
#'  * a list of gene set collections (tmod gene set databases) for which the enrichments have been run
#'  * contrast data frame – to know which genes go up or down, for displaying gene names in color
#'  * if no gene list is given, then using the contrast data frame it is
#'    possible to create a list ordered by p-values. However, since the
#'    gene set database might use a different type of identifiers than the
#'    PrimaryID column of the contrasts data frame, it is necessary to
#'    provide a mapping between the PrimaryIDs and the database gene IDs as well.
#'
#' The sections below discuss these elements.
#'
#'
#' @section tmod gene lists:
#' To display an evidence plot, we need to have an
#' ordered list of genes. This has to be provided from outside, as many
#' different sortings are possible. The parameter tmod_gl is a hierarchical
#' list:
#'  * First level: contrasts
#'  * Second level: gene set databases
#'  * Third level: sorting type
#' 
#' For example, then the `tmod_gl[["Contrast1"]][["tmod"]][["pval"]]` is a
#' vector of gene identifiers which come from the tmod database (in this
#' case, gene names). The vector is named and the names are primary IDs
#' matching the contrast data frames.
#'
#' If this argument is NULL, then the genes will be ordered by p-values
#' from the contrast object provided. However, in this case it is necessary
#' to provide a mapping between the PrimaryIDs of the contrasts and the
#' gene identifiers used by the gene set database.
#'
#' @inheritSection tmodBrowserTableServer Use of tmod database objects
#'
#' @param id ID of the shiny module
#' @param tmod_dbs tmod gene set databases returned by `get_tmod_dbs()`
#' @param tmod_map tmod gene set ID mapping returned by `get_tmod_mapping()`
#' @param tmod_gl tmod gene lists. See details.
#' @param id identifier (same as the one passed to geneBrowserTableUI)
#' @param primary_id name of the column which holds the primary identifiers
#' @param cntr list of contrast results returned by `get_contrasts()`
#' @param annot data frame containing gene annotation 
#' @param gene_id must be a `reactiveValues` object. If not NULL, then
#' clicking on a gene identifier will modify this object (possibly
#' triggering an event in another module).
#' @param gs_id a "reactive values" object (returned by `reactiveValues()`), including 
#' dataset (`ds`), gene set ID (`id`), contrast id (`cntr`), database ID
#' (`db`) and sorting mode (`sort`). If `mod_id` is not `NULL`, these
#' reactive values will be populated, possibly triggering an action in
#' another shiny module.
#' @return returns a reactive value with a selected gene identifier
#' @importFrom shinyBS popify
#' @importFrom tmod evidencePlot getModuleMembers tmodDecideTests tmodSummary tmodPanelPlot
#' @importFrom glue glue
#' @example inst/examples/tmod_browser.R
#' @export
tmodBrowserPlotServer <- function(id, gs_id, tmod_dbs, cntr, tmod_map=NULL, tmod_gl=NULL, annot=NULL, 
                                  primary_id="PrimaryID", gene_id=NULL) {

  stopifnot(!is.null(tmod_gl) || !is.null(tmod_map))

  if(!is.data.frame(annot)) {
    message("tmodBrowserPlotServer: running in multilevel mode")
  } else {
    tmod_dbs <- list(default=tmod_dbs)
    cntr     <- list(default=cntr)
    tmod_map <- list(default=tmod_map)
    tmod_gl  <- list(default=tmod_gl)
    annot    <- list(default=annot)
  }
    
  moduleServer(id, function(input, output, session) {
    message("Launching tmod plot server")

    gene.but <- actionButton("go~%s~%s", label=" \U25B6 ", 
                        onclick=sprintf('Shiny.onInputChange(\"%s-gene_select_button\",  this.id)', id),  
                        class = "btn-primary btn-sm")

    observeEvent(input$gene_select_button, {
      if(!is.null(gene_id)) {
        ids <- strsplit(input$gene_select_button, '~')[[1]]
        gene_id$ds <- ids[2]
        gene_id$id <- ids[3]
      }

    })

    disable("save")

    ## create the evidence plot and display the command line to replicate it
    output$evidencePlot <- renderPlot({
      #message("Rendering plot")
      req(gs_id$id)

      if(is.null(gs_id$id)) { return(NULL) }
      enable("save")

      ds <- gs_id$ds
      .plot_evidence(mod_id=gs_id$id, cntr_id=gs_id$cntr, db_id=gs_id$db, sort_id=gs_id$sort, 
                     cntr=cntr[[ds]], tmod_dbs=tmod_dbs[[ds]], tmod_gl=tmod_gl[[ds]], 
                     tmod_map=tmod_map[[ds]], annot=annot[[ds]], primary_id)
    }, width=800, height=600)


    output$modinfo <- renderText({
      req(gs_id$id)
      if(!isTruthy(gs_id$id)) { return(NULL) }
      ret <- .tmod_browser_mod_info(gs_id$id, gs_id$ds, gs_id$db, gs_id$cntr, gs_id$sort, 
                                    tmod_dbs[[gs_id$ds]], cntr[[gs_id$ds]], tmod_map[[gs_id$ds]])
      return(ret)
    })

    output$moduleGenes <- renderDT({
      req(gs_id$id)
      if(!isTruthy(gs_id$id)) { return(NULL) }
      ds <- gs_id$ds
      .tmod_browser_gene_table(as.character(gene.but), ds,
                                    gs_id$id, gs_id$db, gs_id$cntr, gs_id$sort, 
                                    tmod_dbs[[ds]], cntr[[ds]], tmod_map[[ds]], primary_id)
    })
    


    ## save the plot as PDF
    output$save <- downloadHandler(
      filename = function() {
        req(gs_id$id)
        if(!isTruthy(gs_id$id)) { return(NULL) }
        ret <- sprintf("evidence_plot_%s_%s_%s_%s.pdf", gs_id$ds, gs_id$db, gs_id$cntr, gs_id$id)
        return(ret)
      },
      content = function(file) {
        req(gs_id$id)
        if(!isTruthy(gs_id$id)) { return(NULL) }
        pdf(file=file, width=8, height=5)
        title <- sprintf("%s / %s\nContrast: %s / %s", 
                         gs_id$id, gs_id$db, gs_id$cntr, gs_id$sort)

        ds <- gs_id$ds
        .plot_evidence(mod_id=gs_id$id, cntr_id=gs_id$cntr, db_id=gs_id$db, sort_id=gs_id$sort, 
                     cntr=cntr[[ds]], tmod_dbs=tmod_dbs[[ds]], tmod_gl=tmod_gl[[ds]], 
                     tmod_map=tmod_map[[ds]], annot=annot[[ds]], primary_id)
        dev.off()
      }
    )
  })
}

