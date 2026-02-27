# Normalize covariates to dataset-keyed list form.
# Supports single-data-frame and multi-dataset list inputs.
.heatmap_normalize_covar <- function(covar, datasets) {
  if(is.null(covar)) {
    return(stats::setNames(vector("list", length(datasets)), datasets))
  }

  if(is.data.frame(covar)) {
    if(length(datasets) != 1L) {
      stop("If multiple datasets are present, `covar` must be a named list.")
    }
    return(stats::setNames(list(covar), datasets))
  }

  if(!is.list(covar) || length(covar) == 0L || !all(vapply(covar, is.data.frame, logical(1)))) {
    stop("`covar` must be a data frame or a named list of data frames.")
  }

  if(is.null(names(covar))) {
    names(covar) <- paste0("dataset_", seq_along(covar))
  }
  names(covar) <- as.character(names(covar))
  missing_names <- is.na(names(covar)) | trimws(names(covar)) == ""
  names(covar)[missing_names] <- paste0("dataset_", seq_along(covar))[missing_names]

  if(length(datasets) == 1L && length(covar) == 1L && !datasets[1] %in% names(covar)) {
    names(covar) <- datasets
  }

  missing_ds <- setdiff(datasets, names(covar))
  if(length(missing_ds) > 0L) {
    stop(sprintf("`covar` is missing dataset(s): %s", paste(missing_ds, collapse=", ")))
  }

  covar[datasets]
}

# Draw a ComplexHeatmap object while controlling legend visibility.
.heatmap_draw <- function(hm, show_legend=TRUE) {
  ComplexHeatmap::draw(
    hm,
    show_heatmap_legend=isTRUE(show_legend),
    show_annotation_legend=isTRUE(show_legend)
  )
}

# List covariate columns that can be shown as heatmap annotations.
# Excludes the sample ID column used for joining.
.heatmap_annotation_choices <- function(covar_ds, sample_id_col="SampleID") {
  if(is.null(covar_ds) || !is.data.frame(covar_ds) || ncol(covar_ds) < 1L) {
    return(character(0))
  }

  setdiff(colnames(covar_ds), sample_id_col)
}

# List annotation columns that can be used for heatmap row labels.
# Excludes the identifier column used to match expression row names.
.heatmap_row_label_choices <- function(annot_ds, primary_id="PrimaryID") {
  if(is.null(annot_ds) || !is.data.frame(annot_ds) || ncol(annot_ds) < 1L) {
    return(character(0))
  }

  setdiff(colnames(annot_ds), primary_id)
}

#' Heatmap Module UI
#'
#' @rdname heatmapServer
#' @export
heatmapUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      shiny::tags$div(
        class="well",
        h4("Gene selection options"),
      geneGroupSelectorUI(ns("gene_selector")),
      shiny::tags$div(
        style="margin-top:8px;",
        textOutput(ns("selected_genes_n"), inline=TRUE),
        shiny::tags$span(
          style="margin-left:8px;color:#a94442;font-weight:600;",
          textOutput(ns("selected_genes_warning"), inline=TRUE)
        )
      )
      ),
      shiny::tags$div(
        class="well",
        h4("Heatmap options"),
        gridLayout(
          selectizeInput(
            ns("sel_annot"),
            "Covariates to indicate",
            choices=character(0),
            selected=character(0),
            multiple=TRUE
          ),
          selectizeInput(
            ns("annot_row_col"),
            "Row label column",
            choices=character(0),
            selected="",
            multiple=FALSE
          ),
          figsizeInput(ns("figure_size"), width="100%", selected="800x600"),
          checkboxInput(ns("show_legend"), "Show legend", value=TRUE),
          downloadButton(ns("save"), "Save heatmap to PDF", class="bg-success"),
          .ncol = 2,
          .nrow = 3
        ),
        fluidRow(column(colorPalettesUI(ns("heatmap_color")), width=12))
      ),
      width=3
    ),
    mainPanel(
      plotOutput(ns("heatmap_plot")),
      width=9
    )
  )
}

#' Shiny module for expression heatmaps
#'
#' This module embeds [geneGroupSelectorUI()] / [geneGroupSelectorServer()]
#' on the left side and renders a heatmap for selected genes on the right side.
#'
#' @param id Module identifier (same as passed to [heatmapUI()]).
#' @param annot Annotation data frame (single dataset) or named list of data
#'   frames (multi-dataset mode), as in [geneGroupSelectorServer()].
#' @param exprs Expression matrix/data frame (single dataset) or named list of
#'   matrices/data frames (multi-dataset mode), as in [geneGroupSelectorServer()].
#' @param cntr Optional contrast object, same format as [geneGroupSelectorServer()].
#' @param covar Sample covariate data frame (single dataset) or named list of
#'   data frames (multi-dataset mode).
#' @param palettes Optional reactive expression or `reactiveVal` returning a
#'   palette list forwarded to [plot_heatmap()]. Useful when heatmap colors are
#'   coordinated across modules.
#' @param sample_id_col Name of the sample ID column in `covar`.
#' @param dge_pval_col Optional p-value column name for DGE mode in
#'   [geneGroupSelectorServer()].
#' @param dge_lfc_col Optional log fold-change column name for DGE mode in
#'   [geneGroupSelectorServer()].
#' @param dge_fdr_col Optional adjusted p-value (FDR) column name for DGE mode in
#'   [geneGroupSelectorServer()].
#' @param max_genes Hard limit for the number of genes shown on the heatmap.
#'   If more genes are selected, only the first `max_genes` are displayed.
#' @param primary_id Primary identifier column in `annot` / `cntr`, also used to
#'   map selected genes to row labels in [plot_heatmap()].
#' @param annot_row_col Optional default annotation column used for heatmap row labels.
#'   Can be changed interactively in the module UI.
#'
#' @return A list with reactives: `genes`, `dataset`, and `heatmap`.
#'
#' @examples
#' annot <- data.frame(
#'   PrimaryID=paste0("g", 1:6),
#'   SYMBOL=LETTERS[1:6],
#'   stringsAsFactors=FALSE
#' )
#' exprs <- matrix(
#'   rnorm(24),
#'   nrow=6,
#'   dimnames=list(annot$PrimaryID, paste0("s", 1:4))
#' )
#' cntr <- list(
#'   contrast_a=data.frame(
#'     PrimaryID=annot$PrimaryID,
#'     padj=runif(6),
#'     pvalue=runif(6),
#'     log2FoldChange=rnorm(6),
#'     stringsAsFactors=FALSE
#'   )
#' )
#' covar <- data.frame(
#'   SampleID=colnames(exprs),
#'   group=c("A", "A", "B", "B"),
#'   stringsAsFactors=FALSE
#' )
#'
#' if(interactive()) {
#'   ui <- fluidPage(heatmapUI("hm"))
#'   server <- function(input, output, session) {
#'     heatmapServer("hm", annot=annot, exprs=exprs, cntr=cntr, covar=covar)
#'   }
#'   shinyApp(ui, server)
#' }
#'
#' if(interactive()) {
#'   data(C19)
#'   ui <- fluidPage(heatmapUI("hm"))
#'   server <- function(input, output, session) {
#'     heatmapServer(
#'       "hm",
#'       annot=C19$annotation,
#'       exprs=C19$expression,
#'       cntr=C19$contrasts,
#'       covar=C19$covariates,
#'       sample_id_col="label"
#'     )
#'   }
#'   shinyApp(ui, server)
#' }
#'
#' @export
heatmapServer <- function(id, annot, exprs=NULL, cntr=NULL, covar=NULL,
                          palettes=NULL,
                          sample_id_col="SampleID",
                          primary_id="PrimaryID",
                          dge_pval_col=NULL,
                          dge_lfc_col=NULL,
                          dge_fdr_col=NULL,
                          max_genes=150,
                          annot_row_col=NULL) {
  if(is.null(exprs)) {
    stop("`exprs` must be provided.")
  }

  annot_norm <- .gene_group_normalize_annot(annot)
  datasets <- names(annot_norm)
  exprs_norm <- .gene_group_normalize_exprs(exprs, datasets)
  covar_norm <- .heatmap_normalize_covar(covar, datasets)
  primary_id <- as.character(primary_id)[1]
  max_genes <- suppressWarnings(as.integer(max_genes)[1])
  annot_row_col <- as.character(annot_row_col)[1]

  if(is.na(primary_id) || !nzchar(primary_id)) {
    stop("`primary_id` must be a non-empty column name.")
  }
  if(is.na(max_genes) || max_genes < 1L) {
    stop("`max_genes` must be a positive integer.")
  }
  missing_primary_id <- datasets[!vapply(annot_norm, function(x) primary_id %in% colnames(x), logical(1))]
  if(length(missing_primary_id) > 0L) {
    stop(sprintf(
      "`primary_id` ('%s') not found in annotation for dataset(s): %s",
      primary_id,
      paste(missing_primary_id, collapse=", ")
    ))
  }

  if(!is.null(covar)) {
    missing_sample_col <- datasets[!vapply(covar_norm, function(x) sample_id_col %in% colnames(x), logical(1))]
    if(length(missing_sample_col) > 0L) {
      stop(sprintf(
        "`sample_id_col` ('%s') not found in covariates for dataset(s): %s",
        sample_id_col,
        paste(missing_sample_col, collapse=", ")
      ))
    }
  }

  moduleServer(id, function(input, output, session) {
    selected_ids <- reactiveVal(character(0))
    fig_size <- reactiveValues(width=800, height=600)
    palettes <- palettes %||% reactiveVal(NULL)

    heatmap_col_var <- list(values = list(type="continuous", breaks = c(-2, -1, 0, 1, 2)))
    #heatmap_col <- colorPalettesServer("heatmap_color", heatmap_col_var, compact=TRUE)

    selected_ids_for_heatmap <- reactive({
      ids <- selected_ids()
      if(length(ids) > max_genes) {
        ids[seq_len(max_genes)]
      } else {
        ids
      }
    })

    observeEvent(input$figure_size, {
      size <- sanitize_figsize(input$figure_size, default=c(800, 600))
      fig_size$width <- size$width
      fig_size$height <- size$height
    })

    selector <- geneGroupSelectorServer(
      "gene_selector",
      annot=annot,
      exprs=exprs,
      cntr=cntr,
      primary_id=primary_id,
      dge_pval_col=dge_pval_col,
      dge_lfc_col=dge_lfc_col,
      dge_fdr_col=dge_fdr_col,
      selected_ids=selected_ids
    )

    observe({
      ds <- selector$dataset()
      req(isTruthy(ds))
      req(ds %in% datasets)

      covar_choices <- .heatmap_annotation_choices(covar_norm[[ds]], sample_id_col=sample_id_col)
      isolate({ selected_covar <- intersect(input$sel_annot %||% character(0), covar_choices) })

      updateSelectizeInput(
        session=session,
        inputId="sel_annot",
        choices=covar_choices,
        selected=selected_covar,
        server=TRUE
      )

      row_choices <- .heatmap_row_label_choices(annot_norm[[ds]], primary_id=primary_id)
      isolate({ row_selected <- input$annot_row_col %||% "" })
      if(!nzchar(row_selected) && !is.na(annot_row_col) && nzchar(annot_row_col)) {
        row_selected <- annot_row_col
      }
      if(!row_selected %in% row_choices) {
        row_selected <- ""
      }

      updateSelectizeInput(
        session=session,
        inputId="annot_row_col",
        choices=c("Use primary IDs"="", row_choices),
        selected=row_selected,
        server=TRUE
      )
    })

    output$selected_genes_n <- renderText({
      sprintf("Selected genes: %d", length(selected_ids()))
    })

    output$selected_genes_warning <- renderText({
      n <- length(selected_ids())
      if(n <= max_genes) {
        return("")
      }
      sprintf("Warning: heatmap shows only first %d genes.", max_genes)
    })

    heatmap_obj <- reactive({
      ds <- selector$dataset()
      req(isTruthy(ds))
      req(ds %in% datasets)
      req(length(selected_ids_for_heatmap()) > 0L)

      plot_heatmap(
        exprs=exprs_norm[[ds]],
        genes=selected_ids_for_heatmap(),
        covar=covar_norm[[ds]],
        sample_id_col=sample_id_col,
        annot=annot_norm[[ds]],
        primary_id_col=primary_id,
        annot_row_col=if(isTruthy(input$annot_row_col)) input$annot_row_col else NULL,
        sel_annot=input$sel_annot,
        legend=isTRUE(input$show_legend),
        #col=heatmap_col()$values$pal,
        palettes=palettes()
      )
    })

    observe({ output$heatmap_plot <- renderPlot({
      hm <- heatmap_obj()
      .heatmap_draw(hm, show_legend=isTRUE(input$show_legend))
    }, width=fig_size$width, height=fig_size$height) })

    output$save <- downloadHandler(
      filename = function() {
        ds <- selector$dataset()
        md <- selector$modus()
        sprintf(
          "heatmap_%s_%s.pdf",
          sanitize_filename(ds, "dataset"),
          sanitize_filename(md, "modus")
        )
      },
      content = function(file) {
        hm <- heatmap_obj()
        save_pdf(
          file=file,
          width=fig_size$width / 75,
          height=fig_size$height / 75,
          draw=function() {
            .heatmap_draw(hm, show_legend=isTRUE(input$show_legend))
          }
        )
      }
    )

    return(list(
      genes=reactive(selected_ids()),
      dataset=selector$dataset,
      heatmap=heatmap_obj
    ))
  })
}
