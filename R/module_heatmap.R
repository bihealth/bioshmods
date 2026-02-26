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
  message(".heatmap_annotation_choices: Determining annotation choices from covariate dataset with columns: ", paste(colnames(covar_ds), collapse=", "))
  if(is.null(covar_ds) || !is.data.frame(covar_ds) || ncol(covar_ds) < 1L) {
    return(character(0))
  }

  setdiff(colnames(covar_ds), sample_id_col)
}

#' Heatmap Module UI
#'
#' @rdname heatmapServer
#' @export
heatmapUI <- function(id) {
  ns <- NS(id)

  sidebarLayout(
    sidebarPanel(
      geneGroupSelectorUI(ns("gene_selector")),
      shiny::tags$hr(),
      shiny::tags$div(
        class="well",
        h4("Heatmap options"),
        figsizeInput(ns("figure_size"), width="100%", selected="800x600"),
        selectizeInput(
          ns("sel_annot"),
          "Covariates to indicate",
          choices=character(0),
          selected=character(0),
          multiple=TRUE
        ),
        checkboxInput(ns("show_legend"), "Show legend", value=TRUE)
        , downloadButton(ns("save"), "Save heatmap to PDF", class="bg-success")
      ),
      width=4
    ),
    mainPanel(
      plotOutput(ns("heatmap_plot")),
      width=8
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
#' @param sample_id_col Name of the sample ID column in `covar`.
#' @param primary_id Primary identifier column in `annot` / `cntr`.
#' @param cntr_id_col Identifier column in contrast tables.
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
                          sample_id_col="SampleID",
                          primary_id="PrimaryID", cntr_id_col=primary_id) {
  if(is.null(exprs)) {
    stop("`exprs` must be provided.")
  }

  annot_norm <- .gene_group_normalize_annot(annot)
  datasets <- names(annot_norm)
  exprs_norm <- .gene_group_normalize_exprs(exprs, datasets)
  covar_norm <- .heatmap_normalize_covar(covar, datasets)

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
      cntr_id_col=cntr_id_col,
      selected_ids=selected_ids
    )

    observe({
      ds <- selector$dataset()
      req(isTruthy(ds))
      req(ds %in% datasets)

      covar_choices <- .heatmap_annotation_choices(covar_norm[[ds]], sample_id_col=sample_id_col)
      message("Updating covariate choices for dataset '", ds, "': ", paste(covar_choices, collapse=", "))
      isolate({ selected <- intersect(input$sel_annot %||% character(0), covar_choices) })

      updateSelectizeInput(
        session=session,
        inputId="sel_annot",
        choices=covar_choices,
        selected=selected,
        server=TRUE
      )
    })

    heatmap_obj <- reactive({
      ds <- selector$dataset()
      req(isTruthy(ds))
      req(ds %in% datasets)
      req(length(selected_ids()) > 0L)

      plot_heatmap(
        exprs=exprs_norm[[ds]],
        genes=selected_ids(),
        covar=covar_norm[[ds]],
        sample_id_col=sample_id_col,
        sel_annot=input$sel_annot,
        legend=isTRUE(input$show_legend)
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
