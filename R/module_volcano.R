## Check whether an object is a non-empty list of data frames.
.volcano_is_list_of_data_frames <- function(x) {
  is.list(x) && length(x) > 0L && all(vapply(x, is.data.frame, logical(1)))
}

## Check whether an object is a non-empty list of contrast lists.
.volcano_is_list_of_list_of_data_frames <- function(x) {
  is.list(x) && length(x) > 0L &&
    all(vapply(x, function(ds) .volcano_is_list_of_data_frames(ds), logical(1)))
}

## Collect unique primary IDs across all contrasts in one dataset.
.volcano_cntr_primary_ids <- function(cntr_ds, primary_id) {
  unique(unlist(lapply(cntr_ds, function(df) as.character(df[[primary_id]])), use.names = FALSE))
}

## Build a minimal annotation table from the contrast primary IDs.
.volcano_make_annotation_from_cntr <- function(cntr_ds, primary_id) {
  data.frame(
    stats::setNames(list(.volcano_cntr_primary_ids(cntr_ds, primary_id)), primary_id),
    check.names = FALSE
  )
}

## Normalize and validate contrast and annotation inputs for the module.
.normalize_volcano_inputs <- function(cntr, annot, primary_id, lfc_col, pval_col, annot_show) {
  primary_id <- trimws(as.character(primary_id)[1])
  if(is.na(primary_id) || !nzchar(primary_id)) {
    stop("`primary_id` must be a non-empty character column name.")
  }

  if (.volcano_is_list_of_data_frames(cntr)) {
    cntr <- list(default=cntr)
    if (is.null(annot) || is.data.frame(annot)) {
      annot <- list(default=annot)
    } else if (is.list(annot) && length(annot) == 1L &&
               (is.null(annot[[1]]) || is.data.frame(annot[[1]]))) {
      names(annot) <- "default"
    } else {
      stop("For single-dataset cntr, annot must be NULL, a data frame, or a single-element list.")
    }
  } else if (.volcano_is_list_of_list_of_data_frames(cntr)) {
    if (is.null(names(cntr)) || any(names(cntr) == "")) {
      stop("When cntr is a list of datasets, all datasets must be named.")
    }
    if (is.null(annot)) {
      annot <- stats::setNames(vector("list", length(cntr)), names(cntr))
    } else if (is.data.frame(annot)) {
      stop("When cntr is a list of datasets, annot must be a list of data frames (or NULL).")
    }
  } else {
    stop("cntr must be a list of data frames, or a named list of such lists.")
  }

  if (!is.list(annot)) {
    stop("annot must be NULL, a data frame, or a list of data frames.")
  }

  if (is.null(names(annot))) {
    if (length(annot) != length(cntr)) {
      stop("Unnamed annot list must have the same length as cntr.")
    }
    names(annot) <- names(cntr)
  } else if (any(names(annot) == "")) {
    stop("When annot is a list, all elements must be named.")
  }

  if (!all(names(cntr) %in% names(annot))) {
    stop("annot must contain an entry for each dataset in cntr.")
  }

  annot <- annot[names(cntr)]

  for (ds in names(cntr)) {
    cntr_ds <- cntr[[ds]]

    if (is.null(names(cntr_ds)) || any(names(cntr_ds) == "")) {
      stop(sprintf("All contrasts in dataset '%s' must be named.", ds))
    }

    for (cntr_name in names(cntr_ds)) {
      df <- cntr_ds[[cntr_name]]

      if (!primary_id %in% colnames(df)) {
        stop(sprintf(
          "Contrast '%s' in dataset '%s' is missing '%s'.",
          cntr_name, ds, primary_id
        ))
      }

      if (!all(c(lfc_col, pval_col) %in% colnames(df))) {
        stop(sprintf("Contrast '%s' in dataset '%s' must contain '%s' and '%s'.", cntr_name, ds, lfc_col, pval_col))
      }
    }

    cntr[[ds]] <- cntr_ds

    if (is.null(annot[[ds]])) {
      annot[[ds]] <- .volcano_make_annotation_from_cntr(cntr_ds, primary_id)
    }

    if (!is.data.frame(annot[[ds]])) {
      stop(sprintf("Annotation for dataset '%s' must be a data frame or NULL.", ds))
    }

    if (!primary_id %in% colnames(annot[[ds]])) {
      stop(sprintf("Annotation for dataset '%s' must contain '%s'.", ds, primary_id))
    }

    annot_ids <- as.character(annot[[ds]][[primary_id]])
    for (cntr_name in names(cntr_ds)) {
      cntr_ids <- unique(as.character(cntr_ds[[cntr_name]][[primary_id]]))
      missing_ids <- setdiff(cntr_ids, annot_ids)
      if (length(missing_ids) > 0L) {
        n_show <- min(length(missing_ids), 5L)
        stop(sprintf(
          "Annotation for dataset '%s' is missing %d '%s' values used by contrast '%s' (first %d: %s).",
          ds, length(missing_ids), primary_id, cntr_name, n_show,
          paste(missing_ids[seq_len(n_show)], collapse = ", ")
        ))
      }
    }

    keep_cols <- unique(c(primary_id, annot_show))
    annot[[ds]] <- annot[[ds]] %>% select(any_of(keep_cols))
  }

  list(cntr = cntr, annot = annot)
}

#' @rdname volcanoServer
#' @export
volcanoUI <- function(id, datasets=NULL, lfc_thr=1, pval_thr=.05) {

  datasets <- datasets %||% "default"

  if(length(datasets) == 1L) {
    ds_selector <- hidden(selectInput(NS(id, "dataset"), "Dataset:", datasets, 
                                      selected=datasets))
  } else {
    .ds <- c("_all", datasets)
    names(.ds) <- c("All datasets", datasets)
    ds_selector <- tipify(
                          selectInput(NS(id, "dataset"), "Dataset:", 
                               .ds, selected=.ds[1]),
                          "Choose the dataset to show (use \"all\" to show all data sets", placement="right")
  }

  sidebarLayout(
    sidebarPanel(
      fluidRow(column(width=12, ds_selector)),
      fluidRow(
        column(width=6,
        numericInput(NS(id, "pval_thr"), "P-value threshold:", value=pval_thr, 
                                     min=0, max=1),
               bsTooltip(NS(id, "pval_thr"), "P-value threshold for significant genes")),
        column(width=6,
        numericInput(NS(id, "lfc_thr"), "Log2 FC threshold:", value=lfc_thr),
        bsTooltip(NS(id, "lfc_thr"), "Log2 Fold Change threshold for significant genes")) 
      ),
      fluidRow(column(width=12, 
        tipify(checkboxInput(NS(id, "samescaleX"), "Same X scale for all plots",
          value=TRUE, width="100%"),
               "If checked, the X axis will be identical on all plots"),
                      )),
      fluidRow(column(width=12, 
        tipify(checkboxInput(NS(id, "samescaleY"), "Same Y scale for all plots",
          value=TRUE, width="100%"),
               "If checked, the y axis will be identical on all plots"),
                      )),
      fluidRow(column(width=6,
                      figsizeInput(NS(id, "figure_size"), width="100%"),
                    bsTooltip(NS(id, "figure_size"), "Change the figure size (in pixels, width x height)")),
        column(width=6, numericInput(NS(id, "font_size"), label="Font size", value = 12, 
                                 min=3, step=1, width="100%"),
                    bsTooltip(NS(id, "font_size"), "Change the font size of plot labels"))),
      fluidRow(column(width=12,
                      downloadButton(NS(id, "save"), "Save plot to PDF", class="bg-success"))),
      fluidRow(tableOutput(NS(id, "point_id"))),
    width=2),
    mainPanel(
      column(width=8,
      withSpinner(
                  plotOutput(NS(id, "volcanoPlot"), width="100%", height="100%",
                    hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                    click=NS(id, "plot_click"),
                    brush=NS(id, "plot_brush"))
                  )
      ),
      column(width=4,
        HTML("Click on the button to view an expression profile"),
        tableOutput(NS(id, "sel_genes")),
        uiOutput(NS(id, "show_selected_ui"))
      ), width=10
    )
  )

}



#' Shiny module for displaying volcano plots
#'
#' Shiny module for displaying volcano plots
#' @param id module identifier (same as the one passed to volcanoUI)
#' @param cntr either a named list of data frames, each being the results
#' of differential expression analysis for one contrast, or a list of data
#' sets, each data set being a named list of data frames.
#' @param datasets character vector specifying datasets
#' @param lfc_thr default lfc threshold
#' @param lfc_col,pval_col names of the columns in the contrast data
#' frames which hold the log2 fold changes and p-values, respectively
#' @param pval_thr default p-value threshold
#' @param primary_id name of the primary ID column in contrasts and
#' annotation data frame.
#' @param annot data frame with gene annotations (containing at least the
#' column specified with the `primary_id` parameter) or (if there are
#' multiple data sets) a named list of such data frames. Names of this list
#' must match the names of the `cntr` list.
#' @param gene_id must be a `reactiveValues` object. If not NULL, then
#' clicking on a gene identifier will modify this object (possibly
#' triggering an event in another module).
#' @param selected_ids Optional `reactiveVal()` populated with the primary IDs
#'   from the current plot selection when the user clicks the `Show` button.
#'   If `NULL`, the button is not shown and no external state is updated.
#' @param annot_show which columns from the annotation data frame should be
#' shown when mouse hovers over a gene
#'
#' @examples
#' annot <- data.frame(
#'   PrimaryID=paste0("g", 1:6),
#'   SYMBOL=LETTERS[1:6],
#'   stringsAsFactors=FALSE
#' )
#'
#' cntr <- list(
#'   contrast_a=data.frame(
#'     PrimaryID=annot$PrimaryID,
#'     log2FoldChange=rnorm(6),
#'     padj=runif(6),
#'     stringsAsFactors=FALSE
#'   )
#' )
#'
#' if(interactive()) {
#'   ui <- fluidPage(volcanoUI("volcano"))
#'
#'   server <- function(input, output, session) {
#'     volcanoServer(
#'       "volcano",
#'       cntr = cntr,
#'       annot = annot
#'     )
#'   }
#'
#'   shinyApp(ui, server)
#' }
#'
#' # Example showing how to use shared selected_ids to 
#' link volcano and heatmap modules
#'
#' data(C19)
#'
#' if(interactive()) {
#'   ui <- fluidPage(
#'     tabsetPanel(
#'       tabPanel("Volcano plot", volcanoUI("vol")),
#'       tabPanel("Heatmap", heatmapUI("hm"))
#'       )
#'   )
#'   server <- function(input, output, session) {
#'     selected_ids <- reactiveVal(character(0))
#'     volcanoServer("vol", 
#'                   cntr=C19$contrasts, 
#'                   annot=C19$annotation, 
#'                   selected_ids=selected_ids)
#'     heatmapServer(
#'       "hm",
#'       annot=C19$annotation,
#'       exprs=C19$expression,
#'       cntr=C19$contrasts,
#'       covar=C19$covariates,
#'       sample_id_col="label",
#'       selected_ids=selected_ids
#'     )
#'   }
#'   shinyApp(ui, server)
#' }
#' @export
volcanoServer <- function(id, cntr, lfc_col="log2FoldChange", pval_col="padj", 
                          primary_id="PrimaryID",
                          annot=NULL, gene_id=NULL,
                          selected_ids=NULL,
                          annot_show=c("SYMBOL", "ENTREZID")) {

  normalized <- .normalize_volcano_inputs(cntr, annot, primary_id, lfc_col, pval_col, annot_show)
  cntr <- normalized$cntr
  annot <- normalized$annot

  if(!is.null(selected_ids) && !inherits(selected_ids, "reactiveVal")) {
    stop("`selected_ids` must be NULL or a `reactiveVal()`.")
  }

  if(!"default" %in% names(cntr)) {
    message("volcanoServer: running in multi data set mode")
  }

  df <- .volcano_process_data(cntr, annot, primary_id,
    lfc_col, pval_col)
  df[["Dataset_Contrast"]] <- sprintf("%s\n%s", df[["Dataset"]], df[["Contrast"]])

  moduleServer(id, function(input, output, session) {

    fig_size       <- reactiveValues() ## figure height and width
    hover_genes    <- reactiveVal() ## hover gene selection, shown on the left
    selected_genes <- reactiveVal() ## active gene selection, shown on the right
    dfvar          <- reactiveVal() ## current data frame with genes
    plot_obj       <- reactiveVal()

    output$point_id <- renderTable({
      .df <- hover_genes()
      if(is.null(.df) || nrow(.df) == 0L) {
        return(NULL)
      }

      purrr::map_dfr(unique(.df[["Dataset"]]), function(.ds) {
        .sel <- .df[.df[["Dataset"]] == .ds, c("Dataset", primary_id), drop=FALSE]
        .ann <- annot[[.ds]][ match(.sel[[primary_id]], annot[[.ds]][[primary_id]]), , drop=FALSE ]
        dplyr::bind_cols(.sel["Dataset"], .ann)
      })
    })

    output$sel_genes <- renderTable({
      .df <- selected_genes()
      if(is.null(.df)) { return(NULL) }
      link <- actionButton(NS(id, "gene_id~%s~%s"), label="%s \U25B6 ",
                           onclick=sprintf('Shiny.onInputChange(\"%s-genebutton\",  this.id)', id),
                           class = "btn-primary btn-sm")
      .ds <- .df[["Dataset"]][1]
      .df <- annot[[.ds]][ match(.df[[primary_id]], annot[[.ds]][[primary_id]]), , drop=FALSE ]
      .df[[primary_id]] <- sprintf(as.character(link), .ds, .df[[primary_id]], .df[[primary_id]])
      .df
    }, sanitize.text.function=function(x) x)

    output$show_selected_ui <- renderUI({
      req(!is.null(selected_ids))
      .df <- selected_genes()
      req(!is.null(.df), nrow(.df) > 0L)
      actionButton(NS(id, "show_selected"), "Show", class="btn-default")
    })

    observeEvent(input$genebutton, {
      if(!is.null(gene_id)) {
        ids <- strsplit(input$genebutton, '~')[[1]]
        gene_id$ds <- ids[2]
        gene_id$id <- ids[3]
      }
    })

    observeEvent(input$show_selected, {
      req(!is.null(selected_ids))
      .df <- selected_genes()
      req(!is.null(.df), nrow(.df) > 0L, primary_id %in% colnames(.df))
      ids <- unique(stats::na.omit(as.character(.df[[primary_id]])))
      selected_ids(ids[nzchar(ids)])
    })

    observeEvent(input$figure_size, {
      size <- sanitize_figsize(input$figure_size, default=c(800, 800))
      fig_size$width <- size$width
      fig_size$height <- size$height
    })

    output$save <- downloadHandler(
      filename = function() {
        .ds <- input$dataset
        if(!isTruthy(.ds)) { .ds <- "_all" }
        sprintf("volcano_plot_%s.pdf", sanitize_filename(.ds, "all"))
      },
      content = function(file) {
        req(plot_obj())
        save_pdf(file=file, width=8, height=5, draw=function() {
          print(plot_obj())
        })
      }
    )

    observeEvent(input$plot_hover, {
      .df <- dfvar()
      np <- nearPoints(.df, input$plot_hover, xvar = lfc_col, yvar = "y")
      hover_genes(np[ , c("Dataset", primary_id), drop=FALSE ])
    })

    observeEvent(input$plot_brush, {
      .df <- dfvar()
      np <- brushedPoints(.df, input$plot_brush, xvar = lfc_col, yvar = "y")
      selected_genes(np)
    })

    observeEvent(input$plot_click, {
      .df <- dfvar()
      np <- nearPoints(.df, input$plot_click, xvar = lfc_col, yvar = "y")
      selected_genes(np)
    })


    observe({ output$volcanoPlot <- renderPlot({

      if(input$dataset != "_all") {
        df <- df %>% filter(.data[["Dataset"]] == input$dataset)
      } 

      df <- df %>% mutate(Significant=
                     abs(.data[[lfc_col]]) > input$lfc_thr &
                     .data[[pval_col]] < input$pval_thr)

      ## trickery to fool ggplot in the unholy combination
      ## with nearPoints()
      df <- df %>% mutate(y = -log10(.data[[pval_col]]))

      #yvar <- sprintf("-log10(%s)", pval_col)
      #df[[ yvar ]] <- -log10(df[[pval_col]])


      scales <- ifelse(input$samescaleX, 
                       ifelse(input$samescaleY,
                              "fixed",
                              "free_y"),
                       ifelse(input$samescaleY,
                              "free_x",
                              "free"))

      ## store the data frame for click, hover or brush events
      dfvar(df)
      # print(head(df, 100))

      g <- ggplot(df, aes(x=.data[[lfc_col]], y=.data[["y"]],
                     color   =.data[["Significant"]])) +
        geom_point(alpha=.5) +
        facet_wrap(as.formula(paste('~', "Dataset_Contrast")), scales=scales) +
        scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
                                   theme(text=element_text(size=input$font_size)) +
                                   theme(legend.position="bottom")
      # ggsave("debug_volcano_plot.pdf", g, width=fig_size$width/100, height=fig_size$height/100)
      # saveRDS(g, "debug_volcano_plot.rds")
      plot_obj(g)
      g

      }, width=fig_size$width, 
         height=fig_size$height)

    })
  })

    

}


## create one huge data frame for all contrasts and data sets
.volcano_process_data <- function(cntr, annot, primary_id, lfc_col, pval_col) {


  df <- imap_dfr(cntr, ~ {
    .volcano_process_data_one_ds(.y, .x, annot[[.y]], primary_id, lfc_col, pval_col)
                          })
  return(df)
}


## create one huge data frame for all contrasts
.volcano_process_data_one_ds <- function(ds_id, cntr, annot, primary_id, lfc_col, pval_col) {
  df <- imap_dfr(cntr, ~ {
             stopifnot(!is.null(colnames(.x)))
             stopifnot(all(c(primary_id, lfc_col, pval_col) %in% colnames(.x)))

             .x[["Dataset"]] <- ds_id
             .x[["Contrast"]] <- .y
             return(.x)
             
          })

  selcols <- c(primary_id, lfc_col, pval_col, "Dataset", "Contrast")
  df <- df[ , selcols ]
  return(df)
}
