## Write module-scoped log messages through the shared package logger.
.volcano_log <- function(...) {
  .bioshmods_log(..., .prefix="volcano")
}

## Validate contrast and annotation inputs for the module.
.check_volcano_inputs <- function(cntr, annot, primary_id, lfc_col, pval_col, annot_show) {
  primary_id <- trimws(as.character(primary_id)[1])
  if(is.na(primary_id) || !nzchar(primary_id)) {
    stop("`primary_id` must be a non-empty character column name.")
  }

  cntr <- .check_dataset_contrasts(cntr, "cntr")
  annot <- .check_dataset_data_frames(annot, "annot", datasets=names(cntr))
  annot_full <- annot

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

  list(cntr = cntr, annot = annot, annot_full = annot_full)
}

# Generate readable ggplot code for the current volcano plot.
# The generated code mirrors the renderPlot body without running it.
.volcano_plot_chunk <- function(df, input, lfc_col, pval_col, primary_id,
                                annot_full, env_map=list(df="df")) {
  scales <- ifelse(input$samescaleX,
                   ifelse(input$samescaleY, "fixed", "free_y"),
                   ifelse(input$samescaleY, "free_x", "free"))

  code <- "# Prepare volcano plot data."
  code <- code %+n% sprintf("plot_df <- %s", env_map[["df"]])
  if(!identical(input$dataset, "_all")) {
    code <- code %+n% sprintf("plot_df <- plot_df[plot_df$Dataset == %s, ]", .r_code(input$dataset))
  }
  code <- code %+n% sprintf(
    "plot_df$Significant <- abs(plot_df$%s) > %s & plot_df$%s < %s",
    lfc_col, .r_code(input$lfc_thr), pval_col, .r_code(input$pval_thr)
  )
  code <- code %+n% sprintf("plot_df$y <- -log10(plot_df$%s)", pval_col)
  code <- code %+n% "\n# Pick top labels separately for each contrast."
  code <- code %+n% "label_df <- NULL"

  if(isTRUE(input$show_top_labels)) {
    label_col <- trimws(as.character(input$label_col %||% primary_id)[1])
    if(!nzchar(label_col) || is.na(label_col)) {
      label_col <- primary_id
    }

    code <- code %+n% sprintf("top_n <- %s", .r_code(as.integer(input$top_label_n)[1]))
    code <- code %+n% sprintf("label_df <- plot_df[!is.na(plot_df$%s) & plot_df$%s != \"\", ]", primary_id, primary_id)
    code <- code %+n% sprintf("label_df <- label_df[order(label_df$Dataset_Contrast, -label_df$Significant, -label_df$y, -abs(label_df$%s)), ]", lfc_col)
    code <- code %+n% sprintf("label_df <- do.call(rbind, lapply(split(label_df, label_df$Dataset_Contrast), head, n=top_n))")
    code <- code %+n% sprintf("label_df$label <- as.character(label_df$%s)", primary_id)
    code <- code %+n% sprintf("label_col <- %s", .r_code(label_col))
    code <- code %+n% sprintf("if(label_col != %s && exists(\"annot\")) {", .r_code(primary_id))
    code <- code %+n% "  label_df$label <- vapply(seq_len(nrow(label_df)), function(i) {"
    code <- code %+n% "    dataset <- label_df$Dataset[i]"
    code <- code %+n% sprintf("    id <- label_df$%s[i]", primary_id)
    code <- code %+n% "    annotation <- annot[[dataset]]"
    code <- code %+n% sprintf("    value <- annotation[[label_col]][match(id, annotation[[%s]])]", .r_code(primary_id))
    code <- code %+n% "    value <- as.character(value)[1]"
    code <- code %+n% "    if(is.na(value) || !nzchar(value)) as.character(id) else value"
    code <- code %+n% "  }, character(1))"
    code <- code %+n% "}"
  }

  code <- code %+n% "\n# Draw the volcano plot."
  code <- code %+n% sprintf("g <- ggplot(plot_df, aes(x=%s, y=y, color=Significant)) +", lfc_col)
  code <- code %+n% "  geom_point(alpha=.5) +"
  code <- code %+n% sprintf("  facet_wrap(~ Dataset_Contrast, scales=%s) +", .r_code(scales))
  code <- code %+n% "  scale_color_manual(values=c(\"TRUE\"=\"red\", \"FALSE\"=\"black\")) +"
  code <- code %+n% sprintf("  theme(text=element_text(size=%s)) +", .r_code(input$font_size))
  code <- code %+n% "  theme(legend.position=\"bottom\")"
  code <- code %+n% "if(!is.null(label_df) && nrow(label_df) > 0L) {"
  code <- code %+n% "  g <- g + geom_text(data=label_df, aes(label=label), check_overlap=TRUE,"
  code <- code %+n% sprintf("                     vjust=-0.4, size=max(2.5, %s / 4), show.legend=FALSE)", .r_code(input$font_size))
  code <- code %+n% "}"
  code <- code %+n% "g"

  list(required_pkgs=c("ggplot2"), type="figure", env_map=env_map, code=code)
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
                      checkboxInput(NS(id, "show_top_labels"), "Show top labels", value=FALSE))),
      fluidRow(column(width=12,
                      numericInput(NS(id, "top_label_n"), "Top labels (N)", value=10,
                                   min=1, step=1, width="80%"))),
      fluidRow(column(width=12, uiOutput(NS(id, "label_col_ui")))),
      fluidRow(
        column(width=6, downloadButton(NS(id, "save"), "Save plot to PDF", class="bg-success")),
        column(width=6, .report_add_button(id))
      ),
      fluidRow(tableOutput(NS(id, "point_id"))),
    width=2),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Plot",
          column(width=8,
          withSpinner(
                      plotOutput(NS(id, "volcanoPlot"), width="100%", height="100%",
                        hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                        click=NS(id, "plot_click"),
                        brush=NS(id, "plot_brush"))
                      )
          ),
          column(width=4,
            uiOutput(NS(id, "show_selected_ui")),
            tableOutput(NS(id, "sel_genes"))
          )
        ),
        .report_code_tab(id)
      ),
      width=10
    )
  )

}



#' Shiny module for displaying volcano plots
#'
#' Shiny module for displaying volcano plots
#' @param id module identifier (same as the one passed to volcanoUI)
#' @param cntr named list of dataset contrast lists. For a single dataset, use
#'   `list(default=cntr)`.
#' @param datasets character vector specifying datasets
#' @param lfc_thr default lfc threshold
#' @param lfc_col,pval_col names of the columns in the contrast data
#' frames which hold the log2 fold changes and p-values, respectively
#' @param pval_thr default p-value threshold
#' @param primary_id name of the primary ID column in contrasts and
#' annotation data frame.
#' @param annot named list of dataset annotation data frames. Names must match
#'   `cntr`; for a single dataset, use `list(default=annot)`.
#' @param gene_id must be a `reactiveValues` object. If not NULL, then
#' clicking on a gene identifier will modify this object (possibly
#' triggering an event in another module).
#' @param selection Optional `reactiveValues()` object with an `ids` element
#'   populated with the primary IDs from the current plot selection when the
#'   user clicks the `Show` button. If `NULL`, the button is not shown and no
#'   external state is updated.
#' @param ui_config Optional list configuring UI text. Supported keys:
#'   `show_button_label` for the label shown on the `Show` button.
#' @param report Optional report list. When supplied, clicking `Add to report`
#'        appends the current generated code chunk to `report$chunks`.
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
#'       cntr = list(default=cntr),
#'       annot = list(default=annot)
#'     )
#'   }
#'
#'   shinyApp(ui, server)
#' }
#'
#' # Example showing how to use shared selection to 
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
#'     selection <- reactiveValues(ids=character(0))
#'     volcanoServer("vol", 
#'                   cntr=list(default=C19$contrasts),
#'                   annot=list(default=C19$annotation),
#'                   selection=selection)
#'     heatmapServer(
#'       "hm",
#'       annot=list(default=C19$annotation),
#'       exprs=list(default=C19$expression),
#'       cntr=list(default=C19$contrasts),
#'       covar=list(default=C19$covariates),
#'       sample_id_col="label",
#'       selection=selection
#'     )
#'   }
#'   shinyApp(ui, server)
#' }
#' @export
volcanoServer <- function(id, cntr, lfc_col="log2FoldChange", pval_col="padj", 
                          primary_id="PrimaryID",
                          annot=NULL, gene_id=NULL,
                          selection=NULL,
                          ui_config=NULL,
                          report=NULL,
                          annot_show=c("SYMBOL", "ENTREZID")) {

  checked <- .check_volcano_inputs(cntr, annot, primary_id, lfc_col, pval_col, annot_show)
  cntr <- checked$cntr
  annot <- checked$annot
  annot_full <- checked$annot_full
  ui_config <- .normalize_show_button_ui_config(ui_config)

  .check_selection_reactivevalues(selection, arg_name="selection")

  if(!"default" %in% names(cntr)) {
    .volcano_log("running in multi data set mode.")
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

    output$label_col_ui <- renderUI({
      common_cols <- Reduce(intersect, lapply(annot_full, colnames))
      label_cols <- setdiff(common_cols, primary_id)
      if(length(label_cols) == 0L) {
        return(NULL)
      }

      selected_label_col <- isolate(input$label_col %||% primary_id)
      choices <- stats::setNames(c(primary_id, label_cols), c("Primary ID", label_cols))

      if(!selected_label_col %in% unname(choices)) {
        selected_label_col <- primary_id
      }

      selectInput(
        NS(id, "label_col"),
        "Label column",
        choices=choices,
        selected=selected_label_col,
        width="100%"
      )
    })

    selected_genes_table <- reactive({
      .df <- selected_genes()
      if(is.null(.df) || nrow(.df) == 0L) {
        return(NULL)
      }
      .ds <- .df[["Dataset"]][1]
      .ann <- annot[[.ds]][ match(.df[[primary_id]], annot[[.ds]][[primary_id]]), , drop=FALSE ]
      dplyr::bind_cols(.df["Dataset"], .ann)
    })

    output$sel_genes <- renderTable({
      .df <- selected_genes_table()
      if(is.null(.df)) { return(NULL) }
      link <- actionButton(NS(id, "gene_id~%s~%s"), label="%s \U25B6 ",
                           onclick=sprintf('Shiny.onInputChange(\"%s-genebutton\",  this.id)', id),
                           class = "btn-primary btn-sm")
      .df[[primary_id]] <- sprintf(as.character(link), .df[["Dataset"]], .df[[primary_id]], .df[[primary_id]])
      .df
    }, sanitize.text.function=function(x) x)

    output$show_selected_ui <- renderUI({
      .df <- selected_genes_table()
      req(!is.null(.df), nrow(.df) > 0L)
      buttons <- list(
        downloadButton(NS(id, "export_selected"), "Export to file", class="btn-default")
      )

      if(!is.null(selection)) {
        buttons <- c(
          list(actionButton(NS(id, "show_selected"), ui_config$show_button_label, class="btn-default")),
          buttons
        )
      }

      do.call(
        shiny::tags$div,
        c(
          list(style="display:flex;gap:8px;align-items:center;margin-bottom:8px;"),
          buttons
        )
      )
    })

    observeEvent(input$genebutton, {
      if(!is.null(gene_id)) {
        ids <- strsplit(input$genebutton, '~')[[1]]
        gene_id$ds <- ids[2]
        gene_id$id <- ids[3]
      }
    })

    observeEvent(input$show_selected, {
      req(!is.null(selection))
      .df <- selected_genes()
      req(!is.null(.df), nrow(.df) > 0L, primary_id %in% colnames(.df))
      ids <- unique(stats::na.omit(as.character(.df[[primary_id]])))
      selection$ids <- ids[nzchar(ids)]
    })

    output$export_selected <- downloadHandler(
      filename = function() {
        .ds <- input$dataset
        if(!isTruthy(.ds)) { .ds <- "_all" }
        sprintf("volcano_selected_genes_%s.xlsx", sanitize_filename(.ds, "all"))
      },
      content = function(file) {
        .df <- selected_genes_table()
        req(!is.null(.df), nrow(.df) > 0L)
        writexl::write_xlsx(list(selected_genes=.df), path=file)
      }
    )

    observeEvent(input$figure_size, {
      size <- sanitize_figsize(input$figure_size, default=c(800, 800))
      fig_size$width <- size$width
      fig_size$height <- size$height
    })

    observeEvent(input$show_top_labels, {
      if(isTRUE(input$show_top_labels)) {
        shinyjs::enable("top_label_n")
      } else {
        shinyjs::disable("top_label_n")
      }
    }, ignoreInit=FALSE)

    observeEvent(input$add_to_report, {
      .append_report_chunk(report, input$report_code)
    })

    observeEvent(input$add_to_report_main, {
      .append_report_chunk(report, input$report_code)
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

    ## Keep the report-code tab in sync with the current volcano controls.
    ## This generates text only; the rendered plot still uses the existing path.
    observe({
      req(input$dataset, input$lfc_thr, input$pval_thr)
      code <- .volcano_plot_chunk(
        df=df,
        input=input,
        lfc_col=lfc_col,
        pval_col=pval_col,
        primary_id=primary_id,
        annot_full=annot_full,
        env_map=list(df="df")
      )$code
      updateTextAreaInput(session, "report_code", value=code)
    })

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

      label_df <- NULL
      if(isTRUE(input$show_top_labels)) {
        top_n <- suppressWarnings(as.integer(input$top_label_n)[1])
        if(!is.na(top_n) && top_n > 0L) {
          label_df <- df %>%
            filter(!is.na(.data[[primary_id]]) & .data[[primary_id]] != "") %>%
            arrange(desc(.data[["Significant"]]), desc(.data[["y"]]),
                    desc(abs(.data[[lfc_col]]))) %>%
            distinct(.data[["Dataset_Contrast"]], .data[[primary_id]], .keep_all=TRUE) %>%
            group_by(.data[["Dataset_Contrast"]]) %>%
            slice_head(n=top_n) %>%
            ungroup()

          label_col <- trimws(as.character(input$label_col %||% primary_id)[1])
          if(!nzchar(label_col) || is.na(label_col)) {
            label_col <- primary_id
          }

          label_df[["label"]] <- vapply(seq_len(nrow(label_df)), function(i) {
            .ds <- label_df[["Dataset"]][i]
            .id <- label_df[[primary_id]][i]
            .ann <- annot_full[[.ds]]

            if(is.null(.ann) || !label_col %in% colnames(.ann)) {
              return(as.character(.id))
            }

            .val <- .ann[[label_col]][ match(.id, .ann[[primary_id]]) ]
            .val <- as.character(.val)[1]
            if(is.na(.val) || !nzchar(.val)) {
              as.character(.id)
            } else {
              .val
            }
          }, character(1))
        }
      }

      g <- ggplot(df, aes(x=.data[[lfc_col]], y=.data[["y"]],
                     color   =.data[["Significant"]])) +
        geom_point(alpha=.5) +
        facet_wrap(as.formula(paste('~', "Dataset_Contrast")), scales=scales) +
        scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
                                   theme(text=element_text(size=input$font_size)) +
                                   theme(legend.position="bottom")

      if(!is.null(label_df) && nrow(label_df) > 0L) {
        g <- g + geom_text(
          data=label_df,
          aes(label=.data[["label"]]),
          check_overlap=TRUE,
          vjust=-0.4,
          size=max(2.5, input$font_size / 4),
          show.legend=FALSE
        )
      }

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
