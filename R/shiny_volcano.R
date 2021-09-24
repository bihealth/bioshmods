#' @rdname volcanoServer
#' @export
volcanoUI <- function(id, datasets=NULL, lfc_thr=1, pval_thr=.05) {

  if(is.null(datasets)) { 
    datasets <- "default"
  }

  if(length(datasets) == 1L) {
    ds_selector <- hidden(selectInput(NS(id, "dataset"), "Dataset:", datasets, 
                                      selected=datasets))
  } else {
    .ds <- c("_all", datasets)
    names(.ds) <- c("All datasets", datasets)
    ds_selector <- tipify(
                          selectInput(NS(id, "dataset"), "Dataset:", 
                               .ds, selected=.ds[1]),
                          "Choose the dataset to show (use \"all\" to show all data sets")
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
      fluidRow(tableOutput(NS(id, "point_id"))),
    width=2),
    mainPanel(
      column(width=10,
      withSpinner(
                  plotOutput(NS(id, "volcanoPlot"), width="100%", height="100%",
                    hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                    click=NS(id, "plot_click"),
                    brush=NS(id, "plot_brush"))
                  )
      ),
      column(width=2,
        HTML("Click on the button to view an expression profile"),
        tableOutput(NS(id, "sel_genes"))
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
#' annotation data frame
#' @param annot data frame with gene annotations (containing at least the
#' column specified with the `primary_id` parameter) or (if there are
#' multiple data sets) a named list of such data frames. Names of this list
#' must match the names of the `cntr` list.
#' @param gene_id must be a `reactiveValues` object. If not NULL, then
#' clicking on a gene identifier will modify this object (possibly
#' triggering an event in another module).
#' @export

volcanoServer <- function(id, cntr, lfc_col="log2FoldChange", pval_col="padj", 
                          primary_id="PrimaryID",
                          annot=NULL, gene_id=NULL) {

  if(!is.data.frame(cntr[[1]])) {
    message("volcanoServer: running in multilevel mode")
  } else {
    cntr  <- list(default=cntr)
    annot <- list(annot=cntr)
  }

  stopifnot(all(names(cntr) %in% names(annot)))

  df <- .volcano_process_data(cntr, annot, primary_id,
    lfc_col, pval_col)
  df[["Dataset_Contrast"]] <- sprintf("%s\n%s", df[["Dataset"]], df[["Contrast"]])

  moduleServer(id, function(input, output, session) {

    fig_size       <- reactiveValues() ## figure height and width
    hover_genes    <- reactiveVal() ## hover gene selection, shown on the left
    selected_genes <- reactiveVal() ## active gene selection, shown on the right
    dfvar          <- reactiveVal() ## current data frame with genes

    output$point_id <- renderTable({
      hover_genes()
    })

    output$sel_genes <- renderTable({
      .df <- selected_genes()
      if(is.null(.df)) { return(NULL) }
      link <- actionButton(NS(id, "gene_id~%s~%s"), label="%s \U25B6 ",
                           onclick=sprintf('Shiny.onInputChange(\"%s-genebutton\",  this.id)', id),
                           class = "btn-primary btn-sm")
      .ds <- .df[["Dataset"]][1]
      .df[[primary_id]] <- sprintf(as.character(link), .ds, .df[[primary_id]], .df[[primary_id]])
      .df[ , primary_id, drop=FALSE ]
    }, sanitize.text.function=function(x) x)

    observeEvent(input$genebutton, {
      if(!is.null(gene_id)) {
        ids <- strsplit(input$genebutton, '~')[[1]]
        gene_id$ds <- ids[2]
        gene_id$id <- ids[3]
      }
    })

    observeEvent(input$figure_size, {
        fig_size$width  <- 
          as.numeric(gsub(" *([0-9]+) *x *([0-9]+)", "\\1", input$figure_size))
        fig_size$height <- 
          as.numeric(gsub(" *([0-9]+) *x *([0-9]+)", "\\2", input$figure_size))
    })

    observeEvent(input$plot_hover, {
      .df <- dfvar()
      np <- nearPoints(.df, input$plot_hover)
      hover_genes(np[ , primary_id, drop=FALSE ])
    })

    observeEvent(input$plot_brush, {
      .df <- dfvar()
      np <- brushedPoints(.df, input$plot_brush)
      selected_genes(np)
    })

    observeEvent(input$plot_click, {
      .df <- dfvar()
      np <- nearPoints(.df, input$plot_click)
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
      yvar <- sprintf("-log10(%s)", pval_col)
      df[[ yvar ]] <- -log10(df[[pval_col]])


      scales <- ifelse(input$samescaleX, 
                       ifelse(input$samescaleY,
                              "fixed",
                              "free_y"),
                       ifelse(input$samescaleY,
                              "free_x",
                              "free"))

      ## store the data frame for click, hover or brush events
      dfvar(df)

      ## loads of trickery to get around the dumb decision of using
      ## unquoted vars in ggplot (because typing quotes is SO hard 
      ## so we will make living hell out of an otherwise nice framework)
      ggplot(df, aes_string(x=lfc_col, y=yvar,
                     color   ="Significant",
                     dscon   ="Dataset_Contrast",
                     dataset ="Dataset",
                     contrast="Contrast")) +
        geom_point(alpha=.5) +
        facet_wrap(as.formula(paste('~', "Dataset_Contrast")), scales=scales) +
        scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
                                   theme(text=element_text(size=input$font_size)) +
                                   theme(legend.position="bottom")

      }, width=fig_size$width, 
         height=fig_size$height)

    })
  })

    

}


.volcano_process_data <- function(cntr, annot, primary_id, lfc_col, pval_col) {


  df <- imap_dfr(cntr, ~ {
    .volcano_process_data_one_ds(.y, .x, annot[[.y]], primary_id, lfc_col, pval_col)
                          })


  return(df)
}

.volcano_process_data_one_ds <- function(ds_id, cntr, annot, primary_id, lfc_col, pval_col) {

  df <- imap_dfr(cntr, ~ {
             stopifnot(!is.null(colnames(.x)))

             if(!primary_id %in% colnames(.x) &&
                !is.null(rownames(.x))) {
              warning(".volcano_process_data: primary_id not found in contrast data frame")
              .x[[primary_id]] <- rownames(.x)
             }

             stopifnot(all(c(lfc_col, pval_col) %in% colnames(.x)))

             .x[["Dataset"]] <- ds_id
             .x[["Contrast"]] <- .y
             return(.x)
             
          })

  selcols <- c(primary_id, lfc_col, pval_col, "Dataset", "Contrast")
  df <- df[ , selcols ]
  return(df)
}



