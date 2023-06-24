#' @rdname volcanoServer
#' @export
volcanoUI <- function(id, datasets=NULL, lfc_thr=1, pval_thr=.05) {
  ns <- NS(id)

  if(is.null(datasets)) { datasets <- "default" }

  if(length(datasets) == 1L) {
    ds_selector <- hidden(selectInput(ns("dataset"), "Dataset:", datasets, 
                                      selected=datasets))
  } else {
    .ds <- c("_all", datasets)
    names(.ds) <- c("All datasets", datasets)
    ds_selector <- tipify(
                          selectInput(ns("dataset"), "Dataset:", 
                               .ds, selected=.ds[1]),
                          "Choose the dataset to show (use \"all\" to show all data sets", placement="right")
  }

  sidebarLayout(
    sidebarPanel(
      fluidRow(column(width=12, ds_selector)),
      fluidRow(
        column(width=6,
        numericInput(ns("pval_thr"), "P-value threshold:", value=pval_thr, 
                                     min=0, max=1),
               bsTooltip(ns("pval_thr"), "P-value threshold for significant genes")),
        column(width=6,
        numericInput(ns("lfc_thr"), "Log2 FC threshold:", value=lfc_thr),
        bsTooltip(ns("lfc_thr"), "Log2 Fold Change threshold for significant genes")) 
      ),
      fluidRow(column(width=12, 
        tipify(checkboxInput(ns("samescaleX"), "Same X scale for all plots",
          value=TRUE, width="100%"),
               "If checked, the X axis will be identical on all plots"),
                      )),
      fluidRow(column(width=12, 
        tipify(checkboxInput(ns("samescaleY"), "Same Y scale for all plots",
          value=TRUE, width="100%"),
               "If checked, the y axis will be identical on all plots"),
                      )),
      fluidRow(column(width=6,
                      figsizeInput(ns("figure_size"), width="100%"),
                    bsTooltip(ns("figure_size"), "Change the figure size (in pixels, width x height)")),
        column(width=6, numericInput(ns("font_size"), label="Font size", value = 12, 
                                 min=3, step=1, width="100%"),
                    bsTooltip(ns("font_size"), "Change the font size of plot labels"))),
      fluidRow(tableOutput(ns("point_id"))),
    width=2),
    mainPanel(
      column(width=8,
        fluidRow(
          tipify(downloadButton(ns("volcano_savePDF"), "Save PDF", class="bg-success"), "Save image as PDF"),
          tipify(actionButton(ns("rmd"), "Add to report", class="bg-success"), "Add image to the markdown report")
        ),
        fluidRow(
          withSpinner(
                  plotOutput(ns("volcanoPlot"), width="100%", height="100%",
                    hover=hoverOpts(ns("plot_hover"), delay=50, delayType="throttle"),
                    click=ns("plot_click"),
                    brush=ns("plot_brush"))
                  ))
      ),
      column(width=4,
        HTML("Click on the button to view an expression profile"),
        tableOutput(ns("sel_genes"))
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
#' @param rmd_var a reactive values object which will be used to store the
#'        generated markdown chunks. 
#' @param gene_id must be a `reactiveValues` object. If not NULL, then
#' clicking on a gene identifier will modify this object (possibly
#' triggering an event in another module).
#' @param annot_show which columns from the annotation data frame should be
#' shown when mouse hovers over a gene
#' @export

volcanoServer <- function(id, cntr, annot, lfc_col="log2FoldChange", pval_col="padj", 
                          primary_id="PrimaryID",
                          gene_id=NULL, 
                          annot_show=c("SYMBOL", "ENTREZID"),
                          rmd_var=NULL) {
  parent_frame <- parent.frame()

  varnames <- list(cntr  = deparse(substitute(cntr)),
                   annot = deparse(substitute(annot)))

  if(!is.data.frame(cntr[[1]])) {
    message("volcanoServer: running in multi data set mode")
    mode <- "multi" 
    stopifnot(all(names(cntr) %in% names(annot)))
    annot <- map(annot, ~ .x[ , colnames(.x) %in% c(primary_id, annot_show), drop=FALSE ])
  } else {
    #cntr  <- list(default=cntr)
    #annot <- list(default=annot)
    annot <- annot[ , colnames(annot) %in% c(primary_id, annot_show), drop=FALSE ]
    mode <- "single"
  }



  moduleServer(id, function(input, output, session) {

    code_setup       <- .volcano_generate_preprocess_chunk(id, mode, varnames, primary_id, lfc_col, pval_col, local=FALSE) 
    code_setup_local <- .volcano_generate_preprocess_chunk(id, mode, varnames, primary_id, lfc_col, pval_col) 
    #setup_chunk     <- chunk_generate_rmd(id, code=code_setup, label="initialization", type="setup")
    observe({
      isolate({
        msg("adding setup chunk")
        add_chunk(id=id, type="setup", code=code_setup, rmd_var=rmd_var, title=glue("volcano plots setups for {id}"), label="initialization")
      })
    })
    msg("Running code:\n", code_setup_local)
    eval(parse(text=code_setup_local))
    #print(volcano_df)

    plot_code      <- reactiveVal()
    fig_size       <- reactiveValues() ## figure height and width
    hover_genes    <- reactiveVal() ## hover gene selection, shown on the left
    selected_genes <- reactiveVal() ## active gene selection, shown on the right
    dfvar          <- reactiveVal() ## current data frame with genes

    output$point_id <- renderTable({ hover_genes() })

    output$sel_genes <- renderTable({
      .volcano_selected_genes_table(id, selected_genes(), mode, primary_id, annot)
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

    observeEvent(input$rmd, {
      msg("Add to markdown clicked")
      .code <- plot_code()

      .cap <- .volcano_plot_title(input$dataset, mode)
      add_chunk(id=id, type="plot", code=.code, rmd_var=rmd_var, title=.cap, 
                label="volcano", fig.cap=.cap, fig.width=fig_size$width, 
                fig.height=fig_size$height)
    })

    # need to observe the renderPlot assignment because fig_size is
    # reactive
    observe({ output$volcanoPlot <- renderPlot({
      #print(volcano_df)
      plot_code(.volcano_generate_plot_code(id, input, lfc_col, pval_col))
      msg("plot code: ", plot_code())
      eval(parse(text=plot_code()))

      dfvar(df)
      msg("plotting done")
      g
      }, width=fig_size$width, 
         height=fig_size$height)

    })

    ## Save figure as a PDF
    output$volcano_savePDF <- downloadHandler(
      filename = function() {
        .ds <- .volcano_ds2str(input$dataset, mode)
        ret <- sprintf("volcano_plot_%s_%s.pdf", .ds, id)
        ret <- gsub("[^0-9a-zA-Z_.-]", "", ret)
        return(ret)
      },
      content = function(file) {
        pdf(file=file, width=floor(fig_size$width / 72), height=floor(fig_size$height / 72))
        code <- plot_code()
        print(eval(parse(text=code)))
        dev.off()
      }
    )
 
  })

    

}


.volcano_plot_title <- function(datsets, mode) {
  ds <- .volcano_ds2str(datasets, mode)

  if(mode == "single") {
    return("Volcano plot")
  } else if(ds == "all") {
    return("Volcano plot for all datasets")
  } else if(length(datasets) > 1) {
    return(sprintf("Volcano plot for datasets %s", ds))
  } else {
    return(sprintf("Volcano plot for dataset %s", ds))
  }

}

.volcano_ds2str <- function(datasets, mode) {

  if(mode == "single") {
    return("")
  }

  if(input$dataset == "_all") {
    return("all")
  }

  if(length(input$dataset) > 1) {
    return(paste(input$dataset, collapse=", "))
  }

  return(input$dataset)
}

## code side effects: g (the plot), df (data frame underlying plot)
.volcano_generate_plot_code <- function(id, input, lfc_col, pval_col) {
  ret <- glue("df <- volcano_df")

  ret <- glue('{ret}\ndf <- df %>%')

  if(input$dataset != "_all") {
    ret <- glue('
{ret}\n  filter(Dataset == "{input$dataset}") %>%')
  } 

  ret <- glue('
{ret}
  mutate(Significant=
                      abs({lfc_col}) > {input$lfc_thr} &
                      {pval_col} < {input$pval_thr}) %>%
  mutate(x = {lfc_col},
         y = -log10({pval_col}))')

  ## trickery to fool ggplot in the unholy combination
  ## with nearPoints(). Otherwise it won't work properly


  scales <- ifelse(input$samescaleX, 
                   ifelse(input$samescaleY,
                          "fixed",
                          "free_y"),
                   ifelse(input$samescaleY,
                          "free_x",
                          "free"))

  ret <- glue('
{ret}

g <- ggplot(df, aes(x=x, y=y,
                 color   =Significant,
                 dscon   =Dataset_Contrast,
                 dataset =Dataset,
                 contrast=Contrast)) +
    geom_point(alpha=.5) +
    ylab("-log10({pval_col})") +
    xlab("{lfc_col}") +
    facet_wrap(~ Dataset_Contrast, scales="{scales}") +
    scale_color_manual(values=c("TRUE"="red", "FALSE"="black")) +
                               theme(text=element_text(size={input$font_size})) +
                               theme(legend.position="bottom")
g
')
  return(ret)
}

# if local=TRUE, then the code is for local consumption; if FALSE, it is to
# be evaluated in the parent environment
.volcano_generate_preprocess_chunk <- function(id, mode, varnames, primary_id, lfc_col, pval_col, local=TRUE) {
  cntr_v <- "cntr"
  annot_v <- "annot"
  if(!local) { 
    cntr_v  <- varnames$cntr 
    annot_v <- varnames$annot
  }


  ret <- "\n## preprocess code for volcano plots\nrequire(tidyverse)\n\n"
  ret <- glue("{ret}\nvolcano_df <- ")

  if(mode == "multi") {
    ret <- glue('
{ret}imap_dfr({cntr_v}, ~ {{
  dataset_id <- .y
  imap_dfr(.x, ~ {{
    .x[["Dataset"]]  <- dataset_id
    .x[["Contrast"]] <- .y
    .x
  }})
}})')
  } else {
    ret <- glue('
{ret}imap_dfr({cntr_v}, ~ {{
  .x[["Dataset"]] <- "default"\n  .x[["Contrast"]] <- .y
  .x
}})')
  }

  ret <- glue('
{ret}\n\n## Choose columns to show
selcols <- c("{primary_id}", "{lfc_col}", "{pval_col}", "Dataset", "Contrast")
volcano_df <- volcano_df[ , colnames(volcano_df) %in% selcols]
volcano_df[["Dataset_Contrast"]] <- sprintf("%s\\n%s", 
                                            volcano_df[["Dataset"]], 
                                            volcano_df[["Contrast"]])
  
volcano_annot <- {annot_v}

## end of volcano preprocess code') 

  return(ret)
}


## generates the table on the right containing the details of the selected
## genes
.volcano_selected_genes_table <- function(id, df, mode, primary_id, annot) {
  if(is.null(df)) { 
    msg("selected_genes is NULL")
    return(NULL) 
  }

  msg("sel_genes:")
  print(df)

  link <- actionButton(NS(id, "gene_id~%s~%s"), label="%s \U25B6 ",
                       onclick=sprintf('Shiny.onInputChange(\"%s-genebutton\",  this.id)', id),
                       class = "btn-primary btn-sm")
  .ds <- df[["Dataset"]][1]

  msg("dataset", .ds)
  if(mode == "multi") {
    df <- annot[[.ds]][ match(df[[primary_id]], annot[[.ds]][[primary_id]]), ]
  } else {
    df <- annot[ match(df[[primary_id]], annot[[primary_id]]), ]
  }

  msg("sel_genes:")
  print(df)
  df[[primary_id]] <- sprintf(as.character(link), .ds, df[[primary_id]], df[[primary_id]])
  df
} 
