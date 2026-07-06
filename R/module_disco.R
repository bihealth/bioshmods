.disco_log <- function(...) {
  .bioshmods_log(..., .prefix="disco")
}

#' @rdname discoServer
#' @importFrom shiny textInput
#' @export
discoUI <- function(id, cntr_titles) {

  cntr_titles <- .prep_cntr_titles(cntr_titles)
  cntr_flat   <- unlist(cntr_titles, recursive=FALSE)

  if(!length(cntr_flat) > 1) {
    h4("You need at least two contrasts for this plot")
  } else {

  fluidRow(
    column(width=3,
        fluidRow(
                 column(width=6, selectInput(NS(id, "contrast1"), label = "Contrast 1", 
                             choices = cntr_titles, width="100%")),
        column(width=6, selectInput(NS(id, "contrast2"), label = "Contrast 2", 
                             choices = cntr_titles, selected=cntr_flat[2], width="100%"))),
        fluidRow(uiOutput(NS(id, "matchby"))),
        column(width=12,
        fluidRow(checkboxInput(NS(id, "autoscale"), "Automatic scale", value=TRUE)),
        fluidRow(sliderInput(NS(id, "min"), "Min", min=-150, max=0, value=-100, width="80%")),
        fluidRow(sliderInput(NS(id, "max"), "Max", min=0, max=150, value=100, width="80%")),
        fluidRow(checkboxInput(NS(id, "show_top_labels"), "Show top labels", value=FALSE)),
        fluidRow(numericInput(NS(id, "top_label_n"), "Top labels (N)", value=10, min=1, step=1, width="80%")),
        fluidRow(downloadButton(NS(id, "save"), "Save plot to PDF", class="bg-success")),
        fluidRow(
                 textInput(NS(id, "glabs"), 
                           "Type a comma separated list of gene IDs to label on the plot",
                           width="80%"),
                 actionButton(NS(id, "glabsgo"), label="Update", icon=icon("fa-pen"))
                 ),
        fluidRow(verbatimTextOutput(NS(id, "msg")))
        )
    ),

    column(width=9,
      tabsetPanel(
        tabPanel(
          "Plot",
          withSpinner(plotOutput(NS(id, "discoplot"), 
                     hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                     click=NS(id, "plot_click"),
                     brush=NS(id, "plot_brush"),
                     height="600px")
          ),
          HTML("<br/>Hover to identify genes, click to select, or click & drag to select an area<br/><br/>"),
          fluidRow(tableOutput(NS(id, "point_id")))
        ),
        tabPanel(
          "Genes",
          HTML("Click on the button to view an expression profile"),
          tableOutput(NS(id, "sel_genes"))
        ),
        .report_code_tab(id)
      )
    )
  )
  }
}

# Merge selected point IDs with annotation metadata from both datasets.
# Returns one table with requested columns, suffixed by `_1` and `_2`.
.get_gene_df <- function(pid, selcols, primary_id="PrimaryID", annot1=NULL, annot2=NULL) {

  pid1 <- pid[["__primary_id_1"]]
  pid2 <- pid[["__primary_id_2"]]

  stopifnot(is.data.frame(annot1), is.data.frame(annot2))
  stopifnot(primary_id %in% colnames(annot1), primary_id %in% colnames(annot2))

  cols <- unique(c(primary_id, selcols))
  pid1 <- data.frame(pid1, stringsAsFactors=FALSE)
  pid2 <- data.frame(pid2, stringsAsFactors=FALSE)
  colnames(pid1) <- primary_id
  colnames(pid2) <- primary_id

  pid1 <- pid1 %>%
    left_join(annot1 %>% select(any_of(cols)), by=primary_id) %>%
    select(any_of(cols))

  pid2 <- pid2 %>%
    left_join(annot2 %>% select(any_of(cols)), by=primary_id) %>%
    select(any_of(cols))

  colnames(pid1) <- paste0(colnames(pid1), "_1")
  colnames(pid2) <- paste0(colnames(pid2), "_2")

  return(cbind(pid1, pid2))
}

# Prepare plotted disco points by removing missing values and clipping score range.
# Returns a selection-ready data frame sorted by absolute disco value.
.prep_disco_selection_df <- function(df, lower, upper) {
  ret <- df %>%
    filter(!is.na(.data$log2FoldChange.x) & !is.na(.data$log2FoldChange.y) & !is.na(.data$disco)) %>%
    mutate(disco=ifelse(.data$disco > upper, upper, ifelse(.data$disco < lower, lower, .data$disco))) %>%
    arrange(abs(.data$disco))
  # print(ret)
  ret
}

# Validate `cntr`/`annot` as explicit dataset-keyed lists.
# Ensures required ID columns exist and annotations cover all contrast identifiers.
.check_disco_inputs <- function(cntr, annot, primary_id) {
  cntr <- .check_dataset_contrasts(cntr, "cntr")
  annot <- .check_dataset_data_frames(annot, "annot", datasets=names(cntr))

  for (ds in names(cntr)) {
    cntr_ds <- cntr[[ds]]

    if (!all(vapply(cntr_ds, function(df) primary_id %in% colnames(df), logical(1)))) {
      bad <- names(cntr_ds)[!vapply(cntr_ds, function(df) primary_id %in% colnames(df), logical(1))]
      stop(sprintf(
        "In dataset '%s', all contrasts must contain the '%s' column. Missing in: %s",
        ds, primary_id, paste(bad, collapse = ", ")
      ))
    }

    if (!primary_id %in% colnames(annot[[ds]])) {
      stop(sprintf("Annotation for dataset '%s' must contain the '%s' column.", ds, primary_id))
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
  }

  list(cntr = cntr, annot = annot)
}

#' Shiny Module – disco plots
#'
#' Shiny Module – disco plots
#' @param id identifier of the shiny module (character vector)
#' @param primary_id name of the contrast data frame column with the primary IDs
#' @param cntr named list of dataset contrast lists. Contrast data frames must
#'        have the columns log2FoldChange and pvalue. For a single dataset, use
#'        `list(default=cntr)`.
#' @param annot named list of dataset annotation data frames. Each annotation
#'        data frame must have a column named "PrimaryID". For a single dataset,
#'        use `list(default=annot)`.
#' @param selcols which column in the gene table when genes are selected
#'        from the plot
#' @param cntr_titles character vector containing the IDs of the contrasts
#'        (same as `names(cntr)`).
#' @param gene_id must be a `reactiveValues` object. If not NULL, then
#' clicking on a gene identifier will modify this object (possibly
#' triggering an event in another module).
#' @param report Optional report list. When supplied, clicking `Add to report`
#'        appends the current generated code chunk to `report$chunks`.
#' @return Returns a reactive expression returning the ID of the activated gene
#' @example inst/examples/disco.R
#' @export
discoServer <- function(id, cntr, annot=NULL,
    selcols=c("PrimaryID", "ENTREZ", "SYMBOL"),
    primary_id="PrimaryID", gene_id=NULL,
    report=NULL) {

  checked <- .check_disco_inputs(cntr, annot, primary_id)
  cntr <- checked$cntr
  annot <- checked$annot

  if(!"default" %in% names(cntr)) {
    .disco_log("running in multi dataset mode.")
  }

  link <- actionButton(NS(id, "gene_id~%s~%s"), label="%s \U25B6 ",
                           onclick=sprintf('Shiny.onInputChange(\"%s-genebutton\",  this.id)', id),
                           class = "btn-primary btn-sm")

  moduleServer(id, function(input, output, session) {
    disable("min")
    disable("max")

    disco          <- reactiveVal()
    current_genes  <- reactiveVal()
    selected_genes <- reactiveVal()

    contrast1 <- reactiveVal()
    dataset1  <- reactiveVal()
    contrast2 <- reactiveVal()
    dataset2  <- reactiveVal()
    gene_labs <- reactiveVal()
    plot_obj  <- reactiveVal()
    plot_df   <- reactiveVal()

    env_maps <- lapply(names(cntr), function(ds) {
      lapply(names(cntr[[ds]]), function(contrast) {
        list(
          contrast=sprintf('cntr[["%s"]][["%s"]]', ds, contrast),
          annot=sprintf('annot[["%s"]]', ds)
        )
      }) |> stats::setNames(names(cntr[[ds]]))
    }) |> stats::setNames(names(cntr))

    observeEvent(input$contrast1, {
      contrast1(gsub(".*::", "", input$contrast1))
      dataset1(gsub("::.*", "", input$contrast1))
    })

    observeEvent(input$contrast2, {
      contrast2(gsub(".*::", "", input$contrast2))
      dataset2(gsub("::.*", "", input$contrast2))
    })

    output$matchby <- renderUI({
      tagList(
              column(width=6,
              selectInput(NS(id, "match1"),
                          "Match column 1:",
                          colnames(annot[[dataset1()]]),
                          selected=primary_id, width="100%")),
              column(width=6,
                     selectInput(NS(id, "match2"),
                          "Match column 2:",
                          colnames(annot[[dataset2()]]),
                          selected=primary_id, width="100%"))

              )
    })

    ## genes to indicate on the plot
    observeEvent(input$glabsgo, {
      if(isTruthy(input$glabs)) {
        ids <- strsplit(input$glabs, ",")[[1]]
        ids <- gsub("^[[:space:]]*", "", ids)
        ids <- gsub("[[:space:]]*$", "", ids)
        ids <- sprintf("^(%s)$",
                       paste0(ids, collapse='|'))
        #message("IDS are ", paste(ids, collapse="><"))
        .disco_log("manual label regex=", ids, ".")
        gene_labs(ids)
      } else {
        gene_labs(c())
      }
    })

    ## enable manual color scale
    observeEvent(input$autoscale, { 
      if(input$autoscale) { 
        disable("min") 
        disable("max")
      }
      else { 
        enable("min") 
        enable("max")
      }
    })

    observeEvent(input$show_top_labels, {
      if (isTRUE(input$show_top_labels)) {
        enable("top_label_n")
      } else {
        disable("top_label_n")
      }
    }, ignoreInit = FALSE)

    observeEvent(input$add_to_report, {
      .append_report_chunk(report, input$report_code)
    })

    ## save the disco plot to a PDF file
    output$save <- downloadHandler(
      filename = function() {
        sprintf(
          "disco_plot_%s_%s_vs_%s_%s.pdf",
          sanitize_filename(dataset1(), "dataset1"),
          sanitize_filename(contrast1(), "contrast1"),
          sanitize_filename(dataset2(), "dataset2"),
          sanitize_filename(contrast2(), "contrast2")
        )
      },
      content = function(file) {
        save_pdf(file=file, width=8, height=8, draw=function() {
          g <- plot_obj()
          print(g)
        })
      }
    )

    ## creating the actual plot
    observe({
      req(input$match1)
      req(input$match2)
      .disco_log("calculating plot.")

      .ds <- try(disco_score(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], 
                        annot1=annot[[dataset1()]],
                        annot2=annot[[dataset2()]],
                        primary_id=primary_id,
                        by=c(input$match1, input$match2)))

      disco(.ds)
      if(is(.ds, "try-error")) { 
        .disco_log("plot calculation failed.")
        plot_df(NULL)
       #output$discoplot <- renderPlot({
       #  stop(.ds)
       #})
       #
       # plot_obj(.ds)
        return(NULL) 
      }

      if(input$autoscale) {
        .lower <- -100
        .upper <- 100
      } else {
        .lower <- input$min
        .upper <- input$max
      }

      .plot_df <- .prep_disco_selection_df(.ds, lower=.lower, upper=.upper)
      plot_df(.plot_df)

      if(isTruthy(gene_labs)) { .glabs <- gene_labs() } else { .glabs <- NULL }

      if(input$autoscale) {
        g <- plot_disco(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], 
                        annot1=annot[[dataset1()]], 
                        annot2=annot[[dataset2()]], 
                        disco=disco(),
                        show_top_labels=if (isTRUE(input$show_top_labels)) input$top_label_n else 0,
                        label_sel=.glabs,
                        by=c(input$match1, input$match2))
      } else {
        g <- plot_disco(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], 
                        annot1=annot[[dataset1()]], 
                        annot2=annot[[dataset2()]], 
                        lower=input$min, upper=input$max, 
                        disco=disco(),
                        show_top_labels=if (isTRUE(input$show_top_labels)) input$top_label_n else 0,
                        label_sel=.glabs,
                        by=c(input$match1, input$match2))
      }

      if(dataset1() == "default") {
        .xlab <- contrast1()
      } else {
        .xlab <- paste0(dataset1(), ': ', contrast1())
      }

      if(dataset2() == "default") {
        .ylab <- contrast2()
      } else {
        .ylab <- paste0(dataset2(), ': ', contrast2())
      }

      g <- g + xlab(.xlab) + ylab(.ylab)

      .disco_log("storing plot object.")
      plot_obj(g)
    })

    ## Keep the report-code tab in sync with the current disco plot.
    ## The generated text mirrors the plot logic without executing it.
    observe({
      req(disco())
      req(dataset1(), dataset2(), contrast1(), contrast2(), input$match1, input$match2)

      if(input$autoscale) {
        .lower <- -100
        .upper <- 100
      } else {
        .lower <- input$min
        .upper <- input$max
      }

      .glabs <- if(isTruthy(gene_labs)) gene_labs() else NULL
      code <- plot_disco_chunk(
        cntr[[dataset1()]][[contrast1()]],
        cntr[[dataset2()]][[contrast2()]],
        annot1=annot[[dataset1()]],
        annot2=annot[[dataset2()]],
        lower=.lower,
        upper=.upper,
        show_top_labels=if(isTRUE(input$show_top_labels)) input$top_label_n else 0,
        label_sel=.glabs,
        by=c(input$match1, input$match2),
        primary_id=primary_id,
        env_map=list(
          contrast1=env_maps[[dataset1()]][[contrast1()]]$contrast,
          contrast2=env_maps[[dataset2()]][[contrast2()]]$contrast,
          annot1=env_maps[[dataset1()]][[contrast1()]]$annot,
          annot2=env_maps[[dataset2()]][[contrast2()]]$annot
        )
      )$code

      .xlab <- if(dataset1() == "default") contrast1() else paste0(dataset1(), ": ", contrast1())
      .ylab <- if(dataset2() == "default") contrast2() else paste0(dataset2(), ": ", contrast2())
      code <- code %+n% sprintf("g <- g + xlab(%s) + ylab(%s)", .r_code(.xlab), .r_code(.ylab))

      updateTextAreaInput(session, "report_code", value=code)
    })

    output$discoplot <- renderPlot({
      .disco_log("rendering plot.")
      .ds <- disco()
      if(is(.ds, "try-error")) { stop(.ds) }
      req(plot_obj())
      plot_obj()
    }, width=600, height=600, res=90)
    
    ## React to clicking on the plot: save the current list of genes as a
    ## table on the output, adding buttons for selecting a gene
    output$sel_genes <- renderTable({
      df <- req(selected_genes())
      #df <- isolate(current_genes())

      col_1 <- paste0(primary_id, "_1")
      col_2 <- paste0(primary_id, "_2")

      df[[col_1]] <- sprintf(as.character(link), dataset1(), df[[col_1]], df[[col_1]])
      df[[col_2]] <- sprintf(as.character(link), dataset2(), df[[col_2]], df[[col_2]])
      df
    }, sanitize.text.function=function(x) x)

    observeEvent(input$plot_click, {
      .pdf <- req(plot_df())
      pid <- nearPoints(.pdf, input$plot_click, xvar = "log2FoldChange.x", yvar = "log2FoldChange.y")
      ret <- .get_gene_df(pid, selcols, primary_id, annot[[dataset1()]], annot[[dataset2()]])
      selected_genes(ret)
    })

    ## react to hover over points: enter the close genes into current list
    observeEvent(input$plot_hover, {
      .pdf <- req(plot_df())
      pid <- nearPoints(.pdf, input$plot_hover, xvar = "log2FoldChange.x", yvar = "log2FoldChange.y")
      ret <- .get_gene_df(pid, selcols, primary_id, annot[[dataset1()]], annot[[dataset2()]])
      current_genes(ret)
    })

    ## react to points selected by brush: enter the genes into current list
    observeEvent(input$plot_brush, {
      .pdf <- req(plot_df())
      pid <- brushedPoints(.pdf, input$plot_brush, xvar = "log2FoldChange.x", yvar = "log2FoldChange.y")
      ret <- .get_gene_df(pid, selcols, primary_id, annot[[dataset1()]], annot[[dataset2()]])
      selected_genes(ret)
    })

    ## enter current genes into the output table
    output$point_id <- renderTable({ 
      .cg <- current_genes()
      if(!is.null(.cg) && nrow(.cg) > 0L) {
        return(.cg[1, ])
      } else {
        return(NULL)
      }
    })


    observeEvent(input$genebutton, {
      ids <- strsplit(input$genebutton, '~')[[1]]
      gene_id$ds <- ids[2]
      gene_id$id <- ids[3]
    })

  })
}
