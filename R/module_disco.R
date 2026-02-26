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

    column(width=5, 
      withSpinner(plotOutput(NS(id, "discoplot"), 
                 hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                 click=NS(id, "plot_click"),
                 brush=NS(id, "plot_brush"),
                 height="600px")
      ),
      HTML("<br/>Hover to identify genes, click to select, or click & drag to select an area<br/><br/>"),
      fluidRow(tableOutput(NS(id, "point_id")))
    ),
    column(width=4,
      HTML("Click on the button to view an expression profile"),
      tableOutput(NS(id, "sel_genes"))
      )
  )
  }
}

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

.prep_disco_selection_df <- function(df, lower, upper) {
  ret <- df %>%
    filter(!is.na(.data$log2FoldChange.x) & !is.na(.data$log2FoldChange.y) & !is.na(.data$disco)) %>%
    mutate(disco=ifelse(.data$disco > upper, upper, ifelse(.data$disco < lower, lower, .data$disco))) %>%
    arrange(abs(.data$disco))
  print(ret)
  ret
}

.is_list_of_data_frames <- function(x) {
  is.list(x) && length(x) > 0L && all(vapply(x, is.data.frame, logical(1)))
}

.is_list_of_list_of_data_frames <- function(x) {
  is.list(x) && length(x) > 0L &&
    all(vapply(x, function(ds) .is_list_of_data_frames(ds), logical(1)))
}

.cntr_primary_ids <- function(cntr_ds, primary_id) {
  unique(unlist(lapply(cntr_ds, function(df) as.character(df[[primary_id]])), use.names = FALSE))
}

.make_annotation_from_cntr <- function(cntr_ds, primary_id) {
  data.frame(
    stats::setNames(list(.cntr_primary_ids(cntr_ds, primary_id)), primary_id),
    check.names = FALSE
  )
}

.normalize_disco_inputs <- function(cntr, annot, primary_id) {
  if (.is_list_of_data_frames(cntr)) {
    cntr <- list(default=cntr)
    if (is.null(annot) || is.data.frame(annot)) {
      annot <- list(default=annot)
    } else if (is.list(annot) && length(annot) == 1L &&
               (is.null(annot[[1]]) || is.data.frame(annot[[1]]))) {
      names(annot) <- "default"
    } else {
      stop("For single-dataset cntr, annot must be NULL, a data frame, or a single-element list.")
    }
  } else if (.is_list_of_list_of_data_frames(cntr)) {
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

    if (!all(vapply(cntr_ds, function(df) primary_id %in% colnames(df), logical(1)))) {
      bad <- names(cntr_ds)[!vapply(cntr_ds, function(df) primary_id %in% colnames(df), logical(1))]
      stop(sprintf(
        "In dataset '%s', all contrasts must contain the '%s' column. Missing in: %s",
        ds, primary_id, paste(bad, collapse = ", ")
      ))
    }

    if (is.null(annot[[ds]])) {
      annot[[ds]] <- .make_annotation_from_cntr(cntr_ds, primary_id)
    }

    if (!is.data.frame(annot[[ds]])) {
      stop(sprintf("Annotation for dataset '%s' must be a data frame or NULL.", ds))
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
#' @param cntr list of data frames containing the contrast information.
#'        Data frames must have the columns log2FoldChange and pvalue. Rownames of
#'        the data frames should be IDs.
#' @param annot Annotation data frame. The annotation data frame must have
#'        a column named "PrimaryID" which corresponds to the rownames of the data
#'        frames in the `cntr` list.
#' @param selcols which column in the gene table when genes are selected
#'        from the plot
#' @param cntr_titles character vector containing the IDs of the contrasts
#'        (same as `names(cntr)`).
#' @param gene_id must be a `reactiveValues` object. If not NULL, then
#' clicking on a gene identifier will modify this object (possibly
#' triggering an event in another module).
#' @return Returns a reactive expression returning the ID of the activated gene
#' @example inst/examples/disco.R
#' @export
discoServer <- function(id, cntr, annot=NULL,
    selcols=c("PrimaryID", "ENTREZ", "SYMBOL"),
    primary_id="PrimaryID", gene_id=NULL) {

  normalized <- .normalize_disco_inputs(cntr, annot, primary_id)
  cntr <- normalized$cntr
  annot <- normalized$annot

  if(!"default" %in% names(cntr)) {
    message("discoServer in multi dataset mode")
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
        message("IDS regex: ", ids)
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

    ## save the disco plot to a PDF file
    output$save <- downloadHandler(
      filename = function() {
        sprintf(
          "disco_plot_%s_%s_vs_%s_%s.pdf",
          .sanitize_filename(dataset1(), "dataset1"),
          .sanitize_filename(contrast1(), "contrast1"),
          .sanitize_filename(dataset2(), "dataset2"),
          .sanitize_filename(contrast2(), "contrast2")
        )
      },
      content = function(file) {
        pdf(file=file, width=8, height=8)
        g <- plot_obj()
        print(g)
        dev.off()
      }
    )

    ## creating the actual plot
    observe({
      req(input$match1)
      req(input$match2)
      message("calculating plot")

      .ds <- try(disco_score(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], 
                        annot1=annot[[dataset1()]],
                        annot2=annot[[dataset2()]],
                        primary_id=primary_id,
                        by=c(input$match1, input$match2)))

      disco(.ds)
      if(is(.ds, "try-error")) { 
        message("An error occured")
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

      message("storing plot")
      plot_obj(g)
    })

    output$discoplot <- renderPlot({
      message("rendering plot")
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
