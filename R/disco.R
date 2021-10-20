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

  if(is.null(annot1)) {
    pid1 <- data.frame(pid1)
    colnames(pid1) <- primary_id
  } else {
    pid1 <- annot1 %>% slice(match(pid1, .data[[primary_id]])) %>%
      select(any_of(selcols))
  }

  if(is.null(annot2)) {
    pid2 <- data.frame(pid2)
    colnames(pid2) <- primary_id
  } else {
    pid2 <- annot2 %>% slice(match(pid2, .data[[primary_id]])) %>%
      select(any_of(selcols))
  }

  colnames(pid1) <- paste0(colnames(pid1), "_1")
  colnames(pid2) <- paste0(colnames(pid2), "_2")

  return(cbind(pid1, pid2))
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

  if(!is.data.frame(cntr[[1]])) {
    message("discoServer in multi dataset mode")
  } else {
    cntr  <- list(default=cntr)
    annot <- list(default=annot)
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

    ## save the disco plot to a PDF file
    output$save <- downloadHandler(
      filename = function() {
        ret <- sprintf("disco_plot_%s_%s_vs_%s_%s.pdf", 
                       dataset1(), contrast1(),
                       dataset2(), contrast2())
        ret <- gsub("[^0-9a-zA-Z_.-]", "", ret)
        return(ret)
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
      if(class(.ds) == "try-error") { 
        message("An error occured")
       #output$discoplot <- renderPlot({
       #  stop(.ds)
       #})
       #
       # plot_obj(.ds)
        return(NULL) 
      }

      if(isTruthy(gene_labs)) { .glabs <- gene_labs() } else { .glabs <- NULL }

      if(input$autoscale) {
        g <- plot_disco(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], 
                        annot1=annot[[dataset1()]], 
                        annot2=annot[[dataset2()]], 
                        disco=disco(), label_sel=.glabs,
                        by=c(input$match1, input$match2))
      } else {
        g <- plot_disco(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], 
                        annot1=annot[[dataset1()]], 
                        annot2=annot[[dataset2()]], 
                        lower=input$min, upper=input$max, 
                        disco=disco(), label_sel=.glabs,
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
      if(class(.ds) == "try-error") { stop(.ds) }
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
      selected_genes(current_genes())
    })

    ## react to hover over points: enter the close genes into current list
    observeEvent(input$plot_hover, {
      pid <- disco() %>% nearPoints(input$plot_hover) #%>% pull(.data[[primary_id]])
      ret <- .get_gene_df(pid, selcols, primary_id, annot[[dataset1()]], annot[[dataset2()]])
      current_genes(ret)
    })

    ## react to points selected by brush: enter the genes into current list
    observeEvent(input$plot_brush, {
      pid <- disco() %>% brushedPoints(input$plot_brush) #%>% pull(.data[[primary_id]])
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

