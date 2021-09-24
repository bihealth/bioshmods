#' @rdname discoServer
#' @export
discoUI <- function(id, cntr_titles) {

  cntr_titles <- .prep_cntr_titles(cntr_titles)
  cntr_flat   <- unlist(cntr_titles, recursive=FALSE)

  if(!length(cntr_flat) > 1) {
    h4("You need at least two contrasts for this plot")
  } else {

  fluidRow(
    column(width=1),
    column(width=2,
        fluidRow(selectInput(NS(id, "contrast1"), label = "Contrast 1", 
                             choices = cntr_titles, width="100%"),
        selectInput(NS(id, "contrast2"), label = "Contrast 2", 
                             choices = cntr_titles, selected=cntr_flat[2], width="100%")),
        fluidRow(checkboxInput(NS(id, "autoscale"), "Automatic scale", value=TRUE)),
        fluidRow(sliderInput(NS(id, "min"), "Min", min=-150, max=0, value=-100, width="80%")),
        fluidRow(sliderInput(NS(id, "max"), "Max", min=0, max=150, value=100, width="80%")),
        fluidRow(downloadButton(NS(id, "save"), "Save plot to PDF", class="bg-success")),
        HTML("<br/>Hover to identify genes, click to select, or click & drag to select an area<br/><br/>"),
        fluidRow(verbatimTextOutput(NS(id, "msg"))),
        fluidRow(tableOutput(NS(id, "point_id")))
    ),

    column(width=6, 
      plotOutput(NS(id, "discoplot"), 
                 hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                 click=NS(id, "plot_click"),
                 brush=NS(id, "plot_brush")

      )
    ),
    column(width=2,
      HTML("Click on the button to view an expression profile"),
      tableOutput(NS(id, "sel_genes"))
      )
  )
  }
}

.get_gene_df <- function(pid, selcols, primary_id="PrimaryID", annot=NULL) {

  if(is.null(annot)) {
    ret <- data.frame(PrimaryID=pid)
  } else {
    ret <- annot %>% filter(.data[[primary_id]] %in% pid) %>%
      select(any_of(selcols))
  }
  return(ret)
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
#' @examples
#' if(interactive()) {
#'    cntr1 <- data.frame(log2FoldChange=rnorm(5000),
#'                        pvalue=runif(5000))
#'    rownames(cntr1) <- paste0("ID", 1:5000)
#'    cntr2 <- data.frame(log2FoldChange=cntr1$log2FoldChange + 
#'                                       rnorm(5000),
#'                        pvalue=runif(5000) * cntr1$pvalue)
#'    rownames(cntr2) <- paste0("ID", 1:5000)
#'    cntr <- list("Contrast 1"=cntr1, "Contrast 2"=cntr2)
#'    shinyApp(ui=fluidPage(discoUI("disco", names(cntr))),
#'             server=function(input, output, session) {
#'                discoServer("disco", cntr)
#'             })
#' }
#' @export
discoServer <- function(id, cntr, annot=NULL,
    selcols=c("PrimaryID", "ENTREZ", "SYMBOL"),
    primary_id="PrimaryID", gene_id=NULL) {

  if(!is.data.frame(cntr[[1]])) {
    message("discoServer in multilevel mode")
  } else {
    cntr  <- list(default=cntr)
    annot <- list(default=annot)
  }


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

    observeEvent(input$contrast1, {
      contrast1(gsub(".*::", "", input$contrast1))
      dataset1(gsub("::.*", "", input$contrast1))
    })
    observeEvent(input$contrast2, {
      contrast2(gsub(".*::", "", input$contrast2))
      dataset2(gsub("::.*", "", input$contrast2))
    })


    selcols <- c("PrimaryID", "ENTREZ", "SYMBOL")

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
        if(input$autoscale) {
          g <- plot_disco(cntr[[dataset1()]][[contrast1()]], 
                          cntr[[dataset2()]][[contrast2()]], disco=disco())
        } else {
          g <- plot_disco(cntr[[dataset1()]][[contrast1()]], 
                          cntr[[dataset2()]][[contrast2()]], 
                          lower=input$min, upper=input$max, disco=disco())
        }
        print(g)
        dev.off()
      }
    )

    ## creating the actual plot
    output$discoplot <- renderPlot({
      disco(disco_score(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], by=primary_id))

      if(input$autoscale) {
        g <- plot_disco(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], disco=disco())
      } else {
        g <- plot_disco(cntr[[dataset1()]][[contrast1()]], 
                        cntr[[dataset2()]][[contrast2()]], 
                        lower=input$min, upper=input$max, disco=disco())
      }
      return(g)
    }, width=600, height=600, res=90)
    
    ## React to clicking on the plot: save the current list of genes as a
    ## table on the output, adding buttons for selecting a gene
    output$sel_genes <- renderTable({
      df <- req(selected_genes())
      #df <- isolate(current_genes())
      link <- actionButton(NS(id, "gene_id~%s~%s"), label="%s \U25B6 ",
                           onclick=sprintf('Shiny.onInputChange(\"%s-genebutton\",  this.id)', id),
                           class = "btn-primary btn-sm")
      df[[primary_id]] <- sprintf(as.character(link), dataset1(), df[[primary_id]], df[[primary_id]])
      df
    }, sanitize.text.function=function(x) x)

    observeEvent(input$plot_click, {
      selected_genes(current_genes())
    })

    ## react to hover over points: enter the close genes into current list
    observeEvent(input$plot_hover, {
      pid <- disco() %>% nearPoints(input$plot_hover) %>% pull(.data[[primary_id]])
      ret <- .get_gene_df(pid, selcols, primary_id, annot[[dataset1()]])
      current_genes(ret)
    })

    ## react to points selected by brush: enter the genes into current list
    observeEvent(input$plot_brush, {
      pid <- disco() %>% brushedPoints(input$plot_brush) %>% pull(.data[[primary_id]])
      ret <- .get_gene_df(pid, selcols, primary_id, annot[[dataset1()]])
      selected_genes(ret)
    })

    ## enter current genes into the output table
    output$point_id <- renderTable({ 
      current_genes()
    })


    observeEvent(input$genebutton, {
      ids <- strsplit(input$genebutton, '~')[[1]]
      gene_id$ds <- ids[2]
      gene_id$id <- ids[3]
    })

  })
}

