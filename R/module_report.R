# wrap the chunk code in newlines and chunk bounds
.report_module_wrap_chunk <- function(code) {
  paste0("\n\n```{r}\n", code, "\n```\n")
}

# chunks can be a reactive var, so we need to handle it accordingly
.report_module_read_chunks <- function(chunks) {
  if(is.reactive(chunks)) {
    return(chunks())
  }

  chunks
}

# here we build the rmd from the indiv. chunks
.report_module_build_rmd <- function(template, preamble=NULL, chunks=list()) {
  out_file <- tempfile(fileext=".Rmd")
  file.copy(template, out_file, overwrite=TRUE)

  if(!is.null(preamble) && nzchar(preamble)) {
    cat(.report_module_wrap_chunk(preamble), file=out_file, append=TRUE)
  }

  for(chunk in chunks) {
    cat(.report_module_wrap_chunk(chunk), file=out_file, append=TRUE)
  }

  out_file
}

#' @rdname reportModuleServer
#' @export
reportModuleUI <- function(id) {
  ns <- NS(id)

  fluidRow(
    shinyjs::useShinyjs(),
    column(
      width=3,
      fluidRow(column(
        width=12,
        actionButton(ns("generate_report"), "Generate report", class="btn-primary")
      )),
      fluidRow(column(
        width=12,
        downloadButton(ns("download_rmd"), "Download R Markdown", class="btn-default")
      ))
    ),
    column(
      width=9,
      tabsetPanel(
        tabPanel(
          "Report",
          textAreaInput(
            ns("report_text"),
            label=NULL,
            value="",
            width="100%",
            height="420px"
          )
        ),
        tabPanel(
          "Log",
          textAreaInput(
            ns("messages"),
            label=NULL,
            value="",
            width="100%",
            height="420px"
          )
        )
      )
    )
  )
}

#' Shiny module for generating a simple R Markdown report
#'
#' Shiny module for generating a simple R Markdown report
#'
#' @param id Module identifier.
#' @param chunks List of report chunks, or a reactive expression returning
#'   such a list.
#' @param preamble Optional code string written as the first R chunk.
#' @param template File name of the R Markdown template to copy and extend.
#'
#' @return Invisibly returns the generated report path as a reactive value.
#'
#' @examples
#' template <- tempfile(fileext=".Rmd")
#' writeLines(c("---", "title: Example report", "---", ""), template)
#'
#' if(interactive()) {
#'   chunks <- list(
#'     "plot(1:10)",
#'     "summary(cars)"
#'   )
#'
#'   ui <- fluidPage(reportModuleUI("report"))
#'   server <- function(input, output, session) {
#'     reportModuleServer(
#'       "report",
#'       chunks=chunks,
#'       preamble="library(ggplot2)",
#'       template=template
#'     )
#'   }
#'
#'   shinyApp(ui, server)
#' }
#'
#' if(interactive()) {
#'   data(C19)
#'
#'   gene_id <- reactiveValues(id=C19$annotation$PrimaryID[1], ds="default")
#'   report <- reactiveValues(chunks=list())
#'
#'   ui <- fluidPage(
#'     tabsetPanel(
#'       tabPanel("Gene plot", geneBrowserPlotUI("gplot")),
#'       tabPanel("Report", reportModuleUI("report"))
#'     )
#'   )
#'
#'   server <- function(input, output, session) {
#'     geneBrowserPlotServer(
#'       "gplot",
#'       gene_id=gene_id,
#'       covar=list(default=C19$covariates),
#'       exprs=list(default=C19$expression),
#'       annot=list(default=C19$annotation),
#'       report=report
#'     )
#'
#'     reportModuleServer(
#'       "report",
#'       chunks=reactive(report$chunks),
#'       preamble="data(C19)",
#'       template=template
#'     )
#'   }
#'
#'   shinyApp(ui, server)
#' }
#' @export
reportModuleServer <- function(id, chunks, preamble=NULL, template) {
  moduleServer(id, function(input, output, session) {
    report_file <- reactiveVal(NULL)

    shinyjs::disable("download_rmd")

    observeEvent(input$generate_report, {
      messages <- character(0)
      report_text <- ""

      tryCatch({
        current_chunks <- .report_module_read_chunks(chunks)
        if(is.null(current_chunks)) {
          current_chunks <- list()
        }

        messages <- c(messages, sprintf("Copying template: %s", template))
        out_file <- .report_module_build_rmd(
          template=template,
          preamble=preamble,
          chunks=current_chunks
        )

        report_file(out_file)
        shinyjs::enable("download_rmd")
        report_text <- paste(readLines(out_file, warn=FALSE), collapse="\n")
        messages <- c(
          messages,
          sprintf("Appended preamble: %s", if(is.null(preamble) || !nzchar(preamble)) "no" else "yes"),
          sprintf("Appended chunks: %d", length(current_chunks)),
          sprintf("Report generated: %s", out_file)
        )
      }, error=function(e) {
        report_file(NULL)
        shinyjs::disable("download_rmd")
        report_text <<- ""
        messages <<- c(messages, sprintf("Error: %s", conditionMessage(e)))
      })

      updateTextAreaInput(session, "report_text", value=report_text)
      updateTextAreaInput(session, "messages", value=paste(messages, collapse="\n"))
    })

    output$download_rmd <- downloadHandler(
      filename=function() {
        "report.Rmd"
      },
      content=function(file) {
        req(report_file())
        file.copy(report_file(), file, overwrite=TRUE)
      }
    )

    invisible(report_file)
  })
}
