.compile_report <- function(report_var, template) {
  if(is.null(report_var)) {
    msg("No report_var")
    return("")
  }

  if(is.null(report_var$chunks)) {
    msg("No chunks yet")
    return("")
  }
  msg("Found ", length(report_var$chunks), " chunks\n")

  if(length(report_var$chunks) > 0) {
    enable("Download")
  }

  chunktext <- map_chr(report_var$chunks, ~ .x$chunk)
  chunktext <- paste0(chunktext, collapse="\n\n")

  ret <- paste0(template, "\n\n", chunktext)

  return(ret)
}

#' @rdname reportServer
#' @export
reportUI <- function(id, contrasts=FALSE) {
  ns <- NS(id)
  save_ui <- 
    sidebarPanel(

      # select format type
      selectInput(ns("format"), "Format", choices=c("All (ZIP file)"="zip", 
                                                        "PDF"="pdf", 
                                                        "Rmarkdown"="rmd", 
                                                        "HTML"="html", 
                                                        "Word document"="docx"
      ), selected="zip"),
      tipify(downloadButton(ns("Download"), "Download", class="bg-success"), "Download report"))

    return(sidebarLayout(save_ui,
      mainPanel(
                tabsetPanel(
                            tabPanel("Summary", uiOutput(ns("report_summary"))),
                            tabPanel("Markdown", uiOutput(ns("report_markdown"))),
                            tabPanel("Preview", uiOutput(ns("report_preview")))
                            )
    )))
}


#' Shiny Module – Report generator
#'
#' Shiny Module – generating reports
#'
#' @param id id of the module
#' @param template Rmarkdown or Quarto template to use
#' @param report_var reactive values object holding the report chunks
#' @param generator rmarkdown or quarto
#' @return does not return anything useful
#' @examples
#' if(interactive()) {
#'    ui  <- fluidPage(geneBrowserPlotUI("gplot", FALSE))
#'    serv <- function(input, output, session) {
#'      geneBrowserPlotServer("gplot", list(id="MUZG"), covar, mtx)
#'    }
#'    shinyApp(ui, serv)
#' }
#' @importFrom zip zip
#' @export
reportServer <- function(id, template, report_var, generator="rmarkdown") {
  parent_frame <- parent.frame()


  # start the module server
  moduleServer(id, function(input, output, session) {

    disable("Download")
    rmd_text <- reactiveVal("")

    observe({
      msg("Compiling report")
      rmd_text(.compile_report(report_var, template))
      msg("RMD text: ", rmd_text())
    })

    output$report_summary <- renderUI({
      if(is.null(report_var) || length(report_var$chunks) == 0) return(NULL)
      .generate_report_summary(report_var)
    })

    output$report_markdown <- renderUI({
      if(is.null(report_var) || length(report_var$chunks) == 0) return(NULL)
      tags$pre(rmd_text())
    })

    output$report_preview <- renderUI({
      if(is.null(report_var) || length(report_var$chunks) == 0) return(NULL)
      .generate_report_preview(report_var, parent_frame)
    })

    ## Save figure as a PDF
    output$Download <- downloadHandler(
      filename = function() return(paste0("report.", input$format)),
      content = function(file) .generate_report(file, rmd_text(), input$format)
    )
  
  }) # end of moduleServer
}

## creates a html with the summary of the report contents
.generate_report_summary <- function(report_var) {
  tags$ul(
    map(report_var$chunks, ~ tags$li(HTML(sprintf("%s [%s]: %s",
                                                  .x$module_id, 
                                                  .x$type, 
                                                  .x$title))))
  )
}

## creates a html with the preview of the report contents
.generate_report_preview <- function(report_var, parent_frame) {
  ## run setup code
  map(report_var$chunks, ~ {
    if(.x$type == "setup") {
      eval(parse(text=.x$code), envir=parent_frame)
    } 
  })

  ## generate all plots
  tags$div(
    map(report_var$chunks, ~ {
      if(.x$type == "plot") {
        tags$div(
          tags$h2(paste("Figure:", .x$title)),
          renderPlot({ print(eval(parse(text=.x$code), envir=parent_frame)) },
            width=.x$fig.width, height=.x$fig.height)
        )
      } 
  }))
}


.generate_report <- function(file, contents, type, quiet=FALSE, rmdfile=NULL) {
  msg("Generating report ", file, " of type ", type)
  if(type == "rmd") {
    cat(contents, file=file)
    return()
  }

  if(!quiet) {
    showModal(modalDialog(
      title = "Generating report",
      "Please wait while the report is being generated...",
      footer = NULL,
      easyClose = TRUE,
      size = "m"
    ))
  }

  if(type == "zip") {
    tmpdir <- tempdir()
    msg("Created temporary directory", tmpdir)
    rmdfile <- file.path(tmpdir, "report.rmd")
    cat(contents, file=rmdfile)
    for(fmt in c("pdf", "docx", "html")) {
      .generate_report(file=file.path(tmpdir, paste0("report.", fmt)), 
                       contents=contents, type=fmt, quiet=TRUE, rmdfile=rmdfile)
    }
    msg("Zipping ", tmpdir, " to ", file)
    zip(file, paste0("report.", c("rmd", "pdf", "docx", "html")), root=tmpdir)
    if(!quiet) removeModal()
    return()
  }

  if(is.null(rmdfile)) {
    rmdfile <- tempfile(fileext = ".Rmd")
    cat(contents, file=rmdfile)
  }

  if(type == "docx") type <- "word_document"
  output_format <- match.arg(type, c("pdf_document", "word_document", "html_document"))

  msg("generating report ", rmdfile, " to ", file, " in format ", output_format)
  result <- try(rmarkdown::render(rmdfile, 
                                  output_file=file, 
                                  quiet=TRUE, 
                                  output_format=output_format), 
                silent=TRUE)

  if(class(result) == "try-error") {
    stop("Error generating report: ", result)
    if(!quiet) removeModal()
    return()
  }

  if(!quiet) removeModal()
}
