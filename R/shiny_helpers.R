#' Figure size input
#' @param id ID of the element
#' @param label Label for the input
#' @param choices Which figure size choices (in pixels)
#' @param custom Whether custom fig sizes are allowed
#' @param selected Default figure size
#' @param ... further arguments passed to selectizeInput
#' @export
figsizeInput <- function(id, label="Figure size (w x h)", 
                             choices=c("800x800", 
                                          "600x600",
                                          "600x800",
                                          "800x600",
                                          "1200x600",
                                          "1200x800",
                                          "600x1200",
                                          "800x1200",
                                          "1200x1200"),
                         selected=NULL,
                         custom=TRUE, ...) {

  return(selectizeInput(id, label=label, choices=choices, selected=selected,
                        options=list(create=custom, plugins=list('restore_on_backspace')), ...))

}


# Build a blank report-code tab shared by plot modules.
# The button is intentionally not wired until report storage is defined.
.report_code_tab <- function(id, input_id="report_code", button_id="add_to_report", title="Code") {
  shiny::tabPanel(
    title,
    shiny::fluidRow(
      shiny::br(),
      shiny::textAreaInput(
        shiny::NS(id, input_id),
        label=NULL,
        value="",
        width="100%",
        height="520px"
      ),
      shiny::actionButton(shiny::NS(id, button_id), "Add to report", class="btn-default")
    )
  )
}



