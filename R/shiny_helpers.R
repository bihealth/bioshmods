#' Figure size input
#' @param id ID of the element
#' @param label Label for the input
#' @param choices Which figure size choices (in pixels)
#' @param custom Whether custom fig sizes are allowed
#' @param selected Default figure size
#' @param ... further arguments passed to selectizeInput
#' @export
figsizeInput <- function(id, label="Figure size", 
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




