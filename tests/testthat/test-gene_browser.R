library(shiny)

check_serverrun <- function(cntr, cntr_names=NULL, annot=NULL) {

  expect_error(geneBrowserTableUI(cntr_names))
  expect_s3_class(
                  (ui <- geneBrowserTableUI("gb", cntr_names)),
                  "shiny.tag")
  server <- function(input, output) {
    geneBrowserTableServer("gb", 
                           cntr =cntr,
                           annot=annot)
  }

  expect_type(server, "closure")

  expect_s3_class(
     app <- shinyApp(ui, server),
     "shiny.appobj")
  ## check that server can be started
  testServer(geneBrowserTableServer, args=list(cntr=cntr, annot=annot),
             {

             })
}

test_that("geneBrowserTable can be build", {

  check_serverrun(C19$contrasts, names(C19$contrasts), C19$annotation)
  cntr <- list(testing=C19$contrasts)
  annot <- list(testing=C19$annotation)
  cntr_names <- list(testing=names(C19$contrasts))
  check_serverrun(cntr, cntr_names, annot)

})
