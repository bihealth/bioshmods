## extending the example from tmodBrowserTableServer
if(interactive()) {

  library(shiny)

  ui <- fluidPage(
   fluidRow(tmodBrowserTableUI("tt", names(C19_gs$tmod_res), upset=TRUE)),
   fluidRow(tmodBrowserPlotUI("tp"))
   )

  server <- function(input, output) {
    gs_id <- reactiveValues()
    tmodBrowserTableServer("tt", C19_gs$tmod_res, gs_id = gs_id,
                                 tmod_dbs = C19_gs$tmod_dbs)
    tmodBrowserPlotServer("tp",
             gs_id=gs_id,
             tmod_map=C19_gs$tmod_map,
             tmod_dbs=C19_gs$tmod_dbs,
             tmod_gl=C19_gs$tmod_gl,
             cntr=C19$contrasts,
             annot=C19$annotation)
  }
  runApp(shinyApp(ui, server))
}


