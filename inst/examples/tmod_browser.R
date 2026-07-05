## extending the example from tmodBrowserTableServer
if(interactive()) {

  library(shiny)

  ui <- fluidPage(
   fluidRow(tmodBrowserTableUI("tt", list(default=names(C19_gs$tmod_res)), upset=TRUE)),
   fluidRow(tmodBrowserPlotUI("tp"))
   )

  server <- function(input, output) {
    gs_id <- reactiveValues()
    tmodBrowserTableServer("tt", list(default=C19_gs$tmod_res), gs_id = gs_id,
                                 tmod_dbs = list(default=C19_gs$tmod_dbs))
    tmodBrowserPlotServer("tp",
             gs_id=gs_id,
             tmod_map=list(default=C19_gs$tmod_map),
             tmod_dbs=list(default=C19_gs$tmod_dbs),
             tmod_gl=list(default=C19_gs$tmod_gl),
             cntr=list(default=C19$contrasts),
             annot=list(default=C19$annotation))
  }
  runApp(shinyApp(ui, server))
}

