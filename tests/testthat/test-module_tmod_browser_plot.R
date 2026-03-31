library(shiny)

test_that("tmodBrowserPlotUI contains show placeholder next to save button", {
  ui <- tmodBrowserPlotUI("tp")
  html <- as.character(ui)

  expect_s3_class(ui, "shiny.tag")
  expect_true(grepl("tp-show_selected_ui", html, fixed = TRUE))
  expect_true(grepl("tp-save", html, fixed = TRUE))
})

test_that("tmodBrowserPlotServer updates selection ids on Show button click", {
  gs_id <- reactiveValues(ds = "default", id = "set1", cntr = "contrast_a", db = "db1", sort = "pval")
  selection <- reactiveValues(ids = character(0))

  tmod_dbs <- list(
    db1 = list(
      MODULES2GENES = list(set1 = c("db_g2", "db_g1"))
    )
  )
  cntr <- list(
    contrast_a = data.frame(
      PrimaryID = c("g1", "g2", "g3"),
      pvalue = c(0.02, 0.01, 0.5),
      padj = c(0.03, 0.02, 0.6),
      log2FoldChange = c(1.2, -1.5, 0.1),
      stringsAsFactors = FALSE
    )
  )
  tmod_map <- list(
    dbs = c(db1 = "map1"),
    maps = list(map1 = c(g1 = "db_g1", g2 = "db_g2", g3 = "db_g3"))
  )
  annot <- data.frame(
    PrimaryID = c("g1", "g2", "g3"),
    SYMBOL = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  testServer(
    tmodBrowserPlotServer,
    args = list(
      gs_id = gs_id,
      tmod_dbs = tmod_dbs,
      cntr = cntr,
      tmod_map = tmod_map,
      annot = annot,
      selection = selection,
      ui_config = list(show_button_label = "Send")
    ),
    {
      session$flushReact()

      expect_true(any(grepl("Send", as.character(output$show_selected_ui), fixed = TRUE)))

      session$setInputs(show_selected = 1)
      session$flushReact()
    }
  )

  expect_equal(isolate(selection$ids), c("g2", "g1"))
})

test_that("tmodBrowserPlotServer validates selection and ui_config", {
  gs_id <- reactiveValues(ds = "default", id = "set1", cntr = "contrast_a", db = "db1", sort = "pval")

  tmod_dbs <- list(
    db1 = list(
      MODULES2GENES = list(set1 = c("db_g1"))
    )
  )
  cntr <- list(
    contrast_a = data.frame(
      PrimaryID = "g1",
      pvalue = 0.01,
      padj = 0.02,
      log2FoldChange = 1,
      stringsAsFactors = FALSE
    )
  )
  tmod_map <- list(
    dbs = c(db1 = "map1"),
    maps = list(map1 = c(g1 = "db_g1"))
  )
  annot <- data.frame(PrimaryID = "g1", stringsAsFactors = FALSE)

  expect_error(
    testServer(
      tmodBrowserPlotServer,
      args = list(
        gs_id = gs_id,
        tmod_dbs = tmod_dbs,
        cntr = cntr,
        tmod_map = tmod_map,
        annot = annot,
        selection = character(0)
      ),
      {}
    ),
    "selection"
  )

  expect_error(
    testServer(
      tmodBrowserPlotServer,
      args = list(
        gs_id = gs_id,
        tmod_dbs = tmod_dbs,
        cntr = cntr,
        tmod_map = tmod_map,
        annot = annot,
        ui_config = list(show_button_label = "")
      ),
      {}
    ),
    "show_button_label"
  )
})
