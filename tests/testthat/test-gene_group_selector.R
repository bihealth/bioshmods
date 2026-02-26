library(shiny)

test_that("geneGroupSelectorUI returns controls-only UI", {
  ui <- geneGroupSelectorUI("gsel")
  html <- as.character(ui)

  expect_s3_class(ui, "shiny.tag.list")
  expect_false(grepl("sidebar-layout", html, fixed = TRUE))
  expect_false(grepl("selected_annotation", html, fixed = TRUE))
})

test_that("geneGroupSelectorServer populates external reactiveVal with PrimaryIDs", {
  annot <- data.frame(
    PrimaryID = c("g1", "g2", "g3"),
    SYMBOL = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  selected_ids <- reactiveVal(character(0))

  testServer(
    geneGroupSelectorServer,
    args = list(
      annot = annot,
      selected_ids = selected_ids
    ),
    {
      session$setInputs(modus = "by_name")
      session$setInputs(name_id_col = "SYMBOL")
      session$setInputs(name_list = "B, A, missing")
      session$flushReact()
    }
  )

  expect_equal(isolate(selected_ids()), c("g2", "g1"))

  testServer(
    geneGroupSelectorServer,
    args = list(
      annot = annot,
      selected_ids = selected_ids
    ),
    {
      session$setInputs(modus = "by_name")
      session$setInputs(name_id_col = "PrimaryID")
      session$setInputs(name_list = "")
      session$flushReact()
    }
  )

  expect_equal(isolate(selected_ids()), character(0))
})
