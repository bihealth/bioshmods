library(shiny)

.heatmap_test_inputs <- function(use_label = FALSE) {
  annot <- data.frame(
    PrimaryID = paste0("g", 1:6),
    SYMBOL = LETTERS[1:6],
    stringsAsFactors = FALSE
  )

  set.seed(1)
  exprs <- matrix(
    rnorm(24),
    nrow = 6,
    dimnames = list(annot$PrimaryID, paste0("s", 1:4))
  )

  if(use_label) {
    covar <- data.frame(
      label = colnames(exprs),
      group = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    )
  } else {
    covar <- data.frame(
      SampleID = colnames(exprs),
      group = c("A", "A", "B", "B"),
      stringsAsFactors = FALSE
    )
  }

  list(annot = annot, exprs = exprs, covar = covar)
}

test_that("heatmapUI builds sidebar layout with selector and options", {
  ui <- heatmapUI("hm")
  html <- as.character(ui)

  expect_s3_class(ui, "shiny.tag")
  expect_true(grepl("col-sm-4", html, fixed = TRUE))
  expect_true(grepl("col-sm-8", html, fixed = TRUE))
  expect_true(grepl("hm-gene_selector-modus", html, fixed = TRUE))
  expect_true(grepl("hm-show_legend", html, fixed = TRUE))
  expect_true(grepl("hm-save", html, fixed = TRUE))
  expect_true(grepl("hm-heatmap_plot", html, fixed = TRUE))
})

test_that("heatmapServer returns reactives and computes heatmap for selected genes", {
  x <- .heatmap_test_inputs()

  testServer(
    heatmapServer,
    args = list(
      annot = x$annot,
      exprs = x$exprs,
      covar = x$covar
    ),
    {
      session$setInputs(show_legend = FALSE)
      session$setInputs(`gene_selector-modus` = "by_name")
      session$flushReact()

      session$setInputs(
        `gene_selector-name_id_col` = "PrimaryID",
        `gene_selector-name_list` = "g1,g2"
      )
      session$flushReact()

      ret <- session$returned
      expect_equal(isolate(ret$genes()), c("g1", "g2"))
      expect_equal(isolate(ret$dataset()), "default")
      expect_true(methods::is(isolate(ret$heatmap()), "Heatmap"))
    }
  )
})

test_that("heatmapServer supports custom sample_id_col and validates missing columns", {
  x <- .heatmap_test_inputs(use_label = TRUE)

  testServer(
    heatmapServer,
    args = list(
      annot = x$annot,
      exprs = x$exprs,
      covar = x$covar,
      sample_id_col = "label"
    ),
    {
      session$setInputs(`gene_selector-modus` = "by_name")
      session$flushReact()

      session$setInputs(
        `gene_selector-name_id_col` = "PrimaryID",
        `gene_selector-name_list` = "g3,g4"
      )
      session$flushReact()

      ret <- session$returned
      expect_equal(isolate(ret$genes()), c("g3", "g4"))
      expect_true(methods::is(isolate(ret$heatmap()), "Heatmap"))
    }
  )

  expect_error(
    testServer(
      heatmapServer,
      args = list(
        annot = x$annot,
        exprs = x$exprs,
        covar = x$covar,
        sample_id_col = "missing_col"
      ),
      {}
    ),
    "sample_id_col"
  )
})
