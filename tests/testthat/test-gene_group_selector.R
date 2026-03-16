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
      session$setInputs(name_list = "")
      session$flushReact()
    }
  )

  expect_equal(isolate(selected_ids()), character(0))
})

test_that("geneGroupSelectorServer defaults to expression mode and top 50 by mean", {
  annot <- data.frame(
    PrimaryID = paste0("g", 1:60),
    stringsAsFactors = FALSE
  )

  x1 <- seq_len(60)
  x2 <- x1 + (61 - x1) / 10
  exprs <- rbind(x1, x2)
  exprs <- t(exprs)
  rownames(exprs) <- annot$PrimaryID
  colnames(exprs) <- c("s1", "s2")

  selected_ids <- reactiveVal(character(0))

  testServer(
    geneGroupSelectorServer,
    args = list(
      annot = annot,
      exprs = exprs,
      selected_ids = selected_ids
    ),
    {
      session$flushReact()
      expect_equal(isolate(session$returned$modus()), "by_expression")
    }
  )

  expect_equal(length(isolate(selected_ids())), 50)
  expect_equal(isolate(selected_ids()), paste0("g", 60:11))
})

test_that("geneGroupSelectorServer supports custom mode order and defaults", {
  annot <- data.frame(
    PrimaryID = paste0("g", 1:8),
    stringsAsFactors = FALSE
  )
  exprs <- matrix(
    seq_len(16),
    nrow = 8,
    dimnames = list(annot$PrimaryID, c("s1", "s2"))
  )
  cntr <- list(
    contrast_a = data.frame(
      PrimaryID = annot$PrimaryID,
      pvalue = rep(0.01, 8),
      log2FoldChange = seq_len(8),
      padj = rep(0.05, 8),
      stringsAsFactors = FALSE
    )
  )

  selected_ids <- reactiveVal(character(0))

  testServer(
    geneGroupSelectorServer,
    args = list(
      annot = annot,
      exprs = exprs,
      cntr = cntr,
      mode_order = c("by_name", "by_expression", "by_dge"),
      defaults = list(expr_top_value = 3),
      selected_ids = selected_ids
    ),
    {
      session$flushReact()
      expect_equal(isolate(session$returned$modus()), "by_name")
      session$setInputs(modus = "by_expression")
      session$flushReact()
    }
  )

  expect_equal(isolate(selected_ids()), c("g8", "g7", "g6"))
})

test_that("geneGroupSelectorServer uses server-level DGE column params", {
  annot <- data.frame(
    PrimaryID = paste0("g", 1:4),
    stringsAsFactors = FALSE
  )
  cntr <- list(
    contrast_a = data.frame(
      PrimaryID = annot$PrimaryID,
      PVAL = c(0.001, 0.02, 0.2, 0.03),
      L2FC = c(2, 1, 0.5, 1.5),
      QVAL = c(0.005, 0.02, 0.9, 0.04),
      stringsAsFactors = FALSE
    )
  )

  selected_ids <- reactiveVal(character(0))

  testServer(
    geneGroupSelectorServer,
    args = list(
      annot = annot,
      cntr = cntr,
      dge_pval_col = "PVAL",
      dge_lfc_col = "L2FC",
      dge_fdr_col = "QVAL",
      selected_ids = selected_ids
    ),
    {
      session$setInputs(modus = "by_dge")
      session$setInputs(dge_contrasts = "contrast_a")
      session$setInputs(dge_pval_thr = 0.05)
      session$setInputs(dge_lfc_thr = 0)
      session$flushReact()
    }
  )

  expect_equal(isolate(selected_ids()), c("g1", "g2", "g4"))
})

test_that("geneGroupSelectorServer does not auto-detect DGE columns", {
  annot <- data.frame(
    PrimaryID = paste0("g", 1:4),
    stringsAsFactors = FALSE
  )
  cntr <- list(
    contrast_a = data.frame(
      PrimaryID = annot$PrimaryID,
      pvalue = c(0.001, 0.02, 0.2, 0.03),
      log2FoldChange = c(2, 1, 0.5, 1.5),
      padj = c(0.005, 0.02, 0.9, 0.04),
      stringsAsFactors = FALSE
    )
  )

  selected_ids <- reactiveVal(character(0))

  testServer(
    geneGroupSelectorServer,
    args = list(
      annot = annot,
      cntr = cntr,
      selected_ids = selected_ids
    ),
    {
      session$setInputs(modus = "by_dge")
      session$setInputs(dge_contrasts = "contrast_a")
      session$setInputs(dge_pval_thr = 0.05)
      session$setInputs(dge_lfc_thr = 0)
      session$flushReact()
    }
  )

  expect_equal(isolate(selected_ids()), character(0))
})

test_that("geneGroupSelectorServer synchronizes dataset reactiveVal with UI", {
  annot <- list(
    ds_a = data.frame(
      PrimaryID = c("g1", "g2"),
      SYMBOL = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    ds_b = data.frame(
      PrimaryID = c("h1", "h2"),
      SYMBOL = c("C", "D"),
      stringsAsFactors = FALSE
    )
  )

  dataset <- reactiveVal("ds_b")

  testServer(
    geneGroupSelectorServer,
    args = list(
      annot = annot,
      dataset = dataset
    ),
    {
      session$flushReact()
      expect_equal(isolate(session$returned$dataset()), "ds_b")

      session$setInputs(dataset = "ds_a")
      session$flushReact()
      expect_equal(isolate(session$returned$dataset()), "ds_a")
    }
  )

  expect_equal(isolate(dataset()), "ds_a")
})

test_that("geneGroupSelectorServer validates defaults strictly", {
  annot <- data.frame(
    PrimaryID = c("g1", "g2"),
    stringsAsFactors = FALSE
  )

  expect_error(
    testServer(
      geneGroupSelectorServer,
      args = list(
        annot = annot,
        defaults = list(unknown_default = 1)
      ),
      {
        session$flushReact()
      }
    ),
    "unsupported key"
  )

  expect_error(
    testServer(
      geneGroupSelectorServer,
      args = list(
        annot = annot,
        defaults = list(expr_top_mode = "invalid_mode")
      ),
      {
        session$flushReact()
      }
    ),
    "expr_top_mode"
  )

  expect_error(
    testServer(
      geneGroupSelectorServer,
      args = list(
        annot = annot,
        defaults = list(dge_pval_thr = 2)
      ),
      {
        session$flushReact()
      }
    ),
    "dge_pval_thr"
  )

  expect_error(
    testServer(
      geneGroupSelectorServer,
      args = list(
        annot = annot,
        defaults = list(dge_require_fdr = "yes")
      ),
      {
        session$flushReact()
      }
    ),
    "TRUE or FALSE"
  )
})
