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
  expect_true(grepl("col-sm-3", html, fixed = TRUE))
  expect_true(grepl("col-sm-9", html, fixed = TRUE))
  expect_true(grepl("hm-gene_selector-modus", html, fixed = TRUE))
  expect_true(grepl("hm-selected_genes_n", html, fixed = TRUE))
  expect_true(grepl("hm-selected_genes_warning", html, fixed = TRUE))
  expect_true(grepl("hm-sel_annot", html, fixed = TRUE))
  expect_true(grepl("hm-annot_row_col", html, fixed = TRUE))
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
      session$setInputs(sel_annot = "group")
      session$setInputs(annot_row_col = "SYMBOL")
      session$setInputs(`gene_selector-modus` = "by_name")
      session$flushReact()

      session$setInputs(
        `gene_selector-name_list` = "g1,g2"
      )
      session$flushReact()

      ret <- session$returned
      expect_equal(isolate(ret$genes()), c("g1", "g2"))
      expect_equal(isolate(ret$dataset()), "default")
      expect_true(methods::is(isolate(ret$heatmap()), "Heatmap"))
      expect_false(is.null(isolate(ret$heatmap())@top_annotation))
      expect_equal(isolate(ret$heatmap())@row_names_param$labels, c("A", "B"))
      expect_equal(output$selected_genes_n, "Selected genes: 2")
      expect_equal(output$selected_genes_warning, "")
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

  expect_error(
    testServer(
      heatmapServer,
      args = list(
        annot = x$annot,
        exprs = x$exprs,
        covar = x$covar,
        primary_id = "missing_primary_id_col"
      ),
      {}
    ),
    "primary_id"
  )
})

test_that("heatmapServer forwards DGE column params to gene selector module", {
  x <- .heatmap_test_inputs()
  cntr <- list(
    contrast_a = data.frame(
      PrimaryID = x$annot$PrimaryID,
      PVAL = c(0.001, 0.02, 0.8, 0.03, 0.6, 0.9),
      L2FC = c(2, 1.3, 0.2, 1.1, 0.1, 0.05),
      QVAL = c(0.005, 0.02, 0.9, 0.04, 0.8, 0.95),
      stringsAsFactors = FALSE
    )
  )

  testServer(
    heatmapServer,
    args = list(
      annot = x$annot,
      exprs = x$exprs,
      cntr = cntr,
      covar = x$covar,
      dge_pval_col = "PVAL",
      dge_lfc_col = "L2FC",
      dge_fdr_col = "QVAL"
    ),
    {
      session$setInputs(`gene_selector-modus` = "by_dge")
      session$setInputs(`gene_selector-dge_contrasts` = "contrast_a")
      session$setInputs(`gene_selector-dge_pval_thr` = 0.05)
      session$setInputs(`gene_selector-dge_lfc_thr` = 1)
      session$flushReact()

      ret <- session$returned
      expect_equal(isolate(ret$genes()), c("g1", "g2", "g4"))
      expect_true(methods::is(isolate(ret$heatmap()), "Heatmap"))
    }
  )
})

test_that("heatmapServer limits shown genes and displays warning", {
  n_genes <- 200
  annot <- data.frame(
    PrimaryID = paste0("g", seq_len(n_genes)),
    SYMBOL = paste0("S", seq_len(n_genes)),
    stringsAsFactors = FALSE
  )
  exprs <- matrix(
    rnorm(n_genes * 4),
    nrow = n_genes,
    dimnames = list(annot$PrimaryID, paste0("s", 1:4))
  )
  covar <- data.frame(
    SampleID = colnames(exprs),
    group = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  all_ids <- paste(annot$PrimaryID, collapse = ",")

  testServer(
    heatmapServer,
    args = list(
      annot = annot,
      exprs = exprs,
      covar = covar,
      max_genes = 150
    ),
    {
      session$setInputs(`gene_selector-modus` = "by_name")
      session$setInputs(`gene_selector-name_list` = all_ids)
      session$flushReact()

      ret <- session$returned
      expect_equal(length(isolate(ret$genes())), 200)
      expect_equal(nrow(isolate(ret$heatmap())@matrix), 150)
      expect_match(output$selected_genes_warning, "first 150")
    }
  )
})
