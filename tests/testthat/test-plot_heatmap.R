.plot_heatmap_test_inputs <- function() {
  set.seed(1)
  exprs <- matrix(
    rnorm(24),
    nrow = 6,
    dimnames = list(paste0("g", 1:6), paste0("s", 1:4))
  )
  annot <- data.frame(
    PrimaryID = rownames(exprs),
    SYMBOL = LETTERS[1:6],
    stringsAsFactors = FALSE
  )
  covar <- data.frame(
    SampleID = colnames(exprs),
    group = c("A", "A", "B", "B"),
    sex = c("f", "m", "f", "m"),
    stringsAsFactors = FALSE
  )
  list(exprs = exprs, annot = annot, covar = covar)
}

test_that("plot_heatmap returns Heatmap without annotations by default", {
  x <- .plot_heatmap_test_inputs()

  hm <- plot_heatmap(
    exprs = x$exprs,
    genes = c("g1", "g3", "g5"),
    covar = x$covar
  )

  expect_true(methods::is(hm, "Heatmap"))
  expect_null(hm@top_annotation)
})

test_that("plot_heatmap maps row labels from annot when annot_row_col is available", {
  x <- .plot_heatmap_test_inputs()

  hm <- plot_heatmap(
    exprs = x$exprs,
    genes = c("g1", "g3", "g5"),
    annot = x$annot,
    primary_id_col = "PrimaryID",
    annot_row_col = "SYMBOL"
  )

  expect_equal(hm@row_names_param$labels, c("A", "C", "E"))
})

test_that("plot_heatmap keeps gene IDs when annot_row_col is missing in annot", {
  x <- .plot_heatmap_test_inputs()

  hm <- plot_heatmap(
    exprs = x$exprs,
    genes = c("g2", "g4"),
    annot = x$annot,
    annot_row_col = "NOT_PRESENT"
  )

  expect_equal(hm@row_names_param$labels, c("g2", "g4"))
})

test_that("plot_heatmap uses sel_annot to select annotation bars", {
  x <- .plot_heatmap_test_inputs()

  hm <- plot_heatmap(
    exprs = x$exprs,
    genes = c("g1", "g3", "g5"),
    covar = x$covar,
    sel_annot = "group",
    legend = FALSE
  )

  expect_true(methods::is(hm, "Heatmap"))
  expect_true(methods::is(hm@top_annotation, "HeatmapAnnotation"))
  expect_equal(names(hm@top_annotation@anno_list), "group")
})

test_that("plot_heatmap validates sample_id_col and sel_annot", {
  x <- .plot_heatmap_test_inputs()

  expect_error(
    plot_heatmap(
      exprs = x$exprs,
      genes = "g1",
      covar = x$covar,
      sample_id_col = "missing_id_col"
    ),
    "sample_id_col"
  )

  expect_error(
    plot_heatmap(
      exprs = x$exprs,
      genes = "g1",
      covar = x$covar,
      sel_annot = "missing_annotation"
    ),
    "sel_annot"
  )

  covar_bad <- x$covar
  covar_bad$SampleID <- paste0("x", seq_len(nrow(covar_bad)))
  expect_error(
    plot_heatmap(
      exprs = x$exprs,
      genes = "g1",
      covar = covar_bad,
      sel_annot = "group"
    ),
    "No matching samples"
  )

  expect_error(
    plot_heatmap(
      exprs = x$exprs,
      genes = "g1",
      annot = x$annot[, "SYMBOL", drop = FALSE],
      annot_row_col = "SYMBOL"
    ),
    "primary_id_col"
  )
})
