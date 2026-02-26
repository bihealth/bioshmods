.plot_heatmap_test_inputs <- function() {
  set.seed(1)
  exprs <- matrix(
    rnorm(24),
    nrow = 6,
    dimnames = list(paste0("g", 1:6), paste0("s", 1:4))
  )
  covar <- data.frame(
    SampleID = colnames(exprs),
    group = c("A", "A", "B", "B"),
    sex = c("f", "m", "f", "m"),
    stringsAsFactors = FALSE
  )
  list(exprs = exprs, covar = covar)
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
})
