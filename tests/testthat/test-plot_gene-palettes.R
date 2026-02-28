test_that(".cast_palette_to_ggplot converts discrete palette entries", {
  entry <- list(
    type="categorical",
    levels=c("A", "B"),
    pal=c("A"="#112233", "B"="#445566")
  )

  sc <- bioshmods:::.cast_palette_to_ggplot(entry)

  expect_equal(sc$type, "manual")
  expect_named(sc$values, c("A", "B"))
  expect_equal(unname(sc$values), c("#112233", "#445566"))
})

test_that(".cast_palette_to_ggplot converts continuous palette entries", {
  entry <- list(
    type="continuous",
    breaks=c(0, 1),
    pal=circlize::colorRamp2(c(0, 1), c("#0000FF", "#FF0000"))
  )

  sc <- bioshmods:::.cast_palette_to_ggplot(entry)

  expect_equal(sc$type, "gradientn")
  expect_equal(sc$limits, c(0, 1))
  expect_equal(sc$values, c(0, 1))
  expect_length(sc$colours, 2)
})

test_that(".gene_browser_palette_scale resolves dataset-scoped palettes", {
  pal_all <- list(
    ds_a=list(group=list(type="categorical", levels=c("A", "B"), pal=c("A"="#111111", "B"="#222222"))),
    ds_b=list(group=list(type="categorical", levels=c("A", "B"), pal=c("A"="#333333", "B"="#444444")))
  )

  sc <- bioshmods:::.gene_browser_palette_scale(
    palettes=pal_all,
    dataset="ds_b",
    color_by="group"
  )

  expect_equal(sc$type, "manual")
  expect_equal(unname(sc$values), c("#333333", "#444444"))
})

test_that("plot_gene uses manual colorScale when provided", {
  set.seed(1)
  exprs <- matrix(rnorm(8), nrow=1, dimnames=list("gene1", paste0("s", seq_len(8))))
  covar <- data.frame(
    x=seq_len(8),
    grp=rep(c("A", "B"), each=4),
    stringsAsFactors=FALSE
  )

  scale_spec <- list(type="manual", values=c("A"="#112233", "B"="#445566"))
  g <- plot_gene(
    id="gene1",
    xCovar="x",
    exprs=exprs,
    covar=covar,
    colorBy="grp",
    colorScale=scale_spec
  )

  built <- ggplot2::ggplot_build(g)
  cols <- unique(toupper(built$data[[1]]$colour))
  expect_setequal(cols, toupper(unname(scale_spec$values)))
})
