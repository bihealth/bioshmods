library(shiny)

.volcano_test_contrast <- function(ids, lfc, pval, extra=NULL) {
  x <- data.frame(
    PrimaryID=ids,
    log2FoldChange=lfc,
    padj=pval,
    stringsAsFactors=FALSE
  )
  if(!is.null(extra) && is.list(extra) && length(extra) > 0L) {
    for(nm in names(extra)) {
      x[[nm]] <- extra[[nm]]
    }
  }
  x
}

test_that(".normalize_volcano_inputs keeps explicit default primary_id", {
  cntr <- list(
    contrast_a=data.frame(
      PrimaryID=c("g1", "g2"),
      log2FoldChange=c(1, -1),
      padj=c(0.01, 0.02),
      stringsAsFactors=FALSE
    )
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2"),
    SYMBOL=c("A", "B"),
    stringsAsFactors=FALSE
  )

  normalized <- bioshmods:::.normalize_volcano_inputs(
    cntr=cntr,
    annot=annot,
    primary_id="PrimaryID",
    lfc_col="log2FoldChange",
    pval_col="padj",
    annot_show="SYMBOL"
  )

  expect_named(normalized, c("cntr", "annot"))
  expect_equal(normalized$cntr$default$contrast_a$PrimaryID, c("g1", "g2"))
  expect_equal(normalized$annot$default$PrimaryID, c("g1", "g2"))
})

test_that("volcanoUI contains expected controls and plot outputs", {
  ui <- volcanoUI("vol", datasets=c("ds_a", "ds_b"))
  html <- as.character(ui)

  expect_s3_class(ui, "shiny.tag")
  expect_true(grepl("vol-dataset", html, fixed=TRUE))
  expect_true(grepl("vol-pval_thr", html, fixed=TRUE))
  expect_true(grepl("vol-lfc_thr", html, fixed=TRUE))
  expect_true(grepl("vol-save", html, fixed=TRUE))
  expect_true(grepl("vol-volcanoPlot", html, fixed=TRUE))
  expect_true(grepl("vol-point_id", html, fixed=TRUE))
  expect_true(grepl("vol-sel_genes", html, fixed=TRUE))
})

test_that(".normalize_volcano_inputs errors when primary_id is NULL", {
  cntr <- list(
    contrast_a=data.frame(
      PrimaryID=c("g1", "g2"),
      log2FoldChange=c(1, -1),
      padj=c(0.01, 0.02),
      stringsAsFactors=FALSE
    )
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2"),
    stringsAsFactors=FALSE
  )

  expect_error(
    bioshmods:::.normalize_volcano_inputs(
      cntr=cntr,
      annot=annot,
      primary_id=NULL,
      lfc_col="log2FoldChange",
      pval_col="padj",
      annot_show=character(0)
    ),
    "non-empty character column name"
  )
})

test_that(".normalize_volcano_inputs supports multi-dataset input and trims annot_show columns", {
  cntr <- list(
    ds_a=list(
      c1=.volcano_test_contrast(c("g1", "g2"), c(1, -1), c(0.01, 0.02))
    ),
    ds_b=list(
      c2=.volcano_test_contrast(c("h1", "h2"), c(1.5, -2), c(0.03, 0.04))
    )
  )
  annot <- list(
    ds_a=data.frame(
      PrimaryID=c("g1", "g2"),
      SYMBOL=c("A", "B"),
      EXTRA=c("x", "y"),
      stringsAsFactors=FALSE
    ),
    ds_b=data.frame(
      PrimaryID=c("h1", "h2"),
      SYMBOL=c("C", "D"),
      EXTRA=c("u", "v"),
      stringsAsFactors=FALSE
    )
  )

  normalized <- bioshmods:::.normalize_volcano_inputs(
    cntr=cntr,
    annot=annot,
    primary_id="PrimaryID",
    lfc_col="log2FoldChange",
    pval_col="padj",
    annot_show="SYMBOL"
  )

  expect_named(normalized$cntr, c("ds_a", "ds_b"))
  expect_named(normalized$annot, c("ds_a", "ds_b"))
  expect_named(normalized$annot$ds_a, c("PrimaryID", "SYMBOL"))
  expect_named(normalized$annot$ds_b, c("PrimaryID", "SYMBOL"))
})

test_that(".normalize_volcano_inputs errors when explicit primary_id is missing in contrasts", {
  cntr <- list(
    contrast_a=data.frame(
      log2FoldChange=c(1, -1),
      padj=c(0.01, 0.02),
      stringsAsFactors=FALSE
    )
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2"),
    stringsAsFactors=FALSE
  )

  expect_error(
    bioshmods:::.normalize_volcano_inputs(
      cntr=cntr,
      annot=annot,
      primary_id="PrimaryID",
      lfc_col="log2FoldChange",
      pval_col="padj",
      annot_show=character(0)
    ),
    "missing 'PrimaryID'"
  )
})

test_that(".normalize_volcano_inputs errors when explicit primary_id is missing in annotation", {
  cntr <- list(
    contrast_a=data.frame(
      PrimaryID=c("g1", "g2"),
      log2FoldChange=c(1, -1),
      padj=c(0.01, 0.02),
      stringsAsFactors=FALSE
    )
  )
  annot <- data.frame(
    SYMBOL=c("A", "B"),
    stringsAsFactors=FALSE
  )

  expect_error(
    bioshmods:::.normalize_volcano_inputs(
      cntr=cntr,
      annot=annot,
      primary_id="PrimaryID",
      lfc_col="log2FoldChange",
      pval_col="padj",
      annot_show="SYMBOL"
    ),
    "must contain 'PrimaryID'"
  )
})

test_that(".volcano_process_data_one_ds builds expected long table", {
  cntr_ds <- list(
    c1=.volcano_test_contrast(c("g1", "g2"), c(1, -1), c(0.01, 0.02), extra=list(foo=1:2)),
    c2=.volcano_test_contrast(c("g3"), c(2), c(0.03), extra=list(foo=3))
  )

  df <- bioshmods:::.volcano_process_data_one_ds(
    ds_id="ds_a",
    cntr=cntr_ds,
    annot=data.frame(),
    primary_id="PrimaryID",
    lfc_col="log2FoldChange",
    pval_col="padj"
  )

  expect_named(df, c("PrimaryID", "log2FoldChange", "padj", "Dataset", "Contrast"))
  expect_equal(nrow(df), 3)
  expect_equal(unique(df$Dataset), "ds_a")
  expect_setequal(unique(df$Contrast), c("c1", "c2"))
  expect_equal(df$PrimaryID, c("g1", "g2", "g3"))
})

test_that(".volcano_process_data combines datasets", {
  cntr <- list(
    ds_a=list(c1=.volcano_test_contrast(c("g1"), c(1), c(0.01))),
    ds_b=list(c2=.volcano_test_contrast(c("h1"), c(-1), c(0.02)))
  )
  annot <- list(
    ds_a=data.frame(PrimaryID="g1", stringsAsFactors=FALSE),
    ds_b=data.frame(PrimaryID="h1", stringsAsFactors=FALSE)
  )

  df <- bioshmods:::.volcano_process_data(
    cntr=cntr,
    annot=annot,
    primary_id="PrimaryID",
    lfc_col="log2FoldChange",
    pval_col="padj"
  )

  expect_equal(nrow(df), 2)
  expect_setequal(unique(df$Dataset), c("ds_a", "ds_b"))
  expect_setequal(unique(df$Contrast), c("c1", "c2"))
})

test_that("volcanoServer updates external gene_id on gene button click", {
  cntr <- list(
    contrast_a=.volcano_test_contrast(c("g1", "g2"), c(1, -1), c(0.01, 0.02))
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2"),
    SYMBOL=c("A", "B"),
    stringsAsFactors=FALSE
  )
  gene_id <- reactiveValues(ds=NULL, id=NULL)

  testServer(
    volcanoServer,
    args=list(
      cntr=cntr,
      annot=annot,
      gene_id=gene_id
    ),
    {
      session$setInputs(genebutton="vol-gene_id~default~g2")
      session$flushReact()
    }
  )

  expect_equal(isolate(gene_id$ds), "default")
  expect_equal(isolate(gene_id$id), "g2")
})
