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

  expect_named(normalized, c("cntr", "annot", "annot_full"))
  expect_equal(normalized$cntr$default$contrast_a$PrimaryID, c("g1", "g2"))
  expect_equal(normalized$annot$default$PrimaryID, c("g1", "g2"))
  expect_equal(normalized$annot_full$default$SYMBOL, c("A", "B"))
})

test_that("volcanoUI contains expected controls and plot outputs", {
  ui <- volcanoUI("vol", datasets=c("ds_a", "ds_b"))
  html <- as.character(ui)

  expect_s3_class(ui, "shiny.tag")
  expect_true(grepl("vol-dataset", html, fixed=TRUE))
  expect_true(grepl("vol-pval_thr", html, fixed=TRUE))
  expect_true(grepl("vol-lfc_thr", html, fixed=TRUE))
  expect_true(grepl("vol-save", html, fixed=TRUE))
  expect_true(grepl("vol-show_top_labels", html, fixed=TRUE))
  expect_true(grepl("vol-top_label_n", html, fixed=TRUE))
  expect_true(grepl("vol-label_col_ui", html, fixed=TRUE))
  expect_true(grepl("vol-volcanoPlot", html, fixed=TRUE))
  expect_true(grepl("vol-point_id", html, fixed=TRUE))
  expect_true(grepl("vol-sel_genes", html, fixed=TRUE))
  expect_true(grepl("vol-show_selected_ui", html, fixed=TRUE))
  expect_lt(
    regexpr("vol-show_selected_ui", html, fixed=TRUE)[1],
    regexpr("vol-sel_genes", html, fixed=TRUE)[1]
  )
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
  expect_named(normalized$annot_full, c("ds_a", "ds_b"))
  expect_named(normalized$annot$ds_a, c("PrimaryID", "SYMBOL"))
  expect_named(normalized$annot$ds_b, c("PrimaryID", "SYMBOL"))
  expect_named(normalized$annot_full$ds_a, c("PrimaryID", "SYMBOL", "EXTRA"))
  expect_named(normalized$annot_full$ds_b, c("PrimaryID", "SYMBOL", "EXTRA"))
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

test_that("volcanoServer updates external selection ids on Show button click", {
  cntr <- list(
    contrast_a=.volcano_test_contrast(c("g1", "g2", "g2"), c(1, -1, -1), c(0.01, 0.02, 0.02))
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2"),
    SYMBOL=c("A", "B"),
    stringsAsFactors=FALSE
  )
  selection <- reactiveValues(ids = character(0))

  testServer(
    volcanoServer,
    args=list(
      cntr=cntr,
      annot=annot,
      selection=selection
    ),
    {
      selected_genes(data.frame(
        Dataset=c("default", "default", "default"),
        PrimaryID=c("g2", "g1", "g2"),
        stringsAsFactors=FALSE
      ))
      session$flushReact()

      expect_true(any(grepl("Show", as.character(output$show_selected_ui), fixed=TRUE)))
      expect_true(any(grepl("Export to file", as.character(output$show_selected_ui), fixed=TRUE)))
      expect_true(any(grepl("export_selected", as.character(output$show_selected_ui), fixed=TRUE)))

      session$setInputs(show_selected=1)
      session$flushReact()
    }
  )

  expect_equal(isolate(selection$ids), c("g2", "g1"))
})

test_that("volcanoServer uses configurable show button label", {
  cntr <- list(
    contrast_a=.volcano_test_contrast(c("g1", "g2"), c(1, -1), c(0.01, 0.02))
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2"),
    SYMBOL=c("A", "B"),
    stringsAsFactors=FALSE
  )
  selection <- reactiveValues(ids = character(0))

  testServer(
    volcanoServer,
    args=list(
      cntr=cntr,
      annot=annot,
      selection=selection,
      ui_config=list(show_button_label="Send to heatmap")
    ),
    {
      selected_genes(data.frame(
        Dataset=c("default", "default"),
        PrimaryID=c("g2", "g1"),
        stringsAsFactors=FALSE
      ))
      session$flushReact()

      expect_true(any(grepl("Send to heatmap", as.character(output$show_selected_ui), fixed=TRUE)))
      expect_true(any(grepl("Export to file", as.character(output$show_selected_ui), fixed=TRUE)))
    }
  )
})

test_that("volcanoServer validates selection", {
  cntr <- list(
    contrast_a=.volcano_test_contrast(c("g1"), 1, 0.01)
  )
  annot <- data.frame(
    PrimaryID="g1",
    stringsAsFactors=FALSE
  )

  expect_error(
    testServer(
      volcanoServer,
      args=list(
        cntr=cntr,
        annot=annot,
        selection=character(0)
      ),
      {
      }
    ),
    "selection"
  )

  expect_error(
    testServer(
      volcanoServer,
      args=list(
        cntr=cntr,
        annot=annot,
        ui_config=list(show_button_label="")
      ),
      {
      }
    ),
    "show_button_label"
  )
})

test_that("volcanoServer shows annot_show columns on hover", {
  cntr <- list(
    contrast_a=.volcano_test_contrast(c("g1", "g2"), c(1, -1), c(0.01, 0.02))
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2"),
    SYMBOL=c("A", "B"),
    ENTREZID=c("11", "22"),
    stringsAsFactors=FALSE
  )

  testServer(
    volcanoServer,
    args=list(
      cntr=cntr,
      annot=annot
    ),
    {
      hover_genes(data.frame(Dataset="default", PrimaryID="g2", stringsAsFactors=FALSE))
      session$flushReact()

      html <- output$point_id
      expect_match(html, "Dataset")
      expect_match(html, "PrimaryID")
      expect_match(html, "SYMBOL")
      expect_match(html, "ENTREZID")
      expect_match(html, "default")
      expect_match(html, "g2")
      expect_match(html, "B")
      expect_match(html, "22")
    }
  )
})

test_that("volcanoServer adds text labels when top labels are enabled", {
  cntr <- list(
    contrast_a=.volcano_test_contrast(
      c("g1", "g2", "g3"),
      c(2, -1.5, 0.2),
      c(0.001, 0.01, 0.5)
    )
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2", "g3"),
    stringsAsFactors=FALSE
  )

  testServer(
    volcanoServer,
    args=list(
      cntr=cntr,
      annot=annot
    ),
    {
      session$setInputs(
        .clientdata_output_volcanoPlot_width=800,
        .clientdata_output_volcanoPlot_height=600
      )
      session$setInputs(
        dataset="default",
        lfc_thr=0,
        pval_thr=1,
        samescaleX=TRUE,
        samescaleY=TRUE,
        figure_size="800x600",
        font_size=12,
        show_top_labels=TRUE,
        top_label_n=2
      )
      session$flushReact()
      output$volcanoPlot
      session$flushReact()

      g <- isolate(plot_obj())
      expect_s3_class(g, "ggplot")
      expect_true(any(vapply(g$layers, function(x) inherits(x$geom, "GeomText"), logical(1))))
    }
  )
})

test_that("volcanoServer chooses top labels per contrast", {
  cntr <- list(
    contrast_a=.volcano_test_contrast(
      c("g1", "g2"),
      c(2, 0.5),
      c(0.001, 0.2)
    ),
    contrast_b=.volcano_test_contrast(
      c("g3", "g4"),
      c(-1.8, 0.2),
      c(0.002, 0.4)
    )
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2", "g3", "g4"),
    stringsAsFactors=FALSE
  )

  testServer(
    volcanoServer,
    args=list(
      cntr=cntr,
      annot=annot
    ),
    {
      session$setInputs(
        .clientdata_output_volcanoPlot_width=800,
        .clientdata_output_volcanoPlot_height=600
      )
      session$setInputs(
        dataset="default",
        lfc_thr=0,
        pval_thr=1,
        samescaleX=TRUE,
        samescaleY=TRUE,
        figure_size="800x600",
        font_size=12,
        show_top_labels=TRUE,
        top_label_n=1
      )
      session$flushReact()
      output$volcanoPlot
      session$flushReact()

      g <- isolate(plot_obj())
      text_layer <- g$layers[[which(vapply(g$layers, function(x) inherits(x$geom, "GeomText"), logical(1)))[1]]]
      expect_equal(nrow(text_layer$data), 2)
      expect_setequal(text_layer$data$PrimaryID, c("g1", "g3"))
    }
  )
})

test_that("volcanoServer shows label column selector when annotation columns are available", {
  cntr <- list(
    contrast_a=.volcano_test_contrast(c("g1", "g2"), c(1, -1), c(0.01, 0.02))
  )
  annot <- data.frame(
    PrimaryID=c("g1", "g2"),
    SYMBOL=c("A", "B"),
    ENTREZID=c("11", "22"),
    stringsAsFactors=FALSE
  )

  testServer(
    volcanoServer,
    args=list(
      cntr=cntr,
      annot=annot
    ),
    {
      session$flushReact()
      html <- paste(as.character(output$label_col_ui), collapse="\n")
      expect_match(html, "Label column")
      expect_match(html, "Primary ID")
      expect_match(html, "SYMBOL")
      expect_match(html, "ENTREZID")
    }
  )
})
