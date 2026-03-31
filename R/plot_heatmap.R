# Build a sample-ordered covariate table for heatmap annotations.
# Requires an explicit sample_id_col in covar.
.plot_heatmap_log <- function(...) {
  .bioshmods_log(..., .prefix="plot_heatmap")
}

.plot_heatmap_prepare_covar <- function(covar, sample_ids, sample_id_col="SampleID", sel_annot=NULL) {
  .plot_heatmap_log("prepare_covar called; sample_id_col='", sample_id_col,
                    "', n_sample_ids=", as.character(length(sample_ids)),
                    ", sel_annot={", paste(sel_annot %||% character(0), collapse=","), "}.")
  if(is.null(covar)) {
    .plot_heatmap_log("covar is NULL; no top annotations will be generated.")
    return(NULL)
  }

  covar <- as.data.frame(covar, stringsAsFactors=FALSE)
  .plot_heatmap_log("covar columns={", paste(colnames(covar), collapse=","), "}; nrow=",
                    as.character(nrow(covar)), ".")
  if(!sample_id_col %in% colnames(covar)) {
    stop(sprintf("`sample_id_col` ('%s') not found in `covar`.", sample_id_col))
  }

  sel_annot <- unique(as.character(sel_annot))
  sel_annot <- sel_annot[!is.na(sel_annot) & sel_annot != ""]
  if(length(sel_annot) < 1L) {
    .plot_heatmap_log("no selected annotation columns after cleanup; returning NULL.")
    return(NULL)
  }

  missing_annot <- setdiff(sel_annot, colnames(covar))
  if(length(missing_annot) > 0L) {
    .plot_heatmap_log("missing annotation columns in covar: ", paste(missing_annot, collapse=","))
    stop(sprintf(
      "`sel_annot` column(s) not found in `covar`: %s",
      paste(missing_annot, collapse=", ")
    ))
  }

  covar_ids <- as.character(covar[[sample_id_col]])
  valid <- !is.na(covar_ids) & covar_ids != ""
  covar <- covar[valid, , drop=FALSE]
  covar_ids <- covar_ids[valid]

  if(length(intersect(sample_ids, covar_ids)) < 1L) {
    stop("No matching samples between `exprs` column names and `covar[[sample_id_col]]`.")
  }
  .plot_heatmap_log("matched samples between expression and covar: ",
                    as.character(length(intersect(sample_ids, covar_ids))), ".")

  # Keep first occurrence for duplicated sample IDs.
  dedup <- !duplicated(covar_ids)
  covar <- covar[dedup, , drop=FALSE]
  covar_ids <- covar_ids[dedup]

  ann_cols <- setdiff(sel_annot, sample_id_col)
  if(length(ann_cols) < 1L) {
    .plot_heatmap_log("annotation columns only contained sample_id_col; returning NULL.")
    return(NULL)
  }

  idx <- match(sample_ids, covar_ids)
  ann_df <- covar[idx, ann_cols, drop=FALSE]
  rownames(ann_df) <- sample_ids
  .plot_heatmap_log("prepared annotation data frame with nrow=", as.character(nrow(ann_df)),
                    ", ncol=", as.character(ncol(ann_df)), ", columns={",
                    paste(colnames(ann_df), collapse=","), "}.")
  ann_df
}

# Reorder heatmap columns by selected covariates in the given selector order.
# Keeps ties stable and disables column clustering when a covariate sort is active.
.plot_heatmap_order_columns_by_covar <- function(exprs, ann_df) {
  if(is.null(ann_df) || ncol(ann_df) < 1L || ncol(exprs) < 2L) {
    return(list(exprs=exprs, ann_df=ann_df, sorted=FALSE))
  }

  order_args <- lapply(ann_df, function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  })
  order_args <- c(order_args, list(seq_len(nrow(ann_df)), na.last=TRUE, method="radix"))
  ord <- do.call(order, order_args)

  .plot_heatmap_log("sorting columns by covariates; ordered sample IDs={",
                    paste(rownames(ann_df)[ord], collapse=","), "}.")
  list(
    exprs=exprs[, ord, drop=FALSE],
    ann_df=ann_df[ord, , drop=FALSE],
    sorted=TRUE
  )
}

# Scale expression by gene (row-wise z-score), replacing non-finite values with 0.
.plot_heatmap_scale_rows <- function(x) {
  x <- t(scale(t(x)))
  x[!is.finite(x)] <- 0
  x
}

# Resolve row labels for plotted genes.
# When annot_row_col is available, labels are merged from annot by primary_id_col.
.plot_heatmap_row_labels <- function(gene_ids, annot=NULL, primary_id_col="PrimaryID", annot_row_col=NULL) {
  gene_ids <- as.character(gene_ids)
  row_labels <- gene_ids
  .plot_heatmap_log("row label resolver called; n_genes=", as.character(length(gene_ids)),
                    ", primary_id_col='", primary_id_col,
                    "', annot_row_col='", as.character(annot_row_col)[1], "'.")

  annot_row_col <- as.character(annot_row_col)[1]
  if(is.null(annot) || is.na(annot_row_col) || !nzchar(annot_row_col)) {
    .plot_heatmap_log("annotation or annot_row_col missing; using primary IDs as row labels.")
    return(row_labels)
  }

  annot <- as.data.frame(annot, stringsAsFactors=FALSE)
  if(!annot_row_col %in% colnames(annot)) {
    .plot_heatmap_log("annot_row_col not found in annot; using primary IDs as row labels.")
    return(row_labels)
  }
  if(!primary_id_col %in% colnames(annot)) {
    stop(sprintf(
      "`primary_id_col` ('%s') not found in `annot`, cannot map `annot_row_col`.",
      primary_id_col
    ))
  }

  hm_df <- data.frame(.row_idx=seq_along(gene_ids), stringsAsFactors=FALSE)
  hm_df[[primary_id_col]] <- gene_ids

  annot_df <- annot[, unique(c(primary_id_col, annot_row_col)), drop=FALSE]
  annot_df[[primary_id_col]] <- as.character(annot_df[[primary_id_col]])
  annot_df <- annot_df[!is.na(annot_df[[primary_id_col]]) & annot_df[[primary_id_col]] != "", , drop=FALSE]
  annot_df <- annot_df[!duplicated(annot_df[[primary_id_col]]), , drop=FALSE]

  hm_df <- merge(
    hm_df,
    annot_df,
    by=primary_id_col,
    all.x=TRUE,
    sort=FALSE
  )
  hm_df <- hm_df[order(hm_df$.row_idx), , drop=FALSE]

  lbl <- as.character(hm_df[[annot_row_col]])
  fallback <- as.character(hm_df[[primary_id_col]])
  lbl[is.na(lbl) | lbl == ""] <- fallback[is.na(lbl) | lbl == ""]
  .plot_heatmap_log("row labels resolved from annotation; non-empty labels=",
                    as.character(sum(!is.na(lbl) & lbl != "")), ".")
  lbl
}

#' Plot Expression Heatmap
#'
#' Generate a heatmap object from an expression matrix, a selected gene set,
#' and optional sample covariates.
#'
#' @param exprs Numeric matrix/data frame with genes in rows and samples in columns.
#'   Row names must contain gene identifiers.
#' @param genes Character vector with selected gene identifiers.
#' @param covar Optional sample covariate data frame. Sample identifiers must be in
#'   the column designated by `sample_id_col`.
#' @param sample_id_col Name of the sample ID column in `covar`.
#' @param annot Optional annotation data frame for genes.
#' @param primary_id_col Column in `annot` containing gene identifiers matching
#'   `exprs` row names.
#' @param annot_row_col Optional column in `annot` used as row labels in the heatmap.
#'   If `NULL` or missing in `annot`, row names from `exprs` are shown.
#' @param sel_annot Character vector of column names from `covar` to display as
#'   top annotations. If `NULL` (default), no annotation bars are shown.
#' @param sort_by_covar Logical; if `TRUE` and covariates are selected in
#'   `sel_annot`, samples are ordered by those covariates and column clustering
#'   is disabled.
#' @param legend Logical; whether heatmap and annotation legends should be shown.
#' @param palettes Optional named list of palette definitions. When `col` is
#'   `NULL` and `palettes$values$pal` exists, that value is used for heatmap
#'   colors.
#' @param col Optional heatmap color mapping passed to
#'   `ComplexHeatmap::Heatmap(col=...)`. This overrides `palettes`.
#'
#' @return A `ComplexHeatmap::Heatmap` object.
#'
#' @examples
#' set.seed(1)
#' exprs <- matrix(
#'   rnorm(20),
#'   nrow=5,
#'   dimnames=list(paste0("g", 1:5), paste0("s", 1:4))
#' )
#' covar <- data.frame(
#'   SampleID=colnames(exprs),
#'   group=c("A", "A", "B", "B"),
#'   stringsAsFactors=FALSE
#' )
#' annot <- data.frame(
#'   PrimaryID=rownames(exprs),
#'   SYMBOL=paste0("Gene_", seq_len(nrow(exprs))),
#'   stringsAsFactors=FALSE
#' )
#' hm <- plot_heatmap(
#'   exprs,
#'   genes=c("g1", "g3", "g5"),
#'   covar=covar,
#'   annot=annot,
#'   annot_row_col="SYMBOL",
#'   sel_annot="group"
#' )
#' hm
#'
#' data(C19)
#' hm_c19 <- plot_heatmap(
#'   exprs=C19$expression,
#'   genes=utils::head(C19$annotation$PrimaryID, 30),
#'   covar=C19$covariates,
#'   sample_id_col="label",
#'   annot=C19$annotation,
#'   primary_id_col="PrimaryID",
#'   annot_row_col="SYMBOL",
#'   sel_annot=c("group", "sex"),
#'   col=circlize::colorRamp2(c(-2, 0, 2), c("#2c7bb6", "#ffffbf", "#d7191c")),
#'   legend=FALSE
#' )
#' hm_c19
#'
#' @export
plot_heatmap <- function(exprs, genes, covar=NULL, sample_id_col="SampleID",
                         annot=NULL, primary_id_col="PrimaryID", annot_row_col=NULL,
                         sel_annot=NULL, sort_by_covar=FALSE,
                         legend=TRUE, palettes = NULL, col = NULL) {
  .plot_heatmap_log("plot_heatmap called; exprs dim=", paste(dim(exprs), collapse="x"),
                    ", genes requested=", as.character(length(genes)),
                    ", sample_id_col='", sample_id_col, "'.")
  exprs <- as.matrix(exprs)
  sample_id_col <- as.character(sample_id_col)[1]
  primary_id_col <- as.character(primary_id_col)[1]

  if(!is.numeric(exprs)) {
    stop("`exprs` must be numeric.")
  }
  if(is.null(rownames(exprs)) || any(rownames(exprs) == "")) {
    stop("`exprs` must have non-empty row names with gene identifiers.")
  }
  if(is.null(colnames(exprs)) || any(colnames(exprs) == "")) {
    stop("`exprs` must have non-empty column names with sample identifiers.")
  }
  if(is.na(sample_id_col) || !nzchar(sample_id_col)) {
    stop("`sample_id_col` must be a non-empty column name.")
  }
  if(is.na(primary_id_col) || !nzchar(primary_id_col)) {
    stop("`primary_id_col` must be a non-empty column name.")
  }

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  if(length(genes) < 1L) {
    stop("`genes` is empty.")
  }
  .plot_heatmap_log("genes after cleanup=", as.character(length(genes)), ".")

  genes <- intersect(genes, rownames(exprs))
  if(length(genes) < 1L) {
    stop("None of the selected genes were found in `exprs` row names.")
  }
  .plot_heatmap_log("genes found in expression matrix=", as.character(length(genes)), ".")

  exprs <- exprs[genes, , drop=FALSE]
  exprs <- .plot_heatmap_scale_rows(exprs)
  row_labels <- .plot_heatmap_row_labels(
    gene_ids=rownames(exprs),
    annot=annot,
    primary_id_col=primary_id_col,
    annot_row_col=annot_row_col
  )

  ann_df <- .plot_heatmap_prepare_covar(
    covar,
    colnames(exprs),
    sample_id_col=sample_id_col,
    sel_annot=sel_annot
  )
  sort_res <- .plot_heatmap_order_columns_by_covar(exprs, ann_df)
  if(isTRUE(sort_by_covar) && isTRUE(sort_res$sorted)) {
    exprs <- sort_res$exprs
    ann_df <- sort_res$ann_df
  }
  top_anno <- NULL

  cluster_cols <- !(isTRUE(sort_by_covar) && !is.null(ann_df) && ncol(ann_df) > 0L)

  if(!is.null(ann_df) && ncol(ann_df) > 0L) {
    pal_anno <- NULL
    if(!is.null(palettes)) {
      pal_anno <- lapply(palettes, \(x) x$pal)
    }

    .plot_heatmap_log("building top annotation with columns={",
                      paste(colnames(ann_df), collapse=","), "} and palette keys={",
                      paste(names(pal_anno %||% list()), collapse=","), "}.")
    top_anno <- ComplexHeatmap::HeatmapAnnotation(df=ann_df, show_legend=isTRUE(legend), col=pal_anno)
  } else {
    .plot_heatmap_log("no top annotation columns available; drawing heatmap without annotation bars.")
  }

  .plot_heatmap_log("creating ComplexHeatmap object; legend=", as.character(isTRUE(legend)), ".")
  ComplexHeatmap::Heatmap(
    exprs,
    name="Expression",
    row_labels=row_labels,
    top_annotation=top_anno,
    cluster_rows=TRUE,
    cluster_columns=cluster_cols,
    show_heatmap_legend=isTRUE(legend),
    show_column_names=FALSE,
    col=col,
    show_row_names=TRUE
  )
}
