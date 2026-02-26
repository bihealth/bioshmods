# Build a sample-ordered covariate table for heatmap annotations.
# Requires an explicit sample_id_col in covar.
.plot_heatmap_prepare_covar <- function(covar, sample_ids, sample_id_col="SampleID", sel_annot=NULL) {
  if(is.null(covar)) {
    return(NULL)
  }

  covar <- as.data.frame(covar, stringsAsFactors=FALSE)
  if(!sample_id_col %in% colnames(covar)) {
    stop(sprintf("`sample_id_col` ('%s') not found in `covar`.", sample_id_col))
  }

  sel_annot <- unique(as.character(sel_annot))
  sel_annot <- sel_annot[!is.na(sel_annot) & sel_annot != ""]
  if(length(sel_annot) < 1L) {
    return(NULL)
  }

  missing_annot <- setdiff(sel_annot, colnames(covar))
  if(length(missing_annot) > 0L) {
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

  # Keep first occurrence for duplicated sample IDs.
  dedup <- !duplicated(covar_ids)
  covar <- covar[dedup, , drop=FALSE]
  covar_ids <- covar_ids[dedup]

  ann_cols <- setdiff(sel_annot, sample_id_col)
  if(length(ann_cols) < 1L) {
    return(NULL)
  }

  idx <- match(sample_ids, covar_ids)
  ann_df <- covar[idx, ann_cols, drop=FALSE]
  rownames(ann_df) <- sample_ids
  ann_df
}

# Scale expression by gene (row-wise z-score), replacing non-finite values with 0.
.plot_heatmap_scale_rows <- function(x) {
  x <- t(scale(t(x)))
  x[!is.finite(x)] <- 0
  x
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
#' @param sel_annot Character vector of column names from `covar` to display as
#'   top annotations. If `NULL` (default), no annotation bars are shown.
#' @param legend Logical; whether heatmap and annotation legends should be shown.
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
#' hm <- plot_heatmap(exprs, genes=c("g1", "g3", "g5"), covar=covar, sel_annot="group")
#' hm
#'
#' data(C19)
#' hm_c19 <- plot_heatmap(
#'   exprs=C19$expression,
#'   genes=utils::head(C19$annotation$PrimaryID, 30),
#'   covar=C19$covariates,
#'   sample_id_col="label",
#'   sel_annot=c("group", "sex"),
#'   legend=FALSE
#' )
#' hm_c19
#'
#' @export
plot_heatmap <- function(exprs, genes, covar=NULL, sample_id_col="SampleID", sel_annot=NULL, legend=TRUE) {
  exprs <- as.matrix(exprs)
  sample_id_col <- as.character(sample_id_col)[1]

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

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  if(length(genes) < 1L) {
    stop("`genes` is empty.")
  }

  genes <- intersect(genes, rownames(exprs))
  if(length(genes) < 1L) {
    stop("None of the selected genes were found in `exprs` row names.")
  }

  exprs <- exprs[genes, , drop=FALSE]
  exprs <- .plot_heatmap_scale_rows(exprs)

  ann_df <- .plot_heatmap_prepare_covar(
    covar,
    colnames(exprs),
    sample_id_col=sample_id_col,
    sel_annot=sel_annot
  )
  top_anno <- NULL
  if(!is.null(ann_df) && ncol(ann_df) > 0L) {
    top_anno <- ComplexHeatmap::HeatmapAnnotation(df=ann_df, show_legend=isTRUE(legend))
  }

  ComplexHeatmap::Heatmap(
    exprs,
    name="Expression",
    top_annotation=top_anno,
    cluster_rows=TRUE,
    cluster_columns=TRUE,
    show_heatmap_legend=isTRUE(legend),
    show_column_names=FALSE,
    show_row_names=TRUE
  )
}
