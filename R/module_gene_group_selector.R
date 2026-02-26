# Convert free text into a filesystem-safe filename fragment.
# Replaces disallowed characters and applies a fallback name when empty.
.gene_group_sanitize_filename <- function(x, default="selected_genes") {
  x <- trimws(as.character(x)[1])
  if(is.na(x) || x == "") {
    x <- default
  }
  x <- gsub("[^0-9A-Za-z_.-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if(x == "") {
    x <- default
  }
  x
}

# Parse comma/space separated gene identifiers from text input.
# Returns unique non-empty tokens in their original order.
.gene_group_parse_gene_input <- function(x) {
  x <- x %||% ""
  x <- paste(x, collapse=" ")
  tokens <- unlist(strsplit(x, "[,[:space:]]+"))
  tokens <- trimws(tokens)
  tokens <- tokens[!is.na(tokens) & tokens != ""]
  unique(tokens)
}

# Read a gene list file uploaded through shiny::fileInput.
# Returns parsed identifiers, or character(0) when input/file is missing.
.gene_group_read_gene_file <- function(file_input) {
  if(is.null(file_input) || is.null(file_input$datapath) || !nzchar(file_input$datapath)) {
    return(character(0))
  }
  if(!file.exists(file_input$datapath)) {
    return(character(0))
  }
  txt <- readLines(file_input$datapath, warn=FALSE)
  .gene_group_parse_gene_input(txt)
}

# Convert absolute or percentage "top" requests into an integer k.
# Clamps values to [0, total] and handles invalid inputs safely.
.gene_group_top_k <- function(total, mode, value) {
  if(total <= 0L) {
    return(0L)
  }

  value <- suppressWarnings(as.numeric(value)[1])
  if(is.na(value)) {
    return(0L)
  }

  if(mode == "pct") {
    value <- max(0, min(100, value))
    k <- ceiling(total * value / 100)
  } else {
    k <- floor(value)
  }

  k <- max(0L, min(total, as.integer(k)))
  k
}

# Check whether x is a non-empty list of data frames.
# Used for validating annotation and contrast container structures.
.gene_group_is_named_df_list <- function(x) {
  is.list(x) && length(x) > 0L && all(vapply(x, is.data.frame, logical(1)))
}

# Check whether x is a non-empty list of matrix/data.frame expression objects.
# This accepts both matrix and data.frame expression representations.
.gene_group_is_named_matrix_list <- function(x) {
  is.list(x) && length(x) > 0L &&
    all(vapply(x, function(y) is.matrix(y) || is.data.frame(y), logical(1)))
}

# Check whether one dataset-level contrast object is valid.
# A valid contrast dataset is a non-empty list of data frames.
.gene_group_is_contrast_dataset <- function(x) {
  is.list(x) && length(x) > 0L && all(vapply(x, is.data.frame, logical(1)))
}

# Ensure a list has non-empty names using a prefix for missing entries.
# Keeps list order unchanged while filling invalid names.
.gene_group_ensure_names <- function(x, default_prefix) {
  if(is.null(names(x))) {
    names(x) <- paste0(default_prefix, seq_along(x))
  }
  names(x) <- as.character(names(x))
  names(x)[is.na(names(x)) | trimws(names(x)) == ""] <- paste0(default_prefix, seq_along(x))[is.na(names(x)) | trimws(names(x)) == ""]
  x
}

# Normalize annotation input to a named list keyed by dataset.
# Accepts a single data frame or a list of data frames.
.gene_group_normalize_annot <- function(annot) {
  if(is.data.frame(annot)) {
    return(list(default=annot))
  }
  if(.gene_group_is_named_df_list(annot)) {
    annot <- .gene_group_ensure_names(annot, "dataset_")
    return(annot)
  }
  stop("`annot` must be a data frame or a named list of data frames.")
}

# Normalize expression input to a named list aligned to datasets.
# Supports single matrix/data.frame or per-dataset list input.
.gene_group_normalize_exprs <- function(exprs, datasets) {
  if(is.null(exprs)) {
    return(NULL)
  }

  if(is.matrix(exprs) || is.data.frame(exprs)) {
    if(length(datasets) != 1L) {
      stop("If multiple datasets are present, `exprs` must be a named list.")
    }
    return(stats::setNames(list(exprs), datasets))
  }

  if(!.gene_group_is_named_matrix_list(exprs)) {
    stop("`exprs` must be NULL, a matrix/data frame, or a named list of matrices/data frames.")
  }

  exprs <- .gene_group_ensure_names(exprs, "dataset_")
  if(length(datasets) == 1L && length(exprs) == 1L && !datasets[1] %in% names(exprs)) {
    names(exprs) <- datasets
  }
  missing_ds <- setdiff(datasets, names(exprs))
  if(length(missing_ds) > 0L) {
    stop(sprintf("`exprs` is missing dataset(s): %s", paste(missing_ds, collapse=", ")))
  }

  exprs[datasets]
}

# Normalize contrast input to a named list aligned to datasets.
# Supports single contrast list or per-dataset nested contrast lists.
.gene_group_normalize_cntr <- function(cntr, datasets) {
  if(is.null(cntr)) {
    return(NULL)
  }

  if(.gene_group_is_contrast_dataset(cntr)) {
    if(length(datasets) != 1L) {
      stop("If multiple datasets are present, `cntr` must be a named list of contrast lists.")
    }
    return(stats::setNames(list(cntr), datasets))
  }

  if(!is.list(cntr) || length(cntr) == 0L || !all(vapply(cntr, .gene_group_is_contrast_dataset, logical(1)))) {
    stop("`cntr` must be NULL, a contrast list, or a named list of contrast lists.")
  }

  cntr <- .gene_group_ensure_names(cntr, "dataset_")
  if(length(datasets) == 1L && length(cntr) == 1L && !datasets[1] %in% names(cntr)) {
    names(cntr) <- datasets
  }
  missing_ds <- setdiff(datasets, names(cntr))
  if(length(missing_ds) > 0L) {
    stop(sprintf("`cntr` is missing dataset(s): %s", paste(missing_ds, collapse=", ")))
  }

  cntr[datasets]
}

# Pick the first available preferred numeric column name.
# Falls back to the first candidate column when no preferred names match.
.gene_group_guess_numeric_col <- function(cols, preferred) {
  hit <- preferred[preferred %in% cols]
  if(length(hit) > 0L) {
    return(hit[1])
  }
  if(length(cols) > 0L) {
    cols[1]
  } else {
    NA_character_
  }
}

# Compute numeric columns shared across selected contrast tables.
# Used to build robust DGE UI choices that exist in every selected contrast.
.gene_group_dge_numeric_cols <- function(cntr_ds, selected_contrasts) {
  selected_contrasts <- intersect(selected_contrasts, names(cntr_ds))
  if(length(selected_contrasts) == 0L) {
    selected_contrasts <- names(cntr_ds)
  }
  if(length(selected_contrasts) == 0L) {
    return(character(0))
  }

  numeric_cols <- NULL
  for(nm in selected_contrasts) {
    cn <- names(cntr_ds[[nm]])[vapply(cntr_ds[[nm]], is.numeric, logical(1))]
    if(is.null(numeric_cols)) {
      numeric_cols <- cn
    } else {
      numeric_cols <- intersect(numeric_cols, cn)
    }
  }

  numeric_cols %||% character(0)
}

# Build available selector modes based on provided inputs.
# "By Expression" and "By DGE" are added only when corresponding data exists.
.gene_group_modes <- function(has_exprs, has_cntr) {
  md <- c("By Name"="by_name")
  if(has_exprs) {
    md <- c(md, "By Expression"="by_expression")
  }
  if(has_cntr) {
    md <- c(md, "By DGE"="by_dge")
  }
  md
}

# Create an empty annotation data frame with preserved columns.
# Used as a consistent "no selection" return object.
.gene_group_empty_annot <- function(annot_df) {
  annot_df[0, , drop=FALSE]
}

# Pick preferred column when present, otherwise return the first column.
# Keeps selector logic stable across optional UI choices.
.gene_group_pick_col <- function(df, preferred) {
  if(length(preferred) > 0L && !is.null(preferred) && preferred %in% colnames(df)) {
    return(preferred)
  }
  colnames(df)[1]
}

# Build dataset selector UI.
# Hidden when only one dataset is available.
.gene_group_dataset_input_ui <- function(ns, datasets) {
  if(length(datasets) < 2L) {
    return(hidden(selectizeInput(ns("dataset"), "Dataset", choices=datasets, selected=datasets[1])))
  }
  selectizeInput(ns("dataset"), "Dataset", choices=datasets, selected=datasets[1])
}

# Controls for name-based selection mode.
.gene_group_name_controls_ui <- function(ns, ann_cols, id_default) {
  fluidRow(
    column(
      6,
      selectInput(ns("name_id_col"), "Identifier column", choices=ann_cols, selected=id_default),
      shiny::fileInput(ns("name_file"), "Upload text file with gene IDs", accept=c(".txt", ".csv", ".tsv"))
    ),
    column(
      6,
      shiny::textAreaInput(ns("name_list"), "Gene names / IDs (space or comma separated)", rows=8)
    )
  )
}

# Controls for expression-based selection mode.
.gene_group_expression_controls_ui <- function(ns, ann_cols, expr_id_default, top_mode) {
  top_value_ui <- if(top_mode == "pct") {
    numericInput(ns("expr_top_value"), "Top percentage", value=10, min=0, max=100, step=1)
  } else {
    numericInput(ns("expr_top_value"), "Top N genes", value=100, min=1, step=1)
  }

  fluidRow(
    column(
      6,
      selectInput(ns("expr_id_col"), "Annotation column matching expression IDs", choices=ann_cols, selected=expr_id_default),
      selectInput(ns("expr_metric"), "Ranking metric", choices=c("Top variance"="variance", "Top average expression"="mean"))
    ),
    column(
      6,
      selectInput(ns("expr_top_mode"), "Top in", choices=c("Absolute number"="n", "Percentage"="pct"), selected=top_mode),
      top_value_ui
    )
  )
}

# Compute DGE mode defaults from available contrasts and current selections.
.gene_group_dge_defaults <- function(cntr_ds, selected_contrasts, selected_p_col=NULL, selected_lfc_col=NULL, selected_fdr_col=NULL) {
  cntr_names <- names(cntr_ds)
  sel_cntr <- intersect(selected_contrasts, cntr_names)
  if(length(sel_cntr) == 0L) {
    sel_cntr <- cntr_names[1]
  }

  num_cols <- .gene_group_dge_numeric_cols(cntr_ds, sel_cntr)
  p_col <- .gene_group_guess_numeric_col(num_cols, c("padj", "pvalue", "P.Value", "PValue"))
  lfc_col <- .gene_group_guess_numeric_col(num_cols, c("log2FoldChange", "logFC", "lfc", "LFC"))
  fdr_col <- .gene_group_guess_numeric_col(num_cols, c("padj", "adj.P.Val", "FDR", "qvalue", "q_value"))

  if(isTruthy(selected_p_col) && selected_p_col %in% num_cols) {
    p_col <- selected_p_col
  }
  if(isTruthy(selected_lfc_col) && selected_lfc_col %in% num_cols) {
    lfc_col <- selected_lfc_col
  }
  if(isTruthy(selected_fdr_col) && selected_fdr_col %in% num_cols) {
    fdr_col <- selected_fdr_col
  }

  list(
    contrasts=cntr_names,
    selected_contrasts=sel_cntr,
    numeric_cols=num_cols,
    p_col=p_col,
    lfc_col=lfc_col,
    fdr_col=fdr_col
  )
}

# Controls for DGE-based selection mode.
.gene_group_dge_controls_ui <- function(ns, ann_cols, dge_id_default, dge_defaults, top_mode, rank_mode) {
  top_value_ui <- NULL
  if(top_mode == "n") {
    top_value_ui <- numericInput(ns("dge_top_value"), "Top N genes", value=100, min=1, step=1)
  } else if(top_mode == "pct") {
    top_value_ui <- numericInput(ns("dge_top_value"), "Top percentage", value=10, min=0, max=100, step=1)
  }

  fluidRow(
    column(
      6,
      selectizeInput(ns("dge_contrasts"), "Contrasts", choices=dge_defaults$contrasts, selected=dge_defaults$selected_contrasts, multiple=FALSE),
      selectInput(ns("dge_id_col"), "Annotation column matching contrast IDs", choices=ann_cols, selected=dge_id_default),
      selectInput(ns("dge_pval_col"), "P-value column", choices=dge_defaults$numeric_cols, selected=dge_defaults$p_col),
      numericInput(ns("dge_pval_thr"), "P-value threshold", value=.05, min=0, max=1, step=.001),
      selectInput(ns("dge_fdr_col"), "Adjusted p-value (FDR) column", choices=dge_defaults$numeric_cols, selected=dge_defaults$fdr_col),
      checkboxInput(ns("dge_require_fdr"), "Exclude genes with missing adjusted p-value (FDR)", value=FALSE)
    ),
    column(
      6,
      selectInput(ns("dge_lfc_col"), "logFC column", choices=dge_defaults$numeric_cols, selected=dge_defaults$lfc_col),
      numericInput(ns("dge_lfc_thr"), "Absolute logFC threshold", value=0, min=0, step=.1),
      selectInput(ns("dge_rank_mode"), "Ranking priority", choices=c("Lowest p-value first"="p_first", "Highest absolute logFC first"="lfc_first"), selected=rank_mode),
      selectInput(ns("dge_top_mode"), "Top filter", choices=c("All passing filters"="all", "Top N"="n", "Top %"="pct"), selected=top_mode),
      top_value_ui
    )
  )
}

# Build mode-specific controls for the selector UI.
.gene_group_modus_controls_ui <- function(ns, mode, ann_df, primary_id, expr_ds=NULL, cntr_ds=NULL, input_values=list()) {
  ann_cols <- colnames(ann_df)

  if(mode == "by_name") {
    id_default <- .gene_group_pick_col(ann_df, primary_id)
    return(.gene_group_name_controls_ui(ns, ann_cols, id_default))
  }

  if(mode == "by_expression") {
    if(is.null(expr_ds)) {
      return(p("Expression data not available."))
    }
    expr_id_default <- .gene_group_pick_col(ann_df, primary_id)
    top_mode <- input_values$expr_top_mode %||% "n"
    return(.gene_group_expression_controls_ui(ns, ann_cols, expr_id_default, top_mode))
  }

  if(mode == "by_dge") {
    if(is.null(cntr_ds)) {
      return(p("Contrast data not available."))
    }
    cntr_names <- names(cntr_ds)
    if(length(cntr_names) == 0L) {
      return(p("No contrasts available for this dataset."))
    }

    dge_defaults <- .gene_group_dge_defaults(
      cntr_ds=cntr_ds,
      selected_contrasts=input_values$dge_contrasts,
      selected_p_col=input_values$dge_pval_col,
      selected_lfc_col=input_values$dge_lfc_col,
      selected_fdr_col=input_values$dge_fdr_col
    )
    top_mode <- input_values$dge_top_mode %||% "all"
    rank_mode <- input_values$dge_rank_mode %||% "p_first"
    dge_id_default <- .gene_group_pick_col(ann_df, primary_id)
    return(.gene_group_dge_controls_ui(ns, ann_cols, dge_id_default, dge_defaults, top_mode, rank_mode))
  }

  p("No selection mode available.")
}

# Match query IDs to annotation rows using one annotation column.
.gene_group_match_annot_rows <- function(ann_df, id_col, query_ids) {
  empty_df <- .gene_group_empty_annot(ann_df)
  if(length(query_ids) == 0L) {
    return(empty_df)
  }

  ann_ids <- as.character(ann_df[[id_col]])
  m <- match(query_ids, ann_ids)
  m <- m[!is.na(m)]
  if(length(m) == 0L) {
    return(empty_df)
  }

  ann_df[m, , drop=FALSE]
}

# Select annotation rows from free text and uploaded gene list file.
.gene_group_select_by_name <- function(ann_df, id_col, name_list, name_file) {
  name_ids <- .gene_group_parse_gene_input(name_list)
  file_ids <- .gene_group_read_gene_file(name_file)
  query_ids <- unique(c(name_ids, file_ids))
  .gene_group_match_annot_rows(ann_df, id_col, query_ids)
}

# Select annotation rows based on expression ranking.
.gene_group_select_by_expression <- function(ann_df, expr_ds, id_col, metric="variance", top_mode="n", top_value=100) {
  empty_df <- .gene_group_empty_annot(ann_df)
  if(is.null(expr_ds)) {
    return(empty_df)
  }

  expr_ds <- as.matrix(expr_ds)
  if(nrow(expr_ds) < 1L || is.null(rownames(expr_ds))) {
    return(empty_df)
  }

  score <- if(metric == "mean") {
    rowMeans(expr_ds, na.rm=TRUE)
  } else {
    apply(expr_ds, 1, stats::var, na.rm=TRUE)
  }
  score <- score[!is.na(score)]
  if(length(score) == 0L) {
    return(empty_df)
  }

  k <- .gene_group_top_k(length(score), top_mode, top_value)
  if(k < 1L) {
    return(empty_df)
  }

  ids <- names(sort(score, decreasing=TRUE))[seq_len(k)]
  .gene_group_match_annot_rows(ann_df, id_col, ids)
}

# Build a compact per-gene DGE summary from selected contrasts.
.gene_group_collect_dge_rows <- function(cntr_ds, contrasts, cntr_id_col, p_col, lfc_col, fdr_col, require_fdr, p_thr, lfc_thr) {
  dge_rows <- lapply(contrasts, function(nm) {
    df <- cntr_ds[[nm]]
    required_cols <- c(cntr_id_col, p_col, lfc_col)
    if(require_fdr) {
      required_cols <- c(required_cols, fdr_col)
    }
    if(!all(required_cols %in% colnames(df))) {
      return(NULL)
    }

    ids <- as.character(df[[cntr_id_col]])
    pv <- suppressWarnings(as.numeric(df[[p_col]]))
    lfc <- suppressWarnings(as.numeric(df[[lfc_col]]))
    fdr <- if(require_fdr) suppressWarnings(as.numeric(df[[fdr_col]])) else rep(NA_real_, nrow(df))

    keep <- !is.na(ids) & ids != "" & !is.na(pv) & !is.na(lfc)
    if(require_fdr) {
      keep <- keep & !is.na(fdr)
    }
    keep <- keep & pv <= p_thr & abs(lfc) >= lfc_thr
    if(!any(keep)) {
      return(NULL)
    }

    data.frame(
      id=ids[keep],
      pvalue=pv[keep],
      abs_lfc=abs(lfc[keep]),
      stringsAsFactors=FALSE
    )
  })

  dge_rows <- dge_rows[!vapply(dge_rows, is.null, logical(1))]
  if(length(dge_rows) == 0L) {
    return(NULL)
  }

  do.call(rbind, dge_rows)
}

# Rank genes by strongest DGE evidence across contrasts.
.gene_group_rank_dge_rows <- function(dge_rows, rank_mode="p_first") {
  split_ids <- split(dge_rows, dge_rows$id)
  rank_df <- data.frame(
    id=names(split_ids),
    min_p=vapply(split_ids, function(x) min(x$pvalue, na.rm=TRUE), numeric(1)),
    max_abs_lfc=vapply(split_ids, function(x) max(x$abs_lfc, na.rm=TRUE), numeric(1)),
    stringsAsFactors=FALSE
  )

  if(rank_mode == "lfc_first") {
    rank_df[order(-rank_df$max_abs_lfc, rank_df$min_p, rank_df$id), , drop=FALSE]
  } else {
    rank_df[order(rank_df$min_p, -rank_df$max_abs_lfc, rank_df$id), , drop=FALSE]
  }
}

# Select annotation rows based on DGE filters/ranking.
.gene_group_select_by_dge <- function(ann_df, cntr_ds, id_col, cntr_id_col,
                                      contrasts, p_col, lfc_col, fdr_col,
                                      require_fdr=FALSE, p_thr=.05, lfc_thr=0,
                                      rank_mode="p_first", top_mode="all", top_value=100) {
  empty_df <- .gene_group_empty_annot(ann_df)
  if(is.null(cntr_ds) || length(cntr_ds) < 1L) {
    return(empty_df)
  }

  contrasts <- intersect(contrasts, names(cntr_ds))
  if(length(contrasts) == 0L) {
    return(empty_df)
  }

  if(!isTruthy(p_col) || !isTruthy(lfc_col) || (isTRUE(require_fdr) && !isTruthy(fdr_col))) {
    return(empty_df)
  }

  p_thr <- suppressWarnings(as.numeric(p_thr)[1])
  lfc_thr <- suppressWarnings(as.numeric(lfc_thr)[1])
  if(is.na(p_thr) || is.na(lfc_thr)) {
    return(empty_df)
  }

  dge_rows <- .gene_group_collect_dge_rows(
    cntr_ds=cntr_ds,
    contrasts=contrasts,
    cntr_id_col=cntr_id_col,
    p_col=p_col,
    lfc_col=lfc_col,
    fdr_col=fdr_col,
    require_fdr=isTRUE(require_fdr),
    p_thr=p_thr,
    lfc_thr=lfc_thr
  )
  if(is.null(dge_rows) || nrow(dge_rows) == 0L) {
    return(empty_df)
  }

  rank_df <- .gene_group_rank_dge_rows(dge_rows, rank_mode=rank_mode)
  if(top_mode %in% c("n", "pct")) {
    k <- .gene_group_top_k(nrow(rank_df), top_mode, top_value)
    if(k < 1L) {
      return(empty_df)
    }
    rank_df <- rank_df[seq_len(k), , drop=FALSE]
  }

  .gene_group_match_annot_rows(ann_df, id_col, rank_df$id)
}

# Dispatch selection mode to the corresponding pure selection helper.
.gene_group_selected_annotation <- function(mode, ann_df, cntr_id_col,
                                            expr_ds=NULL, cntr_ds=NULL, input_values=list()) {
  if(mode == "by_name") {
    id_col <- .gene_group_pick_col(ann_df, input_values$name_id_col)
    return(.gene_group_select_by_name(ann_df, id_col, input_values$name_list, input_values$name_file))
  }

  if(mode == "by_expression") {
    id_col <- .gene_group_pick_col(ann_df, input_values$expr_id_col)
    metric <- input_values$expr_metric %||% "variance"
    top_mode <- input_values$expr_top_mode %||% "n"
    top_value <- input_values$expr_top_value %||% if(top_mode == "pct") 10 else 100
    return(.gene_group_select_by_expression(ann_df, expr_ds, id_col, metric, top_mode, top_value))
  }

  if(mode == "by_dge") {
    id_col <- .gene_group_pick_col(ann_df, input_values$dge_id_col)
    rank_mode <- input_values$dge_rank_mode %||% "p_first"
    top_mode <- input_values$dge_top_mode %||% "all"
    top_value <- input_values$dge_top_value %||% if(top_mode == "pct") 10 else 100
    return(.gene_group_select_by_dge(
      ann_df=ann_df,
      cntr_ds=cntr_ds,
      id_col=id_col,
      cntr_id_col=cntr_id_col,
      contrasts=input_values$dge_contrasts,
      p_col=input_values$dge_pval_col,
      lfc_col=input_values$dge_lfc_col,
      fdr_col=input_values$dge_fdr_col,
      require_fdr=input_values$dge_require_fdr,
      p_thr=input_values$dge_pval_thr,
      lfc_thr=input_values$dge_lfc_thr,
      rank_mode=rank_mode,
      top_mode=top_mode,
      top_value=top_value
    ))
  }

  .gene_group_empty_annot(ann_df)
}

#' @rdname geneGroupSelectorServer
#' @export
# Build controls-only UI for selecting gene groups.
# Intended to be placed inside a parent sidebarPanel.
geneGroupSelectorUI <- function(id) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("dataset_ui")),
    selectizeInput(ns("modus"), "Modus", choices=character(0), selected=NULL),
    uiOutput(ns("modus_controls")),
    br(),
    h4("Current selection"),
    textOutput(ns("selected_count")),
    textOutput(ns("selected_preview")),
    downloadButton(ns("save_genes"), "Export selected PrimaryIDs")
  )
}

#' Shiny module for selecting groups of genes
#'
#' The module supports three selection modi:
#' 1) by name (`data_type = by_name`) via free text and/or file upload;
#' 2) by expression statistics (top variance or top average expression);
#' 3) by differential expression results from contrast tables.
#'
#' In DGE mode, users can choose ranking priority (`lowest p-value first`
#' or `highest absolute logFC first`) and optionally exclude rows with
#' missing adjusted p-values (FDR).
#'
#' @param id module identifier (same as passed to [geneGroupSelectorUI()])
#' @param annot annotation data frame (minimum input), or named list of data frames
#'   for multiple datasets
#' @param exprs optional expression matrix/data frame, or named list of matrices/data
#'   frames for multiple datasets
#' @param cntr optional contrasts list (named list of data frames), or named list of
#'   such lists for multiple datasets
#' @param primary_id default annotation column used as the primary gene identifier
#' @param cntr_id_col contrast table column containing gene identifiers
#' @param selected_ids optional `reactiveVal()` that will be updated
#'   with the currently selected PrimaryIDs (column `primary_id`) as a
#'   character vector
#' @return A list with reactives: `genes`, `annotation`, `dataset`, and `modus`.
#' @examples
#' ## Minimal example
#' annot <- data.frame(
#'   PrimaryID=paste0("g", 1:5),
#'   SYMBOL=LETTERS[1:5]
#' )
#' exprs <- matrix(
#'   rnorm(20),
#'   nrow=5,
#'   dimnames=list(annot$PrimaryID, paste0("s", 1:4))
#' )
#' cntr <- list(
#'   contrast_a=data.frame(
#'     PrimaryID=annot$PrimaryID,
#'     padj=c(0.01, 0.20, 0.03, NA, 0.50),
#'     pvalue=c(0.001, 0.10, 0.02, 0.04, 0.80),
#'     log2FoldChange=c(2.1, -0.2, 1.3, 0.8, -3.2)
#'   )
#' )
#' if(interactive()) {
#'   ui <- fluidPage(
#'     sidebarLayout(
#'       sidebarPanel(geneGroupSelectorUI("gsel")),
#'       mainPanel(tableOutput("selected_ids"))
#'     )
#'   )
#'   server <- function(input, output, session) {
#'     selected_ids <- reactiveVal(character())
#'
#'     geneGroupSelectorServer(
#'       "gsel",
#'       annot=annot,
#'       exprs=exprs,
#'       cntr=cntr,
#'       selected_ids=selected_ids
#'     )
#'
#'     output$selected_ids <- renderTable({
#'       data.frame(PrimaryID=selected_ids(), stringsAsFactors=FALSE)
#'     })
#'   }
#'   shinyApp(ui, server)
#' }
#'
#' ## Example with C19 dataset
#' data(C19)
#' if(interactive()) {
#'   ui <- fluidPage(
#'     sidebarLayout(
#'       sidebarPanel(geneGroupSelectorUI("gsel")),
#'       mainPanel(tableOutput("selected_ids"))
#'     )
#'   )
#'   server <- function(input, output, session) {
#'     selected_ids <- reactiveVal(character())
#'
#'     geneGroupSelectorServer(
#'       "gsel",
#'       annot=C19$annotation,
#'       exprs=C19$expression,
#'       cntr=C19$contrasts,
#'       selected_ids=selected_ids
#'     )
#'
#'     output$selected_ids <- renderTable({
#'       data.frame(PrimaryID=selected_ids(), stringsAsFactors=FALSE)
#'     })
#'   }
#'   shinyApp(ui, server)
#' }
#' @export
# Register server logic for mode-specific gene selection and filtering.
# Returns reactives for selected PrimaryIDs, annotation subset, dataset, and mode.
geneGroupSelectorServer <- function(id, annot, exprs=NULL, cntr=NULL,
                                    primary_id="PrimaryID",
                                    cntr_id_col=primary_id,
                                    selected_ids=NULL) {
  annot <- .gene_group_normalize_annot(annot)
  for(ds in names(annot)) {
    if(!primary_id %in% colnames(annot[[ds]])) {
      stop(sprintf("`primary_id` ('%s') not found in annotation for dataset '%s'.", primary_id, ds))
    }
  }

  if(!is.null(selected_ids)) {
    if(!inherits(selected_ids, "reactiveVal")) {
      stop("`selected_ids` must be NULL or a `reactiveVal()`.")
    }
  }

  datasets <- names(annot)
  exprs <- .gene_group_normalize_exprs(exprs, datasets)
  cntr <- .gene_group_normalize_cntr(cntr, datasets)
  modes <- .gene_group_modes(!is.null(exprs), !is.null(cntr))

  moduleServer(id, function(input, output, session) {
    output$dataset_ui <- renderUI({
      .gene_group_dataset_input_ui(session$ns, datasets)
    })

    observeEvent(TRUE, {
      shiny::updateSelectizeInput(
        session=session,
        inputId="modus",
        choices=modes,
        selected=modes[1],
        server=TRUE
      )
    }, once=TRUE)

    dataset_selected <- reactive({
      if(length(datasets) == 1L) {
        return(datasets[1])
      }
      .ds <- input$dataset
      if(!isTruthy(.ds) || !.ds %in% datasets) {
        return(datasets[1])
      }
      .ds
    })

    output$modus_controls <- renderUI({
      ds <- dataset_selected()
      mode <- input$modus %||% modes[1]
      .gene_group_modus_controls_ui(
        ns=session$ns,
        mode=mode,
        ann_df=annot[[ds]],
        primary_id=primary_id,
        expr_ds=if(is.null(exprs)) NULL else exprs[[ds]],
        cntr_ds=if(is.null(cntr)) NULL else cntr[[ds]],
        input_values=list(
          expr_top_mode=input$expr_top_mode,
          dge_contrasts=input$dge_contrasts,
          dge_pval_col=input$dge_pval_col,
          dge_lfc_col=input$dge_lfc_col,
          dge_fdr_col=input$dge_fdr_col,
          dge_top_mode=input$dge_top_mode,
          dge_rank_mode=input$dge_rank_mode
        )
      )
    })

    selected_annotation <- reactive({
      ds <- dataset_selected()
      mode <- input$modus %||% modes[1]
      .gene_group_selected_annotation(
        mode=mode,
        ann_df=annot[[ds]],
        cntr_id_col=cntr_id_col,
        expr_ds=if(is.null(exprs)) NULL else exprs[[ds]],
        cntr_ds=if(is.null(cntr)) NULL else cntr[[ds]],
        input_values=list(
          name_id_col=input$name_id_col,
          name_list=input$name_list,
          name_file=input$name_file,
          expr_id_col=input$expr_id_col,
          expr_metric=input$expr_metric,
          expr_top_mode=input$expr_top_mode,
          expr_top_value=input$expr_top_value,
          dge_contrasts=input$dge_contrasts,
          dge_id_col=input$dge_id_col,
          dge_pval_col=input$dge_pval_col,
          dge_lfc_col=input$dge_lfc_col,
          dge_fdr_col=input$dge_fdr_col,
          dge_require_fdr=input$dge_require_fdr,
          dge_pval_thr=input$dge_pval_thr,
          dge_lfc_thr=input$dge_lfc_thr,
          dge_rank_mode=input$dge_rank_mode,
          dge_top_mode=input$dge_top_mode,
          dge_top_value=input$dge_top_value
        )
      )
    })

    selected_primary_ids <- reactive({
      sel <- selected_annotation()
      if(nrow(sel) == 0L) {
        return(character(0))
      }

      ids <- as.character(sel[[primary_id]])
      ids <- ids[!is.na(ids) & ids != ""]
      unique(ids)
    })

    observe({
      if(!is.null(selected_ids)) {
        selected_ids(selected_primary_ids())
      }
    })

    output$selected_count <- renderText({
      sprintf("Selected genes: %d", length(selected_primary_ids()))
    })

    output$selected_preview <- renderText({
      ids <- selected_primary_ids()
      if(length(ids) == 0L) {
        return("First PrimaryIDs: none")
      }
      preview <- utils::head(ids, 10)
      suffix <- if(length(ids) > 10L) ", ..." else ""
      sprintf("First PrimaryIDs: %s%s", paste(preview, collapse=", "), suffix)
    })

    output$save_genes <- downloadHandler(
      filename = function() {
        ds <- dataset_selected()
        md <- input$modus %||% modes[1]
        sprintf(
          "selected_genes_%s_%s.txt",
          .gene_group_sanitize_filename(ds, "dataset"),
          .gene_group_sanitize_filename(md, "modus")
        )
      },
      content = function(file) {
        writeLines(selected_primary_ids(), con=file)
      }
    )

    return(list(
      genes=selected_primary_ids,
      annotation=selected_annotation,
      dataset=dataset_selected,
      modus=reactive(input$modus %||% modes[1])
    ))
  })
}
