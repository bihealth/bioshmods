## ---------------------------------------------
## General Utilities
## ---------------------------------------------

# Write module-scoped log messages through the shared package logger.
# Keeps the selector log prefix consistent across helper functions.
.gene_group_log <- function(...) {
  .bioshmods_log(..., .prefix="gene_group_selector")
}

# Parse comma- or space-separated gene identifiers from text input.
# Returns unique non-empty tokens while preserving their input order.
.gene_group_parse_gene_input <- function(x) {
  x <- x %||% ""
  x <- paste(x, collapse=" ")
  tokens <- unlist(strsplit(x, "[,[:space:]]+"))
  tokens <- trimws(tokens)
  tokens <- tokens[!is.na(tokens) & tokens != ""]
  tokens <- unique(tokens)
  .gene_group_log("parsed gene input tokens n=", as.character(length(tokens)), ".")
  tokens
}

# Read a gene list file uploaded through shiny::fileInput.
# Returns parsed identifiers, or character(0) when input/file is missing.
.gene_group_read_gene_file <- function(file_input) {
  if(is.null(file_input) || is.null(file_input$datapath) || !nzchar(file_input$datapath)) {
    .gene_group_log("gene file input missing; returning empty list.")
    return(character(0))
  }
  if(!file.exists(file_input$datapath)) {
    .gene_group_log("gene file path does not exist: ", file_input$datapath)
    return(character(0))
  }
  .gene_group_log("reading gene file from: ", file_input$datapath)
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


## ---------------------------------------------
## Sanity checks
## ---------------------------------------------

# Check that the specified primary ID column exists in all annotation data frames.
.gene_group_check_primary_id <- function(annot, primary_id) {
  for(ds in names(annot)) {
    if(!primary_id %in% colnames(annot[[ds]])) {
      stop(sprintf("`primary_id` ('%s') not found in annotation for dataset '%s'.", primary_id, ds))
    }
  }
}


# Check that "selected_ids" is either NULL or a reactiveVal object.
.gene_group_check_selected_ids <- function(selected_ids) {
  if(!is.null(selected_ids)) {
    if(!inherits(selected_ids, "reactiveVal")) {
      stop("`selected_ids` must be NULL or a `reactiveVal()`.")
    }
  }
}

# Check that "dataset" is either NULL or a reactiveVal object.
.gene_group_check_dataset <- function(dataset) {
  if(!is.null(dataset) && !inherits(dataset, "reactiveVal")) {
    stop("`dataset` must be NULL or a `reactiveVal()`.")
  }
}



## ---------------------------------------------
## Input Normalization
## ---------------------------------------------

# Normalize annotation input to a named list keyed by dataset.
# Accepts a single data frame or a list of data frames.
.gene_group_normalize_annot <- function(annot) {
  .gene_group_log("normalize annot input class=", paste(class(annot), collapse="/"), ".")
  if(is.data.frame(annot)) {
    .gene_group_log("annot provided as data.frame; assigning dataset 'default'.")
    return(list(default=annot))
  }
  if(.gene_group_is_named_df_list(annot)) {
    annot <- .gene_group_ensure_names(annot, "dataset_")
    .gene_group_log("annot datasets={", paste(names(annot), collapse=","), "}.")
    return(annot)
  }
  stop("`annot` must be a data frame or a named list of data frames.")
}

# Normalize expression input to a named list aligned to datasets.
# Supports single matrix/data.frame or per-dataset list input.
.gene_group_normalize_exprs <- function(exprs, datasets) {
  .gene_group_log("normalize exprs input class=", paste(class(exprs), collapse="/"),
                  "; expected datasets={", paste(datasets, collapse=","), "}.")
  if(is.null(exprs)) {
    .gene_group_log("exprs is NULL.")
    return(NULL)
  }

  if(is.matrix(exprs) || is.data.frame(exprs)) {
    if(length(datasets) != 1L) {
      stop("If multiple datasets are present, `exprs` must be a named list.")
    }
    .gene_group_log("exprs provided as matrix/data.frame for single dataset '", datasets, "'.")
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
    .gene_group_log("exprs missing datasets: ", paste(missing_ds, collapse=","))
    stop(sprintf("`exprs` is missing dataset(s): %s", paste(missing_ds, collapse=", ")))
  }

  .gene_group_log("normalized exprs datasets={", paste(datasets, collapse=","), "}.")
  exprs[datasets]
}

# Normalize contrast input to a named list aligned to datasets.
# Supports single contrast list or per-dataset nested contrast lists.
.gene_group_normalize_cntr <- function(cntr, datasets) {
  .gene_group_log("normalize cntr input class=", paste(class(cntr), collapse="/"),
                  "; expected datasets={", paste(datasets, collapse=","), "}.")
  if(is.null(cntr)) {
    .gene_group_log("cntr is NULL.")
    return(NULL)
  }

  if(.gene_group_is_contrast_dataset(cntr)) {
    if(length(datasets) != 1L) {
      stop("If multiple datasets are present, `cntr` must be a named list of contrast lists.")
    }
    .gene_group_log("cntr provided as contrast list for single dataset '", datasets, "'.")
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
    .gene_group_log("cntr missing datasets: ", paste(missing_ds, collapse=","))
    stop(sprintf("`cntr` is missing dataset(s): %s", paste(missing_ds, collapse=", ")))
  }

  .gene_group_log("normalized cntr datasets={", paste(datasets, collapse=","), "}.")
  cntr[datasets]
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

# Standard selection mode identifiers.
# Used for validating mode order customization.
.gene_group_allowed_modes <- function() {
  c("by_expression", "by_dge", "by_name")
}

# Normalize user-provided mode ordering.
# Keeps only known mode ids, removes duplicates, and appends missing modes.
.gene_group_normalize_mode_order <- function(mode_order) {
  allowed <- .gene_group_allowed_modes()
  if(is.null(mode_order)) {
    return(allowed)
  }

  if(!is.character(mode_order)) {
    stop("`mode_order` must be NULL or a character vector.")
  }

  mode_order <- trimws(as.character(mode_order))
  mode_order <- mode_order[!is.na(mode_order) & mode_order != ""]
  mode_order <- unique(mode_order[mode_order %in% allowed])
  c(mode_order, setdiff(allowed, mode_order))
}

# Default values used to initialize and evaluate selector controls.
# Can be overridden via `geneGroupSelectorServer(..., defaults=...)`.
.gene_group_default_settings <- function() {
  list(
    expr_metric="mean",
    expr_top_mode="n",
    expr_top_value=50,
    dge_top_mode="all",
    dge_top_value=100,
    dge_rank_mode="p_first",
    dge_pval_thr=.05,
    dge_lfc_thr=0,
    dge_require_fdr=TRUE
  )
}

# Normalize optional selector defaults supplied by users.
# Invalid values raise explicit errors.
.gene_group_normalize_defaults <- function(defaults) {
  base <- .gene_group_default_settings()
  if(is.null(defaults)) {
    return(base)
  }

  if(!is.list(defaults)) {
    stop("`defaults` must be NULL or a list.")
  }

  extra_keys <- setdiff(names(defaults), names(base))
  if(length(extra_keys) > 0L) {
    stop(sprintf(
      "`defaults` contains unsupported key(s): %s",
      paste(extra_keys, collapse=", ")
    ))
  }

  x <- utils::modifyList(base, defaults)

  x$expr_metric <- as.character(x$expr_metric)[1]
  if(is.na(x$expr_metric) || !x$expr_metric %in% c("variance", "mean")) {
    stop("`defaults$expr_metric` must be one of: 'variance', 'mean'.")
  }

  x$expr_top_mode <- as.character(x$expr_top_mode)[1]
  if(is.na(x$expr_top_mode) || !x$expr_top_mode %in% c("n", "pct")) {
    stop("`defaults$expr_top_mode` must be one of: 'n', 'pct'.")
  }

  x$expr_top_value <- suppressWarnings(as.numeric(x$expr_top_value)[1])
  if(!is.finite(x$expr_top_value)) {
    stop("`defaults$expr_top_value` must be a finite number.")
  }
  if(x$expr_top_mode == "pct") {
    if(x$expr_top_value < 0 || x$expr_top_value > 100) {
      stop("`defaults$expr_top_value` must be in [0, 100] when `expr_top_mode='pct'`.")
    }
  } else {
    if(x$expr_top_value < 1 || abs(x$expr_top_value - round(x$expr_top_value)) > .Machine$double.eps^0.5) {
      stop("`defaults$expr_top_value` must be an integer >= 1 when `expr_top_mode='n'`.")
    }
    x$expr_top_value <- as.integer(round(x$expr_top_value))
  }

  x$dge_top_mode <- as.character(x$dge_top_mode)[1]
  if(is.na(x$dge_top_mode) || !x$dge_top_mode %in% c("all", "n", "pct")) {
    stop("`defaults$dge_top_mode` must be one of: 'all', 'n', 'pct'.")
  }

  x$dge_top_value <- suppressWarnings(as.numeric(x$dge_top_value)[1])
  if(!is.finite(x$dge_top_value)) {
    stop("`defaults$dge_top_value` must be a finite number.")
  }
  if(x$dge_top_mode == "pct") {
    if(x$dge_top_value < 0 || x$dge_top_value > 100) {
      stop("`defaults$dge_top_value` must be in [0, 100] when `dge_top_mode='pct'`.")
    }
  } else if(x$dge_top_mode == "n") {
    if(x$dge_top_value < 1 || abs(x$dge_top_value - round(x$dge_top_value)) > .Machine$double.eps^0.5) {
      stop("`defaults$dge_top_value` must be an integer >= 1 when `dge_top_mode='n'`.")
    }
    x$dge_top_value <- as.integer(round(x$dge_top_value))
  }

  x$dge_rank_mode <- as.character(x$dge_rank_mode)[1]
  if(is.na(x$dge_rank_mode) || !x$dge_rank_mode %in% c("p_first", "lfc_first")) {
    stop("`defaults$dge_rank_mode` must be one of: 'p_first', 'lfc_first'.")
  }

  x$dge_pval_thr <- suppressWarnings(as.numeric(x$dge_pval_thr)[1])
  if(!is.finite(x$dge_pval_thr) || x$dge_pval_thr < 0 || x$dge_pval_thr > 1) {
    stop("`defaults$dge_pval_thr` must be a number in [0, 1].")
  }

  x$dge_lfc_thr <- suppressWarnings(as.numeric(x$dge_lfc_thr)[1])
  if(!is.finite(x$dge_lfc_thr) || x$dge_lfc_thr < 0) {
    stop("`defaults$dge_lfc_thr` must be a finite number >= 0.")
  }

  if(!is.logical(x$dge_require_fdr) || length(x$dge_require_fdr) != 1L || is.na(x$dge_require_fdr)) {
    stop("`defaults$dge_require_fdr` must be TRUE or FALSE.")
  }

  x$dge_require_fdr <- as.logical(x$dge_require_fdr)
  x
}

# Resolve mode-specific values from the current parameter state.
# Uses live inputs first and falls back to defaults from the selector config.
.gene_group_resolve_mode_inputs <- function(param_state, selector_config) {
  mode <- param_state$mode
  input_values <- param_state$input_values
  defaults <- selector_config$defaults

  if(mode == "by_expression") {
    return(list(
      metric=input_values$expr_metric %||% defaults$expr_metric,
      top_mode=input_values$expr_top_mode %||% defaults$expr_top_mode,
      top_value=input_values$expr_top_value %||% defaults$expr_top_value
    ))
  }

  if(mode == "by_dge") {
    return(list(
      top_mode=input_values$dge_top_mode %||% defaults$dge_top_mode,
      top_value=input_values$dge_top_value %||% defaults$dge_top_value,
      rank_mode=input_values$dge_rank_mode %||% defaults$dge_rank_mode,
      p_thr=input_values$dge_pval_thr %||% defaults$dge_pval_thr,
      lfc_thr=input_values$dge_lfc_thr %||% defaults$dge_lfc_thr,
      require_fdr=input_values$dge_require_fdr %||% defaults$dge_require_fdr
    ))
  }

  list()
}

# Build available selector modes based on provided inputs.
# "By Expression" and "By DGE" are added only when corresponding data exists.
.gene_group_modes <- function(has_exprs, has_cntr, mode_order=.gene_group_allowed_modes()) {
  mode_labels <- c(
    by_expression="By Expression",
    by_dge="By DGE",
    by_name="By Name"
  )

  available <- c("by_name")
  if(has_exprs) {
    available <- c(available, "by_expression")
  }
  if(has_cntr) {
    available <- c(available, "by_dge")
  }

  ordered <- .gene_group_normalize_mode_order(mode_order)
  ordered <- c(intersect(ordered, available), setdiff(available, ordered))
  .gene_group_log("available modes={", paste(available, collapse=","),
                  "}; ordered modes={", paste(ordered, collapse=","), "}.")
  stats::setNames(ordered, mode_labels[ordered])
}







## ---------------------------------------------
## UI Builders
## ---------------------------------------------

# Build dataset selector UI.
# Hidden when only one dataset is available.
.gene_group_dataset_input_ui <- function(ns, datasets, selected=NULL) {
  if(!isTruthy(selected) || !selected %in% datasets) {
    selected <- datasets[1]
  }

  # hide the datasdet selector if we only have one dataset
  if(length(datasets) < 2L) {
    .gene_group_log("dataset UI hidden (single dataset='", selected, "').")
    return(hidden(selectizeInput(ns("dataset"), "Dataset", choices=datasets, selected=selected)))
  }

  .gene_group_log("dataset UI shown; datasets={", paste(datasets, collapse=","),
                  "}, selected='", selected, "'.")
  selectizeInput(ns("dataset"), "Dataset", choices=datasets, selected=selected)
}

# Controls for name-based selection mode.
.gene_group_name_controls_ui <- function(ns, ann_cols, selected_col, name_list=NULL) {
  if(!isTruthy(selected_col) || !selected_col %in% ann_cols) {
    selected_col <- ann_cols[1]
  }
  name_list <- name_list %||% ""

  fluidRow(
    column(
      12,
      selectInput(
        ns("name_id_col"),
        "Annotation column for names/IDs",
        choices=ann_cols,
        selected=selected_col
      ),
      shiny::fileInput(ns("name_file"), "Upload text file with gene IDs", accept=c(".txt", ".csv", ".tsv"))
    ),
    column(
      12,
      shiny::textAreaInput(
        ns("name_list"),
        "Gene names / IDs (space or comma separated)",
        rows=8,
        value=name_list
      )
    )
  )
}

# Build the controls for expression-based gene selection.
# Renders metric and top-k inputs for the expression mode UI.
.gene_group_expression_controls_ui <- function(ns, metric, top_mode, top_value) {
  top_value <- suppressWarnings(as.numeric(top_value)[1])
  if(!is.finite(top_value)) {
    top_value <- if(top_mode == "pct") 10 else 100
  }

  top_value_ui <- if(top_mode == "pct") {
    numericInput(
      ns("expr_top_value"),
      "",
      value=max(0, min(100, top_value)),
      min=0,
      max=100,
      step=1
    )
  } else {
    numericInput(ns("expr_top_value"), "", value=max(1, floor(top_value)), min=1, step=1)
  }

  shiny::tags$div(
    fluidRow(
      column(
        12,
        selectInput(
          ns("expr_metric"),
          "Ranking by",
          choices=c("Variance"="variance", "Average expression"="mean"),
          selected=metric
        )
      )
    ), fluidRow(
      column(
        8,
        selectInput(ns("expr_top_mode"), "", choices=c("Absolute number:"="n", "Percentage:"="pct"), selected=top_mode)
      ),
      column(
        4,
        top_value_ui
      )
    )
  )
}

# Resolve the DGE columns used by the selector UI and logic.
# Restricts choices to columns stored in the selector config.
.gene_group_dge_defaults <- function(cntr_ds, param_state, selector_config) {
  selected_contrasts <- param_state$input_values$dge_contrasts
  p_col <- selector_config$dge$pval_col
  lfc_col <- selector_config$dge$lfc_col
  fdr_col <- selector_config$dge$fdr_col
  cntr_names <- names(cntr_ds)
  sel_cntr <- intersect(selected_contrasts, cntr_names)
  if(length(sel_cntr) == 0L) {
    sel_cntr <- cntr_names[1]
  }

  num_cols <- .gene_group_dge_numeric_cols(cntr_ds, sel_cntr)
  resolved_p_col <- NULL
  resolved_lfc_col <- NULL
  resolved_fdr_col <- NULL

  if(isTruthy(p_col) && p_col %in% num_cols) {
    resolved_p_col <- p_col
  }
  if(isTruthy(lfc_col) && lfc_col %in% num_cols) {
    resolved_lfc_col <- lfc_col
  }
  if(isTruthy(fdr_col) && fdr_col %in% num_cols) {
    resolved_fdr_col <- fdr_col
  }

  list(
    contrasts=cntr_names,
    selected_contrasts=sel_cntr,
    p_col=resolved_p_col,
    lfc_col=resolved_lfc_col,
    fdr_col=resolved_fdr_col
  )
}

# Build the controls for DGE-based gene selection.
# Renders thresholds, ranking, and top-k inputs for DGE mode.
.gene_group_dge_controls_ui <- function(ns, dge_defaults, top_mode, top_value,
                                        rank_mode, p_thr, lfc_thr, require_fdr) {
  top_value <- suppressWarnings(as.numeric(top_value)[1])
  if(!is.finite(top_value)) {
    top_value <- if(top_mode == "pct") 10 else 100
  }

  p_thr <- suppressWarnings(as.numeric(p_thr)[1])
  if(!is.finite(p_thr)) {
    p_thr <- .05
  }

  lfc_thr <- suppressWarnings(as.numeric(lfc_thr)[1])
  if(!is.finite(lfc_thr)) {
    lfc_thr <- 0
  }

  top_value_ui <- NULL
  if(top_mode == "n") {
    top_value_ui <- numericInput(ns("dge_top_value"), "Top N genes", value=max(1, floor(top_value)), min=1, step=1)
  } else if(top_mode == "pct") {
    top_value_ui <- numericInput(
      ns("dge_top_value"),
      "Top percentage",
      value=max(0, min(100, top_value)),
      min=0,
      max=100,
      step=1
    )
  }

  shiny::tags$div(
    fluidRow(
      column(
        12,
        selectizeInput(ns("dge_contrasts"), "Contrast", choices=dge_defaults$contrasts, selected=dge_defaults$selected_contrasts, multiple=FALSE),
      )
    ),
    fluidRow(
      column(
        6,
        numericInput(ns("dge_pval_thr"), "P-value threshold", value=max(0, min(1, p_thr)), min=0, max=1, step=.001),
        selectInput(ns("dge_rank_mode"), "Ranking priority", choices=c("Lowest p-value first"="p_first", "Highest absolute logFC first"="lfc_first"), selected=rank_mode),
        checkboxInput(ns("dge_require_fdr"), "Exclude genes with missing adjusted p-value (FDR)", value=isTRUE(require_fdr))
      ),
      column(
        6,
        numericInput(ns("dge_lfc_thr"), "Absolute logFC threshold", value=max(0, lfc_thr), min=0, step=.1),
        selectInput(ns("dge_top_mode"), "Top filter", choices=c("All passing filters"="all", "Top N"="n", "Top %"="pct"), selected=top_mode),
        top_value_ui
      )
    )
  )
}

# Build the mode-specific control panel for the selector.
# Dispatches the UI rendering based on the current parameter state.
.gene_group_modus_controls_ui <- function(ns, ann_df=NULL, expr_ds=NULL, cntr_ds=NULL,
                                          param_state, selector_config) {
  mode <- param_state$mode
  input_values <- param_state$input_values

  if(mode == "by_name") {
    ann_cols <- colnames(ann_df)
    selected_col <- input_values$name_id_col %||% "PrimaryID"
    if(!selected_col %in% ann_cols && "PrimaryID" %in% ann_cols) {
      selected_col <- "PrimaryID"
    }
    return(.gene_group_name_controls_ui(
      ns,
      ann_cols,
      selected_col,
      name_list=input_values$name_list
    ))
  }

  if(mode == "by_expression") {
    if(is.null(expr_ds)) {
      return(p("Expression data not available."))
    }
    resolved <- .gene_group_resolve_mode_inputs(param_state, selector_config)
    return(.gene_group_expression_controls_ui(ns, resolved$metric, resolved$top_mode, resolved$top_value))
  }

  if(mode == "by_dge") {
    if(is.null(cntr_ds)) {
      return(p("Contrast data not available."))
    }
    cntr_names <- names(cntr_ds)
    if(length(cntr_names) == 0L) {
      return(p("No contrasts available for this dataset."))
    }

    dge_defaults <- .gene_group_dge_defaults(cntr_ds, param_state, selector_config)
    resolved <- .gene_group_resolve_mode_inputs(param_state, selector_config)
    return(.gene_group_dge_controls_ui(
      ns=ns,
      dge_defaults=dge_defaults,
      top_mode=resolved$top_mode,
      top_value=resolved$top_value,
      rank_mode=resolved$rank_mode,
      p_thr=resolved$p_thr,
      lfc_thr=resolved$lfc_thr,
      require_fdr=resolved$require_fdr
    ))
  }

  p("No selection mode available.")
}







## ---------------------------------------------
## Selection Logic
## ---------------------------------------------

# Create an empty annotation data frame with preserved columns.
# Keeps "no selection" results compatible with downstream code.
.gene_group_empty_annot <- function(annot_df) {
  annot_df[0, , drop=FALSE]
}

# Match query IDs against one annotation column.
# Preserves query order and drops IDs that are not present.
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

# Select annotation rows from typed and uploaded gene IDs.
# Combines both sources into one ordered list before matching.
.gene_group_select_by_name <- function(ann_df, id_col, name_list, name_file) {
  name_ids <- .gene_group_parse_gene_input(name_list)
  file_ids <- .gene_group_read_gene_file(name_file)
  query_ids <- unique(c(name_ids, file_ids))
  .gene_group_log("select_by_name using column='", id_col, "'; query_ids n=",
                  as.character(length(query_ids)), ".")
  .gene_group_match_annot_rows(ann_df, id_col, query_ids)
}

# Select annotation rows by expression-derived ranking.
# Uses the chosen metric and top-k rule to keep the highest-ranked genes.
.gene_group_select_by_expression <- function(ann_df, expr_ds, primary_id, metric="variance", top_mode="n", top_value=100) {
  empty_df <- .gene_group_empty_annot(ann_df)
  if(is.null(expr_ds)) {
    .gene_group_log("select_by_expression skipped: expr_ds is NULL.")
    return(empty_df)
  }

  expr_ds <- as.matrix(expr_ds)
  if(nrow(expr_ds) < 1L || is.null(rownames(expr_ds))) {
    .gene_group_log("select_by_expression skipped: expression has no rows or rownames.")
    return(empty_df)
  }

  score <- if(metric == "mean") {
    rowMeans(expr_ds, na.rm=TRUE)
  } else {
    apply(expr_ds, 1, stats::var, na.rm=TRUE)
  }
  score <- score[!is.na(score)]
  if(length(score) == 0L) {
    .gene_group_log("select_by_expression skipped: no finite scores.")
    return(empty_df)
  }

  k <- .gene_group_top_k(length(score), top_mode, top_value)
  if(k < 1L) {
    .gene_group_log("select_by_expression skipped: top-k resolved to 0.")
    return(empty_df)
  }

  ids <- names(sort(score, decreasing=TRUE))[seq_len(k)]
  .gene_group_log("select_by_expression metric='", metric, "', top_mode='", top_mode,
                  "', top_value=", as.character(top_value), ", selected n=",
                  as.character(length(ids)), ".")
  .gene_group_match_annot_rows(ann_df, primary_id, ids)
}

# Build a compact per-gene DGE summary across selected contrasts.
# Keeps the values needed for later ranking and threshold filtering.
.gene_group_collect_dge_rows <- function(cntr_ds, contrasts, primary_id, p_col, lfc_col, fdr_col, require_fdr, p_thr, lfc_thr) {
  dge_rows <- lapply(contrasts, function(nm) {
    df <- cntr_ds[[nm]]
    required_cols <- c(primary_id, p_col, lfc_col)
    if(require_fdr) {
      required_cols <- c(required_cols, fdr_col)
    }
    if(!all(required_cols %in% colnames(df))) {
      return(NULL)
    }

    ids <- as.character(df[[primary_id]])
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

# Rank genes by their strongest DGE evidence across contrasts.
# Supports ranking by p-value first or by absolute log fold change first.
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

# Select annotation rows from DGE results after filtering and ranking.
# Applies thresholding first and then an optional top-k restriction.
.gene_group_select_by_dge <- function(ann_df, cntr_ds, primary_id,
                                      contrasts, p_col, lfc_col, fdr_col,
                                      require_fdr=FALSE, p_thr=.05, lfc_thr=0,
                                      rank_mode="p_first", top_mode="all", top_value=100) {
  empty_df <- .gene_group_empty_annot(ann_df)
  if(is.null(cntr_ds) || length(cntr_ds) < 1L) {
    .gene_group_log("select_by_dge skipped: cntr_ds is NULL/empty.")
    return(empty_df)
  }

  contrasts <- intersect(contrasts, names(cntr_ds))
  if(length(contrasts) == 0L) {
    .gene_group_log("select_by_dge skipped: no matching contrasts selected.")
    return(empty_df)
  }

  if(!isTruthy(p_col) || !isTruthy(lfc_col) || (isTRUE(require_fdr) && !isTruthy(fdr_col))) {
    .gene_group_log("select_by_dge skipped: invalid DGE column configuration ",
                    "(p_col='", as.character(p_col), "', lfc_col='", as.character(lfc_col),
                    "', fdr_col='", as.character(fdr_col), "', require_fdr=",
                    as.character(isTRUE(require_fdr)), ").")
    return(empty_df)
  }

  p_thr <- suppressWarnings(as.numeric(p_thr)[1])
  lfc_thr <- suppressWarnings(as.numeric(lfc_thr)[1])
  if(is.na(p_thr) || is.na(lfc_thr)) {
    .gene_group_log("select_by_dge skipped: invalid numeric thresholds.")
    return(empty_df)
  }

  dge_rows <- .gene_group_collect_dge_rows(
    cntr_ds=cntr_ds,
    contrasts=contrasts,
    primary_id=primary_id,
    p_col=p_col,
    lfc_col=lfc_col,
    fdr_col=fdr_col,
    require_fdr=isTRUE(require_fdr),
    p_thr=p_thr,
    lfc_thr=lfc_thr
  )
  if(is.null(dge_rows) || nrow(dge_rows) == 0L) {
    .gene_group_log("select_by_dge: no rows passed filters.")
    return(empty_df)
  }

  rank_df <- .gene_group_rank_dge_rows(dge_rows, rank_mode=rank_mode)
  if(top_mode %in% c("n", "pct")) {
    k <- .gene_group_top_k(nrow(rank_df), top_mode, top_value)
    if(k < 1L) {
      .gene_group_log("select_by_dge: top filter reduced selection to 0.")
      return(empty_df)
    }
    rank_df <- rank_df[seq_len(k), , drop=FALSE]
  }

  .gene_group_log(
    "select_by_dge contrasts={", paste(contrasts, collapse=","), "}, rank_mode='", rank_mode,
    "', top_mode='", top_mode, "', top_value=", as.character(top_value),
    ", thresholds(p=", as.character(p_thr), ", lfc=", as.character(lfc_thr),
    "), selected n=", as.character(nrow(rank_df)), "."
  )
  .gene_group_match_annot_rows(ann_df, primary_id, rank_df$id)
}

# Dispatch the active mode to the matching selection helper.
# Returns one annotation subset from the current state and selector config.
.gene_group_selected_annotation <- function(ann_df, expr_ds=NULL, cntr_ds=NULL,
                                            param_state, selector_config) {
  mode <- param_state$mode
  input_values <- param_state$input_values
  primary_id <- selector_config$primary_id

  .gene_group_log("selected_annotation dispatch mode='", mode, "'.")
  if(mode == "by_name") {
    name_id_col <- input_values$name_id_col %||% primary_id
    name_id_col <- as.character(name_id_col)[1]
    if(!name_id_col %in% colnames(ann_df)) {
      name_id_col <- primary_id
    }
    return(.gene_group_select_by_name(ann_df, name_id_col, input_values$name_list, input_values$name_file))
  }

  if(mode == "by_expression") {
    resolved <- .gene_group_resolve_mode_inputs(param_state, selector_config)
    return(.gene_group_select_by_expression(ann_df, expr_ds, primary_id, resolved$metric, resolved$top_mode, resolved$top_value))
  }

  if(mode == "by_dge") {
    resolved <- .gene_group_resolve_mode_inputs(param_state, selector_config)
    dge_defaults <- .gene_group_dge_defaults(cntr_ds, param_state, selector_config)
    return(.gene_group_select_by_dge(
      ann_df=ann_df,
      cntr_ds=cntr_ds,
      primary_id=primary_id,
      contrasts=input_values$dge_contrasts,
      p_col=dge_defaults$p_col,
      lfc_col=dge_defaults$lfc_col,
      fdr_col=dge_defaults$fdr_col,
      require_fdr=resolved$require_fdr,
      p_thr=resolved$p_thr,
      lfc_thr=resolved$lfc_thr,
      rank_mode=resolved$rank_mode,
      top_mode=resolved$top_mode,
      top_value=resolved$top_value
    ))
  }

  .gene_group_empty_annot(ann_df)
}







## ---------------------------------------------
## Server State Helpers
## ---------------------------------------------

# Capture the current selector inputs as a plain list.
# Keeps server code from rebuilding the same input bundle inline.
.gene_group_input_values <- function(input) {
  list(
    name_id_col=input$name_id_col,
    name_list=input$name_list,
    name_file=input$name_file,
    expr_metric=input$expr_metric,
    expr_top_mode=input$expr_top_mode,
    expr_top_value=input$expr_top_value,
    dge_contrasts=input$dge_contrasts,
    dge_require_fdr=input$dge_require_fdr,
    dge_pval_thr=input$dge_pval_thr,
    dge_lfc_thr=input$dge_lfc_thr,
    dge_top_mode=input$dge_top_mode,
    dge_top_value=input$dge_top_value,
    dge_rank_mode=input$dge_rank_mode
  )
}

# Build the static selector configuration shared across helper calls.
# Groups defaults, ID settings, and DGE column mappings in one object.
.gene_group_selector_config <- function(primary_id, defaults, dge_pval_col, dge_lfc_col, dge_fdr_col) {
  list(
    primary_id=primary_id,
    defaults=defaults,
    dge=list(
      pval_col=dge_pval_col,
      lfc_col=dge_lfc_col,
      fdr_col=dge_fdr_col
    )
  )
}

# Override name-based inputs with an externally supplied ID selection.
# Clears the uploaded file so the external IDs fully define the query.
.gene_group_apply_external_name_state <- function(input_values, external_name_state) {
  if(is.null(external_name_state)) {
    return(input_values)
  }

  input_values$name_id_col <- external_name_state$name_id_col
  input_values$name_list <- external_name_state$name_list
  input_values$name_file <- NULL
  input_values
}

# Convert selected IDs into the stored external by-name override state.
# Uses the primary ID column and newline-separated IDs for the text area.
.gene_group_external_name_state <- function(primary_id, ids) {
  list(
    name_id_col=primary_id,
    name_list=paste(as.character(ids), collapse="\n")
  )
}

# Build the dynamic selector parameter state from reactive inputs.
# Combines dataset, active mode, and the effective control values.
.gene_group_param_state <- function(dataset, mode, input_values) {
  list(
    dataset=dataset,
    mode=mode,
    input_values=input_values
  )
}







## ---------------------------------------------
## Public UI
## ---------------------------------------------

#' @rdname geneGroupSelectorServer
#' @export
# Build controls-only UI for selecting gene groups.
# Intended to be placed inside a parent sidebarPanel.
geneGroupSelectorUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(6, uiOutput(ns("dataset_ui"))),
      column(6, selectizeInput(ns("modus"), "Modus", choices=character(0), selected=NULL))
    ),
    shiny::tags$hr(),
    uiOutput(ns("modus_controls")),
    br(),
    textOutput(ns("selected_count")),
    downloadButton(ns("save_genes"), "Export selected PrimaryIDs")
  )
}







## ---------------------------------------------
## Public Server
## ---------------------------------------------

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
#' @param dge_pval_col optional p-value column name in contrast tables for DGE mode.
#' @param dge_lfc_col optional log fold-change column name in contrast tables for
#'   DGE mode.
#' @param dge_fdr_col optional adjusted p-value (FDR) column name in contrast tables
#'   for DGE mode.
#' @param mode_order order of available selection modes. Values must come from
#'   `c("by_expression", "by_dge", "by_name")`. Unavailable modes are skipped.
#'   Default is expression -> DGE -> by name.
#' @param defaults optional list of default selector parameters. Supported keys:
#'   `expr_metric`, `expr_top_mode`, `expr_top_value`, `dge_top_mode`,
#'   `dge_top_value`, `dge_rank_mode`, `dge_pval_thr`, `dge_lfc_thr`,
#'   and `dge_require_fdr`. Invalid or unsupported values raise an error.
#' @param selected_ids optional `reactiveVal()` that will be updated
#'   with the currently selected PrimaryIDs (column `primary_id`) as a
#'   character vector. When it changes externally, the selector switches to
#'   `by_name` mode and mirrors those IDs in the text field.
#' @param dataset optional `reactiveVal()` that stores the currently selected
#'   dataset name. If `NULL`, it is initialized internally with `reactiveVal()`.
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
                                    mode_order=c("by_expression", "by_dge", "by_name"),
                                    primary_id="PrimaryID",
                                    dge_pval_col=NULL,
                                    dge_lfc_col=NULL,
                                    dge_fdr_col=NULL,
                                    defaults=NULL,
                                    selected_ids=NULL,
                                    dataset=NULL) {
  .gene_group_log("geneGroupSelectorServer init id='", id, "'.")

  # after this, annot is a list with data frames (different datasets)
  # or a single one called "default"
  annot <- .gene_group_normalize_annot(annot)

  # sanity checks; stop with informative errors if the input is not as expected
  .gene_group_check_primary_id(annot, primary_id)
  .gene_group_check_selected_ids(selected_ids)
  .gene_group_check_dataset(dataset)

  defaults <- .gene_group_normalize_defaults(defaults)
  mode_order <- .gene_group_normalize_mode_order(mode_order)
  dge_pval_col <- trimws(as.character(dge_pval_col %||% "")[1])
  dge_lfc_col <- trimws(as.character(dge_lfc_col %||% "")[1])
  dge_fdr_col <- trimws(as.character(dge_fdr_col %||% "")[1])

  if(!nzchar(dge_pval_col)) {
    dge_pval_col <- NULL
  }

  if(!nzchar(dge_lfc_col)) {
    dge_lfc_col <- NULL
  }

  if(!nzchar(dge_fdr_col)) {
    dge_fdr_col <- NULL
  }

  datasets <- names(annot)
  dataset <- dataset %||% reactiveVal()
  exprs <- .gene_group_normalize_exprs(exprs, datasets)
  cntr <- .gene_group_normalize_cntr(cntr, datasets)
  selector_config <- .gene_group_selector_config(
    primary_id=primary_id,
    defaults=defaults,
    dge_pval_col=dge_pval_col,
    dge_lfc_col=dge_lfc_col,
    dge_fdr_col=dge_fdr_col
  )
  modes <- .gene_group_modes(!is.null(exprs), !is.null(cntr), mode_order=mode_order)

  .gene_group_log("server datasets={", paste(datasets, collapse=","), "}; modes={",
                  paste(unname(modes), collapse=","), "}.")

  ## -----------------------------------------------------
  ## Reactive and observer definitions
  ## -----------------------------------------------------
  moduleServer(id, function(input, output, session) {
    .gene_group_log("moduleServer started for id='", id, "'.")
    default_mode <- unname(modes[1])

    # state vars.
    # external_name_state reacts to external changes to selected_ids
    external_name_state <- reactiveVal(NULL)
    pushing_selected_ids <- reactiveVal(FALSE)
    last_published_ids <- reactiveVal(character(0))

    # dynamically render the dataset selection input based on the provided annotation datasets
    # b/c dataset() is a reactiveVal.
    output$dataset_ui <- renderUI({
      .gene_group_log("rendering dataset UI with dataset()='", as.character(dataset()), "'.")
      .gene_group_dataset_input_ui(session$ns, datasets, selected=dataset())
    })

    .gene_group_log("initializing modus choices={", paste(unname(modes), collapse=","), "}.")
    updateSelectizeInput(
      session=session,
      inputId="modus",
      choices=modes,
      selected=default_mode,
      server=TRUE
    )

    # Keep the shared dataset reactiveVal on a valid dataset name.
    # Repairs invalid external or initial values before the rest of the module uses them.
    observe({
      .ds <- isolate({ dataset() })
      if(!isTruthy(.ds) || !.ds %in% datasets) {
        .gene_group_log("dataset reactive invalid ('", as.character(.ds),
                        "'); resetting to '", datasets[1], "'.")
        dataset(datasets[1])
      }
    })

    # Mirror changes from the dataset input widget into the shared dataset reactiveVal.
    # Normalizes invalid UI values to the first available dataset before publishing them.
    observeEvent(input$dataset, {
      .ds <- input$dataset
      if(!isTruthy(.ds) || !.ds %in% datasets) {
        .ds <- datasets[1]
      }
      .gene_group_log("dataset input changed to '", as.character(input$dataset),
                      "'; effective dataset='", .ds, "'.")
      isolate({ 
        if(!identical(dataset(), .ds)) {
          dataset(.ds)
        }
      })
    }, ignoreInit=FALSE)

    # Derive the effective selector mode from the mode widget, with a by-name fallback.
    # External selected IDs only influence mode while the widget is transiently blank.
    current_mode <- reactive({
      mode_input <- trimws(as.character(input$modus %||% "")[1])
      if(!nzchar(mode_input)) {
        if(!is.null(external_name_state())) {
          return("by_name")
        }
        return(default_mode)
      }
      if(!mode_input %in% unname(modes)) {
        return(default_mode)
      }
      mode_input
    })

    # Build the current dynamic parameter state used by UI rendering and selection logic.
    # Combines dataset, effective mode, and the effective control input values in one object.
    param_state <- reactive({
      .gene_group_param_state(
        dataset=dataset(),
        mode=current_mode(),
        input_values=.gene_group_apply_external_name_state(
          .gene_group_input_values(input),
          external_name_state()
        )
      )
    })

    # prepare the UI controls for the active mode
    output$modus_controls <- renderUI({
      ds <- dataset()
      mode <- current_mode()
      input_values <- isolate(.gene_group_apply_external_name_state(
        .gene_group_input_values(input),
        external_name_state()
      ))
      state <- .gene_group_param_state(
        dataset=ds,
        mode=mode,
        input_values=input_values
      )
      .gene_group_log("rendering modus controls for dataset='", ds, "', mode='", mode, "'.")

      .gene_group_modus_controls_ui(
        ns=session$ns,
        ann_df=annot[[ds]],
        expr_ds=if(is.null(exprs)) NULL else exprs[[ds]],
        cntr_ds=if(is.null(cntr)) NULL else cntr[[ds]],
        param_state=state,
        selector_config=selector_config
      )
    })

    # Compute the selected annotation rows for the current dataset and mode.
    # Dispatches to the appropriate by-name, expression, or DGE selection helper.
    selected_annotation <- reactive({
      state <- param_state()
      ds <- state$dataset
      out <- .gene_group_selected_annotation(
        ann_df=annot[[ds]],
        expr_ds=if(is.null(exprs)) NULL else exprs[[ds]],
        cntr_ds=if(is.null(cntr)) NULL else cntr[[ds]],
        param_state=state,
        selector_config=selector_config
      )
      .gene_group_log("selected_annotation dataset='", ds, "', mode='", state$mode,
                      "' -> nrows=", as.character(nrow(out)), ".")
      out
    })

    # Extract the unique primary IDs from the selected annotation rows.
    # This is the canonical gene list shared with downstream outputs and modules.
    selected_primary_ids <- reactive({
      sel <- selected_annotation()
      if(nrow(sel) == 0L) {
        .gene_group_log("selected_primary_ids: no rows selected.")
        return(character(0))
      }

      ids <- as.character(sel[[primary_id]])
      ids <- ids[!is.na(ids) & ids != ""]
      ids <- unique(ids)
      .gene_group_log("selected_primary_ids n=", as.character(length(ids)), ".")
      ids
    })

    # Drop the external by-name override when the user intentionally changes mode.
    observeEvent(input$modus, {
      mode_input <- trimws(as.character(input$modus %||% "")[1])
      if(nzchar(mode_input) && !identical(mode_input, "by_name")) {
        external_name_state(NULL)
      }
    }, ignoreInit=TRUE)

    # Drop the external by-name override when the ID-column selector is changed manually.
    # This returns control of by-name settings to the live UI once the user edits them.
    observeEvent(input$name_id_col, {
      state <- external_name_state()
      if(!is.null(state) &&
         !identical(input$name_id_col %||% "", state$name_id_col %||% "")) {
        external_name_state(NULL)
      }
    }, ignoreInit=TRUE)

    # Drop the external by-name override when the gene text area is edited manually.
    # This lets free-text editing take precedence over an imported external selection.
    observeEvent(input$name_list, {
      state <- external_name_state()
      if(!is.null(state) &&
         !identical(input$name_list %||% "", state$name_list %||% "")) {
        external_name_state(NULL)
      }
    }, ignoreInit=TRUE)

    # if selected_ids is provided, sync it with the selector's current selection
    if(!is.null(selected_ids)) {
      # Import externally supplied selected IDs into the selector as a by-name override.
      # Ignores the module's own last publication so bidirectional syncing does not loop.
      observeEvent(selected_ids(), {
        ids <- as.character(selected_ids() %||% character(0))

        if(isTRUE(pushing_selected_ids()) &&
           identical(ids, last_published_ids())) {
          pushing_selected_ids(FALSE)
          return()
        }

        pushing_selected_ids(FALSE)
        external_name_state(.gene_group_external_name_state(selector_config$primary_id, ids))

        if(!identical(trimws(as.character(isolate(input$modus %||% default_mode))[1]), "by_name")) {
          shiny::updateSelectizeInput(
            session=session,
            inputId="modus",
            choices=modes,
            selected="by_name",
            server=TRUE
          )
        }
      }, ignoreInit=TRUE)
    }

    # Publish the current selected primary IDs to the external selected_ids reactiveVal.
    # Skips unchanged values and marks self-originated updates to prevent feedback loops.
    observe({
      if(!is.null(selected_ids)) {
        ids <- selected_primary_ids()
        .gene_group_log("updating external selected_ids reactiveVal with n=",
                        as.character(length(ids)), ".")
        if(!identical(isolate(selected_ids()), ids)) {
          pushing_selected_ids(TRUE)
          last_published_ids(ids)
          selected_ids(ids)
        }
      }
    })

    output$selected_count <- renderText({
      sprintf("Selected genes: %d", length(selected_primary_ids()))
    })

    output$save_genes <- downloadHandler(
      filename = function() {
        ds <- dataset()
        md <- current_mode()
        .gene_group_log("preparing gene export filename for dataset='", ds, "', mode='", md, "'.")
        sprintf(
          "selected_genes_%s_%s.txt",
          sanitize_filename(ds, "dataset"),
          sanitize_filename(md, "modus")
        )
      },
      content = function(file) {
        .gene_group_log("writing selected genes to file='", file, "'; n=",
                        as.character(length(selected_primary_ids())), ".")
        writeLines(selected_primary_ids(), con=file)
      }
    )

    return(list(
      genes=selected_primary_ids,
      annotation=selected_annotation,
      dataset=dataset,
      modus=current_mode
    ))
  })
}
