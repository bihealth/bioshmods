# Check whether x is a non-empty list whose elements are all data frames.
# Used to classify data payloads for export behavior.
.file_export_is_df_list <- function(x) {
  is.list(x) && length(x) > 0L && all(vapply(x, is.data.frame, logical(1)))
}

# Check whether x is a list of non-empty data-frame lists.
# This identifies the "list of lists of data frames" case.
.file_export_is_nested_df_list <- function(x) {
  is.list(x) && length(x) > 0L && all(vapply(x, .file_export_is_df_list, logical(1)))
}

# Classify the data payload into single_df, df_list, nested_df_list, or other.
# The class decides which UI controls and export format are used.
.file_export_item_structure <- function(data_obj) {
  if(is.data.frame(data_obj)) {
    return("single_df")
  }
  if(.file_export_is_df_list(data_obj)) {
    return("df_list")
  }
  if(.file_export_is_nested_df_list(data_obj)) {
    return("nested_df_list")
  }
  "other"
}

# Sanitize a text label into a filesystem-safe base filename.
# Falls back to default when the input is missing or empty.
.file_export_safe_name <- function(x, default="file") {
  x <- trimws(as.character(x)[1])

  if(is.na(x) || x == "") {
    return(default)
  }

  x <- gsub("[^0-9A-Za-z_.-]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_+|_+$", "", x)

  if(x == "") {
    x <- default
  }

  x
}

# Normalize candidate Excel sheet names to valid unique names.
# Enforces Excel constraints such as max length and illegal characters.
.file_export_sheet_names <- function(x, prefix="Sheet") {
  x <- as.character(x)
  n <- length(x)
  if(n == 0L) {
    stop("Cannot create an Excel file from an empty list.")
  }

  x[is.na(x)] <- ""
  x <- trimws(x)

  defaults <- paste0(prefix, seq_len(n))
  x[x == ""] <- defaults[x == ""]
  x <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", x)
  x <- gsub("^'+|'+$", "", x)
  x[nchar(x) == 0L] <- defaults[nchar(x) == 0L]
  x <- substr(x, 1L, 31L)

  out <- character(n)
  used <- character(0)
  for(i in seq_len(n)) {
    base <- x[[i]]
    if(!nzchar(base)) {
      base <- defaults[[i]]
    }
    candidate <- base
    j <- 1L
    while(candidate %in% used) {
      suffix <- paste0("_", j)
      max_base <- 31L - nchar(suffix)
      if(max_base < 1L) {
        max_base <- 1L
      }
      candidate <- paste0(substr(base, 1L, max_base), suffix)
      j <- j + 1L
    }
    out[[i]] <- candidate
    used <- c(used, candidate)
  }

  out
}

# Convert a single data frame or list of data frames into named sheets.
# Output is a named list ready for writexl::write_xlsx().
.file_export_as_sheet_list <- function(data_obj, fallback_sheet="Sheet1") {

  if(is.data.frame(data_obj)) {
    sheet_name <- .file_export_sheet_names(fallback_sheet, prefix="Sheet")
    return(setNames(list(data_obj), sheet_name))
  }

  if(.file_export_is_df_list(data_obj)) {
    nm <- names(data_obj)
    if(is.null(nm)) {
      nm <- rep("", length(data_obj))
    }
    nm <- as.character(nm)
    nm[is.na(nm)] <- ""
    nm <- trimws(nm)

    empty <- nm == ""
    if(any(empty)) {
      base <- .file_export_safe_name(fallback_sheet, default="Sheet")
      nm[empty] <- paste0(base, "_", seq_len(length(data_obj))[empty])
    }

    nm <- .file_export_sheet_names(nm, prefix=.file_export_safe_name(fallback_sheet, "Sheet"))
    return(setNames(data_obj, nm))
  }

  stop("Expected a data frame or a list of data frames.")
}

# Build stable, non-empty labels for nested list selection controls.
# Names are made unique to avoid ambiguous selectInput choices.
.file_export_choice_names <- function(data_obj) {
  nm <- names(data_obj)
  if(is.null(nm)) {
    nm <- rep("", length(data_obj))
  }
  nm <- as.character(nm)
  nm[is.na(nm)] <- ""
  nm <- trimws(nm)
  default_names <- paste0("Selection ", seq_len(length(data_obj)))
  nm[nm == ""] <- default_names[nm == ""]
  make.unique(nm, sep="_")
}

# Convert a help object into UI markup for display in one export row.
# Accepts shiny tags directly or plain text converted to HTML.
.file_export_help_ui <- function(help_obj) {
  if(is.null(help_obj)) {
    return(NULL)
  }

  if(inherits(help_obj, "shiny.tag") || inherits(help_obj, "shiny.tag.list")) {
    return(help_obj)
  }

  txt <- paste(as.character(help_obj), collapse="<br/>")
  if(!nzchar(trimws(txt))) {
    return(NULL)
  }

  p(HTML(txt), class="text-muted", style="margin:0;")
}

## Validate and normalize a single user-provided object definition.
.file_export_normalize_single_object <- function(obj, i) {

  if(!is.list(obj)) {
    stop(sprintf("objects[[%d]] must be a list.", i))
  }

  data_obj <- obj[["data"]]
  if(is.null(data_obj)) {
    stop(sprintf("objects[[%d]] must contain `data`.", i))
  }

  title       <- obj[["title"]] %||% sprintf("Object %d", i)
  description <- obj[["description"]] %||% ""
  data_type   <- obj[["data_type"]] %||% "rds"
  data_type   <- match.arg(data_type, c("dataframes", "rds"))

  structure_type <- .file_export_item_structure(data_obj)

  if(data_type == "dataframes" && structure_type == "other") {
    cat("Invalid data structure:\n-----------------------\n")
    print(str(data_obj))
    stop(sprintf(
      "objects[[%d]] with data_type='dataframes' must be a data frame, list of data frames, or list of lists of data frames.",
      i
    ))
  }

  if(data_type == "rds") {
    structure_type <- "rds"
  }

  list(
    data=data_obj,
    title=title,
    description=description,
    data_type=data_type,
    structure_type=structure_type,
    help=obj[["help"]]
  )
}


# Validate and normalize user-provided object definitions.
# Returns a standardized list used by both UI and server logic.
.file_export_normalize_objects <- function(objects) {
  if(!is.list(objects) || length(objects) < 1L) {
    stop("`objects` must be a non-empty list.")
  }

  lapply(seq_along(objects), function(i) {
    obj <- objects[[i]]
    .file_export_normalize_single_object(obj, i)
  })
}

# Build one UI row for one export object definition.
# Row includes title/description/help plus context-specific controls.
.file_export_line_ui <- function(id, spec, idx) {
  save_id <- sprintf("save_%d", idx)
  select_id <- sprintf("select_%d", idx)
  save_all_id <- sprintf("save_all_%d", idx)

  info_ui <- tagList(
    p(strong(spec$title), " ", spec$description, style="margin:0 0 4px 0;")
  )

  help_ui <- .file_export_help_ui(spec$help)
  if(!is.null(help_ui)) {
    info_ui <- tagList(info_ui, help_ui)
  }

  controls_ui <- NULL
  if(spec$data_type == "dataframes" && spec$structure_type == "nested_df_list") {
    choices <- .file_export_choice_names(spec$data)
    controls_ui <- tagList(
      selectInput(NS(id, select_id), label=NULL, choices=choices, selected=choices[1], width="100%"),
      downloadButton(NS(id, save_id), "save as"),
      HTML("&nbsp;"),
      downloadButton(NS(id, save_all_id), "save all")
    )
  } else {
    controls_ui <- downloadButton(NS(id, save_id), "save as")
  }

  fluidRow(
    column(8, info_ui),
    column(4, controls_ui),
    style="margin-bottom: 12px;"
  )
}

# Write a data frame payload to an Excel file (one or many sheets).
# Sheet names are normalized to ensure valid workbook output.
.file_export_write_xlsx <- function(data_obj, file, fallback_sheet="Sheet1") {
  sheets <- .file_export_as_sheet_list(data_obj, fallback_sheet=fallback_sheet)
  writexl::write_xlsx(sheets, path=file)
}

# Write all nested data-frame groups into separate xlsx files and zip them.
# Each nested entry becomes one workbook inside the final archive.
.file_export_write_zip_all <- function(nested_item, file) {
  choices <- .file_export_choice_names(nested_item)
  # Create a temporary staging directory for xlsx files before zipping.
  tmp_dir <- tempfile("file_export_zip_")
  dir.create(tmp_dir, recursive=TRUE)
  on.exit(unlink(tmp_dir, recursive=TRUE), add=TRUE)

  base_names <- vapply(seq_along(nested_item), function(i) {
    .file_export_safe_name(choices[[i]], default=sprintf("selection_%d", i))
  }, character(1))
  base_names <- make.unique(base_names, sep="_")
  xlsx_names <- paste0(base_names, ".xlsx")

  # Materialize one workbook per nested group using normalized filenames.
  for(i in seq_along(nested_item)) {
    .file_export_write_xlsx(
      nested_item[[i]],
      file.path(tmp_dir, xlsx_names[[i]]),
      fallback_sheet=choices[[i]]
    )
  }

  zip::zipr(zipfile=file, files=xlsx_names, root=tmp_dir)
}

#' @rdname fileExportServer
#' @export
# Build the static UI rows for all export object definitions.
# The output is a list of row components, one per object.
fileExportUI <- function(id, objects) {
  # Normalize once so rendering and server logic share one schema.
  specs <- .file_export_normalize_objects(objects)
  tagList(lapply(seq_along(specs), function(i) .file_export_line_ui(id, specs[[i]], i)))
}

#' Shiny module for exporting tabular and generic objects
#'
#' Shiny module for exporting a list of objects as Excel (`.xlsx`), RDS
#' (`.rds`), or ZIP archives (`.zip`) depending on each object's
#' specification.
#'
#' Each entry in `objects` must be a list with:
#' - `data`: data payload
#' - `title`: short title shown in the UI
#' - `description`: descriptive text shown after the title
#' - `data_type`: either `"dataframes"` or `"rds"`
#' - `help`: optional help text/tag shown below description
#'
#' For `data_type = "dataframes"`:
#' - `data` can be a single data frame, a list of data frames, or a list of
#'   lists of data frames.
#' - single/list: one `"save as"` button exports one Excel file (one sheet
#'   per data frame).
#' - list of lists: UI includes a selection dropdown, `"save as"` for the
#'   selected list, and `"save all"` for a ZIP of Excel files.
#'
#' For `data_type = "rds"`:
#' - `"save as"` exports the object as RDS.
#'
#' @param id module identifier (same as the one passed to [fileExportUI()])
#' @param objects list of object definitions as described above
#' @return `fileExportUI` returns UI elements. `fileExportServer` returns
#'   nothing useful.
#' @examples
#' data(C19)
#' if(interactive()) {
#'   export_objects <- list(
#'     list(
#'       data = C19$covariates,
#'       title = "Covariates",
#'       description = "Sample metadata table",
#'       data_type = "dataframes",
#'       help = "Saved as a single-sheet Excel file"
#'     ),
#'     list(
#'       data = list(
#'         dataset_a = C19$contrasts[1:2],
#'         dataset_b = C19$contrasts[1:2]
#'       ),
#'       title = "Contrasts",
#'       description = "Choose one contrast set or export all",
#'       data_type = "dataframes",
#'       help = "Save one set as xlsx or all sets as a zip archive"
#'     ),
#'     list(
#'       data = C19,
#'       title = "Full C19 object",
#'       description = "Complete data object",
#'       data_type = "rds",
#'       help = "Saved as RDS"
#'     )
#'   )
#'
#'   ui <- fluidPage(fileExportUI("fexp", export_objects))
#'   server <- function(input, output, session) {
#'     fileExportServer("fexp", export_objects)
#'   }
#'   shinyApp(ui, server)
#' }
#' @export
# Register download handlers for each static export object in the module.
# Handler type depends on data_type and payload structure classification.
fileExportServer <- function(id, objects) {
  # Normalize once at startup; this module currently uses static objects.
  specs <- .file_export_normalize_objects(objects)

  moduleServer(id, function(input, output, session) {
    for(i in seq_along(specs)) {
      local({
        idx <- i
        spec <- specs[[idx]]
        # Build a stable base filename from the user-visible object title.
        base_name <- .file_export_safe_name(spec$title, default=sprintf("object_%d", idx))

        # Build per-row output/input ids for download and nested selection controls.
        save_id <- sprintf("save_%d", idx)
        select_id <- sprintf("select_%d", idx)
        save_all_id <- sprintf("save_all_%d", idx)

        if(spec$data_type == "rds") {
          output[[save_id]] <- downloadHandler(
            filename = function() {
              sprintf("%s.rds", base_name)
            },
            content = function(file) {
              saveRDS(spec$data, file=file)
            }
          )
          return(invisible(NULL))
        }

        if(spec$structure_type == "nested_df_list") {
          # Choice labels map UI selection to nested dataframe groups.
          choices <- .file_export_choice_names(spec$data)

          output[[save_id]] <- downloadHandler(
            filename = function() {
              selected <- input[[select_id]]
              if(is.null(selected) || is.na(selected) || !selected %in% choices) {
                selected <- choices[[1]]
              }
              sprintf("%s_%s.xlsx", base_name, .file_export_safe_name(selected, "selection"))
            },
            content = function(file) {
              selected <- input[[select_id]]
              if(is.null(selected) || is.na(selected) || !selected %in% choices) {
                selected <- choices[[1]]
              }
              # Resolve the selected label back to the nested list index.
              selected_idx <- match(selected, choices)
              .file_export_write_xlsx(
                spec$data[[selected_idx]],
                file=file,
                fallback_sheet=selected
              )
            }
          )

          output[[save_all_id]] <- downloadHandler(
            filename = function() {
              sprintf("%s_all.zip", base_name)
            },
            content = function(file) {
              .file_export_write_zip_all(spec$data, file=file)
            }
          )

          return(invisible(NULL))
        }

        output[[save_id]] <- downloadHandler(
          filename = function() {
            sprintf("%s.xlsx", base_name)
          },
          content = function(file) {
            .file_export_write_xlsx(
              spec$data,
              file=file,
              fallback_sheet=spec$title
            )
          }
        )
      })
    }
  })
}
