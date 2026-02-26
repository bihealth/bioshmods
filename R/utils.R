`%||%` <- function(x, y) if(is.null(x)) y else x

#' Sanitize Text for Filenames
#'
#' Normalize free text into a filesystem-safe filename fragment by replacing
#' unsupported characters and applying a fallback when needed.
#'
#' @param x Character value to sanitize.
#' @param default Fallback used when `x` is missing, empty, or sanitizes to an
#'   empty string.
#'
#' @return A length-one character vector containing the sanitized filename
#'   fragment.
#'
#' @examples
#' sanitize_filename("Dataset A / Contrast 1")
#' sanitize_filename("", default = "untitled")
#'
#' @export
sanitize_filename <- function(x, default="file") {
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

#' Parse and Sanitize Figure Size Strings
#'
#' Parse figure size input in `"WIDTHxHEIGHT"` format and return a validated
#' size in pixels. Invalid inputs fall back to defaults and dimensions are
#' clamped to a safe range.
#'
#' @param x Character value in `"WIDTHxHEIGHT"` format.
#' @param default Numeric vector of length 2 giving fallback width and height
#'   in pixels.
#' @param min_px Minimum allowed width/height in pixels.
#' @param max_px Maximum allowed width/height in pixels.
#'
#' @return A list with numeric `width` and `height`.
#'
#' @examples
#' sanitize_figsize("1200x800")
#' sanitize_figsize("bad input", default = c(640, 480))
#'
#' @export
sanitize_figsize <- function(x, default=c(800, 600), min_px=100, max_px=6000) {
  default <- suppressWarnings(as.numeric(default))
  if(length(default) < 2L || any(!is.finite(default[1:2]))) {
    default <- c(800, 600)
  } else {
    default <- default[1:2]
  }

  default <- round(default)
  default <- pmax(min_px, pmin(max_px, default))

  x <- trimws(as.character(x %||% "")[1])
  if(is.na(x) || !nzchar(x)) {
    return(list(width=default[[1]], height=default[[2]]))
  }

  m <- regexec("^([0-9]+)\\s*[xX]\\s*([0-9]+)$", x)
  parts <- regmatches(x, m)[[1]]
  if(length(parts) != 3L) {
    return(list(width=default[[1]], height=default[[2]]))
  }

  width <- suppressWarnings(as.numeric(parts[[2]]))
  height <- suppressWarnings(as.numeric(parts[[3]]))

  if(!is.finite(width) || !is.finite(height)) {
    return(list(width=default[[1]], height=default[[2]]))
  }

  width <- round(width)
  height <- round(height)
  width <- pmax(min_px, pmin(max_px, width))
  height <- pmax(min_px, pmin(max_px, height))

  list(width=width, height=height)
}

#' Save Output to PDF Safely
#'
#' Open a PDF graphics device, execute drawing code, and always close the
#' device on exit.
#'
#' @param file Output file path.
#' @param width,height PDF width and height in inches.
#' @param draw Function or expression executed while the PDF device is open.
#'
#' @return Invisibly returns `NULL`.
#'
#' @examples
#' tf <- tempfile(fileext = ".pdf")
#' save_pdf(tf, width = 4, height = 3, draw = function() {
#'   plot(1:5, 1:5, main = "Example PDF")
#' })
#' file.exists(tf)
#'
#' @export
save_pdf <- function(file, width=8, height=5, draw) {
  width <- suppressWarnings(as.numeric(width)[1])
  height <- suppressWarnings(as.numeric(height)[1])

  if(!is.finite(width) || width <= 0) {
    width <- 8
  }

  if(!is.finite(height) || height <= 0) {
    height <- 5
  }

  grDevices::pdf(file=file, width=width, height=height)
  on.exit(grDevices::dev.off(), add=TRUE)

  if(is.function(draw)) {
    draw()
  } else {
    force(draw)
  }

  invisible(NULL)
}

#' Create a Grid of Shiny Columns
#'
#' Build a fixed grid as a `div` containing `.nrow` `fluidRow()` blocks,
#' each with `.ncol` `column()` elements.
#'
#' @param ... UI elements to place in the grid.
#' @param .ncol Number of columns per row.
#' @param .nrow Number of rows.
#' @param .byrow Logical; fill elements row-wise (`TRUE`) or column-wise (`FALSE`).
#' @param .colwidths Optional bootstrap column widths. Must have length `1` or
#'   `.ncol`, and sum to at most `12`. If `NULL`, defaults to
#'   `floor(12 / .nrow)` for each column.
#'
#' @return A `shiny.tag` (`div`) containing the requested grid.
#'
#' @examples
#' gridLayout("A", "B", .ncol=2, .nrow=1, .colwidths=6)
#' gridLayout("A", "B", "C", .ncol=2, .nrow=2, .byrow=FALSE, .colwidths=c(4, 8))
#'
#' @export
gridLayout <- function(..., .ncol=1, .nrow=1, .byrow=TRUE, .colwidths=NULL) {
  .ncol <- suppressWarnings(as.integer(.ncol)[1])
  .nrow <- suppressWarnings(as.integer(.nrow)[1])

  if(is.na(.ncol) || .ncol < 1L) {
    stop("`.ncol` must be a positive integer.")
  }
  if(is.na(.nrow) || .nrow < 1L) {
    stop("`.nrow` must be a positive integer.")
  }

  if(!is.logical(.byrow) || length(.byrow) != 1L || is.na(.byrow)) {
    stop("`.byrow` must be TRUE or FALSE.")
  }

  items <- list(...)
  capacity <- .ncol * .nrow
  if(length(items) > capacity) {
    stop("Number of UI elements in `...` cannot exceed `.ncol * .nrow`.")
  }

  if(is.null(.colwidths)) {
    .colwidths <- floor(12 / .ncol)
  }

  message(".colwidths: ", paste(.colwidths, collapse=", "))

  .colwidths <- suppressWarnings(as.integer(.colwidths))
  if(length(.colwidths) < 1L || any(is.na(.colwidths))) {
    stop("`.colwidths` must be numeric/integer.")
  }
  if(length(.colwidths) == 1L) {
    .colwidths <- rep(.colwidths[1], .ncol)
  } else if(length(.colwidths) != .ncol) {
    stop("`.colwidths` must have length 1 or `.ncol`.")
  }
  if(any(.colwidths < 1L)) {
    stop("All values in `.colwidths` must be >= 1.")
  }
  if(sum(.colwidths) > 12L) {
    stop("Sum of `.colwidths` must not be larger than 12.")
  }

  padded <- c(items, rep(list(NULL), capacity - length(items)))
  idx <- matrix(seq_len(capacity), nrow=.nrow, ncol=.ncol, byrow=isTRUE(.byrow))

  rows <- lapply(seq_len(.nrow), function(r) {
    cols <- lapply(seq_len(.ncol), function(c) {
      shiny::column(.colwidths[c], padded[[idx[r, c]]])
    })
    do.call(shiny::fluidRow, cols)
  })

  do.call(shiny::tags$div, rows)
}
