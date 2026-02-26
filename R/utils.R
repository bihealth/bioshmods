`%||%` <- function(x, y) if(is.null(x)) y else x

# Normalize text into a filesystem-safe filename fragment.
# Uses `default` when input is missing or sanitizes to empty.
.sanitize_filename <- function(x, default="file") {
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

# Parse and sanitize a "WIDTHxHEIGHT" size string (pixels).
# Invalid input falls back to `default`, and values are clamped.
.sanitize_figsize <- function(x, default=c(800, 600), min_px=100, max_px=6000) {
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

# Open a PDF device, run drawing code, and always close the device.
.save_pdf <- function(file, width=8, height=5, draw) {
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
