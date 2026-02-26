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
