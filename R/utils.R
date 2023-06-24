msg <- function(...) {
  #print(str(...))
  cat(..., "\n", file=stderr(), sep=" ")
}

# some NA values are displayed as "N/A" in the interface
is_na <- function(x) { return(is.null(x) || is.na(x) || x == "N/A") }

add_chunk <- function(id, rmd_var, type, code, title, label, ...) {
  msg("add_chunk: adding chunk of type", type, "by", id, "label", label)
  if(is.null(rmd_var)) {
    warning("add_chunk: rmd_var is NULL")
    return(NULL)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d-%H:%M:%S")
  additional <- list(...)
  chunk <- c(list(code=code, label=label, type=type, module_id=id, title=title, timestamp=timestamp), additional)

  rmd <- chunk_generate_rmd(chunk)
  chunk$chunk <- rmd

  if(is.null(rmd_var$chunks)) { rmd_var$chunks <- list() }


  msg("adding chunk\n", rmd)
  rmd_var$chunks <- c(rmd_var$chunks, list(chunk))
}

# generates a markdown chunk based on code and type, including label and
# title if provided
chunk_generate_rmd <- function(chunk, dpi=72) {
  id <- chunk$module_id
  msg(paste("generate_rmd_chunk: generating chunk of type", chunk$type, "by", id))

  if(is_na(id)) {
    stop("generate_rmd_chunk: id is NA")
    return(NULL)
  }

  if(is_na(chunk$code)) {
    stop("generate_rmd_chunk: code is NA")
    return(NULL)
  }

  if(is_na(chunk$type)) {
    stop("generate_rmd_chunk: type is NA")
    return(NULL)
  }

  if(is_na(chunk$label)) {
    stop("generate_rmd_chunk: label is NA")
    return(NULL)
  }

  if(chunk$type == "plot" && (is_na(chunk$fig.width) || is_na(chunk$fig.height)
                              || (is_na(chunk$title) && is_na(chunk$fig.cap)))) {
    stop("generate_rmd_chunk: type == plot and fig.width or fig.height is NA")
    return(NULL)
  }

  rand <- paste0(chunk$timestamp, '-', floor(runif(1, 1000, 9999)))

  ret <- "```{r "
  ret <- glue('{ret}{id}-{chunk$type}-{chunk$label}-{rand}')
  if(chunk$type == "plot") {
    if(is_na(chunk$fig.cap)) { chunk$fig.cap <- title }
    ret <- glue('{ret}, fig.width={floor(chunk$fig.width/dpi)}, fig.height={floor(chunk$fig.height/dpi)}')
  }

  ret <- glue('{ret}}}\n#| label: {id}-{chunk$type}-{chunk$label}-{rand}\n')

  if(chunk$type == "plot") {
    ret <- glue('{ret}\n#| fig-width: {floor(chunk$fig.width/dpi)}\n')
    ret <- glue('{ret}\n#| fig-height: {floor(chunk$fig.height/dpi)}\n')
    ret <- glue('{ret}\n#| fig-cap: {chunk$fig.cap}\n')
  }

  ret <- glue('{ret}\n{chunk$code}\n```\n')

  msg("generated rmd chunk:\n", ret)
  return(ret)
}


