# Resolve a report template path from the installed package or the source tree.
# This keeps template loading working during development and after installation.
.find_report_template <- function(template_name) {
  file_name <- if(grepl("\\.ya?ml$", template_name)) template_name else paste0(template_name, ".yaml")

  inst_path <- system.file("report_templates", file_name, package="bioshmods")
  if(nzchar(inst_path)) {
    return(inst_path)
  }

  search_dir <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
  while(TRUE) {
    dev_path <- file.path(search_dir, "inst", "report_templates", file_name)
    if(file.exists(dev_path)) {
      return(dev_path)
    }

    parent_dir <- dirname(search_dir)
    if(identical(parent_dir, search_dir)) {
      break
    }
    search_dir <- parent_dir
  }

  stop(sprintf("Template '%s' not found.", template_name))
}

# Match one object against the simple type labels used in the YAML templates.
# XXX -- add richer checks later if the template language grows.
.matches_template_type <- function(x, type_spec) {
  type_spec <- unlist(type_spec, use.names=FALSE)

  any(vapply(type_spec, function(.type) {
    switch(
      as.character(.type),
      any = TRUE,
      character = is.character(x),
      logical = is.logical(x),
      numeric = is.numeric(x),
      integer = is.integer(x),
      matrix = is.matrix(x),
      data.frame = is.data.frame(x),
      list = is.list(x),
      environment = is.environment(x),
      FALSE
    )
  }, logical(1)))
}

# Parse one YAML template into a plain list used by the reporting helpers.
# The structure is kept intentionally simple so the first iteration is easy to inspect.
parse_template <- function(template_name) {
  template_path <- .find_report_template(template_name)
  template <- yaml::read_yaml(template_path)

  template$required <- template$required %||% list()
  template$optional <- template$optional %||% list()
  template$template_path <- template_path
  template
}

# Check that required mappings exist, point to real objects, and match the simple type spec.
# Optional mappings are only checked when they are provided.
check_mappings <- function(env, template, mappings) {
  mappings <- mappings %||% list()

  for(.var in names(template$required)) {
    mapped_name <- mappings[[.var]]
    if(is.null(mapped_name)) {
      stop(sprintf("Missing mapping for required variable '%s'.", .var))
    }
    if(!exists(mapped_name, envir=env, inherits=FALSE)) {
      stop(sprintf("Mapped object '%s' for variable '%s' does not exist in `env`.", mapped_name, .var))
    }
    if(!.matches_template_type(env[[mapped_name]], template$required[[.var]]$type)) {
      stop(sprintf("Mapped object '%s' does not match the required type for '%s'.", mapped_name, .var))
    }
  }

  for(.var in intersect(names(template$optional), names(mappings))) {
    mapped_name <- mappings[[.var]]
    if(!exists(mapped_name, envir=env, inherits=FALSE)) {
      stop(sprintf("Mapped object '%s' for optional variable '%s' does not exist in `env`.", mapped_name, .var))
    }
    if(!.matches_template_type(env[[mapped_name]], template$optional[[.var]]$type)) {
      stop(sprintf("Mapped object '%s' does not match the required type for '%s'.", mapped_name, .var))
    }
  }

  invisible(TRUE)
}

# Turn one default specification into executable R code for the generated chunk.
# Supports either a literal default value or a literal code expression.
.template_default_code <- function(spec) {
  if(!is.null(spec$default_code)) {
    return(as.character(spec$default_code)[1])
  }

  paste(capture.output(dput(spec$default)), collapse="\n")
}

# Turn a scalar value into plain ASCII R code.
# This avoids locale-dependent quote characters in generated chunks.
.as_r_code <- function(x) {
  paste(capture.output(dput(x)), collapse="\n")
}

# Build the assignment preamble for a generated chunk.
# Required variables always come from `env`; optional ones can fall back to defaults.
.template_assignment_lines <- function(template, mappings) {
  mappings <- mappings %||% list()
  lines <- character(0)

  for(.var in names(template$required)) {
    lines <- c(lines, sprintf("%s <- env[[%s]]", .var, .as_r_code(mappings[[.var]])))
  }

  for(.var in names(template$optional)) {
    if(!is.null(mappings[[.var]])) {
      lines <- c(lines, sprintf("%s <- env[[%s]]", .var, .as_r_code(mappings[[.var]])))
    } else {
      lines <- c(lines, sprintf("%s <- %s", .var, .template_default_code(template$optional[[.var]])))
    }
  }

  lines
}

# Generate one executable text chunk from a YAML template and an environment mapping.
# The returned string is meant to be shown, evaluated in Shiny, or embedded into reports.
generate_chunk <- function(template_name, env, mappings) {
  template <- parse_template(template_name)
  check_mappings(env, template, mappings)

  chunk_builder <- get(template$chunk_builder, mode="function")
  lines <- c(
    "# Generated code chunk",
    .template_assignment_lines(template, mappings),
    "",
    chunk_builder(template)
  )

  paste(lines, collapse="\n")
}
