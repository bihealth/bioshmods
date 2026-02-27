# Validate and normalize variable palette specification.
.color_palettes_validate_variables <- function(variables) {
  if(!is.list(variables) || length(variables) < 1L) {
    stop("`variables` must be a non-empty named list.")
  }

  var_names <- names(variables)
  if(is.null(var_names)) {
    stop("`variables` must be named.")
  }

  var_names <- trimws(as.character(var_names))
  if(any(is.na(var_names) | var_names == "")) {
    stop("All entries in `variables` must have non-empty names.")
  }

  out <- stats::setNames(vector("list", length(variables)), var_names)

  for(i in seq_along(variables)) {
    spec <- variables[[i]]
    var_name <- var_names[i]

    if(!is.list(spec)) {
      stop(sprintf("`variables[['%s']]` must be a list.", var_name))
    }

    type <- tolower(trimws(as.character(spec$type)[1]))
    if(is.na(type) || !type %in% c("categorical", "ordinal", "continuous")) {
      stop(sprintf(
        "`variables[['%s']]$type` must be one of: categorical, ordinal, continuous.",
        var_name
      ))
    }

    if(type %in% c("categorical", "ordinal")) {
      lev <- as.character(spec$levels)
      lev <- lev[!is.na(lev) & lev != ""]
      lev <- unique(lev)

      if(length(lev) < 1L) {
        stop(sprintf(
          "`variables[['%s']]` must define non-empty `levels` for type '%s'.",
          var_name, type
        ))
      }

      out[[i]] <- list(type=type, levels=lev)
    } else {
      br <- suppressWarnings(as.numeric(spec$breaks))
      br <- br[is.finite(br)]
      br <- unique(br)
      br <- sort(br)

      if(length(br) < 2L) {
        stop(sprintf(
          "`variables[['%s']]` must define numeric `breaks` with at least 2 values for type 'continuous'.",
          var_name
        ))
      }

      out[[i]] <- list(type=type, breaks=br)
    }
  }

  out
}

# Palette selectors per variable type.
.color_palettes_discrete_choices <- function() {
  c(
    "Brewer: Set2"="brewer:Set2",
    "Brewer: Dark2"="brewer:Dark2",
    "Brewer: Paired"="brewer:Paired",
    "Brewer: Set3"="brewer:Set3",
    "Brewer: Accent"="brewer:Accent",
    "Viridis: viridis"="viridis:viridis",
    "Viridis: cividis"="viridis:cividis",
    "Viridis: magma"="viridis:magma",
    "Viridis: plasma"="viridis:plasma",
    "Viridis: inferno"="viridis:inferno",
    "Viridis: turbo"="viridis:turbo"
  )
}

.color_palettes_ordinal_choices <- function() {
  c(
    "Brewer: YlGnBu"="brewer:YlGnBu",
    "Brewer: Blues"="brewer:Blues",
    "Brewer: Greens"="brewer:Greens",
    "Brewer: Oranges"="brewer:Oranges",
    "Brewer: Purples"="brewer:Purples",
    "Brewer: Greys"="brewer:Greys",
    "Viridis: viridis"="viridis:viridis",
    "Viridis: cividis"="viridis:cividis",
    "Viridis: magma"="viridis:magma",
    "Viridis: plasma"="viridis:plasma",
    "Viridis: inferno"="viridis:inferno",
    "Viridis: turbo"="viridis:turbo"
  )
}

.color_palettes_continuous_choices <- function() {
  c(
    "Viridis: viridis"="viridis:viridis",
    "Viridis: cividis"="viridis:cividis",
    "Viridis: magma"="viridis:magma",
    "Viridis: plasma"="viridis:plasma",
    "Viridis: inferno"="viridis:inferno",
    "Viridis: turbo"="viridis:turbo",
    "Brewer: RdBu"="brewer:RdBu",
    "Brewer: Spectral"="brewer:Spectral",
    "Brewer: PiYG"="brewer:PiYG",
    "Brewer: YlOrBr"="brewer:YlOrBr"
  )
}

.color_palettes_palette_choices <- function(type) {
  switch(
    type,
    categorical=.color_palettes_discrete_choices(),
    ordinal=.color_palettes_ordinal_choices(),
    continuous=.color_palettes_continuous_choices(),
    character(0)
  )
}

.color_palettes_default_palette <- function(type) {
  switch(
    type,
    categorical="brewer:Set2",
    ordinal="brewer:YlGnBu",
    continuous="viridis:viridis",
    "viridis:viridis"
  )
}

# Generate n colors from a Brewer or viridis-like palette selector.
.color_palettes_make_colors <- function(palette_id, n) {
  n <- suppressWarnings(as.integer(n)[1])
  if(is.na(n) || n < 1L) {
    stop("`n` must be a positive integer.")
  }

  palette_id <- as.character(palette_id)[1]
  if(is.na(palette_id) || !nzchar(palette_id)) {
    stop("`palette_id` must be a non-empty string.")
  }

  if(grepl("^brewer:", palette_id)) {
    pal_name <- sub("^brewer:", "", palette_id)
    info <- RColorBrewer::brewer.pal.info
    if(!pal_name %in% rownames(info)) {
      stop(sprintf("Unknown Brewer palette '%s'.", pal_name))
    }

    max_col <- info[pal_name, "maxcolors"]
    min_col <- 3L

    if(n <= max_col) {
      n_use <- max(min_col, n)
      cols <- RColorBrewer::brewer.pal(n_use, pal_name)[seq_len(n)]
    } else {
      base_cols <- RColorBrewer::brewer.pal(max_col, pal_name)
      cols <- grDevices::colorRampPalette(base_cols)(n)
    }

    return(cols)
  }

  if(grepl("^viridis:", palette_id)) {
    opt <- sub("^viridis:", "", palette_id)
    if(opt == "turbo") {
      return(viridisLite::turbo(n))
    }
    return(viridisLite::viridis(n, option=opt))
  }

  stop(sprintf("Unknown palette selector '%s'.", palette_id))
}

# Build one stored palette object for one variable spec.
.color_palettes_build_entry <- function(spec, palette_id, reverse=FALSE, shift=0L) {
  stopifnot(is.list(spec), is.character(spec$type), length(spec$type) == 1L)

  if(spec$type %in% c("categorical", "ordinal")) {
    lev <- spec$levels
    cols <- .color_palettes_make_colors(palette_id, length(lev))
    if(spec$type == "categorical") {
      cols <- .color_palettes_rotate(cols, shift=shift)
    }
    names(cols) <- lev
    return(list(type=spec$type, levels=lev, pal=cols))
  }

  br <- spec$breaks
  cols <- .color_palettes_make_colors(palette_id, length(br))
  if(isTRUE(reverse)) {
    cols <- rev(cols)
  }
  pal_fn <- circlize::colorRamp2(br, cols)
  list(type=spec$type, breaks=br, pal=pal_fn)
}

# Circularly rotate colors by shift positions.
.color_palettes_rotate <- function(cols, shift=0L) {
  n <- length(cols)
  if(n < 2L) {
    return(cols)
  }

  shift <- suppressWarnings(as.integer(shift)[1])
  if(is.na(shift)) {
    shift <- 0L
  }

  shift <- shift %% n
  if(shift == 0L) {
    return(cols)
  }

  c(cols[(shift + 1L):n], cols[seq_len(shift)])
}

# Pick readable text color for a colored tile.
.color_palettes_text_color <- function(bg_color) {
  rgb <- grDevices::col2rgb(bg_color)[, 1] / 255
  lum <- (0.2126 * rgb[1]) + (0.7152 * rgb[2]) + (0.0722 * rgb[3])
  if(lum < 0.5) "#FFFFFF" else "#1A1A1A"
}

# Preview for categorical / ordinal palettes.
.color_palettes_preview_discrete <- function(levels, cols) {
  tiles <- lapply(seq_along(levels), function(i) {
    this_col <- cols[i]
    this_txt <- .color_palettes_text_color(this_col)
    tags$div(
      style=sprintf(
        paste(
          "display:inline-block;",
          "margin:2px;",
          "padding:4px 8px;",
          "border-radius:3px;",
          "border:1px solid rgba(0,0,0,0.15);",
          "font-size:12px;",
          "background:%s;",
          "color:%s;"
        ),
        this_col,
        this_txt
      ),
      levels[i]
    )
  })

  do.call(tags$div, tiles)
}

# Preview for continuous palettes.
.color_palettes_preview_continuous <- function(breaks, cols) {
  labels <- lapply(seq_along(breaks), function(i) {
    tags$span(
      style="display:inline-block;min-width:36px;text-align:center;font-size:11px;color:#555;",
      format(breaks[i], trim=TRUE, scientific=FALSE)
    )
  })

  tags$div(
    tags$div(
      style=sprintf(
        paste(
          "height:16px;",
          "border:1px solid rgba(0,0,0,0.2);",
          "border-radius:3px;",
          "background:linear-gradient(to right, %s);"
        ),
        paste(cols, collapse=", ")
      )
    ),
    tags$div(
      style="display:flex;justify-content:space-between;margin-top:4px;gap:4px;",
      labels
    )
  )
}

# Build palette preview UI for one stored palette object.
.color_palettes_preview_ui <- function(entry) {
  stopifnot(is.list(entry), is.character(entry$type), length(entry$type) == 1L)

  if(entry$type %in% c("categorical", "ordinal")) {
    cols <- unname(entry$pal)
    lev <- names(entry$pal)
    return(.color_palettes_preview_discrete(levels=lev, cols=cols))
  }

  breaks <- entry$breaks
  cols <- .color_palettes_make_colors("viridis:viridis", length(breaks))
  if(is.function(entry$pal)) {
    cols <- entry$pal(breaks)
  }
  .color_palettes_preview_continuous(breaks=breaks, cols=cols)
}

# Stable internal input/output IDs per variable row.
.color_palettes_row_ids <- function(n) {
  paste0("var_", seq_len(n))
}

# Build UI for one variable palette row.
.color_palettes_row_ui <- function(ns, var_name, spec, row_id, compact=FALSE) {
  palette_input_id <- paste0("palette_", row_id)
  preview_output_id <- paste0("preview_", row_id)
  reverse_input_id <- paste0("reverse_", row_id)
  cycle_minus_id <- paste0("cycle_minus_", row_id)
  cycle_plus_id <- paste0("cycle_plus_", row_id)

  extra_controls <- NULL

  if(spec$type == "continuous") {
    extra_controls <- checkboxInput(
      ns(reverse_input_id),
      "Reverse palette",
      value=FALSE
    )
  } else if(spec$type == "categorical") {
    extra_controls <- tags$div(
      style="display:flex;align-items:center;gap:6px;margin-top:8px;",
      tags$span(style="font-size:12px;color:#666;", "Cycle colors:"),
      actionButton(ns(cycle_minus_id), "-", width="36px"),
      actionButton(ns(cycle_plus_id), "+", width="36px")
    )
  }

  description <- tags$div(
        tags$span(style="font-weight:600;", var_name),
        tags$span(
          style="font-size:11px;color:#666;margin-bottom:8px;",
          sprintf("Type: %s", spec$type)
        ))

  pal_input <- selectInput(
          ns(palette_input_id),
          "",
          choices=.color_palettes_palette_choices(spec$type),
          selected=.color_palettes_default_palette(spec$type),
          width="100%"
        )

  preview <- uiOutput(ns(preview_output_id))

  if(compact) {
    return(list(pal_input, extra_controls))
  } else {
    return(list(description, pal_input, preview, extra_controls))
  }
}

# Build the full palette controls UI for all variables.
.color_palettes_ui <- function(ns, var_names, variables, row_ids, compact=FALSE) {
  rows <- lapply(seq_along(var_names), function(i) {
    .color_palettes_row_ui(
      ns=ns,
      var_name=var_names[i],
      spec=variables[[i]],
      row_id=row_ids[i],
      compact=compact
    )
  })

  rows <- unlist(rows, recursive=FALSE)
  if(compact) {
    rows$.ncol <- 2L
    rows$.nrow <- length(var_names)
    rows$.colwidths <- c(9, 3)
  } else {
    rows$.ncol <- 4L
    rows$.nrow <- length(var_names)
    rows$.colwidths <- c(1, 5, 3, 3)
  }
  do.call(gridLayout, rows)
}

#' Color Palettes Module UI
#'
#' @param id Shiny module id.
#'
#' @return A UI container for dynamic variable-specific palette controls.
#'
#' @rdname colorPalettesServer
#' @export
colorPalettesUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("rows"))
}

#' Color Palettes Module Server
#'
#' Create palette controls for each declared variable and return the selected
#' palettes as a reactive expression.
#'
#' @param id Shiny module id (same as passed to [colorPalettesUI()]).
#' @param variables Named list of variable specifications. Each entry must
#'   contain `type` (`"categorical"`, `"ordinal"`, or `"continuous"`).
#'   For discrete types, provide `levels`. For continuous types, provide `breaks`.
#' @param palettes Optional reactiveVal-like placeholder kept for compatibility.
#'   If supplied, it must be a function; the module returns its own reactive
#'   palette specification and does not mutate external state.
#' @param compact Logical; whether to use a more compact layout for palette rows.
#'
#' @return A reactive expression returning the current palette list.
#'
#' @examples
#' if(interactive()) {
#'   variables <- list(
#'     sex=list(type="categorical", levels=c("female", "male")),
#'     stage=list(type="ordinal", levels=c("I", "II", "III", "IV")),
#'     expression=list(type="continuous", breaks=c(-2, -1, 0, 1, 2))
#'   )
#'
#'   ui <- fluidPage(
#'     colorPalettesUI("pal"),
#'     h4("Current palette keys"),
#'     verbatimTextOutput("pal_keys")
#'   )
#'
#'   server <- function(input, output, session) {
#'     palettes <- colorPalettesServer("pal", variables=variables)
#'
#'     output$pal_keys <- renderPrint(names(palettes()))
#'   }
#'
#'   shinyApp(ui, server)
#' }
#'
#' @export
colorPalettesServer <- function(id, variables, palettes = NULL, compact = FALSE) {
  variables <- .color_palettes_validate_variables(variables)

  if(is.null(palettes)) {
    palettes <- reactiveVal(list())
  }

  if(!is.function(palettes)) {
    stop("`palettes` must be a reactiveVal-like function.")
  }

  var_names <- names(variables)
  row_ids <- .color_palettes_row_ids(length(var_names))

  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    shift_state <- reactiveVal(stats::setNames(rep(0L, length(var_names)), var_names))

    output$rows <- renderUI({
      .color_palettes_ui(
          ns=ns,
          var_names=var_names,
          variables=variables,
          row_ids=row_ids,
          compact=compact
          )

      #do.call(tagList, rows)
    })

    for(i in seq_along(var_names)) {
      if(variables[[i]]$type != "categorical") {
        next
      }

      local({
        idx <- i
        var_name <- var_names[idx]
        cycle_minus_id <- paste0("cycle_minus_", row_ids[idx])
        cycle_plus_id <- paste0("cycle_plus_", row_ids[idx])

        observeEvent(input[[cycle_plus_id]], {
          x <- shift_state()
          x[[var_name]] <- x[[var_name]] + 1L
          shift_state(x)
        }, ignoreInit=TRUE)

        observeEvent(input[[cycle_minus_id]], {
          x <- shift_state()
          x[[var_name]] <- x[[var_name]] - 1L
          shift_state(x)
        }, ignoreInit=TRUE)
      })
    }

    palettes <- reactive({
      out <- stats::setNames(vector("list", length(var_names)), var_names)

      for(i in seq_along(var_names)) {
        var_name <- var_names[i]
        spec <- variables[[i]]
        palette_input_id <- paste0("palette_", row_ids[i])
        reverse_input_id <- paste0("reverse_", row_ids[i])
        selected_palette <- input[[palette_input_id]] %||% .color_palettes_default_palette(spec$type)
        out[[i]] <- .color_palettes_build_entry(
          spec=spec,
          palette_id=selected_palette,
          reverse=isTRUE(input[[reverse_input_id]]),
          shift=shift_state()[[var_name]] %||% 0L
        )
      }

      out
    })

    for(i in seq_along(var_names)) {
      local({
        idx <- i
        var_name <- var_names[idx]
        spec <- variables[[idx]]
        palette_input_id <- paste0("palette_", row_ids[idx])
        reverse_input_id <- paste0("reverse_", row_ids[idx])
        preview_output_id <- paste0("preview_", row_ids[idx])

        output[[preview_output_id]] <- renderUI({
          selected_palette <- input[[palette_input_id]] %||% .color_palettes_default_palette(spec$type)
          entry <- .color_palettes_build_entry(
            spec=spec,
            palette_id=selected_palette,
            reverse=isTRUE(input[[reverse_input_id]]),
            shift=shift_state()[[var_name]] %||% 0L
          )
          .color_palettes_preview_ui(entry)
        })
      })
    }

    return(palettes)
  })
}
