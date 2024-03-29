.cleanup_ids <- function(ids) {
  ids <- gsub("_", " ", ids)

  #return(ids)

  ids <- strsplit(ids, " ")
  min.l <- min(sapply(ids, length))

  max_prefix <- 0
  for(i in 1:min.l) {

    if(length(unique(sapply(ids, function(x) x[i]))) == 1L) {
      max_prefix <- i
    }
  }

  if(max_prefix > 0) {
    ids <- lapply(ids, function(x) x[ -1:-max_prefix ])
  }

  ids <- unlist(lapply(ids, paste, collapse=" "))

  if(!any(grepl("[a-z]", ids))) {
    ids <- tolower(ids)
    substr(ids, 1, 1) <- toupper(substr(ids, 1, 1))
  }

  ids
}


#' Create a tmod panel plot using ggplot
#'
#' Create a tmod panel plot using ggplot
#' @param res a list with tmod results (each element of the list is a data
#' frame returned by a tmod test function)
#' @param pie a list with summaries of significantly DE genes by gene set.
#' Each a element of the list is a matrix returned by tmodDecideTests.
#' @param auc_thr gene sets enrichments with AUC below `auc_thr` will not be shown
#' @param q_thr gene sets enrichments with q (adjusted P value) above `q_thr` will not be shown
#' @param filter_row_q Gene sets will only be shown if at least in one
#' contrast q is smaller than `filter_row_q`
#' @param filter_row_auc Gene sets will only be shown if at least in one
#' contrast AUC is larger than `filter_row_q`
#' @param q_cutoff Any q value below `q_cutoff` will be considered equal to
#' `q_cutoff`
#' @param label_angle The angle at which column labels are shown
#' @param cleanup attempt to clean up gene set titls
#' @param add_ids add IDs of gene sets to their titles on the plot
#' @importFrom tidyr pivot_longer pivot_wider everything
#' @importFrom tibble rownames_to_column
#' @importFrom purrr flatten flatten_chr
#' @export
gg_panelplot <- function(res, pie, auc_thr=.5, q_thr=.05,
                         filter_row_q=.01,
                         filter_row_auc=.65,
                         q_cutoff=1e-12,
                         label_angle=45, cleanup=TRUE, add_ids=TRUE) {

  label_angle=as.numeric(label_angle)
  resS <- tmodSummary(res) 

  if(any(resS$Title != resS$ID)) {
    modnames <- resS$Title
    if(cleanup) {
      modnames <- .cleanup_ids(modnames)
    }
    if(add_ids) {
      modnames <- paste(resS$ID, modnames)
    }
  } else {
    modnames <- resS$ID
  }


  names(modnames) <- resS$ID

  resS_l <- resS %>%
    pivot_longer(starts_with(c("AUC", "q")), 
                 names_to=c("Param", "Contrast"), 
                 names_sep="\\.", 
                 values_to="Value") %>% 
    pivot_wider(all_of(c("ID", "Title", "Contrast")), 
                names_from="Param", 
                values_from="Value") %>%
    mutate("q" = ifelse(is.na(.data[["q"]]), 1, .data[["q"]])) %>%
    mutate("AUC" = ifelse(is.na(.data[["AUC"]]), .5, .data[["AUC"]])) %>%
    filter(.data[["AUC"]] > auc_thr & .data[["q"]] < q_thr) %>%
    mutate(q = ifelse(.data[["q"]] < q_cutoff, q_cutoff, .data[["q"]])) %>%
    mutate(alpha = -log10(.data[["q"]]) / -log10(q_cutoff))

  selMod <- resS$ID
  if(!is.na(filter_row_auc)) {
    .s <- resS_l %>% filter(.data[["AUC"]] > filter_row_auc) %>% pull("ID")
    selMod <- intersect(selMod, .s)
  }

  if(!is.na(filter_row_q)) {
    .s <- resS_l %>% filter(.data[["q"]] < filter_row_q) %>% pull("ID")
    selMod <- intersect(selMod, .s)
  }

  resS_l <- resS_l %>% filter(.data[["ID"]] %in% selMod)

  pie <- imap(pie, ~ { 
               colnames(.x) <- paste0(.y, '.', colnames(.x))
               .x %>% as.data.frame() %>% rownames_to_column("ID")
                })
  pieS <- Reduce(function(x, y) merge(x, y, all=TRUE), pie) %>%
    pivot_longer(-1, names_to=c("Contrast", "Direction"), 
                 names_sep="\\.", 
                 values_to="Number")

  df <- merge(resS_l, pieS, by=c("ID", "Contrast"), all.x=TRUE) %>%
    group_by(paste(.data[["ID"]], .data[["Contrast"]])) %>%
    mutate(Tot=sum(.data[["Number"]])) %>%
    ungroup() %>%
    mutate(Number = .data[["Number"]] * .data[["AUC"]] / .data[["Tot"]]) %>%
    mutate(Direction = factor(.data[["Direction"]], levels=c("up", "N", "down"))) 

  df <- df %>% group_by(paste0(.data[["Contrast"]], .data[["Direction"]])) %>%
    slice(match(resS$ID, .data[["ID"]])) %>%
    ungroup()

  df$Contrast <- factor(df$Contrast, levels=names(res))

  colors <- c("red", "grey", "blue")
  names(colors) <- levels(df$Direction)

  minq <- max(-log10(df[["q"]]))


  ggplot(df, aes(x=.data[["ID"]], y=.data[["Number"]], 
                 fill=.data[["Direction"]],
                 contrast=.data[["Contrast"]],
                 id=.data[["ID"]],
                 alpha=-log10(.data[["q"]]))) + 
    facet_wrap(~ .data[["Contrast"]], nrow=1, drop=FALSE) + 
    geom_bar(stat="identity") + 
    coord_flip() +
    scale_fill_manual(values=colors) +
    scale_x_discrete(breaks=names(modnames), labels=modnames) +
    theme(strip.text.x = element_text(angle=label_angle), strip.background=element_blank(),
          axis.text.x=element_text(angle=90), axis.title.y=element_blank()) +
    scale_y_continuous(breaks=c(0, .5, 1), labels=c("0", ".5", "1"), limits=c(0, 1)) +
    guides(alpha = guide_legend(override.aes = list(fill = "grey"))) +
    lims(alpha = c(0, minq)) +
    ylab("Effect size")
    
}

mp <- function(x) {
  message(paste(x, collapse=", "))
}

mf <- function(x, ...) {
  message(sprintf(x, ...))
}

#' @importFrom shinycssloaders withSpinner
#' @rdname tmodPanelPlotServer
#' @export
tmodPanelPlotUI <- function(id, datasets=NULL) {

  ttip <- list(

    contrast_select=paste("Select which contrasts are to be shown on the panel plot"),
    gene_pval=paste("Determines which genes are considered significant. Significant genes show as colored",
                    "fragments on the plot."),
    filter_auc=paste("Filter the gene sets by setting a minimal AUC threshold on the maximum AUC",
                     "in any of the contrasts. In other words, remove rows in which no test achieves at",
                     "least the AUC threshold."),
    filter_pval=paste("Filter the gene sets by setting a maximum p-value threshold on the minimum p value",
                     "in any of the contrasts. In other words, remove rows in which no test achieves ",
                     "the p-value threshold."),

    database=paste("Different gene set databases can be used for enrichment analysis.",
                   "Examples include KEGG, the pathways from Kyoto Encyclopaedia of Genes and",
                   "genomes, GO - gene ontology, Hallmark - 50 predefined gene sets from the",
                   "MSigDB and tmod, the transcritpional modules included in the tmod package."),

    sorting_order=paste("Enrichment analysis starts with a sorted list of genes. Depending on ",
                   "the details of the analysis, several approaches to sorting can be used.",
                   "By default, the genes are ordered by the p-value.")
    )

  if(is.null(datasets)) {
    datasets <- "default"
  }

  
  if(length(datasets) == 1L) {
    .inp_dataset <- hidden(selectInput(NS(id, "dataset"), label="Dataset", choices=datasets, width="100%"))
  } else {
    datasets <- c("_all", datasets)
    names(datasets) <- c("All datasets", datasets[-1])
    .inp_dataset <- fluidRow(selectInput(NS(id, "dataset"), label="Dataset", choices=datasets, width="100%"))
  }

    sidebarLayout(
      sidebarPanel(
        .inp_dataset,
        fluidRow(popify(uiOutput(NS(id, "contrast_select_field")),
                           "Select contrast", ttip$contrast_select, placement="right")),

        fluidRow(
        column(
           fluidRow(popify(uiOutput(NS(id, "db_field")),
                    "Gene set database to be shown", ttip$database, placement="right")),
           fluidRow(popify(uiOutput(NS(id, "sort_field")),
                           "Sorting order for te enrichment", ttip$sorting_order, placement="right")),
                    #selectInput(NS(id, "sort"),      label="Sorting",  choices=sorting,     width="100%"),
           fluidRow(popify(numericInput(NS(id, "gene_pval"), label="P-value significance threshold for genes", 
                                 value = 0.05, min=0, step=.01, width="100%"),
                    "P-value significance threshold for genes", ttip$gene_pval)),
           fluidRow(popify(numericInput(NS(id, "gene_lfc"), label="L2FC significance threshold for genes", 
                                 value = 0.5, min=0, step=.1, width="100%"),
                    "L2FC significance threshold for genes", ttip$gene_pval)),
           width=5),
        column(
           fluidRow(numericInput(NS(id, "font_size"), label="Font size", value = 12, 
                                 min=3, step=1, width="100%"),
                    bsTooltip(NS(id, "font_size"), "Change the font size of plot labels")),
           fluidRow(figsizeInput(NS(id, "figure_size"), width="100%"),
                bsTooltip(NS(id, "figure_size"), 
                  "Change the figure size (in pixels), width x height. Press backspace to enter your own sizes.", placement="right")),
           fluidRow(selectizeInput(NS(id, "label_angle"), label="Contrast label", 
                                choices=c(Vertical=90,
                                          Slanted=45, 
                                          Horizontal=0),
                                 width="100%"),
                bsTooltip(NS(id, "label_angle"), 
                  "How the contrast label should be displayed on the image.", placement="right")),
            
           fluidRow(numericInput(NS(id, "filter_auc"),  label="Filter by AUC (per row)",  value=0.5,
                                 min=.1, max=1, step=.05, width="100%"),
                    bsTooltip(NS(id, "filter_auc"), ttip$filter_auc)),
           fluidRow(numericInput(NS(id, "filter_pval"),  label="Filter by p-value (per row)",  
                                 value=0.05, min=0, max=1, step=0.01, width="100%"),
                    bsTooltip(NS(id, "filter_pval"), ttip$filter_auc)),
           fluidRow(downloadButton(NS(id, "save"), "Save plot to PDF", class="bg-success")),
           width=5, offset=1)),
        width=3),
      mainPanel(
          fluidRow( textOutput(NS(id, "hover_pos"))),
          column(width=12,
            withSpinner(plotOutput(NS(id, "panelPlot"), 
                                   hover=hoverOpts(NS(id, "plot_hover"), delay=50, delayType="throttle"),
                                   click=NS(id, "plot_click"),
                                   height="100%")),
           ), 
                width=9)
           



      )

 
}



#' Shiny module displaying tmod panel plots
#'
#' Shiny module displaying tmod panel plots
#' 
#' Tmod results, mapping and databases.
#'
#' The `tmod_res` object is a nested list with the following levels:
#'
#'  * top level are the contrasts. `names(tmod_res)` must be (set)equal to 
#'    `names(cntr)`.
#'  * next level are the names of tmod databases. `names(tmod_res[[1]])`
#'  must be (set)equal to `names(tmod_dbs)`
#*  * third level is the sorting or some other search parameter
#*  * lowest level is a data frame containing the actual result for the
#*    given contrast, database and sorting.
#' @param id Module ID
#' @param tmod_map mapping between the PrimaryIDs from the contrasts and
#' the gene IDs from the gene set databases.
#' @param cntr list of data frames with the results of differential
#' expression analysis. Rownames must correspond to the 'PrimaryID' column
#' of data frame annot.
#' @param annot data frame containing at least the column 'PrimaryID'
#' @param tmod_res list of tmod gene set enrichment analysis
#' results. See Details.
#' @param gs_id a "reactive values" object (returned by `reactiveValues()`), including 
#' dataset (`ds`), gene set ID (`id`), contrast id (`cntr`), database ID
#' (`db`) and sorting mode (`sort`). If `mod_id` is not `NULL`, these
#' reactive values will be populated, possibly triggering an action in
#' another shiny module.
#' @param tmod_dbs named list of tmod databases, see Details.
#' @param datasets if there are multiple data sets, this character vector
#'        defines what they are to show an approppriate selector in the UI
#' @return Returns a reactive value which is a list with elements
#' `contrast` and `id`.
#' @importFrom shiny observe selectizeInput
#' @importFrom shinyBS bsTooltip addTooltip
#' @importFrom purrr map_lgl
#' @export
tmodPanelPlotServer <- function(id, cntr, tmod_res, tmod_dbs, tmod_map, gs_id=NULL, annot=NULL) {

  if(!is.data.frame(cntr[[1]])) {
    message("tmodPanelPlotServer: cntr[[1]] is not a data frame, assuming multilevel mode")
  } else {
    cntr=list(default=cntr)
    tmod_res=list(default=tmod_res)
    tmod_dbs=list(default=tmod_dbs)
    tmod_map=list(default=tmod_map)
    annot=list(default=annot)
  }

  ds_ids       <- names(cntr)  # dataset ids
  is_single_ds <- length(ds_ids) == 1L

  ## in this module, we use labels which combine the dataset and the
  ## contrast ID, because the representation is not hierarchical, but flat
  ## this is used to map the row clicked by the user to the actual contrast
  ## name
  ## the unique identifiers are of form 'dataset: contrast'
  if(is_single_ds) {
    ds_label_map <- unlist(imap(cntr, ~ rep(.y, length(.x))))
    names(ds_label_map) <- names(cntr_label_map) <- 
                             cntr_label_map <- names(cntr[[1]])
  } else {
    cntr_label_map <- unlist(map(cntr, names))
    names(cntr_label_map) <- unlist(imap(cntr, ~ {
                                  paste0(.y, ': ', names(.x))
              }))
    ds_label_map <- unlist(imap(cntr, ~ rep(.y, length(.x))))
    names(ds_label_map) <- names(cntr_label_map)

    cntr <- imap(cntr, ~ {
                   names(.x) <- paste0(.y, ': ', names(.x))
                   .x
    })
    tmod_res <- imap(tmod_res, ~ {
                       names(.x) <- paste0(.y, ': ', names(.x))
                       .x
    })
  }

	moduleServer(id, function(input, output, session) {
    message("Launching tmod panelplot server")
    disable("save")

    observeEvent(input$dataset, {
      if(!is_single_ds) {
        .ds <- input$dataset
      } else {
        .ds <- ds_ids[1]
      }

      .ds_1 <- .ds
      if(.ds_1 == "_all") {
        .ds_1 <- ds_ids[1]
      }

      output$db_field <- renderUI({
        dbs <- names(tmod_res[[.ds_1]][[1]])
        selectInput(NS(id, "db"),  label="Database", choices=dbs, width="100%")
      })

      output$sort_field <- renderUI({
        sorting <- names(tmod_res[[.ds_1]][[1]][[1]])
        selectInput(NS(id, "sort"),  label="Sorting", choices=sorting, width="100%")
      })

      output$contrast_select_field <- renderUI({
        #if(!is_single_ds) {
        if(.ds == "_all") {
          print("ds is _all")
          choices <- c("_all", flatten_chr(map(cntr, names)))
        } else {
          print("ds is not _all")
          choices <- c("_all", names(tmod_res[[.ds]]))
        }
        names(choices) <- c("All contrasts", choices[-1])
        print("printing choices:")
        print(choices)

        selectInput(NS(id, "contrast_select"), label="Select contrasts", 
                      multiple=TRUE, choices=choices, selected="_all")

      })

    })


   ## Save figure as a PDF
   output$save <- downloadHandler(
     filename = function() {
       req(res())
       ret <- sprintf("tmod_panel_plot_%s_%s.pdf",
                      input$db, input$sort)
       ret <- gsub("[^0-9a-zA-Z_.-]", "", ret)
       return(ret)
     },
     content = function(file) {
       req(res())

      .pie <- pie()
      if(is.null(.pie)) { return(NULL) }

      .res <- flatten(res())

      if(! "_all" %in% input$contrast_select) {
        sel <- input$contrast_select
        sel <- sel[ sel %in% names(.res) ]
        .res <- .res[ sel ]
      }

      mf("Saving to file %s", file)
      pdf(file=file, width=fig_size$width / 75, height=fig_size$height / 75)
      g <- gg_panelplot(.res, pie=.pie, 
                     filter_row_auc=input$filter_auc,
                     filter_row_q=input$filter_pval,
                     label_angle=input$label_angle) + 
                                   theme(text=element_text(size=input$font_size))

       print(g)
       dev.off()
       message("returning")
     }
   )

    output$hover_pos <- renderText({

      if(is.null(input$plot_hover)) { return("Hover over the plot to identify the gene sets.") }

      contrast <- cntr_label_map[ input$plot_hover$panelvar1 ]
      dataset  <-   ds_label_map[ input$plot_hover$panelvar1 ]
      id       <- unlist(input$plot_hover$domain$discrete_limits$y)[
                                                        round(input$plot_hover$y) ]

      if(!is_single_ds) {
        ret <- sprintf("Dataset %s, Contrast %s, ID %s. Click to view in tmod browser panel.", 
              dataset, contrast, id)
      } else {
        ret <- sprintf("Contrast %s, ID %s. Click to view in tmod browser panel.", 
              contrast, id)
      }
    })


    ## prepare the list of tmod results to display
    ## a reactive expression called whenever results are required
    res <- reactive({
      if(!(isTruthy(input$db) && isTruthy(input$sort)
           && isTruthy(input$contrast_select))) { return(NULL) }

      if(input$dataset == "_all") {
        .datasets <- ds_ids
      } else {
        .datasets <- input$dataset
        warning(".datasets changed")
        warning(paste(.datasets))
      }

      # select results which correspond to the selected data sets and 
      # selected sorting order 
      ret <- imap(tmod_res[.datasets], ~ {
             .ds <- .y
             .res <- .x
             ret <- map(.res, ~ .x[[input$db]][[input$sort]]) 

             if(! "_all" %in% input$contrast_select) {
               if(any(names(ret) %in% input$contrast_select)) {
                 ret <- ret[ names(ret) %in% input$contrast_select ]
               } else {
                 ret <- NULL
               }
             }

             ret

      })

      ret <- ret[ !map_lgl(ret, is.null) ]

      return(ret)
    })

    ## Important: we make the pie for *all* datasets
    ## the returned result is flattened!
    pie <- reactive({

      # following is needed for cases when the UI is not ready yet

      if(!(isTruthy(input$db))) { return(NULL) }

      res <- res()

      if(!(isTruthy(res))) { return(NULL) }
      
      .make_pie(res(), dbname=input$db, 
                     cntr=cntr, 
                     tmod_dbs=tmod_dbs,
                     tmod_map=tmod_map,
                     gene_pval=input$gene_pval,
                     gene_lfc =input$gene_lfc)
    })

    ## record the details o the selected gene set when the user clicks on
    ## the plot. The gs_id reactive variable will hold the results.
    observeEvent(input$plot_click, {
      gs_id$db    <- input$db
      gs_id$sort  <- input$sort
      gs_id$cntr  <- cntr_label_map[ input$plot_click$panelvar1 ]
      gs_id$ds    <- ds_label_map[ input$plot_click$panelvar1 ]
      gs_id$id    <- unlist(input$plot_click$domain$discrete_limits$y)[
                      round(input$plot_click$y) ]
      gs_id$click <-  gs_id$click + 1
    })

    fig_size <- reactiveValues()

    ## record the details of the figure size selected by the user
    observeEvent(input$figure_size, {
      fig_size$width <- 
        as.numeric(gsub(" *([0-9]+) *x *([0-9]+)", "\\1", input$figure_size))

      fig_size$height <- 
        as.numeric(gsub(" *([0-9]+) *x *([0-9]+)", "\\2", input$figure_size))
    })

    observe({ output$panelPlot <- renderPlot({
      enable("save")
      #browser()

      .pie <- pie()
      if(is.null(.pie)) { 
        warning("pie() is NULL")
        return(NULL)
      }

      .res <- flatten(res())
      if(! "_all" %in% input$contrast_select) {
        sel <- input$contrast_select
        sel <- sel[ sel %in% names(.res) ]
        .res <- .res[ sel ]
      }


      #print(names(.res))

      g <- gg_panelplot(.res, pie=.pie, 
                     filter_row_auc=input$filter_auc,
                     filter_row_q=input$filter_pval,
                     label_angle=input$label_angle) + 
                                   theme(text=element_text(size=input$font_size))


      g
    }, width=fig_size$width, height=fig_size$height) })

	})

}




.make_pie <- function(tmod_res, dbname, cntr=NULL, tmod_dbs=NULL, tmod_map=NULL, 
                           gene_pval=0.05, gene_lfc=0.5) {

  if(is.null(cntr)) { return(NULL) }
  #message("Making \U1F967 pie")

  names(datasets) <- datasets <- names(tmod_res)

  pie <- map(datasets, ~ {
               ds <- .x
               .make_pie_ds(tmod_res[[ds]], dbname=dbname, 
                     cntr=cntr[[ds]], 
                     tmod_db_obj=tmod_dbs[[ds]][[dbname]],
                     tmod_map=tmod_map[[ds]],
                     gene_pval=gene_pval,
                     gene_lfc =gene_lfc)
  })

  return(flatten(pie))
}

## make pie for a single dataset
.make_pie_ds <- function(res, dataset, dbname, cntr=NULL, tmod_db_obj=NULL, tmod_map=NULL, 
                           gene_pval=0.05, gene_lfc=0.5, ...) {

  stopifnot(all(names(res) %in% names(cntr)))
  stopifnot(!is.null(tmod_db_obj) && !is.null(dbname) && !is.null(tmod_map))

  tmod_s <- tmodSummary(res)
  dbobj <- tmod_db_obj[ tmod_s[["ID"]] ]
  if(is(dbobj, "tmodGS")) {
    genes_s <- dbobj[["gv"]]
  } else {
    genes_s <- unique(unlist(dbobj[["MODULES2GENES"]]))
  }

  mp_id <- tmod_map$dbs[[dbname]]
  mp <- tmod_map$maps[[mp_id]]

  genes_sel <- unique(unlist(map(cntr, ~ .x[["PrimaryID"]])))
  genes_sel <- genes_sel[ mp[ genes_sel ] %in% genes_s ]

  ## XXXX TODO: clean up these hard encoded PrimaryIDs
  lfcs  <- map_dfc(cntr, ~ .x[ match(genes_sel, .x[["PrimaryID"]]), ][["log2FoldChange"]])
  pvals <- map_dfc(cntr, ~ .x[ match(genes_sel, .x[["PrimaryID"]]), ][["padj"]])

  ## XXX this is a workaround for a bug in tmod; in new versions it will
  ## not be necessary.
  lfcs[ is.na(lfcs) ] <- 0
  pvals[ is.na(pvals) ] <- 1

  pie <- tmodDecideTests(g = mp[ genes_sel ], lfc = lfcs, pval = pvals, mset=dbobj,
    lfc.thr = gene_lfc, pval.thr = gene_pval)

  return(pie)
}



.evidence_plot <- function(res, pie=NULL, ...) {

  tmodPanelPlot(res, pie=pie, grid="b", ...)
}

