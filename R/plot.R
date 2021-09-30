

#' Map disco scores to colors
#'
#' Map disco scores to colors
#' @param x disco score
#' @param lower,upper lower and upper boundaries of the score
#' @param int how many colors in the palette
#' @param alpha transparency (string with two hexadecimal digits)
#' @return character vector of colors corresponding to x
#' @importFrom methods is
#' @export
disco_color_scale <- function(x, lower=-100, upper=100, int=255, alpha="66") {
  x[ is.na(x) ] <- 0
  x[ x > upper ] <- upper
  x[ x < lower ] <- lower

  pal <- colorRampPalette(c("blue", "grey", "red"))(int)
  pal <- paste0(pal, alpha)

  int <- findInterval(x, seq(lower, upper, length.out=int + 1), all.inside=TRUE)
  pal[int]
}






#' Create a disco plot
#'
#' Make a concordance / discordance plots for two contrasts
#'
#' A concordance / discordance plot is simply a log~2~ FC vs log~2~ FC
#' comparison between two contrasts. 
#' Disco score is a heuristic score 
#' reflecting the strength of similarity or dissimilarity between two
#' log~2~ FC values weighted by the corresponding p-values and given by the
#' formula log2FC.x * log2FC.y * (-log10(pval.x) - log10(pval.y)).
#'
#' To manually select which labels to show, use a truncated `annot` object
#' – remove all rows you don't want to show, and set the show_top_labels
#' parameter to `Inf`.
#' @param contrast1,contrast2 data frames with rownames corresponding to
#'        IDs (they don't need to be in the same order) and columns `log2FoldChange`
#'        and `pvalue`.
#' @param lower,upper lower and upper boundaries for coloring of the score
#' @param show_top_labels sort the genes by descending absolute disco score and show top N labels
#' @param annot annotation object returned by `get_annot()` or any other data
#'        frame with columns "PrimaryID" and "SYMBOL"
#' @param top_labels_both should top labels from both negative and positive
#'        disco scores be shown, or only for the absolute top, whether only negative
#'        or only positive or both?
#' @param alpha transparency
#' @param disco result of `disco_score`; if provided, it will be used to
#'        avoid unnecessary computations
#' @param label_col name of the column in `annot` which should be used for
#'        labels
#' @param by column by which the contrast data frames should be merged
#'       (passed to `merge`). Default: merge by row names
#' @param primary_id the name which should be assigned to the identifier
#'        column which results from the merge
#' @param label_sel identifiers of labels to be shown on the plot
#' @examples
#' ## Generate example data
#' c1 <- data.frame(log2FoldChange=rnorm(5000, sd=2))
#' c1$pvalue <- pnorm(abs(c1$log2FoldChange), sd=2, lower.tail=FALSE)
#' c2 <- data.frame(log2FoldChange=c1$log2FoldChange + rnorm(5000, sd=3))
#' c2$pvalue <- pnorm(abs(c2$log2FoldChange), sd=2, lower.tail=FALSE)
#'
#' ## Example disco plot
#' plot_disco(c1, c2)
#' @return a ggplot object (plot)
#' @import ggplot2 ggrepel
#' @importFrom stats cor
#' @export
plot_disco <- function(contrast1, contrast2, lower=-100, upper=100,
  show_top_labels=0, top_labels_both=TRUE, annot=NULL, alpha=.5, disco=NULL, by=0,
  primary_id="PrimaryID", label_col="SYMBOL", label_sel=NULL) {

  if(is.null(disco)) {
    cc <- disco_score(contrast1, contrast2, by=by, primary_id=primary_id)
  } else {
    cc <- disco
  }

  cc$col <- disco_color_scale(cc$disco, lower=lower, upper=upper)
  cc <- cc[ order(-abs(cc$disco)), ]

  if(show_top_labels > 0) {
    cc$label <- ""
    if(!is.null(annot) && all(c(primary_id, label_col) %in% colnames(annot))) {
      cc$label <- as.character(annot[[label_col]])[ match(cc[[primary_id]], annot[[primary_id]]) ]
      sel <- is.na(cc$label)
      cc$label[sel] <- cc[[primary_id]][sel]
      cc$label[is.na(cc$label)] <- ""
    } else if(is.null(cc$label <- cc[[primary_id]])) {
      cc$label <- ""
    }
    if(top_labels_both) {
      o1 <- which(cc$disco > 0)[ 1:show_top_labels ]
      o2 <- which(cc$disco < 0)[ 1:show_top_labels ]
      cc$label[ ! 1:nrow(cc) %in% c(o1, o2) ] <- ""

    } else {
      cc$label[ 1:nrow(cc) > show_top_labels ] <- ""
    }

    cc$label[is.na(cc$label)] <- ""
  }

  if(!is.null(label_sel) && length(label_sel) > 0L) {
    cc$label_sel <- ""
    if(!is.null(annot) && all(c(primary_id, label_col) %in% colnames(annot))) {
      cc$label_sel <- as.character(annot[[label_col]])[ match(cc[[primary_id]], annot[[primary_id]]) ]
      sel <- is.na(cc$label_sel)
      cc$label_sel[sel] <- cc[[primary_id]][sel]
      cc$label_sel[is.na(cc$label_sel)] <- ""
    } else if(is.null(cc$label_sel <- cc[[primary_id]])) {
      cc$label_sel <- ""
    }

    sel <- cc$label_sel %in% label_sel
    cc$label_sel[!sel] <- NA
  }


  cc <- cc %>% filter(!is.na(.data$log2FoldChange.x) & !is.na(.data$log2FoldChange.y) & !is.na(.data$disco)) %>%
    mutate(disco=ifelse(.data$disco > upper, upper, ifelse(.data$disco < lower, lower, .data$disco))) %>%
    arrange(abs(.data$disco))


  g <- ggplot(cc, aes_string(x="log2FoldChange.x", y="log2FoldChange.y")) +
    geom_point(aes(color=.data$disco), alpha=alpha) + 
    scale_color_gradient2(low="blue", mid="grey", high="red") + 
    theme(legend.position="none") +
    geom_hline(aes(yintercept=0), color="grey") +
    geom_vline(aes(xintercept=0), color="grey") +
    geom_abline(aes(slope=1, intercept=0), color="grey")

  g <- g + ggtitle(
      sprintf("Correlation between log2FC\nr=%.2f rho=%.2f",
        cor(cc$log2FoldChange.x, cc$log2FoldChange.y, use="p"),
        cor(cc$log2FoldChange.x, cc$log2FoldChange.y, method="s", use="p")))

  if(show_top_labels > 0) {
    g <- g + geom_label_repel(aes(label=.data$label, color=.data$disco,
      force_pull=1.5, max.overlaps=Inf))
  }

  if(!is.null(label_sel) && length(label_sel) > 0L) {
    g <- g + geom_label(aes(label=.data$label_sel, color=.data$disco))
  }

  return(g)
}


#' Calculate the disco score
#'
#' Merge two contrasts and calculate the disco score
#'
#' Disco score is a heuristic score 
#' reflecting the strength of similarityor dissimilarity between two
#' log~2~ FC values weighted by the corresponding p-values and given by the
#' formula log2FC.x * log2FC.y * (-log10(pval.x) - log10(pval.y)).
#' @param contrast1,contrast2 data frames with rownames corresponding to
#'        IDs (they don't need to be in the same order) and columns `log2FoldChange`
#'        and `pvalue`.
#' @param minp minimum p-value
#' @param by column by which the contrast data frames should be merged
#'       (passed to `merge`). Default: merge by row names
#' @param primary_id the name which should be assigned to the identifier
#'        column which results from the merge
#' @return a merged data frame containing column "disco.score"
#' @importFrom methods as
#' @export
disco_score <- function(contrast1, contrast2, minp=1e-16, by=0, primary_id="PrimaryID") {

  if(length(by) > 1) {
    cc <- merge(as(contrast1, "data.frame"), as(contrast2, "data.frame"), by.x=by[1], by.y=by[2])
  } else {
    cc <- merge(as(contrast1, "data.frame"), as(contrast2, "data.frame"), by=by)
  }
  colnames(cc)[1] <- "PrimaryID"
  rownames(cc) <- cc[,1]
  cc$disco <- with(cc, log2FoldChange.x * log2FoldChange.y * (-log10(pvalue.x + minp) -log10(pvalue.y  + minp)))
  cc$disco[ is.na(cc$disco) ] <- 0
  return(cc)
}



#' Create a PCA plot with plotly
#'
#' Generate a PCA 3D or 2D plot with plotly
#' @param mtx PCA matrix (typically `prcomp(...)$x`)
#' @param covariate_data data frame with covariates to display (character, factor or continuous)
#' @param threeD if FALSE, a 2D plot will be produced
#' @param cov_default default covariate names (must be names of columns of `covariate_data`)
#' @return the plotly object
#' @import purrr tibble dplyr RColorBrewer
#' @importFrom plotly plot_ly add_trace layout plotly_build
#' @importFrom grDevices colorRampPalette palette
#' @export
plot_ly_pca <- function(mtx, covariate_data, threeD=TRUE, cov_default=NULL) {

	numplots <- floor(ncol(mtx)/2)
	numplots <- min(numplots, 6)

	df <- data.frame(
		label=rownames(covariate_data),
		mtx[, paste0("PC", 1:(2*numplots))],
		stringsAsFactors=FALSE
	)

	## selecting covariates to use: not all unique and also not only one value
	sel <- apply(covariate_data, 2, function(x) {
  	.u <- unique(x)
  	length(.u) > 1 && length(.u) < nrow(covariate_data)
	})

	covariates <- union(cov_default, colnames(covariate_data)[sel])
	covariate_data <- covariate_data[ , covariates, drop=FALSE ]

	## we are sure that the order is the same
	df <- cbind(covariate_data, df)

	## which covariates are numerical?
	cov_numeric <- covariates[ map_lgl(covariates, ~ is.numeric(df[[.]])) ]

  ## make sure the remaining covariates are factors
  df <- df %>% mutate_at(setdiff(covariates, cov_numeric), factor)

  ## make the PCA data more terse
  df <- df %>% mutate_at(vars(starts_with("PC")), signif, digits=4)

  ## auto select default covariates if cov_default is missing
  if(is.null(cov_default)) {
    if(length(covariates) > 0) {
      cov_default <- covariates[1]
    }
  }

  plotly_ids <- paste0("ID", 1:nrow(df))

## First, assign palettes to factor covariates
## We need to set up coloring manually, because plotly is stupid
## use palettes from RColorBrewer, preferring the colorblind ones
  palettes <- brewer.pal.info %>% rownames_to_column("palette") %>% filter(.data$category == "qual") %>% arrange(-.data$colorblind)
  palette_names <- palettes[["palette"]]

  covariates_fac <- covariates[ sapply(covariates, function(cov) is.factor(df[[cov]])) ]
  if(nrow(palettes) < length(covariates_fac)) {
    palette_names <- rep(palette_names, ceiling(length(covariates_fac)/nrow(palettes)))
  }

  covariates_pal_names <- set_names(palette_names[1:length(covariates_fac)], covariates_fac)

## get colors for each of the factor covariates
  covariates_col <- lapply(set_names(covariates_fac), function(cov) {
    pal.name <- covariates_pal_names[cov]
    n.col <- length(levels(df[[cov]]))
    max.col <- subset(palettes, palette == pal.name)$maxcolors[1]
    .pal <- brewer.pal(n=min(max.col, max(3, n.col)), name=pal.name)
    if(max.col < n.col) {
      warning(sprintf("covariate: %s more levels (%d) than colors (%d) in palette (%s)",
        cov, n.col, max.col, pal.name))
      .pal <- rep(.pal, ceiling(n.col / max.col))
    }
    names(.pal) <- levels(df[[cov]])

    set_names(.pal[ as.character(df[[cov]]) ], plotly_ids)
  })

  cont_pal <- colorRampPalette(c("#440154FF", "#21908CFF", "#FDE725FF"))(32)
  covariates_col_cont <- lapply(set_names(setdiff(covariates, covariates_fac)), function(cov) {
    x <- df[[cov]]
    col_i <- seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=33)
    i <- findInterval(x, col_i, all.inside=TRUE)
    set_names(cont_pal[i], plotly_ids)
  })

  covariates_col <- c(covariates_col, covariates_col_cont)


## Same for symbols
  symbols <- c("circle", "square", "diamond", "cross", "x", "triangle",
               "pentagon", "hexagram", "star", "diamond", "hourglass",
               "bowtie", "asterisk", "hash")

  cov_symbols <- lapply(set_names(covariates), function(cov) {
    lev <- levels(df[[cov]])
    .s <- symbols
    if(length(lev) > length(.s)) {
      # expand symbol list if necessary. Moot, because if there are more
      # than three symbols it is pointless anyway
      .s <- rep(.s, ceiling(length(lev) / length(.s)))
    }

    .s <- .s[1:length(lev)]
    names(.s) <- lev
    return(.s)
  })

  plotly_ids <- paste0("ID", 1:nrow(df))

## create hovertext as a combination of covariates
  hovertext <- apply(map_dfc(set_names(covariates), ~ sprintf("%s: %s", ., df[[.]])), 1, paste, collapse="\n")
  hovertext <- paste(hovertext, sprintf("plotlyID: %s", plotly_ids), sep="\n")

## set up default colors and symbols

  def_symbol <- def_symbols <- def_color <- def_colors <- NULL
  if(length(cov_default) > 0) {
    def_color <- df[[ cov_default[1] ]]
    if(is.factor(def_color)) {
      levs <- levels(def_color)
      matches <- match(levs, df[[ cov_default[1] ]])
      def_colors <- covariates_col[[cov_default[1]]][ matches ]
      names(def_colors) <- levs
    } else {
      def_colors <- covariates_col[[cov_default[1]]]
    }
  }



  if(length(cov_default) > 1) {
    def_symbol <- df[[ cov_default[2] ]]
    def_symbols <- cov_symbols[[ cov_default[2] ]]
  } 

  if(threeD) {
    numplots <- floor(numplots * 2/3)
  }

## create the initial plotly object. We need it to know how plotly cuts up the data into segments
  plotly.args <- list(data=df, type="scatter", mode="markers",
               x=df$PC1, y=df$PC2,
               hovertext=hovertext,
               ids=plotly_ids,
               marker=list(size=10),
               visible=TRUE,
               color= def_color, colors= def_colors,
               symbol=def_symbol, symbols=def_symbols)

  if(threeD) {
    plotly.args$z <- df$PC3
    plotly.args$type <- "scatter3d"
    plotly.args$marker <- list(size=5)
  }

  p <- do.call("plot_ly", plotly.args)


## if there are more than two components to be shown, we need to add the data
  if(numplots > 1) {
    if(threeD) {
      for(i in 2:numplots) {
        pc1 <- paste0("PC", i * 3 - 2)
        pc2 <- paste0("PC", i * 3 - 1)
        pc3 <- paste0("PC", i * 3)
        p <- p %>% add_trace(data=df, type="scatter3d", mode="markers", x=df[[pc1]], y=df[[pc2]], z=df[[pc3]], visible=FALSE)
      }
    } else {
      for(i in 2:numplots) {
        pc1 <- paste0("PC", i * 2 - 1)
        pc2 <- paste0("PC", i * 2)
        p <- p %>% add_trace(data=df, type="scatter", mode="markers", x=df[[pc1]], y=df[[pc2]], visible=FALSE)
      }
    }
  }


## here we build the plotly object and dissect it to figure out
## how plotly cut up the data. Note that plotly is efficient and that if
## factor levels overlap between factors (i.e. factors are colinear) then
## there will be fewer blocks than simply length(levels(f1)) * length(levels(f2))
## Also, we need the IDs in each block to be able to match the correct values
  p.obj <- plotly_build(p)$x$data
  nblocks <- length(p.obj) / numplots
  p.obj <- p.obj[ 1:nblocks ]
  p.ids <- map(p.obj, ~ .$ids)

## here we build the menu for updating symbol style
  sym_update <- lapply(covariates_fac, function(cov) {
    sym_pal <- cov_symbols[[cov]]

    ## map levels to symbol
    syms <- lapply(p.ids, function(ids) {
      ord <- match(ids, paste0("ID", 1:nrow(df)))
      var <- df[[cov]][ ord ]
      sym_pal[ as.character(var) ]
    })

    list(
          method="restyle",
          label=sprintf("Symbol by %s", cov),
          args=list(list(
              marker.symbol=syms
              ))
        )
  })

## set up update menus for colors
  col_update <- lapply(covariates, function(cov) {
    cols <- lapply(p.ids, function(ids) {
      covariates_col[[cov]][ids]
    })
    list(
          method="restyle",
          label=sprintf("Color by %s", cov),
          args=list(list(
              marker.color=cols
              ))
        )
  })

## prepare the buttons to change the coordinates
  if(threeD) {
    pc_update <- lapply(1:numplots, function(i) {
      pc1 <- paste0("PC", i * 3 - 2)
      pc2 <- paste0("PC", i * 3 - 1)
      pc3 <- paste0("PC", i * 3)
      list(
        method="update",
        label=sprintf("PC%d/%d/%d", i * 3 - 2, i * 3 - 1, i * 3),
        args=list(
          list(visible=rep((1:numplots) == i, each=nblocks)),
          list(scene=list(
              xaxis=list(title=pc1),
              yaxis=list(title=pc2),
              zaxis=list(title=pc3)
            ))
        ))
    })
  } else {
    pc_update <- lapply(1:numplots, function(i) {
      pc1 <- paste0("PC", i * 2 - 1)
      pc2 <- paste0("PC", i * 2)
      list(
        method="restyle",
        label=paste(pc1, "vs", pc2),
        args=list(
          list(visible=rep((1:numplots) == i, each=nblocks))
       ))
    })
  }


  menus <- list(
    list(
      type="buttons",
      direction="right",
      y=1.1,
      yanchor="top",
      x=0.1,
      xanchor="left",
      buttons=pc_update
    ),
    list(
      yanchor="bottom",
      y=0,
      x=-.1,
      xanchor="right",
      direction="up",
      buttons=col_update
    ),
    list(
      yanchor="bottom",
      y=.1,
      xanchor="right",
      x=-.1,
      direction="up",
      buttons=sym_update
    )
  )

  scene <- list(
    xaxis=list(title="PC1"),
    yaxis=list(title="PC2")
  )

  if(threeD) {
    scene$zaxis <- list(title="PC3")
  }

  p <- p %>% layout(updatemenus=menus,
    scene=scene,
    legend=list(
      xanchor="right",
      x=-.1
    ))

  return(p)

}


#' Show gene expression in relation to a covariate
#'
#' Show gene expression in relation to a covariate
#'
#' @param id PrimaryID of the gene (usually ENSEMBL ID)
#' @param xCovar the x covariate – column name from the covariate table
#' @param exprs gene expression matrix to show on the y axis; rownames must
#'        be PrimaryIDs. If NULL, the rld object from the pipeline is used.
#' @param annot_symb_col name of the column in the annot data frame which should be added to the title of the plot.
#' @param annot_id_col name of the column in the annot data frame which corresponds to the rownames of the expression matrix. 
#' @param annot annotation data frame (as returned by the get_annot()
#'        function). If empty, it will be loaded.
#' @param covar the covariate data frame containing the column `xCovar`
#' @param groupBy name of the covariate column by which to group and connect by lines the data points 
#' @param symbolBy name of the covariate column by which to select point symbols
#' @param colorBy name of the covariate column by which to color the data
#' @param trellisBy name of the covariate column for use in a trellis (multipanel) plot
#' @return a ggplot2 object
#' @import ggplot2 
#' @export
plot_gene <- function(id, xCovar, exprs, covar, annot=NULL, 
                               annot_id_col="PrimaryID",
                               annot_symb_col="SYMBOL",
                               groupBy = NA, colorBy = NA, symbolBy = NA,
                               trellisBy=NA) {

  df <- data.frame(covar, Expression=exprs[id, ])
  if(!is.null(annot)) {
    title <- sprintf("%s (%s)", id, annot[ match(id, annot[[annot_id_col]]), ][[annot_symb_col]])
  } else {
    title <- id
  }

  if(!is.na(colorBy)) {
    if(!is.na(groupBy)) {
      g <- ggplot(df, aes(x=.data[[xCovar]], y=.data[["Expression"]], group=.data[[groupBy]], color=.data[[colorBy]]))
    } else {
      g <- ggplot(df, aes(x=.data[[xCovar]], y=.data[["Expression"]], color=.data[[colorBy]]))
    }
  } else {
    if(!is.na(groupBy)) {
      g <- ggplot(df, aes(x=.data[[xCovar]], y=.data[["Expression"]], group=.data[[groupBy]]))
    } else {
      g <- ggplot(df, aes(x=.data[[xCovar]], y=.data[["Expression"]]))
    }
  }

  if(!is.numeric(df[[xCovar]]) && is.na(groupBy)) {
    g <- g + geom_boxplot() + geom_jitter(size=3, alpha=.5, width=.1)
  } else {
    if(!is.na(symbolBy)) {
      g <- g + geom_point(aes(shape=.data[[symbolBy]], size=3))
    } else {
      g <- g + geom_point(size=3)
    }
  }

  if(!is.na(groupBy)) {
    g <- g + geom_line()
  }

  if(!is.na(trellisBy)) {
    g <- g + facet_wrap(trellisBy)
  }


  g  <- g + ggtitle(title)

  return(g)
}



