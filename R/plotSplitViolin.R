#' plotSplitViolin
#'
#' @param dataList raw data from metabolome study
#' @param grouping.variable grouping variable to plot..need to be binary or multinomial
#' @param markers list of markers to plots
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import ggpubr
#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
#' @import ggpubr
#' @import rstatix
plotSplitViolin <- function(dataList=dataList,
                            grouping.variable = grouping.variable,
                            markers= markers) {

  violin.plot <- list()

  ## load imputed data matrix
  ##----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  ## subset metadata
  ##----------------------------------------------------------------
  imputed.data <- imputed.data[, colnames(imputed.data) %in% markers,
                               drop = FALSE]

  ## load metadata
  ##----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## subset metadata
  ##----------------------------------------------------------------
  select.columns <- grouping.variable
  metadata.data <- metadata.data[, colnames(metadata.data) %in% select.columns,
                                 drop = FALSE]

  ## define factors
  ##----------------------------------------------------------------
  for (c in colnames(metadata.data)) {
    if (mode(metadata.data[[c]]) %in% c("character", "factor")) {
      metadata.data[[c]] <- as.factor(metadata.data[[c]])
    } else if (mode(metadata.data[[c]]) == "difftime") {
      metadata.data[[c]] <- as.numeric(metadata.data[[c]])
    } else
      metadata.data[[c]] <- as.numeric(metadata.data [[c]])
  }

  ## merge Data
  ##----------------------------------------------------------------
  data <- merge(metadata.data,imputed.data,by=0) %>%
    column_to_rownames("Row.names")

  grouping.factors <- unique(data[[grouping.variable]])

  for (i in grouping.factors) {
    data.subset <- data
    ## define groups
    data.subset[[grouping.variable]] <- as.factor(ifelse(data.subset[[grouping.variable]] %in% i,
                                                         i, "rest"))

    data.subset[[grouping.variable]] <- relevel(data.subset[[grouping.variable]], ref = i)

    ## gather data
    gData <- gather(data.subset, key = variable, value = value,-grouping.variable)

    name <- paste(i, "vs", "rest")

    p <- gData %>%
      mutate(metabolite = str_wrap(metabolite, width = 20)) %>%
      ggplot(aes(x=variable,
                 y=value, fill=.data[[grouping.variable]])) +
      .splitViolin() +
      stat_compare_means(aes(group = .data[[grouping.variable]]), label = "p.signif", paired = TRUE) +
      scale_fill_manual(values = brewer.pal(2, "Set1")[1:2]) +
      theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(
          size = 11,
          #face = "bold",
          colour = "black"
        ),
        axis.title = element_text(size = 12, face = "bold")
      ) +
      ggtitle(name) +
      xlab("") +
      ylab("relative changes") +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))

    violin.plot[[name]] <- p

  }

  return(plot=violin.plot)

}



#' GeomSplitViolin
#'
#' @param self auxillary function
#' @param data data
#' @param ...
#' @param draw_quantiles whether to draw quantile
#'
#' @import plyr
#' @import ggplot2
#' @import scales
#' @import grid
GeomSplitViolin <- ggproto("GeomSplitViolin",
                           GeomViolin,
                           draw_group = function(self,
                                                 data, ...,
                                                 draw_quantiles = NULL) {
                             data <- transform(data,
                                               xminv = x - violinwidth * (x - xmin),
                                               xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- arrange(transform(data,
                                                          x = if (grp %% 2 == 1) xminv
                                                          else xmaxv),
                                                if (grp %% 2 == 1) y
                                                else -y)
                             newdata <- rbind(newdata[1, ],
                                              newdata,
                                              newdata[nrow(newdata), ],
                                              newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1,
                                       nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0),
                                         all(draw_quantiles <= 1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)),
                                                  setdiff(names(data),
                                                          c("x", "y")),
                                                  drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin",
                                                grid::grobTree(GeomPolygon$draw_panel(newdata, ...),
                                                               quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin",
                                                GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

#' .splitViolin
#'
#' @param mapping ggplot mapping
#' @param data ggplot data
#' @param stat stats
#' @param position position
#' @param ... extensions
#' @param draw_quantiles TRUE/FALSE
#' @param trim TRUE/FALSE
#' @param scale scale type
#' @param na.rm TRUE/FALSE
#' @param show.legend TRUE/FALSE
#' @param inherit.aes TRUE/FALSE
#'
#' @import ggplot2

.splitViolin <- function(mapping = NULL,
                         data = NULL,
                         stat = "ydensity",
                         position = "identity", ...,
                         draw_quantiles = NULL,
                         trim = TRUE,
                         scale = "area",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  layer(data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomSplitViolin,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(trim = trim,
                      scale = scale,
                      draw_quantiles = draw_quantiles,
                      na.rm = na.rm, ...))
}

