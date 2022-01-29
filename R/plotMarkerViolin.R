#' plotMarkerViolin
#'
#' @param dataList raw data list of metabolome data
#' @param grouping.variable grouping variable to plot..need to be binary or multinomial
#' @param markers list of markers to plots
#' @param n.markers number of markers to plots
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import ggpubr
#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
plotMarkerViolin <- function(dataList=dataList,
                             grouping.variable = grouping.variable,
                             markers=markers,
                             n.markers = 5) {

  if (length(markers) < 5) {
    n.markers <- length(markers)
  } else
    n.markers <- n.markers

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

  ## gather data
  gData <- gather(data, key = metabolite, value = value,-grouping.variable)

  p <- ggplot(gData,
              aes(x=.data[[grouping.variable]],
                  y=value,
                  fill=.data[[grouping.variable]])) +
      geom_violin(alpha = 0.5) +
    geom_jitter(shape = 21,
                size = 0.2,
                height = 0,
                width = 0.2,
                alpha = 0.2,
                show.legend = FALSE) +
      scale_fill_manual(values = brewer.pal(length(unique(gData[[grouping.variable]])), "Set1")) +
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
      ggtitle(i) +
      xlab("") +
      ylab("Relative qunatification") +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
    facet_wrap(~metabolite, ncol = 1) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 1))


  return(p)

}