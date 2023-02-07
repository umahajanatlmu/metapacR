#' pcaPlot
#'
#' Principle component analysis
#'
#' @param data metabolome dta
#' @param group stratifying variable
#' @param drop.grouping.var drop grouping variables not necessary for PCA
#' @param ... ggplot extensions
#'
#' @import tidyverse
#' @import utils
#' @import ggplot2
#' @import RColorBrewer
#' @import ggpubr
#' @import factoextra
#' @import stats
#' @import graphics
#' @import grDevices
pcaPlot <- function (data = data,
                 group = group, drop.grouping.var = NULL, ...) {
  ## grouping variable
  dataGroup <- data[, colnames(data) %in% group, drop = FALSE]


  colnames(dataGroup) <- group

  var <- setdiff(colnames(data), c(group, drop.grouping.var))
  ## select numeric data, filter NaN and NAs
  dataNumeric <- data[, colnames(data) %in% var] %>%
    select_if(., is.numeric) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    select_if(~!all(is.na(.))) %>%
    as.data.frame()
  ## perform PCA
  pca <- prcomp(dataNumeric)
  ## plot PCA
  plot_dat <- as.data.frame(pca$x)
  plot_dat <- merge(plot_dat, dataGroup, by = 0)

  ## PCA componants
  pc1var <- round(summary(pca)$importance[2, 1] * 100, digits = 2)
  pc2var <- round(summary(pca)$importance[2, 2] * 100, digits = 2)

  ## plot
  p <- ggplot(plot_dat, aes(x = PC1,
                            y = PC2,
                            fill = .data[[group]],
                            color = .data[[group]])) +
    stat_ellipse(type = "norm",
                 lty = "dotted",
                 show.legend = FALSE) +
    geom_point(
      size = 3,
      shape = 21,
      color = "grey20",
      alpha = 0.8
    ) +
    xlab(paste0("PC1: ", pc1var, "%")) +
    ylab(paste0("PC2: ", pc2var, "%")) +
    labs(fill = group) +
    theme_bw() +
    theme(
      axis.line = element_line(size = 0.75),
      axis.text = element_text(
        size = 11,
        face = "bold",
        colour = "black"
      ),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    scale_fill_manual(values = brewer.pal(length(unique(.data[[group]])),
                                          "Set1")) +
    scale_color_manual(values = brewer.pal(length(unique(.data[[group]])),
                                           "Set1")) +
    theme(legend.position = "right")

  return(list(plot=p,
         screeplot=screeplot(pca),
         results=pca))

}
