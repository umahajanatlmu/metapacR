#' oplsPlot
#'
#' Partial Least-Squares Discriminant Analysis (PLS-DA) is a multivariate dimensionality-reduction tool.
#'
#' @param data metabolome data
#' @param group stratifying variable
#' @param drop.grouping.var drop grouping varibales not necessary
#' @param ... ggplot extensions
#'
#' @import tidyverse
#' @import ropls
#' @import ggplot2
#' @import RColorBrewer
#' @import ggpubr
#' @import utils
#' @import stats
#'
#' @return
#' @export
#'
#' @examples
#' oplsPlot(data=data, group="group")
oplsPlot <- function(data = data,
                     group = group, drop.grouping.var = NULL, ...) {
  dataGroup <- data[, colnames(data) %in% group]

  var <- setdiff(colnames(data), c(group, drop.grouping.var))

  ## select numeric data, filter NaN and NAs
  dataNumeric <- data[, colnames(data) %in% var] %>%
    select_if(., is.numeric) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    select_if(~!all(is.na(.))) %>%
    as.data.frame()
  ## perform opls
  oplsSummary <-  opls(dataNumeric,
                       scaleC = "none",
                       y = dataGroup,
                       log10L = FALSE,
                       append = FALSE)
  ## plot OPLS
  oplsSummaryPlot <- data.frame(oplsSummary@scoreMN) %>%
    tibble::rownames_to_column(var="group")

  N <- nrow(oplsSummaryPlot)
  pscores <- oplsSummaryPlot[["p1"]]
  oscores <- oplsSummaryPlot[["p2"]]
  hotFisN <- (N - 1) * 2 * (N^2 - 1) / (N^2 * (N - 2)) * qf(0.95, 2, N - 2)

  ## plot
  p <-  ggplot(oplsSummaryPlot, aes(x = p1, y = p2, fill = dataGroup)) +
    gg_circle(
      rx = sqrt(var(pscores) * hotFisN),
      ry = sqrt(var(oscores) * hotFisN),
      xc = 0,
      yc = 0,
      size = 1,
      color = "black",
      fill = "white") +
    geom_hline(yintercept = 0, size = 1, color = "black") +
    geom_vline(xintercept = 0, size = 1, color = "black") +
    geom_point(
      size = 3,
      shape = 21,
      color = "black",
      alpha = 0.75) +
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
    scale_fill_manual(values = brewer.pal(length(unique(dataGroup)),
                                          "Set1")) +
    theme(legend.position = "right") +
    ggtitle(oplsSummary@typeC) +
    xlab(paste("t1 =", oplsSummary@modelDF[1,1]*100, "%", "[between groups variation] ", sprintf('\u2192'))) +
    ylab(paste("t2 =", oplsSummary@modelDF[2,1]*100, "%", "[within groups variation] ",sprintf('\u2192'))) +
    theme(panel.grid.major = element_line(colour = "grey90")) +
    annotate(geom = 'text', label = paste("R2X =",oplsSummary@summaryDF[[1]]),
             x = min(oplsSummaryPlot$p1),
             y = max(oplsSummaryPlot$p2),
             vjust = 1,
             size = 4) +
    annotate(geom = 'text', label = paste("R2Y =",oplsSummary@summaryDF[[2]]),
             x = min(oplsSummaryPlot$p1),
             y = 0.9*max(oplsSummaryPlot$p2),
             vjust = 1,
             size = 4) +
    annotate(geom = 'text', label = paste("RMSEE =",oplsSummary@summaryDF[[4]]),
             x = min(oplsSummaryPlot$p1),
             y = 0.8*max(oplsSummaryPlot$p2),
             vjust = 1,
             size = 4) +
    labs(fill = "Disease diagnosis") +
    theme(panel.background = element_rect(fill = "grey80"))

  return(list(plot=p,
              results=oplsSummary))
}


# Function to plot Hotelling's T-squared ellipse
# Adapted from https://github.com/tyrannomark/bldR/blob/master/R/L2017.R
# GPL-3 license
gg_circle <- function(rx, ry, xc, yc, color = "black", fill = NA, ...) {
  x <- xc + rx * cos(seq(0, pi, length.out = 100))
  ymax <- yc + ry * sin(seq(0, pi, length.out = 100))
  ymin <- yc + ry * sin(seq(0, -pi, length.out = 100))
  annotate(
    "ribbon",
    x = x, ymin = ymin, ymax = ymax,
    color = color, fill = fill, ...
  )
}
