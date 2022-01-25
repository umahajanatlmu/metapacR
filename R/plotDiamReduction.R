#' plotDiamReduction
#'
#' @param dataList raw data list of metabolome data
#' @param results results from diamntionality reduction
#' @param grouping.variables list of grouping variables
#' @param dist.variables list of metabolites
#' @param diam.method method of diamentionality reduction, pca, opls, tsne, tsne_pca, umap
#'
#' @import tidyverse
#' @import ropls
#' @import umap
#' @import Rtsne
#' @import factoextra
#' @import ggplot2
#' @import RColorBrewer
#' @import ggpubr
#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
#' @import patchwork
#' @import scales
#'
plotDiamReduction <- function(dataList = dataList,
                              results = results,
                              diam.method = "pca",
                              grouping.variables = NULL,
                              dist.variables = NULL) {
  ## load imputed data matrix
  ##----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  ## load metadata
  ##----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## add kmean clusters to metadata
  ##----------------------------------------------------------------
  cl <- results$kmeans$cluster
  metadata.data$clusters <- as.factor(cl)

  ## select row for which metadata is available
  ##----------------------------------------------------------------
  imputed.data <- imputed.data[rownames(metadata.data),]

  ## subset imputed data for markers
  ##----------------------------------------------------------------
  select.markers <- dist.variables
  imputed.data <- imputed.data[, colnames(imputed.data) %in% select.markers, drop = FALSE]

  ## subset metadata
  ##----------------------------------------------------------------
  select.columns <- c(grouping.variables, "clusters")
  metadata.data <- metadata.data[, colnames(metadata.data) %in% select.columns, drop = FALSE]

  ## group <- plotting.variable

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

  ## number of groups
  ## nVar <- length(unique(metadata.data[[group]]))

  ## plot diamentionality reductions
  plot <- results[[diam.method]]

  plot.list.group <- list()
  plot.list.dist <- list()

  if (diam.method == "pca") {
    plot.dat <- as.data.frame(plot$x)
    plot.dat <- plot.dat[,1:2]
    plot.dat <- merge(plot.dat, metadata.data, by = 0) %>%
      column_to_rownames("Row.names")
    plot.dat <- merge(plot.dat, imputed.data, by = 0)


    ## PCA componants
    pc1var <- round(summary(plot)$importance[2, 1] * 100, digits = 2)
    pc2var <- round(summary(plot)$importance[2, 2] * 100, digits = 2)

    for (pl in select.columns) {

      nVar <- length(unique(plot.dat[[pl]]))
      ## plot
      p_diam <- ggplot(plot.dat, aes(x = PC1,
                                        y = PC2,
                                        color = .data[[pl]])) +
        geom_point(
          size = 3,
          alpha = 0.8
        ) +
        xlab(paste0("PC1: ", pc1var, "%")) +
        ylab(paste0("PC2: ", pc2var, "%")) +
        labs(color = pl) +
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
        scale_color_manual(values = brewer.pal(nVar,
                                               "Set1")) +
        theme(legend.position = "right") +
        ggtitle(pl)
      ## add to list
      plot.list.group[[pl]] <- p_diam
    }

    for (dist in dist.variables) {
      ## plot
      p_dist <- ggplot(plot.dat, aes(x = PC1,
                                        y = PC2,
                                        color = .data[[dist]])) +
        geom_point(
          size = 3,
          alpha = 0.8
        ) +
        xlab(paste0("PC1: ", pc1var, "%")) +
        ylab(paste0("PC2: ", pc2var, "%")) +
        labs(color = "") +
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
        scale_color_gradientn(colours = brewer.pal(8, "Greens"),
                              limits = c(-3,3),
                              oob = scales::squish,
                              name = '') +
        theme(legend.position = "right") +
        ggtitle(dist) +
        guides(colour = guide_colourbar(barwidth = unit(0.4, "cm"),
                                        ticks.colour = "black",
                                         frame.colour = "black"))
      ## add to list
      plot.list.dist[[dist]] <- p_dist
    }

  } else if  (diam.method == "opls") {

    ## plot opls
    plot.dat <- data.frame(plot@scoreMN)
    plot.dat <- merge(plot.dat, metadata.data, by = 0) %>%
      column_to_rownames("Row.names")
    plot.dat <- merge(plot.dat, imputed.data, by = 0)

    N <- nrow(plot.dat)
    pscores <- plot.dat[["p1"]]
    oscores <- plot.dat[["p2"]]
    hotFisN <- (N - 1) * 2 * (N^2 - 1) / (N^2 * (N - 2)) * qf(0.95, 2, N - 2)

    for (pl in select.columns) {

      nVar <- length(unique(plot.dat[[pl]]))
      ## plot
      ggplot(plot.dat, aes(x = p1, y = p2, color = .data[[pl]])) +
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
          alpha = 0.75) +
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
        scale_color_manual(values = brewer.pal(nVar,"Set1")) +
        theme(legend.position = "right") +
        ggtitle(paste(plot@typeC, ":", pl)) +
        xlab(paste("t1 =", plot@modelDF[1,1]*100, "%",
                   "[between groups variation] ", sprintf('\u2192'))) +
        ylab(paste("t2 =", plot@modelDF[2,1]*100, "%",
                   "[within groups variation] ",sprintf('\u2192'))) +
        theme(panel.grid.major = element_line(colour = "grey93")) +
        annotate(geom = 'text', label = paste("R2X =",plot@summaryDF[[1]]),
                 x = min(plot.dat$p1),
                 y = max(plot.dat$p2),
                 vjust = 1,
                 size = 4) +
        annotate(geom = 'text', label = paste("R2Y =",plot@summaryDF[[2]]),
                 x = min(plot.dat$p1),
                 y = 0.9*max(plot.dat$p2),
                 vjust = 1,
                 size = 4) +
        annotate(geom = 'text', label = paste("RMSEE =",plot@summaryDF[[4]]),
                 x = min(plot.dat$p1),
                 y = 0.8*max(plot.dat$p2),
                 vjust = 1,
                 size = 4) +
        labs(color = pl) +
        theme(panel.background = element_rect(fill = "grey90"))
      ## add to list
      plot.list.group[[pl]] <- p_diam
    }

    for (dist in dist.variables) {
      ## plot
      p_dist <- ggplot(plot.dat, aes(x = p1, y = p2, color = .data[[dist]])) +
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
          alpha = 0.75) +
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
        scale_color_gradientn(colours = brewer.pal(8, "Greens"),
                              limits = c(-3,3),
                              oob = scales::squish,
                              name = '') +
        theme(legend.position = "right") +
        ggtitle(dist) +
        xlab(paste("t1 =", plot@modelDF[1,1]*100, "%",
                   "[between groups variation] ", sprintf('\u2192'))) +
        ylab(paste("t2 =", plot@modelDF[2,1]*100, "%",
                   "[within groups variation] ",sprintf('\u2192'))) +
        theme(panel.grid.major = element_line(colour = "grey93")) +
        annotate(geom = 'text', label = paste("R2X =",plot@summaryDF[[1]]),
                 x = min(plot.dat$p1),
                 y = max(plot.dat$p2),
                 vjust = 1,
                 size = 4) +
        annotate(geom = 'text', label = paste("R2Y =",plot@summaryDF[[2]]),
                 x = min(plot.dat$p1),
                 y = 0.9*max(plot.dat$p2),
                 vjust = 1,
                 size = 4) +
        annotate(geom = 'text', label = paste("RMSEE =",plot@summaryDF[[4]]),
                 x = min(plot.dat$p1),
                 y = 0.8*max(plot.dat$p2),
                 vjust = 1,
                 size = 4) +
        labs(color = "") +
        theme(panel.background = element_rect(fill = "grey90")) +
        guides(colour = guide_colourbar(barwidth = unit(0.3, "cm"),
                                        ticks.colour = "black",
                                        frame.colour = "black"))
      ## add to list
      plot.list.dist[[dist]] <- p_dist
    }

  } else if (diam.method == "umap") {

    plot.dat<- as.data.frame(plot$layout)
    plot.dat <- merge(plot.dat, metadata.data, by = 0) %>%
      column_to_rownames("Row.names")
    plot.dat <- merge(plot.dat, imputed.data, by = 0)

    for (pl in select.columns) {

      nVar <- length(unique(plot.dat[[pl]]))
      ## plot
      p_diam <- ggplot(plot.dat, aes(x = V1,
                                     y = V2,
                                     color = .data[[pl]])) +
        geom_point(
          size = 3,
          alpha = 0.8
        ) +
        xlab("UMAP1") +
        ylab("UMAP2") +
        labs(color = pl) +
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
        scale_color_manual(values = brewer.pal(nVar,
                                               "Set1")) +
        theme(legend.position = "right") +
        ggtitle(pl)
      ## add to list
      plot.list.group[[pl]] <- p_diam
    }

    for (dist in dist.variables) {
      ## plot
      p_dist <- ggplot(plot.dat, aes(x = V1,
                                     y = V2,
                                     color = .data[[dist]])) +
        geom_point(
          size = 3,
          alpha = 0.8
        ) +
        xlab("UMAP1") +
        ylab("UMAP2") +
        labs(color = "") +
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
        scale_color_gradientn(colours = brewer.pal(8, "Greens"),
                              limits = c(-3,3),
                              oob = scales::squish,
                              name = '') +
        theme(legend.position = "right") +
        ggtitle(dist) +
        guides(colour = guide_colourbar(barwidth = unit(0.4, "cm"),
                                        ticks.colour = "black",
                                        frame.colour = "black"))
      ## add to list
      plot.list.dist[[dist]] <- p_dist
    }

  }

  ## plot rtsne
  plot.dat<- as.data.frame(plot$Y)
  rownames(plot.dat) <- rownames(imputed.data)
  plot.dat <- merge(plot.dat, metadata.data, by = 0) %>%
    column_to_rownames("Row.names")
  plot.dat <- merge(plot.dat, imputed.data, by = 0)

  for (pl in select.columns) {

    nVar <- length(unique(plot.dat[[pl]]))
    ## plot
    p_diam <- ggplot(plot.dat, aes(x = V1,
                                   y = V2,
                                   color = .data[[pl]])) +
      geom_point(
        size = 3,
        alpha = 0.8
      ) +
      xlab("tSNE1") +
      ylab("tSNE2") +
      labs(color = pl) +
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
      scale_color_manual(values = brewer.pal(nVar,
                                             "Set1")) +
      theme(legend.position = "right") +
      ggtitle(pl)
    ## add to list
    plot.list.group[[pl]] <- p_diam
  }

  for (dist in dist.variables) {
    ## plot
    p_dist <- ggplot(plot.dat, aes(x = V1,
                                   y = V2,
                                   color = .data[[dist]])) +
      geom_point(
        size = 3,
        alpha = 0.8
      ) +
      xlab("tSNE1") +
      ylab("tSNE2") +
      labs(color = "") +
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
      scale_color_gradientn(colours = brewer.pal(8, "Greens"),
                            limits = c(-3,3),
                            oob = scales::squish,
                            name = '') +
      theme(legend.position = "right") +
      ggtitle(dist) +
      guides(colour = guide_colourbar(barwidth = unit(0.4, "cm"),
                                      ticks.colour = "black",
                                      frame.colour = "black"))
    ## add to list
    plot.list.dist[[dist]] <- p_dist
  }

  ## print grouping plots
  if (length(plot.list.group) > 3) {
    ncol = 3
  } else
    ncol = length(plot.list.group)

  p1 <- wrap_plots(plot.list.group, ncol = ncol)

  ## print distribution plots
  if (length(plot.list.dist) > 3) {
    ncol = 3
  } else
    ncol = length(plot.list.dist)

  p2 <- wrap_plots(plot.list.dist, ncol = ncol)

  return(list(group.plot=p1,
         distribution.plot=p2))
}
