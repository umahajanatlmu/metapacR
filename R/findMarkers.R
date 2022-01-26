#' findMarkers
#'
#' method to select classifer markers by AUC.
#'
#' @param results anova results obtained from normalizeDat.binary function
#' @param dataList raw data list of metabolome data
#' @param group grouping variable, need to be one variable
#' @param p.value.cutoff p-value cutoff value
#' @param auc.threshould auc cutoff threshold
#' @param nmarkers number of markers to select
#' @param fold.changes.cutoff  cutoff for fold changes cutoff, need ony higher cutoff.
#' @param rank.plot whether to plot rank plot
#' @param dot.plot whether to plot dot plot
#' @param heatmap whether to plot heatmap
#'
#'
#' @import pROC
#' @import ggdendra
#' @import ggtree
#' @import RColorBrewer
#' @import patchwork
#' @import pheatmap
#' @import pROC
#' @import stats
#' @import graphics
#' @import grDevices
#' @import ggplot2
#' @import tidyverse
#'
findMarkers <- function (results = results,
                         dataList = dataList,
                         group = group,
                         p.value.cutoff = 0.05,
                         auc.threshould = 0.60,
                         fold.changes.cutoff = 1.5,
                         nmarkers=5,
                         rank.plot=TRUE,
                         dot.plot=TRUE,
                         heatmap=TRUE
                         ) {
  options(warn=-1) ## supress all warning

  ## load imputed data matrix
  ##----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  ## load metadata
  ##----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## subset metadata
  ##----------------------------------------------------------------
  select.columns <- group
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

  roc.results <- data.frame()
  for (m in unique(results$contrast)) {
    results.subset <- results[results$contrast %in% m, ]
    ## get comparison group
    compGroup <- sub(" - .*", "", m)

    ## p-values cutoff
    results.subset <- results.subset %>%
      filter(adj.P.Val < p.value.cutoff)
    ## subset data
    data.subset <- data[,colnames(data) %in% c(group, results.subset$Metabolite)]
    data.subset[[group]] <- ifelse(data.subset[[group]] == compGroup, compGroup,"aaaa")
    ## perform roc analysis per metabolite
    roc.results.subset <- data.frame()
    loop <- 1

    for (r in colnames(data.subset)[!colnames(data.subset) %in% group]) {
      rocObj <- roc(data.subset[[group]],
                    data.subset[[r]])

      roc.results.subset[loop,"metabolite"] <- r
      roc.results.subset[loop,"auc"] <- as.numeric(rocObj$auc)
      roc.results.subset[loop,"contrast"] <- m

      loop <- loop + 1
    }
    roc.results <- rbind(roc.results.subset, roc.results)
  }

  roc.results.cutoff <- roc.results %>%
    group_by(contrast) %>%
    arrange(desc(auc))

  if(rank.plot = TRUE) {
    ## plot distribution
    plot.list <- list()

    for (p in unique(roc.results.cutoff$contrast)) {

      dist <- roc.results.cutoff %>%
        filter(contrast == p) %>%
        arrange(desc(auc)) %>%
        slice(1:nmarkers) %>%
        mutate(metabolite = str_wrap(metabolite, width = 20),
               color = ifelse(.$auc > auc.threshould, "high", "low")) %>%
        ggplot(aes(x=reorder(metabolite, -auc),
                   y=auc,
                   fill=color)) +
        # geom_segment(aes(x=reorder(metabolite, -auc),
        #                  xend=reorder(metabolite, -auc),
        #                  y=0,
        #                  yend=auc),
        #              size = 0.5,
        #              lty = "dotted") +
        geom_point(shape=21,
                   size=4,
                   color="black",
                   show.legend = FALSE) +
        xlab("") +
        ylab("AUC") +
        geom_hline(yintercept = auc.threshould) +
        theme_bw() +
        geom_text(aes(label=metabolite), angle = 90, nudge_y = -0.1, hjust=1) +
        scale_fill_manual(values = brewer.pal(2, "Set1")) +
        theme(
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.text = element_text(
            size = 11,
            #face = "bold",
            colour = "black"
          ),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 12, face = "bold")
        ) +
        ylim(c(0,1)) +
        ggtitle(p)

      plot.list[[p]] <- dist
    }
    ## print distribution plots
    if (length(plot.list) > 5) {
      ncol = 5
    } else {
      ncol = length(plot.list)
    }

    dist.p <- wrap_plots(plot.list, ncol = ncol)

  } else
    dist.p <- NULL

  if (dot.plot = TRUE) {

    ## select markers from auc results
    selected.metabolites <- roc.results.cutoff %>%
      group_by(contrast) %>%
      arrange(desc(auc)) %>%
      slice(1:nmarkers) %>%
      ungroup() %>%
      select(metabolite)

    ## row dendrogram
    mat.row <- results %>%
      filter(Metabolite %in% selected.metabolites$metabolite) %>%
      mutate(contrast = sub(" - .*", "", contrast)) %>%
      filter(adj.P.Val < p.value.cutoff, logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)) %>%
      select(contrast, Metabolite, logFC) %>%
      pivot_wider(names_from = contrast, values_from = logFC) %>%
      column_to_rownames("Metabolite") %>%
      as.matrix() %>%
      dist() %>%
      replace(is.na(.), 0)

    ## create clusters
    clust.row <- hclust(mat.row)

    ## create dendrogram
    ddgram.row <- as.dendrogram(clust.row)

    ## plot  dendrogram
    p.row <- ggtree(ddgram.row)

    ## cloumn dendrogram
    mat.col <- results %>%
      filter(Metabolite %in% selected.metabolites$metabolite) %>%
      mutate(contrast = sub(" - .*", "", contrast)) %>%
      filter(adj.P.Val < p.value.cutoff, logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)) %>%
      select(contrast, Metabolite, logFC) %>%
      pivot_wider(names_from = Metabolite, values_from = logFC) %>%
      column_to_rownames("contrast") %>%
      as.matrix() %>%
      dist() %>%
      replace(is.na(.), 0)

    ## create clusters
    clust.col <- hclust(mat.col)

    ## create dendrogram
    ddgram.col <- as.dendrogram(clust.col)

    p.col <- ggtree(ddgram.col) + layout_dendrogram()

    p.dot <- results %>%
      mutate(contrast = sub(" - .*", "", contrast)) %>%
      filter(Metabolite %in% selected.metabolites$metabolite) %>%
      mutate(Metabolite = factor(Metabolite,
                                 levels = clust$labels[clust$order])) %>%
      mutate(contrast = factor(contrast,
                               levels = clust.col$labels[clust.col$order])) %>%
      filter(adj.P.Val < p.value.cutoff,
             logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)) %>%
      ggplot(aes(x=contrast,
                 y = Metabolite,
                 color = logFC,
                 size = t.ratio)) +
      geom_point() +
      theme_bw() +
      scale_colour_gradientn(colours = rev(brewer.pal(10, "RdYlBu")),
                             limits = c(0,3),
                             oob = scales::squish,
                             name = 'log10(ratio)') +
      guides(colour = guide_colourbar(barwidth = unit(0.3, "cm"),
                                      ticks.colour = "black",
                                      frame.colour = "black")) +
      labs(x="", y="", size="t-ratio") +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5)) +
      scale_y_discrete(position = "right")

    # p <- plot_grid(p1, NULL, p2, nrow = 1, rel_widths = c(0.5,-0.05,2), align = 'h')

    p.dot = plot_spacer() +
      plot_spacer() +
      p.col +
      plot_spacer() +
      plot_spacer() +
      plot_spacer() +
      p.row +
      plot_spacer() +
      p.dot +
      plot_layout(ncol = 3,
                  widths = c(2, -0.2, 4),
                  heights = c(1, 0.1, 4))
  } else
    p.dot <- NULL

  if (heatmap = TRUE) {
    ## select markers from auc results
    selected.metabolites <- roc.results.cutoff %>%
      group_by(contrast) %>%
      arrange(desc(auc)) %>%
      slice(1:nmarkers) %>%
      filter(auc > auc.threshould) %>%
      ungroup() %>%
      select(metabolite)

    ## subset data for selected metabolites
    data.heatmap <- data[, colnames(data) %in% c(group,
                                                 selected.metabolites$metabolite)]
    n.data.heatmap <- data.heatmap %>%
      select(-group) %>%
      as.matrix()

    c.data.heatmap <- data.heatmap %>%
      select(group) %>%
      rename(group = group)

    panel.col <-
      brewer.pal(length(sort(unique(c.data.heatmap[["group"]]))), "Dark2")
    names(panel.col) <- sort(unique(c.data.heatmap[["group"]]))

    combo.cols <-
      list(group = panel.col)

    p.heatmap <- pheatmap(n.data.heatmap,
                  scale = "row",
                  color = colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(100),
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  show_rownames = FALSE,
                  show_colnames = FALSE,
                  annotation_row = c.data.heatmap,
                  annotation_colors = combo.cols)
  } else
    p.heatmap <- NULL

  met.list <- list()

  for (ll in unique(roc.results$contrast)) {
    roc.results.cutoff.list <- roc.results %>%
      filter(contrast == ll) %>%
      filter(auc > auc.threshould) %>%
      arrange(desc(auc)) %>%
      select(metabolite)

    met.list[[ll]] <- roc.results.cutoff.list$metabolite
  }

  return(list(
    raw.results=roc.results,
    metabolite.rank.plot=dist.p,
    dot.plot=p.dot,
    heatmap=p.heatmap,
    marker.metabolites=met.list
  ))

  options(warn=0) ## reset all warnings
}

