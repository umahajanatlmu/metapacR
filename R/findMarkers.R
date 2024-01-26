#' @title findMarkers
#'
#' @description method to select classifer markers by AUC.
#'
#' @param results anova results obtained from normalizeDat.binary function
#' @param dataList raw data list of metabolome data
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param species species to use "hsa" or "mmu"
#' @param group grouping variable, need to be one variable
#' @param p.value.cutoff p-value cutoff value
#' @param auc.threshould auc cutoff threshold
#' @param nmarkers number of markers to select
#' @param fold.changes.cutoff  cutoff for fold changes cutoff, need ony higher cutoff.
#' @param rank.plot whether to plot rank plot
#' @param dot.plot whether to plot dot plot, only if more than 1 comparisons
#' @param heatmap whether to plot heatmap
#'
#' @importFrom pROC roc
#' @importFrom ggtree ggtree layout_dendrogram
#' @importFrom RColorBrewer brewer.pal
#' @importFrom patchwork plot_spacer plot_layout wrap_plots
#' @importFrom pheatmap pheatmap
#' @import stats
#' @import graphics
#' @import grDevices
#' @import tidyverse
#' @importFrom scales squish
#'
#' @return findMarkers analyses results in list object.
#'   The object contains the following:\itemize{
#'     \item raw.results all the results of roc analysis
#'     \item metabolite.rank.plot ranking plot
#'     \item dot.plot dot.plot of markers
#'     \item heatmap heatmap distribution of the markers
#'     \item marker.metabolites list of marker metabolites
#'   }
#'
#' @export

findMarkers <- function(results,
                        dataList,
                        group,
                        data.type = c("MH", "Metabolon", "Others"),
                        species = c("hsa", "mmu"),
                        p.value.cutoff = 0.05,
                        auc.threshould = 0.60,
                        fold.changes.cutoff = 1.5,
                        nmarkers = 5,
                        rank.plot = TRUE,
                        dot.plot = TRUE,
                        heatmap = TRUE) {
  options(warn = -1) ## supress all warning

  stopifnot(inherits(dataList, "list"))
  validObject(dataList)
  species <- match.arg(species, c("hsa", "mmu"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  stopifnot(inherits(results, "data.frame"))
  validObject(results)

  if (is.null(group)) {
    stop("group variable is missing")
  } else if (length(group) != 1) {
    stop("multiple group variables available....provide only one group variable")
  }

  ## add metabolite annotation
  ## ----------------------------------------------------------------
  if (data.type == "Metabolon" && species %in% c("hsa", "mmu")) {
    data("chemicalMetadata")
    metabolite.class <- force(chemicalMetadata)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    results[["MetaboliteClass"]] <-
      metabolite.class[["SUPER_PATHWAY"]][match(results[["Metabolite"]], metabolite.class[["MET_CHEM_NO"]])]

    results <- results %>%
      full_join(metabolite.class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      rename(c("MetaboliteName" = "CHEMICAL_NAME"))
  }

  if (data.type == "MH" && species == "hsa") {

    data("chemicalMetadata_MH")
    metabolite.class <- force(chemicalMetadata_MH)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    results <- results %>%
      full_join(metabolite.class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      rename(
        c(
          "MetaboliteClass" = "ONTOLOGY1_NAME",
          "lipidClass" = "ONTOLOGY2_NAME",
          "MetaboliteName" = "METABOLITE_NAME"
        )
      )
  }

  if (data.type == "MH" && species == "mmu") {
    data("chemicalMetadata_MH_mmu")
    metabolite.class <- force(chemicalMetadata_MH_mmu)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    results <- results %>%
      full_join(metabolite.class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      rename(
        c(
          "MetaboliteClass" = "ONTOLOGY1_NAME",
          "lipidClass" = "ONTOLOGY2_NAME",
          "MetaboliteName" = "METABOLITE_NAME"
        )
      )
  }

  if (data.type == "Others") {
    metabolite.class <- Other_metadata

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    results <- results %>%
      full_join(metabolite.class, by = "Metabolite") %>%
      rename(
        c(
          "MetaboliteClass" = "Ontology_Class",
          "lipidClass" = "Ontology_Subclass",
          "MetaboliteName" = "Metabolite_Name"
        )
      )
  }

  ## load imputed data matrix
  ## ----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  ## load metadata
  ## ----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## subset metadata
  ## ----------------------------------------------------------------
  select.columns <- group
  metadata.data <- metadata.data[, colnames(metadata.data) %in% select.columns,
                                 drop = FALSE
  ]

  ## define factors
  ## ----------------------------------------------------------------
  for (c in colnames(metadata.data)) {
    if (mode(metadata.data[[c]]) %in% c("character", "factor")) {
      metadata.data[[c]] <- as.factor(metadata.data[[c]])
    } else if (mode(metadata.data[[c]]) == "difftime") {
      metadata.data[[c]] <- as.numeric(metadata.data[[c]])
    } else {
      metadata.data[[c]] <- as.numeric(metadata.data[[c]])
    }
  }

  ## merge Data
  ## ----------------------------------------------------------------
  data <- merge(metadata.data, imputed.data, by = 0) %>%
    column_to_rownames("Row.names")

  roc.results <- data.frame()
  for (m in unique(na.omit(results$contrast))) {
    results.subset <- results[results$contrast %in% m, ]
    ## get comparison group
    compGroup <- sub(" - .*", "", m)
    compGroup <- gsub(
      "^[(]", "",
      gsub("[)]$", "", compGroup)
    )


    ## p-values cutoff
    results.subset <- results.subset[results.subset$adj.P.Val < p.value.cutoff, ]
    ## subset data
    data.subset <- data[, colnames(data) %in% c(group, results.subset$Metabolite)]
    colnames(data.subset)[!colnames(data.subset) %in% group] <- results.subset$MetaboliteName[match(colnames(data.subset)[!colnames(data.subset) %in% group], results.subset$Metabolite)]
    data.subset[[group]] <- ifelse(data.subset[[group]] == compGroup, compGroup, "aaaa")

    data.subset <- data.subset[, !is.na(colnames(data.subset))]

    ## perform roc analysis per metabolite
    roc.results.subset <- data.frame()
    loop <- 1

    for (r in colnames(data.subset)[!colnames(data.subset) %in% group]) {
      rocObj <- pROC::roc(
        data.subset[[group]],
        data.subset[[r]]
      )

      roc.results.subset[loop, "metabolite"] <- r
      roc.results.subset[loop, "auc"] <- as.numeric(rocObj$auc)
      roc.results.subset[loop, "contrast"] <- m

      loop <- loop + 1
    }
    roc.results <- rbind(roc.results.subset, roc.results)
  }

  roc.results.cutoff <- roc.results %>%
    group_by(contrast) %>%
    arrange(desc(auc)) %>%
    ungroup()

  if (rank.plot == TRUE) {
    ## plot distribution
    plot.list <- list()

    for (p in unique(roc.results.cutoff$contrast)) {
      dist <- roc.results.cutoff %>%
        dplyr::filter(contrast %in% p) %>%
        arrange(desc(auc)) %>%
        slice(1:nmarkers) %>%
        mutate(
          metabolite = str_wrap(metabolite, width = 20),
          color = ifelse(.$auc > auc.threshould, "high", "low")
        ) %>%
        ggplot(aes(
          x = reorder(metabolite, -auc),
          y = auc,
          fill = color
        )) +
        # geom_segment(aes(x=reorder(metabolite, -auc),
        #                  xend=reorder(metabolite, -auc),
        #                  y=0,
        #                  yend=auc),
        #              size = 0.5,
        #              lty = "dotted") +
        geom_point(
          shape = 21,
          size = 4,
          color = "black",
          show.legend = FALSE
        ) +
        xlab("") +
        ylab("AUC") +
        geom_hline(yintercept = auc.threshould) +
        theme_bw() +
        geom_text(aes(label = metabolite), angle = 90, nudge_y = -0.1, hjust = 1) +
        scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set1")[1:2]) +
        theme(
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          axis.text = element_text(
            size = 11,
            # face = "bold",
            colour = "black"
          ),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 12, face = "bold")
        ) +
        ylim(c(0, 1)) +
        ggtitle(p)

      plot.list[[p]] <- dist
    }
    ## print distribution plots
    if (length(plot.list) > 5) {
      ncol <- 5
    } else {
      ncol <- length(plot.list)
    }

    dist.p <- patchwork::wrap_plots(plot.list, ncol = ncol)
  } else {
    dist.p <- NULL
  }



  if (dot.plot == TRUE) {
    ## select markers from auc results
    selected.metabolites <- roc.results.cutoff %>%
      group_by(contrast) %>%
      arrange(desc(auc)) %>%
      slice(1:nmarkers) %>%
      ungroup() %>%
      dplyr::select(metabolite)

    ## row dendrogram
    mat.row <- results %>%
      dplyr::filter(MetaboliteName %in% selected.metabolites$metabolite) %>%
      mutate(contrast = sub(" - .*", "", contrast)) %>%
      dplyr::filter(adj.P.Val < p.value.cutoff, logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)) %>%
      dplyr::select(contrast, MetaboliteName, logFC) %>%
      pivot_wider(names_from = contrast, values_from = logFC) %>%
      column_to_rownames("MetaboliteName") %>%
      as.matrix() %>%
      dist() %>%
      replace(is.na(.), 0)

    ## create clusters
    clust.row <- hclust(mat.row)

    ## create dendrogram
    ddgram.row <- as.dendrogram(clust.row)

    ## plot  dendrogram
    p.row <- ggtree::ggtree(ddgram.row)

    if (length(unique(na.omit(results$contrast))) < 2) {
      p.dot <- results %>%
        mutate(contrast = sub(" - .*", "", contrast)) %>%
        dplyr::filter(MetaboliteName %in% selected.metabolites$metabolite) %>%
        mutate(MetaboliteName = factor(MetaboliteName,
                                       levels = clust.row$labels[clust.row$order]
        )) %>%
        dplyr::filter(
          adj.P.Val < p.value.cutoff,
          logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)
        ) %>%
        ggplot(aes(
          x = contrast,
          y = MetaboliteName,
          color = logFC,
          size = t.ratio
        )) +
        geom_point() +
        theme_bw() +
        scale_colour_gradientn(
          colours = rev(RColorBrewer::brewer.pal(10, "RdYlBu")),
          limits = c(0, 3),
          oob = scales::squish,
          name = "log10(ratio)"
        ) +
        guides(colour = guide_colourbar(
          barwidth = unit(0.3, "cm"),
          ticks.colour = "black",
          frame.colour = "black"
        )) +
        labs(x = "", y = "", size = "t-ratio") +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        )) +
        scale_y_discrete(position = "right")

      p.dot <-

        patchwork::plot_spacer() +
        patchwork::plot_spacer() +
        p.row +
        p.dot +
        patchwork::plot_layout(
          ncol = 2,
          widths = c(2, 4),
          heights = c(1, 4)
        )
    } else {
      ## cloumn dendrogram
      mat.col <- results %>%
        dplyr::filter(MetaboliteName %in% selected.metabolites$metabolite) %>%
        mutate(contrast = sub(" - .*", "", contrast)) %>%
        dplyr::filter(adj.P.Val < p.value.cutoff, logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)) %>%
        dplyr::select(contrast, MetaboliteName, logFC) %>%
        pivot_wider(names_from = MetaboliteName, values_from = logFC) %>%
        column_to_rownames("contrast") %>%
        as.matrix() %>%
        dist() %>%
        replace(is.na(.), 0)

      ## create clusters
      clust.col <- hclust(mat.col)

      ## create dendrogram
      ddgram.col <- as.dendrogram(clust.col)

      p.col <- ggtree::ggtree(ddgram.col) + ggtree::layout_dendrogram()

      p.dot <- results %>%
        mutate(contrast = sub(" - .*", "", contrast)) %>%
        dplyr::filter(MetaboliteName %in% selected.metabolites$metabolite) %>%
        mutate(MetaboliteName = factor(MetaboliteName,
                                       levels = clust.row$labels[clust.row$order]
        )) %>%
        mutate(contrast = factor(contrast,
                                 levels = clust.col$labels[clust.col$order]
        )) %>%
        dplyr::filter(
          adj.P.Val < p.value.cutoff,
          logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)
        ) %>%
        ggplot(aes(
          x = contrast,
          y = MetaboliteName,
          color = logFC,
          size = t.ratio
        )) +
        geom_point() +
        theme_bw() +
        scale_colour_gradientn(
          colours = rev(RColorBrewer::brewer.pal(10, "RdYlBu")),
          limits = c(0, 3),
          oob = scales::squish,
          name = "log10(ratio)"
        ) +
        guides(colour = guide_colourbar(
          barwidth = unit(0.3, "cm"),
          ticks.colour = "black",
          frame.colour = "black"
        )) +
        labs(x = "", y = "", size = "t-ratio") +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        )) +
        scale_y_discrete(position = "right")

      p.dot <- patchwork::plot_spacer() +
        patchwork::plot_spacer() +
        p.col +
        patchwork::plot_spacer() +
        patchwork::plot_spacer() +
        patchwork::plot_spacer() +
        p.row +
        patchwork::plot_spacer() +
        p.dot +
        patchwork::plot_layout(
          ncol = 3,
          widths = c(2, -0.2, 4),
          heights = c(1, 0.1, 4)
        )
    }



  } else {
    p.dot <- NULL
  }

  if (heatmap == TRUE) {
    ## select markers from auc results
    selected.metabolites <- roc.results.cutoff %>%
      group_by(contrast) %>%
      arrange(desc(auc)) %>%
      slice(1:nmarkers) %>%
      dplyr::filter(auc > auc.threshould) %>%
      ungroup() %>%
      dplyr::select(metabolite)

    ## subset data for selected metabolites

    colnames(data)[!colnames(data) %in% group] <- results$MetaboliteName[match(colnames(data)[!colnames(data) %in% group], results$Metabolite)]

    data.heatmap <- data[, colnames(data) %in% c(
      group,
      selected.metabolites$metabolite
    )]
    n.data.heatmap <- data.heatmap %>%
      dplyr::select(-any_of(group)) %>%
      as.matrix()

    c.data.heatmap <- data.heatmap %>%
      dplyr::select(group) %>%
      rename(group = group)

    panel.col <-
      RColorBrewer::brewer.pal(length(sort(unique(c.data.heatmap[["group"]]))), "Dark2")
    names(panel.col) <- sort(unique(c.data.heatmap[["group"]]))

    combo.cols <-
      list(group = panel.col)

    p.heatmap <- pheatmap::pheatmap(n.data.heatmap,
                                    scale = "row",
                                    color = colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdYlBu")))(100),
                                    cluster_rows = TRUE,
                                    cluster_cols = TRUE,
                                    show_rownames = FALSE,
                                    show_colnames = TRUE,
                                    annotation_row = c.data.heatmap,
                                    annotation_colors = combo.cols
    )
  } else {
    p.heatmap <- NULL
  }

  met.list <- list()

  for (ll in unique(roc.results$contrast)) {
    roc.results.cutoff.list <- roc.results %>%
      dplyr::filter(contrast == ll) %>%
      dplyr::filter(auc > auc.threshould) %>%
      arrange(desc(auc)) %>%
      dplyr::select(metabolite)

    met.list[[ll]] <- roc.results.cutoff.list$metabolite
  }

  return(list(
    raw.results = roc.results,
    metabolite.rank.plot = dist.p,
    dot.plot = p.dot,
    heatmap = p.heatmap,
    marker.metabolites = met.list
  ))

  options(warn = 0) ## reset all warnings
}
