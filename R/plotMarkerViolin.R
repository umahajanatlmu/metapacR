#' @title plotMarkerViolin
#'
#' @description plot violin plot for selected markers
#'
#' @param dataList raw data list of metabolome data
#' @param grouping.variable grouping variable to plot..need to be binary or multinomial
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param species species to use "hsa" or "mmu"
#' @param markers list of markers to plots. List can be supplemented from findMarkers function.
#' @param n.markers number of markers to plots
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
#'
#' @return markers plots.
#'
#' @export

plotMarkerViolin <- function(dataList,
                             grouping.variable,
                             data.type = c("MH", "Metabolon", "Others"),
                             species = c("hsa", "mmu"),
                             markers = markers,
                             n.markers = 5) {
  stopifnot(inherits(dataList, "list"))
  validObject(dataList)
  species <- match.arg(species, c("hsa", "mmu"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  if (is.null(grouping.variable)) {
    stop("group variable is missing....please provide atleast one grouping variable")
  } else if (length(grouping.variable) != 1) {
    stop("multiple group variables available....provide only one group variable")
  }

  if (length(markers) < 5) {
    n.markers <- length(markers)
  } else {
    n.markers <- n.markers
  }

  ## add metabolite annotation
  ## ----------------------------------------------------------------
  if (data.type == "Metabolon" && species %in% c("hsa", "mmu")) {
    data("chemicalMetadata")
    metabolite.class <- force(chemicalMetadata)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character)) %>%
      rename(c("Metabolite" = "MET_CHEM_NO",
               "MetaboliteName" = "CHEMICAL_NAME"))
  }

  if (data.type == "MH" && species == "hsa") {

    data("chemicalMetadata_MH")
    metabolite.class <- force(chemicalMetadata_MH)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character)) %>%
      rename(
        c("MetaboliteClass" = "ONTOLOGY1_NAME",
          "lipidClass" = "ONTOLOGY2_NAME",
          "MetaboliteName" = "METABOLITE_NAME"
        )
      )
  }

  if (data.type == "MH" && species == "mmu") {
    data("chemicalMetadata_MH_mmu")
    metabolite.class <- force(chemicalMetadata_MH_mmu) %>%
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
      mutate(across(everything(), as.character)) %>%
      rename(
        c(
          "MetaboliteClass" = "Ontology_Class",
          "lipidClass" = "Ontology_Subclass",
          "MetaboliteName" = "Metabolite_Name"
        )
      )
  }

  violin.plot <- list()

  ## load imputed data matrix
  ## ----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  colnames(imputed.data) <- metabolite.class$MetaboliteName[match(colnames(imputed.data), metabolite.class$Metabolite)]


  ## subset metadata
  ## ----------------------------------------------------------------
  imputed.data <- imputed.data[, colnames(imputed.data) %in% markers[1:n.markers],
                               drop = FALSE
  ]

  ## load metadata
  ## ----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## subset metadata
  ## ----------------------------------------------------------------
  select.columns <- grouping.variable
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

  ## gather data
  gData <- gather(data, key = metabolite, value = value, -grouping.variable)

  p <- ggplot(
    gData,
    aes(
      x = .data[[grouping.variable]],
      y = value,
      fill = .data[[grouping.variable]]
    )
  ) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width= 0.1) +
    geom_jitter(
      shape = 21,
      size = 0.2,
      height = 0,
      width = 0.2,
      alpha = 0.2,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(gData[[grouping.variable]])), "Set1")) +
    theme_bw() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      axis.text = element_text(
        size = 11,
        # face = "bold",
        colour = "black"
      ),
      axis.title = element_text(size = 12, face = "bold")
    ) +
    ggtitle(paste("identified markers:", n.markers)) +
    xlab("") +
    ylab("Relative qunatification") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_wrap(~metabolite, ncol = 1, scales = "free_y",strip.position = "right") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    guides(fill = guide_legend(nrow = 1))


  return(p)
}
