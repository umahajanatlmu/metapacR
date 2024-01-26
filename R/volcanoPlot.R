#' @title volcanoPlot
#'
#' @description plot volcano plots
#'
#' @param species species to use "hsa" or "mmu"
#' @param data fold changes data
#' @param path saving path
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param save either "pdf", "svg" ,"png" or "none"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi  dpi only applicable for png
#' @param Other_metadata dataframe with metadata....it must have  columns: Metabolite, Metabolite_Name, Ontology_Class, Ontology_Subclass
#'
#' @import tidyverse
#' @importFrom here here
#' @importFrom ggrepel geom_text_repel
#' @importFrom RColorBrewer brewer.pal
#' @import graphics
#' @import grDevices
#'
#' @return plots in save object in defined path
#'
#' @export

volcanoPlot <- function(data,
                        path = NULL,
                        species = c("hsa", "mmu"),
                        save = c("pdf", "svg", "png", "none"),
                        data.type = c("MH", "Metabolon", "Others"),
                        fig.width = 12,
                        fig.height = 9,
                        dpi = 300,
                        Other_metadata = NULL) {
  stopifnot(inherits(data, "data.frame"))
  validObject(data)

  species <- match.arg(species,c("hsa", "mmu"))
  save <- match.arg(save, c("pdf", "svg", "png", "none"))

  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  if (data.type == "Others") {
    stopifnot(inherits(Other_metadata, "data.frame"))
    validObject(Other_metadata)
  }

  if (is.null(path)) {
    path <- here::here()
    ifelse(!dir.exists(file.path(paste0(path), "results")),
           dir.create(file.path(paste0(path), "results")),
           FALSE
    )
    path <- paste(path, "results", sep = "/")
  } else {
    path <- path
  }


  if (data.type == "Metabolon" && species %in% c("hsa", "mmu")) {
    data("chemicalMetadata")
    metabolite.class <- force(chemicalMetadata)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    data[["MetaboliteClass"]] <- metabolite.class[["SUPER_PATHWAY"]][match(
      data[["Metabolite"]], metabolite.class[["MET_CHEM_NO"]]
    )]
    data <- data %>%
      full_join(metabolite.class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      rename(c(
        "MetaboliteName" = "CHEMICAL_NAME"
      ))
  }

  if (data.type == "MH" && species == "hsa") {
    data("chemicalMetadata_MH")
    metabolite.class <- force(chemicalMetadata_MH)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    data <- data %>%
      full_join(metabolite.class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      rename(c(
        "MetaboliteClass" = "ONTOLOGY1_NAME",
        "lipidClass" = "ONTOLOGY2_NAME",
        "MetaboliteName" = "METABOLITE_NAME"
      ))
  }

  if (data.type == "MH" && species == "mmu") {
    data("chemicalMetadata_MH_mmu")
    metabolite.class <- force(chemicalMetadata_MH_mmu)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    data <- data %>%
      full_join(metabolite.class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      rename(c(
        "MetaboliteClass" = "ONTOLOGY1_NAME",
        "lipidClass" = "ONTOLOGY2_NAME",
        "MetaboliteName" = "METABOLITE_NAME"
      ))
  }

  if (data.type == "Others") {
    metabolite.class <- Other_metadata

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    data <- data %>%
      full_join(metabolite.class, by = "Metabolite") %>%
      rename(c(
        "MetaboliteClass" = "Ontology_Class",
        "lipidClass" = "Ontology_Subclass",
        "MetaboliteName" = "Metabolite_Name"
      ))
  }

  ## prepare volcano data
  datVolcano <- data %>%
    drop_na(MetaboliteClass) %>%
    mutate(foldChanges = log2(logFC))

  groups <- unique(datVolcano$contrast)

  ## colors
  colorsOntologyOne <-
    data.frame(
      MetaboliteClass = unique(datVolcano$MetaboliteClass),
      color = RColorBrewer::brewer.pal(length(unique(
        datVolcano$MetaboliteClass
      )), "Paired")
    )
  ## match colors
  matchColumnColors <-
    match(datVolcano$MetaboliteClass,
          colorsOntologyOne$MetaboliteClass,
          nomatch = 0
    )
  datVolcano$color <- c("")
  datVolcano$color[datVolcano$MetaboliteClass %in%
                     colorsOntologyOne$MetaboliteClass] <-
    as.character(colorsOntologyOne$color)[matchColumnColors]

  ## plot volcano plots

  for (i in unique(na.omit(groups))) {
    filteredData <- datVolcano[datVolcano$contrast %in% i, ]

    ## plot
    p <-
      filteredData %>%
      ggplot(aes(foldChanges, -log10(adj.P.Val))) +
      geom_point(
        aes(fill = MetaboliteClass),
        size = 3,
        color = "black",
        pch = 21,
        alpha = 0.5
      ) +
      ggtitle(i) +
      ggrepel::geom_text_repel(
        data = head(filteredData, 5),
        aes(label = Metabolite),
        min.segment.length = 0
      ) +
      ggrepel::geom_text_repel(
        data = tail(filteredData, 5),
        aes(label = MetaboliteName),
        min.segment.length = 0
      ) +
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
      theme(legend.text.align = 0) +
      scale_fill_manual(values = unique(filteredData$color)) +
      labs(
        fill = "Metabololites category",
        x = "Relative Abundance",
        y = "p value (-log10)"
      )

    if (save == "none") {
      ## print
      print(p)
    }

    ## save plots
    if (save != "none") {
      ggsave(
        filename = paste(path, "/volcanoPlots_", paste0(i, ".", save), sep = ""),
        plot = p,
        width = fig.width,
        height = fig.height,
        dpi = dpi
      )
    }
  }
}
