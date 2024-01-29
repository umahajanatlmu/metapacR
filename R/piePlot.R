#' @title piePlot
#'
#' @description plot class distribution using pie charts
#'
#' @param species species to use "hsa" or "mmu"
#' @param data anova results obtained from normalizeDat.binary or normalizeDat function
#' @param path saving path
#' @param cutoff significance cutoff
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param lipid.class wheather to plot lipid or not
#' @param save either "pdf", "svg", "png" or "none"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi  dpi only applicable for png
#' @param Other_metadata dataframe with metadata....it must have  columns: Metabolite, Metabolite_Name, Ontology_Class, Ontology_Subclass
#'
#' @import tidyverse
#' @importFrom here here
#' @importFrom RColorBrewer brewer.pal
#' @import graphics
#' @import grDevices
#'
#' @return output in save format in defined path
#'
#' @export

piePlot <- function(data,
                    path = NULL,
                    species = c("hsa", "mmu"),
                    cutoff = 0.05,
                    data.type = c("MH", "Metabolon", "Others"),
                    lipid.class = TRUE,
                    save = c("pdf", "svg", "png", "none"),
                    Other_metadata = NULL,
                    fig.width = 12,
                    fig.height = 9,
                    dpi = 300) {
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

  ## prepare data
  datPie <- data %>%
    dplyr::filter(adj.P.Val < cutoff) %>%
    drop_na(MetaboliteClass) %>%
    group_by(contrast, MetaboliteClass) %>%
    summarise(Freq = length(MetaboliteClass)) %>%
    arrange(MetaboliteClass) %>%
    ungroup()

  groups <- unique(datPie$contrast)

  ## colors
  colorsOntologyOne <-
    data.frame(
      MetaboliteClass = unique(datPie$MetaboliteClass),
      color = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(
        unique(datPie$MetaboliteClass)
      ))
    )
  ## match colors
  matchColumnColors <-
    match(datPie$MetaboliteClass,
          colorsOntologyOne$MetaboliteClass,
          nomatch = 0
    )
  datPie$color <- c("")
  datPie$color[datPie$MetaboliteClass %in%
                 colorsOntologyOne$MetaboliteClass] <-
    as.character(colorsOntologyOne$color)[matchColumnColors]

  ## plot Pie plots

  for (i in unique(na.omit(groups))) {
    filteredData <- datPie[datPie$contrast %in% i, ]

    ## plot
    p <- filteredData %>%
      ggplot(aes(x = "", y = Freq)) +
      geom_bar(
        aes(fill = MetaboliteClass),
        width = 0.1,
        stat = "identity",
        color = "gray50"
      ) +
      coord_polar("y", start = 0) +
      scale_fill_manual(values = filteredData$color) +
      theme_void() +
      labs(title = paste(
        "Distribution of Metabolites:",
        i
      )) +
      theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      ) +
      guides(fill = guide_legend(ncol = 3)) +
      theme(plot.margin = unit(
        c(1, 1, 1, 1),
        "lines"
      ))
    if (save == "none") {
      ## print
      print(p)
    }

    else if (save != "none") {
      ggsave(
        filename = paste(path, "/piePlots_", paste0(i, ".", save), sep = ""),
        plot = p,
        width = fig.width,
        height = fig.height,
        dpi = dpi
      )
    }
  }

  ## lipid class
  if (isTRUE(lipid.class)) {

    ## prepare distibution data
    if (data.type == "Metabolon") {
      datPie.lipid <- data %>%
        dplyr::filter(adj.P.Val < cutoff) %>%
        drop_na(MetaboliteClass) %>%
        dplyr::filter(MetaboliteClass == "Complex lipids") %>%
        mutate(lipidClass = ifelse(grepl("^TAG", Metabolite),
                                   gsub("TAG.*", "TAG", Metabolite),
                                   gsub("[(].*", "", Metabolite)
        )) %>%
        group_by(contrast, lipidClass) %>%
        summarise(Freq = length(lipidClass)) %>%
        arrange(lipid.class) %>%
        ungroup()
    } else {
      datPie.lipid <- data %>%
        dplyr::filter(adj.P.Val < cutoff) %>%
        drop_na(MetaboliteClass) %>%
        dplyr::filter(grepl("Complex lipids", MetaboliteClass)) %>%
        group_by(contrast, lipidClass) %>%
        summarise(Freq = length(lipidClass)) %>%
        arrange(lipidClass) %>%
        ungroup()
    }


    groups <- unique(datPie.lipid$contrast)

    ## colors
    colorsOntologyOne <-
      data.frame(
        lipidClass = unique(datPie.lipid$lipidClass),
        color = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(
          unique(datPie.lipid$lipidClass)
        ))
      )
    ## match colors
    matchColumnColors <-
      match(datPie.lipid$lipidClass,
            colorsOntologyOne$lipidClass,
            nomatch = 0
      )
    datPie.lipid$color <- c("")
    datPie.lipid$color[datPie.lipid$lipidClass %in%
                         colorsOntologyOne$lipidClass] <-
      as.character(colorsOntologyOne$color)[matchColumnColors]

    ## plot Pie plots

    for (i in unique(na.omit(groups))) {
      filteredData <- datPie.lipid[datPie.lipid$contrast %in% i, ]

      ## plot
      p1 <- filteredData %>%
        ggplot(aes(x = "", y = Freq)) +
        geom_bar(
          aes(fill = lipidClass),
          width = 0.1,
          stat = "identity",
          color = "gray50"
        ) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = filteredData$color) +
        theme_void() +
        labs(title = paste(
          "Distribution of Metabolites:",
          i
        ),
        subtitle = " Complex lipids") +
        theme(
          legend.title = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)
        ) +
        guides(fill = guide_legend(ncol = 3)) +
        theme(plot.margin = unit(
          c(1, 1, 1, 1),
          "lines"
        ))
      if (save == "none") {
        ## print
        print(p1)
      }

      ## save plots
      if (save != "none") {
        ggsave(
          filename = paste(path, "/pieplot_lipids_", paste0(i, ".", save), sep = ""),
          plot = p,
          width = fig.width,
          height = fig.height,
          dpi = dpi
        )
      }
    }

  }
}
