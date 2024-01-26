#' @title distributionPlot
#'
#' @description plot distribution of all metabolites with fold changes and p-values.
#'
#' @param species species to use "hsa" or "mmu"
#' @param data fold changes data
#' @param path saving path
#' @param cutoff cutoff of p value
#' @param lipid.class TRUE/FALSE to plot lipid classes
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param save either "pdf", "svg" , "png", "none
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
#' @export

distributionPlot <- function(data,
                             path = NULL,
                             species = c("hsa", "mmu"),
                             cutoff = 0.01,
                             lipid.class = TRUE,
                             data.type = c("MH", "Metabolon", "Others"),
                             save = c("pdf", "svg", "png", "none"),
                             fig.width = 12,
                             fig.height = 9,
                             dpi = 300,
                             Other_metadata = NULL) {
  stopifnot(inherits(data, "data.frame"))
  validObject(data)

  species <- match.arg(species,c("hsa", "mmu"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  if (data.type == "Others") {
    stopifnot(inherits(Other_metadata, "data.frame"))
    validObject(Other_metadata)
  }

  save <- match.arg(save, c("pdf", "svg", "png", "none"))

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

  # metabolite.class <- readRDS("inst/extdata/ref/Chemical_annotations.rds")
  # use_data(metabolite.class, overwrite = TRUE)
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

  ## prepare distibution data
  dat <- data %>%
    drop_na(MetaboliteClass) %>%
    mutate(foldChanges = log2(logFC))

  groups <- unique(dat$contrast)

  ## colors
  colorsOntologyOne <-
    data.frame(
      MetaboliteClass = unique(dat$MetaboliteClass),
      color = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(
        unique(dat$MetaboliteClass)
      ))
    )
  ## match colors
  matchColumnColors <-
    match(dat$MetaboliteClass,
          colorsOntologyOne$MetaboliteClass,
          nomatch = 0
    )
  dat$color <- c("")
  dat$color[dat$MetaboliteClass %in%
              colorsOntologyOne$MetaboliteClass] <-
    as.character(colorsOntologyOne$color)[matchColumnColors]

  ## plot distribution plots

  for (i in unique(na.omit(groups))) {
    filteredData <- dat[dat$contrast %in% i, ]

    ## plot
    p1 <- ggplot(
      filteredData,
      aes(
        x = Metabolite,
        y = -log10(adj.P.Val),
        fill = MetaboliteClass
      )
    ) +
      geom_point(aes(size = logFC),
                 shape = 21,
                 alpha = 0.5,
                 show.legend = FALSE
      ) +
      facet_grid(. ~ MetaboliteClass, scales="free_x") +
      scale_fill_manual(values = unique(filteredData$color)) +
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
      theme(axis.ticks.x = element_blank()) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      ) +
      geom_hline(yintercept = -log10(cutoff), linetype = "dotted") +
      geom_hline(yintercept = -log10(cutoff * 5), linetype = "dotted") +
      ggrepel::geom_text_repel(
        aes(label = ifelse(filteredData$adj.P.Val < cutoff,
                           filteredData$MetaboliteName, NA
        )),
        min.segment.length = 0,
        size = 2,
        show.legend = FALSE
      ) +
      ggtitle(i) +
      xlab("") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    if (save == "none") {
      ## print
      print(p1)
    }
    ## save plots
    else if (save != "none") {
      ggsave(
        filename = paste(path, "/distributionPlots_", paste0(i, ".", save), sep = ""),
        plot = p1,
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
      dat <- dat %>%
        drop_na(MetaboliteClass) %>%
        mutate(foldChanges = log2(logFC)) %>%
        dplyr::filter(MetaboliteClass == "Complex lipids") %>%
        mutate(lipidClass = ifelse(grepl("^TAG", Metabolite),
                                   gsub("TAG.*", "TAG", Metabolite),
                                   gsub("[(].*", "", Metabolite)
        ))
    } else {
      dat <- dat %>%
        drop_na(MetaboliteClass) %>%
        mutate(foldChanges = log2(logFC)) %>%
        dplyr::filter(grepl("Complex lipids", MetaboliteClass))
    }

    groups <- unique(dat$contrast)

    ## colors
    colorsOntologyOne <-
      data.frame(
        lipidClass = unique(dat$lipidClass),
        color = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(
          unique(dat$lipidClass)
        ))
      )
    ## match colors
    matchColumnColors <-
      match(dat$lipidClass,
            colorsOntologyOne$lipidClass,
            nomatch = 0
      )
    dat$color <- c("")
    dat$color[dat$lipidClass %in%
                colorsOntologyOne$lipidClass] <-
      as.character(colorsOntologyOne$color)[matchColumnColors]

    ## plot distribution plots

    for (i in unique(na.omit(groups))) {
      filteredData <- dat[dat$contrast %in% i, ]

      ## plot
      p <- ggplot(
        filteredData,
        aes(
          x = Metabolite,
          y = -log10(adj.P.Val),
          fill = lipidClass
        )
      ) +
        geom_point(aes(size = logFC),
                   shape = 21,
                   alpha = 0.5,
                   show.legend = FALSE
        ) +
        facet_grid(. ~ lipidClass, scales="free_x") +
        scale_fill_manual(values = unique(filteredData$color)) +
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
        theme(axis.ticks.x = element_blank()) +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        geom_hline(yintercept = -log10(cutoff), linetype = "dotted") +
        geom_hline(yintercept = -log10(cutoff * 5), linetype = "dotted") +
        ggrepel::geom_text_repel(
          aes(label = ifelse(filteredData$adj.P.Val < cutoff,
                             filteredData$MetaboliteName, NA
          )),
          min.segment.length = 0,
          size = 2,
          show.legend = FALSE
        ) +
        ggtitle(i) +
        xlab("") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      if (save == "none") {
        ## print
        print(p)
      }
      ## save plots
      else if (save != "none") {
        ggsave(
          filename = paste(path, "/distributionPlots_lipids_", paste0(i, ".", save), sep = ""),
          plot = p,
          width = fig.width,
          height = fig.height,
          dpi = dpi
        )
      }
    }
  }
}
