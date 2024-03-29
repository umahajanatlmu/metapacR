#' @title plotMetaboliteAlteration
#'
#' @description compute and plot distribution of signficantly altered metabolites per class.
#'
#' @param species species to use "hsa" or "mmu"
#' @param data fold changes data
#' @param path saving path
#' @param cutoff significance cutoff
#' @param lipid.class wheather to plot lipid or not
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
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
#' @return save object in defined path.
#'
#' @export

plotMetaboliteAlteration <- function(data,
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

  if (data.type == "Metabolon" && species %in% c("mmu","hsa")) {
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
    dplyr::filter(adj.P.Val < cutoff) %>%
    drop_na(MetaboliteClass) %>%
    mutate(foldChanges = log2(logFC)) %>%
    select(contrast, MetaboliteClass, logFC) %>%
    mutate(trend = ifelse(log2(logFC) > 1, "up",
                          ifelse(log2(logFC) <= -1, "down", "unchanged")
    )) %>%
    dplyr::filter(trend != "unchanged") %>%
    group_by(contrast, MetaboliteClass, trend) %>%
    summarise(FreqClass = length(trend)) %>%
    mutate(FreqClass = ifelse(trend == "down", -1 * FreqClass, FreqClass)) %>%
    ungroup() %>%
    arrange(MetaboliteClass)

  groups <- unique(na.omit(dat$contrast))


  ## plot distribution plots

  for (i in groups) {
    filteredData <- dat[dat$contrast %in% i, ]

    ## plot
    p <- ggplot(filteredData, aes(
      x = MetaboliteClass,
      y = FreqClass,
      fill = trend
    )) +
      geom_bar(
        stat = "identity",
        color = "black",
        size = 0.5
      ) +
      scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
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
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      labs(
        x = NULL, y = "Number of metabolites with sign. changes",
        title = paste("Alteration in metabolites:", i)
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
      theme(legend.position = "none") +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 0.95
      ))
    if (save == "none") {
      ## print
      print(p)
    }
    ## save plots
    if (save != "none") {
      sjPlot::save_plot(
        filename = paste(path, "/metabolite_alterationPlots_", paste0(i, ".", save), sep = ""),
        fig = p,
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
      dat <- data %>%
        drop_na(MetaboliteClass) %>%
        mutate(foldChanges = log2(logFC)) %>%
        dplyr::filter(MetaboliteClass == "Complex lipids") %>%
        mutate(lipidClass = ifelse(grepl("^TAG", Metabolite),
                                   gsub("TAG.*", "TAG", Metabolite),
                                   gsub("[(].*", "", Metabolite)
        )) %>%
        select(contrast, lipidClass, logFC) %>%
        mutate(trend = ifelse(log2(logFC) > 1, "up",
                              ifelse(log2(logFC) <= -1, "down", "unchanged")
        )) %>%
        dplyr::filter(trend != "unchanged") %>%
        group_by(contrast, lipidClass, trend) %>%
        summarise(FreqClass = length(trend)) %>%
        mutate(FreqClass = ifelse(trend == "down", -1 * FreqClass, FreqClass)) %>%
        ungroup() %>%
        arrange(lipidClass)
    } else {
      dat <- data %>%
        drop_na(MetaboliteClass) %>%
        mutate(foldChanges = log2(logFC)) %>%
        dplyr::filter(grepl("Complex lipids", MetaboliteClass)) %>%
        select(contrast, lipidClass, logFC) %>%
        mutate(trend = ifelse(log2(logFC) > 1, "up",
                              ifelse(log2(logFC) <= -1, "down", "unchanged")
        )) %>%
        dplyr::filter(trend != "unchanged") %>%
        group_by(contrast, lipidClass, trend) %>%
        summarise(FreqClass = length(trend)) %>%
        mutate(FreqClass = ifelse(trend == "down", -1 * FreqClass, FreqClass)) %>%
        ungroup() %>%
        arrange(lipidClass)
    }

    groups <- unique(na.omit(dat$contrast))

    ## plot distribution plots

    for (i in groups) {
      filteredData <- dat[dat$contrast %in% i, ]

      ## plot
      ## plot
      p1 <- ggplot(filteredData, aes(
        x = lipidClass,
        y = FreqClass,
        fill = trend
      )) +
        geom_bar(
          stat = "identity",
          color = "black",
          size = 0.5
        ) +
        scale_fill_manual(values = c("#377eb8", "#e41a1c")) +
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
        theme(
          panel.border = element_rect(colour = "black", fill = NA, size = 1),
          axis.text = element_text(
            size = 11,
            # face = "bold",
            colour = "black"
          ),
          axis.title = element_text(size = 12, face = "bold")
        ) +
        geom_hline(yintercept = 0, linetype = "solid", color = "black") +
        labs(
          x = NULL, y = "Number of metabolites with sign. changes",
          title = paste("Alteration in metabolites:", i),
          subtitle = "Complex lipids"
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
        theme(legend.position = "none") +
        theme(axis.ticks.x = element_blank()) +
        theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 0.95
        ))
      if (save == "none") {
        ## print
        print(p1)
      }
      ## save plots
      if (save != "none") {
        ggsave(
          filename = paste(path, "/metabolite_alterationPlots_lipids_", paste0(i, ".", save), sep = ""),
          plot = p1,
          width = fig.width,
          height = fig.height,
          dpi = dpi
        )
      }
    }
  }
}
