#' @title lipidChainLengthDistribution
#'
#' @description compute distribution of chain length per class
#'
#' @param species species to use "hsa" or "mmu"
#' @param results fold changes data
#' @param p.value.cutoff cutoff of p-values to be used
#' @param fold.changes.cutoff higher fold changes cutoff to be used
#' @param path saving path
#' @param save either "pdf", "svg" or "png"
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi dpi only applicable for png
#' @param Other_metadata dataframe with metadata....it must have  columns: Metabolite, Metabolite_Name, Ontology_Class, Ontology_Subclass
#'
#' @import tidyverse
#' @importFrom here here
#' @importFrom sjPlot save_plot
#' @importFrom RColorBrewer brewer.pal
#' @import graphics
#' @import grDevices
#' @importFrom scales squish
#'
#' @return plot in save object in defined path.
#'
#' @export

lipidChainLengthDistribution <- function(results,
                                         species = c("hsa", "mmu"),
                                         p.value.cutoff = 0.05,
                                         fold.changes.cutoff = 1.5,
                                         path = NULL,
                                         save = c("pdf", "svg", "png"),
                                         data.type = c("MH", "Metabolon", "Others"),
                                         fig.width = 12,
                                         fig.height = 9,
                                         dpi = 300,
                                         Other_metadata = NULL) {
  stopifnot(inherits(results, "data.frame"))
  validObject(results)

  species <- match.arg(species,c("hsa", "mmu"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  if (data.type == "Others") {
    stopifnot(inherits(Other_metadata, "data.frame"))
    validObject(Other_metadata)
  }

  save <- match.arg(save, c("pdf", "svg", "png"))

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

  if (save == "pdf") {
    pdf(paste(path, "chainLengthDistribution.pdf", sep = "/"),
      paper = "a4r",
      onefile = TRUE
    )
  } else if (save != "pdf") {
    dir.create(paste(here(), "chainLengthDistribution", sep = "/"))
  }


  ## load annotation file
  if (data.type == "Metabolon" && species == "hsa") {
    data("chemicalMetadata")
    metabolite.class <- force(chemicalMetadata)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character))

    ## define metabolites
    results[["MetaboliteClass"]] <- metabolite.class[["SUPER_PATHWAY"]][match(
      results[["Metabolite"]], metabolite.class[["MET_CHEM_NO"]]
    )]
    results <- results %>%
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
    results <- results %>%
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
    results <- results %>%
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
    results <- results %>%
      full_join(metabolite.class, by = "Metabolite") %>%
      rename(c(
        "MetaboliteClass" = "Ontology_Class",
        "lipidClass" = "Ontology_Subclass",
        "MetaboliteName" = "Metabolite_Name"
      ))
  }

  ## subset chain lengths

  if (data.type == "Metabolon") {
    ## subset chain lengths
    results <- results %>%
      dplyr::filter(adj.P.Val < p.value.cutoff, logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)) %>%
      dplyr::filter(MetaboliteClass == "Complex lipids") %>%
      mutate(newMet = MetaboliteName) %>%
      mutate(newMet = gsub("O-|P-", "", newMet)) %>%
      separate(newMet, c("lipid.class", "fatty.acid"), "[()]|-") %>%
      mutate(fatty.acid = gsub(".*/", "", fatty.acid)) %>%
      mutate(lipid.class = case_when(
        str_detect(lipid.class, "TAG") ~ "TAG",
        TRUE ~ lipid.class
      )) %>%
      mutate(fatty.acid = gsub("FA", "", fatty.acid)) %>%
      select(contrast, logFC, adj.P.Val, lipid.class, fatty.acid) %>%
      group_by(contrast, lipid.class, fatty.acid) %>%
      summarise_all(funs(mean)) %>%
      arrange(fatty.acid) %>%
      mutate(new.fatty.acid = fatty.acid) %>%
      separate(new.fatty.acid, c("chain.length", "saturation"), ":") %>%
      mutate(chain.length = as.numeric(as.character(chain.length))) %>%
      mutate(saturation.class = ifelse(saturation == "0", "saturated",
        ifelse(saturation == "1", "mono-unsaturated",
          "poly-unsaturated"
        )
      ))
  } else {
    results <- results %>%
      dplyr::filter(adj.P.Val < p.value.cutoff, logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)) %>%
      dplyr::filter(grepl("Complex lipids", MetaboliteClass)) %>%
      mutate(newMet = MetaboliteName) %>%
      mutate(newMet = gsub("O-|P-", "", newMet)) %>%
      separate(newMet, c("lipid.class", "fatty.acid"), "[()]|-") %>%
      mutate(fatty.acid = gsub(".*/", "", fatty.acid)) %>%
      mutate(lipid.class = case_when(
        str_detect(lipid.class, "TAG") ~ "TAG",
        TRUE ~ lipid.class
      )) %>%
      mutate(fatty.acid = gsub("FA", "", fatty.acid)) %>%
      select(contrast, logFC, adj.P.Val, lipid.class, fatty.acid) %>%
      group_by(contrast, lipid.class, fatty.acid) %>%
      summarise_all(funs(mean)) %>%
      arrange(fatty.acid) %>%
      mutate(new.fatty.acid = fatty.acid) %>%
      separate(new.fatty.acid, c("chain.length", "saturation"), ":") %>%
      mutate(chain.length = as.numeric(as.character(chain.length))) %>%
      mutate(saturation.class = ifelse(saturation == "0", "saturated",
        ifelse(saturation == "1", "mono-unsaturated",
          "poly-unsaturated"
        )
      ))
  }

  groups <- unique(results$contrast)
  for (i in groups) {
    p <- results %>%
      dplyr::filter(contrast == i) %>%
      ggplot(aes(
        x = lipid.class,
        y = fatty.acid,
        color = log2(logFC),
        size = -log10(adj.P.Val),
        shape = saturation.class
      )) +
      geom_point() +
      theme_bw() +
      theme(
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1
        ),
        axis.text = element_text(
          size = 11,
          # face = "bold",
          colour = "black"
        ),
        axis.title = element_text(size = 12, face = "bold")
      ) +
      scale_colour_gradientn(
        colours = rev(RColorBrewer::brewer.pal(10, "RdYlBu")),
        limits = c(-2, 2),
        oob = scales::squish,
        name = "fold changes"
      ) +
      guides(colour = guide_colourbar(
        barwidth = unit(0.3, "cm"),
        ticks.colour = "black",
        frame.colour = "black"
      )) +
      labs(
        title = i,
        x = "",
        y = "fatty acid chain length",
        size = "-log10(P-value)",
        shape = "saturation"
      ) +
      theme(axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ))

    if (save == "pdf") {
      ## print
      print(p)
    }
    ## save plots
    if (save != "pdf") {
      sjPlot::save_plot(
        filename = paste(here(), "chainLengthDistribution", paste0(groups[i], ".", save), sep = "/"),
        fig = p,
        width = fig.width,
        height = fig.height,
        dpi = dpi
      )
    }
  }
  if (save == "pdf") {
    dev.off()
  }
}
