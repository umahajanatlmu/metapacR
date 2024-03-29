#' @title lipidChainLengthCorrelation
#'
#' @description compute correlation of chain length with abundance.
#'
#' @param species species to use "hsa" or "mmu"
#' @param results fold changes data
#' @param path saving path
#' @param save either "pdf", "svg" , "png" or "none"
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi dpi only applicable for png
#' @param Other_metadata dataframe with metadata....it must have  columns: Metabolite, Metabolite_Name, Ontology_Class, Ontology_Subclass
#'
#' @import tidyverse
#' @importFrom here here
#' @importFrom RColorBrewer brewer.pal
#' @import graphics
#' @import grDevices
#' @importFrom ggpmisc stat_fit_glance
#'
#' @return plot as save object in defined path.
#'
#' @export

lipidChainLengthCorrelation <- function(results,
                                        path = NULL,
                                        species = c("hsa", "mmu"),
                                        save = c("pdf", "svg", "png", "none"),
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

  ## load annotation file
  if (data.type == "Metabolon" && species %in% c("hsa","mmu")) {
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
    res <- results %>%
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
      separate(fatty.acid, c("chain.length", "saturation"), ":") %>%
      mutate(chain.length = as.numeric(as.character(chain.length))) %>%
      mutate(saturation.class = ifelse(saturation == "0", "saturated",
                                       ifelse(saturation == "1", "mono-unsaturated",
                                              "poly-unsaturated"
                                       )
      ))
  } else {
    res <- results %>%
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
      separate(fatty.acid, c("chain.length", "saturation"), ":") %>%
      mutate(chain.length = as.numeric(as.character(chain.length))) %>%
      mutate(saturation.class = ifelse(saturation == "0", "saturated",
                                       ifelse(saturation == "1", "mono-unsaturated",
                                              "poly-unsaturated"
                                       )
      ))
  }

  ## plot
  groups <- unique(res$contrast)

  for (i in groups) {
    res.subset <- res %>%
      dplyr::filter(contrast == i)

    p <- ggplot(
      res.subset,
      aes(
        x = chain.length,
        y = t.ratio,
        color = saturation.class
      )
    ) +
      geom_hline(yintercept = 2, linetype = "dotted") +
      geom_hline(yintercept = -2, linetype = "dotted") +
      geom_point(alpha = 0.7) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        show.legend = FALSE
      ) +
      ggpmisc::stat_fit_glance(
        method = "lm",
        label.x = "left",
        label.y = "bottom",
        method.args = list(formula = y ~ x),
        aes(label = sprintf(
          'R^2~"="~%.3f~~italic(p)~"="~%.2f',
          after_stat(..r.squared..),
          stat(..p.value..)
        )),
        parse = TRUE,
        size = 4
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
      theme(axis.ticks.x = element_blank()) +
      ggtitle(i) +
      xlab("Fatty acid chain length") +
      ylab("t-ratio") +
      facet_wrap(~lipid.class, ncol = 7) +
      scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set1")) +
      theme(legend.position = "bottom") +
      labs(color = "Saturation status")

    ## save plots
    if (save != "none") {
      ggsave(
        filename = paste(
          path,
          paste0(i, "_chainLengthCorrelation", ".", save),
          sep = "/"
        ),
        plot = p,
        width = fig.width,
        height = fig.height,
        dpi = dpi
      )
    } else if (save == "none") {
      print(p)
    }
  }

}
