#' @title diffAbundanceScore
#'
#' @description The differential abundance (DA) score captures the tendency for a pathway t
#' o have increased levels of metabolites, relative to a control group.
#' The score is calculating by first applying a non-parametric differential
#' abundance test (in this study, Benjamini-Hochberg corrected Mann-Whitney
#' U-tests) to all metabolites in a pathway. Then, after determining which
#' metabolites are significantly increased/decreased in abundance, the
#' differential abundance score is defined as:
#' $$
#' DAS = {n(Metabolites_{up}) - n(Metabolites_{down})}{n(Metabolites_{up}) + n(Metabolites_{down})}
#' $$
#' Thus, the DA score varies from -1 to 1. A score of -1 indicates that all
#' metabolites in a pathway decreased in abundance, while a score of 1
#' indicates that all metabolites increased.
#'
#' @param species species to use "hsa" or "mmu"
#' @param ref.path saving path
#' @param results fold changes results
#' @param p.value.cutoff cutoff value of p.value
#' @param fold.changes.cutoff higher cutoff value of fold changes
#' @param common.mets minimum number of common metaolites
#' @param save either "pdf", "svg", "png", "none"
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi dpi only applicable for png
#' @param Other_metadata dataframe with metadata....it must have  columns: Metabolite, Metabolite_Name, Ontology_Class, Ontology_Subclass.
#, KEGG
#'
#' @import tidyverse
#' @importFrom here here
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales rescale
#' @import graphics
#' @import grDevices
#' @importFrom KEGGREST keggGet keggLink
#' @import stats
#' @importFrom reshape2 dcast
#' @importFrom scales squish
#'
#' @return output in defined path
#'
#' @export

diffAbundanceScore <- function(species = c("hsa", "mmu"),
                               ref.path = NULL,
                               results,
                               p.value.cutoff = 0.05,
                               fold.changes.cutoff = 1.5,
                               common.mets = 15,
                               save = c("pdf", "svg", "png", "none"),
                               data.type = c("MH", "Metabolon", "Others"),
                               fig.width = 12,
                               fig.height = 9,
                               dpi = 300,
                               Other_metadata = NULL) {
  stopifnot(inherits(results, "data.frame"))
  validObject(results)

  species <- match.arg(species,c("hsa", "mmu"))
  save <- match.arg(save,c("pdf", "svg", "png", "none"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  if (is.null(ref.path)) {
    ref.path <- here::here()
    ifelse(!dir.exists(file.path(paste0(ref.path), "results")),
           dir.create(file.path(paste0(ref.path), "results")),
           FALSE
    )
    path <- paste(ref.path, "results", sep = "/")
  } else {
    ref.path <- ref.path
  }

  if (data.type == "Others") {
    stopifnot(inherits(Other_metadata, "data.frame"))
    validObject(Other_metadata)
  }

  ## load annotation file
  if (data.type == "Metabolon" && species %in% c("hsa", "mmu")) {
    data("chemicalMetadata")
    chemicalMetadata <- force(chemicalMetadata)

    chemicalMetadata <- chemicalMetadata %>%
      mutate(across(everything(), as.character))

    ## define metabolite classes
    columnToSelect <- c("SUPER_PATHWAY", "CHEMICAL_NAME", "KEGG", "MET_CHEM_NO")
    metabolite_class <- chemicalMetadata %>%
      dplyr::select(any_of(columnToSelect)) %>%
      separate_rows(KEGG, sep=",") %>%
      mutate(KEGG=trimws(KEGG))

    ## load enriched data
    pathDat <- results %>%
      dplyr::filter(adj.P.Val < p.value.cutoff) %>%
      left_join(metabolite_class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      dplyr::rename(keggID =KEGG) %>%
      drop_na(keggID)
  }

  if (data.type == "MH" && species == "hsa") {
    data("chemicalMetadata_MH")
    chemicalMetadata <- force(chemicalMetadata_MH)

    chemicalMetadata <- chemicalMetadata %>%
      mutate(across(everything(), as.character))

    ## define metabolite classes
    columnToSelect <- c("ONTOLOGY1_NAME", "METABOLITE_NAME", "KEGG", "MET_CHEM_NO")
    metabolite_class <- chemicalMetadata %>%
      dplyr::select(any_of(columnToSelect)) %>%
      separate_rows(KEGG, sep=",") %>%
      mutate(KEGG=trimws(KEGG))

    ## load enriched data
    pathDat <- results %>%
      dplyr::filter(adj.P.Val < p.value.cutoff) %>%
      left_join(metabolite_class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      dplyr::rename(keggID =KEGG) %>%
      drop_na(keggID)
  }

  if (data.type == "MH" && species == "mmu") {
    data("chemicalMetadata_MH_mmu")
    chemicalMetadata <- force(chemicalMetadata_MH_mmu)

    chemicalMetadata <- chemicalMetadata %>%
      mutate(across(everything(), as.character))

    ## define metabolite classes
    columnToSelect <- c("ONTOLOGY1_NAME", "METABOLITE_NAME", "KEGG_ID", "MET_CHEM_NO")
    metabolite_class <- chemicalMetadata %>%
      dplyr::select(any_of(columnToSelect)) %>%
      mutate(KEGG = KEGG_ID) %>%
      dplyr::select(-KEGG_ID) %>%
      separate_rows(KEGG, sep=",") %>%
      mutate(KEGG=trimws(KEGG))

    ## load enriched data
    pathDat <- results %>%
      dplyr::filter(adj.P.Val < p.value.cutoff) %>%
      left_join(metabolite_class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      dplyr::rename(keggID =KEGG) %>%
      drop_na(keggID)
  }

  if (data.type == "Others") {
    chemicalMetadata <- Other_metadata

    chemicalMetadata <- chemicalMetadata %>%
      mutate(across(everything(), as.character))

    columnToSelect <- c("Metabolite", "Metabolite_Name", "Ontology_Class", "Ontology_Subclass", "KEGG")
    metabolite_class <- chemicalMetadata %>%
      dplyr::select(any_of(columnToSelect)) %>%
      rename(c(
        "MET_CHEM_NO" = "Metabolite",
        "METABOLITE_NAME" = "Metabolite_Name",
        "ONTOLOGY1_NAME" = "Ontology_Class"
      )) %>%
      separate_rows(KEGG, sep=",") %>%
      mutate(KEGG=trimws(KEGG))

    ## load enriched data
    pathDat <- results %>%
      dplyr::filter(adj.P.Val < p.value.cutoff) %>%
      left_join(metabolite_class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
      dplyr::rename(keggID =KEGG) %>%
      drop_na(keggID)
  }

  keggTest <- KEGGREST::keggLink("pathway", unique(pathDat$keggID))
  keggTest <- data.frame(keggID=names(keggTest), keggPath=keggTest, row.names = NULL)
  keggTest$keggID <- gsub("^.*?:", "", keggTest$keggID)
  keggTest$keggPath <- gsub("path:map", "", keggTest$keggPath)

  # identify reference pathways wrt species
  keggRef <- KEGGREST::keggLink("pathway", species)
  keggRef <- data.frame(geneID=names(keggRef), keggPath=keggRef, row.names = NULL)
  # rownamesRef <- names(keggRef)  ## gene ID
  keggRef$keggPath <- gsub(paste0("path:", species), "", keggRef$keggPath)

  ## select pathway dataset
  pathDatDas <- keggRef %>%
    inner_join(keggTest, by ="keggPath") %>%
    inner_join(pathDat, by = "keggID") %>%
    group_by(contrast, keggPath) %>%
    mutate(direction = ifelse(log2(logFC) > fold.changes.cutoff, "positive",
                              ifelse(log2(logFC) < -fold.changes.cutoff, "negative",
                                     "nochange"
                              )
    )) %>%
    dplyr::filter(direction != "nochange") %>%
    dplyr::select(contrast, keggPath, direction) %>%
    drop_na() %>%
    mutate(count = n()) %>%
    reshape2::dcast(contrast + keggPath ~ direction, value.var = "count") %>%
    replace(is.na(.), 0) %>%
    ##  more than 15 metabolites in common
    dplyr::filter((positive + negative) > common.mets) %>%
    mutate(das = (positive - negative) /
             (positive + negative)) %>%
    ungroup() %>%
    mutate(keggPath = paste0("hsa", keggPath))

  ## create pathway annotations

  ## empty dataframe
  annotation <- data.frame()
  for (i in seq_along(unique(pathDatDas$keggPath))) {
    ## select kegg pathway name
    keggPath <- KEGGREST::keggGet(pathDatDas$keggPath[i])
    ## add name to annotations
    annotation[i, "keggPath"] <- unname(keggPath[[1]]$ENTRY)
    annotation[i, "keggName"] <- unname(gsub("\\ - .*", "", keggPath[[1]]$NAME))
    ## create class
    if (is.null(keggPath[[1]]$CLASS)) {
      annotation[i, "class"] <- unname(gsub("\\ - .*", "", keggPath[[1]]$NAME))
    } else {
      annotation[i, "class"] <- unname(gsub("^.*?; ", "", keggPath[[1]]$CLASS))
    }
  }
  ## join annotation data-frame with pathway data
  pathDatDas <- pathDatDas %>%
    inner_join(annotation, by = "keggPath")

  # ## add ceramides and sphigomylins subclass
  # subclass <- keggTest %>%
  #   mutate(subclass = ifelse(keggID %in% "C00195", "Cer",
  #                            ifelse(keggID %in% "C00550",
  #                                   "SMs", NA))) %>%
  #   mutate(keggPath = paste0(species, keggPath)) %>%
  #   select(keggPath, subclass) %>%
  #   unique()
  # ##----------------------------------------------------------------
  # ##       join class/subclass data-frame with pathway data       --
  # ##----------------------------------------------------------------
  # pathDatDas <-pathDatDas %>%
  #   full_join(subclass, by="keggPath") %>%
  #   mutate(colors = ifelse(class == "Lipid metabolism", "#e41a1c",
  #                          ifelse(class == "Amino acid metabolism", "#377eb8", "gray"))) %>%
  #   drop_na(.,any_of(c("Subtype","keggName")))

  ## differential abundance score

  ## select groups
  groups <- na.omit(unique(pathDatDas$contrast))

  for (j in seq_along(groups)) {
    ## filtered dataset
    datFiltered <- pathDatDas[pathDatDas$contrast %in% groups[j], ] %>%
      arrange(keggName) %>%
      mutate(size = (positive + negative),
             size = scales::rescale(size, to = c(0, 100)))
    ## plot
    p <- datFiltered %>%
      drop_na(keggName) %>%
      ggplot(aes(x = keggName, y = das)) +
      geom_segment(
        aes(xend = keggName, yend = 0, fill = das),
        width = 0.15, alpha = 0.5
      ) +
      geom_point(aes(
        size = size,
        fill = das
      ), shape = 21) +
      geom_hline(yintercept = 0, lty = 1) +
      labs(
        fill = "DAS",
        size = "Pathway size",
        x = "",
        y = "Differential Abundance score",
        title = paste0(groups[j])
      ) +
      # geom_text(aes(label=subclass,
      #               fontface= "bold"),
      #           size=3,
      #           hjust=-1) +
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
      scale_fill_gradientn(
        colours = rev(RColorBrewer::brewer.pal(10, "RdYlBu")),
        limits = c(-2, 2),
        oob = scales::squish,
        name = "fold changes"
      ) +
      guides(fill = guide_colourbar(
        barwidth = unit(0.3, "cm"),
        ticks.colour = "black",
        frame.colour = "black"
      )) +
      theme(axis.text.x = element_text(size = 8)) +
      coord_flip()

    if (save == "none") {
      print(p)
    }
    ## save results
    if (save != "none") {
      ggsave(
        paste(ref.path, "/diffAbundenceScore_", paste0("das_", groups[j], ".", save), sep = ""),
        plot = p,
        width = 22,
        height = 12,
        dpi = 300
      )
    }
  }
}

