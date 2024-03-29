#' @title enrichmentScore.KEGG
#'
#' @description enrichmentScore calualtion based on enrichement pathway analysis.
#'
#' @param species species to use "hsa" or "mmu"
#' @param ref.path saving path
#' @param results fold changes results
#' @param p.value.cutoff cutoff value of p.value
#' @param fold.changes.cutoff higher cutoff value of fold changes
#' @param save either "pdf", "svg", "png" or "none"
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
#' @import graphics
#' @import grDevices
#' @importFrom KEGGREST keggLink keggGet
#' @import stats
#'
#' @return data.frame onject with enrichement.results and plots as save object in defined path.
#'
#' @export

enrichmentScore.KEGG <- function(species = c("hsa", "mmu"),
                                 ref.path = NULL,
                                 results,
                                 p.value.cutoff = 0.05,
                                 fold.changes.cutoff = 1.5,
                                 save = c("pdf", "svg", "png", "none"),
                                 data.type = c("MH", "Metabolon", "Others"),
                                 fig.width = 12,
                                 fig.height = 9,
                                 dpi = 300,
                                 Other_metadata = NULL) {

  options(warn = -1) ## supress all warning

  stopifnot(inherits(results, "data.frame"))
  validObject(results)

  species <- match.arg(species, c("hsa", "mmu"))
  save <- match.arg(save, c("pdf", "svg", "png", "none"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))


  enrichment.results <- data.frame()
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
  }

  if (file.exists(paste(ref.path, "keggDB.rds", sep = "/")) == FALSE) {
    ## reference kegg dataset building
    keggReferences <- KEGGREST::keggLink("pathway", species)
    referencesPathway <- unique(keggReferences[1:length(keggReferences)])
    referencesPathway <- gsub("path:", "", referencesPathway)

    ## empty data-frame
    keggReferenceDB <- data.frame()

    ## reference KEGG List

    for (rp in seq_along(referencesPathway)) {
      ## selected pathway
      try({
        pathway <- referencesPathway[rp]
        pathwayName <- KEGGREST::keggGet(pathway)[[1]][["NAME"]]

        ## list reference compounds
        listReferenceComp <- suppressWarnings(KEGGREST::keggGet(pathway)[[1]][["COMPOUND"]])
        compoundID <- names(listReferenceComp)
        nCompound <- length(compoundID)

        ## build dataset
        keggReferenceDB[rp, "pathway"] <- pathway
        keggReferenceDB[rp, "pathwayName"] <- pathwayName
        keggReferenceDB[rp, "compoundID"] <- paste(compoundID, collapse = ",")
        keggReferenceDB[rp, "nCompound"] <- nCompound
      })
    }

    ## filter non-NA compounds
    keggReferenceDB <- keggReferenceDB[keggReferenceDB$nCompound != 0, ]

    ## save reference kegg dataset
    keggReferenceDB <- saveRDS(
      keggReferenceDB,
      paste(ref.path, "keggDB.rds", sep = "/")
    )
  } else {
    ## read reference kegg dataset
    keggReferenceDB <- readRDS(paste(ref.path, "keggDB.rds", sep = "/"))
  }

  ## tidy reference list
  keggReferenceDBList <- keggReferenceDB %>%
    dplyr::select(pathway, pathwayName, compoundID) %>%
    separate_rows(compoundID)

  ## length of compounds
  nCompoundkeggReferenceDBList <- length(unique(keggReferenceDBList$compoundID))


  ## duplicate metabolite id listing in test dataset
  metaboliteOccurances <- data.frame(table(
    metabolite_class$KEGG
  ))

  colnames(metaboliteOccurances) <- c("keggID", "Freq")

  dupliMetaboliteOccurances <- metaboliteOccurances[
    metaboliteOccurances$Freq > 1,
  ]


  ## number of kegg listed metabolites

  nKeggIdList <- length(na.omit(metabolite_class$KEGG))
  dupliMetaboliteOccurances$Percent <- dupliMetaboliteOccurances$Freq /
    nKeggIdList * 100

  ##  balance reference dataset for duplicate IDs in test dataset

  for (jj in 1:nrow(keggReferenceDB)) {
    (keggReferenceDB[jj, "nCompoundAdjusted"] <- keggReferenceDB[jj, "nCompound"])

    for (iii in 1:nrow(dupliMetaboliteOccurances)) {
      keggID <- dupliMetaboliteOccurances[iii, "keggID"]

      if (keggID %in% unlist(strsplit(keggReferenceDB[jj, "compoundID"], ","))) {
        count <- length(unlist(
          strsplit(
            keggReferenceDB[jj, "compoundID"], ","
          )
        )) +
          ((nCompoundkeggReferenceDBList * dupliMetaboliteOccurances[iii, "Percent"]) /
             100)
        keggReferenceDB[jj, "nCompoundAdjusted"] <- as.numeric(round(count))
      }
    }
  }
  ## balance unique reference dataset count
  nCompoundkeggReferenceDBListAdjusted <- nCompoundkeggReferenceDBList
  for (ii in 1:nrow(dupliMetaboliteOccurances)) {
    nCompoundkeggReferenceDBListAdjusted <- round(nCompoundkeggReferenceDBListAdjusted + ((
      nCompoundkeggReferenceDBListAdjusted * dupliMetaboliteOccurances[ii, "Percent"]) /
        100))
  }

  ## load enriched data
  enrichDat <- results %>%
    dplyr::filter(adj.P.Val < p.value.cutoff) %>%
    full_join(metabolite_class, by = c("Metabolite" = "MET_CHEM_NO")) %>%
    dplyr::rename("keggID" ="KEGG")

  ## define direction
  enrichDat$direction <- ifelse(log2(enrichDat$logFC) > fold.changes.cutoff, "up",
                                ifelse(log2(enrichDat$logFC) < -fold.changes.cutoff, "down",
                                       "nochange"
                                )
  )

  ## define groups
  groups <- na.omit(unique(enrichDat$contrast))

  ## enrichement pathway

  for (i in seq_along(groups)) {
    ## create empty data-frame
    enrichTableUp <- data.frame()
    enrichTableDown <- data.frame()
    ## filtered data by group
    enrichDatFiltered <- enrichDat[enrichDat$contrast %in% groups[i], ]
    ## remove missing KeggIDs
    enrichDatFiltered <- enrichDatFiltered[!is.na(enrichDatFiltered$keggID), ]
    ## subset data by direction
    enrichDatFilteredUp <- enrichDatFiltered[enrichDatFiltered$direction == "up", ]
    enrichDatFilteredDown <- enrichDatFiltered[enrichDatFiltered$direction == "down", ]
    ## define reference data
    keggReferenceDBDirection <- keggReferenceDB
    ## get adjusted list
    keggReferenceDBDirection$nCompoundkeggReferenceDBListAdjusted <- nCompoundkeggReferenceDBListAdjusted
    ## select non-negative data
    if (nrow(enrichDatFilteredUp) != 0 && nrow(enrichDatFilteredDown) != 0) {
      ## list direction data
      dfList <- list(enrichDatFilteredUp, enrichDatFilteredDown)
      for (list in dfList) {
        ## calculate expected metabolite count
        listedMetabolite <- unique(list$keggID)
        ## length of metabolite
        keggReferenceDBDirection[["nListedMetabolite"]] <-
          length(listedMetabolite)
        ## create expected compound
        keggReferenceDBDirection <- keggReferenceDBDirection %>%
          mutate(expectedCompound = (nCompoundAdjusted * nListedMetabolite) / nCompoundkeggReferenceDBListAdjusted)
        ## create empty data frame
        keggReferenceDBDirection[["nCompoundFC"]] <- 0
        ## count number of metabolites
        for (j in seq_along(listedMetabolite)) {
          ## select metabolite
          metabolite <- listedMetabolite[j]
          ## create keggReferenceDBDirection
          for (k in 1:nrow(keggReferenceDBDirection)) {
            if (metabolite %in% unlist(strsplit(keggReferenceDBDirection[k, "compoundID"], ","))) {
              keggReferenceDBDirection[k, "nCompoundFC"] <- keggReferenceDBDirection[k, "nCompoundFC"] + 1
            }
          }
        }
        ## calculate enrichment
        for (l in 1:nrow(keggReferenceDBDirection)) {
          ## define matrix
          a <- keggReferenceDBDirection[l, "nCompoundFC"]
          b <- keggReferenceDBDirection[l, "nCompoundAdjusted"] - a
          c <- keggReferenceDBDirection[l, "nListedMetabolite"] - a
          d <- keggReferenceDBDirection[l, "nCompoundkeggReferenceDBListAdjusted"] - a - b - c
          ## create matrix
          fischerMatrix <- matrix(c(a, b, c, d), nrow = 2)
          ## perform fischer exact
          pValue <- fisher.test(fischerMatrix, alternative = "greater")$p.value
          ## adjust p-values
          pValueAdjusted <- p.adjust(pValue, method = "BH")
          ## calculate enrichment
          pathwayEnrichment <- keggReferenceDBDirection[l, "nCompoundFC"] /
            keggReferenceDBDirection[l, "expectedCompound"]
          ## save results as as up or down-regulated pathway
          if (length(na.omit(list$direction)) == 0) {
            next
          }
          if (unique(list$direction) == "up") {
            enrichTableUp[l, "pathway"] <- keggReferenceDBDirection[l, "pathway"]
            enrichTableUp[l, "pathwayName"] <- gsub("\\ - .*", "", keggReferenceDBDirection[l, "pathwayName"])
            enrichTableUp[l, "enrichment"] <- pathwayEnrichment
            enrichTableUp[l, "pValue"] <- pValue
            enrichTableUp[l, "adjPValueFDR"] <- pValueAdjusted
            enrichTableUp[l, "direction"] <- "up"
          } else if (unique(list$direction) == "down") {
            enrichTableDown[l, "pathway"] <- keggReferenceDBDirection[l, "pathway"]
            enrichTableDown[l, "pathwayName"] <- gsub("\\ - .*", "", keggReferenceDBDirection[l, "pathwayName"])
            enrichTableDown[l, "enrichment"] <- pathwayEnrichment
            enrichTableDown[l, "pValue"] <- pValue
            enrichTableDown[l, "adjPValueFDR"] <- pValueAdjusted
            enrichTableDown[l, "direction"] <- "down"
          }
        }
      }
    } else if (nrow(enrichDatFilteredUp) != 0) {
      dfList <- enrichDatFilteredUp
      ## calculate expected metabolite count
      listedMetabolite <- unique(dfList$keggID)
      ## length of metabolite
      keggReferenceDBDirection[["nListedMetabolite"]] <-
        length(listedMetabolite)
      ## create expected compound
      keggReferenceDBDirection <- keggReferenceDBDirection %>%
        mutate(expectedCompound = (nCompoundAdjusted * nListedMetabolite) / nCompoundkeggReferenceDBListAdjusted)
      ## create empty data frame
      keggReferenceDBDirection[["nCompoundFC"]] <- 0
      ## count number of metabolites
      for (j in seq_along(listedMetabolite)) {
        ## select metabolite
        metabolite <- listedMetabolite[j]
        ## create keggReferenceDBDirection
        for (k in 1:nrow(keggReferenceDBDirection)) {
          if (metabolite %in% unlist(strsplit(keggReferenceDBDirection[k, "compoundID"], ","))) {
            keggReferenceDBDirection[k, "nCompoundFC"] <- keggReferenceDBDirection[k, "nCompoundFC"] + 1
          }
        }
      }
      ## calculate enrichment
      for (l in 1:nrow(keggReferenceDBDirection)) {
        ## define matrix
        a <- keggReferenceDBDirection[l, "nCompoundFC"]
        b <- keggReferenceDBDirection[l, "nCompoundAdjusted"] - a
        c <- keggReferenceDBDirection[l, "nListedMetabolite"] - a
        d <- keggReferenceDBDirection[l, "nCompoundkeggReferenceDBListAdjusted"] - a - b - c
        ## create matrix
        fischerMatrix <- matrix(c(a, b, c, d), nrow = 2)
        ## perform fischer exact
        pValue <- fisher.test(fischerMatrix, alternative = "greater")$p.value
        ## adjust p-values
        pValueAdjusted <- p.adjust(pValue, method = "BH")
        ## calculate enrichment
        pathwayEnrichment <- keggReferenceDBDirection[l, "nCompoundFC"] /
          keggReferenceDBDirection[l, "expectedCompound"]
        ## save results as as up-regulated pathways
        enrichTableUp[l, "pathway"] <- keggReferenceDBDirection[l, "pathway"]
        enrichTableUp[l, "pathwayName"] <- gsub("\\ - .*", "", keggReferenceDBDirection[l, "pathwayName"])
        enrichTableUp[l, "enrichment"] <- pathwayEnrichment
        enrichTableUp[l, "pValue"] <- pValue
        enrichTableUp[l, "adjPValueFDR"] <- pValueAdjusted
        enrichTableUp[l, "direction"] <- "up"
      }
    } else if (nrow(enrichDatFilteredDown) != 0) {
      dfList <- enrichDatFilteredDown
      ## calculate expected metabolite count
      listedMetabolite <- unique(dfList$keggID)
      ## length of metabolite
      keggReferenceDBDirection[["nListedMetabolite"]] <-
        length(listedMetabolite)
      ## create expected compound
      keggReferenceDBDirection <- keggReferenceDBDirection %>%
        mutate(expectedCompound = (nCompoundAdjusted * nListedMetabolite) / nCompoundkeggReferenceDBListAdjusted)
      ## create empty data frame
      keggReferenceDBDirection[["nCompoundFC"]] <- 0
      ## count number of metabolites
      for (j in seq_along(listedMetabolite)) {
        ## define metabolite
        metabolite <- listedMetabolite[j]
        ## create keggReferenceDBDirection
        for (k in 1:nrow(keggReferenceDBDirection)) {
          if (metabolite %in% unlist(strsplit(keggReferenceDBDirection[k, "compoundID"], ","))) {
            keggReferenceDBDirection[k, "nCompoundFC"] <- keggReferenceDBDirection[k, "nCompoundFC"] + 1
          }
        }
      }
      ## calculation of enrichment
      for (l in 1:nrow(keggReferenceDBDirection)) {
        ## select metabolite
        a <- keggReferenceDBDirection[l, "nCompoundFC"]
        b <- keggReferenceDBDirection[l, "nCompoundAdjusted"] - a
        c <- keggReferenceDBDirection[l, "nListedMetabolite"] - a
        d <- keggReferenceDBDirection[l, "nCompoundkeggReferenceDBListAdjusted"] - a - b - c
        ## create matrix
        fischerMatrix <- matrix(c(a, b, c, d), nrow = 2)
        ## perform fischer exact
        pValue <- fisher.test(fischerMatrix, alternative = "greater")$p.value
        ## adjust p-values
        pValueAdjusted <- p.adjust(pValue, method = "BH")
        ## calculate enrichement
        pathwayEnrichment <- keggReferenceDBDirection[l, "nCompoundFC"] /
          keggReferenceDBDirection[l, "expectedCompound"]
        ## save results as as up-regulated pathways
        enrichTableDown[l, "pathway"] <- keggReferenceDBDirection[l, "pathway"]
        enrichTableDown[l, "pathwayName"] <- gsub("\\ - .*", "", keggReferenceDBDirection[l, "pathwayName"])
        enrichTableDown[l, "enrichment"] <- pathwayEnrichment
        enrichTableDown[l, "pValue"] <- pValue
        enrichTableDown[l, "adjPValueFDR"] <- pValueAdjusted
        enrichTableDown[l, "direction"] <- "down"
      }
    }
    ## merge data
    enrichTable <- rbind(enrichTableUp, enrichTableDown)

    if (nrow(enrichTable) == 0) {
      next
    } else if (nrow(enrichTable[enrichTable$pValue < p.value.cutoff, ]) != 0) {
      ## tidy enrichTable
      enrichTable <- enrichTable %>%
        dplyr::filter(pValue < p.value.cutoff) %>%
        mutate(enrichment = ifelse(direction == "down", -enrichment, enrichment))

      enrichTable$pathway <- paste0(
        enrichTable$pathway, "-",
        enrichTable$direction
      )
      enrichTable$contrast <- groups[i]
      ## save enrichTable
      enrichment.results <- bind_rows(enrichment.results, enrichTable)
      ## plot
      p <- enrichTable %>%
        mutate(pathway = gsub("\\-.*", "", pathway)) %>%
        dplyr::select(pathway, pathwayName, enrichment, direction, contrast) %>%
        group_by(pathway) %>%
        summarise(n = n(), across()) %>%
        mutate(enrichment = case_when(n == 2 ~ (enrichment[match("up", direction)] + enrichment[match("down", direction)]), TRUE ~ enrichment)) %>%
        dplyr::select(-direction) %>%
        distinct() %>%
        mutate(direction = case_when(
          enrichment > 0 ~ "up",
          enrichment < 0 ~ "down",
          enrichment == 0 ~ "nochange"
        )) %>%
        arrange(desc(enrichment)) %>%
        ungroup() %>%
        ggplot(aes(x = reorder(pathwayName, enrichment, FUN = sum), y = enrichment, group = direction, fill = direction)) +
        geom_bar(stat = "identity", color = "black", size = 0.25) +
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
        coord_flip() +
        theme(legend.position = "none") +
        xlab("KEGG pathways") +
        ylab("Enrichment score") +
        ggtitle(paste0("Enriched pathways:", groups[i])) +
        scale_fill_manual(values = c(up = "#e41a1c", down = "#377eb8"))
      if (save == "none") {
        ## print
        print(p)
      }
      ## save results
      if (save != "none") {
        ggsave(
          paste(ref.path, "/KEGG.enrichmentScore_", paste0("barplot_", groups[i], ".", save), sep = ""),
          plot = p,
          width = fig.width,
          height = fig.height,
          dpi = dpi
        )
      }
    }
  }

  return(results = enrichment.results)
  options(warn = 0) ## reset all warnings
}
