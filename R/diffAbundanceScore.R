#' diffAbundanceScore
#'
#' The differential abundance (DA) score captures the tendency for a pathway t
#' o have increased levels of metabolites, relative to a control group.
#' The score is calculating by first applying a non-parametric differential
#' abundance test (in this study, Benjamini-Hochberg corrected Mann-Whitney
#' U-tests) to all metabolites in a pathway. Then, after determining which
#' metabolites are significantly increased/decreased in abundance, the
#' differential abundance score is defined as:
#' $$
#' DAS = \frac{n(Metabolites_{up}) - n(Metabolites_{down})}{n(Metabolites_{up}) + n(Metabolites_{down})}
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
#' @param save either "pdf", "svg" or "png"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi dpi only applicable for png
#'
#' @import tidyverse
#' @import here
#' @import ggplot2
#' @import ggrepel
#' @import RColorBrewer
#' @import ggpubr
#' @import graphics
#' @import grDevices
#' @import sjPlot
#' @import KEGGREST
#' @import stats
#' @import reshape2

diffAbundanceScore <- function(species=species,
                               ref.path=NULL,
                               results = results,
                               p.value.cutoff = 0.05,
                               fold.changes.cutoff = 1.5,
                               common.mets =15,
                               save = "pdf",
                               fig.width = 12,
                               fig.height = 9,
                               dpi = 300) {
  enrichment.results <- data.frame()
  if(is.null(ref.path)) {
    ref.path = here()
  } else
    ref.path = ref.path

  if (save == "pdf"){
    pdf(paste(ref.path, "diffAbundenceScore.pdf", sep = "/"),
        onefile = TRUE)
  } else if (save != "pdf") {
    dir.create(paste(ref.path, "diffAbundenceScore", sep = "/"))
  }

  ## load annotation file
  chemicalMetadata <- system.file("inst/extdata/ref",
                                  "Chemical_annotations.csv",
                                  package="metapacR")

  ## define metabolite classes
  metabolite_class <- chemicalMetadata[, colnames(chemicalMetadata) %in% c("SUPER_PATHWAY", "CHEMICAL_NAME","KEGG")]

  ## load enriched data
  pathDat <- results %>%
    filter(adj.P.Val < p.value.cutoff)
  # metabolite ID to kegg ID
  # match ID
  matchColumnID <-
    match(pathDat$Metabolite, metabolite_class$CHEMICAL_NAME)
  # add names
  pathDat$keggID[pathDat$Metabolite %in% metabolite_class$CHEMICAL_NAME] <-
    metabolite_class$KEGG[matchColumnID]
  ## define direction
  pathDat$direction <- ifelse(log2(pathDat$logFC) > fold.changes.cutoff, "up",
                              ifelse(log2(pathDat$logFC) < -fold.changes.cutoff, "down",
                                     "nochange"))

  ## identify pathways associated with obtained metabolites

  keggTest <-keggLink("pathway",sort(unique(pathDat$keggID)))
  rownames <- names(keggTest)
  keggTest <- as.data.frame(gsub("path:map","",keggTest))
  rownames<- gsub("^.*?:","",rownames)
  keggTest <- cbind(rownames, keggTest)
  names(keggTest) <- c("keggID", "keggPath")

  ## identify reference pathways wrt species

  keggRef <- keggLink("pathway", species)
  # rownamesRef <- names(keggRef)  ## gene ID
  keggRef <- as.data.frame(gsub(paste0("path:",species),"",keggRef))
  # rownamesRef <- gsub("^.*?:","",rownamesRef) ## gene ID
  # keggRef <- cbind(rownamesRef, keggRef) ## gene ID
  names(keggRef) <- c("keggPath")

  ## Merge two datasets

  keggDf <- keggTest[keggTest$keggPath %in% keggRef$keggPath, ]
  keggDf$keggPath <- paste0(species, keggDf$keggPath)

  ## select pathway dataset
  pathDatDas <- keggDf %>%
    full_join(pathDat, by="keggID") %>%
    group_by(contrast, keggPath) %>%
    mutate(direction = ifelse(log2(logFC) > fold.changes.cutoff, "positive",
                              ifelse(log2(logFC) < -fold.changes.cutoff, "negative",
                                     "nochange"))) %>%
    filter(direction != "nochange") %>%
    select(contrast, keggPath, direction) %>%
    drop_na() %>%
    mutate(count = n()) %>%
    dcast(contrast + keggPath ~ direction, value.var="count") %>%
    replace(is.na(.),0) %>%
    ##  more than 15 metabolites in common
    filter((positive + negative) > common.mets) %>%
    mutate(das=(positive - negative) /
             (positive + negative))

  ## create pathway annotations

  ## empty dataframe
  annotation <- data.frame()
  for (i in seq_along(unique(pathDatDas$keggPath))) {
    ## select kegg pathway name
    keggPath <- keggGet(pathDatDas$keggPath[i])
    ## add name to annotations
    annotation[i,"keggPath"] <- unname(keggPath[[1]]$ENTRY)
    annotation[i,"keggName"] <- unname(gsub("\\ - .*","",keggPath[[1]]$NAME))
    ## create class
    if (is.null(keggPath[[1]]$CLASS)) {
      annotation[i,"class"] <- unname(gsub("\\ - .*","",keggPath[[1]]$NAME))
    } else
      annotation[i,"class"] <- unname(gsub("^.*?; ","",keggPath[[1]]$CLASS))

  }
  ## join annotation data-frame with pathway data
  pathDatDas <-pathDatDas %>%
    full_join(annotation, by="keggPath")

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
  groups <- unique(pathDatDas$contrast)
  for (j in seq_along(groups)) {
    ## filtered dataset
    datFiltered <- pathDatDas[pathDatDas$contrast %in% groups[j],]  %>%
      arrange(keggName)
    ## plot
    p <- datFiltered %>%
      drop_na(keggName) %>%
      ggplot(aes(x=keggName, y=das)) +
      geom_segment(
        aes(xend=keggName, yend=0, fill= das),
        width=0.15, alpha = 0.5) +
      geom_point(aes(size = (positive + negative),
                     fill = das), shape=21) +
      geom_hline(yintercept=0, lty=1) +
      labs(fill= "DAS",
           size= "Pathway size",
           x= "",
           y= "Differential Abundance score",
           title = paste0(groups[j])) +
      # geom_text(aes(label=subclass,
      #               fontface= "bold"),
      #           size=3,
      #           hjust=-1) +
      theme_bw() +
      theme(
        panel.border = element_rect(colour = "black",
                                    fill=NA,
                                    size=1),
        axis.text = element_text(
          size = 11,
          #face = "bold",
          colour = "black"
        ),
        axis.title = element_text(size = 12, face = "bold")
      ) +
      scale_fill_gradientn(colours = rev(brewer.pal(10, "RdYlBu")),
                           limits = c(-2,2),
                           oob = scales::squish,
                           name = 'fold changes') +
      guides(fill = guide_colourbar(barwidth = unit(0.3, "cm"),
                                    ticks.colour = "black",
                                    frame.colour = "black")) +
      theme(axis.text.x = element_text(size=8)) +
      coord_flip()
    if (save == "pdf") {
      print(p)
    }
    ## save results
    if (save != "pdf") {
      save_plot(
        paste(ref.path, "diffAbundenceScore",paste0("das_",groups[j],".",save), sep = "/"),
        fig = p,
        width = 22,
        height = 12,
        dpi = 300
      )
    }

  }
  if (save == "pdf") {
    dev.off()
  }

}

}
