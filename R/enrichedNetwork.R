#' enrichedNetwork
#'
#' @param species species to use "hsa" or "mmu"
#' @param ref.path saving path
#' @param results fold changes results
#' @param p.value.cutoff cutoff value of p.value
#' @param fold.changes.cutoff higher cutoff value of fold changes
#' @param network.method "diffusion" as method of choice, for other method refer to FELLA.
#' @param legend.pathway number of pathways to map on plot
#' @param save either "pdf", "svg" or "png"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param ... plot extensions
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
#' @import FELLA
#' @import igraph
#' @import scales
enrichedNetwork <- function(species=species,
                            ref.path=NULL,
                            results = results,
                            p.value.cutoff = 0.05,
                            fold.changes.cutoff = 1.5,
                            network.method = "diffusion",
                            legend.pathway = 3,
                            save = "pdf",
                            fig.width = 12,
                            fig.height = 9,
                            ...) {

  enrichment.results <- data.frame()

  plot.pathway.impact <- list()

  if(is.null(ref.path)) {
    ref.path = here()
  } else
    ref.path = ref.path

  if (save == "pdf"){
    pdf(paste(ref.path, "enrichmentNetwork.pdf", sep = "/"), paper = "a4r",
        onefile = TRUE)
  } else if (save != "pdf") {
    dir.create(paste(ref.path, "enrichmentNetwork", sep = "/"))
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

  ## create pathway fella data
  pathDatFella <- keggDf %>%
    full_join(pathDat, by = "keggID", all = FALSE) %>%
    mutate(keggPath = gsub(species, "", keggPath),
           keggID = gsub(" ", "", keggID)) %>%
    drop_na(contrast)


  ## filter overview pathways
  graph <- buildGraphFromKEGGREST(organism = species,
                                  filter.path = NULL)

  ## Cannot be overwritten
  tmpdir <- paste0(tempdir())
  unlink(tmpdir, recursive = TRUE)

  ## build data from graph
  buildDataFromGraph(
    keggdata.graph = graph,
    databaseDir = tmpdir,
    internalDir = FALSE,
    matrices = network.method,
    normality = "diffusion",
    niter = 250
  )
  if (species == "hsa") {
    alias2entrez <- as.list(org.Hs.egSYMBOL2EG)
  } else if (species == "mmu") {
    alias2entrez <- as.list(org.Mm.egSYMBOL2EG)
  }

  ## get associated pathways
  entrez2ec <- keggLink("enzyme", species)
  entrez2path <- keggLink("pathway", species)

  ## perform FELLA
  fellaData <- loadKEGGdata(databaseDir = tmpdir,
                            internalDir = FALSE,
                            loadMatrix = network.method)

  ## get FELLA community

  ## compound
  idCpd <- getCom(fellaData, level = 5,
                  format = "id") %>% names
  ## reaction
  idRx <- getCom(fellaData, level = 4,
                 format = "id") %>% names
  ## enzymes
  idEc <- getCom(fellaData, level = 3,
                 format = "id") %>% names

  groups <- unique(pathDatFella$contrast)

  for (i in seq_along(groups)) {
    ## filtered data
    pathDatFellaFiltered <- pathDatFella[pathDatFella$contrast %in% groups[i], ]
    ## define compunds
    cpd <- unique(pathDatFellaFiltered$keggID)
    ## perform enricment analysis
    analysis <- FELLA::enrich(
      compounds = cpd,
      data = fellaData,
      method = network.method,
      approx = "normality"
    )

    ## create fellaplot function

    ## fella plot
    # net <- FELLA::plot(
    #   analysis,
    #   method = "diffusion",
    #   data = fellaData,
    #   nlimit = 250,
    #   plotLegend = TRUE,
    #   vertex.label.cex = 1,
    #   vertex.label.degree = pi,
    #   rescale = TRUE
    # )

    ## generateResultsGraph

    g <- generateResultsGraph(object = analysis,
                              data = fellaData,
                              method = network.method)
    ## define undirected
    unionGraphUndir <- as.undirected(g, mode = "collapse")
    ## create igraph obj
    clp <- cluster_edge_betweenness(unionGraphUndir)

    ## obtain matrix properties
    hubscore <- hub.score(g)$vector
    authscore <- authority.score(g)$vector
    eigenvalue <- eigen_centrality(g)$vector
    graph.strength <- graph.strength(g)
    centrality <- degree(g) # degree centality
    ## generate matrix table
    gDf <- as.data.frame(
      list(
        Hubscore = hubscore,
        Authscore = authscore,
        Eigen = eigenvalue,
        strength = graph.strength,
        centrality = centrality
      ),
      stringsAsFactors = FALSE
    )
    ## define rownames
    gDf$keggPath <- row.names(gDf)
    ## add p scores statistics
    pscores <- getPscores(object = analysis,
                          method = network.method)
    ## generateResultsTable
    table <- generateResultsTable(object = analysis,
                                  data = fellaData,
                                  method = network.method)
    ## p.adjust holm
    table$holm <- p.adjust(table$p.score, method = "holm")
    ## p.adjust fdr
    table$FDR <- p.adjust(table$p.score, method = "fdr")
    ## merge tables
    colnames(table)[1] <- "keggPath"
    table <- merge(table, gDf, by = "keggPath", all = FALSE)
    ## pathway impact
    table.plot <- subset(table, Entry.type %in% c("pathway", "reaction", "module" ))
    table.plot$KEGG.name <- gsub(" -.*", "", table.plot$KEGG.name)

    ## lattice of annotated networks
    clust <- data.frame(cl = clp$membership)
    rownames(clust) <- names(V(unionGraphUndir))
    clust$desc <- table.plot$KEGG.name[match(rownames(clust), table.plot$keggPath)]

    ## coalesce cluster names
    clust <- clust %>%
      drop_na() %>%
      group_by(cl) %>%
      summarise(desc = paste(desc, collapse=":"))

    n <- legend.pathway
    pat <- paste0('^([^:]+(?::[^:]+){',n-1,'}).*')
    clust$desc <- sub(pat, '\\1', clust$desc)
    clust$desc <- gsub(":", "\n", clust$desc)

    colors <- colorRampPalette(brewer.pal(9,"Set1"))(length(clust$cl))
    colors <- paste0(colors, "50")
    Group <- gl(n= length(clust$cl), 1,
                labels = clust$desc[clust$cl])
    mark.col <- colors[Group]

    if (save == "svg") {
      svg(
        paste(ref.path, "enrichmentNetwork",paste0("enrichmentNetwork_",groups[i],".",save), sep = "/"),
        width = fig.width,
        height = fig.height
      )
    } else if (save == "png") {
      png(
        paste(ref.path, "enrichmentNetwork",paste0("enrichmentNetwork_",groups[i],".",save), sep = "/"),
        width = fig.width,
        height = fig.height
      )
    }

    ## plot grid
    par(mfrow=c(1,2))

    ## plot
    plot(clp,
         unionGraphUndir,
         alpha = 0.5,
         mark.border="black",
         vertex.size = (V(unionGraphUndir)$input + 0.75) * 5,
         vertex.label = NA,
         mark.col = mark.col,
         main = groups[i]
    )
    plot(NA)
    legend("center",
           bty = "n",
           cex=0.7,
           legend=levels(Group),
           fill=colors)
    ## print
    dev.off()

    ## plot pathway impact
    ## plot
    p <- ggplot(table.plot,
                aes(x = centrality,
                    y = -log10(p.score))) +
      geom_point(aes(size = centrality,
                     fill = -log10(p.score)),
                 color = "black",
                 pch = 21) +
      ggtitle(paste0("Pathway impact:", groups[i])) +
      theme(legend.text.align = 0) +
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
      )  +
      scale_fill_gradientn(colours = rev(brewer.pal(10, "RdYlBu")),
                           limits = c(0,8),
                           oob = scales::squish,
                           name = 'fold changes') +
      guides(fill = guide_colourbar(barwidth = unit(0.3, "cm"),
                                    ticks.colour = "black",
                                    frame.colour = "black")) +
      labs(x = "Pathway Impact",
           y = "p value (-log10)",
           fill = "p value",
           size = "Pathway size") +
      geom_text_repel(data = subset(table.plot,
                                    centrality >= quantile(centrality, 0.75)[[1]]),
                      aes(label = KEGG.name))

    ## print plot
    plot.pathway.impact[[groups[i]]] <- p

    ## save enrichTable
    enrichment.results <- bind_rows(enrichment.results, table.plot)

  }
  return(list(results=enrichment.results,
              plot.pathway.impact = plot.pathway.impact))
}
