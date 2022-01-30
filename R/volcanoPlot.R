#' @title volcanoPlot
#'
#' @description plot volcano plots
#'
#' @param data fold changes data
#' @param path saving path
#' @param save either "pdf", "svg" or "png"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi  dpi only applicable for png
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
#'
#' @return plots in save object in defined path

volcanoPlot <- function (data,
                         path = NULL,
                         save= c("pdf", "svg","png"),
                         fig.width = 12,
                         fig.height = 9,
                         dpi = 300) {

  stopifnot(inherits(data, "data.frame"))
  validObject(data)

  save <- match.arg(save)

  if(is.null(path)) {
    path = here()
    ifelse(!dir.exists(file.path(paste0(path), "results")),
           dir.create(file.path(paste0(path), "results")),
           FALSE)
    path = paste(path,"results", sep = "/")
  } else
    path = path

  if (save == "pdf"){
  pdf(paste(path, "volcanoPlots.pdf", sep = "/"),
      onefile = TRUE)
  } else if (save != "pdf") {
    dir.create(paste(here(), "volcanoPlots", sep = "/"))
  }

  metabolite.class <- system.file("inst/extdata/ref",
                                  "Chemical_annotations_lipids.csv",
                                  package="metapacR")


  ## define metabolites
  data[["MetaboliteClass"]] <- metabolite.class[["SUPER_PATHWAY"]][match(
    data[["Metabolite"]], metabolite.class[["CHEMICAL_NAME"]])]

  ## prepare volcano data
  datVolcano <- data %>%
    drop_na(MetaboliteClass) %>%
    mutate(foldChanges = log2(logFC))

  groups <- unique(datVolcano$contrast)

  ## colors
  colorsOntologyOne <-
    data.frame(
      MetaboliteClass = unique(datVolcano$MetaboliteClass),
      color = brewer.pal(length(unique(
        datVolcano$MetaboliteClass
      )), "Paired")
    )
  ## match colors
  matchColumnColors <-
    match(datVolcano$MetaboliteClass,
          colorsOntologyOne$MetaboliteClass,
          nomatch = 0)
  datVolcano$color <- c("")
  datVolcano$color[datVolcano$MetaboliteClass %in%
                     colorsOntologyOne$MetaboliteClass] <-
    as.character(colorsOntologyOne$color)[matchColumnColors]

  ## plot volcano plots

  for (i in seq_along(groups)) {

    filteredData <- datVolcano[datVolcano$contrast %in% groups[i],]

    ## plot
    p <-
      filteredData %>%
      ggplot(aes(foldChanges, -log10(adj.P.Val))) +
      geom_point(
        aes(fill = MetaboliteClass),
        size = 3,
        color = "black",
        pch = 21
      ) +
      ggtitle(groups[i]) +
      geom_text_repel(data = head(filteredData, 5),
                      aes(label = Metabolite),
                      min.segment.length = 0) +
      geom_text_repel(data = tail(filteredData, 5),
                      aes(label = Metabolite),
                      min.segment.length = 0) +
      theme_bw() +
      theme(
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.text = element_text(
          size = 11,
          #face = "bold",
          colour = "black"
        ),
        axis.title = element_text(size = 12, face = "bold")
      ) +
      theme(legend.text.align = 0) +
      scale_fill_manual(values = unique(filteredData$color)) +
      labs(fill = "Metabololites category",
           x = "Relative Abundance",
           y = "p value (-log10)")

    if (save == "pdf") {
      ## print
      print(p)
    }

    ## save plots
    if (save != "pdf") {
      save_plot(filename = paste(here(), "volcanoPlots", paste0(groups[i], ".", save), sep = "/"),
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

