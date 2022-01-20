#' piePlot
#'
#' @param data fold changes data
#' @param path saving path
#' @param cutoff significance cutoff
#' @param save either "pdf", "svg" or "png"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi  dpi only applicable for png
#' @param ...
#'
#' @import tidyverse
#' @import here
#' @import ggplot2
#' @import ggrepel
#' @import RColorBrewer
#' @import ggpubr
#' @import graphics
#' @import grDevices
#' @import readxl
#' @import sjPlot
piePlot <- function (data = data,
                     path = NULL,
                     cutoff = 0.05,
                     save = "pdf",
                     fig.width = 12,
                     fig.height = 9,
                     dpi = 300,
                     ...) {
  if(is.null(path)) {
    path = here()
  } else
    path = path
  if (save == "pdf"){
  pdf(paste(path, "piePlots.pdf", sep = "/"),
      onefile = TRUE)
} else if (save != "pdf") {
    dir.create(paste(here(), "piePlots", sep = "/"))
  }
  metabolite.class <- readxl::read_excel("./inst/extdata/ref/Chemical_annotations_lipids.xlsx")

  ## define metabolites
  data[["MetaboliteClass"]] <- metabolite.class[["SUPER_PATHWAY"]][match(
    data[["Metabolite"]], metabolite.class[["CHEMICAL_NAME"]])]

  ## prepare data
  datPie <- data %>%
    filter(adj.P.Val < cutoff) %>%
    drop_na(MetaboliteClass) %>%
    group_by(contrast, MetaboliteClass) %>%
    summarise(Freq = length(MetaboliteClass)) %>%
    arrange(MetaboliteClass)

  groups <- unique(datPie$contrast)

  ## colors
  colorsOntologyOne <-
    data.frame(
      MetaboliteClass = unique(datPie$MetaboliteClass),
      color = brewer.pal(length(unique(
        datPie$MetaboliteClass
      )), "Paired")
    )
  ## match colors
  matchColumnColors <-
    match(datPie$MetaboliteClass,
          colorsOntologyOne$MetaboliteClass,
          nomatch = 0)
  datPie$color <- c("")
  datPie$color[datPie$MetaboliteClass %in%
                     colorsOntologyOne$MetaboliteClass] <-
    as.character(colorsOntologyOne$color)[matchColumnColors]

  ## plot Pie plots

  for (i in seq_along(groups)) {

    filteredData <- datPie[datPie$contrast %in% groups[i],]

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
      labs(title = paste("Distribution of Metabolites:",
                         groups[i])) +
      theme(
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
      ) +
      guides(fill = guide_legend(ncol = 3)) +
      theme(plot.margin = unit(c(1, 1, 1, 1),
                               "lines"))
    ## print
    print(p)

    if (save != "pdf") {
      save_plot(filename = paste(here(), "piePlots", paste0(groups[i], ".", save), sep = "/"),
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

