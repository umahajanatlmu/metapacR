#' distributionPlot
#'
#' @param data fold changes data
#' @param path saving path
#' @param cutoff cutoff of p value
#' @param lipid.class TRUE/FALSE to plot lipid classes
#' @param save either "pdf", "svg" or "png"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi  dpi only applicable for png
#' @param ...   ggplot extensions
#'
#' @import tidyverse
#' @import here
#' @import ggplot2
#' @import ggrepel
#' @import RColorBrewer
#' @import ggpubr
#' @import graphics
#' @import grDevices
distributionPlot <- function (data = data,
                              path = NULL,
                              cutoff = 0.01,
                              lipid.class = TRUE,
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
    pdf(paste(path, "distributionPlots.pdf", sep = "/"),
        paper = "a4r",
        onefile = TRUE)
  } else if (save != "pdf") {
    dir.create(paste(here(), "distributionPlots", sep = "/"))
  }

  metabolite.class <- system.file("inst/extdata/ref", "Chemical_annotations_lipids.csv", package="metapacR")

  ## define metabolites
  data[["MetaboliteClass"]] <- metabolite.class[["SUPER_PATHWAY"]][match(
    data[["Metabolite"]], metabolite.class[["CHEMICAL_NAME"]])]

  ## prepare distibution data
  dat <- data %>%
    drop_na(MetaboliteClass) %>%
    mutate(foldChanges = log2(logFC))

  groups <- unique(dat$contrast)

  ## colors
  colorsOntologyOne <-
    data.frame(
      MetaboliteClass = unique(dat$MetaboliteClass),
      color = brewer.pal(length(unique(
        dat$MetaboliteClass
      )), "Paired")
    )
  ## match colors
  matchColumnColors <-
    match(dat$MetaboliteClass,
          colorsOntologyOne$MetaboliteClass,
          nomatch = 0)
  dat$color <- c("")
  dat$color[dat$MetaboliteClass %in%
              colorsOntologyOne$MetaboliteClass] <-
    as.character(colorsOntologyOne$color)[matchColumnColors]

  ## plot distribution plots

  for (i in seq_along(groups)) {

    filteredData <- dat[dat$contrast %in% groups[i],]

    ## plot
    p <- ggplot(filteredData,
                aes(x = Metabolite,
                    y = -log10(adj.P.Val),
                    fill = MetaboliteClass)) +
      geom_point(aes(size = logFC),
                 shape = 21,
                 alpha = 0.5,
                 show.legend = FALSE)  +
      facet_grid(. ~ MetaboliteClass) +
      scale_fill_manual(values = unique(filteredData$color)) +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line = element_line(size = 0.75),
            axis.title = element_text(size = 11, face = "bold")) +
      geom_hline(yintercept = -log10(cutoff), linetype='dotted') +
      geom_hline(yintercept = -log10(cutoff*5), linetype='dotted') +
      geom_text_repel(aes(label = ifelse(filteredData$adj.P.Val < cutoff,
                                         filteredData$Metabolite, NA)),
                      min.segment.length = 0,
                      size = 2,
                      show.legend = FALSE) +
      ggtitle(groups[i]) +
      xlab("")
    ## print
    print(p)

    ## save plots
    if (save != "pdf") {
      save_plot(filename = paste(here(), "distributionPlots", paste0(groups[i], ".", save), sep = "/"),
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
  ## lipid class
  if (isTRUE(lipid.class)) {
    if (save == "pdf"){
      pdf(paste(path, "distributionPlots_lipids.pdf", sep = "/"),
          paper = "a4r",
          onefile = TRUE)
    } else if (save != "pdf") {
      dir.create(paste(here(), "distributionPlots_lipids", sep = "/"))
    }
    ## prepare distibution data
    dat <- dat %>%
      drop_na(MetaboliteClass) %>%
      mutate(foldChanges = log2(logFC)) %>%
      filter(MetaboliteClass == "Complex lipids") %>%
      mutate(lipidClass = ifelse(grepl("^TAG", Metabolite),
                                 gsub("TAG.*", "TAG", Metabolite),
                                 gsub("[(].*", "", Metabolite)))

    groups <- unique(dat$contrast)

    ## colors
    colorsOntologyOne <-
      data.frame(
        lipidClass = unique(dat$lipidClass),
        color = brewer.pal(length(unique(
          dat$lipidClass
        )), "Paired")
      )
    ## match colors
    matchColumnColors <-
      match(dat$lipidClass,
            colorsOntologyOne$lipidClass,
            nomatch = 0)
    dat$color <- c("")
    dat$color[dat$lipidClass %in%
                colorsOntologyOne$lipidClass] <-
      as.character(colorsOntologyOne$color)[matchColumnColors]

    ## plot distribution plots

    for (i in seq_along(groups)) {

      filteredData <- dat[dat$contrast %in% groups[i],]

      ## plot
      p <- ggplot(filteredData,
                  aes(x = Metabolite,
                      y = -log10(adj.P.Val),
                      fill = lipidClass)) +
        geom_point(aes(size = logFC),
                   shape = 21,
                   alpha = 0.5,
                   show.legend = FALSE)  +
        facet_grid(. ~ lipidClass) +
        scale_fill_manual(values = unique(filteredData$color)) +
        theme_classic() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line = element_line(size = 0.75),
              axis.title = element_text(size = 11, face = "bold")) +
        geom_hline(yintercept = -log10(cutoff), linetype='dotted') +
        geom_hline(yintercept = -log10(cutoff*5), linetype='dotted') +
        geom_text_repel(aes(label = ifelse(filteredData$adj.P.Val < cutoff,
                                           filteredData$Metabolite, NA)),
                        min.segment.length = 0,
                        size = 2,
                        show.legend = FALSE) +
        ggtitle(groups[i]) +
        xlab("")
      ## print
      print(p)

      ## save plots
      if (save != "pdf") {
        save_plot(filename = paste(here(), "distributionPlots_lipids", paste0(groups[i], ".", save), sep = "/"),
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
}
