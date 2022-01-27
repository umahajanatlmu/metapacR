#' lipidChainLengthDistribution
#'
#' @param results fold changes data
#' @param p.value.cutoff cutoff of p-values to be used
#' @param fold.changes.cutoff higher fold changes cutoff to be used
#' @param path saving path
#' @param save either "pdf", "svg" or "png"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi dpi only applicable for png
#' @param ... ggplot extension
#'
#' @import tidyverse
#' @import here
#' @import ggplot2
#' @import ggrepel
#' @import sjPlot
#' @import RColorBrewer
#' @import ggpubr
#' @import graphics
#' @import grDevices
lipidChainLengthDistribution <- function (results=results,
                                          p.value.cutoff=0.05,
                                          fold.changes.cutoff=1.5,
                                          path = NULL,
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
    pdf(paste(path, "chainLengthDistribution.pdf", sep = "/"),
        paper = "a4r",
        onefile = TRUE)
  } else if (save != "pdf") {
    dir.create(paste(here(), "chainLengthDistribution", sep = "/"))
  }


  ## load annotation file
  chemicalMetadata <- system.file("inst/extdata/ref",
                                  "Chemical_annotations.csv",
                                  package="metapacR")

  ## annotate for class
  results$class <- chemicalMetadata$SUPER_PATHWAY[match(results$Metabolite,
                                                        chemicalMetadata$CHEMICAL_NAME)]
  ## subset chain lengths
  results <- results %>%
    filter(adj.P.Val < p.value.cutoff, logFC > fold.changes.cutoff | logFC < (fold.changes.cutoff - 1)) %>%
    filter(class == "Complex lipids") %>%
    mutate(newMet = Metabolite) %>%
    mutate(newMet = gsub("O-|P-", "", newMet)) %>%
    separate(newMet, c("lipid.class", "fatty.acid"), "[()]|-") %>%
    mutate(fatty.acid=gsub(".*/", "", fatty.acid)) %>%
    mutate(lipid.class = case_when(str_detect(lipid.class, "TAG") ~ "TAG",
                                   TRUE ~ lipid.class)) %>%
    mutate(fatty.acid = gsub("FA", "", fatty.acid)) %>%
    select(contrast,logFC, adj.P.Val, lipid.class, fatty.acid) %>%
    group_by(contrast, lipid.class, fatty.acid) %>%
    summarise_all(funs(mean)) %>%
    arrange(fatty.acid) %>%
    mutate(new.fatty.acid = fatty.acid) %>%
    separate(new.fatty.acid, c("chain.length","saturation"), ":") %>%
    mutate(chain.length = as.numeric(as.character(chain.length))) %>%
    mutate(saturation.class = ifelse(saturation == "0", "saturated",
                                     ifelse(saturation == "1", "mono-unsaturated",
                                            "poly-unsaturated")))

  groups <- unique(results$contrast)
  for (i in groups) {

    p <- results %>%
      filter(contrast== i) %>%
      ggplot(aes(x=lipid.class,
                 y = fatty.acid,
                 color = log2(logFC),
                 size = -log10(adj.P.Val),
                 shape = saturation.class)) +
      geom_point() +
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
      scale_colour_gradientn(colours = rev(brewer.pal(10, "RdYlBu")),
                             limits = c(-2,2),
                             oob = scales::squish,
                             name = 'fold changes') +
      guides(colour = guide_colourbar(barwidth = unit(0.3, "cm"),
                                      ticks.colour = "black",
                                      frame.colour = "black")) +
      labs(title=i,
           x="",
           y="fatty acid chain length",
           size="-log10(P-value)",
           shape="saturation") +
      theme(axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = 0.5))

    if (save == "pdf") {
      ## print
      print(p)
    }
    ## save plots
    if (save != "pdf") {
      save_plot(filename = paste(here(), "chainLengthDistribution", paste0(groups[i], ".", save), sep = "/"),
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
