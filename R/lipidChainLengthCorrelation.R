#' @title lipidChainLengthCorrelation
#'
#' @description compute correlation of chain length with abundance.
#'
#' @param results fold changes data
#' @param path saving path
#' @param save either "pdf", "svg" or "png"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param dpi dpi only applicable for png
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
#' @import ggpmisc
#'
#' @return plot as save object in defined path.
#'
#' @export

lipidChainLengthCorrelation <- function (results,
                                         path = NULL,
                                         save = c("pdf", "svg", "png"),
                                         fig.width = 12,
                                         fig.height = 9,
                                         dpi = 300) {

  stopifnot(inherits(results, "data.frame"))
  validObject(results)

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
    pdf(paste(path, "chainLengthCorrelation.pdf", sep = "/"),
        paper = "a4r",
        onefile = TRUE)
  } else if (save != "pdf") {
    dir.create(paste(here(), "chainLengthCorrelations", sep = "/"))
  }


  ## load annotation file
  chemicalMetadata <- system.file("inst/extdata/ref",
                                  "Chemical_annotations.csv",
                                  package="metapacR")

  ## annotate for class
  results$class <- chemicalMetadata$SUPER_PATHWAY[match(results$Metabolite,
                                                    chemicalMetadata$CHEMICAL_NAME)]
  ## subset chain lengths
  res <- results %>%
    filter(class == "Complex lipids") %>%
    mutate(newMet = Metabolite) %>%
    mutate(newMet = gsub("O-|P-", "", newMet)) %>%
    separate(newMet, c("lipid.class", "fatty.acid"), "[()]|-") %>%
    mutate(fatty.acid=gsub(".*/", "", fatty.acid)) %>%
    mutate(lipid.class = case_when(str_detect(lipid.class, "TAG") ~ "TAG",
                                   TRUE ~ lipid.class)) %>%
    mutate(fatty.acid = gsub("FA", "", fatty.acid)) %>%
    separate(fatty.acid, c("chain.length","saturation"), ":") %>%
    mutate(chain.length = as.numeric(as.character(chain.length))) %>%
    mutate(saturation.class = ifelse(saturation == "0", "saturated",
                                     ifelse(saturation == "1", "mono-unsaturated",
                                            "poly-unsaturated")))

  ## plot
  groups <- unique(res$contrast)

  for (i in groups) {

    res.subset <- res %>%
      filter(contrast == i)

    p <- ggplot(res.subset,
                aes(x=chain.length,
                    y=t.ratio,
                    color=saturation.class)) +
      geom_hline(yintercept = 2, linetype="dotted") +
      geom_hline(yintercept = -2, linetype="dotted") +
      geom_point(alpha=0.7) +
      geom_smooth(method = "lm",
                  formula = y ~ x,
                  se = FALSE,
                  show.legend = FALSE) +
      stat_fit_glance(method = "lm",
                      label.x="left",
                      label.y="bottom",
                      method.args = list(formula = y ~ x),
                      aes(label = sprintf('R^2~"="~%.3f~~italic(p)~"="~%.2f',
                                          stat(..r.squared..),
                                          stat(..p.value..))
                      ),
                      parse = TRUE,
                      size = 4) +
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
      theme(axis.ticks.x = element_blank()) +
      ggtitle(i) +
      xlab("Fatty acid chain length") +
      ylab("t-ratio") +
      facet_wrap(~lipid.class, ncol=7) +
      scale_color_manual(values = brewer.pal(3, "Set1")) +
      theme(legend.position = "bottom") +
      labs(color="Saturation status")

    if (save == "pdf") {
    ## print
    print(p)
    }
    ## save plots
    if (save != "pdf") {
      save_plot(filename = paste(here(), "chainLengthCorrelations", paste0(groups[i], ".", save), sep = "/"),
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
