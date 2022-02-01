#' @title plotMetaboliteAlteration
#'
#' @description compute and plot distribution of signficantly altered metabolites per class.
#'
#' @param data fold changes data
#' @param path saving path
#' @param cutoff significance cutoff
#' @param lipid.class wheather to plot lipid or not
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
#' @import usethis
#'
#' @return save object in defined path.
#'
#' @export

plotMetaboliteAlteration <- function (data,
                              path = NULL,
                              cutoff = 0.01,
                              lipid.class = TRUE,
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
    pdf(paste(path, "metabolite_alterationPlots.pdf", sep = "/"),
        paper = "a4r",
        onefile = TRUE)
  } else if (save != "pdf") {
    dir.create(paste(here(), "metabolite_alterationPlots", sep = "/"))
  }

  metabolite.class <- readRDS("inst/extdata/ref/Chemical_annotations.rds")
  use_data(metabolite.class, overwrite = TRUE)

  ## define metabolites
  data[["MetaboliteClass"]] <- metabolite.class[["SUPER_PATHWAY"]][match(
    data[["Metabolite"]], metabolite.class[["CHEMICAL_NAME"]])]

  ## prepare distibution data
  dat <- data %>%
    filter(adj.P.Val < cutoff) %>%
    drop_na(MetaboliteClass) %>%
    mutate(foldChanges = log2(logFC)) %>%
    select(contrast, MetaboliteClass, logFC) %>%
    mutate(trend = ifelse(log2(logFC) > 1, "up",
                          ifelse(log2(logFC) <= -1, "down", "unchanged"))) %>%
    filter(trend != "unchanged") %>%
    group_by(contrast, MetaboliteClass, trend) %>%
    summarise(FreqClass = length(trend)) %>%
    mutate(FreqClass = ifelse(trend == "down", -1*FreqClass, FreqClass)) %>%
    ungroup() %>%
    arrange(MetaboliteClass)

  groups <- unique(dat$contrast)


  ## plot distribution plots

  for (i in groups) {

    filteredData <- dat[dat$contrast %in% i,]

    ## plot
    p <- ggplot(filteredData,aes(x= MetaboliteClass,
                                        y= FreqClass,
                                        fill = trend)) +
      geom_bar(stat = "identity",
               color = "black",
               size = 0.5) +
      scale_fill_manual(values=c("#377eb8", "#e41a1c")) +
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
      geom_hline(yintercept=0, linetype="solid", color = "black") +
      labs(x=NULL, y="Number of metabolites with sign. changes",
           title=paste("Alteration in metabolites:", i)) +
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
      theme(legend.position='none') +
      theme(axis.ticks.x = element_blank()) +
      theme(axis.text.x = element_text(angle = 90,
                                       vjust = 0.5,
                                       hjust = 0.95))
    if (save == "pdf") {
      ## print
      print(p)
    }
    ## save plots
    if (save != "pdf") {
      save_plot(filename = paste(here(), "metabolite_alterationPlots", paste0(groups[i], ".", save), sep = "/"),
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
      pdf(paste(path, "metabolite_alterationPlots_lipids.pdf", sep = "/"),
          paper = "a4r",
          onefile = TRUE)
    } else if (save != "pdf") {
      dir.create(paste(here(), "metabolite_alterationPlots_lipids", sep = "/"))
    }
    ## prepare distibution data
    dat <- dat %>%
      drop_na(MetaboliteClass) %>%
      mutate(foldChanges = log2(logFC)) %>%
      filter(MetaboliteClass == "Complex lipids") %>%
      mutate(lipidClass = ifelse(grepl("^TAG", Metabolite),
                                 gsub("TAG.*", "TAG", Metabolite),
                                 gsub("[(].*", "", Metabolite))) %>%
      select(contrast, lipidClass, logFC) %>%
      mutate(trend = ifelse(log2(logFC) > 1, "up",
                            ifelse(log2(logFC) <= -1, "down", "unchanged"))) %>%
      filter(trend != "unchanged") %>%
      group_by(contrast, lipidClass, trend) %>%
      summarise(FreqClass = length(trend)) %>%
      mutate(FreqClass = ifelse(trend == "down", -1*FreqClass, FreqClass)) %>%
      ungroup() %>%
      arrange(lipidClass)

    groups <- unique(dat$contrast)

    ## plot distribution plots

    for (i in groups) {

      filteredData <- dat[dat$contrast %in% i,]

      ## plot
      ## plot
      p <- ggplot(filteredData,aes(x= lipidClass,
                                   y= FreqClass,
                                   fill = trend)) +
        geom_bar(stat = "identity",
                 color = "black",
                 size = 0.5) +
        scale_fill_manual(values=c("#377eb8", "#e41a1c")) +
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
        theme(
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.text = element_text(
            size = 11,
            #face = "bold",
            colour = "black"
          ),
          axis.title = element_text(size = 12, face = "bold")
        ) +
        geom_hline(yintercept=0, linetype="solid", color = "black") +
        labs(x=NULL, y="Number of metabolites with sign. changes",
             title=paste("Alteration in metabolites:", i)) +
        ggplot_theme +
        theme(legend.position='none') +
        theme(axis.ticks.x = element_blank()) +
        theme(axis.text.x = element_text(angle = 90,
                                         vjust = 0.5,
                                         hjust = 0.95))
      if (save == "pdf") {
        ## print
        print(p)
      }
      ## save plots
      if (save != "pdf") {
        save_plot(filename = paste(here(), "metabolite_alterationPlots_lipids", paste0(groups[i], ".", save), sep = "/"),
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
