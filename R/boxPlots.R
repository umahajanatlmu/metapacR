#' boxPlots
#'
#' @param data normalized data
#' @param group grouping variables
#' @param drop.var drop un-necessary variables
#' @param path saving path
#'
#' @import tidyverse
#' @import here
#' @import ggplot2
#' @import RColorBrewer
#' @import ggpubr
#' @import rstatix
#' @import stats
#'
#' @return
#' @export
#'
#' @examples
#' boxPlots(data=data, group="var", drop.var = NULL, path = NULL
boxPlots <- function(data = data,
                     group = group,
                     drop.var = NULL,
                     path = NULL) {
  if(is.null(path)) {
    path = here()
  } else
    path = path

  pdf(paste(path, "boxplots.pdf", sep = "/"),
      paper= "a4r",
      onefile = TRUE)

  ## drop varibales
  data <- data[, !colnames(data) %in% drop.var]

  ## convert to characters
  data[[group]] <- as.character(data[[group]])

  ## convert to numeric
  data <- data %>%
    select_if(names(.) == paste0(group) | sapply(., is.numeric)) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    select_if(~!all(is.na(.)))

  ## progrss bar
  niter <- ncol(data) - 1
  pb <- txtProgressBar(min = 0,
                       max = niter,
                       style = 3,
                       char = "=")

  for(j in 1:niter) {

    for (i in colnames(data)) {

      if ( i == group) {
        next
      }
      ## stats
      plot.Dat <- data %>%
        select(i, group)

      formula <-
        as.formula(paste0("`", i,"`", "~", paste0(group)))

      stat <- plot.Dat %>%
        pairwise_t_test(formula, p.adjust.method = "BH") %>%
        add_xy_position(
          x = paste0(group),
          fun = "max",
          dodge = 1,
          step.increase = 0.15
        ) %>%
        filter(p.adj.signif != "ns")

      ## plot
      p <- ggplot(plot.Dat,
                  aes(x = .data[[group]],
                      y = .data[[i]],
                      fill = .data[[group]])) +
        geom_boxplot(width= 0.1) +
        geom_violin(alpha = 0.5,
                    show.legend = FALSE) +
        geom_jitter(
          shape = 21,
          size = 0.5,
          width = 0.2,
          color = "black",
          alpha = 0.3,
          show.legend = TRUE
        ) +
        scale_fill_manual(values = brewer.pal(length(unique(
          plot.Dat[[group]])), "Set1")) +
        theme_bw() +
        theme(
          axis.line = element_line(size = 0.75),
          axis.text = element_text(
            size = 11,
            face = "bold",
            colour = "black"
          ),
          axis.title = element_text(size = 12, face = "bold")
        ) +
        theme(axis.ticks.x = element_blank()) +
        theme(legend.position = "bottom",
              legend.box="vertical") +
        stat_pvalue_manual(
          stat,
          label = "p.adj.signif",
          tip.length = 0.01,
          color = "gray50",
          size = 6,
          bracket.size = 0.8,
          inherit.aes = FALSE
        ) +
        xlab("") +
        labs( fill = paste0(group)) +
        theme(axis.text.x = element_blank()) +
        ylab("median scaled log10 [absolute concentration]") +
        ggtitle(i) +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
      ## print
      print(p)
    }
    Sys.sleep(0.01)
    setTxtProgressBar(pb, j)
  }
  dev.off()
  close(pb)
}
