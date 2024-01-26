#' @title plotSplitViolin
#'
#' @description plot split violin plot for comparison of two groups.
#'
#' @param dataList raw data from metabolome study
#' @param grouping.variable grouping variable to plot..need to be binary or multinomial variable
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param species species to use "hsa" or "mmu"
#' @param markers list of markers to plots
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
#' @importFrom ggpubr stat_compare_means
#'
#' @return violin plot
#'
#' @export

plotSplitViolin <- function(dataList,
                            grouping.variable = NULL,
                            data.type = c("MH", "Metabolon", "Others"),
                            species = c("hsa", "mmu"),
                            markers = NULL) {
  stopifnot(inherits(dataList, "list"))
  validObject(dataList)
  species <- match.arg(species, c("hsa", "mmu"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  if (is.null(grouping.variable)) {
    stop("group variable is missing....please provide atleast one grouping variable")
  } else if (length(grouping.variable) != 1) {
    stop("multiple group variables available....provide only one group variable")
  }

  if (is.null(markers)) {
    stop("disribution variable is missing....please provide atleast one marker")
  }

  ## add metabolite annotation
  ## ----------------------------------------------------------------
  if (data.type == "Metabolon" && species %in% c("hsa", "mmu")) {
    data("chemicalMetadata")
    metabolite.class <- force(chemicalMetadata)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character)) %>%
      rename(c("Metabolite" = "MET_CHEM_NO",
               "MetaboliteName" = "CHEMICAL_NAME"))
  }

  if (data.type == "MH" && species == "hsa") {

    data("chemicalMetadata_MH")
    metabolite.class <- force(chemicalMetadata_MH)

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character)) %>%
      rename(
        c("MetaboliteClass" = "ONTOLOGY1_NAME",
          "lipidClass" = "ONTOLOGY2_NAME",
          "MetaboliteName" = "METABOLITE_NAME"
        )
      )
  }

  if (data.type == "MH" && species == "mmu") {
    data("chemicalMetadata_MH_mmu")
    metabolite.class <- force(chemicalMetadata_MH_mmu) %>%
      rename(
        c(
          "MetaboliteClass" = "ONTOLOGY1_NAME",
          "lipidClass" = "ONTOLOGY2_NAME",
          "MetaboliteName" = "METABOLITE_NAME"
        )
      )
  }

  if (data.type == "Others") {
    metabolite.class <- Other_metadata

    metabolite.class <- metabolite.class %>%
      mutate(across(everything(), as.character)) %>%
      rename(
        c(
          "MetaboliteClass" = "Ontology_Class",
          "lipidClass" = "Ontology_Subclass",
          "MetaboliteName" = "Metabolite_Name"
        )
      )
  }


  violin.plot <- list()

  ## load imputed data matrix
  ## ----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  colnames(imputed.data) <- metabolite.class$MetaboliteName[match(colnames(imputed.data), metabolite.class$Metabolite)]

  ## subset metadata
  ## ----------------------------------------------------------------
  imputed.data <- imputed.data[, colnames(imputed.data) %in% markers,
                               drop = FALSE
  ]

  ## load metadata
  ## ----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## subset metadata
  ## ----------------------------------------------------------------
  select.columns <- grouping.variable
  row_names <- rownames(metadata.data)
  metadata.data <- metadata.data[, colnames(metadata.data) %in% select.columns,
                                 drop = FALSE
  ]

  ## define factors
  ## ----------------------------------------------------------------
  for (c in colnames(metadata.data)) {
    if (mode(metadata.data[[c]]) %in% c("character", "factor")) {
      metadata.data[[c]] <- as.factor(metadata.data[[c]])
    } else if (mode(metadata.data[[c]]) == "difftime") {
      metadata.data[[c]] <- as.numeric(metadata.data[[c]])
    } else {
      metadata.data[[c]] <- as.numeric(metadata.data[[c]])
    }
  }

  rownames(metadata.data) <- row_names

  ## merge Data
  ## ----------------------------------------------------------------
  data <- bind_cols(metadata.data, imputed.data)

  grouping.factors <- sort(unique(data[[grouping.variable]]))

  if (length(grouping.factors) == 2) {
    data.subset <- data
    ## define groups
    data.subset[[grouping.variable]] <- as.factor(data.subset[[grouping.variable]])

    #data.subset[[grouping.variable]] <- relevel(data.subset[[grouping.variable]], ref = grouping.factors[1])

    ## gather data
    gData <- gather(data.subset, key = metabolite, value = value, -grouping.variable)

    name <- paste(grouping.factors[1], "vs", grouping.factors[2])

    p <- gData %>%
      mutate(metabolite = str_wrap(metabolite, width = 20)) %>%
      ggplot(aes(
        x = metabolite,
        y = value, fill = .data[[grouping.variable]]
      )) +
      .splitViolin() +
      ggpubr::stat_compare_means(aes(group = .data[[grouping.variable]]), label = "p.signif", paired = TRUE) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(2, "Set1")[1:2]) +
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
      ggtitle(name) +
      xlab("") +
      ylab("relative changes") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

    violin.plot <- p

  } else {

    for (i in grouping.factors) {
      data.subset <- data
      ## define groups
      data.subset[[grouping.variable]] <- as.factor(ifelse(data.subset[[grouping.variable]] %in% i,
                                                           i, "rest"
      ))

      data.subset[[grouping.variable]] <- relevel(data.subset[[grouping.variable]], ref = i)

      ## gather data
      gData <- gather(data.subset, key = metabolite, value = value, -grouping.variable)

      name <- paste(i, "vs", "rest")

      p <- gData %>%
        mutate(metabolite = str_wrap(metabolite, width = 20)) %>%
        ggplot(aes(
          x = metabolite,
          y = value, fill = .data[[grouping.variable]]
        )) +
        .splitViolin() +
        ggpubr::stat_compare_means(aes(group = .data[[grouping.variable]]), label = "p.signif", paired = TRUE) +
        scale_fill_manual(values = RColorBrewer::brewer.pal(2, "Set1")[1:2]) +
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
        ggtitle(name) +
        xlab("") +
        ylab("relative changes") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

      violin.plot[[name]] <- p
    }
  }

  return(plot = violin.plot)
}


.GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin",
                                     ggplot2::GeomViolin,
                                     draw_group = function(self,
                                                           data, ...,
                                                           draw_quantiles = NULL) {
                                       data <- transform(data,
                                                         xminv = x - violinwidth * (x - xmin),
                                                         xmaxv = x + violinwidth * (xmax - x)
                                       )
                                       grp <- data[1, "group"]
                                       newdata <- arrange(
                                         transform(data,
                                                   x = if (grp %% 2 == 1) {
                                                     xminv
                                                   } else {
                                                     xmaxv
                                                   }
                                         ),
                                         if (grp %% 2 == 1) {
                                           y
                                         } else {
                                           -y
                                         }
                                       )
                                       newdata <- rbind(
                                         newdata[1, ],
                                         newdata,
                                         newdata[nrow(newdata), ],
                                         newdata[1, ]
                                       )
                                       newdata[c(
                                         1, nrow(newdata) - 1,
                                         nrow(newdata)
                                       ), "x"] <- round(newdata[1, "x"])

                                       if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                         stopifnot(
                                           all(draw_quantiles >= 0),
                                           all(draw_quantiles <= 1)
                                         )
                                         quantiles <- .create_quantile_segment_frame(data, draw_quantiles)
                                         aesthetics <- data[rep(1, nrow(quantiles)),
                                                            setdiff(
                                                              names(data),
                                                              c("x", "y")
                                                            ),
                                                            drop = FALSE
                                         ]
                                         aesthetics$alpha <- rep(1, nrow(quantiles))
                                         both <- cbind(quantiles, aesthetics)
                                         quantile_grob <- GeomPath$draw_panel(both, ...)
                                         .ggname(
                                           "geom_split_violin",
                                           grid::grobTree(
                                             GeomPolygon$draw_panel(newdata, ...),
                                             quantile_grob
                                           )
                                         )
                                       } else {
                                         .ggname(
                                           "geom_split_violin",
                                           GeomPolygon$draw_panel(newdata, ...)
                                         )
                                       }
                                     }
)


.splitViolin <- function(mapping = NULL,
                         data = NULL,
                         stat = "ydensity",
                         position = "identity", ...,
                         draw_quantiles = NULL,
                         trim = TRUE,
                         scale = "area",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = .GeomSplitViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      draw_quantiles = draw_quantiles,
      na.rm = na.rm, ...
    )
  )
}


.create_quantile_segment_frame <- function(data, draw_quantiles) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles) # these are all the y-values for quantiles

  # Get the violin bounds for the requested quantiles.
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)

  # We have two rows per segment drawn. Each segment gets its own group.
  new_data_frame(list(
    x = interleave(violin.xminvs, violin.xmaxvs),
    y = rep(ys, each = 2),
    group = rep(ys, each = 2)
  ))
}



# geom_rangeframe is adapted from ggthemes::geom_rangeframe, but it uses the panel_scales
# to compute the endpoints of the lines rather than the data (as ggthemes::geom_rangeframe does)
.ggname <- function(prefix, grob) {
  # copy of ggthemes:::ggname
  grob$name <- grid::grobName(grob, prefix)
  grob
}
