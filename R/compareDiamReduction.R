#' @title  compareDiamReduction
#'
#' @description function to compare different diamentionality reduction methods, namely, PCA, OPLS, TSNE and UMAP.
#'
#' @param dataList raw metabolome data list from imputeTransformScale function.It need to have imputed.matrix and metadata.
#' @param plotting.variable plotting grouping variable..should be 1
#' @param crossvalI number of cross-validation segments
#'
#' @import tidyverse
#' @importFrom ropls opls
#' @importFrom umap umap
#' @importFrom Rtsne Rtsne
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr ggarrange
#' @import utils
#' @import stats
#' @import graphics
#' @import grDevices
#'
#' @return Multivariate analyses results in list object.
#'   The object contains the following:\itemize{
#'     \item plot comparative diamemtionality reduction plot
#'     \item pca S4 object of pca results
#'     \item opls S4 object of opls results
#'     \item rtsne S4 object of tsne results
#'     \item rtsme.pca S4 object of tsne + pca results
#'     \item kmeans S4 object of kmeans results
#'     \item umap S4 object of umap results
#'   }
#'
#' @export

compareDiamReduction <- function(dataList,
                                 plotting.variable = NULL,
                                 crossvalI = 7) {
  stopifnot(inherits(dataList, "list"))
  validObject(dataList)


  if (is.null(plotting.variable)) {
    stop("plotting variable is missing")
  } else if (length(plotting.variable) != 1) {
    stop("multiple plotting variables available....provide only one plotting variable")
  }


  ## load imputed data matrix
  ## ----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  ## load metadata
  ## ----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## select row for which metadata is available
  ## ----------------------------------------------------------------
  imputed.data <- imputed.data[rownames(metadata.data), ]

  ## subset metadata
  ## ----------------------------------------------------------------
  select.columns <- plotting.variable
  row_names <- rownames(metadata.data)
  metadata.data <- metadata.data[, colnames(metadata.data) %in% select.columns, drop = FALSE]

  group <- plotting.variable

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


  ## number of groups
  nVar <- length(unique(metadata.data[[group]]))

  # perform pca
  ## ----------------------------------------------------------------
  ## select numeric data, filter NaN and NAs
  dataNumeric <- imputed.data %>%
    select_if(., is.numeric) %>%
    mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
    select_if(~ !all(is.na(.))) %>%
    as.data.frame()

  # perform Kmeans clustering
  ## ----------------------------------------------------------------
  print(paste("Performing kmeans for", nVar, "clusters...."))

  kmeans <- kmeans(dataNumeric,
                   centers = nVar,
                   nstart = 100
  )

  print("Performing PCA ....")
  ## perform PCA
  pca <- prcomp(dataNumeric)
  ## plot PCA
  plot.dat.pca <- as.data.frame(pca$x)
  plot.dat.pca <- bind_cols(plot.dat.pca, metadata.data)

  ## PCA componants
  pc1var <- round(summary(pca)$importance[2, 1] * 100, digits = 2)
  pc2var <- round(summary(pca)$importance[2, 2] * 100, digits = 2)

  ## plot
  p_pca <- ggplot(plot.dat.pca, aes(
    x = PC1,
    y = PC2,
    color = .data[[group]]
  )) +
    geom_point(
      size = 3,
      alpha = 0.8
    ) +
    xlab(paste0("PC1: ", pc1var, "%")) +
    ylab(paste0("PC2: ", pc2var, "%")) +
    labs(color = group) +
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
    scale_color_manual(values = RColorBrewer::brewer.pal(
      nVar,
      "Set1"
    )) +
    theme(legend.position = "right") +
    ggtitle("PCA")

  print("Performing OPLS ....")
  # perform opls
  ## ----------------------------------------------------------------
  ## perform opls
  opls <- ropls::opls(dataNumeric,
                      scaleC = "none",
                      y = metadata.data[[group]],
                      log10L = FALSE,
                      crossvalI = crossvalI
  )
  ## plot opls
  plot.dat.opls <- data.frame(opls@scoreMN)
  plot.dat.opls <- bind_cols(plot.dat.opls, metadata.data)

  N <- nrow(plot.dat.opls)
  pscores <- plot.dat.opls[["p1"]]
  oscores <- plot.dat.opls[["p2"]]
  hotFisN <- (N - 1) * 2 * (N^2 - 1) / (N^2 * (N - 2)) * qf(0.95, 2, N - 2)

  ## plot
  p_opls <- ggplot(plot.dat.opls, aes(x = p1, y = p2, color = .data[[group]])) +
    gg_circle(
      rx = sqrt(var(pscores) * hotFisN),
      ry = sqrt(var(oscores) * hotFisN),
      xc = 0,
      yc = 0,
      size = 1,
      color = "black",
      fill = "white"
    ) +
    geom_hline(yintercept = 0, size = 1, color = "black") +
    geom_vline(xintercept = 0, size = 1, color = "black") +
    geom_point(
      size = 3,
      alpha = 0.75
    ) +
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
    scale_color_manual(values = RColorBrewer::brewer.pal(nVar, "Set1")) +
    theme(legend.position = "right") +
    ggtitle(opls@typeC) +
    xlab(paste("t1 =", opls@modelDF[1, 1] * 100, "%", "[between groups variation] ", sprintf("\u2192"))) +
    ylab(paste("t2 =", opls@modelDF[2, 1] * 100, "%", "[within groups variation] ", sprintf("\u2192"))) +
    theme(panel.grid.major = element_line(colour = "grey93")) +
    annotate(
      geom = "text", label = paste("R2X =", opls@summaryDF[[1]]),
      x = min(plot.dat.opls$p1),
      y = max(plot.dat.opls$p2),
      vjust = 1,
      size = 4
    ) +
    annotate(
      geom = "text", label = paste("R2Y =", opls@summaryDF[[2]]),
      x = min(plot.dat.opls$p1),
      y = 0.9 * max(plot.dat.opls$p2),
      vjust = 1,
      size = 4
    ) +
    annotate(
      geom = "text", label = paste("RMSEE =", opls@summaryDF[[4]]),
      x = min(plot.dat.opls$p1),
      y = 0.8 * max(plot.dat.opls$p2),
      vjust = 1,
      size = 4
    ) +
    labs(color = group) +
    theme(panel.background = element_rect(fill = "grey90"))

  print("Performing t-SNE....")
  # perform t-SNE plus PCA
  ## ----------------------------------------------------------------
  rtsne.pca <- Rtsne::Rtsne(dataNumeric,
                            perplexity = floor((nrow(dataNumeric)-1)/3),
                            dims = 2)

  ## plot rtsne
  plot.dat.rtsne.pca <- as.data.frame(rtsne.pca$Y)
  rownames(plot.dat.rtsne.pca) <- rownames(dataNumeric)
  plot.dat.rtsne.pca <- bind_cols(plot.dat.rtsne.pca, metadata.data)

  ## plot
  p_rtsne.pca <- ggplot(plot.dat.rtsne.pca, aes(
    x = V1,
    y = V2,
    color = .data[[group]]
  )) +
    geom_point(
      size = 3,
      alpha = 0.8
    ) +
    xlab("tSNE1") +
    ylab("tSNE2") +
    labs(color = group) +
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
    scale_color_manual(values = RColorBrewer::brewer.pal(
      nVar,
      "Set1"
    )) +
    theme(legend.position = "right") +
    ggtitle("PCA + tSNE")

  print("Performing t-SNE and PCA ....")
  # perform t-SNE w/o PCA
  ## ----------------------------------------------------------------
  rtsne <- Rtsne::Rtsne(dataNumeric,
                        perplexity = floor((nrow(dataNumeric)-1)/3),
                        dims = 2,
                        pca=FALSE)

  ## plot rtsne
  plot.dat.rtsne <- as.data.frame(rtsne$Y)
  rownames(plot.dat.rtsne) <- rownames(dataNumeric)
  plot.dat.rtsne <- bind_cols(plot.dat.rtsne, metadata.data)

  ## plot
  p_rtsne <- ggplot(plot.dat.rtsne, aes(
    x = V1,
    y = V2,
    color = .data[[group]]
  )) +
    geom_point(
      size = 3,
      alpha = 0.8
    ) +
    xlab("tSNE1") +
    ylab("tSNE2") +
    labs(color = group) +
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
    scale_color_manual(values = RColorBrewer::brewer.pal(
      nVar,
      "Set1"
    )) +
    theme(legend.position = "right") +
    ggtitle("tSNE")

  print("Performing UMAP....")
  # perform umap
  ## ----------------------------------------------------------------
  n_neighbors <- min(nrow(dataNumeric) - 1, 10)
  umap <- umap::umap(dataNumeric, n_neighbors = n_neighbors)

  ## plot rtsne
  plot.dat.umap <- as.data.frame(umap$layout)
  plot.dat.umap <- bind_cols(plot.dat.umap, metadata.data)

  p_umap <- ggplot(plot.dat.umap, aes(
    x = V1,
    y = V2,
    color = .data[[group]]
  )) +
    geom_point(
      size = 3,
      alpha = 0.8
    ) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    labs(color = group) +
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
    scale_color_manual(values = RColorBrewer::brewer.pal(
      nVar,
      "Set1"
    )) +
    theme(legend.position = "right") +
    ggtitle("UMAP")

  p <- ggpubr::ggarrange(p_pca,
                         p_opls,
                         p_rtsne,
                         p_rtsne.pca,
                         p_umap,
                         ncol = 3,
                         nrow = 2,
                         common.legend = TRUE,
                         legend = "bottom"
  )

  return(list(
    plot = p,
    pca = pca,
    opls = opls,
    tsne = rtsne,
    tsne_pca = rtsne.pca,
    kmeans = kmeans,
    umap = p_umap
  ))
}


# Function to plot Hotelling's T-squared ellipse
# Adapted from https://github.com/tyrannomark/bldR/blob/master/R/L2017.R
# GPL-3 license

#' gg_circle
#'
#' plot opls circle
#'
#' @param rx x axis cordinates on r
#' @param ry y axis cordinates on r
#' @param xc x axis cordinate on c
#' @param yc y axis cordinate on c
#' @param color circle border color
#' @param fill circle fill
#' @param ... ggplot extensions
gg_circle <- function(rx, ry, xc, yc, color = "black", fill = NA, ...) {
  x <- xc + rx * cos(seq(0, pi, length.out = 100))
  ymax <- yc + ry * sin(seq(0, pi, length.out = 100))
  ymin <- yc + ry * sin(seq(0, -pi, length.out = 100))
  annotate(
    "ribbon",
    x = x, ymin = ymin, ymax = ymax,
    color = color, fill = fill, ...
  )
}
