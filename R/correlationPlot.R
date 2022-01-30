#' @title correlationPlot
#'
#' @description function to plot intragroup correlation plot.
#'
#' @param dataList normalized data
#' @param h cluster size
#' @param save either "pdf", "svg" or "png"
#' @param fig.width plot width not applicable for pdf
#' @param fig.height plot height not applicable for pdf
#' @param path  saving path
#'
#' @import tidyverse
#' @import Hmisc
#' @import pheatmap
#' @import graphics
#' @import grDevices
#' @import here
#'
#' @return correlation matrix as a list and pdf output in defined path folder.

correlationPlot <- function (dataList,
                             h = 3,
                             path = NULL,
                             save= c("pdf", "svg","png"),
                             fig.width = 12,
                             fig.height = 9) {

  stopifnot(inherits(dataList, "list"))
  validObject(dataList)
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
    pdf(paste(path, "correlationPlots.pdf", sep = "/"),
        onefile = TRUE)
  } else if (save != "pdf") {
    dir.create(paste(here(), "correlationPlots.pdf", sep = "/"))
  }

  ## load imputed data matrix
  ##----------------------------------------------------------------
  data <- dataList[["imputed.matrix"]]

  ## convert to numeric
  data <- data %>%
    select_if(., is.numeric) %>%
    mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    select_if(~!all(is.na(.)))

  ## make correlation matrix
  corrMat <- rcorr(as.matrix(data))

  ## define matrix colors
  macolor = colorRampPalette(c("navyblue", "white", "red"))(100)

  ## define matrix
  M <- corrMat$r
  p_mat <- corrMat$p
  row_names <- rownames(M)

  ## plot pheatmap
  p <- pheatmap(M,
                color = macolor,
                silent = TRUE)

  tree_cut <- cutree(p$tree_row, h = h)

  tc <- data.frame(tip = names(tree_cut),
                   clust_membership = as.character(unname(
                     tree_cut)))
  row.names(tc) <- tc$tip
  tc <- tc['clust_membership']

  ## plot heatmap
  p <- pheatmap(M,
                color = rev(macolor),
                clustering_method = "complete",
                annotation_row = tc,
                annotation_col = tc,
                fontsize_row = 0.8,
                fontsize_col = 0.8)

  if (save == "svg") {
    svg(paste(path, "correlationPlots.svg", sep = "/"), height = fig.height, width=fig.width)
  } else if (save == "png") {
    svg(paste(path, "correlationPlots.png", sep = "/"), height = fig.height, width=fig.width)
  }

  ## print
  print(p)

  ## dev off
  dev.off()

  return(result=corrMat)

}
