#' correlationPlot
#'
#' @param data normalized data
#' @param drop.var grouping variables to drop
#' @param ...
#'
#' @import tidyverse
#' @import Hmisc
#' @import pheatmap
#' @import dichromat
#'
#' @return
#' @export
#'
#' @examples
#' correlationPlot(data=data, cluster_size=6, drop.var=NULL)
correlationPlot <- function (data = data,
                             drop.var = NULL, ...) {
  ## drop varibales
  data <- data[, !colnames(data) %in% drop.var]

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

  tree_cut <- cutree(p$tree_row)

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

  return(list(plot = p,
              result=corrMat))

}
