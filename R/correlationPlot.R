#' @title correlationPlot
#'
#' @description function to plot intragroup correlation plot.
#'
#' @param dataList normalized data
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param species species to use "hsa" or "mmu"
#' @param group grouping variables
#'
#' @import tidyverse
#' @importFrom Hmisc rcorr
#' @importFrom psych corr.test
#' @importFrom ggnewscale new_scale_fill
#' @importFrom reshape2 melt
#' @import graphics
#' @import grDevices
#' @import here
#'
#' @return correlation matrix as a list and pdf output in defined path folder.
#'
#' @export

correlationPlot <- function(dataList,
                            data.type = c("MH", "Metabolon", "Others"),
                            species = c("hsa", "mmu"),
                            group = NULL,
                            var.imp = NULL) {

  options(warn = -1) ## supress all warning

  stopifnot(inherits(dataList, "list"))
  validObject(dataList)

  species <- match.arg(species, c("hsa", "mmu"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  if (is.null(group)) {
    stop("group variable is missing")
  } else if (length(group) != 1) {
    stop("multiple group variables available....provide only one group variable")
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


  ## load imputed data matrix
  ## ----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  colnames(imputed.data) <- metabolite.class$MetaboliteName[match(colnames(imputed.data), metabolite.class$Metabolite)]

  ## load metadata
  ## ----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## subset metadata
  ## ----------------------------------------------------------------
  select.columns <- group
  row_names <- rownames(metadata.data)
  metadata.data <- metadata.data[, colnames(metadata.data) %in% select.columns, drop = FALSE]

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

  ## convert to characters
  data[[group]] <- as.character(data[[group]])

  data <- data[, !is.na(colnames(data))]

  ## unique combinations
  if (length(unique(data[[group]])) == 2) {

    corr_mat_list <- list()

    for (i in unique(data[[group]])) {
      temp <- data[data[[group]] == i, !colnames(data) %in% group]

      temp <- temp[, colSums(is.na(temp)) != nrow(temp)]

      #temp <- t(temp)
      corr_mat <- psych::corr.test(as.matrix(temp),
                                   y = NULL,
                                   use = "complete.obs",
                                   method="kendall",
                                   adjust="fdr",
                                   alpha=.05,
                                   ci=FALSE,
                                   normal=FALSE)
      ## save corr mat
      corr_mat_list[[i]] <- corr_mat
    }
  }
  else if (length(unique(data[[group]])) > 2) {

    if (is.null(var.imp)) {
      stop("provide one var.imp")
    }

    data[[group]] <- ifelse(data[[group]] == var.imp, var.imp, "rest")


    corr_mat_list <- list()

    for (i in unique(data[[group]])) {
      temp <- data[data[[group]] == i, !colnames(data) %in% group]
      #temp <- t(temp)
      corr_mat <- psych::corr.test(as.matrix(temp),
                                   y = NULL,
                                   use = "complete.obs",
                                   method="kendall",
                                   adjust="fdr",
                                   alpha=.05,
                                   ci=FALSE,
                                   minlength=5,
                                   normal=FALSE)
      ## save corr mat
      corr_mat_list[[i]] <- corr_mat
    }


  }

  return(list(results = corr_mat_list))

  options(warn = 0) ## reset all warnings

}

