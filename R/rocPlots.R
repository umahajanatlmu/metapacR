#' @title rocPlots
#'
#' @description compute and plot roc curves
#'
#' @param dataList normalized data
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param species species to use "hsa" or "mmu"
#' @param group grouping variables
#' @param path saving path
#' @param var.imp  disease for which ROC should be compared---should be one.
#'
#' @import tidyverse
#' @importFrom here here
#' @importFrom RColorBrewer brewer.pal
#' @import stats
#' @importFrom pROC roc coords ggroc
#' @import utils
#' @import graphics
#' @import grDevices
#'
#' @return summaryROC results table. rocPlots in save oblect in defined path.
#'
#' @export

rocPlots <- function(dataList,
                     data.type = c("MH", "Metabolon", "Others"),
                     species = c("hsa", "mmu"),
                     group,
                     path = NULL,
                     var.imp = "PDAC") {
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

  if (is.null(var.imp)) {
    stop("variable importance is missing")
  } else if (length(var.imp) != 1) {
    stop("variable importance variables available....provide only one variable importance variable")
  }

  if (is.null(path)) {
    path <- here::here()
    ifelse(!dir.exists(file.path(paste0(path), "results")),
           dir.create(file.path(paste0(path), "results")),
           FALSE
    )
    path <- paste(path, "results", sep = "/")
  } else {
    path <- path
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

  ## save as pdf
  ## ----------------------------------------------------------------
  pdf(paste(path, "rocCurves.pdf", sep = "/"),
      paper = "a4r",
      onefile = TRUE
  )

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

  ## merge Data
  ## ----------------------------------------------------------------
  data <- merge(metadata.data, imputed.data, by = 0) %>%
    column_to_rownames("Row.names")

  ## convert to characters
  data[[group]] <- as.character(data[[group]])
  data <- data[, !is.na(colnames(data))]

  ## convert to numeric
  data <- data %>%
    select_if(names(.) == paste0(group) | sapply(., is.numeric)) %>%
    mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
    select_if(~ !all(is.na(.)))

  ## unique combinations
  uniqueComparisons <- combn(unique(data[[group]]), 2)

  summaryROC <- data.frame()

  ## progrss bar
  niter <- ncol(data) - 1
  pb <- txtProgressBar(
    min = 0,
    max = niter,
    style = 3,
    char = "="
  )

  for (i in colnames(data)) {
    if (class(data[[i]]) %in% c("character", "factor")) {
      next
    } else if (sum(is.nan(data[[i]])) == nrow(data)) {
      next
    }
    ## plot per comparisons
    for (j in 1:ncol(uniqueComparisons)) {
      uniqueComp <- as.vector(uniqueComparisons[, j])

      dataSubset <-
        data[data[[group]] %in% uniqueComp, ]

      if (uniqueComp[1] == var.imp) {
        dataSubset[["shortCode"]] <-
          ifelse(dataSubset[[group]] == uniqueComp[1], 0, 1)
      } else {
        dataSubset[["shortCode"]] <-
          ifelse(dataSubset[[group]] == uniqueComp[2], 1, 0)
      }
      ## roc
      rocObj <- pROC::roc(
        dataSubset[["shortCode"]],
        dataSubset[[i]]
      )

      ## Get the best threshold
      bestObj <-
        pROC::coords(rocObj, "best", ret = "all", transpose = FALSE)

      rownames(bestObj)[1] <- paste(i,
                                    "_",
                                    uniqueComp[1],
                                    "--",
                                    uniqueComp[2],
                                    sep = ""
      )

      ## add results
      summaryROC <- rbind(summaryROC, bestObj)

      ## plot
      ## ----------------------------------------------------------------
      p <- pROC::ggroc(rocObj, size = 2, color = "#377eb8") +
        geom_segment(
          aes(
            x = 1,
            xend = 0,
            y = 0,
            yend = 1
          ),
          color = "grey",
          linetype = "dashed"
        ) +
        geom_point(
          inherit.aes = FALSE,
          data = bestObj,
          aes(specificity, sensitivity),
          pch = 21,
          size = 5,
          fill = "#e41a1c"
        ) +
        annotate(
          "text",
          x = 0.5,
          y = 0.20,
          size = 4,
          label = paste0("AUC = ", round(rocObj$auc, 3))
        ) +
        annotate(
          "text",
          x = 0.5,
          y = 0.12,
          size = 4,
          label = paste0("Sensitivity = ", round(bestObj$sensitivity, 3))
        ) +
        annotate(
          "text",
          x = 0.5,
          y = 0.04,
          size = 4,
          label = paste0("Specificity = ", round(bestObj$specificity, 3))
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
        ggtitle(paste(
          i,
          ":",
          "\n",
          uniqueComp[1],
          "--",
          uniqueComp[2],
          sep = ""
        )) +
        theme(plot.title = element_text(
          size = 12,
          face = "bold",
          hjust = 0.5
        ))
      ## print
      print(p)

    }
    Sys.sleep(0.01)
    setTxtProgressBar(pb, which(colnames(data) == i))
  }

  dev.off()

  close(pb)
  options(warn = 0) ## reset all warning
  summaryROC <- summaryROC %>%
    rownames_to_column("metabolite") %>%
    extract(metabolite, c("metabolite", "comparison"), "(.*)_([^_]+$)")
  return(summaryROC)
}
