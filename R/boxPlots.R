
#' @title boxPlots
#'
#' @description function to plot box and violin plot from normalized metabolome data
#'
#' @param dataList raw metabolome data list from imputeTransformScale function.It need to have imputed.matrix and metadata.
#' @param data.type  select platform used, c("MH", "Metabolon", "Others")
#' @param species species to use "hsa" or "mmu"
#' @param group plotting grouping variable..should be 1
#' @param path saving path
#'
#' @import tidyverse
#' @importFrom here here
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr stat_pvalue_manual
#' @importFrom rstatix pairwise_t_test add_xy_position
#' @import stats
#' @import graphics
#' @import grDevices
#'
#' @return returns compiled pdf file.
#'
#' @export

boxPlots <- function(dataList,
                     data.type = c("MH", "Metabolon", "Others"),
                     species = c("hsa", "mmu"),
                     group = NULL,
                     path = NULL) {

  stopifnot(inherits(dataList, "list"))
  validObject(dataList)

  species <- match.arg(species, c("hsa", "mmu"))
  data.type <- match.arg(data.type, c("MH", "Metabolon", "Others"))

  if (is.null(group)) {
    stop("grouping variable is missing")
  } else if (length(group) != 1) {
    stop("multiple grouping variables available....provide only one grouping variable")
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

  ## save as pdf
  ## ----------------------------------------------------------------
  pdf(paste(path, "boxplots.pdf", sep = "/"),
      paper = "a4r",
      onefile = TRUE
  )

  ## convert to characters
  data[[group]] <- as.character(data[[group]])

  data <- data[, !is.na(colnames(data))]

  ## convert to numeric
  data <- data %>%
    select_if(names(.) == paste0(group) | sapply(., is.numeric)) %>%
    mutate_all(~ ifelse(is.nan(.), NA, .)) %>%
    select_if(~ !all(is.na(.)))

  ## progrss bar
  niter <- ncol(data) - 1
  pb <- txtProgressBar(
    min = 0,
    max = niter,
    style = 3,
    char = "="
  )


  for (i in colnames(data)) {
    if (i == group) {
      next
    }
    ## stats
    plot.Dat <- data %>%
      dplyr::select(any_of(c(i, group)))

    formula <-
      as.formula(paste0("`", i, "`", "~", paste0(group)))

    stat <- plot.Dat %>%
      rstatix::pairwise_t_test(formula, p.adjust.method = "BH") %>%
      rstatix::add_xy_position(
        x = paste0(group),
        fun = "max",
        dodge = 1,
        step.increase = 0.15
      ) %>%
      dplyr::filter(p.adj.signif != "ns")

    ## plot
    p <- ggplot(
      plot.Dat,
      aes(
        x = .data[[group]],
        y = .data[[i]],
        fill = .data[[group]]
      )
    ) +
      geom_boxplot(width = 0.1) +
      geom_violin(
        alpha = 0.5,
        show.legend = FALSE
      ) +
      geom_jitter(
        shape = 21,
        size = 0.5,
        width = 0.2,
        color = "black",
        alpha = 0.3,
        show.legend = TRUE
      ) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(length(unique(
        plot.Dat[[group]]
      )), "Set1")) +
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
      theme(axis.ticks.x = element_blank()) +
      theme(
        legend.position = "bottom",
        legend.box = "vertical"
      ) +
      ggpubr::stat_pvalue_manual(
        stat,
        label = "p.adj.signif",
        tip.length = 0.01,
        color = "gray50",
        size = 6,
        bracket.size = 0.8,
        inherit.aes = FALSE
      ) +
      xlab("") +
      labs(fill = paste0(group)) +
      theme(axis.text.x = element_blank()) +
      ylab("median scaled log10 [absolute concentration]") +
      ggtitle(i) +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
    ## print
    print(p)

    Sys.sleep(0.01)

    setTxtProgressBar(pb, which(colnames(data) == i))
  }
  dev.off()

  close(pb)
}
