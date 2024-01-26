#' @title normalizeDat.binary
#'
#' @description normalize data using fixed effect as well as mixed effects models for comparison of binary groups ie group vs rest.
#'
#' @param dataList raw metabolome data list from imputeTransformScale function.It need to have imputed.matrix and metadata.
#' @param confounders list of confounders, NULL for mixed effect model
#' @param stratifier classifier variable of interest, Should not be in confounders. Should be one length 1.
#' @param fix.effect equation of fixed effect, NULL for fixed effect model
#' @param random.effect equation of random effects, NULL for fixed effect model
#'
#' @import tidyverse
#' @import utils
#' @importFrom emmeans emmeans ref_grid
#' @importFrom nlme lme
#' @import stats
#' @import graphics
#'
#' @return Analyses results in list object.
#'   The object contains the following:\itemize{
#'     \item fitted.value fitted values for each variable
#'     \item summaryFC anova results
#'   }
#'
#' @export

normalizeDat.binary <- function(dataList,
                                confounders = NULL,
                                stratifier,
                                fix.effect = NULL,
                                random.effect = NULL) {
  stopifnot(inherits(dataList, "list"))
  validObject(dataList)

  if (is.null(stratifier)) {
    stop("stratifier variable is missing")
  } else if (length(stratifier) != 1) {
    stop("multiple stratifier variables available....provide only one stratifier variable")
  }

  ## load imputed data matrix
  ## ----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

  ## load metadata
  ## ----------------------------------------------------------------
  metadata.data <- dataList[["metadata"]]

  ## subset metadata
  ## ----------------------------------------------------------------
  select.columns <- c(confounders, stratifier)
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


  ## define met, empty data-frames
  ## -----------------------------------------------------------------
  met <- setdiff(colnames(data), c(confounders, stratifier))
  normData <- data
  summaryFC <- data.frame()

  ## define column order
  ## -----------------------------------------------------------------
  column.order <- c(
    "contrast",
    "Metabolite",
    "logFC",
    "OR",
    "estimate",
    "SE",
    "df",
    "t.ratio",
    "p.value",
    "adj.P.Val"
  )

  ## unique groups
  groups <- unique(data[[stratifier]])

  pb <- txtProgressBar(
    min = 0,
    max = length(groups),
    style = 3,
    char = "="
  )


  for (g in groups) {
    data[["bin.stratifier"]] <- ifelse(data[[stratifier]] == g, g, "rest")

    data[["bin.stratifier"]] <- as.factor(data[["bin.stratifier"]])

    data[["bin.stratifier"]] <- relevel(data[["bin.stratifier"]], ref = g)


    for (i in met) {
      response <- i

      if (class(data[[response]]) %in% c("character", "factor")) {
        next
      } else if (sum(is.nan(data[[response]])) == nrow(data) |
        sum(is.na(data[[response]])) == nrow(data)) {
        next
      }

      if (is.null(fix.effect) && is.null(random.effect)) {
        ## dependent variables
        ## -----------------------------------------------------------------
        if (is.null(confounders)) {
          indepVars <- "bin.stratifier"
        } else {
          indepVars <-
            paste(
              c(
                confounders,
                "bin.stratifier"
              ),
              collapse = " + ",
              sep = ""
            )
        }
        ## formula
        ## -----------------------------------------------------------------
        formula <-
          as.formula(paste(paste("`", response, "`", "~", sep = ""), indepVars))

        ## perform modeling
        ## -----------------------------------------------------------------
        model <- lm(formula,
          data = data,
          na.action = na.exclude
        )

        ## perform modeling for annova
        ## -----------------------------------------------------------------
        anova.grid <- ref_grid(model)
        anova.emmeans <- emmeans::emmeans(model, "bin.stratifier")
        anova.model <- pairs(anova.emmeans)
        anova.results <- anova.model %>%
          as.data.frame() %>%
          mutate(
            Metabolite = response,
            logFC = 10^estimate,
            OR = exp(estimate)
          )

        ## perform FDR correction
        ## -----------------------------------------------------------------
        anova.results.fdr <- update(anova.model, adjust = "BH") %>%
          as.data.frame() %>%
          select(contrast, p.value)

        colnames(anova.results.fdr)[colnames(anova.results.fdr) == "p.value"] <- "adj.P.Val"

        anova.results.fdr <- anova.results.fdr %>%
          full_join(anova.results, by = "contrast") %>%
          select(all_of(column.order))

        ## fold changes results
        ## -----------------------------------------------------------------
        summaryFC <- bind_rows(summaryFC, anova.results.fdr)

        ## fit the model
        ## -----------------------------------------------------------------
        normData[, response] <- fitted.values(model)
      } else if (!is.null(fix.effect)) {
        print("fix effect variable is missing")
        break
      } else if (!is.null(random.effect)) {
        print("fix effect variable is missing")
        break
      } else {
        ## formula
        ## -----------------------------------------------------------------
        formula <-
          as.formula(paste(paste("`", response, "`", "~", sep = ""), fix.effect))

        ## perform modeling
        ## -----------------------------------------------------------------\
        model <- nlme::lme(formula,
          random = random.effect,
          data = data,
          na.action = na.exclude
        )

        ## perform modeling for annova
        ## -----------------------------------------------------------------
        anova.grid <- emmeans::ref_grid(model)
        anova.emmeans <- emmeans::emmeans(model, stratifier)
        anova.model <- pairs(anova.emmeans)
        anova.results <- anova.model %>%
          as.data.frame() %>%
          mutate(
            Metabolite = response,
            logFC = 10^estimate,
            OR = exp(estimate)
          )

        ## perform FDR correction
        ## -----------------------------------------------------------------
        anova.results.fdr <- update(anova.model, adjust = "BH") %>%
          as.data.frame() %>%
          select(contrast, p.value)

        colnames(anova.results.fdr)[colnames(anova.results.fdr) == "p.value"] <- "adj.P.Val"

        anova.results.fdr <- anova.results.fdr %>%
          full_join(anova.results, by = "contrast") %>%
          select(all_of(column.order))

        ## fold changes results
        ## -----------------------------------------------------------------
        summaryFC <- bind_rows(summaryFC, anova.results.fdr)

        ## fit the model
        ## -----------------------------------------------------------------
        normData[, response] <- fitted.values(model)
      }
    }

    Sys.sleep(0.01)
    setTxtProgressBar(pb, which(groups == g))
  }

  close(pb)

  return(list(
    fitted.values = normData,
    summaryFC = summaryFC
  ))
}
