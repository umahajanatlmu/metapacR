#' normalizeDat
#'
#' Function to normalize data using fixed effect as well as mixed effects models.
#'
#' @param data metabolome raw data
#' @param confounders list of confounders, NULL for mixed effect model
#' @param stratifier classifier variable of interest
#' @param fix.effect equation of fixed effect, NULL for fixed effect model
#' @param random.effect equation of random effects, NULL for fixed effect model
#'
#' @import tidyverse
#' @import utils
#' @import emmeans
#' @import nlme
#' @import stats
normalizeDat <- function (data = data,
                          confounders = NULL,
                          stratifier = stratifier,
                          fix.effect = NULL,
                          random.effect = NULL) {
  ## define met, empty data-frames
  ##-----------------------------------------------------------------
  met <- setdiff(colnames(data), confounders)
  normData <- data
  summaryFC <- data.frame()

  ## define column order
  ##-----------------------------------------------------------------
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

  niter <- ncol(data)
  pb <- txtProgressBar(min = 0,
                       max = niter,
                       style = 3,
                       char = "=")

  for(j in 1:niter) {

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
      ##-----------------------------------------------------------------
      if (is.null(confounders)) {
        indepVars <-
          paste0(stratifier)
      } else
        indepVars <-
          paste(c(confounders,
                stratifier),
                collapse = " + ",
                sep = "")
      ## formula
      ##-----------------------------------------------------------------
      formula <-
        as.formula(paste(paste("`", response, "`", "~", sep = ""), indepVars))

      ## perform modeling
      ##-----------------------------------------------------------------
      model <- lm(formula,
                  data = data,
                  na.action = na.exclude)

      ## perform modeling for annova
      ##-----------------------------------------------------------------
      anova.grid <- ref_grid(model)
      anova.emmeans <- emmeans(model, stratifier)
      anova.model <- pairs(anova.emmeans)
      anova.results <- anova.model %>%
        as.data.frame() %>%
        mutate(Metabolite = response,
               logFC = 10 ^ estimate,
               OR = exp(estimate))

      ## perform FDR correction
      ##-----------------------------------------------------------------
      anova.results.fdr <- update(anova.model, adjust = "BH") %>%
        as.data.frame() %>%
        select(contrast, p.value) %>%
        rename(adj.P.Val = p.value) %>%
        full_join(anova.results, by = "contrast") %>%
        select(all_of(column.order))

      ## fold changes results
      ##-----------------------------------------------------------------
      summaryFC <- bind_rows(summaryFC, anova.results.fdr)

      ## fit the model
      ##-----------------------------------------------------------------
      normData[, response] <- fitted.values(model)
    } else if (!is.null(fix.effect)) {
      print("fix effect variable is missing")
      break
    } else if (!is.null(random.effect)) {
      print("fix effect variable is missing")
      break
    } else {
      ## formula
      ##-----------------------------------------------------------------
      formula <-
        as.formula(paste(paste("`", response, "`", "~", sep = ""), fix.effect))

      ## perform modeling
      ##-----------------------------------------------------------------\
      model <- lme(formula,
                   random = random.effect,
                   data = data,
                   na.action = na.exclude)

      ## perform modeling for annova
      ##-----------------------------------------------------------------
      anova.grid <- ref_grid(model)
      anova.emmeans <- emmeans(model, stratifier)
      anova.model <- pairs(anova.emmeans)
      anova.results <- anova.model %>%
        as.data.frame() %>%
        mutate(Metabolite = response,
               logFC = 10 ^ estimate,
               OR = exp(estimate))

      ## perform FDR correction
      ##-----------------------------------------------------------------
      anova.results.fdr <- update(anova.model, adjust = "BH") %>%
        as.data.frame() %>%
        select(contrast, p.value) %>%
        rename(adj.P.Val = p.value) %>%
        full_join(anova.results, by = "contrast") %>%
        select(all_of(column.order))

      ## fold changes results
      ##-----------------------------------------------------------------
      summaryFC <- bind_rows(summaryFC, anova.results.fdr)

      ## fit the model
      ##-----------------------------------------------------------------
      normData[, response] <- fitted.values(model)

    }

  }
    Sys.sleep(0.01)
    setTxtProgressBar(pb, j)
  }

  close(pb)

  return(list(fitted.values = normData,
              summaryFC = summaryFC))

}
