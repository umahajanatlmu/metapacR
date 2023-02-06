#' @title leveneStat
#'
#' @description computes Levene's stats for homogeneity of variance across groups
#'
#' @param group grouping variable
#' @param dataList  metabolome raw data expDataList
#' @param location mean or median, median is default
#' @param trim.alpha specify the percentage of trimmed mean
#' @param bootstrap whether to perform bootstrapping TRUE/FALSE
#' @param num.bootstrap number of bootstraping steps
#' @param kruskal.test whether to perform kruskal test TRUE/FALSE
#'
#' @import tidyverse
#' @import utils
#' @importFrom car leveneTest
#' @import stats
#'
#' @return result table
#'
#' @export

leveneStat <- function(group,
                       dataList,
                       location = c("mean", "median"),
                       trim.alpha = 0.25,
                       bootstrap = TRUE,
                       num.bootstrap = 1000,
                       kruskal.test = TRUE) {
  stopifnot(inherits(dataList, "list"))
  validObject(dataList)

  location <- match.arg(location)

  if (is.null(group)) {
    stop("group variable is missing")
  } else if (length(group) != 1) {
    stop("multiple group variables available....provide only one group variable")
  }

  # load imputed data matrix
  ## ----------------------------------------------------------------
  imputed.data <- dataList[["imputed.matrix"]]

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

  ## select variables

  var <- setdiff(colnames(data), group)

  result <- data.frame()

  for (y in var) {
    if (class(data[[y]]) != "numeric") {
      next
    } else {
      y.var <- na.omit(as.numeric(as.vector(data[[y]])))
    }
    group.var <- data[[group]][!is.na(data[[y]])]
    nn <- length(y.var)
    lower <- ceiling(nn * trim.alpha) + 1
    upper <- floor(nn * (1 - trim.alpha))

    if (lower > upper) {
      print(paste("frac.trim.alpha value is too large", y))
    } else {
      leveneMet <- car::leveneTest(
        y = y.var,
        group = group.var,
        location = location,
        trim.alpha = trim.alpha,
        bootstrap = bootstrap,
        num.bootstrap = num.bootstrap,
        kruskal.test = kruskal.test
      )
    }

    leveneMet <- leveneMet %>%
      as.data.frame() %>%
      drop_na() %>%
      mutate(Metabolite = y) %>%
      dplyr::select(Metabolite, `F value`, `Pr(>F)`)

    row.names(leveneMet) <- NULL

    result <- bind_rows(result, leveneMet)
  }
  return(result)
}
