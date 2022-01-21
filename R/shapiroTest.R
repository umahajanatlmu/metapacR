#' shapiroTest
#'
#' Shapiro-Wilk normality test
#'
#' @param data metabolome raw data expDataList
#'
#' @import tidyverse
#' @import utils
#' @import stats
shapiroTest <- function(dataList=dataList) {

  data <- dataList[["imputed.matrix"]]

  result <- data.frame()

  for (y in var) {
    if (class(data[[y]]) != "numeric") {
      next
    } else if (length(na.omit(data[[y]])) < 5 ||
      sum(is.nan(data[[y]])) == nrow(data)) {
    next
  } else
    shapiroMet <- shapiro.test(as.numeric(as.vector(data[[y]])))

  shapiroMet <- do.call(cbind, shapiroMet) %>%
    as.data.frame() %>%
    mutate(Metabolite = y) %>%
    select(-data.name, -method) %>%
    remove_rownames()

  result <- bind_rows(result, shapiroMet)

  }

  return(result)

}
