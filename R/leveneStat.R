#' leveneStat
#'
#' Computes Levene's test for homogeneity of variance across groups.
#' @param group grouing variable
#' @param data  dataframe
#' @param location mean or median, median is default
#' @param trim.alpha specify the % of trimmed mean
#' @param bootstrap whether to perform bootstrapping TRUE/FALSE
#' @param num.bootstrap bootstraping steps
#' @param kruskal.test whether to perform kruskal test TRUE/FALSE
#'
#' @import tidyverse
#' @import utils
#' @import car
#' @import stats
#'
#' @return
#' @export
#'
#' @examples
#' leveneStat(group = "group", data = data, trim.alpha = 0.25)
#'
leveneStat <- function(group = group,
                       data = data,
                       location = "median",
                       trim.alpha = 0.25,
                       bootstrap = TRUE,
                       num.bootstrap = 1000,
                       kruskal.test = TRUE) {

  var <- setdiff(colnames(data), group)

  result <- data.frame()

  for (y in var) {
    if (class(data[[y]]) != "numeric") {
      next
    } else
  y.var <- na.omit(as.numeric(as.vector(data[[y]])))
  group.var <- data[[group]][!is.na(data[[y]])]
  nn <- length(y.var)
  lower <- ceiling(nn * trim.alpha) + 1
  upper <- floor(nn * (1 - trim.alpha))

  if (lower > upper) {
    print(paste("frac.trim.alpha value is too large", y))
  } else
    leveneMet <- leveneTest(
      y = y.var,
      group = group.var,
      location = "median",
      trim.alpha = trim.alpha,
      bootstrap = bootstrap,
      num.bootstrap = num.bootstrap,
      kruskal.test = kruskal.test
    )

  leveneMet <- leveneMet %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(Metabolite = y,
           Group = group)

  result <- bind_rows(result, leveneMet)
}
  return(result)
}

