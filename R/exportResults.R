#' exportResults
#'
#' export results with precolored formatting
#'
#' @param data log fold changes data generated from normalizeDat function
#' @param path saving path of the data
#' @param filename name of the file
#'
#' @import tidyverse
#' @importFrom here here
#' @importFrom condformat condformat rule_fill_discrete rule_text_bold condformat2excel
#'
#' @return color coded result table in excel formant
#'
#' @export

exportResults <- function(data, path = NULL, filename = NULL) {
  stopifnot(inherits(data, "data.frame"))
  validObject(dataList)

  ## create result tables groups
  groups <- unique(data$contrast)

  ## progrss bar
  niter <- length(groups)
  pb <- txtProgressBar(
    min = 0,
    max = niter,
    style = 3,
    char = "="
  )

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

  if (is.null(filename)) {
    filename <- "results"
  }
  ## create keys
  keys <- data.frame(
    legends = c(
      "ratio > 1 and p < 0.01",
      "ratio > 1 and p < 0.05",
      "ratio > 1 and p < 0.1",
      "ratio < 1 and p < 0.01",
      "ratio < 1 and p < 0.05",
      "ratio < 1 and p < 0.1",
      "ratio < 0.5 or ratio > 2 = bold",
      "0.05 <= p < 0.10",
      "0.01 <= p <0.05",
      "p < 0.01"
    )
  )
  ## format keys
  keys <- condformat::condformat(keys) %>%
    condformat::rule_fill_discrete(
      legends,
      colours = c(
        "ratio > 1 and p < 0.01" = "#CC0033",
        "ratio > 1 and p < 0.05" = "#FF3366",
        "ratio > 1 and p < 0.1" = "#FF6666",
        "ratio < 1 and p < 0.01" = "#3333FF",
        "ratio < 1 and p < 0.05" = "#3399FF",
        "ratio < 1 and p < 0.1" = "#99CCFF",
        "0.05 <= p < 0.10" = "#CCCCCC",
        "0.01 <= p <0.05" = "#FFFF99",
        "p < 0.01" = "#FFFF66"
      )
    )
  # export legends
  condformat::condformat2excel(keys,
    paste(path,
      "results.xlsx",
      sep = "/"
    ),
    sheet_name = "legends"
  )


  for (i in groups) {
    ## subset dataframe
    csv <- data[data$contrast %in% i, ]


    ## file name
    ws <- gsub(" ", "_", i)

    if (nchar(ws) > 31) {
      ws <- substr(ws, start = 1, stop = 31)
    }
    ## export results
    csvFormat <- condformat::condformat(csv) %>%
      condformat::rule_fill_discrete(
        logFC,
        expression = ifelse(
          logFC > 1 &
            adj.P.Val < 0.1 & adj.P.Val >= 0.05,
          "high",
          ifelse(
            logFC > 1 &
              adj.P.Val < 0.05 & adj.P.Val >= 0.01,
            "vhigh",
            ifelse(
              logFC > 1 & adj.P.Val < 0.01,
              "vvhigh",
              ifelse(
                logFC < 1 & adj.P.Val < 0.1 & adj.P.Val >= 0.05,
                "low",
                ifelse(
                  logFC < 1 &
                    adj.P.Val < 0.05 & adj.P.Val >= 0.01,
                  "vlow",
                  ifelse(logFC < 1 &
                    adj.P.Val < 0.1, "vvlow", "")
                )
              )
            )
          )
        ),
        colours = c(
          "high" = "#FF6666",
          "vhigh" = "#FF3366",
          "vvhigh" = "#CC0033",
          "low" = "#99CCFF",
          "vlow" = "#3399FF",
          "vvlow" = "#3333FF"
        )
      ) %>%
      condformat::rule_text_bold(logFC,
        expression = logFC < 0.5 |
          logFC > 2
      ) %>%
      condformat::rule_fill_discrete(
        adj.P.Val,
        expression = ifelse(
          adj.P.Val < 0.1 & adj.P.Val >= 0.05,
          "*",
          ifelse(
            adj.P.Val < 0.05 & adj.P.Val >= 0.01,
            "**",
            ifelse(adj.P.Val < 0.01, "***", "")
          )
        ),
        colours = c(
          "*" = "#CCCCCC",
          "**" = "#FFFF99",
          "***" = "#FFFF66"
        )
      )

    ## compile results
    condformat::condformat2excel(
      csvFormat,
      paste(
        path,
        paste0(filename, ".xlsx"),
        sep = "/"
      ),
      sheet_name = ws
    )

    Sys.sleep(0.01)
    setTxtProgressBar(pb, which(groups == i))
  }

  close(pb)
  cat("\n")
  return(paste("results saved in", path, "folder"))
}
