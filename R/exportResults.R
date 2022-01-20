#' exportResults
#'
#'export results with precolored formatting
#'
#' @param data log fold changes data generated from normalizeDat function
#' @param path saving path of the data
#'
#' @import tidyverse
#' @import here
#' @import condformat
exportResults <- function(data, path=NULL) {
  ## progrss bar
  niter <- length(unique(data$contrast))
  pb <- txtProgressBar(min = 0,
                       max = niter,
                       style = 3,
                       char = "=")

  for(j in 1:niter) {

    if(is.null(path)) {
      path = here()
    } else
      path = path
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
    keys <- condformat(keys) %>%
      rule_fill_discrete(
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
    condformat2excel(keys,
                     paste(path,
                           "results.xlsx", sep = "/"),
                     sheet_name = "legends")

    ## create result tables
    groups <- unique(data$contrast)

    for (i in groups) {
      ## subset dataframe
      csv <- data[data$contrast %in% i, ]


      ## file name
      ws <- gsub(" ", "_", i)

      if (nchar(ws) > 31) {
        ws <- substr(ws, start = 1, stop = 31)
      }
      ## export results
      csvFormat <- condformat(csv) %>%
        rule_fill_discrete(
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
        rule_text_bold(logFC,
                       expression = logFC < 0.5 |
                         logFC > 2) %>%
        rule_fill_discrete(
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
      condformat2excel(
        csvFormat,
        paste(
          path,
          "results.xlsx",
          sep = "/"
        ),
        sheet_name = ws
      )


    }
    Sys.sleep(0.01)
    setTxtProgressBar(pb, j)
  }

  close(pb)
  cat("\n")
  return(paste("results saved in", path, "folder"))
}
