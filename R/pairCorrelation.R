#' pairCorrelation
#'
#' This function plot pair correlation
#' @param data dataframe to load
#' @keywords install, load cran and bioconductor packages 
#' @export
#' @examples
#' pairCorrelation(x)
#' 

pairCorrelation <- function (data, method, cex.labels) {
   pairs(data, 
        upper.panel = function(x, y, digits = 2, cex.cor, ...) {
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(0, 1, 0, 1))
          # correlation coefficient
          r <- cor(x, y, method=method)
          txt <- format(c(r, 0.123456789), digits = digits)[1]
          txt <- paste("r = ", txt, sep = "")
          text(0.5, 0.6, txt, cex = cex.labels)
          # p-value calculation
          p <- cor.test(x, y, method=method)$p.value
          txt2 <- format(c(p, 0.123456789), digits = digits)[1]
          txt2 <- paste("p = ", txt2, sep = "")
          text(0.5, 0.4, txt2, cex = cex.labels)
        },
        diag.panel = function(x, ...) {
          usr <- par("usr"); on.exit(par(usr))
          par(usr = c(usr[1:2], 0, 1.5) )
          h <- hist(x, plot = FALSE)
          breaks <- h$breaks; nB <- length(breaks)
          y <- h$counts; y <- y/max(y)
          rect(breaks[-nB], 0, breaks[-1], y, col = "dodgerblue", ...)
        }, 
        cex.labels = cex.labels,
        lower.panel = panel.smooth)
}
