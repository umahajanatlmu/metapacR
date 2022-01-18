#' installScriptLibs
#'
#' function to install and load required listed packages
#' @param packages list of required packages
#'
#' @import utils
#' @import BiocManager
#'
#' @return
#' @export
#'
#' @examples
#' pkgs <- c("ggplot2", "ggrepel", "arsenal", "limma")
#' installScriptLibs(packages=pkgs)
#'
installScriptLibs <- function(packages) {
  options(warn=-1) ## supress all warning
  for (i in packages) {
    if (i %in% .packages(all.available = TRUE)) {
      ## load library
      library(i, character.only = TRUE)
    } else if (i %in% names(available.packages()[,1])) {
      print(paste("installing:",i))
      ## install CRAN packages
      install.packages(i, character.only = TRUE)
      ## load library
      library(i, character.only = TRUE)
    } else if (i %in% BiocManager::available(i)) {
      print(paste("installing:",i))
      ## install Bioconductor packages
      BiocManager::install(i)
      ## load library
      library(i, character.only = TRUE)
    } else
      print(paste("package", i, "not found in CRAN or Bioconductor"))
  }
  options(warn=0) ## reset all warnings
}
