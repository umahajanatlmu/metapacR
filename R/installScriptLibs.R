#' @title installScriptLibs
#'
#' @description function to install and load listed R-packages from CRAN and Bioconductor
#'
#' @param packages list of required packages
#'
#' @import utils
#' @importFrom BiocManager install
#' @importFrom BiocManager available
#'
#' @export

installScriptLibs <- function(packages) {

  stopifnot(inherits(packages, "character"))
  validObject(packages)

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
    } else if (i %in% available(i)) {
      print(paste("installing:",i))
      ## install Bioconductor packages
      install(i)
      ## load library
      library(i, character.only = TRUE)
    } else
      print(paste("package", i, "not found in CRAN or Bioconductor"))
  }
  options(warn=0) ## reset all warnings
}
