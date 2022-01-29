#' @title banner
#'
#' @description function to print tidy labels and separators
#'
#' @param txt text to print in banner
#' @param symbol character separator (eg. "-", "=", "*", ....)
#'
#' @import utils
#'
#' @return tidy labels.
#'
#' @export

banner <-function(txt,
                  symbol ="-") {

  stopifnot(inherits(txt, "character"))
  validObject(txt)

  nchar <- 64
  ## head tail
  headTail <- strrep(symbol, nchar)
  hash <- paste0("##")
  charEnd <- strrep(symbol,2)
  headTailSent <- paste0(hash, headTail)
  ## text
  textChar <- nchar(txt)

  if (textChar > nchar) {
    txt <- substring(txt, first = 1, last = 60)
    textChar <- nchar(txt)
  }
  space <- " "
  centering <- (nchar - 2) - textChar
  spacesBeforeAfter <- strrep(space, centering/2)

  textSent <- paste0(hash,
                     spacesBeforeAfter,
                     txt,
                     spacesBeforeAfter,
                     charEnd)

  return(cat(paste0(headTailSent ,
                    "\n",textSent ,
                    "\n",headTailSent,
                    "\n")))
}
