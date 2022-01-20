#' banner
#'
#' print banner as header with separator
#'
#' @param txt text to print in banner
#' @param Char character seperator (eg. "-", "=", "*")
#'
#' @import utils
banner <-function(txt, Char ="-") {
  nchar <- 64
  ## head tail
  headTail <- strrep(Char, nchar)
  hash <- paste0("##")
  charEnd <- strrep(Char,2)
  headTailSent <- paste0(hash, headTail)
  ## text
  textChar <- nchar(txt)
  if (textChar > nchar) {
    txt <- substring(txt, first = 1, last = 60)
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
