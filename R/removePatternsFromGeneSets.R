#' Remove entries matching a vector of patterns from a list of vectors
#'
#' @param geneSets a list of vectors (e.g. gene sets)
#' @param patterns a vector of patterns for grepl()
#' @param verbose a boolean whether to show debugging output
#'
#' @return the gene sets (list of vectors) without entries matching patterns
#'
#' @examples
#' tmp <- list(a=c("abc", "def", "deg", "ghi"), b=c("def", "deg", "ijk"))
#' removePatternsFromGeneSets(tmp, "def", TRUE)
#'
#' @concept geneSetPathwayAnalysis
#' @export
removePatternsFromGeneSets <- function(geneSets, patterns, verbose=FALSE) {
  tmpGeneSets <- lapply(geneSets, function(x) {
    for(pattern in patterns) {
      tmp <- !grepl(pattern, x)

      if(verbose) {
        cat("REMOVED: ", length(which(!tmp)), "\n")
      }

      x <- x[tmp]
    }

    return(x)
  })

  return(tmpGeneSets)
}
