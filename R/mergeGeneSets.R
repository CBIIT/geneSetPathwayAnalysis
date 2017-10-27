#' Merge Gene Sets
#'
#' @param geneSets a list of vectors (e.g. gene sets)
#' @param name the name for the merged gene set
#' @param idConversion convert IDs from Uniprot to HGNC (Default: FALSE)
#' @param returnCounts changes return to a list where one entry is "counts" (Default: FALSE)
#'
#' @return a list of vectors with one entry (i.e. a gene set) with the given
#'   name
#'
#' @examples
#' tmp <- list(a=c("abc", "def", "deg", "ghi"), b=c("def", "deg", "ijk"))
#' mergeGeneSets(tmp, "y")
#'
#' @concept geneSetPathwayAnalysis
#' @export
mergeGeneSets <- function(geneSets, name="merged", idConversion=FALSE, returnCounts=FALSE, forceCache=FALSE) {
  if(forceCache) {
    require(simpleRCache)
    setCacheRootPath()

    convertIdsCached <- addMemoization(convertIds)
  }

  tmpGeneSet <- list()
  t1 <- unname(unlist(geneSets))

  # From paxtoolsr
  if(idConversion) {
    t1 <- convertIds(t1)
  }

  t2 <- unique(t1)

  tmpGeneSet[[name]] <- sort(t2[!is.na(t2)])

  if(returnCounts) {
    results <- list(counts=sort(table(t1), decreasing = TRUE), geneSet=tmpGeneSet)
  } else {
    results <- tmpGeneSet
  }

  return(results)
}
