#' Make a Tree List
#'
#' @param df a data.frame
#' @param root the entry in the data.frame that serves at the tree root node
#'
#' @return a list
#'
#' @concept geneSetPathwayAnalysis
#' @export
makeTreeList <- function(df, root=df[1,1]) {
  if(is.factor(root)) root <- as.character(root)
  r <- list(name=root)
  children = df[df[,1]==root,2]
  if(is.factor(children)) children <- as.character(children)
  if(length(children)>0) {
    for(i in 1:length(children)) {
      child <- children[i]
      r[[child]] <- makeTreeList(df=df, root=child)
    }
  }
  r
}
