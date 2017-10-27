#--------------------------------------------------------------------------------------------------
# PREPARE AND SAVE GENE SETS LIST
#--------------------------------------------------------------------------------------------------
library(stringr)
geneAnnotTab <- read.table(file="inst/extdata/cellminer_gene_annot.txt", header=TRUE, sep="\t",
  stringsAsFactors=FALSE, comment.char="", quote="", na.strings = "-")

geneAnnotTab <- geneAnnotTab[!is.na(geneAnnotTab$ANNOT), ]
stopifnot(all(!duplicated(geneAnnotTab$GENE_NAME)))

geneAnnotTab <- geneAnnotTab[order(geneAnnotTab$GENE_NAME), ]
rownames(geneAnnotTab) <- geneAnnotTab$GENE_NAME

trimAnnot <- function(annotStr){
  tmp <- c(str_split(annotStr, pattern = ";"), recursive = TRUE)
  trimmedAnnot <- tmp[1]
  if (length(tmp) > 1){
    trimmedAnnot <- paste0(trimmedAnnot, ";", tmp[2])
  }
  if (str_sub(trimmedAnnot, start = 1, end = 1) == "\""){
    trimmedAnnot <- str_sub(trimmedAnnot, start = 2)
  }
  return(trimmedAnnot)
}

geneAnnotTab$SHORT_ANNOT <- vapply(geneAnnotTab$ANNOT, trimAnnot, character(1))

save(geneAnnotTab, file = "data/geneAnnotTab.RData")
#--------------------------------------------------------------------------------------------------
