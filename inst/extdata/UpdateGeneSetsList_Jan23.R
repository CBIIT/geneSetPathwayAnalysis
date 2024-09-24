#--------------------------------------------------------------------------------------------------
# Update GENE SETS LIST made by Sudhir and DTB team on April 11, 2022 plus
# cell surface gene set
#--------------------------------------------------------------------------------------------------

gsets = read.csv("inst/extdata/230131.Genelist.csv", check.names = F, stringsAsFactors = F)
dim(gsets) # 9566   24 ## removed genes without hugo names
colnames(gsets)[1]= "ZZgene"
oldcol = colnames(gsets)

gsets = gsets[,sort(oldcol)]
## check
length(gsets$ZZgene) # 9566
length(unique(gsets$ZZgene)) # 9566
nbzero = apply(gsets,1,function(x) length(which(x==0)))
table(nbzero)

# 16   17   18   19   20   21   22   23
# 1    5   10   90  518 1650 6389  903
## 903 genes to discard with all zeros
##

geneSets <- list()
geneSets[["All Gene Sets"]] <- ""

for (k in 1:(ncol(gsets)-1)) {
  idxk = which(gsets[,k]==1)

  geneSets[[colnames(gsets)[k]]] <- gsets$ZZgene[idxk]

}

geneSets[["All Gene Sets"]] <- setdiff(sort(unique(c(geneSets, recursive = TRUE))), "")

geneSets <- lapply(geneSets, FUN = unique)

stopifnot(all(sort(vapply(geneSets, length, integer(1)))) > 0)
## 8663 genes
save(geneSets, file = "data/geneSets.RData")

# gs = geneSetPathwayAnalysis::geneSets
#
## stop here --------------------------------------------------

