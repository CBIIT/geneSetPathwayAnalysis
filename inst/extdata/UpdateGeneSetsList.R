#--------------------------------------------------------------------------------------------------
# Update GENE SETS LIST mage by Sudhir and DTB team on April 11, 2022, ONLY HUGO NAMES
#--------------------------------------------------------------------------------------------------

gsets = read.csv("inst/extdata/220411.Genelist.csv", check.names = F, stringsAsFactors = F)
dim(gsets) # 6780 -23 ## removed genes without hugo names
colnames(gsets)[1]= "ZZgene"
oldcol = colnames(gsets)

gsets = gsets[,sort(oldcol)]
## check
length(gsets$ZZgene) # 6780
length(unique(gsets$ZZgene)) # 6780
nbzero = apply(gsets,1,function(x) length(which(x==0)))
table(nbzero)

# 16   17   18   19   20   21   22
#  6    9   33  224 1226 3749 1533
# so additional 1533 genes with all zeros >> 6780-1533 >> 5247 common genes
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

save(geneSets, file = "data/geneSets.RData")

# gs = geneSetPathwayAnalysis::geneSets
#
## stop here --------------------------------------------------

