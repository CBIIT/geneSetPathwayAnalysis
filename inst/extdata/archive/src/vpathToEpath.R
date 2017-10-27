# lapply(all_shortest_paths(g, aIdx, bIdx)$res, function(x) {
#     cat("IT: ", get.edge.attribute(g, "interactionType", E(g, path=x)), "\n")
#     E(g, path=x)
# })

## Get Network
library(paxtoolsr)
library(rcellminer)
library(igraph)

sif <- downloadPc2("PathwayCommons.8.All.EXTENDED_BINARY_SIF.hgnc.txt.gz")
#sif <- downloadPc2("PathwayCommons.8.drugbank.BINARY_SIF.hgnc.txt.sif.gz")

x <- getSifInteractionCategories()
excludeInteractions <- c("in-complex-with")
sifX <- filterSif(sif$edges, setdiff(c(x$BetweenProteins, x$BetweenProteinSmallMolecule), excludeInteractions))

ignoreDb <- "CTD"
t1 <- lapply(sifX$INTERACTION_DATA_SOURCE, function(x) {
    if(ignoreDb %in% x && length(x) == 1) {
        return(FALSE)
    } else {
        return(TRUE)
    }
})

t2 <- unlist(t1)

sifY <- sifX[t2, ]

g1 <- loadSifInIgraph(sifY, directed=TRUE)
g <- loadSifInIgraph(sifY, directed=FALSE)

idA <- "GEN1"
idB <- "CCND1"

aIdx <- match(idA, V(g)$name)
bIdx <- match(idB, V(g)$name)    

mode <- "all"
weights <- NULL

targetCor <- "609699"

db <- list()
db[["exp"]] <- getAllFeatureData(rcellminerData::molData)[["exp"]]
db[["act"]] <- exprs(getAct(rcellminerData::drugData))

m1 <- getShortestPathSif(g1, idA, idB, mode=mode, weights=weights, filterFun=filterFun, db, targetCor)

# s2 <- all_shortest_paths(g, aIdx, bIdx, mode=mode, weights=weights)
# s2$res
# 
# s2$res[[1]]$name
# 
# vpaths <- s2$res
# 
# i <- 1
# j <- 1 

filterFun <- function(vpaths, db, targetCor) {
    results <- NULL
    pathLength <- length(vpaths[[1]]$name)

    for(i in 1:length(vpaths)) {
        path <- vpaths[[i]]$name

        tmpResults <- 0

        if(any(grepl("^CHEBI", path))) {
            next
        }

        for(j in 1:length(path)) {
            id <- path[j]

            #cat("ID: ", id, "\n")

            if(id %in% rownames(db[["exp"]])) {
                x <- db[["act"]][targetCor,]
                y <- db[["exp"]][id,]
                t1 <- cor.test(x, y)

                tmpResults <- tmpResults + abs(t1$estimate)
            }
        }

        results <- c(results, tmpResults)
    }

    n1 <- lapply(vpaths, function(x) {
        paste(names(x), collapse=":")
    })

    n1 <- unlist(n1)

    cat("CORR: ", paste(results/pathLength, collapse=", "), "\n")
    cat("PATHS: ", paste(n1, collapse=", "), "\n")

    vpath <- vpaths[[which.max(results)]]

    return(vpath)
}

# y <- filterFun(vpaths, db, targetCor)
# e <- E(g, path=y) # Works for undirected
# e
# e$interactionType
# 
# gX <- make_ring(10, directed=TRUE)
# V(gX)$name <- letters[1:10]
# yX <- all_shortest_paths(gX, 1, 6)$res[[1]]
# E(gX, path=yX)
# 
# are.connected(g1, "RB1", "GEN1") # T
# are.connected(g1, "TP53", "RB1") # T
# 
# are.connected(g1, "GEN1", "RB1") # F
# are.connected(g1, "RB1", "TP53") # F
# 
# # Evidence it has issues if the path is not ordered corrected
# E(g1, P=c("RB1", "GEN1"))
# #+ 1/381041 edge (vertex names):
# #    [1] RB1->GEN1
# E(g1, P=c("GEN1", "RB1"))
# #Error in .Call("R_igraph_es_pairs", graph, as.igraph.vs(graph, P) - 1,  : 
# #                   At type_indexededgelist.c:1173 : Cannot get edge id, no such edge, Invalid value
# 
# are.connected(g1, "RELA", "GEN1") # T
# are.connected(g1, "TP53", "RELA") # F
# 
# are.connected(g1, "GEN1", "RELA") # F
# are.connected(g1, "RELA", "TP53") # T
# 
# r1 <- E(g1, P=c("RB1", "GEN1"))
# r2 <- E(g1, P=c("RELA", "TP53"))
# 
# E(g1)[c(r1, r2)]
# 
# v1 <- y$name
# v1 <- s2$res[[5]]$name
# e1 <- NULL 
# 
# for(i in 1:(length(v1)-1)) {
#     if(are.connected(g1, v1[i], v1[i+1])) {
#         r1 <- E(g1, P=c(v1[i], v1[i+1]))
#     } else {
#         r1 <- E(g1, P=c(v1[i+1], v1[i]))
#     }
#     
#     e1 <- c(e1, r1)
# }
# 
# E(g1)[e1]
# E(g1)[e1]$interactionType
# 
# m1 <- getShortestPathSif(g1, idA, idB, mode="all", weights=NULL, filterFun=filterFun, db, targetCor)
