library(RColorBrewer)
library(paxtoolsr)

##
lov <- readGmt(file.path(.lmp, "gene_set_pathway_analysis", "data", "lmp_gene_sets.txt"))
q <- c("LMP_DDR", "LMP_Apoptosis", "LMP_Oncogenes", "LMP_Tumor Suppressors", "LMP_ABC Transporters", "LMP_SLC")

lov <- lov[q]

r1 <- split(rep(names(lov), lengths(lov)), unlist(lov))
r2 <- lapply(r1, paste, collapse=";")

##
enAnnot <- read.table(file.path("inst", "extdata", "en_annot_table.txt"), sep="\t", header=TRUE, stringsAsFactors = FALSE)
nscs <- unique(enAnnot$nsc)
i <- 1
nscs[i]
cors <- enAnnot[enAnnot$nsc == nscs[i], "respCor"]
genes <- enAnnot[enAnnot$nsc == nscs[i], "gene"]

b <- 100
a <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(b)
cut(runif(10), seq(-1, 1, length.out = (b+1)), a)

nsc <- nscs[i]
idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]

idBs <- genes

netR4 <- NULL

for(i in 1:length(idBs)) {
    idB <- idBs[i]
    cat("idB", idB, "\n")
    tryCatch({
        yT1 <- getShortestPathSif(g, idA[1], idB, "out")
        netR4 <- rbind(netR4, yT1)   
    }, error = function(e) {
        cat("ERROR\n")
    })
}    

g1 <- loadSifInIgraph(netR4)

g1 <- addAttributeList(g1, "lmp", r2)

col <- rep("pink", length(V(g1)))
col[!is.na(vertex_attr(g1, "lmp"))] <- "lightblue"

g1 <- set_vertex_attr(g1, "color", value=col)

plot(g1, layout = layout.fruchterman.reingold, vertex.size=10, edge.arrow.size=0.25)

# vertex_attr(g1, "lmp")
# vertex_attr(g1, "color")

graphFile <- "~/Downloads/graphFile.graphml"
write_graph(g1, graphFile, "graphml")

###

nsc <- "681239"
idA <- chebiNscMap[chebiNscMap$NSC == nsc, "PARTICIPANT_A"]
idB <- "ACIN1"

aIdx <- match(idA, V(g)$name)
bIdx <- match(idB, V(g)$name)    

s1 <- shortest_paths(g, aIdx, bIdx, mode="out", output="vpath")
s1$vpath
cat("PATH LENGTH: ", length(s1$vpath[[1]]), "\n")
s2 <- all_shortest_paths(g, aIdx, bIdx, mode="out")
s2$res
cat("ALL PATHS: ", length(s2$res), "\n")

y <- lapply(s2$res, function(x) { 
    length(which(x$name %in% unique(unlist(lov))))    
})

tAS1 <- as.vector(unlist(y))
idx <- which(tAS1 == max(tAS1))
s2$res[idx]





