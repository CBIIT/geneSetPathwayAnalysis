library(igraph)
library(sqldf)
library(RCytoscape)

# LOAD FUNCTIONS ----
source(file.path(.lmp, "gene_set_pathway_analysis", "src", "getPathwayOverrepresentation.R"))
source(file.path(.lmp, "ElasticNetExps", "src", "ElasticNet.R"))
source(file.path(.lmp, "Utils", "src", "sourceDir.R"))
sourceDir(path = file.path(.lmp, "Utils", "src"))

# LOAD DATA ----
load(file.path(.lmp, "RData", "lmpdb.Rdata"))

#interactions <- read.table(file.path(.lmp, "gene_set_pathway_analysis", "data", "PC_sif_102813.gz"), sep=" ", as.is=TRUE)
#interactions <- read.table(file.path(.lmp, "gene_set_pathway_analysis", "data", "Pathway\ Commons.5.All.BINARY_SIF.tsv.gz"), sep="\t", quote="", as.is=TRUE)
filename <- file.path(.lmp, "gene_set_pathway_analysis", "data", "pc5_extended_sif_072414_edges.zip")
interactions <- read.table(unz(filename, "edges.tsv"), sep="\t", quote="", as.is=TRUE, header=TRUE)

load(file.path(Sys.getenv("HOME"), "nsc_290193.Rdata"))

# ADD CORRELATED PREDICTORS ----
elNetResults <- addPredictorCorrelatedFeatures(elNetResults, elNetMolDataNCI60, corThreshold = 0.8)

# GENERATE SEED SET ----
#seedSet <- substr(names(elNetResults$predictorWts), 4, 100)
seedSet <- substr(elNetResults$predictorsAndCorFeatures, 4, 100)

userRemoveSet <- NULL
#userRemoveSet <- c("SAR1B", "SEC24A")

# FILTER NETWORK ----
# NOTE: Use with PC_sif_102813.gz
# filteredInteractions <- sqldf("SELECT V1, V2, V3
# 	FROM interactions
# 	WHERE V1 != V3
# 	AND V2 != 'GENERIC_OF'
# 	AND V2 != 'IN_SAME_COMPONENT'
# 	AND V2 != 'INTERACTS_WITH'")

# NOTE: Use with Pathway\ Commons.5.All.BINARY_SIF.tsv.gz
filteredInteractions <- sqldf("SELECT V1, V2, V3
	FROM interactions
	WHERE V1 != V3
	AND V2 != 'neighbor-of'
	AND V2 != 'chemical-affects'
	AND V2 != 'in-complex-with'  
  AND V2 != 'interacts-with'
	AND V2 != 'controls-transport-of-chemical'
	AND V2 != 'consumption-controlled-by'
	AND V2 != 'used-to-produce'
	AND V2 != 'controls-production-of'
	AND V2 != 'controls-phosphorylation-of'")

# Include all interactions
#filteredInteractions <- interactions

# GENERATE IGRAPH GRAPHS ----
# Used directed, otherwise finding information on interactions becomes harder 
gWithType <- graph.edgelist(as.matrix(filteredInteractions[, c(1, 3)]), directed=TRUE)
E(gWithType)$type <- filteredInteractions[, 2]

# Remove duplicates for these generic interactions (interactions without type)
filteredInteractionsGeneric <- filteredInteractions[,c(1,3)]
filteredInteractionsGeneric <- filteredInteractionsGeneric[which(!duplicated(filteredInteractionsGeneric)),]
g <- graph.edgelist(as.matrix(filteredInteractionsGeneric), directed=FALSE)

V(g)[which(V(g)$name %in% seedSet)]$type <- "seed"
V(g)[which(!(V(g)$name %in% seedSet))]$type <- "linker"

# INITIALIZE PVALS ----
pvals <- data.frame(name=V(g)$name, 
										type=V(g)$type, 
										pval=rep(NA, vcount(g)), 
										pvalAdj=rep(NA, vcount(g)), 
										globalDegree=degree(g),
										neighborsInSeedSet=rep(NA, vcount(g)),
										seedSetLength=rep(NA, vcount(g)),
										genesNotInSeedSet=rep(NA, vcount(g)))

seedSetInNetwork <- V(g)[which(V(g)$type == "seed")]$name

# PERFORM LINKER SEARCH ----
pb <- txtProgressBar(min=1, max=vcount(g), style=3)

for (i in 1:vcount(g)) {
	setTxtProgressBar(pb, i)
	
	if(V(g)[i]$type == "linker") {
		neighborSet <- neighbors(g, V(g)[i], mode="all")
		
		# Neighbors returns multiple entries if there are multiple edges to neighbors
		neighborSet <- unique(neighborSet)

		#  Count number of neighbors for this linker in the seedSetInNetwork
		counter <- length(which(V(g)[neighborSet]$name %in% seedSetInNetwork))
		
		if(counter > 1) {
			globalDegree <- degree(g)[V(g)[i]]
			
			neighborsInSeedSet <- length(which(V(g)[neighborSet]$name %in% seedSetInNetwork)) 
			seedSetLength <- sum(V(g)$type == "seed")
			genesNotInSeedSet <- length(which(!(V(g)$name %in% seedSetInNetwork)))
			#neigborsNotInSeedSet <- length(which(!(V(g)[neighborSet]$name %in% seedSetInNetwork)))
			
			pval <- phyper(neighborsInSeedSet-1, 
										 seedSetLength, 
										 genesNotInSeedSet, 
										 globalDegree,
										 lower.tail=FALSE)
			
			V(g)[i]$pval <- pval
			
			pvals[i,"pval"] <- pval
			pvals[i,"neighborsInSeedSet"] <- neighborsInSeedSet
			pvals[i,"seedSetLength"] <- seedSetLength
			pvals[i,"genesNotInSeedSet"] <- genesNotInSeedSet
		} else {
			# Linkers without p-values are not connected multiple seeds
			pvals[i,"pval"] <- NA			
		}
	}
}	

# ADJUST PVALS ----
idx <- which(!(is.na(pvals[,"pval"])), arr.ind=TRUE)
pvals[idx,"pvalAdj"] <- p.adjust(pvals[idx,"pval"], "fdr")

# REMOVE INSIGNIFICANT LINKERS ----
# NA is NULL in sqldf
removeSet <- sqldf("SELECT *
	FROM pvals
	WHERE type = 'linker' 
	AND ((pvalAdj > 0.05) OR (pvalAdj IS NULL))")

removeSetVector <- c(userRemoveSet, as.vector(removeSet$name))
g2 <- delete.vertices(g, V(g)[removeSetVector])

# EXPORT GENERIC INTERACTIONS WITHOUT TYPE ----
write.graph(g2, file="~/Desktop/linkerNetwork.txt", format="ncol")

# EXPORT GENERIC INTERACTIONS WITH TYPE ----
# Iterate through all the edges in the generic g2 network and grab the edges
# with corresponding vertices from gWithType. Add these edges and their types
g3 <- graph.empty(0, directed=TRUE) + vertices(V(g2)$name)

for(i in 1:ecount(g2)) {
	v1 <- get.edgelist(g2)[i,1]
	v2 <- get.edgelist(g2)[i,2]
		
	tmpType <- E(gWithType)[v1 %--% v2]$type
	tmpEdges <- get.edges(gWithType, E(gWithType)[v1 %--% v2])
	
	for(j in 1:nrow(tmpEdges)) {
		g3 <- g3 + edge(V(gWithType)[tmpEdges[j,1]]$name, 
										V(gWithType)[tmpEdges[j,2]]$name, 
										type=tmpType[j])
	}
}

outputSif <- data.frame(participant1=get.edgelist(g3)[,1], 
												interactionType=E(g3)$type, 
												participant2=get.edgelist(g3)[,2])
write.table(outputSif, file="~/Desktop/linkerNetworkwithTypes.txt", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

# Display added in linkers 
V(g2)$name[which(!(V(g2)$name %in% seedSet))]

# PERFORM ENRICHMENT ANALYSIS ----
geneSets <- list(geneSet1=V(g2)$name)
gseaResults <- getPathwayOverrepresentation(geneSets, pathwayList="all", totalGenesNum=nrow(elNetMolDataNCI60$exp))

sortedGseaResults <- gseaResults[with(gseaResults, order(pval_adjust)), ]

write.table(sortedGseaResults, file="~/Desktop/gseaResults.txt", quote=FALSE, sep="\t", row.names=FALSE)

sort(seedSet[seedSet %in% V(g2)$name])

# FIND MODULES ----
ebc <- edge.betweenness.community(g2, directed=F)

topMod <- names(sort(table(ebc$membership), decreasing=TRUE)[1])
topModIdx <- which(ebc$membership == topMod)
topModGenes <- V(g2)$name[topModIdx]

# GET DRUG GO RESULTS ----
goResults <- getHgncDrugGoResults(V(g2)$name)
tmp <- goResults[which(!(is.na(goResults[,3]))), ]
sortedGoResults <- tmp[with(tmp, order(Gene)), ]
goResultsTable <- table(sortedGoResults[,3])

# VISUALIZE NETWORK ----
tryCatch({
	tmp <- elNetResults$predictorWts
	names(tmp) <- substr(names(elNetResults$predictorWts), 4, 100)

	gCw <- new("graphNEL", edgemode = "directed")

	gCw <- initNodeAttribute(gCw, "weight", "numeric", 0)
	gCw <- initNodeAttribute(gCw, "type", "char", "NA")
	gCw <- initEdgeAttribute(gCw, "edgeType", "char", "undefined")
	
	for(i in 1:vcount(g2)) {
		gCw <- graph::addNode(V(g2)[i]$name, gCw)	
		
		if(V(g2)[i]$name %in% names(tmp)) {
			nodeData(gCw, V(g2)[i]$name, "weight") <- as.numeric(tmp[V(g2)[i]$name])
		}
		
		nodeData(gCw, V(g2)[i]$name, "type") <- V(g2)[i]$type
	}
	
	for(i in 1:ecount(g2)) {
		v1 <- get.edgelist(g2)[i,1]
		v2 <- get.edgelist(g2)[i,2]
		
		tmpType <- E(gWithType)[v1 %--% v2]$type
		tmpEdges <- get.edges(gWithType, E(gWithType)[v1 %--% v2])
		
		for(j in 1:nrow(tmpEdges)) {
			gCw <- graph::addEdge(V(gWithType)[tmpEdges[j,1]]$name, 
														V(gWithType)[tmpEdges[j,2]]$name, 
														gCw)
			
			edgeData(gCw, 
							 V(gWithType)[tmpEdges[j,1]]$name, 
							 V(gWithType)[tmpEdges[j,2]]$name, 
							 "edgeType") <- tmpType[j]
		}
	}
	
	cw <- new.CytoscapeWindow('linkerNetwork', graph=gCw, overwriteWindow=TRUE)
	
	displayGraph(cw)
	layoutNetwork(cw, layout.name="force-directed")
	redraw(cw)
}, error=function(e) {
	cat("ERROR: ", e$message, "\n")
})
