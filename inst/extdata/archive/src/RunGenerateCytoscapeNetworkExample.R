# NOTE: RCytoscape properties: http://rcytoscape.systemsbiology.net/versions/current/index.html
library(RCytoscape)

source(file.path(.lmp, "gene_set_pathway_analysis", "src", "readGmtPathwayList.R"))

# LOAD GENE LIST -----------------------------------------
interactions <- read.table(file.path(.lmp, "gene_set_pathway_analysis", "data", "pc_sif_050814.gz"), sep="\t", as.is=TRUE, quote="")

pathways <- readGmtPathwayList(file.path(.lmp, "gene_set_pathway_analysis", "data", "c2.cp.v4.0.symbols.gmt"))
genes <- pathways[["REACTOME_SIGNALLING_TO_ERKS"]]

tmp <- interactions[which(interactions$V1 %in% genes), ]
filteredNetwork <- tmp[which(tmp$V3 %in% genes), ]

# FILTER INTERACTIONS ---------------------------------------
# The interaction types to remove 
notTheseIntTypes <- c("GENERIC_OF", "IN_SAME_COMPONENT")
filteredNetwork <- filteredNetwork[!(filteredNetwork$V2 %in% notTheseIntTypes),]

# Remove self interactions 
filteredNetwork <- filteredNetwork[filteredNetwork$V1 != filteredNetwork$V3,]

# APPEND NODES/EDGES TO NETWORK ----------------------------------------------
g <- new("graphNEL", edgemode="directed")

genesToBeAdded <- unique(c(as.vector(filteredNetwork$V1), as.vector(filteredNetwork$V3)))

for(gene in genesToBeAdded) {
	g <- graph::addNode(gene, g)
}

for(i in 1:nrow(filteredNetwork)) {
	g <- graph::addEdge(filteredNetwork[i,1], filteredNetwork[i,3], g)
	cat(filteredNetwork[i,1], " ", filteredNetwork[i,2], " ", filteredNetwork[i,3], "\n") 
} 

# INITIALIZE CYTOSCAPE GRAPH -----------------------------------------------------
cw <- new.CytoscapeWindow("example", graph=g, overwriteWindow=TRUE)
setDefaultBackgroundColor (cw, "#0000FF")
setDefaultNodeColor(cw, "#FFFFFF")
setDefaultNodeBorderColor(cw, "#000000")
setDefaultNodeBorderWidth(cw, 5)
setDefaultEdgeColor (cw, "#00FF00")

# ADD NODE/EDGE ATTRIBUTES ------------------------------------------------------
cg <- cw@graph 
cg <- initNodeAttribute(graph=cg, attribute.name="moleculeType", attribute.type="char", default.value="undefined")
cg <- initEdgeAttribute(graph=cg, attribute.name="edgeType", attribute.type="char", default.value="unspecified")

for(i in 1:length(genesToBeAdded)) {	
	if(genesToBeAdded[i] %in% genes) {
		nodeData(cg, genesToBeAdded[i], "moleculeType") <- "inPathway"
	} else {
		nodeData(cg, genesToBeAdded[i], "moleculeType") <- "outPathway"	
	}
}

for(i in 1:nrow(filteredNetwork)) {
	edgeData(cg, filteredNetwork[i,1], filteredNetwork[i,3], 'edgeType') <- filteredNetwork[i,2]
} 

# DISPLAY NETWORK IN CYTOSCAPE ---------------------------------------------------
cw <- setGraph(cw, cg)

# MODIFY PROPERTIES BASED ON ATTRIBUTES ----------------------------------------
displayGraph(cw)

edgeTypeValues <- c("REACTS_WITH", "INTERACTS_WITH", "STATE_CHANGE", "CO_CONTROL", "SEQUENTIAL_CATALYSIS")
edgeColors <- c("#FF0000", "#FFFF00", "#00FFFF", "#FF00FF", "#0000FF")

nodeTypeValues <- c("inPathway", "outPathway")
nodeColors <- c("#FF0000", "#00FF00")

setNodeColorRule(cw, "moleculeType", nodeTypeValues, nodeColors, "lookup")
setEdgeColorRule(cw, "edgeType", edgeTypeValues, edgeColors, mode="lookup")

# LAYOUT NETWORK ---------------------------------------------------------------
layoutNetwork(cw, layout.name="force-directed")
redraw(cw)

fitContent(cw)

# SAVE NETWORK TO IMAGE ----------------------------------------------------------------
imageType <- "pdf" 
saveImage (cw, "example.pdf", "pdf", 2.0) 
