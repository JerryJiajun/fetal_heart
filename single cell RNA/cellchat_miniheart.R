# load the library
library(CellChat)
library(patchwork)
library(Seurat)

# prepare input data for CellChat analysis
data_input <- GetAssayData(heart,layer = "data",assay = "RNA")


meta <- data.frame(labels = Idents(heart),row.names = names(Idents(heart)))

unique(meta$labels)

# Create a CellChat object
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# simply use the default CellChatDB
# CellChatDB_use <- subsetDB(CellChatDB,search = "Secreted Signaling")
# CellChatDB in human contains 1,939 validated molecular interactions
# including 61.8% of paracrine/autocrine signaling interactions
# 21.7% of extracellular matrix (ECM)-receptor interactions 
# and 16.5% of cell-cell contact interactions

CellChatDB_use <- CellChatDB

cellchat@DB <- CellChatDB_use

# Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


# Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.8, margin = 0.4)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 0.8, margin = 0.4)

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")

netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

netVisual_bubble(cellchat,remove.isolate = FALSE)

mat <- cellchat@net$weight
par(mfrow = c(3,4))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 10,height = 15,font.size = 4)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 10,height = 15,font.size = 4)
ht1 + ht2

# We can also show all the significant interactions (L-R pairs) from some cell groups to other cell groups using netVisual_bubble.
netVisual_bubble(cellchat, sources.use = -8, targets.use = 8, remove.isolate = FALSE,grid.on = FALSE) + coord_flip()
netVisual_bubble(cellchat, sources.use = -1, targets.use = 1, remove.isolate = FALSE,grid.on = FALSE) + coord_flip()

netVisual_bubble(cellchat, sources.use = -8, targets.use = 8, remove.isolate = FALSE,grid.on = TRUE)
netVisual_bubble(cellchat, sources.use = -1, targets.use = 1, remove.isolate = FALSE,grid.on = TRUE)

netVisual_bubble(cellchat, sources.use = 8, targets.use = 1, remove.isolate = FALSE,grid.on = TRUE)

netVisual_bubble(cellchat, sources.use = 6, targets.use = 8, remove.isolate = FALSE,grid.on = TRUE)

netVisual_bubble(cellchat, sources.use = 8, targets.use = -8, remove.isolate = FALSE,grid.on = TRUE)

# Plot the signaling gene expression distribution using violin/dot plot
pathways.show <- "NRG"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "EPHA"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "EPHB"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "NCAM"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "BMP"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "CDH"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

plotGeneExpression(cellchat, signaling = "NRG")
plotGeneExpression(cellchat, signaling = "EPHA")
plotGeneExpression(cellchat, signaling = "NRG",enriched.only = FALSE)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("NRG"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("NRG", "EPHA","EPHB","NCAM","CDH","BMP","FGF","CADM"),pattern = "incoming")

library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")


selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "incoming")

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")


saveRDS(cellchat, file = "cellchat_mini_heart.rds")


# subset cellchat to see only interactions between neuron and SAN
# Subset the communication matrix
df.comm <- subsetCommunication(cellchat)
# Rank by communication strength (you can change the threshold as needed)
df.comm_Neuron_SAN <- df.comm[df.comm$source == 'Neuronal' & df.comm$target == 'SAN',]
df.comm_Neuron_SAN <- df.comm_Neuron_SAN[order(-df.comm_Neuron_SAN$prob), ]  # Sort by strength (highest first)
# Get the top N pairs (e.g., top 20)
top_interactions <- head(df.comm_Neuron_SAN, 5)

netVisual_bubble(cellchat,
                 sources.use = 6,  # From cell types (sources)
                 targets.use = 8,    # To cell types (targets)
                 signaling = top_interactions$pathway_name,
                 remove.isolate = TRUE,
                 grid.on = FALSE)



# subset cellchat to see only interactions between SAN and ACM
# Subset the communication matrix
df.comm <- subsetCommunication(cellchat)
# Rank by communication strength (you can change the threshold as needed)
df.comm_SAN_ACM <- df.comm[df.comm$source == 'SAN' & df.comm$target == 'ACM',]
df.comm_SAN_ACM <- df.comm_SAN_ACM[order(-df.comm_SAN_ACM$prob), ]  # Sort by strength (highest first)
# Get the top N pairs (e.g., top 20)
top_interactions <- head(df.comm_SAN_ACM, 5)

netVisual_bubble(cellchat,
                 sources.use = 8,  # From cell types (sources)
                 targets.use = 1,    # To cell types (targets)
                 signaling = top_interactions$pathway_name,
                 remove.isolate = TRUE,
                 grid.on = FALSE)



# subset cellchat to see only interactions between ACM and VCM
# Subset the communication matrix
df.comm <- subsetCommunication(cellchat)
# Rank by communication strength (you can change the threshold as needed)
df.comm_ACM_VCM <- df.comm[df.comm$source == 'ACM' & df.comm$target == 'VCM',]
df.comm_ACM_VCM <- df.comm_ACM_VCM[order(-df.comm_ACM_VCM$prob), ]  # Sort by strength (highest first)
# Get the top N pairs (e.g., top 20)
top_interactions <- head(df.comm_ACM_VCM, 5)

netVisual_bubble(cellchat,
                 sources.use = 1,  # From cell types (sources)
                 targets.use = 9,    # To cell types (targets)
                 signaling = top_interactions$pathway_name,
                 remove.isolate = TRUE,
                 grid.on = FALSE)
