# load the library
library(CellChat)
library(patchwork)
# options(stringsAsFactors = FALSE)
library(Seurat)

# load the data
fetal_heart <- readRDS("fetal_heart_normal_spatial_annotated(006).rds")

SpatialDimPlot(fetal_heart,label =  TRUE, label.size = 3)

# prepare input data for CellChat analysis
data_input <- GetAssayData(fetal_heart,layer = "data",assay = "SCT")

meta <- data.frame(labels = Idents(fetal_heart),row.names = names(Idents(fetal_heart)))

unique(meta$labels)

# load spatial imaging information
spatial_locs <- Seurat::GetTissueCoordinates(fetal_heart,scale = NULL, cols = c("imagerow","imagecol"))

scale_factors <- jsonlite::fromJSON(txt = file.path("Fetal_heart006/outs/spatial",'scalefactors_json.json'))

scale_factors <- list(spot.diameter = 65, 
                      spot = scale_factors$spot_diameter_fullres,
                      fiducial = scale_factors$fiducial_diameter_fullres,
                      hires = scale_factors$tissue_hires_scalef,
                      lowres = scale_factors$tissue_lowres_scalef)


# Create a CellChat object
cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels",
                           datatype = "spatial",coordinates = spatial_locs[,c(1,2)] ,scale.factors = scale_factors)

cellchat <- addMeta(cellchat, meta = meta)

cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group


# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)


# simply use the default CellChatDB
# CellChatDB_use <- subsetDB(CellChatDB,search = "Secreted Signaling")

CellChatDB_use <- CellChatDB

cellchat@DB <- CellChatDB_use


# Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)

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

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.8, margin = 0.2)

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 0.8, margin = 0.2)

df.net <- subsetCommunication(cellchat, sources.use = 5, targets.use = 2)

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

# We can also show all the significant interactions (L-R pairs) from some cell groups to other cell groups using netVisual_bubble.
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,3,4,5,6,7), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(1,3,4,5,6,7), targets.use = 2 , remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 3, targets.use = 1 , remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 2, targets.use = 1 , remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1,2) , remove.isolate = FALSE)

# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 10,height = 20,font.size = 4)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 10,height = 20,font.size = 4)
ht1 + ht2

# Plot the signaling gene expression distribution using violin/dot plot
pathways.show <- "COLLAGEN"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "LAMININ"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

# MDK-NCL
pathways.show <- "MK"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "EPHB"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "EPHA"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "AGRN"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "HSPG"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "RELN"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "PTN"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "CDH"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "NCAM"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

pathways.show <- "VWF"
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)

plotGeneExpression(cellchat, signaling = "NRG")
plotGeneExpression(cellchat, signaling = "EPHA")
plotGeneExpression(cellchat, signaling = "AGRN",enriched.only = FALSE)


# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")

netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Blues")

netVisual_bubble(cellchat,remove.isolate = FALSE)

# Bubble plot
netVisual_bubble(cellchat,sources.use = 5,targets.use = 2,remove.isolate = FALSE)
# Chord diagram
netVisual_chord_gene(cellchat, sources.use = 5, targets.use = 2, lab.cex = 0.5,legend.pos.y = 30)

# Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat,signaling = "SEMA3")

# Visualize each signaling pathway using different plots
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)

spatialFeaturePlot(cellchat, features = c("COL1A1","ITGB1"))

#you can also intuitive way to visualize the main senders (sources) and receivers (targets) in two-dimensional space using scatter plots.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MK", "PTN"))

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("NPNT", "RELN", "CADM", "GRN","DESMOSOME"))

library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
# Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 3.
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")

# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")


nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")

# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")

saveRDS(cellchat, file = "fetal_heart_spatial(006)_cellchat.rds")


# subset cellchat to see only interactions between neuron and SAN
# Subset the communication matrix
df.comm <- subsetCommunication(cellchat)
# Rank by communication strength (you can change the threshold as needed)
df.comm_Neuron_SAN <- df.comm[df.comm$source == 'Neuron' & df.comm$target == 'SAN',]
df.comm_Neuron_SAN <- df.comm_Neuron_SAN[df.comm_Neuron_SAN$pathway_name != 'COLLAGEN' & df.comm_Neuron_SAN$pathway_name != 'LAMININ' & df.comm_Neuron_SAN$pathway_name != 'SEMA3',]
df.comm_Neuron_SAN <- df.comm_Neuron_SAN[order(-df.comm_Neuron_SAN$prob), ]  # Sort by strength (highest first)
# Get the top N pairs (e.g., top 20)
top_interactions <- head(df.comm_Neuron_SAN, 5)

netVisual_bubble(cellchat,
                 sources.use = 3,  # From cell types (sources)
                 targets.use = 1,    # To cell types (targets)
                 signaling = top_interactions$pathway_name,
                 remove.isolate = TRUE,
                 grid.on = FALSE)  + coord_flip()

ggplot(df.comm_Neuron_SAN[1:10,], aes(x=prob,y=reorder(interaction_name,prob),color=prob,size=-log10(pval+1))) + 
  geom_point()+
  scale_color_gradient(low = "red", high = "blue") + 
  theme_bw() + 
  xlab('Communication Prob')+
  ylab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


df.comm_SAN_CM <- df.comm[df.comm$source == 'SAN' & df.comm$target == 'CM',]
df.comm_SAN_CM <- df.comm_SAN_CM[df.comm_SAN_CM$pathway_name != 'COLLAGEN' & df.comm_SAN_CM$pathway_name != 'LAMININ' & df.comm_SAN_CM$pathway_name != 'SEMA3',]
df.comm_SAN_CM <- df.comm_SAN_CM[order(-df.comm_SAN_CM$prob), ]  # Sort by strength (highest first)
# Get the top N pairs (e.g., top 20)
top_interactions <- head(df.comm_SAN_CM, 5)

netVisual_bubble(cellchat,
                 sources.use = 1,  # From cell types (sources)
                 targets.use = 2,    # To cell types (targets)
                 signaling = top_interactions$pathway_name,
                 remove.isolate = TRUE,
                 grid.on = FALSE)  + coord_flip()

