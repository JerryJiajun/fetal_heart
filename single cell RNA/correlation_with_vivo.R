# Load your Seurat object
library(Seurat)
library(UCell)
library(pheatmap)
library(Biobase)
library(patchwork)
library(ggplot2)

sn_markers <- read.csv(file = "markers_of_each_cluster_RNA.csv",row.names = 1)

sn_markers <- sn_markers[-which(grepl("^MT-",rownames(sn_markers))),]

ACM <- sn_markers[sn_markers$cluster == 'ACM',]
Endothelial <- sn_markers[sn_markers$cluster == 'Endothelial',]
Epicardial <- sn_markers[sn_markers$cluster == 'Epicardial',]
Fibroblast <- sn_markers[sn_markers$cluster == 'Fibroblast',]
Neuronal <- sn_markers[sn_markers$cluster == 'Neuronal_1',]
SAN <- sn_markers[sn_markers$cluster == 'SAN',]
VCM <- sn_markers[sn_markers$cluster == 'VCM',]

gene.signatures <- list(
  ACM = setdiff(ACM$gene[1:100],union(SAN$gene[1:100],VCM$gene[1:100])),
  Endothelial = Endothelial$gene[1:40],
  Epicardial = Epicardial$gene[1:40],
  Fibroblast = Fibroblast$gene[1:40],
  Neuronal = setdiff(Neuronal$gene[1:40],union(ACM$gene[1:40],Epicardial$gene[1:40])),
  SAN = setdiff(SAN$gene[1:100],union(ACM$gene[1:100],VCM$gene[1:100])),
  VCM = setdiff(VCM$gene[1:100],union(SAN$gene[1:100],ACM$gene[1:100]))
)

# Compute signature scores per cell
# With UCell (fast & robust):
my_obj <- AddModuleScore_UCell(my_obj, features = gene.signatures, name = names(gene.signatures))

# This creates meta.data columns like "ACMACM", "VCMVCM", "SANSAN"
head(my_obj@meta.data)

# Visualize similarity
# UMAP plots
FeaturePlot(my_obj, features = "SANSAN", reduction = "umap") + labs(x = 'UMAP1',y = 'UMAP2',title = "SAN Module")

FeaturePlot(my_obj, features = "ACMACM", reduction = "umap") + labs(x = 'UMAP1',y = 'UMAP2',title = "ACM Module")

FeaturePlot(my_obj, features = "VCMVCM", reduction = "umap") + labs(x = 'UMAP1',y = 'UMAP2',title = "VCM Module")

FeaturePlot(my_obj, features = "NeuronalNeuronal", reduction = "umap") + labs(x = 'UMAP1',y = 'UMAP2',title = "Neuronal Module")

# Violin plots by cluster
VlnPlot(my_obj, features = c("SANSAN"), group.by = "celltype",pt.size = 0) + labs(x = '',y = '',title = "SAN Module")

VlnPlot(my_obj, features = c("ACMACM"), group.by = "celltype",pt.size = 0) + labs(x = '',y = '',title = "ACM Module")

VlnPlot(my_obj, features = c("VCMVCM"), group.by = "celltype",pt.size = 0) + labs(x = '',y = '',title = "VCM Module")

VlnPlot(my_obj, features = c("NeuronalNeuronal"), group.by = "celltype",pt.size = 0) + labs(x = '',y = '',title = "Neuronal Module")

# Heatmap of cluster-average scores
# Calculate average per cluster
cluster_scores <- AggregateExpression(my_obj, features = unlist(gene.signatures),
                                      assays = "RNA", group.by = "celltype", return.seurat = FALSE)

# Instead of raw expression, use metadata scores
cluster_means <- aggregate(my_obj@meta.data[, c("ACMACM","EndothelialEndothelial","EpicardialEpicardial","FibroblastFibroblast","NeuronalNeuronal", "SANSAN", "VCMVCM")],
                           by = list(cluster = Idents(my_obj)), mean)

rownames(cluster_means) <- cluster_means$cluster
cluster_means <- cluster_means[,-1]

# Plot heatmap
pheatmap(cluster_means, cluster_rows = FALSE, cluster_cols = FALSE, scale = "column")

custom_labels <- c("ACM", "Endothelial", "Epicardial", "Neuronal","SAN","VCM")
pheatmap(cluster_means[c(1,2,3,6,8,9),c(1,2,3,5,6,7)], cluster_rows = FALSE, cluster_cols = FALSE, scale = "column",labels_col = custom_labels)


