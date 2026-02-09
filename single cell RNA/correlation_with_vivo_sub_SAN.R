# Load your Seurat object
library(Seurat)
library(UCell)
library(pheatmap)
library(Biobase)
library(patchwork)
library(ggplot2)

SAN_subtype_genes <- read.csv(file = "fetal_SAN_subtype_enriched_genes.csv",row.names = 1)

SAN_subtype_genes <- SAN_subtype_genes[-which(grepl("^MT-",rownames(SAN_subtype_genes))),]

head <- SAN_subtype_genes[SAN_subtype_genes$cluster == 'SAN-head',]
tail <- SAN_subtype_genes[SAN_subtype_genes$cluster == 'SAN-tail',]
TZ <- SAN_subtype_genes[SAN_subtype_genes$cluster == 'SAN-TZ',]

head$gene <- ifelse(
  head$avg_log2FC > 0,
  paste0(head$gene, "+"),
  ifelse(
    head$avg_log2FC < 0,
    paste0(head$gene, "-"),
    head$gene
  )
)

tail$gene <- ifelse(
  tail$avg_log2FC > 0,
  paste0(tail$gene, "+"),
  ifelse(
    tail$avg_log2FC < 0,
    paste0(tail$gene, "-"),
    tail$gene
  )
)

TZ$gene <- ifelse(
  TZ$avg_log2FC > 0,
  paste0(TZ$gene, "+"),
  ifelse(
    TZ$avg_log2FC < 0,
    paste0(TZ$gene, "-"),
    TZ$gene
  )
)

gene.signatures <- list(
  head = head$gene,
  tail = tail$gene,
  TZ = TZ$gene
)

# Compute signature scores per cell
# With UCell (fast & robust):
my_obj <- AddModuleScore_UCell(my_obj, features = gene.signatures, name = names(gene.signatures))

# This creates meta.data columns like "ACMACM", "VCMVCM", "SANSAN"
head(my_obj@meta.data)

# Visualize similarity
# UMAP plots
FeaturePlot(my_obj, features = "headhead", reduction = "umap") + labs(x = 'UMAP1',y = 'UMAP2',title = "SAN-head Module")

FeaturePlot(my_obj, features = "tailtail", reduction = "umap") + labs(x = 'UMAP1',y = 'UMAP2',title = "SAN-tail Module")

FeaturePlot(my_obj, features = "TZTZ", reduction = "umap") + labs(x = 'UMAP1',y = 'UMAP2',title = "SAN-TZ Module")

# Violin plots by cluster
VlnPlot(my_obj, features = c("headhead"), group.by = "celltype",pt.size = 0) + labs(x = '',y = '',title = "SAN-head Module")

VlnPlot(my_obj, features = c("tailtail"), group.by = "celltype",pt.size = 0) + labs(x = '',y = '',title = "SAN-tail Module")

VlnPlot(my_obj, features = c("TZTZ"), group.by = "celltype",pt.size = 0) + labs(x = '',y = '',title = "SAN-TZ Module")

# Heatmap of cluster-average scores
# Calculate average per cluster
cluster_scores <- AggregateExpression(my_obj, features = unlist(gene.signatures),
                                      assays = "RNA", group.by = "celltype", return.seurat = FALSE)

# Instead of raw expression, use metadata scores
cluster_means <- aggregate(my_obj@meta.data[, c("headhead","tailtail","TZTZ")],
                           by = list(cluster = Idents(my_obj)), mean)

rownames(cluster_means) <- cluster_means$cluster
cluster_means <- cluster_means[,-1]

# Plot heatmap
pheatmap(cluster_means, cluster_rows = FALSE, cluster_cols = FALSE, scale = "column")



