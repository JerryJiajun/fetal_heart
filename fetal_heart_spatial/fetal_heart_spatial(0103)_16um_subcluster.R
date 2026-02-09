library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

fetal_heart <- readRDS("~/spatial_fetal_heart_new/fetal_heart_0103_16um_annotated.rds")

DimPlot(fetal_heart,reduction = "umap.016um",label = TRUE, repel = TRUE) + labs(x = "UMAP1",y = "UMAP2")
DotPlot(fetal_heart,features = c("NPPA","CDH5","DCN","SPP1","PLP1","HBG1","SHOX2","MYH11")) + RotatedAxis() + labs(x = '',y='')
SpatialDimPlot(fetal_heart,label = TRUE,label.size = 3,repel = TRUE)
levels(fetal_heart) <- c("ACM","EC","FB","Macrophages","Neuronal",'RBC',"SAN","SMC")

mycols <- c(
  "orange",
  "darkgreen",
  "#F781BF", 
  "#FF7F00", 
  "#984EA3",
  "#AEC6CF",
  "#E41A1C",
  "#A65628" 
)
SpatialDimPlot(fetal_heart,label = TRUE,label.size = 3,repel = TRUE) + scale_fill_manual(values = mycols) + theme_void()

sub_SAN <- subset(fetal_heart,subset = celltype == 'SAN')

sub_SAN <- NormalizeData(sub_SAN)
sub_SAN <- FindVariableFeatures(sub_SAN)
sub_SAN <- ScaleData(sub_SAN)
sub_SAN <- RunPCA(sub_SAN, reduction.name = "pca.016um")
sub_SAN <- FindNeighbors(sub_SAN, reduction = "pca.016um", dims = 1:20)
sub_SAN <- FindClusters(sub_SAN, resolution = 1, cluster.name = "seurat_cluster.016um")
sub_SAN <- RunUMAP(sub_SAN, reduction = "pca.016um", reduction.name = "umap.016um", dims = 1:20)

DimPlot(sub_SAN,reduction = "umap.016um",label = TRUE, repel = TRUE) + labs(x = "UMAP1",y = "UMAP2")
DotPlot(sub_SAN,features = c("TBX18","TBX3","SHOX2","HCN4","NPPA","NKX2-5")) + RotatedAxis() + labs(x = '',y='')
SpatialDimPlot(sub_SAN,label = TRUE,label.size = 3,repel = TRUE)

saveRDS(sub_SAN,file = 'fetal_heart_0103_16um_sub_SAN.rds')

sub_SAN <- RenameIdents(sub_SAN, '3' = 'SAN-head','2' = 'SAN-tail','1' = 'SAN-TZ','0' = 'SAN-TZ')
sub_SAN$subtype <- Idents(fetal_heart)
DimPlot(sub_SAN,reduction = "umap.016um",label = TRUE, repel = TRUE) + labs(x = "UMAP1",y = "UMAP2")
DotPlot(sub_SAN,features = c("SHOX2","TBX18","NKX2-5","NPPA")) + RotatedAxis() + labs(x = '',y='')
SpatialDimPlot(sub_SAN,label = TRUE,label.size = 3,repel = TRUE)

DotPlot(sub_SAN,features = c('MST1','AMOT','SAV1','LATS1','MOB1A')) + RotatedAxis() + labs(x = '',y='')

SAN_subtype_genes <- FindAllMarkers(sub_SAN,logfc.threshold = 0.1,min.pct = 0.1)

write.csv(SAN_subtype_genes,file = "fetal_SAN_subtype_enriched_genes.csv")
