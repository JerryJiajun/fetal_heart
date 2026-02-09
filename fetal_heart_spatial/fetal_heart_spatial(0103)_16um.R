library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# load the spatial data
data.dir <- "2025-0103_heart/outs/binned_outputs/square_016um/"

fetal_heart <- Load10X_Spatial(
  data.dir = data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial.016um",
  slice = "tissue_hires_image.png"
)

VlnPlot(fetal_heart,features = "nCount_Spatial.016um",pt.size = 0.1) + NoLegend()
VlnPlot(fetal_heart,features = "nFeature_Spatial.016um",pt.size = 0.1) + NoLegend()

SpatialFeaturePlot(fetal_heart,features = "nCount_Spatial.016um") + theme(legend.position = "right")
SpatialFeaturePlot(fetal_heart,features = "nFeature_Spatial.016um") + theme(legend.position = "right")

fetal_heart <- subset(fetal_heart,subset = nCount_Spatial.016um > 100 )

fetal_heart <- NormalizeData(fetal_heart)

SpatialFeaturePlot(fetal_heart, features = "SHOX2",alpha = c(0.1,1),pt.size.factor = 5) + ggtitle("Fetal heart SHOX2 expression (8um)")

fetal_heart <- FindVariableFeatures(fetal_heart)
fetal_heart <- ScaleData(fetal_heart)
fetal_heart <- RunPCA(fetal_heart, reduction.name = "pca.016um")
fetal_heart <- FindNeighbors(fetal_heart, reduction = "pca.016um", dims = 1:20)
fetal_heart <- FindClusters(fetal_heart, resolution = 1, cluster.name = "seurat_cluster.016um")

fetal_heart <- RunUMAP(fetal_heart, reduction = "pca.016um", reduction.name = "umap.016um", dims = 1:20)

SpatialDimPlot(fetal_heart,cells.highlight = CellsByIdentities(object = fetal_heart,idents = 7),pt.size.factor = 3,facet.highlight = TRUE)

DimPlot(fetal_heart, reduction = "umap.016um", group.by = "seurat_cluster.016um", label = TRUE, repel = T) + labs(x = 'UMAP1',y = 'UMAP2',title = '') + NoLegend()

SpatialDimPlot(fetal_heart, group.by = "seurat_cluster.016um", label = FALSE, pt.size.factor = 1.5) + theme(legend.position = "right")

saveRDS(fetal_heart,file = 'fetal_heart_0103_16um.rds')

# Smooth Muscle Cells (SMCs), cluster 1,4,10,16,21,23
SpatialFeaturePlot(fetal_heart,features = c("MYH11"),alpha = c(0.1,1),pt.size.factor = 3)
FeaturePlot(fetal_heart,features = c("MYH11")) + labs(x = 'UMAP1', y = 'UMAP2')
DotPlot(fetal_heart,features = c("ACTA2","MYH11","TAGLN","CNN1")) + RotatedAxis()

# Endothelial Cells (ECs), cluster 7,9,20
DotPlot(fetal_heart,features = c("CDH5","PECAM1","VWF","KDR")) + RotatedAxis()

# Blood Endothelial Cells (bEC) cluster 9,20

# Lymphatic Endothelial Cells (lEC) cluster 7

# activated macrophages, cluster 6
DotPlot(fetal_heart,features = c("F13A1","SPP1","C1QA")) + RotatedAxis()

# Lymphatic endothelial cells, cluster 18
DotPlot(fetal_heart,features = c("CCL21","LYVE1")) + RotatedAxis()

# Neuronal Cells, cluster 12, 22
DotPlot(fetal_heart,features = c("PLP1","SOX10","PHOX2B")) + RotatedAxis()

# Fibroblasts, cluster 5,11,19
FeaturePlot(fetal_heart,features = c("DCN")) + labs(x = 'UMAP1', y = 'UMAP2')
DotPlot(fetal_heart,features = c("COL1A1","COL3A1","DCN","LUM")) + RotatedAxis()

# Red Blood Cells (RBCs), cluster 14,20
DotPlot(fetal_heart,features = c("HBG1","HBG2","HBA2")) + RotatedAxis()


# add annotations
fetal_heart <- RenameIdents(fetal_heart, '1' = 'SMC','4' = 'SMC','10' = 'SMC','16' = 'SMC','21' = 'SMC','23' = 'SMC')
fetal_heart <- RenameIdents(fetal_heart, '9' = 'EC','20' = 'EC')
fetal_heart <- RenameIdents(fetal_heart, '5' = 'FB','11' = 'FB','18' = 'FB','19' = 'FB')
fetal_heart <- RenameIdents(fetal_heart, '6' = 'Macrophages')
fetal_heart <- RenameIdents(fetal_heart, '12' = 'Neuronal','22' = 'Neuronal')
fetal_heart <- RenameIdents(fetal_heart, '14' = 'RBC','25' = 'RBC')
fetal_heart <- RenameIdents(fetal_heart, '0' = 'ACM','2' = 'ACM','3' = 'ACM','7' = 'ACM','8' = 'ACM','13' = 'ACM','17' = 'ACM','24' = 'ACM')
fetal_heart <- RenameIdents(fetal_heart, '15' = 'SAN')

fetal_heart$celltype <- Idents(fetal_heart)
levels(fetal_heart) <- c("ACM","EC","FB","Macrophages","Neuronal",'RBC',"SAN","SMC")
DimPlot(fetal_heart,reduction = "umap.016um",label = TRUE, repel = TRUE) + labs(x = "UMAP1",y = "UMAP2")
DotPlot(fetal_heart,features = c("MYL7","VWF","DCN","SPP1","PLP1","HBG1","SHOX2","MYH11")) + RotatedAxis() + labs(x = '',y='')
SpatialDimPlot(fetal_heart,label = TRUE,label.size = 3,repel = TRUE)
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

saveRDS(fetal_heart,file = 'fetal_heart_0103_16um_annotated.rds')


SpatialFeaturePlot(fetal_heart,features = c("NPPA"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("NPPA")) + labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("NPPA"),pt.size = 0) + labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("MYL7"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("MYL7")) + labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("MYL7"),pt.size = 0) + labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("CDH5"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("CDH5"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("CDH5"),pt.size = 0)+ labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("VWF"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("VWF"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("VWF"),pt.size = 0)+ labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("DCN"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("DCN"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("DCN"),pt.size = 0)+ labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("CCL21"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("CCL21"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("CCL21"),pt.size = 0)+ labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("SPP1"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("SPP1"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("SPP1"),pt.size = 0)+ labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("PLP1"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("PLP1"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("PLP1"),pt.size = 0)+ labs(x='')

DotPlot(fetal_heart,features = c("STMN2","PRPH","PHOX2B","RET","PHOX2A","CHAT","TH")) + RotatedAxis() + labs(x = '',y='')
SpatialFeaturePlot(fetal_heart,features = c("STMN2"),alpha = c(0.1,1),pt.size.factor = 5)
SpatialFeaturePlot(fetal_heart,features = c("PRPH"),alpha = c(0.1,1),pt.size.factor = 5)
SpatialFeaturePlot(fetal_heart,features = c("PHOX2B"),alpha = c(0.1,1),pt.size.factor = 5)
SpatialFeaturePlot(fetal_heart,features = c("RET"),alpha = c(0.1,1),pt.size.factor = 5)
SpatialFeaturePlot(fetal_heart,features = c("PHOX2A"),alpha = c(0.1,1),pt.size.factor = 5)
SpatialFeaturePlot(fetal_heart,features = c("CHAT"),alpha = c(0.1,1),pt.size.factor = 5)
SpatialFeaturePlot(fetal_heart,features = c("TH"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("STMN2"))+ labs(x='UMAP1',y='UMAP2')
FeaturePlot(fetal_heart,features = c("PRPH"))+ labs(x='UMAP1',y='UMAP2')
FeaturePlot(fetal_heart,features = c("PHOX2B"))+ labs(x='UMAP1',y='UMAP2')
FeaturePlot(fetal_heart,features = c("TH"))+ labs(x='UMAP1',y='UMAP2')
FeaturePlot(fetal_heart,features = c("DBH"))+ labs(x='UMAP1',y='UMAP2')
FeaturePlot(fetal_heart,features = c("RET"))+ labs(x='UMAP1',y='UMAP2')

SpatialFeaturePlot(fetal_heart,features = c("SHOX2"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("SHOX2"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("SHOX2"),pt.size = 0)+ labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("NRG1"),alpha = c(0.1,1),pt.size.factor = 5)
SpatialFeaturePlot(fetal_heart,features = c("ERBB2"),alpha = c(0.1,1),pt.size.factor = 5)
SpatialFeaturePlot(fetal_heart,features = c("ERBB4"),alpha = c(0.1,1),pt.size.factor = 5)
DotPlot(fetal_heart,features = c("NRG1","ERBB4")) + RotatedAxis() + coord_flip()
FeaturePlot(fetal_heart,features = c("NRG1"))+ labs(x='UMAP1',y='UMAP2')
FeaturePlot(fetal_heart,features = c("ERBB4"))+ labs(x='UMAP1',y='UMAP2')

SpatialFeaturePlot(fetal_heart,features = c("TEAD1"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("TEAD1"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("TEAD1"),pt.size = 0)+ labs(x='')
DotPlot(fetal_heart,features = c("ERBB4","TEAD1","YAP1")) + RotatedAxis() + coord_flip() + labs(x = '',y = '')

SpatialFeaturePlot(fetal_heart,features = c("YAP1"),alpha = c(0.1,1),pt.size.factor = 5)

SpatialFeaturePlot(fetal_heart,features = c("WWTR1"),alpha = c(0.1,1),pt.size.factor = 5)


DotPlot(fetal_heart,features = c("SHOX2","TBX3","CACNA1G","TBX18","HCN4","HCN1","ISL1")) + RotatedAxis() + labs(x = '',y = '')

DotPlot(fetal_heart,features = c("SHOX2","TBX3","CACNA1G","TBX18","HCN4","HCN1","ISL1","CACNA2D2")) + RotatedAxis() + labs(x = '',y = '')

SpatialFeaturePlot(fetal_heart,features = c("TBX3"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("TBX3"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("TBX3"),pt.size = 0)+ labs(x='')
SpatialFeaturePlot(fetal_heart,features = c("TBX18"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("TBX18"))+ labs(x='UMAP1',y='UMAP2')
SpatialFeaturePlot(fetal_heart,features = c("HCN4"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("HCN4"))+ labs(x='UMAP1',y='UMAP2')
SpatialFeaturePlot(fetal_heart,features = c("CACNA1G"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("CACNA1G"))+ labs(x='UMAP1',y='UMAP2')
SpatialFeaturePlot(fetal_heart,features = c("ISL1"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("ISL1"))+ labs(x='UMAP1',y='UMAP2')
SpatialFeaturePlot(fetal_heart,features = c("CACNA2D2"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("CACNA2D2"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("CACNA2D2"),pt.size = 0)+ labs(x='')


SpatialFeaturePlot(fetal_heart,features = c("MYH11"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("MYH11"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("MYH11"),pt.size = 0)+ labs(x='')

SAN_genes <- FindMarkers(fetal_heart,ident.1 = 'SAN',logfc.threshold = 0.25,min.pct = 0.1)
SAN_enriched_genes <- SAN_genes[SAN_genes$avg_log2FC > 0 & SAN_genes$p_val_adj < 1e-5,]
write.csv(SAN_enriched_genes,file = "SAN_enriched_genes.csv")

ACM_genes <- FindMarkers(fetal_heart,ident.1 = 'ACM',logfc.threshold = 0.25,min.pct = 0.1)
ACM_enriched_genes <- ACM_genes[ACM_genes$avg_log2FC > 0 & ACM_genes$p_val_adj < 1e-5,]
write.csv(ACM_enriched_genes,file = "ACM_enriched_genes.csv")

Neuronal_genes <- FindMarkers(fetal_heart,ident.1 = 'Neuronal',logfc.threshold = 0.25,min.pct = 0.1)
Neuronal_enriched_genes <- Neuronal_genes[Neuronal_genes$avg_log2FC > 0 & Neuronal_genes$p_val_adj < 1e-5,]
write.csv(Neuronal_enriched_genes,file = "Neuronal_enriched_genes.csv")

SpatialFeaturePlot(fetal_heart,features = c("RELN"),alpha = c(0.1,1),pt.size.factor = 5)
FeaturePlot(fetal_heart,features = c("RELN"))+ labs(x='UMAP1',y='UMAP2')
VlnPlot(fetal_heart,features = c("TUBB3"),pt.size = 0)+ labs(x='')

SpatialFeaturePlot(fetal_heart,features = c("TUBB3"),alpha = c(0.1,1),pt.size.factor = 5)

SpatialFeaturePlot(fetal_heart,features = c("UCHL1"),alpha = c(0.1,1),pt.size.factor = 5)

SpatialFeaturePlot(fetal_heart,features = c("NEFL"),alpha = c(0.1,1),pt.size.factor = 5)

SpatialFeaturePlot(fetal_heart,features = c("MPZ"),alpha = c(0.1,1),pt.size.factor = 5)

DotPlot(fetal_heart,features = c("TUBB3","UCHL1","NEFL","RELN","PRPH","PHOX2B","RET","PHOX2A","CHAT","TH","MPZ")) + RotatedAxis() + labs(x = '',y='')
