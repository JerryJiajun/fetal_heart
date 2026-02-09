library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)

# load the spatial data
data.dir <- "Fetal_heart006/outs/"
fetal_heart_normal <- Load10X_Spatial(
  data.dir = data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "tissue_hires_image.png",
)

# data preprocessing
plot1 <- VlnPlot(fetal_heart_normal,features = "nCount_Spatial",pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(fetal_heart_normal,features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1,plot2)

fetal_heart_normal <- SCTransform(fetal_heart_normal,assay = "Spatial",verbose = FALSE)

# gene expression visualization
SpatialFeaturePlot(fetal_heart_normal,features = c("SHOX2","TNNT2","MYH11"))

SpatialFeaturePlot(fetal_heart_normal,features = c("SHOX2","TNNT2","MYH11"),alpha = c(0.1,1))

# pt.size.factor- This will scale the size of the spots. Default is 1.6
p1 <- SpatialFeaturePlot(fetal_heart_normal,features = "SHOX2",pt.size.factor = 1)
# alpha - minimum and maximum transparency. Default is c(1, 1)
# Try setting to alpha c(0.1, 1), to downweight the transparency of points with lower expression
p2 <- SpatialFeaturePlot(fetal_heart_normal,features = "SHOX2", alpha = c(0.1,1))

p1 + p2


# Dimensionality reduction, clustering and visualization
fetal_heart_normal <- RunPCA(fetal_heart_normal,assay = "SCT",verbose = FALSE)

ElbowPlot(fetal_heart_normal)

fetal_heart_normal <- FindNeighbors(fetal_heart_normal,reduction = "pca",dims = 1:20)

fetal_heart_normal <- FindClusters(fetal_heart_normal,verbose = FALSE,resolution = 0.5)

fetal_heart_normal <- RunUMAP(fetal_heart_normal,reduction = "pca",dim = 1:20)

p1 <- DimPlot(fetal_heart_normal,reduction = "umap",label = TRUE)

p2 <- SpatialDimPlot(fetal_heart_normal,label = TRUE,label.size = 4)

p1 + p2

SpatialDimPlot(fetal_heart_normal,cells.highlight = CellsByIdentities(object = fetal_heart_normal,idents = c(0,2,3,5,7,8)),facet.highlight = TRUE,ncol = 3)

SpatialFeaturePlot(fetal_heart_normal,features = "SHOX2",interactive = TRUE)

LinkedDimPlot(fetal_heart_normal)

# Identification of Spatially Variable Features
SpatialFeaturePlot(fetal_heart_normal,features = c("SHOX2","TBX3","HCN4"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("SHOX2","TBX3","HCN4"),ncol = 3)
DotPlot(subset(fetal_heart_normal,idents = c(0,2,3,5,7,8)),features = c("SHOX2","TBX3","TBX18")) + RotatedAxis()
SpatialDimPlot(fetal_heart_normal,cells.highlight = CellsByIdentities(object = fetal_heart_normal,idents = 8),facet.highlight = TRUE,ncol = 1)

SpatialFeaturePlot(fetal_heart_normal,features = c("TNNT2","MYL7"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("TNNT2","MYL7"))
DotPlot(fetal_heart_normal,features = c("TNNT2","MYL7")) + RotatedAxis()
SpatialDimPlot(fetal_heart_normal,cells.highlight = CellsByIdentities(object = fetal_heart_normal,idents = c(0,2,3,5,7)),facet.highlight = TRUE,ncol = 5)

SpatialFeaturePlot(fetal_heart_normal,features = c("NRXN1","PLP1","NRXN3"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("NRXN1","PLP1","NRXN3"),ncol = 3)
DotPlot(fetal_heart_normal,features = c("NRXN1","PLP1","NRXN3")) + RotatedAxis()
SpatialDimPlot(fetal_heart_normal,cells.highlight = CellsByIdentities(object = fetal_heart_normal,idents = 11),facet.highlight = TRUE,ncol = 1)

SpatialFeaturePlot(fetal_heart_normal,features = c("CD3E","IL7R","CD40LG"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CD3E","IL7R","CD40LG"),ncol = 3)
DotPlot(fetal_heart_normal,features = c("CD3E","IL7R","CD40LG")) + RotatedAxis()
SpatialDimPlot(fetal_heart_normal,cells.highlight = CellsByIdentities(object = fetal_heart_normal,idents = 13),facet.highlight = TRUE,ncol = 1)

SpatialFeaturePlot(fetal_heart_normal,features = c("CD163","CD14","C1QA"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CD163","CD14","C1QA"),ncol = 3)
DotPlot(fetal_heart_normal,features = c("CD163","CD14","C1QA")) + RotatedAxis()
SpatialDimPlot(fetal_heart_normal,cells.highlight = CellsByIdentities(object = fetal_heart_normal,idents = c(1,4)),facet.highlight = TRUE,ncol = 2)

SpatialFeaturePlot(fetal_heart_normal,features = c("ACTA2","MYH11"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("ACTA2","MYH11"))
DotPlot(fetal_heart_normal,features = c("ACTA2","MYH11")) + RotatedAxis()
SpatialDimPlot(fetal_heart_normal,cells.highlight = CellsByIdentities(object = fetal_heart_normal,idents = c(6,9,10,12)),facet.highlight = TRUE,ncol = 4)

# Erythroid cells (Red Blood Cell)
SpatialFeaturePlot(fetal_heart_normal,features = c("HBB","HBG1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("HBB","HBG1"))
DotPlot(fetal_heart_normal,features = c("HBB","HBG1")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("WT1","BNC1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("WT1","BNC1"))
DotPlot(fetal_heart_normal,features = c("WT1","BNC1")) + RotatedAxis()
SpatialDimPlot(fetal_heart_normal,cells.highlight = CellsByIdentities(object = fetal_heart_normal,idents = c(4,7)),facet.highlight = TRUE,ncol = 2)

saveRDS(fetal_heart_normal,file = "fetal_heart_normal_spatial(006).rds")

# Identification of Spatially Variable Features
# find markers for every cluster compared to all remaining cells, report only the positive ones
fetal_heart_normal <- PrepSCTFindMarkers(fetal_heart_normal)
SAN.markers <- FindAllMarkers(fetal_heart_normal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SAN.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top20 <- SAN.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(SAN, features = top20$gene) + NoLegend()
# add annotations
fetal_heart_normal <- RenameIdents(fetal_heart_normal, '5' = 'Lymphoid')
fetal_heart_normal <- RenameIdents(fetal_heart_normal, '3' = 'VSM','6' = 'VSM','8' = 'VSM')
fetal_heart_normal <- RenameIdents(fetal_heart_normal, '1' = 'Neuron')
fetal_heart_normal <- RenameIdents(fetal_heart_normal, '0' = 'CM','2' = 'CM','4' = 'CM','9' = 'CM')
fetal_heart_normal <- RenameIdents(fetal_heart_normal, '7' = 'SAN')


fetal_heart_normal$celltype <- Idents(fetal_heart_normal)

DimPlot(fetal_heart_normal,reduction = "umap",label = TRUE)

SpatialDimPlot(fetal_heart_normal,label = TRUE,label.size = 5)

saveRDS(fetal_heart_normal,file = "fetal_heart_normal_spatial_annotated(006).rds")

SpatialFeaturePlot(fetal_heart_normal,features = c("TNNT2","MYL7"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("TNNT2","MYL7"),ncol = 2)
DotPlot(fetal_heart_normal,features = c("TNNT2","MYL7")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("SHOX2","TBX3","TBX18"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("SHOX2","TBX3","TBX18"))
DotPlot(fetal_heart_normal,features = c("SHOX2","TBX3","TBX18")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("NRXN1","PLP1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("NRXN1","PLP1"))
DotPlot(fetal_heart_normal,features = c("NRXN1","PLP1")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("SCARA5","ANGPTL1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("SCARA5","ANGPTL1"),ncol = 2)
DotPlot(fetal_heart_normal,features = c("SCARA5","ANGPTL1")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("CD3E","IL7R"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CD3E","IL7R"))
DotPlot(fetal_heart_normal,features = c("CD3E","IL7R")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("CD163","CD14"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CD163","CD14"))
DotPlot(fetal_heart_normal,features = c("CD163","CD14")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("CCL14","SPP1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CCL14","SPP1"))
DotPlot(fetal_heart_normal,features = c("CCL14","SPP1")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("ACTA2","MYH11"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("ACTA2","MYH11"))
DotPlot(fetal_heart_normal,features = c("ACTA2","MYH11")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("WT1","BNC1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("WT1","BNC1"))
DotPlot(fetal_heart_normal,features = c("WT1","BNC1")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("HBB","HBG1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("HBB","HBG1"))
DotPlot(fetal_heart_normal,features = c("HBB","HBG1")) + RotatedAxis()

SpatialFeaturePlot(fetal_heart_normal,features = c("PLP1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("PLP1"))
VlnPlot(fetal_heart_normal,features = c("PLP1"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("WT1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("WT1"))
VlnPlot(fetal_heart_normal,features = c("WT1"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("CCL14"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CCL14"))
VlnPlot(fetal_heart_normal,features = c("CCL14"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("CD3E"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CD3E"))
VlnPlot(fetal_heart_normal,features = c("CD3E"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("MYH11"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("MYH11"))
VlnPlot(fetal_heart_normal,features = c("MYH11"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("SHOX2"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("SHOX2"))
VlnPlot(fetal_heart_normal,features = c("SHOX2"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("ANGPT1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("ANGPT1"))
VlnPlot(fetal_heart_normal,features = c("ANGPT1"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("GJA1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("GJA1"))
VlnPlot(fetal_heart_normal,features = c("GJA1"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("GJC1"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("GJC1"))
VlnPlot(fetal_heart_normal,features = c("GJC1"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("TBX3"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("TBX3"))
VlnPlot(fetal_heart_normal,features = c("TBX3"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("CACNA1D"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CACNA1D"))
VlnPlot(fetal_heart_normal,features = c("CACNA1D"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("CCL21"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CCL21"))
VlnPlot(fetal_heart_normal,features = c("CCL21"),pt.size = 0)

levels(fetal_heart_normal) <- c("VSM","SAN","Neuron","Lymphoid","CM")
DotPlot(fetal_heart_normal,features = c("MYH11","SHOX2","PLP1",'CCL21',"ANGPT1")) + RotatedAxis()

DotPlot(fetal_heart_normal,features = c("SHOX2","TBX3","CACNA1D")) + RotatedAxis()

# General Autonomic Neuron Markers--PHOX2A,PHOX2B,TH,DBH
DotPlot(fetal_heart_normal,features = c("TH","DBH","PHOX2A","PHOX2B")) + RotatedAxis()
SpatialFeaturePlot(fetal_heart_normal,features = c("STMN2"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("STMN2"))
VlnPlot(fetal_heart_normal,features = c("STMN2"),pt.size = 0)

SpatialFeaturePlot(fetal_heart_normal,features = c("PRPH"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("PRPH"))
VlnPlot(fetal_heart_normal,features = c("PRPH"),pt.size = 0)

# Sypathetic Neurons--TH, DBH, PRPH, NTPK1
DotPlot(fetal_heart_normal,features = c("TH","DBH","PRPH")) + RotatedAxis()
SpatialFeaturePlot(fetal_heart_normal,features = c("TH"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("TH"))
VlnPlot(fetal_heart_normal,features = c("TH"),pt.size = 0)
SpatialFeaturePlot(fetal_heart_normal,features = c("DBH"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("DBH"))

# Parasympathetic Neurons--CHAT, PHOX2B, RET
DotPlot(fetal_heart_normal,features = c("CHAT","PHOX2B","RET")) + RotatedAxis()
SpatialFeaturePlot(fetal_heart_normal,features = c("CHAT"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("CHAT"))
SpatialFeaturePlot(fetal_heart_normal,features = c("PHOX2B"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("PHOX2B"))
SpatialFeaturePlot(fetal_heart_normal,features = c("RET"),alpha = c(0.1,1))
FeaturePlot(fetal_heart_normal,features = c("RET"))

DotPlot(fetal_heart_normal,features = c("TH","DBH","PRPH","CHAT","PHOX2B","RET")) + RotatedAxis()
