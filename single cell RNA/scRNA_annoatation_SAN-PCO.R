library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)

SAN_PCO.data <- Read10X(data.dir = "SAN-PCO/filtered_feature_bc_matrix/")

SAN_PCO <- CreateSeuratObject(counts = SAN_PCO.data, project = "SAN_PCO", min.cells = 3, min.features = 200)

SAN_PCO


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
SAN_PCO[["percent.mt"]] <- PercentageFeatureSet(SAN_PCO, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(SAN_PCO@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(SAN_PCO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(SAN_PCO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SAN_PCO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

SAN_PCO <- subset(SAN_PCO, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)

SAN_PCO <- NormalizeData(SAN_PCO, normalization.method = "LogNormalize", scale.factor = 10000)

SAN_PCO <- FindVariableFeatures(SAN_PCO, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SAN_PCO), 10)

plot3 <- VariableFeaturePlot(SAN_PCO)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

# Scaling the data
#'regress out' heterogeneity associated with mitochondrial contamination
SAN_PCO <- ScaleData(SAN_PCO, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
SAN_PCO <- RunPCA(SAN_PCO, features = VariableFeatures(object = SAN_PCO))

# Examine and visualize PCA results a few different ways
print(SAN_PCO[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(SAN_PCO, dims = 1:2, reduction = "pca")

DimPlot(SAN_PCO, reduction = "pca")

# DimHeatmap(SAN_PCO, dims = 1, cells = 500, balanced = TRUE)

# DimHeatmap(SAN_PCO, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset
# SAN_PCO <- JackStraw(SAN_PCO, num.replicate = 100)
# SAN_PCO <- ScoreJackStraw(SAN_PCO, dims = 1:20)

# JackStrawPlot(SAN_PCO, dims = 1:20,ymax = 0.8)

ElbowPlot(SAN_PCO)

#Cluster the cells
SAN_PCO <- FindNeighbors(SAN_PCO, dims = 1:20)
SAN_PCO <- FindClusters(SAN_PCO, resolution = 0.5)
# resolution can adjust from 0.2 to 1.5

head(Idents(SAN_PCO), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
SAN_PCO <- RunUMAP(SAN_PCO, dims = 1:20)

DimPlot(SAN_PCO, reduction = "umap")

SAN_PCO <- RunTSNE(SAN_PCO, dims = 1:20)

DimPlot(SAN_PCO, reduction = "tsne")

UMAPPlot(SAN_PCO,label=TRUE,repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

TSNEPlot(SAN_PCO,label=TRUE) + labs(x = 'TSNE1',y = 'TSNE2')

# find all markers of cluster 6,5
cluster0_vs_all <- FindMarkers(SAN_PCO, ident.1 = 0, min.pct = 0.1)
head(cluster0_vs_all, n = 20)

cluster0_vs_7 <- FindMarkers(SAN_PCO, ident.1 = 0,ident.2 = 7, min.pct = 0.1)
head(cluster0_vs_7, n = 20)


# find markers for every cluster compared to all remaining cells, report only the positive ones
SAN_PCO.markers <- FindAllMarkers(SAN_PCO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SAN_PCO.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- SAN_PCO.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(SAN_PCO, features = top20$gene) + NoLegend()
write.csv(SAN_PCO.markers, file = "ACM_VCM_SAN_cluster_markers.csv")

#find each cluster characters
VlnPlot(SAN_PCO, features = c("SHOX2"),pt.size = 0)
FeaturePlot(SAN_PCO,features = c("SHOX2"))
FeaturePlot(SAN_PCO,features = c("SHOX2"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("MYH6"))
FeaturePlot(SAN_PCO,features = c("MYH6"))
FeaturePlot(SAN_PCO,features = c("MYH6"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("HCN1"))
FeaturePlot(SAN_PCO,features = c("HCN1"))
FeaturePlot(SAN_PCO,features = c("HCN1"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("HCN4"))
FeaturePlot(SAN_PCO,features = c("HCN4"))
FeaturePlot(SAN_PCO,features = c("HCN4"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("CACNA1D"))
FeaturePlot(SAN_PCO,features = c("CACNA1D"))
FeaturePlot(SAN_PCO,features = c("CACNA1D"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("CACNA1G"))
FeaturePlot(SAN_PCO,features = c("CACNA1G"))
FeaturePlot(SAN_PCO,features = c("CACNA1G"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("GJC1"))
FeaturePlot(SAN_PCO,features = c("GJC1"))
FeaturePlot(SAN_PCO,features = c("GJC1"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("ISL1"))
FeaturePlot(SAN_PCO,features = c("ISL1"))
FeaturePlot(SAN_PCO,features = c("ISL1"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("TBX3"))
FeaturePlot(SAN_PCO,features = c("TBX3"))
FeaturePlot(SAN_PCO,features = c("TBX3"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("TBX5"))
FeaturePlot(SAN_PCO,features = c("TBX5"))
FeaturePlot(SAN_PCO,features = c("TBX5"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("TBX18"))
FeaturePlot(SAN_PCO,features = c("TBX18"))
FeaturePlot(SAN_PCO,features = c("TBX18"),reduction = "tsne")

VlnPlot(SAN_PCO, features = c("NKX2-5"))
FeaturePlot(SAN_PCO,features = c("NKX2-5"))
FeaturePlot(SAN_PCO,features = c("NKX2-5"),reduction = "tsne")

table(Idents(SAN_PCO))

# Atrial cluster 0,1: NPPA,MYL7,NR2F2,KCNA5,GJA5
VlnPlot(SAN_PCO, features = c("NPPA"),pt.size = 0)
FeaturePlot(SAN_PCO,features = c("NPPA"))
DotPlot(SAN_PCO,features = c("NPPA","NR2F2","GJA5")) + RotatedAxis()

# Fibroblast cluster 2,5: POSTN,DCN,ACTA2,COL1A1,PDGFRA
VlnPlot(SAN_PCO, features = c("COL3A1"),pt.size = 0)
FeaturePlot(SAN_PCO,features = c("DCN"))
DotPlot(SAN_PCO,features = c("DCN","POSTN","COL1A1","PDGFRA","COL3A1")) + RotatedAxis()

# cluster 3 -- SAN pacemaker cells
VlnPlot(SAN_PCO, features = c("SHOX2"),pt.size = 0)
FeaturePlot(SAN_PCO,features = c("SHOX2"))
DotPlot(SAN_PCO,features = c("SHOX2","HCN1")) + RotatedAxis()
DotPlot(SAN_PCO,features = c("SHOX2","HCN4","HCN1","TBX5","CACNA1G")) + RotatedAxis()

# cluster 4,10,13-- Epithelial cells (EPCAM,CDH1)
VlnPlot(SAN_PCO, features = c("EPCAM"),pt.size = 0)
FeaturePlot(SAN_PCO,features = c("EPCAM"))
FeaturePlot(SAN_PCO,features = c("EPCAM"),reduction = "tsne")
DotPlot(SAN_PCO,features = c("EPCAM","CDH1")) + RotatedAxis()

# cluster 6 -- Endothelial cells (PECAM1,CDH5,VWF)
VlnPlot(SAN_PCO, features = c("PECAM1"),pt.size = 0)
DotPlot(SAN_PCO,features = c("PECAM1","CDH5")) + RotatedAxis()
FeaturePlot(SAN_PCO,features = c("CDH5"))
FeaturePlot(SAN_PCO,features = c("CDH5"),reduction = "tsne")
DotPlot(SAN_PCO,features = c("PECAM1","CDH5")) + RotatedAxis()

# cluster 8-- Proliferating cells (TOP2A,CDK1,CENPF)
# TOP2A Highly expressed in proliferating cardiomyocytes and cardiac progenitor cells
# CDK1 Found in cycling cardiomyocytes and cardiac stem/progenitor cells
# CENPF Expressed in mitotic cells, including dividing cardiomyocytes and progenitor cells
VlnPlot(SAN_PCO, features = c("TOP2A"),pt.size = 0)
FeaturePlot(SAN_PCO,features = c("TOP2A"))
DotPlot(SAN_PCO,features = c("TOP2A","CDK1","CENPF")) + RotatedAxis()

# cluster 14--Cardiac Neuron (STMN2,SST,TAC3,NTRK3,ELAVL4)
# STMN2 Highly expressed in neuronal-like cells of the cardiac conduction system
# SST Found in intrinsic cardiac neurons and conduction system cells, regulating parasympathetic signaling and cardiac autonomic function
# TAC3 Expressed in neuronal and conduction cells, potentially involved in neurotransmitter signaling within the heart
# NTRK3 Found in neuronal-like cells and conduction system cells, playing a role in neurotrophic signaling and cardiac development
# ELAVL4 Expressed in neuronal-like conduction system cells, regulating mRNA stability and synaptic plasticity
DotPlot(SAN_PCO,features = c("STMN2","SST","TAC3","NTRK3","ELAVL4")) + RotatedAxis()

DotPlot(SAN_PCO,features = c("TH","PRPH","CHAT","PHOX2B","RET")) + RotatedAxis()

VlnPlot(SAN_PCO, features = c("STMN2"),pt.size = 0) 
FeaturePlot(SAN_PCO,features = c("STMN2")) + labs(x = 'UMAP1',y = 'UMAP2')

# cluster 11--neural crest
VlnPlot(SAN_PCO, features = c("SOX2"),pt.size = 0) 
FeaturePlot(SAN_PCO,features = c("SOX2")) + labs(x = 'UMAP1',y = 'UMAP2')

# cluster 12 -- epicardial (WT1,TBX18) WT1, TBX18, SEMA3D, ALDH1A2, GATA5, TCF21
VlnPlot(SAN_PCO, features = c("WT1"),pt.size = 0)
FeaturePlot(SAN_PCO,features = c("WT1"))
FeaturePlot(SAN_PCO,features = c("WT1"),reduction = "tsne")

DotPlot(SAN_PCO,features = c("WT1","BNC1")) + RotatedAxis()

saveRDS(SAN_PCO, file = "SAN-PCO.rds")

# add annotations
SAN_PCO <- RenameIdents(SAN_PCO, '3' = 'SAN')

SAN_PCO <- RenameIdents(SAN_PCO, '0' = 'ACM','1' = 'ACM')

SAN_PCO <- RenameIdents(SAN_PCO, '4' = 'Epithelial','7' = 'Epithelial','10' = 'Epithelial','13' = 'Epithelial')

SAN_PCO <- RenameIdents(SAN_PCO, '2' = 'FB','5' = 'FB','9' = 'FB')

SAN_PCO <- RenameIdents(SAN_PCO, '6' = 'Endothelial')

SAN_PCO <- RenameIdents(SAN_PCO, '12' = 'Epicardial')

SAN_PCO <- RenameIdents(SAN_PCO, '14' = 'Neuronal')

SAN_PCO <- RenameIdents(SAN_PCO, '11' = 'Neural_Crest')

SAN_PCO <- RenameIdents(SAN_PCO, '8' = 'Proliferating')

SAN_PCO$celltype <- Idents(SAN_PCO)

levels(SAN_PCO) <- c("ACM","Endothelial","Epicardial","Epithelial","FB",
                     "Neuronal","Neural_Crest","Proliferating","SAN")

UMAPPlot(SAN_PCO,label=TRUE, repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

DotPlot(SAN_PCO,features = c("NPPA","CDH5","WT1","EPCAM","POSTN","STMN2","SOX2","TOP2A","SHOX2")) + RotatedAxis()

DotPlot(SAN_PCO,features = c("NPPA","CDH5","WT1","EPCAM","COL1A1","STMN2","SOX2","TOP2A","SHOX2")) + RotatedAxis() + labs(x = '',y = '')

VlnPlot(SAN_PCO,features = "NPPA",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("NPPA")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "CDH5",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("CDH5")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "WT1",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("WT1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "EPCAM",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("EPCAM")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "COL1A1",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("COL1A1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "STMN2",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("STMN2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "SOX2",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("SOX2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "TOP2A",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("TOP2A")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "SHOX2",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("SHOX2")) + labs(x = 'UMAP1',y = 'UMAP2')

VlnPlot(SAN_PCO,features = "PECAM1",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("PECAM1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "PHOX2B",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("PHOX2B")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "BNC1",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("BNC1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "ISL1",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("ISL1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "POSTN",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("POSTN")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "DCN",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("DCN")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(SAN_PCO,features = "NR2F2",pt.size = 0)
FeaturePlot(SAN_PCO,features = c("NR2F2")) + labs(x = 'UMAP1',y = 'UMAP2')

DotPlot(SAN_PCO,features = c("SHOX2","HCN1","HCN4","CACNA1D","CACNA1G")) + RotatedAxis() + labs(x='',y = '')
VlnPlot(SAN_PCO,features = "SHOX2",pt.size = 0)
FeaturePlot(SAN_PCO,features = "SHOX2")
VlnPlot(SAN_PCO,features = "HCN1",pt.size = 0)
FeaturePlot(SAN_PCO,features = "HCN1")
VlnPlot(SAN_PCO,features = "CACNA1F",pt.size = 0)
FeaturePlot(SAN_PCO,features = "CACNA1D")
VlnPlot(SAN_PCO,features = "TBX5",pt.size = 0)
FeaturePlot(SAN_PCO,features = "TBX5")

DotPlot(SAN_PCO,features = c("TH","PRPH","CHAT","PHOX2B","RET","PHOX2A")) + RotatedAxis()

DotPlot(subset(SAN_PCO,subset = celltype %in% c('ACM','SAN')),features = c("NPPA","GJA5","OPCML","CSMD1","BMPER","ACOXL","SHOX2","ISL1","HCN1","HCN4","TBX5","TBX18")) + RotatedAxis() + labs(x='',y='')

saveRDS(SAN_PCO, file = "SAN-PCO_annotated.rds")
