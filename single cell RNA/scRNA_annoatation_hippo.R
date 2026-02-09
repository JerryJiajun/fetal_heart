library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

Control.data <- Read10X(data.dir = "SAN_Control/filtered_feature_bc_matrix/")

Control <- CreateSeuratObject(counts = Control.data, project = "Control", min.cells = 3, min.features = 200)

Control

Forskolin.data <- Read10X(data.dir = "SAN_Forskolin/filtered_feature_bc_matrix/")

Forskolin <- CreateSeuratObject(counts = Forskolin.data, project = "Forskolin", min.cells = 3, min.features = 200)

Forskolin

XMU.data <- Read10X(data.dir = "SAN_XMU/filtered_feature_bc_matrix/")

XMU <- CreateSeuratObject(counts = XMU.data, project = "XMU", min.cells = 3, min.features = 200)

XMU

Combine <- merge(x = Control, y = c(Forskolin,XMU) ,add.cell.ids = c("Control","Forskolin","XMU"),project = "combine")


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Combine[["percent.mt"]] <- PercentageFeatureSet(Combine, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(Combine@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(Combine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Combine, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Combine, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Combine <- subset(Combine, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 20)

Combine <- NormalizeData(Combine, normalization.method = "LogNormalize", scale.factor = 10000)

Combine <- FindVariableFeatures(Combine, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Combine), 10)

plot3 <- VariableFeaturePlot(Combine)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

#Scaling the data
#'regress out' heterogeneity associated with mitochondrial contamination
Combine <- ScaleData(Combine, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
Combine <- RunPCA(Combine, features = VariableFeatures(object = Combine))

# Examine and visualize PCA results a few different ways
print(Combine[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Combine, dims = 1:2, reduction = "pca")

DimPlot(Combine, reduction = "pca")

# DimHeatmap(Combine, dims = 1, cells = 500, balanced = TRUE)

# DimHeatmap(Combine, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset
# Combine <- JackStraw(Combine, num.replicate = 100)
# Combine <- ScoreJackStraw(Combine, dims = 1:20)

# JackStrawPlot(Combine, dims = 1:20,ymax = 0.8)

ElbowPlot(Combine)

# join layers
Combine <- JoinLayers(Combine)

#Cluster the cells
Combine <- FindNeighbors(Combine, dims = 1:20)
Combine <- FindClusters(Combine, resolution = 0.5)
# resolution can adjust from 0.2 to 1.5

head(Idents(Combine), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
Combine <- RunUMAP(Combine, dims = 1:20)

DimPlot(Combine, reduction = "umap")

Combine <- RunTSNE(Combine, dims = 1:20)

DimPlot(Combine, reduction = "tsne")

UMAPPlot(Combine,label=TRUE,repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
UMAPPlot(Combine,label=TRUE,split.by = 'orig.ident',repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
UMAPPlot(Combine,group.by = 'orig.ident') + labs(x = 'UMAP1',y = 'UMAP2')

TSNEPlot(Combine,label=TRUE) + labs(x = 'TSNE1',y = 'TSNE2')

# find all markers of cluster 6,5
cluster6_vs_all <- FindMarkers(Combine, ident.1 = 6, min.pct = 0.1)
head(cluster6_vs_all, n = 20)

cluster0_vs_7 <- FindMarkers(Combine, ident.1 = 0,ident.2 = 7, min.pct = 0.1)
head(cluster0_vs_7, n = 20)


# find markers for every cluster compared to all remaining cells, report only the positive ones
Combine.markers <- FindAllMarkers(Combine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combine.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- Combine.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(Combine, features = top15$gene) + NoLegend()

#find each cluster characters
VlnPlot(Combine, features = c("SHOX2"),pt.size = 0)
FeaturePlot(Combine,features = c("SHOX2"))
FeaturePlot(Combine,features = c("SHOX2"),reduction = "tsne")

VlnPlot(Combine, features = c("MYH6"))
FeaturePlot(Combine,features = c("MYH6"))
FeaturePlot(Combine,features = c("MYH6"),reduction = "tsne")

VlnPlot(Combine, features = c("HCN1"))
FeaturePlot(Combine,features = c("HCN1"))
FeaturePlot(Combine,features = c("HCN1"),reduction = "tsne")

VlnPlot(Combine, features = c("HCN4"))
FeaturePlot(Combine,features = c("HCN4"))
FeaturePlot(Combine,features = c("HCN4"),reduction = "tsne")

VlnPlot(Combine, features = c("CACNA1D"))
FeaturePlot(Combine,features = c("CACNA1D"))
FeaturePlot(Combine,features = c("CACNA1D"),reduction = "tsne")

VlnPlot(Combine, features = c("CACNA1G"))
FeaturePlot(Combine,features = c("CACNA1G"))
FeaturePlot(Combine,features = c("CACNA1G"),reduction = "tsne")

VlnPlot(Combine, features = c("GJC1"))
FeaturePlot(Combine,features = c("GJC1"))
FeaturePlot(Combine,features = c("GJC1"),reduction = "tsne")

VlnPlot(Combine, features = c("ISL1"))
FeaturePlot(Combine,features = c("ISL1"))
FeaturePlot(Combine,features = c("ISL1"),reduction = "tsne")

VlnPlot(Combine, features = c("TBX3"))
FeaturePlot(Combine,features = c("TBX3"))
FeaturePlot(Combine,features = c("TBX3"),reduction = "tsne")

VlnPlot(Combine, features = c("TBX5"))
FeaturePlot(Combine,features = c("TBX5"))
FeaturePlot(Combine,features = c("TBX5"),reduction = "tsne")

VlnPlot(Combine, features = c("TBX18"))
FeaturePlot(Combine,features = c("TBX18"))
FeaturePlot(Combine,features = c("TBX18"),reduction = "tsne")

VlnPlot(Combine, features = c("NKX2-5"))
FeaturePlot(Combine,features = c("NKX2-5"))
FeaturePlot(Combine,features = c("NKX2-5"),reduction = "tsne")

table(Idents(Combine))

# Fibroblast cluster 2 : POSTN,DCN,ACTA2,COL1A1,PDGFRA
VlnPlot(Combine, features = c("DCN"),pt.size = 0)
FeaturePlot(Combine,features = c("DCN"))+ labs(x = 'UMAP1',y = 'UMAP2')
DotPlot(Combine,features = c("DCN","POSTN","COL1A1","PDGFRA","COL3A1")) + RotatedAxis()

# cluster 0,3,5,6 -- pacemaker cells
VlnPlot(Combine, features = c("SHOX2"),pt.size = 0)
FeaturePlot(Combine,features = c("SHOX2"))+ labs(x = 'UMAP1',y = 'UMAP2')
DotPlot(Combine,features = c("SHOX2","HCN1")) + RotatedAxis()
DotPlot(Combine,features = c("SHOX2","HCN4","HCN1","TBX5","CACNA1G")) + RotatedAxis()

# cluster 1,4,7-- Epithelial cells (EPCAM,CDH1)
VlnPlot(Combine, features = c("EPCAM"),pt.size = 0)
FeaturePlot(Combine,features = c("EPCAM"))+ labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(Combine,features = c("EPCAM"),reduction = "tsne")
DotPlot(Combine,features = c("EPCAM","CDH1")) + RotatedAxis()

# cluster  -- Endothelial cells (PECAM1,CDH5,VWF)
VlnPlot(Combine, features = c("PECAM1"),pt.size = 0)
DotPlot(Combine,features = c("PECAM1","CDH5")) + RotatedAxis()
FeaturePlot(Combine,features = c("CDH5"))+ labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(Combine,features = c("CDH5"),reduction = "tsne")
DotPlot(Combine,features = c("PECAM1","CDH5")) + RotatedAxis()


# cluster 9-- Proliferating cells (TOP2A,CDK1,CENPF)
# TOP2A Highly expressed in proliferating cardiomyocytes and cardiac progenitor cells
# CDK1 Found in cycling cardiomyocytes and cardiac stem/progenitor cells
# CENPF Expressed in mitotic cells, including dividing cardiomyocytes and progenitor cells
VlnPlot(Combine, features = c("TOP2A"),pt.size = 0)
FeaturePlot(Combine,features = c("TOP2A"))+ labs(x = 'UMAP1',y = 'UMAP2')
DotPlot(Combine,features = c("TOP2A","CDK1","CENPF")) + RotatedAxis()

# cluster 2--Cardiac Neuron (STMN2,SST,TAC3,NTRK3,ELAVL4)
# STMN2 Highly expressed in neuronal-like cells of the cardiac conduction system
# SST Found in intrinsic cardiac neurons and conduction system cells, regulating parasympathetic signaling and cardiac autonomic function
# TAC3 Expressed in neuronal and conduction cells, potentially involved in neurotransmitter signaling within the heart
# NTRK3 Found in neuronal-like cells and conduction system cells, playing a role in neurotrophic signaling and cardiac development
# ELAVL4 Expressed in neuronal-like conduction system cells, regulating mRNA stability and synaptic plasticity
DotPlot(Combine,features = c("STMN2","PRPH")) + RotatedAxis()

DotPlot(Combine,features = c("TH","PRPH","CHAT","PHOX2B","RET")) + RotatedAxis()

VlnPlot(Combine, features = c("STMN2"),pt.size = 0) 
FeaturePlot(Combine,features = c("STMN2")) + labs(x = 'UMAP1',y = 'UMAP2')

# cluster 10 -- epicardial (WT1,TBX18) WT1, TBX18, SEMA3D, ALDH1A2, GATA5, TCF21
VlnPlot(Combine, features = c("WT1"),pt.size = 0)
FeaturePlot(Combine,features = c("WT1"))
FeaturePlot(Combine,features = c("WT1"),reduction = "tsne")
DotPlot(Combine,features = c("WT1","BNC1")) + RotatedAxis()

saveRDS(Combine, file = "CTL_forskolin_XMU.rds")


# add annotations
Combine <- RenameIdents(Combine, '1' = 'Epithelial','4' = 'Epithelial','7' = 'Epithelial',
                        '9' = 'Proliferating',
                        '8' = 'Epicardial','10' = 'Epicardial',
                        '2' = 'Neuronal',
                        '0' = 'SAN','3' = 'SAN','5' = 'SAN','6' = 'SAN','11' = 'SAN')

Combine$celltype <- Idents(Combine)

levels(Combine) <- c("ACM","Endothelial","Epicardial","Epicardial_progenitor","Epithelial","FB",
                     "Neuronal","Neural_Crest","Proliferating","SAN","SAN_progenitor","VCM")

saveRDS(Combine, file = "CTL_forskolin_XMU_annotated.rds")

UMAPPlot(Combine,label=TRUE, repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

UMAPPlot(Combine,label=TRUE,split.by = "orig.ident",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

UMAPPlot(Combine,group.by = 'orig.ident') + labs(x = 'UMAP1',y = 'UMAP2')

DotPlot(Combine,features = c("EPCAM","TOP2A","WT1","STMN2","SHOX2")) + RotatedAxis()

VlnPlot(Combine,features = "STMN2",pt.size = 0)
FeaturePlot(Combine,features = c("STMN2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "PHOX2B",pt.size = 0)
FeaturePlot(Combine,features = c("PHOX2B")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "CDH5",pt.size = 0)
FeaturePlot(Combine,features = c("CDH5")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "PECAM1",pt.size = 0)
FeaturePlot(Combine,features = c("PECAM1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "WT1",pt.size = 0)
FeaturePlot(Combine,features = c("WT1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "BNC1",pt.size = 0)
FeaturePlot(Combine,features = c("BNC1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "EPCAM",pt.size = 0)
FeaturePlot(Combine,features = c("EPCAM")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "SOX2",pt.size = 0)
FeaturePlot(Combine,features = c("SOX2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "SHOX2",pt.size = 0)
FeaturePlot(Combine,features = c("SHOX2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "ISL1",pt.size = 0)
FeaturePlot(Combine,features = c("ISL1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "COL3A1",pt.size = 0)
FeaturePlot(Combine,features = c("COL3A1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "POSTN",pt.size = 0)
FeaturePlot(Combine,features = c("POSTN")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "DCN",pt.size = 0)
FeaturePlot(Combine,features = c("DCN")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "TOP2A",pt.size = 0)
FeaturePlot(Combine,features = c("TOP2A")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "MYL2",pt.size = 0)
FeaturePlot(Combine,features = c("MYL2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "HEY2",pt.size = 0)
FeaturePlot(Combine,features = c("HEY2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "NPPA",pt.size = 0)
FeaturePlot(Combine,features = c("NPPA")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(Combine,features = "NR2F2",pt.size = 0)
FeaturePlot(Combine,features = c("NR2F2")) + labs(x = 'UMAP1',y = 'UMAP2')

DotPlot(Combine,features = c("SHOX2","HCN1","TBX3","TBX18")) + RotatedAxis()
VlnPlot(Combine,features = "SHOX2",pt.size = 0)
FeaturePlot(Combine,features = "SHOX2")
VlnPlot(Combine,features = "HCN1",pt.size = 0)
FeaturePlot(Combine,features = "HCN1")
VlnPlot(Combine,features = "CACNA1F",pt.size = 0)
FeaturePlot(Combine,features = "CACNA1D")
VlnPlot(Combine,features = "TBX5",pt.size = 0)
FeaturePlot(Combine,features = "TBX5")

DotPlot(Combine,features = c("MYH7","MYL2","MYL3","IRX4","HEY2","ANKRD1","FHL2","MYOZ2","CSRP3")) + RotatedAxis()
DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN_progenitor','SAN')),features = c("MYH7","MYL2","IRX4","HEY2","ANKRD1")) + RotatedAxis()

DotPlot(Combine,features = c("NPPA","GJA5","NR2F2","OPCML","CSMD1","BMPER","ACOXL")) + RotatedAxis()
DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN_progenitor','SAN')),features = c("GJA5","OPCML","CSMD1","BMPER","ACOXL")) + RotatedAxis()

DotPlot(Combine,features = c("SHOX2","ISL1","TBX18","HCN1")) + RotatedAxis()
DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN_progenitor','SAN')),features = c("SHOX2","ISL1","TBX18","HCN1")) + RotatedAxis()

# subset only SHOX2 positive SAN cells
# SAN <- subset(Combine, subset = seurat_clusters %in% c(0,3,5,6) & WT1 == 0 )
SAN <- subset(Combine, subset = celltype == 'SAN')

#Cluster the cells
SAN <- FindNeighbors(SAN, dims = 1:20)
SAN <- FindClusters(SAN, resolution = 0.8)
# resolution can adjust from 0.2 to 1.5
#Run non-linear dimensional reduction (UMAP/tSNE)
SAN <- RunUMAP(SAN, dims = 1:20)

UMAPPlot(SAN,label=TRUE, repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

UMAPPlot(SAN,label=TRUE,split.by = "orig.ident",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

VlnPlot(SAN,features = "SHOX2",pt.size = 0)
VlnPlot(SAN,features = "TBX18",pt.size = 0)
VlnPlot(SAN,features = "NPPA",pt.size = 0)
VlnPlot(SAN,features = "HCN4",pt.size = 0)
VlnPlot(SAN,features = "NKX2-5",pt.size = 0)

FeaturePlot(SAN,features = "SHOX2", repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "TBX18", repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "TBX3", repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "HCN4",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "NKX2-5",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "NPPA",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

SAN <- RenameIdents(SAN, '2'='SAN-head','8' = 'SAN-tail',
                    '5' = 'SAN-tail','6'='SAN-tail','7' = 'SAN-tail',
                    '0' = 'SAN-TZ','1' = 'SAN-TZ','3' = 'SAN-TZ','4' = 'SAN-TZ','9' = 'SAN-TZ')

SAN$celltype <- Idents(SAN)

saveRDS(SAN, file = "CTL_forskolin_XMU_SAN_subcluster.rds")

UMAPPlot(SAN,label=TRUE, repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

UMAPPlot(SAN,label=TRUE,split.by = "orig.ident",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

DotPlot(SAN,features = c("SHOX2","TBX18","NKX2-5","NPPA")) + RotatedAxis() + labs(x = '',y = '')

VlnPlot(SAN,features = "SHOX2",pt.size = 0)
VlnPlot(SAN,features = "TBX18",pt.size = 0)
VlnPlot(SAN,features = "NKX2-5",pt.size = 0)
VlnPlot(SAN,features = "NPPA",pt.size = 0)
VlnPlot(SAN,features = "HCN4",pt.size = 0)

FeaturePlot(SAN,features = "SHOX2", repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "TBX18", repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "TBX3", repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "HCN4",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "NKX2-5",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(SAN,features = "NPPA",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')
# find out the cellular compositions in each sample
meta <- SAN@meta.data
meta$origin <- meta$orig.ident
SAN@meta.data <- meta

celltype <- SAN$celltype
origin <- SAN$origin

celltype_table <- table(SAN$origin,SAN$celltype)
celltype_prop_table <- prop.table(celltype_table,margin = 1)

sample <- rep(rownames(celltype_prop_table),each = dim(celltype_prop_table)[2])
cell <- rep(colnames(celltype_prop_table),times = dim(celltype_prop_table)[1])
ratio <- c(t(celltype_prop_table))

combine_data <- data.frame(sample,cell,ratio)

combine_data$sample <- factor(combine_data$sample,levels = c('Control','Forskolin','XMU'))
# Stacked + percent
ggplot(combine_data, aes(fill=cell, y=ratio, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()

