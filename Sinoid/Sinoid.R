library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
library(corrplot)


Sinoid.data <- Read10X(data.dir = "Epicardial_perturb/STARsolo/filtered/")

Sinoid <- CreateSeuratObject(counts = Sinoid.data, project = "Sinoid", min.cells = 3, min.features = 200)

Sinoid


# Lets examine a few genes in the first thirty cells
Sinoid.data[c("SHOX2", "TBX3", "TBX5"), 1:30]

dense.size <- object.size(as.matrix(Sinoid.data))

dense.size

sparse.size <- object.size(Sinoid.data)

sparse.size

dense.size/sparse.size


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Sinoid[["percent.mt"]] <- PercentageFeatureSet(Sinoid, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(Sinoid@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(Sinoid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Sinoid, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Sinoid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

Sinoid <- subset(Sinoid, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 15)

cell_barcodes <- rownames(Sinoid@meta.data)

names(cell_barcodes) <- colnames(Sinoid)

Sinoid <- AddMetaData(
  object = Sinoid,
  metadata = cell_barcodes,
  col.name = 'cell.barcodes'
)

Sinoid <- NormalizeData(Sinoid, normalization.method = "LogNormalize", scale.factor = 10000)

Sinoid <- FindVariableFeatures(Sinoid, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Sinoid), 10)

plot3 <- VariableFeaturePlot(Sinoid)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

#Scaling the data
all.genes <- rownames(Sinoid)
Sinoid <- ScaleData(Sinoid, features = all.genes)

#'regress out' heterogeneity associated with mitochondrial contamination
Sinoid <- ScaleData(Sinoid, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
Sinoid <- RunPCA(Sinoid, features = VariableFeatures(object = Sinoid))

# Examine and visualize PCA results a few different ways
print(Sinoid[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(Sinoid, dims = 1:2, reduction = "pca")

DimPlot(Sinoid, reduction = "pca")

DimHeatmap(Sinoid, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(Sinoid, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset
Sinoid <- JackStraw(Sinoid, num.replicate = 100)
Sinoid <- ScoreJackStraw(Sinoid, dims = 1:20)

JackStrawPlot(Sinoid, dims = 1:20,ymax = 0.8)

ElbowPlot(Sinoid)

#Cluster the cells
Sinoid <- FindNeighbors(Sinoid, dims = 1:20)
Sinoid <- FindClusters(Sinoid, resolution = 0.5)
# resolution can adjust from 0.2 to 1.5

head(Idents(Sinoid), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
Sinoid <- RunUMAP(Sinoid, dims = 1:20)

DimPlot(Sinoid, reduction = "umap")

Sinoid <- RunTSNE(Sinoid, dims = 1:20)

DimPlot(Sinoid, reduction = "tsne")

UMAPPlot(Sinoid,label=TRUE)

TSNEPlot(Sinoid,label=TRUE)

# find all markers of cluster 6,5
cluster6_vs_all <- FindMarkers(Sinoid, ident.1 = 6, min.pct = 0.1)
head(cluster6_vs_all, n = 20)

cluster5_vs_all <- FindMarkers(Sinoid, ident.1 = 5, min.pct = 0.1)
head(cluster5_vs_all, n = 20)

cluster2_vs_all <- FindMarkers(Sinoid, ident.1 = 2, min.pct = 0.1)
head(cluster2_vs_all, n = 20)

# find markers for every cluster compared to all remaining cells, report only the positive ones
Sinoid.markers <- FindAllMarkers(Sinoid, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Sinoid.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- Sinoid.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(Sinoid, features = top15$gene) + NoLegend()

#find each cluster characters
VlnPlot(Sinoid, features = c("SHOX2"))
FeaturePlot(Sinoid,features = c("SHOX2"))
FeaturePlot(Sinoid,features = c("SHOX2"),reduction = "tsne")

VlnPlot(Sinoid, features = c("MYH6"))
FeaturePlot(Sinoid,features = c("MYH6"))
FeaturePlot(Sinoid,features = c("MYH6"),reduction = "tsne")

VlnPlot(Sinoid, features = c("HCN1"))
FeaturePlot(Sinoid,features = c("HCN1"))
FeaturePlot(Sinoid,features = c("HCN1"),reduction = "tsne")

VlnPlot(Sinoid, features = c("HCN4"))
FeaturePlot(Sinoid,features = c("HCN4"))
FeaturePlot(Sinoid,features = c("HCN4"),reduction = "tsne")

VlnPlot(Sinoid, features = c("CACNA1D"))
FeaturePlot(Sinoid,features = c("CACNA1D"))
FeaturePlot(Sinoid,features = c("CACNA1D"),reduction = "tsne")

VlnPlot(Sinoid, features = c("CACNA1G"))
FeaturePlot(Sinoid,features = c("CACNA1G"))
FeaturePlot(Sinoid,features = c("CACNA1G"),reduction = "tsne")

VlnPlot(Sinoid, features = c("GJC1"))
FeaturePlot(Sinoid,features = c("GJC1"))
FeaturePlot(Sinoid,features = c("GJC1"),reduction = "tsne")

VlnPlot(Sinoid, features = c("ISL1"))
FeaturePlot(Sinoid,features = c("ISL1"))
FeaturePlot(Sinoid,features = c("ISL1"),reduction = "tsne")

VlnPlot(Sinoid, features = c("TBX3"))
FeaturePlot(Sinoid,features = c("TBX3"))
FeaturePlot(Sinoid,features = c("TBX3"),reduction = "tsne")

VlnPlot(Sinoid, features = c("TBX5"))
FeaturePlot(Sinoid,features = c("TBX5"))
FeaturePlot(Sinoid,features = c("TBX5"),reduction = "tsne")

VlnPlot(Sinoid, features = c("TBX18"))
FeaturePlot(Sinoid,features = c("TBX18"))
FeaturePlot(Sinoid,features = c("TBX18"),reduction = "tsne")

VlnPlot(Sinoid, features = c("NKX2-5"))
FeaturePlot(Sinoid,features = c("NKX2-5"))
FeaturePlot(Sinoid,features = c("NKX2-5"),reduction = "tsne")

table(Idents(Sinoid))

# cluster -- Vascular Smooth Muscle Cells (VSMCs): TAGLN, MYLK, ACTA2, CCN2,MYH11
VlnPlot(Sinoid, features = c("TAGLN"))
FeaturePlot(Sinoid,features = c("MYLK"))
FeaturePlot(Sinoid,features = c("TAGLN"))
FeaturePlot(Sinoid,features = c("ACTA2"))
DotPlot(Sinoid,features = c("TAGLN","MYLK")) + RotatedAxis()

# cluster -- Epithelial: EPCAM,CD24
VlnPlot(Sinoid, features = c("EPCAM"))
FeaturePlot(Sinoid,features = c("EPCAM"))
DotPlot(Sinoid,features = c("EPCAM","CD24")) + RotatedAxis()

# cluster  -- endothelial cells (PECAM1,CDH5,VWF)
DotPlot(Sinoid,features = c("PECAM1","CDH5")) + RotatedAxis()
FeaturePlot(Sinoid,features = c("PECAM1"),reduction = "tsne")
FeaturePlot(Sinoid,features = c("CDH5"),reduction = "tsne")
DotPlot(Sinoid,features = c("PECAM1","CDH5")) + RotatedAxis()

# cluster -- Myoblasts/Fibroblast (DLK1,DCN,LUN)
VlnPlot(Sinoid, features = c("DLK1"))
FeaturePlot(Sinoid,features = c("DLK1"))
FeaturePlot(Sinoid,features = c("DLK1"),reduction = "tsne")
DotPlot(Sinoid,features = c("DLK1","DDR2","PDGFRB")) + RotatedAxis()

# cluster - Proliferating cells (TOP2A,CDK1,CENPF)
# TOP2A Highly expressed in proliferating cardiomyocytes and cardiac progenitor cells
# CDK1 Found in cycling cardiomyocytes and cardiac stem/progenitor cells
# CENPF Expressed in mitotic cells, including dividing cardiomyocytes and progenitor cells
VlnPlot(Sinoid, features = c("TOP2A"))
FeaturePlot(Sinoid,features = c("TOP2A"))
DotPlot(Sinoid,features = c("TOP2A","CDK1","CENPF")) + RotatedAxis()


# cluster  -- epicardial (WT1,TCF21,TBX18)
# WT1 (Wilms Tumor 1): A transcription factor widely recognized as a marker for epicardial cells
# TCF21 
# TBX18 Expressed in the proepicardial organ during embryonic development, TBX18 contributes to the formation of the epicardium
VlnPlot(Sinoid, features = c("WT1"))
FeaturePlot(Sinoid,features = c("WT1"))
FeaturePlot(Sinoid,features = c("WT1"),reduction = "tsne")
DotPlot(Sinoid,features = c("WT1","TCF21","TBX18")) + RotatedAxis()

# cluster  --Cardiac Neuron (STMN2,SST,TAC3,NTRK3,ELAVL4)
# STMN2 Highly expressed in neuronal-like cells of the cardiac conduction system
# SST Found in intrinsic cardiac neurons and conduction system cells, regulating parasympathetic signaling and cardiac autonomic function
# TAC3 Expressed in neuronal and conduction cells, potentially involved in neurotransmitter signaling within the heart
# NTRK3 Found in neuronal-like cells and conduction system cells, playing a role in neurotrophic signaling and cardiac development
# ELAVL4 Expressed in neuronal-like conduction system cells, regulating mRNA stability and synaptic plasticity
VlnPlot(Sinoid, features = c("STMN2"))
FeaturePlot(Sinoid,features = c("STMN2"))
DotPlot(Sinoid,features = c("STMN2")) + RotatedAxis()

# cluster  SAN
DotPlot(Sinoid,features = c("SHOX2","HCN1","CACNA1G","CACNA1D","ISL1","HCN4","TBX5")) + RotatedAxis()

saveRDS(Sinoid, file = "Sinoid.rds")


# add annotations
Sinoid <- RenameIdents(Sinoid, '0' = 'SAN','1' = 'SAN','2' = 'SAN','8' = 'SAN','10' = 'SAN')

Sinoid <- RenameIdents(Sinoid, '6' = 'Epithelial')

Sinoid <- RenameIdents(Sinoid, '9' = 'Endothelial')

Sinoid <- RenameIdents(Sinoid, '3' = 'Epicardial')

Sinoid <- RenameIdents(Sinoid, '5' = 'Cardiac Neuron')

Sinoid <- RenameIdents(Sinoid, '4' = 'Proliferating Cells') 

Sinoid$celltype <- Idents(Sinoid)

UMAPPlot(Sinoid,label=FALSE)

UMAPPlot(Sinoid,label=TRUE)

TSNEPlot(Sinoid,label=TRUE)
levels(Sinoid) <- c("SAN","Proliferating Cells","Cardiac Neuron","Epithelial","Epicardial","Endothelial")
DotPlot(Sinoid,features = c("HCN4","TOP2A","STMN2","EPCAM","WT1","CDH5")) + RotatedAxis()
levels(Sinoid) <- c("Endothelial","Epicardial","Epithelial","Cardiac Neuron","Proliferating Cells","SAN")
VlnPlot(Sinoid, features = c("STMN2"),pt.size = 0.5)
FeaturePlot(Sinoid,features = c("STMN2"))
VlnPlot(Sinoid,features = "TOP2A",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("TOP2A"))
VlnPlot(Sinoid,features = "WT1",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("WT1"))
VlnPlot(Sinoid,features = "CDH5",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("CDH5"))
VlnPlot(Sinoid,features = "EPCAM",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("EPCAM"))
VlnPlot(Sinoid,features = "HCN4",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("HCN4"))

DotPlot(Sinoid,features = c("TBX5","CACNA1D","CACNA1G","ISL1","SHOX2","HCN1")) + RotatedAxis()
VlnPlot(Sinoid,features = "SHOX2",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("SHOX2"))
VlnPlot(Sinoid,features = "HCN1",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("HCN1"))
VlnPlot(Sinoid,features = "CACNA1G",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("CACNA1G"))
VlnPlot(Sinoid,features = "CACNA1D",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("CACNA1D"))
VlnPlot(Sinoid,features = "ISL1",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("ISL1"))
VlnPlot(Sinoid,features = "GJC1",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("GJC1"))
VlnPlot(Sinoid,features = "TBX5",pt.size = 0.5)
FeaturePlot(Sinoid,features = c("TBX5"))

DotPlot(Sinoid,features = c("STMN2","DLK1","ROBO2","BST2")) + RotatedAxis()
VlnPlot(Sinoid, features = c("STMN2"),pt.size = 0.5)
VlnPlot(Sinoid, features = c("DLK1"),pt.size = 0.5)
VlnPlot(Sinoid, features = c("ROBO2"),pt.size = 0.5)
VlnPlot(Sinoid, features = c("BST2"),pt.size = 0.5)
FeaturePlot(Sinoid,features = c("STMN2"))
FeaturePlot(Sinoid,features = c("DLK1"))
FeaturePlot(Sinoid,features = c("ROBO2"))
FeaturePlot(Sinoid,features = c("BST2"))

saveRDS(Sinoid, file = "Sinoid_annotated.rds")


# Subcluster SAN population
SAN <- subset(Sinoid,subset = celltype == 'SAN')

SAN <- FindNeighbors(SAN, dims = 1:20)
SAN <- FindClusters(SAN, resolution = 0.5)
# resolution can adjust from 0.2 to 1.5

head(Idents(SAN), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
SAN <- RunUMAP(SAN, dims = 1:20)

UMAPPlot(SAN,label=TRUE)

dotplot <- DotPlot(SAN,features = c("SHOX2","TBX18","NKX2-5","NPPA")) + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  c("SAN-head","SAN-tail","SAN-TZ"))

FeaturePlot(SAN,features = c("SHOX2","TBX18","NKX2-5","NPPA"))
FeaturePlot(SAN,features = "HCN4")
VlnPlot(SAN,features = "SHOX2")
FeaturePlot(SAN,features = "SHOX2")
FeaturePlot(SAN,features = "TBX18")
FeaturePlot(SAN,features = "NKX2-5")
FeaturePlot(SAN,features = "NPPA")

# add annotations
SAN <- RenameIdents(SAN, '1' = 'SAN-tail','4' = 'SAN-tail','5' = 'SAN-tail')
SAN <- RenameIdents(SAN, '3' = 'SAN-head')
SAN <- RenameIdents(SAN, '0' = 'SAN-TZ','2' = 'SAN-TZ')

SAN$celltype <- Idents(SAN)

UMAPPlot(SAN,label=TRUE)

VlnPlot(SAN,features = "HCN4")
VlnPlot(SAN,features = "SHOX2")
VlnPlot(SAN,features = "TBX18")
VlnPlot(SAN,features = "NKX2-5")
VlnPlot(SAN,features = "NPPA")

saveRDS(SAN, file = "Sinoid_SAN_Subcluster.rds")

DotPlot(SAN,features = c("L3MBTL4","CALD1","ZNF385B","DAPK2","NBPF24","CPNE5","PTPRK","CHST11","KCNH7","MYOM2","SDK1","TRPM3")) + RotatedAxis()

DotPlot(SAN,features = c("TBX5","TEAD1","MEF2A","TBX20","MEF2D","ESRRB","NFIA","GATA4","GATA6","MEF2C","NR2F2","TGIF2","MITF","MGA","ESRRA")) + RotatedAxis()
FeaturePlot(SAN,features = "NFIA")
FeaturePlot(SAN,features = "TBX5")
FeaturePlot(SAN,features = "GATA6")


head_vs_tail_TZ <- FindMarkers(SAN,ident.1 = 'SAN-head',ident.2 = c('SAN-tail','SAN-TZ'), min.pct = 0.1)
write.csv(head_vs_tail_TZ,file = 'SAN_head_vs_tail_TZ_DEA.csv')

head_vs_tail <- FindMarkers(SAN,ident.1 = 'SAN-head',ident.2 = 'SAN-tail', min.pct = 0.1)
write.csv(head_vs_tail,file = 'SAN_head_vs_tail_DEA.csv')

tail_vs_head_TZ <- FindMarkers(SAN,ident.1 = 'SAN-tail',ident.2 = c('SAN-TZ','SAN-head'), min.pct = 0.1)
write.csv(tail_vs_head_TZ,file = 'SAN_tail_vs_TZ_head_DEA.csv')

TZ_vs_head_tail <- FindMarkers(SAN,ident.1 = 'SAN-TZ',ident.2 = c('SAN-head','SAN-tail'), min.pct = 0.1)
write.csv(TZ_vs_head_tail,file = 'SAN_TZ_vs_head_tail_DEA.csv')

tail_vs_TZ <- FindMarkers(SAN,ident.1 = 'SAN-tail',ident.2 = 'SAN-TZ', min.pct = 0.1)
write.csv(tail_vs_TZ,file = 'SAN_tail_vs_TZ_DEA.csv')

# SAN-head vs SAN-tail
EnhancedVolcano(head_vs_tail,
                lab = rownames(head_vs_tail),
                x = "avg_log2FC",
                y = "p_val",
                pCutoff = 1e-50,
                FCcutoff = 1,
                xlim = c(-4,4),
                ylim = c(0,350),
                pointSize = 1,
                labSize = 4,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle = NULL)

# SAN-tail vs SAN-TZ
EnhancedVolcano(tail_vs_TZ,
                lab = rownames(tail_vs_TZ),
                x = "avg_log2FC",
                y = "p_val",
                pCutoff = 1e-25,
                FCcutoff = 1,
                xlim = c(-4,4),
                ylim = c(0,150),
                pointSize = 1,
                labSize = 4,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle = NULL)



