library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
library(corrplot)

coculture.data <- Read10X(data.dir = "Coculture_with_intron/Cellranger/")

coculture <- CreateSeuratObject(counts = coculture.data, project = "coculture", min.cells = 3, min.features = 200)

coculture


mixture.data <- Read10X(data.dir = "mixed_with_intron/Cellranger/")

mixture <- CreateSeuratObject(counts = mixture.data, project = "mixture", min.cells = 3, min.features = 200)

mixture


Combine <-merge(x = coculture, y = mixture ,add.cell.ids = c("coculture","mixture"),project = "combine")

Combine


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

Combine <- subset(Combine, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 10)

Combine <- NormalizeData(Combine, normalization.method = "LogNormalize", scale.factor = 10000)

Combine <- FindVariableFeatures(Combine, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Combine), 10)

plot3 <- VariableFeaturePlot(Combine)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

#Scaling the data
all.genes <- rownames(Combine)
Combine <- ScaleData(Combine, features = all.genes)

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

Combine <- JoinLayers(Combine)

# Determine the 'dimensionality' of the dataset
# Combine <- JackStraw(Combine, num.replicate = 100)
# Combine <- ScoreJackStraw(Combine, dims = 1:20)

# JackStrawPlot(Combine, dims = 1:20,ymax = 0.8)

ElbowPlot(Combine)

#Cluster the cells
Combine <- FindNeighbors(Combine, dims = 1:18)
Combine <- FindClusters(Combine, resolution = 0.5)
# resolution can adjust from 0.2 to 1.5

head(Idents(Combine), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
Combine <- RunUMAP(Combine, dims = 1:18)

DimPlot(Combine, reduction = "umap")

Combine <- RunTSNE(Combine, dims = 1:18)

DimPlot(Combine, reduction = "tsne")

UMAPPlot(Combine,label=TRUE)
UMAPPlot(Combine,label=TRUE,split.by = 'orig.ident')
UMAPPlot(Combine,group.by = 'orig.ident')

TSNEPlot(Combine,label=TRUE)

# find all markers of cluster 6,5
cluster0_vs_all <- FindMarkers(Combine, ident.1 = 0, min.pct = 0.1)
head(cluster0_vs_all, n = 20)

cluster0_vs_7 <- FindMarkers(Combine, ident.1 = 0,ident.2 = 7, min.pct = 0.1)
head(cluster0_vs_7, n = 20)


# find markers for every cluster compared to all remaining cells, report only the positive ones
Combine.markers <- FindAllMarkers(Combine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combine.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- Combine.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(Combine, features = top15$gene) + NoLegend()

#find each cluster characters
VlnPlot(Combine, features = c("SHOX2"))
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

# cluster 2--LMOD2 (Leiomodin 2): Essential for actin filament formation in cardiac muscle. Highly expressed in ventricular and atrial cardiomyocytes.
# cluster 2--MYH7 (Myosin Heavy Chain 7, β-MHC): The dominant myosin heavy chain isoform in adult ventricular cardiomyocytes (especially in slow-contracting or failing hearts). High in fetal, hypertrophic, and failing hearts
# cluster 3--KCNQ5, Encodes a potassium channel involved in electrical excitability and action potential regulation, Expressed in ventricular and atrial cardiomyocytes
# cluster 5--NPPB (Natriuretic Peptide B, also known as BNP), Highly expressed in ventricular cardiomyocytes, particularly in response to ventricular stretch, hypertrophy, or heart failure.
VlnPlot(Combine, features = c("MYH7"))
FeaturePlot(Combine,features = c("MYH7"))
FeaturePlot(Combine,features = c("MYH7"))
DotPlot(Combine,features = c("MYH7","NPPB")) + RotatedAxis()

# cluster 1 and 13 -- SAN pacemaker cells
VlnPlot(Combine, features = c("SHOX2"))
FeaturePlot(Combine,features = c("SHOX2"))
DotPlot(Combine,features = c("SHOX2","HCN1")) + RotatedAxis()
DotPlot(Combine,features = c("SHOX2","HCN4","HCN1","TBX5","CACNA1D")) + RotatedAxis()

# cluster 4,9,10-- Epithelial cells (EPCAM,CDH1)
VlnPlot(Combine, features = c("EPCAM"))
FeaturePlot(Combine,features = c("EPCAM"))
FeaturePlot(Combine,features = c("EPCAM"),reduction = "tsne")
DotPlot(Combine,features = c("EPCAM","CDH1")) + RotatedAxis()


# cluster 6 -- Endothelial cells (PECAM1,CDH5,VWF)
DotPlot(Combine,features = c("PECAM1","CDH5")) + RotatedAxis()
FeaturePlot(Combine,features = c("PECAM1"),reduction = "tsne")
FeaturePlot(Combine,features = c("CDH5"),reduction = "tsne")
DotPlot(Combine,features = c("PECAM1","CDH5")) + RotatedAxis()

# cluster 0, 12-- Myofibroblasts (COL1A1,COL1A2,COL3A1,PDGFRB)
# COL1A1,COL1A2,COL3A1 Highly expressed in myofibroblasts key components of the extracellular matrix (ECM),
VlnPlot(Combine, features = c("COL1A1"))
FeaturePlot(Combine,features = c("COL1A1"))
FeaturePlot(Combine,features = c("COL1A1"),reduction = "tsne")
DotPlot(Combine,features = c("COL1A1","COL3A1","COL1A2")) + RotatedAxis()

# cluster 7-- Proliferating cells (TOP2A,CDK1,CENPF)
# TOP2A Highly expressed in proliferating cardiomyocytes and cardiac progenitor cells
# CDK1 Found in cycling cardiomyocytes and cardiac stem/progenitor cells
# CENPF Expressed in mitotic cells, including dividing cardiomyocytes and progenitor cells
VlnPlot(Combine, features = c("TOP2A"))
FeaturePlot(Combine,features = c("TOP2A"))
DotPlot(Combine,features = c("TOP2A","CDK1","CENPF")) + RotatedAxis()

# cluster 11--Cardiac Fibroblasts, Mesenchymal Cells
# ZFHX4, TWIST1, and TWIST2 are strongly linked to fibrosis and tissue remodeling, particularly in response to cardiac injury or stress
VlnPlot(Combine, features = c("ZFHX4"))
FeaturePlot(Combine,features = c("ZFHX4"))
DotPlot(Combine,features = c("ZFHX4","TWIST1","TWIST2")) + RotatedAxis()

# cluster 8--Cardiac Neuron (UNC5C,DCLK1,NPAS3,FOXP2)
# UNC5C Expressed in intrinsic cardiac neurons
# DCLK1 Found in neuronal-like cells within the cardiac conduction system
# NPAS3 Expressed in neuronal-like and conduction system cells
# FOXP2 Found in neuronal-like cardiac cells, possibly autonomic ganglia
VlnPlot(Combine, features = c("DCLK1"))
FeaturePlot(Combine,features = c("DCLK1"))
DotPlot(Combine,features = c("UNC5C","DCLK1","NPAS3","FOXP2")) + RotatedAxis()

# cluster 14--Cardiac Neuron (STMN2,SST,TAC3,NTRK3,ELAVL4)
# STMN2 Highly expressed in neuronal-like cells of the cardiac conduction system
# SST Found in intrinsic cardiac neurons and conduction system cells, regulating parasympathetic signaling and cardiac autonomic function
# TAC3 Expressed in neuronal and conduction cells, potentially involved in neurotransmitter signaling within the heart
# NTRK3 Found in neuronal-like cells and conduction system cells, playing a role in neurotrophic signaling and cardiac development
# ELAVL4 Expressed in neuronal-like conduction system cells, regulating mRNA stability and synaptic plasticity
VlnPlot(Combine, features = c("STMN2"))
FeaturePlot(Combine,features = c("STMN2"))
DotPlot(Combine,features = c("STMN2","SST","TAC3","NTRK3","ELAVL4")) + RotatedAxis()

DotPlot(Combine,features = c("TH","PRPH","CHAT","PHOX2B","RET")) + RotatedAxis()
VlnPlot(Combine, features = c("PHOX2B"),pt.size = 0)
FeaturePlot(Combine,features = c("PHOX2B"))
VlnPlot(Combine, features = c("PRPH"),pt.size = 0)
FeaturePlot(Combine,features = c("PRPH"))

# cluster not detected -- epicardial (WT1,TBX18) WT1, TBX18, SEMA3D, ALDH1A2, GATA5, TCF21
VlnPlot(Combine, features = c("WT1"))
FeaturePlot(Combine,features = c("WT1"))
FeaturePlot(Combine,features = c("WT1"),reduction = "tsne")

VlnPlot(Combine, features = c("BNC1"))
FeaturePlot(Combine,features = c("BNC1"))
FeaturePlot(Combine,features = c("BNC1"),reduction = "tsne")

DotPlot(Combine,features = c("WT1")) + RotatedAxis()

saveRDS(Combine, file = "Combine.rds")


# add annotations
Combine <- RenameIdents(Combine, '0' = 'Myofibroblasts','8' = 'Myofibroblasts','12' = 'Myofibroblasts')

Combine <- RenameIdents(Combine, '1' = 'SAN','13' = 'SAN')

Combine <- RenameIdents(Combine, '2' = 'CM','3' = 'CM','5' = 'CM')

Combine <- RenameIdents(Combine, '4' = 'Epithelial','9' = 'Epithelial','10' = 'Epithelial')

Combine <- RenameIdents(Combine, '6' = 'Endothelial')

Combine <- RenameIdents(Combine, '7' = 'Proliferating Cells')

Combine <- RenameIdents(Combine, '11' = 'Mesenchymal Cells')

Combine <- RenameIdents(Combine, '14' = 'Cardiac Neuron')

Combine$celltype <- Idents(Combine)

UMAPPlot(Combine,label=FALSE)

UMAPPlot(Combine,label=TRUE)

UMAPPlot(Combine,label=TRUE,split.by = "orig.ident")

DotPlot(Combine,features = c("STMN2","SST","TAC3","NTRK3","ELAVL4",
                             "ZFHX4","TWIST1","TWIST2",
                             "UNC5C","DCLK1","NPAS3","FOXP2",
                             "TOP2A","CDK1","CENPF",
                             "PECAM1","CDH5",
                             "EPCAM","CDH1",
                             "MYH7","NPPB","MYL2",
                             "SHOX2","HCN1","TBX5",
                             "COL1A1","COL3A1","COL1A2")) + RotatedAxis()
DotPlot(Combine,features = c("STMN2","ZFHX4","TOP2A","PECAM1","EPCAM","MYH7","HCN4","DLK1")) + RotatedAxis()
DotPlot(Combine,features = c("STMN2","ZFHX4","TOP2A","PECAM1","EPCAM","MYH7","HCN4","COL1A1")) + RotatedAxis()
levels(Combine) <- c("SAN","Proliferating Cells","Cardiac Neuron","Myofibroblasts","Mesenchymal Cells","Epithelial","Endothelial","CM")
DotPlot(Combine,features = c("HCN4","TOP2A","STMN2","COL1A1","ZFHX4","EPCAM","PECAM1","MYH7")) + RotatedAxis()

levels(Combine) <- c("CM","Endothelial","Epithelial","Mesenchymal Cells","Myofibroblasts","Cardiac Neuron","Proliferating Cells","SAN")
VlnPlot(Combine,features = "STMN2",pt.size = 0)
FeaturePlot(Combine,features = "STMN2")
VlnPlot(Combine,features = "ZFHX4",pt.size = 0)
FeaturePlot(Combine,features = "ZFHX4")
VlnPlot(Combine,features = "TOP2A",pt.size = 0)
FeaturePlot(Combine,features = "TOP2A")
VlnPlot(Combine,features = "PECAM1",pt.size = 0)
FeaturePlot(Combine,features = "PECAM1")
VlnPlot(Combine,features = "EPCAM",pt.size = 0)
FeaturePlot(Combine,features = "EPCAM")
VlnPlot(Combine,features = "MYH7",pt.size = 0)
FeaturePlot(Combine,features = "MYH7")
VlnPlot(Combine,features = "SHOX2",pt.size = 0)
FeaturePlot(Combine,features = "SHOX2")
VlnPlot(Combine,features = "HCN4",pt.size = 0)
FeaturePlot(Combine,features = "HCN4")
VlnPlot(Combine,features = "DLK1",pt.size = 0)
FeaturePlot(Combine,features = "DLK1")
VlnPlot(Combine,features = "COL1A1",pt.size = 0)
FeaturePlot(Combine,features = "COL1A1")

DotPlot(Combine,features = c("SHOX2","HCN4","HCN1","CACNA1D","TBX5")) + RotatedAxis()
VlnPlot(Combine,features = "SHOX2",pt.size = 0)
FeaturePlot(Combine,features = "SHOX2")
VlnPlot(Combine,features = "HCN1",pt.size = 0)
FeaturePlot(Combine,features = "HCN1")
VlnPlot(Combine,features = "CACNA1D",pt.size = 0)
FeaturePlot(Combine,features = "CACNA1D")
VlnPlot(Combine,features = "TBX5",pt.size = 0)
FeaturePlot(Combine,features = "TBX5")

DotPlot(Combine,features = c("MYH7","CSRP3","LMOD2","TNNI3K","ANKRD1")) + RotatedAxis()
VlnPlot(Combine,features = "MYH7",pt.size = 0)
FeaturePlot(Combine,features = "MYH7")
VlnPlot(Combine,features = "CSRP3",pt.size = 0)
FeaturePlot(Combine,features = "CSRP3")
VlnPlot(Combine,features = "LMOD2",pt.size = 0)
FeaturePlot(Combine,features = "LMOD2")
VlnPlot(Combine,features = "TNNI3K",pt.size = 0)
FeaturePlot(Combine,features = "TNNI3K")
VlnPlot(Combine,features = "ANKRD1",pt.size = 0)
FeaturePlot(Combine,features = "ANKRD1")


saveRDS(Combine, file = "Combine_annotated.rds")

# find out the differentially expressed genes between mixed and coculture in SAN cluster
celltype <- Combine$celltype
origin <- Combine$orig.ident

Idents(Combine) <- paste(origin,celltype,sep = '_')
levels(Combine)

CM_coculture_vs_mixture <- FindMarkers(Combine, ident.1 = 'coculture_CM', ident.2 = 'mixture_CM')

SAN_coculture_vs_mixture <- FindMarkers(Combine, ident.1 = 'coculture_SAN', ident.2 = 'mixture_SAN')

EC_coculture_vs_mixture <- FindMarkers(Combine, ident.1 = 'coculture_Endothelial', ident.2 = 'mixture_Endothelial')

write.csv(CM_coculture_vs_mixture,file = "CM_coculture_vs_mixture.csv")

write.csv(SAN_coculture_vs_mixture,file = "SAN_coculture_vs_mixture.csv")

write.csv(EC_coculture_vs_mixture,file = "EC_coculture_vs_mixture.csv")

# remove mitochodria gene
CM_coculture_vs_mixture <- CM_coculture_vs_mixture[-which(grepl("^MT-",rownames(CM_coculture_vs_mixture))),]

SAN_coculture_vs_mixture <- SAN_coculture_vs_mixture[-which(grepl("^MT-",rownames(SAN_coculture_vs_mixture))),]

EC_coculture_vs_mixture <- EC_coculture_vs_mixture[-which(grepl("^MT-",rownames(EC_coculture_vs_mixture))),]

EnhancedVolcano(CM_coculture_vs_mixture,
                lab = rownames(CM_coculture_vs_mixture),
                selectLab = rownames(CM_coculture_vs_mixture[(CM_coculture_vs_mixture$avg_log2FC > 0.5 | CM_coculture_vs_mixture$avg_log2FC < -1) & CM_coculture_vs_mixture$p_val_adj < 1e-25, ]),
                x = "avg_log2FC",
                y = "p_val",
                pCutoff = 1e-25,
                FCcutoff = 0.5,
                xlim = c(-5,4),
                ylim = c(0,180),
                pointSize = 1,
                labSize = 5,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle = NULL)

EnhancedVolcano(SAN_coculture_vs_mixture,
                lab = rownames(SAN_coculture_vs_mixture),
                x = "avg_log2FC",
                y = "p_val",
                pCutoff = 1e-10,
                FCcutoff = 0.5,
                xlim = c(-4,4),
                ylim = c(0,120),
                pointSize = 1,
                labSize = 4,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle = NULL)

EnhancedVolcano(EC_coculture_vs_mixture,
                lab = rownames(EC_coculture_vs_mixture),
                x = "avg_log2FC",
                y = "p_val",
                pCutoff = 1e-10,
                FCcutoff = 0.5,
                xlim = c(-4,4),
                ylim = c(0,50),
                pointSize = 1,
                labSize = 4,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle = NULL)

# subclustering of CM + SAN
SAN_CM <- subset(Combine,celltype == 'SAN'| celltype == 'CM')
#Cluster the cells
SAN_CM <- FindNeighbors(SAN_CM, dims = 1:18)
SAN_CM <- FindClusters(SAN_CM, resolution = 1)
# resolution can adjust from 0.2 to 1.5
head(Idents(SAN_CM), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
SAN_CM <- RunUMAP(SAN_CM, dims = 1:18)

SAN_CM <- RunTSNE(SAN_CM, dims = 1:18)

UMAPPlot(SAN_CM,label=TRUE)

UMAPPlot(SAN_CM,label=TRUE,split.by = "orig.ident")

saveRDS(SAN_CM, file = "SAN_CM_subpopulation.rds")

FeaturePlot(SAN_CM,features = 'SHOX2')
VlnPlot(SAN_CM,features = 'SHOX2')

FeaturePlot(SAN_CM,features = 'NKX2-5')
VlnPlot(SAN_CM,features = 'NKX2-5')

FeaturePlot(SAN_CM,features = 'NPPA')
VlnPlot(SAN_CM,features = 'NPPA')

DotPlot(SAN_CM,features = c("SHOX2","NKX2-5","NPPA","TBX18","ISL1")) + RotatedAxis()

SAN_CM <- RenameIdents(SAN_CM,'5' = 'SAN-head','11' = 'SAN-head','12' = 'SAN-head',
                       '1' = 'SAN-tail','3' = 'SAN-tail','9' = 'SAN-TZ',
                       '0' = 'CM','2' = 'CM', '4' = 'CM', '6' = 'CM', '7' = 'CM', '8' = 'CM', '10' = 'CM')

UMAPPlot(SAN_CM,label=TRUE)

DotPlot(SAN_CM,features = c("HCN4","SHOX2","NKX2-5","NPPA","TBX18","ISL1")) + RotatedAxis()

UMAPPlot(SAN_CM,label=TRUE,split.by = 'orig.ident')

saveRDS(SAN_CM, file = "SAN_CM_subpopulation_annotated.rds")


# check the subpopulation composition between cocuture and mixture
sub_ratio <- table(SAN_CM$orig.ident,Idents(SAN_CM))

sub_ratio <- sub_ratio[,1:3]

sub_ratio <- t(prop.table(sub_ratio,margin = 1))


cluster <- rep(rownames(sub_ratio),each = dim(sub_ratio)[2])
group <- rep(colnames(sub_ratio),times = dim(sub_ratio)[1])
ratio <- c(t(sub_ratio))

combine_data <- data.frame(cluster,group,ratio)

# Stacked + percent
dist_plot<- ggplot(combine_data, aes(fill=cluster, y=ratio, x=group)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dist_plot + RotatedAxis()


