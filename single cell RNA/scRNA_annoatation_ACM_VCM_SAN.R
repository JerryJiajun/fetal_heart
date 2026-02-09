library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)

SAN.data <- Read10X(data.dir = "SAN_Control/filtered_feature_bc_matrix/")

SAN <- CreateSeuratObject(counts = SAN.data, project = "SAN", min.cells = 3, min.features = 200)

SAN

ACM.data <- Read10X(data.dir = "ACM/filtered_feature_bc_matrix/")

ACM <- CreateSeuratObject(counts = ACM.data, project = "ACM", min.cells = 3, min.features = 200)

ACM

VCM.data <- Read10X(data.dir = "VCM/filtered_feature_bc_matrix/")

VCM <- CreateSeuratObject(counts = VCM.data, project = "VCM", min.cells = 3, min.features = 200)

VCM

SAN_PCO.data <- Read10X(data.dir = "SAN-PCO/filtered_feature_bc_matrix/")

SAN_PCO <- CreateSeuratObject(counts = SAN_PCO.data, project = "SAN_PCO", min.cells = 3, min.features = 200)

SAN_PCO

Combine <- merge(x = SAN_PCO, y = c(SAN,ACM,VCM) ,add.cell.ids = c("SAN_PCO","SAN","ACM","VCM"),project = "combine")


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

Combine <- subset(Combine, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 20)

Combine <- NormalizeData(Combine, normalization.method = "LogNormalize", scale.factor = 10000)

Combine <- FindVariableFeatures(Combine, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Combine), 10)

plot3 <- VariableFeaturePlot(Combine)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

# Joinlayers
Combine <- JoinLayers(Combine)

# Scaling the data
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
cluster0_vs_all <- FindMarkers(Combine, ident.1 = 0, min.pct = 0.1)
head(cluster0_vs_all, n = 20)

cluster0_vs_7 <- FindMarkers(Combine, ident.1 = 0,ident.2 = 7, min.pct = 0.1)
head(cluster0_vs_7, n = 20)


# find markers for every cluster compared to all remaining cells, report only the positive ones
Combine.markers <- FindAllMarkers(Combine, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combine.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- Combine.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(Combine, features = top20$gene) + NoLegend()
write.csv(Combine.markers, file = "ACM_VCM_SAN_cluster_markers.csv")

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

# Atrial cluster 2,12: NPPA,MYL7,NR2F2,KCNA5,GJA5
VlnPlot(Combine, features = c("NPPA"),pt.size = 0)
FeaturePlot(Combine,features = c("NPPA"))
DotPlot(Combine,features = c("NPPA","NR2F2","GJA5")) + RotatedAxis()

# Ventricle cluster 6,13: MYL2,IRX4,HEY2,SCN5A,GJA1
VlnPlot(Combine, features = c("MYL2"),pt.size = 0)
FeaturePlot(Combine,features = c("MYL2"))
DotPlot(Combine,features = c("MYL2","IRX4","HEY2","SCN5A")) + RotatedAxis()

# Fibroblast cluster 3,5: POSTN,DCN,ACTA2,COL1A1,PDGFRA
FeaturePlot(Combine,features = c("DCN"))
DotPlot(Combine,features = c("DCN","POSTN","COL1A1","PDGFRA","COL3A1")) + RotatedAxis()

# cluster 5 -- Combine pacemaker cells
VlnPlot(Combine, features = c("SHOX2"),pt.size = 0)
FeaturePlot(Combine,features = c("SHOX2"))
DotPlot(Combine,features = c("SHOX2","HCN1")) + RotatedAxis()
DotPlot(Combine,features = c("SHOX2","HCN4","HCN1","TBX5","CACNA1G")) + RotatedAxis()

# cluster 4,10,13,19,17-- Epithelial cells (EPCAM,CDH1)
VlnPlot(Combine, features = c("EPCAM"),pt.size = 0)
FeaturePlot(Combine,features = c("EPCAM"))
FeaturePlot(Combine,features = c("EPCAM"),reduction = "tsne")
DotPlot(Combine,features = c("EPCAM","CDH1")) + RotatedAxis()

# cluster 8 -- Endothelial cells (PECAM1,CDH5,VWF)
VlnPlot(Combine, features = c("PECAM1"),pt.size = 0)
DotPlot(Combine,features = c("PECAM1","CDH5")) + RotatedAxis()
FeaturePlot(Combine,features = c("CDH5"))
FeaturePlot(Combine,features = c("CDH5"),reduction = "tsne")
DotPlot(Combine,features = c("PECAM1","CDH5")) + RotatedAxis()

# cluster 11-- Proliferating cells (TOP2A,CDK1,CENPF)
# TOP2A Highly expressed in proliferating cardiomyocytes and cardiac progenitor cells
# CDK1 Found in cycling cardiomyocytes and cardiac stem/progenitor cells
# CENPF Expressed in mitotic cells, including dividing cardiomyocytes and progenitor cells
VlnPlot(Combine, features = c("TOP2A"),pt.size = 0)
FeaturePlot(Combine,features = c("TOP2A"))
DotPlot(Combine,features = c("TOP2A","CDK1","CENPF")) + RotatedAxis()

# cluster 21--Cardiac Neuron (STMN2,SST,TAC3,NTRK3,ELAVL4)
# STMN2 Highly expressed in neuronal-like cells of the cardiac conduction system
# SST Found in intrinsic cardiac neurons and conduction system cells, regulating parasympathetic signaling and cardiac autonomic function
# TAC3 Expressed in neuronal and conduction cells, potentially involved in neurotransmitter signaling within the heart
# NTRK3 Found in neuronal-like cells and conduction system cells, playing a role in neurotrophic signaling and cardiac development
# ELAVL4 Expressed in neuronal-like conduction system cells, regulating mRNA stability and synaptic plasticity
DotPlot(Combine,features = c("STMN2","SST","TAC3","NTRK3","ELAVL4")) + RotatedAxis()

DotPlot(Combine,features = c("TH","PRPH","CHAT","PHOX2B","RET")) + RotatedAxis()

VlnPlot(Combine, features = c("STMN2"),pt.size = 0) 
FeaturePlot(Combine,features = c("STMN2")) + labs(x = 'UMAP1',y = 'UMAP2')

# cluster 18--neural crest
VlnPlot(Combine, features = c("SOX2"),pt.size = 0) 
FeaturePlot(Combine,features = c("SOX2")) + labs(x = 'UMAP1',y = 'UMAP2')

# cluster 16 -- epicardial (WT1,TBX18) WT1, TBX18, SEMA3D, ALDH1A2, GATA5, TCF21
VlnPlot(Combine, features = c("WT1"),pt.size = 0)
FeaturePlot(Combine,features = c("WT1"))
FeaturePlot(Combine,features = c("WT1"),reduction = "tsne")

DotPlot(Combine,features = c("WT1","BNC1")) + RotatedAxis()

saveRDS(Combine, file = "ACM_VCM_SAN_SAN-PCO.rds")

# add annotations
Combine <- RenameIdents(Combine, '5' = 'SAN')

Combine <- RenameIdents(Combine, '0' = 'ACM','1' = 'ACM')

Combine <- RenameIdents(Combine, '7' = 'VCM','12' = 'VCM')

Combine <- RenameIdents(Combine, '4' = 'Epithelial','10' = 'Epithelial','13' = 'Epithelial','17' = 'Epithelial')

Combine <- RenameIdents(Combine, '2' = 'FB_1','14' = 'FB_1','15' = 'FB_1')

Combine <- RenameIdents(Combine, '3' = 'FB_2','6' = 'FB_2','9' = 'FB_2','20' = 'FB_2')

Combine <- RenameIdents(Combine, '8' = 'Endothelial')

Combine <- RenameIdents(Combine, '16' = 'Epicardial')

Combine <- RenameIdents(Combine, '21' = 'Neuronal')

Combine <- RenameIdents(Combine, '18' = 'Neural_Crest')

Combine <- RenameIdents(Combine, '11' = 'Proliferating','19' = 'Proliferating')


Combine$celltype <- Idents(Combine)

levels(Combine) <- c("ACM","Endothelial","Epicardial","Epithelial","FB_1","FB_2",
                     "Neuronal","Neural_Crest","Proliferating","SAN","VCM")

UMAPPlot(Combine,label=TRUE, repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

Combine$orig.ident <- factor(Combine$orig.ident,levels = c('ACM','VCM','SAN','SAN_PCO'))

UMAPPlot(Combine,label=TRUE,split.by = "orig.ident",repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

UMAPPlot(Combine,group.by = 'orig.ident') + labs(x = 'UMAP1',y = 'UMAP2')

DotPlot(Combine,features = c("NPPA","CDH5","WT1","EPCAM","POSTN","PDZRN4","STMN2","SOX2","TOP2A","SHOX2","MYL2")) + RotatedAxis()

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
VlnPlot(Combine,features = "PDZRN4",pt.size = 0)
FeaturePlot(Combine,features = c("PDZRN4")) + labs(x = 'UMAP1',y = 'UMAP2')
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
DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN')),features = c("MYH7","MYL2","MYL3","IRX4","HEY2","ANKRD1","FHL2","MYOZ2","CSRP3")) + RotatedAxis()

DotPlot(Combine,features = c("NPPA","GJA5","NR2F2","OPCML","CSMD1","BMPER","ACOXL")) + RotatedAxis()
DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN')),features = c("NPPA","GJA5","NR2F2","OPCML","CSMD1","BMPER","ACOXL")) + RotatedAxis()

DotPlot(Combine,features = c("SHOX2","ISL1","TBX5","HCN1","CACNA1G",'TBX18')) + RotatedAxis()
DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN')),features = c("SHOX2","ISL1","TBX18","HCN1","CACNA1G",'TBX5')) + RotatedAxis()

DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN')),features = c("NPPA","OPCML","CSMD1","SHOX2","ISL1","HCN1","HEY2","MYL2","ANKRD1")) + RotatedAxis()

DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN')),features = c("NPPA","OPCML","CSMD1","GJA5","ACOXL","SHOX2","ISL1","HCN1","HCN4","TBX5","HEY2","MYL2","ANKRD1","MYOZ2","IRX4")) + RotatedAxis() + labs(x='',y='')

DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN')),features = c("GJA5","KCNJ2","KCNQ1","CACNA1C","HCN1","CACNA1D","SCN5A","GJA1")) + RotatedAxis() + labs(x='',y='')


# Differential expression of ion channel related genes underlies the distinct electrophysiological properties of sinoatrial node (SAN), atrial, and ventricular cardiomyocytes. 
# SAN cells exhibit high HCN1 and CACNA1D expression, supporting funny current–driven diastolic depolarization and Ca²⁺-dependent action potential upstroke in the setting of low SCN5A 
# while reduced GJA1/GJA5 expression limits electrical coupling and slows conduction. 
# In atrial cardiomyocytes, elevated SCN5A and GJA5/GJA1 promote rapid depolarization and efficient intercellular conduction, with KCNH2 facilitating rapid repolarization and short action potentials. 
# Ventricular cardiomyocytes are distinguished by high SCN5A and GJA1 expression for synchronized activation, dominant CACNA1C supporting a prolonged plateau phase, and strong KCNH2 stabilizing the resting membrane potential

DotPlot(subset(Combine,subset = celltype %in% c('ACM','VCM','SAN')),features = c("HCN1","CACNA1D","CACNA1C","SCN5A","GJA1","GJA5","KCNJ2","KCNH2")) + RotatedAxis() + labs(x='',y='')

DotPlot(Combine,features = c("TH","PRPH","CHAT","PHOX2B","RET","PHOX2A")) + RotatedAxis()

saveRDS(Combine, file = "ACM_VCM_SAN_SAN-PCO_annotated.rds")

# find out the cellular compositions in each sample
meta <- Combine@meta.data
meta$origin <- meta$orig.ident
Combine@meta.data <- meta

celltype <- Combine$celltype
origin <- Combine$origin

celltype_table <- table(Combine$origin,Combine$celltype)
celltype_prop_table <- prop.table(celltype_table,margin = 1)

sample <- rep(rownames(celltype_prop_table),each = dim(celltype_prop_table)[2])
cell <- rep(colnames(celltype_prop_table),times = dim(celltype_prop_table)[1])
ratio <- c(t(celltype_prop_table))

combine_data <- data.frame(sample,cell,ratio)

combine_data$sample <- factor(combine_data$sample,levels = c('ACM','VCM','SAN','SAN_PCO'))
combine_data$cell <- factor(combine_data$cell,levels = c("ACM","Endothelial","Epicardial","Epithelial","FB_1","FB_2",
                                                         "Neuronal","Neural_Crest","Proliferating","SAN","VCM"))

# Stacked + percent
ggplot(combine_data, aes(fill=cell, y=ratio, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()

# find out the differentially expressed genes between mixed and coculture in Combine cluster
celltype <- Combine$celltype
origin <- Combine$orig.ident

Idents(Combine) <- paste(origin,celltype,sep = '_')
levels(Combine)

Combine$group <- paste(origin,celltype,sep = '_')

ACM_coculture_DEA <- FindMarkers(Combine, ident.1 = 'SAN_PCO_ACM', ident.2 = 'ACM_ACM')

SAN_coculture_DEA <- FindMarkers(Combine, ident.1 = c('SAN_PCO_SAN'), ident.2 = c('SAN_SAN'))

# remove mitochodria gene
ACM_coculture_DEA <- ACM_coculture_DEA[-which(grepl("^MT-",rownames(ACM_coculture_DEA))),]
write.csv(ACM_coculture_DEA,file = "ACM_coculture_vs_ACM_DEA.csv")

SAN_coculture_DEA <- SAN_coculture_DEA[-which(grepl("^MT-",rownames(SAN_coculture_DEA))),]

DotPlot(subset(Combine,subset = group %in% c('ACM_ACM','SAN_PCO_ACM')),features = c("STXBP1","SNAP25","SV2B","RAB3A","RAB3C","CPLX1")) + RotatedAxis()

DotPlot(subset(Combine,subset = group %in% c('ACM_ACM','SAN_PCO_ACM')),features = c("APLN","APLNR","GNAQ","GNG2")) + RotatedAxis()

DotPlot(subset(Combine,subset = group %in% c('ACM_ACM','SAN_PCO_ACM')),features = c("STXBP1","SNAP25","SV2B","RAB3A","RAB3C","CPLX1","APLN","APLNR","GNAQ","GNG2")) + RotatedAxis() + labs(x='',y='')

DotPlot(subset(Combine,subset = group %in% c('SAN_Neuronal','SAN_PCO_Neuronal')),features = c("TH","PRPH","CHAT","PHOX2B","RET","PHOX2A")) + RotatedAxis() + labs(x='',y='')

DotPlot(subset(Combine,subset = group %in% c('SAN_Endothelial','SAN_PCO_Endothelial')),features = c("PECAM1","CDH5","VWF","LYVE1","FLT1")) + RotatedAxis() + labs(x='',y='')

DotPlot(subset(Combine,subset = group %in% c('SAN_SAN','SAN_PCO_SAN')),features = c("SHOX2","ISL1","HCN1","CACNA1G",'TBX5')) + RotatedAxis() + labs(x='',y='')

DotPlot(subset(Combine,subset = group %in% c('SAN_Epithelial','SAN_PCO_Epithelial')),features = c("EPCAM","CDH1","CLDN1","FOXA1","KRT5")) + RotatedAxis() + labs(x='',y='')

DotPlot(subset(Combine,subset = group %in% c('SAN_Epicardial','SAN_PCO_Epicardial')),features = c("WT1","TBX18","TCF21","MSLN",'SEMA3D')) + RotatedAxis() + labs(x='',y='')

EnhancedVolcano(ACM_coculture_DEA,
                lab = rownames(ACM_coculture_DEA),
                x = "avg_log2FC",
                y = "p_val",
                pCutoff = 1e-50,
                FCcutoff = 0.5,
                xlim = c(-2,2),
                ylim = c(0,300),
                pointSize = 1,
                labSize = 4,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle = NULL)

EnhancedVolcano(SAN_coculture_DEA,
                lab = rownames(SAN_coculture_DEA),
                x = "avg_log2FC",
                y = "p_val",
                pCutoff = 1e-50,
                FCcutoff = 0.5,
                xlim = c(-10,5),
                ylim = c(0,400),
                pointSize = 1,
                labSize = 4,
                colAlpha = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                title = NULL,
                subtitle = NULL)


