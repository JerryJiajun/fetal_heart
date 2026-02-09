library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

heart.data <- Read10X(data.dir = "mini_heart/filtered_feature_bc_matrix/")

heart <- CreateSeuratObject(counts = heart.data, project = "heart", min.cells = 3, min.features = 200)

heart


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(heart@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

heart <- subset(heart, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 10)

heart <- NormalizeData(heart, normalization.method = "LogNormalize", scale.factor = 10000)

heart <- FindVariableFeatures(heart, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(heart), 10)

plot3 <- VariableFeaturePlot(heart)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

#Scaling the data
all.genes <- rownames(heart)
heart <- ScaleData(heart, features = all.genes)

#'regress out' heterogeneity associated with mitochondrial contamination
heart <- ScaleData(heart, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
heart <- RunPCA(heart, features = VariableFeatures(object = heart))

# Examine and visualize PCA results a few different ways
print(heart[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(heart, dims = 1:2, reduction = "pca")

DimPlot(heart, reduction = "pca")

# DimHeatmap(heart, dims = 1, cells = 500, balanced = TRUE)

# DimHeatmap(heart, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset
# heart <- JackStraw(heart, num.replicate = 100)
# heart <- ScoreJackStraw(heart, dims = 1:20)

# JackStrawPlot(heart, dims = 1:20,ymax = 0.8)

ElbowPlot(heart)

#Cluster the cells
heart <- FindNeighbors(heart, dims = 1:20)
heart <- FindClusters(heart, resolution =0.8)
# resolution can adjust from 0.2 to 1.5

head(Idents(heart), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
heart <- RunUMAP(heart, dims = 1:20)

DimPlot(heart, reduction = "umap")

heart <- RunTSNE(heart, dims = 1:20)

DimPlot(heart, reduction = "tsne")

UMAPPlot(heart,label=TRUE,repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

TSNEPlot(heart,label=TRUE) + labs(x = 'TSNE1',y = 'TSNE2')

# find all markers of cluster 6,5
cluster5_vs_all <- FindMarkers(heart, ident.1 = 5, min.pct = 0.1)
head(cluster5_vs_all, n = 20)

cluster0_vs_7 <- FindMarkers(heart, ident.1 = 0,ident.2 = 7, min.pct = 0.1)
head(cluster0_vs_7, n = 20)


# find markers for every cluster compared to all remaining cells, report only the positive ones
heart.markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
heart.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- heart.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(heart, features = top15$gene) + NoLegend()

#find each cluster characters
VlnPlot(heart, features = c("SHOX2"),pt.size = 0)
FeaturePlot(heart,features = c("SHOX2"))
FeaturePlot(heart,features = c("SHOX2"),reduction = "tsne")

VlnPlot(heart, features = c("MYH6"))
FeaturePlot(heart,features = c("MYH6"))
FeaturePlot(heart,features = c("MYH6"),reduction = "tsne")

VlnPlot(heart, features = c("HCN1"))
FeaturePlot(heart,features = c("HCN1"))
FeaturePlot(heart,features = c("HCN1"),reduction = "tsne")

VlnPlot(heart, features = c("HCN4"))
FeaturePlot(heart,features = c("HCN4"))
FeaturePlot(heart,features = c("HCN4"),reduction = "tsne")

VlnPlot(heart, features = c("CACNA1D"))
FeaturePlot(heart,features = c("CACNA1D"))
FeaturePlot(heart,features = c("CACNA1D"),reduction = "tsne")

VlnPlot(heart, features = c("CACNA1G"))
FeaturePlot(heart,features = c("CACNA1G"))
FeaturePlot(heart,features = c("CACNA1G"),reduction = "tsne")

VlnPlot(heart, features = c("GJC1"))
FeaturePlot(heart,features = c("GJC1"))
FeaturePlot(heart,features = c("GJC1"),reduction = "tsne")

VlnPlot(heart, features = c("ISL1"))
FeaturePlot(heart,features = c("ISL1"))
FeaturePlot(heart,features = c("ISL1"),reduction = "tsne")

VlnPlot(heart, features = c("TBX3"))
FeaturePlot(heart,features = c("TBX3"))
FeaturePlot(heart,features = c("TBX3"),reduction = "tsne")

VlnPlot(heart, features = c("TBX5"))
FeaturePlot(heart,features = c("TBX5"))
FeaturePlot(heart,features = c("TBX5"),reduction = "tsne")

VlnPlot(heart, features = c("TBX18"))
FeaturePlot(heart,features = c("TBX18"))
FeaturePlot(heart,features = c("TBX18"),reduction = "tsne")

VlnPlot(heart, features = c("NKX2-5"))
FeaturePlot(heart,features = c("NKX2-5"))
FeaturePlot(heart,features = c("NKX2-5"),reduction = "tsne")

table(Idents(heart))

# Atrial cluster 1,9,13: NPPA,MYL7,NR2F2,KCNA5,GJA5
VlnPlot(heart, features = c("NPPA"),pt.size = 0)
FeaturePlot(heart,features = c("NPPA"))+ labs(x = 'UMAP1',y = 'UMAP2')
DotPlot(heart,features = c("NPPA","NR2F2","GJA5")) + RotatedAxis()

# Ventricle cluster 0,2,3: MYL2,IRX4,HEY2,SCN5A,GJA1
VlnPlot(heart, features = c("HEY2"),pt.size = 0)
VlnPlot(heart, features = c("HAND1"),pt.size = 0)
FeaturePlot(heart,features = c("HEY2"))+ labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(heart,features = c("HAND1"))+ labs(x = 'UMAP1',y = 'UMAP2')
DotPlot(heart,features = c("MYL2","IRX4","HEY2","HAND1")) + RotatedAxis()

# Fibroblast cluster 4: POSTN,DCN,ACTA2,COL1A1,PDGFRA
FeaturePlot(heart,features = c("DCN"))+ labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = c("POSTN"),pt.size = 0)
DotPlot(heart,features = c("DCN","POSTN","COL1A1","PDGFRA","COL3A1")) + RotatedAxis()

# cluster 6.13 -- heart pacemaker cells
VlnPlot(heart, features = c("SHOX2"),pt.size = 0)
FeaturePlot(heart,features = c("SHOX2"))+ labs(x = 'UMAP1',y = 'UMAP2')
DotPlot(heart,features = c("SHOX2","HCN1")) + RotatedAxis()
DotPlot(heart,features = c("SHOX2","HCN4","HCN1","TBX5","CACNA1G")) + RotatedAxis()

# cluster 8-- Epithelial cells (EPCAM,CDH1)
VlnPlot(heart, features = c("EPCAM"),pt.size = 0)
FeaturePlot(heart,features = c("EPCAM"))
FeaturePlot(heart,features = c("EPCAM"),reduction = "tsne")
DotPlot(heart,features = c("EPCAM","CDH1")) + RotatedAxis()

# cluster 11 -- Endothelial cells (PECAM1,CDH5,VWF)
VlnPlot(heart, features = c("PECAM1"),pt.size = 0)
DotPlot(heart,features = c("PECAM1","CDH5")) + RotatedAxis()
FeaturePlot(heart,features = c("CDH5"))+ labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(heart,features = c("CDH5"),reduction = "tsne")
DotPlot(heart,features = c("PECAM1","CDH5")) + RotatedAxis()

# cluster 10-- Proliferating cells (TOP2A,CDK1,CENPF)
# TOP2A Highly expressed in proliferating cardiomyocytes and cardiac progenitor cells
# CDK1 Found in cycling cardiomyocytes and cardiac stem/progenitor cells
# CENPF Expressed in mitotic cells, including dividing cardiomyocytes and progenitor cells
VlnPlot(heart, features = c("TOP2A"),pt.size = 0)
FeaturePlot(heart,features = c("TOP2A"))+ labs(x = 'UMAP1',y = 'UMAP2')
DotPlot(heart,features = c("TOP2A","CDK1","CENPF")) + RotatedAxis()

# cluster 4,12--Cardiac Neuron (STMN2,SST,TAC3,NTRK3,ELAVL4)
# STMN2 Highly expressed in neuronal-like cells of the cardiac conduction system
# SST Found in intrinsic cardiac neurons and conduction system cells, regulating parasympathetic signaling and cardiac autonomic function
# TAC3 Expressed in neuronal and conduction cells, potentially involved in neurotransmitter signaling within the heart
# NTRK3 Found in neuronal-like cells and conduction system cells, playing a role in neurotrophic signaling and cardiac development
# ELAVL4 Expressed in neuronal-like conduction system cells, regulating mRNA stability and synaptic plasticity
DotPlot(heart,features = c("STMN2","SST","TAC3","NTRK3","ELAVL4")) + RotatedAxis()

DotPlot(heart,features = c("PRPH","RET","PHOX2A")) + RotatedAxis()

VlnPlot(heart, features = c("STMN2"),pt.size = 0) 
FeaturePlot(heart,features = c("STMN2")) + labs(x = 'UMAP1',y = 'UMAP2')

# cluster 7 -- epicardial (WT1,TBX18) WT1, TBX18, SEMA3D, ALDH1A2, GATA5, TCF21
VlnPlot(heart, features = c("WT1"),pt.size = 0)
FeaturePlot(heart,features = c("WT1"))+ labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(heart,features = c("WT1"),reduction = "tsne")

DotPlot(heart,features = c("WT1","BNC1")) + RotatedAxis()

VlnPlot(heart, features = c("SOX2"),pt.size = 0) 
FeaturePlot(heart,features = c("SOX2"))

saveRDS(heart, file = "mini_heart.rds")

# add annotations
heart <- RenameIdents(heart, '6' = 'SAN','13' = 'SAN')

heart <- RenameIdents(heart, '10' = 'Proliferating')

heart <- RenameIdents(heart, '0' = 'VCM','2' = 'VCM','3' = 'VCM')

heart <- RenameIdents(heart, '1' = 'ACM','5' = 'ACM','9' = 'ACM','15' = 'ACM')

heart <- RenameIdents(heart, '8' = 'Epithelial')

heart <- RenameIdents(heart, '11' = 'Endothelial')

heart <- RenameIdents(heart, '7' = 'Epicardial')

heart <- RenameIdents(heart, '14' = 'FB')

heart <- RenameIdents(heart, '4' = 'Neuronal','12' = 'Neuronal')

heart$celltype <- Idents(heart)

levels(heart) <- c("ACM","Endothelial","Epicardial","Epithelial","FB",
                   "Neuronal","Proliferating","SAN","VCM")

UMAPPlot(heart,label=TRUE, repel = TRUE) + labs(x = 'UMAP1',y = 'UMAP2')

DotPlot(heart,features = c("NPPA","CDH5","WT1","EPCAM","DCN","STMN2","TOP2A","SHOX2","HEY2")) + RotatedAxis()

DotPlot(heart,features = c("NPPA","CDH5","WT1","EPCAM","DCN","STMN2","SOX2","TOP2A","SHOX2","HEY2")) + RotatedAxis()

VlnPlot(heart,features = "STMN2",pt.size = 0)
FeaturePlot(heart,features = c("STMN2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "CDH5",pt.size = 0)
FeaturePlot(heart,features = c("CDH5")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "PECAM1",pt.size = 0)
FeaturePlot(heart,features = c("PECAM1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "WT1",pt.size = 0)
FeaturePlot(heart,features = c("WT1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "BNC1",pt.size = 0)
FeaturePlot(heart,features = c("BNC1")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "EPCAM",pt.size = 0)
FeaturePlot(heart,features = c("EPCAM")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "SHOX2",pt.size = 0)
FeaturePlot(heart,features = c("SHOX2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "HCN4",pt.size = 0)
FeaturePlot(heart,features = c("HCN4")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "POSTN",pt.size = 0)
FeaturePlot(heart,features = c("POSTN")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "DCN",pt.size = 0)
FeaturePlot(heart,features = c("DCN")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "TOP2A",pt.size = 0)
FeaturePlot(heart,features = c("TOP2A")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "MYL2",pt.size = 0)
FeaturePlot(heart,features = c("MYL2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "HEY2",pt.size = 0)
FeaturePlot(heart,features = c("HEY2")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "NPPA",pt.size = 0)
FeaturePlot(heart,features = c("NPPA")) + labs(x = 'UMAP1',y = 'UMAP2')
VlnPlot(heart,features = "MYL7",pt.size = 0)
FeaturePlot(heart,features = c("MYL7")) + labs(x = 'UMAP1',y = 'UMAP2')

DotPlot(heart,features = c("SHOX2","HCN1","ISL1","TBX18","HCN4")) + RotatedAxis()
VlnPlot(heart,features = "SHOX2",pt.size = 0)
FeaturePlot(heart,features = "SHOX2")
VlnPlot(heart,features = "HCN1",pt.size = 0)
FeaturePlot(heart,features = "HCN1")
VlnPlot(heart,features = "CACNA1F",pt.size = 0)
FeaturePlot(heart,features = "CACNA1D")
VlnPlot(heart,features = "TBX5",pt.size = 0)
FeaturePlot(heart,features = "TBX5")

DotPlot(heart,features = c("MYH7","MYL2","MYL3","IRX4","HEY2","ANKRD1","FHL2","MYOZ2","CSRP3")) + RotatedAxis()
DotPlot(subset(heart,subset = celltype %in% c('ACM','VCM','SAN')),features = c("MYL2","IRX4","HEY2","ANKRD1")) + RotatedAxis()

DotPlot(heart,features = c("NPPA","GJA5","NR2F2","OPCML","CSMD1","BMPER","ACOXL")) + RotatedAxis()
DotPlot(subset(heart,subset = celltype %in% c('ACM','VCM','SAN')),features = c("GJA5","OPCML","CSMD1","BMPER","ACOXL")) + RotatedAxis()

DotPlot(heart,features = c("SHOX2","ISL1","TBX18","HCN1","HCN4","TBX5")) + RotatedAxis()
DotPlot(subset(heart,subset = celltype %in% c('ACM','VCM','SAN')),features = c("SHOX2","ISL1","TBX18","HCN1","HCN4","TBX5")) + RotatedAxis()

DotPlot(subset(heart,subset = celltype %in% c('ACM','VCM','SAN')),features = c("MYH7","HEY2","ANKRD1","CSMD1","BMPER","ACOXL","SHOX2","ISL1","TBX18","HCN1","HCN4","TBX5")) + RotatedAxis()

saveRDS(heart, file = "mini_heart_annotated.rds")

# find out the cellular compositions in each sample
meta <- heart@meta.data
meta$origin <- meta$orig.ident
heart@meta.data <- meta

celltype <- heart$celltype
origin <- heart$origin

celltype_table <- table(heart$origin,heart$celltype)
celltype_prop_table <- prop.table(celltype_table,margin = 1)

sample <- rep(rownames(celltype_prop_table),each = dim(celltype_prop_table)[2])
cell <- rep(colnames(celltype_prop_table),times = dim(celltype_prop_table)[1])
ratio <- c(t(celltype_prop_table))

heart_data <- data.frame(sample,cell,ratio)

heart_data$sample <- factor(heart_data$sample,levels = c('ACM','VCM','SAN','SAN_PCO'))
# Stacked + percent
ggplot(heart_data, aes(fill=cell, y=ratio, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()


# find out the differentially expressed genes between mixed and coculture in heart cluster
celltype <- heart$celltype
origin <- heart$orig.ident

Idents(heart) <- paste(origin,celltype,sep = '_')
levels(heart)

ACM_coculture_DEA <- FindMarkers(heart, ident.1 = 'SAN_PCO_ACM', ident.2 = 'ACM_ACM')

SAN_coculture_DEA <- FindMarkers(heart, ident.1 = c('SAN_PCO_SAN'), ident.2 = c('SAN_SAN'))

# remove mitochodria gene
ACM_coculture_DEA <- ACM_coculture_DEA[-which(grepl("^MT-",rownames(ACM_coculture_DEA))),]

SAN_coculture_DEA <- SAN_coculture_DEA[-which(grepl("^MT-",rownames(SAN_coculture_DEA))),]

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



# subclustering of CM + heart
heart_CM <- subset(heart,celltype == 'heart'| celltype == 'CM')
#Cluster the cells
heart_CM <- FindNeighbors(heart_CM, dims = 1:18)
heart_CM <- FindClusters(heart_CM, resolution = 1)
# resolution can adjust from 0.2 to 1.5
head(Idents(heart_CM), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
heart_CM <- RunUMAP(heart_CM, dims = 1:18)

heart_CM <- RunTSNE(heart_CM, dims = 1:18)

UMAPPlot(heart_CM,label=TRUE)

UMAPPlot(heart_CM,label=TRUE,split.by = "orig.ident")

saveRDS(heart_CM, file = "heart_CM_subpopulation.rds")

FeaturePlot(heart_CM,features = 'SHOX2')
VlnPlot(heart_CM,features = 'SHOX2')

FeaturePlot(heart_CM,features = 'NKX2-5')
VlnPlot(heart_CM,features = 'NKX2-5')

FeaturePlot(heart_CM,features = 'NPPA')
VlnPlot(heart_CM,features = 'NPPA')

DotPlot(heart_CM,features = c("SHOX2","NKX2-5","NPPA","TBX18","ISL1")) + RotatedAxis()

heart_CM <- RenameIdents(heart_CM,'5' = 'heart-head','11' = 'heart-head','12' = 'heart-head',
                           '1' = 'heart-tail','3' = 'heart-tail','9' = 'heart-TZ',
                           '0' = 'CM','2' = 'CM', '4' = 'CM', '6' = 'CM', '7' = 'CM', '8' = 'CM', '10' = 'CM')

UMAPPlot(heart_CM,label=TRUE)

DotPlot(heart_CM,features = c("HCN4","SHOX2","NKX2-5","NPPA","TBX18","ISL1")) + RotatedAxis()

UMAPPlot(heart_CM,label=TRUE,split.by = 'orig.ident')

saveRDS(heart_CM, file = "heart_CM_subpopulation_annotated.rds")


# check the subpopulation composition between cocuture and mixture
sub_ratio <- table(heart_CM$orig.ident,Idents(heart_CM))

sub_ratio <- sub_ratio[,1:3]

sub_ratio <- t(prop.table(sub_ratio,margin = 1))


cluster <- rep(rownames(sub_ratio),each = dim(sub_ratio)[2])
group <- rep(colnames(sub_ratio),times = dim(sub_ratio)[1])
ratio <- c(t(sub_ratio))

heart_data <- data.frame(cluster,group,ratio)

# Stacked + percent
dist_plot<- ggplot(heart_data, aes(fill=cluster, y=ratio, x=group)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dist_plot + RotatedAxis()


