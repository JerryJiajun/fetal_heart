library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(EnhancedVolcano)
library(RColorBrewer)
library(pheatmap)
library(corrplot)

SAN_perturb.data <- Read10X(data.dir = "SAN_perturb/STARsolo/filtered/")

SAN_perturb <- CreateSeuratObject(counts = SAN_perturb.data, project = "SAN_perturb", min.cells = 3, min.features = 200)

SAN_perturb


# Lets examine a few genes in the first thirty cells
SAN_perturb.data[c("SHOX2", "TBX3", "TBX5"), 1:30]

dense.size <- object.size(as.matrix(SAN_perturb.data))

dense.size

sparse.size <- object.size(SAN_perturb.data)

sparse.size

dense.size/sparse.size


# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
SAN_perturb[["percent.mt"]] <- PercentageFeatureSet(SAN_perturb, pattern = "^MT-")

# Show QC metrics for the first 5 cells
head(SAN_perturb@meta.data, 5)


# Visualize QC metrics as a violin plot
VlnPlot(SAN_perturb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(SAN_perturb, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SAN_perturb, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

SAN_perturb <- subset(SAN_perturb, subset = nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 20)

cell_barcodes <- rownames(SAN_perturb@meta.data)

names(cell_barcodes) <- colnames(SAN_perturb)

SAN_perturb <- AddMetaData(
  object = SAN_perturb,
  metadata = cell_barcodes,
  col.name = 'cell.barcodes'
)

SAN_perturb <- NormalizeData(SAN_perturb, normalization.method = "LogNormalize", scale.factor = 10000)

SAN_perturb <- FindVariableFeatures(SAN_perturb, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SAN_perturb), 10)

plot3 <- VariableFeaturePlot(SAN_perturb)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4

#Scaling the data
all.genes <- rownames(SAN_perturb)
SAN_perturb <- ScaleData(SAN_perturb, features = all.genes)

#'regress out' heterogeneity associated with mitochondrial contamination
SAN_perturb <- ScaleData(SAN_perturb, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction
SAN_perturb <- RunPCA(SAN_perturb, features = VariableFeatures(object = SAN_perturb))

# Examine and visualize PCA results a few different ways
print(SAN_perturb[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(SAN_perturb, dims = 1:2, reduction = "pca")

DimPlot(SAN_perturb, reduction = "pca")

DimHeatmap(SAN_perturb, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(SAN_perturb, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the 'dimensionality' of the dataset
SAN_perturb <- JackStraw(SAN_perturb, num.replicate = 100)
SAN_perturb <- ScoreJackStraw(SAN_perturb, dims = 1:20)

JackStrawPlot(SAN_perturb, dims = 1:20,ymax = 0.8)

ElbowPlot(SAN_perturb)

#Cluster the cells
SAN_perturb <- FindNeighbors(SAN_perturb, dims = 1:20)
SAN_perturb <- FindClusters(SAN_perturb, resolution = 0.5)
# resolution can adjust from 0.2 to 1.5

head(Idents(SAN_perturb), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
SAN_perturb <- RunUMAP(SAN_perturb, dims = 1:20)

DimPlot(SAN_perturb, reduction = "umap")

SAN_perturb <- RunTSNE(SAN_perturb, dims = 1:20)

DimPlot(SAN_perturb, reduction = "tsne")

UMAPPlot(SAN_perturb,label=TRUE)

TSNEPlot(SAN_perturb,label=TRUE)

# find all markers of cluster 6,5
cluster6_vs_all <- FindMarkers(SAN_perturb, ident.1 = 6, min.pct = 0.1)
head(cluster6_vs_all, n = 20)

cluster5_vs_all <- FindMarkers(SAN_perturb, ident.1 = 5, min.pct = 0.1)
head(cluster5_vs_all, n = 20)

cluster2_vs_all <- FindMarkers(SAN_perturb, ident.1 = 2, min.pct = 0.1)
head(cluster2_vs_all, n = 20)

# find markers for every cluster compared to all remaining cells, report only the positive ones
SAN_perturb.markers <- FindAllMarkers(SAN_perturb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SAN_perturb.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top15 <- SAN_perturb.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
DoHeatmap(SAN_perturb, features = top15$gene) + NoLegend()

#find each cluster characters
VlnPlot(SAN_perturb, features = c("SHOX2"))
FeaturePlot(SAN_perturb,features = c("SHOX2"))
FeaturePlot(SAN_perturb,features = c("SHOX2"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("MYH6"))
FeaturePlot(SAN_perturb,features = c("MYH6"))
FeaturePlot(SAN_perturb,features = c("MYH6"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("HCN1"))
FeaturePlot(SAN_perturb,features = c("HCN1"))
FeaturePlot(SAN_perturb,features = c("HCN1"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("HCN4"))
FeaturePlot(SAN_perturb,features = c("HCN4"))
FeaturePlot(SAN_perturb,features = c("HCN4"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("CACNA1D"))
FeaturePlot(SAN_perturb,features = c("CACNA1D"))
FeaturePlot(SAN_perturb,features = c("CACNA1D"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("CACNA1G"))
FeaturePlot(SAN_perturb,features = c("CACNA1G"))
FeaturePlot(SAN_perturb,features = c("CACNA1G"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("GJC1"))
FeaturePlot(SAN_perturb,features = c("GJC1"))
FeaturePlot(SAN_perturb,features = c("GJC1"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("ISL1"))
FeaturePlot(SAN_perturb,features = c("ISL1"))
FeaturePlot(SAN_perturb,features = c("ISL1"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("TBX3"))
FeaturePlot(SAN_perturb,features = c("TBX3"))
FeaturePlot(SAN_perturb,features = c("TBX3"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("TBX5"))
FeaturePlot(SAN_perturb,features = c("TBX5"))
FeaturePlot(SAN_perturb,features = c("TBX5"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("TBX18"))
FeaturePlot(SAN_perturb,features = c("TBX18"))
FeaturePlot(SAN_perturb,features = c("TBX18"),reduction = "tsne")

VlnPlot(SAN_perturb, features = c("NKX2-5"))
FeaturePlot(SAN_perturb,features = c("NKX2-5"))
FeaturePlot(SAN_perturb,features = c("NKX2-5"),reduction = "tsne")

table(Idents(SAN_perturb))

VlnPlot(SAN_perturb, features = c("EGFP"))
FeaturePlot(SAN_perturb,features = c("EGFP"))

VlnPlot(SAN_perturb, features = c("CAS9"))
FeaturePlot(SAN_perturb,features = c("CAS9"))

VlnPlot(SAN_perturb, features = c("MCHERRY"))
FeaturePlot(SAN_perturb,features = c("MCHERRY"))

# cluster 3-- Vascular Smooth Muscle Cells (VSMCs): TAGLN, MYLK, ACTA2, CCN2
VlnPlot(SAN_perturb, features = c("TAGLN"))
FeaturePlot(SAN_perturb,features = c("MYLK"))
FeaturePlot(SAN_perturb,features = c("TAGLN"))
FeaturePlot(SAN_perturb,features = c("ACTA2"))
DotPlot(SAN_perturb,features = c("TAGLN","MYLK")) + RotatedAxis()

# cluster 5-- Epithelial: EPCAM,CD24
VlnPlot(SAN_perturb, features = c("EPCAM"))
FeaturePlot(SAN_perturb,features = c("EPCAM"))
DotPlot(SAN_perturb,features = c("EPCAM","CD24")) + RotatedAxis()

# cluster 8 -- endothelial cells (PECAM1,CDH5,VWF)
DotPlot(SAN_perturb,features = c("PECAM1","CDH5")) + RotatedAxis()
FeaturePlot(SAN_perturb,features = c("PECAM1"),reduction = "tsne")
FeaturePlot(SAN_perturb,features = c("CDH5"),reduction = "tsne")
DotPlot(SAN_perturb,features = c("PECAM1","CDH5")) + RotatedAxis()

# cluster 8-- Myoblasts/Fibroblast (DLK1,DCN,LUN)
VlnPlot(SAN_perturb, features = c("DLK1"))
FeaturePlot(SAN_perturb,features = c("DLK1"))
FeaturePlot(SAN_perturb,features = c("DLK1"),reduction = "tsne")
DotPlot(SAN_perturb,features = c("DLK1","DDR2","PDGFRB")) + RotatedAxis()

# cluster 6-- Proliferating cells (TOP2A,CDK1,CENPF)
# TOP2A Highly expressed in proliferating cardiomyocytes and cardiac progenitor cells
# CDK1 Found in cycling cardiomyocytes and cardiac stem/progenitor cells
# CENPF Expressed in mitotic cells, including dividing cardiomyocytes and progenitor cells
VlnPlot(SAN_perturb, features = c("TOP2A"))
FeaturePlot(SAN_perturb,features = c("TOP2A"))
DotPlot(SAN_perturb,features = c("TOP2A","CDK1","CENPF")) + RotatedAxis()


# cluster 4 -- epicardial (WT1,TCF21,TBX18)
# WT1 (Wilms Tumor 1): A transcription factor widely recognized as a marker for epicardial cells
# TCF21 
# TBX18 Expressed in the proepicardial organ during embryonic development, TBX18 contributes to the formation of the epicardium
VlnPlot(SAN_perturb, features = c("WT1"))
FeaturePlot(SAN_perturb,features = c("WT1"))
FeaturePlot(SAN_perturb,features = c("WT1"),reduction = "tsne")
DotPlot(SAN_perturb,features = c("WT1","TCF21","TBX18")) + RotatedAxis()

# cluster 1,2 SAN normal, cluster 3,6 SAN low ion channel, cluster 0 SAN low SHOX2
DotPlot(SAN_perturb,features = c("SHOX2","HCN4","HCN1","CACNA1D")) + RotatedAxis()

DotPlot(SAN_perturb,features = c("SHOX2","HCN4","CACNA1D","WT1","BNC1","CD24","EPCAM","PECAM1","CDH5")) + RotatedAxis()

DotPlot(SAN_perturb,features = c("WT1","BNC1","CD24","EPCAM","PECAM1","CDH5","DDR2")) + RotatedAxis()

DotPlot(SAN_perturb,features = c("HCN1","HCN4","CACNA1D","CACNA1G","TBX5","ISL1","SHOX2","NKX2-5")) + RotatedAxis()

DotPlot(SAN_perturb,features = c("SHOX2","HCN1","HCN4","CACNA1D","ISL1","GJC1"),group.by = "orig.ident") + RotatedAxis()

DotPlot(SAN_perturb,c("SHOX2","HCN4",input_perturb_target$target_gene[1:5])) + RotatedAxis()

saveRDS(SAN_perturb, file = "SAN_perturb.rds")


# check the fetal heart identified genes expression in hESC derived SAN
specific_reference <- read.csv(file = "fetal_heart_markers.csv",header = TRUE,row.names = 1)
SAN_specific_reference <- specific_reference %>%
  filter(cluster == '10',avg_log2FC > 1,pct.1 > 0.1) 
DotPlot(SAN_perturb,features = SAN_specific_reference$gene) + RotatedAxis()

DoHeatmap(subset(SAN_perturb,downsample = 500),features = SAN_specific_reference$gene,label = TRUE)


# add annotations
SAN_perturb <- RenameIdents(SAN_perturb, '0' = 'SAN','1' = 'SAN','2' = 'SAN','3' = 'SAN')

SAN_perturb <- RenameIdents(SAN_perturb, '5' = 'Epithelial')

SAN_perturb <- RenameIdents(SAN_perturb, '8' = 'Endothelial')

SAN_perturb <- RenameIdents(SAN_perturb, '4' = 'Epicardial')

SAN_perturb <- RenameIdents(SAN_perturb, '6' = 'Proliferating Cells') 

SAN_perturb <- RenameIdents(SAN_perturb, '7' = 'Progenitors')


SAN_perturb$celltype <- Idents(SAN_perturb)

UMAPPlot(SAN_perturb,label=FALSE)

UMAPPlot(SAN_perturb,label=TRUE)

TSNEPlot(SAN_perturb,label=TRUE)

DotPlot(SAN_perturb,features = c("ISL1","TOP2A","WT1","CDH5","EPCAM","HCN4","HCN1","SHOX2")) + RotatedAxis()
VlnPlot(SAN_perturb,features = "TOP2A")
VlnPlot(SAN_perturb,features = "WT1")
VlnPlot(SAN_perturb,features = "CDH5")
VlnPlot(SAN_perturb,features = "EPCAM")
VlnPlot(SAN_perturb,features = "HCN4")
VlnPlot(SAN_perturb,features = "SHOX2")
FeaturePlot(SAN_perturb,features = "HCN4")
FeaturePlot(SAN_perturb,features = "SHOX2")
FeaturePlot(SAN_perturb,features = "NPPA")

saveRDS(SAN_perturb, file = "SAN_perturb_annotated.rds")

# Subcluster SAN population
SAN <- subset(SAN_perturb,subset = celltype == 'SAN')

SAN <- FindNeighbors(SAN, dims = 1:20)
SAN <- FindClusters(SAN, resolution = 0.5)
# resolution can adjust from 0.2 to 1.5

head(Idents(SAN), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
SAN <- RunUMAP(SAN, dims = 1:20)

UMAPPlot(SAN,label=TRUE)

FeaturePlot(SAN,features = "SHOX2")
# add annotations
SAN <- RenameIdents(SAN, '2' = 'SAN-tail','5' = 'SAN-tail')
SAN <- RenameIdents(SAN, '1' = 'SAN-head')
SAN <- RenameIdents(SAN, '0' = 'SAN-TZ','3' = 'SAN-TZ','4' = 'SAN-TZ','6' = 'SAN-TZ')

SAN$celltype <- Idents(SAN)

UMAPPlot(SAN,label=TRUE)

VlnPlot(SAN,features = "HCN4")
VlnPlot(SAN,features = "SHOX2")
VlnPlot(SAN,features = "TBX18")
VlnPlot(SAN,features = "NKX2-5")
VlnPlot(SAN,features = "NPPA")

saveRDS(SAN, file = "Sinoid_SAN_Subcluster.rds")

# subselect the perturb cells
library(dplyr)
perturb_infor <- read.csv(file = "SAN_CBC_GBC_reference.csv",row.names = 1)
perturb_infor <- perturb_infor[,c(2,8)]
perturb_infor <- perturb_infor %>% group_by(cell_barcode) %>% count(target) %>% filter(n > 4)
UniqueCBC <- perturb_infor %>%
  group_by(cell_barcode) %>%
  summarise(counts = n()) %>%
  filter(counts ==1) 
perturb_infor <- perturb_infor[perturb_infor$cell_barcode %in% UniqueCBC$cell_barcode,]

SAN_perturb_GBC <- subset(SAN, subset = cell.barcodes  %in% perturb_infor$cell_barcode)

perturb_target <- perturb_infor$target
names(perturb_target) <- perturb_infor$cell_barcode

SAN_perturb_GBC <- AddMetaData(
  object = SAN_perturb_GBC,
  metadata = perturb_target,
  col.name = 'perturb_target'
)

'%notin%' <- Negate('%in%')

saveRDS(SAN_perturb_GBC, file = "SAN_perturb_GBC.rds")

UMAPPlot(SAN_perturb_GBC,label=FALSE)

TSNEPlot(SAN_perturb_GBC,label=FALSE)

UMAPPlot(SAN_perturb_GBC,group.by = 'perturb_target')

VlnPlot(SAN_perturb_GBC,features = c("SHOX2"),group.by = 'perturb_target')

DoHeatmap(SAN_perturb_GBC,features = c("SHOX2"),group.by = 'perturb_target',angle = 90,size = 2)

DotPlot(SAN_perturb_GBC,features = c("SHOX2","HCN4","HCN1","TBX5","CACNA1G","GJC1","CACNA1D"),group.by = 'perturb_target') + coord_flip() + RotatedAxis()

DotPlot(SAN_perturb_GBC,features = c("SHOX2","TBX18","NKX2-5","NPPA","TBX3","GJC1","GJA1"),group.by = 'perturb_target') + coord_flip() + RotatedAxis()

DotPlot(SAN_perturb_GBC,features = sort(unique(perturb_target)),group.by = 'perturb_target') + RotatedAxis()


perturb_target_ratio <- table(SAN_perturb_GBC$perturb_target,SAN_perturb_GBC$celltype)
perturb_target_ratio <- t(prop.table(perturb_target_ratio,margin = 1))

write.csv(perturb_target_ratio,file = "SAN_perturb_target_ratio.csv",row.names = TRUE)

SAN_dist <- read.csv(file = "SAN_perturb_target_ratio.csv",row.names = 1)

cluster <- rep(rownames(SAN_dist),each = dim(SAN_dist)[2])
perturb_target <- rep(colnames(SAN_dist),times = dim(SAN_dist)[1])
ratio <- c(t(SAN_dist))

combine_data <- data.frame(cluster,perturb_target,ratio)

select_perturb_target <- table(SAN_perturb_GBC$perturb_target)

select_perturb_target <- select_perturb_target[select_perturb_target >= 10]

combine_data_select <- combine_data %>%
  filter(perturb_target %in% names(select_perturb_target))

# Stacked + percent
dist_plot<- ggplot(combine_data, aes(fill=cluster, y=ratio, x=perturb_target)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dist_plot + RotatedAxis()

# Stacked + percent
combine_data_select$perturb_target <- factor(combine_data_select$perturb_target,levels = c('ZNF385B','L3MBTL4','CALD1','ARHGAP24','RELN','PLEKHA7',
                                                                                           'TACC2','TRDN','ATP2A2','CHRM2','EPHA4','NPPA','CMYA5',
                                                                                           'DAPK2','NBPF24','PTPRK','CPNE5','MAST4',
                                                                                           'CHST11','KCNH7','MYOM2','SDK1','TRPM3','SFRP1',
                                                                                           'PDE4D','PLCB4','PRKAG2','RBM24','ADCY5',
                                                                                           'AFF3','FBXO32','CORIN','DPF3','MLIP','TBX5','NonTarget'
                                                                                           ))

combine_data_select$cluster <- factor(combine_data_select$cluster,levels = c('SAN-head','SAN-tail','SAN-TZ'))
combine_data_select$perturb_target <- factor(combine_data_select$perturb_target,levels = c('L3MBTL4','CALD1','ARHGAP24','PLEKHA7','TACC2','TRDN','EPHA4',
                                                                                           'ATP2A2','CHRM2','CORIN','NPPA',"FBXO32",'SFRP1','NonTarget',
                                                                                           'TBX5','PDE4D','DPF3','AFF3','PTPRK','MLIP','CPNE5','CHST11','ADCY5'
))


dist_plot<- ggplot(combine_data_select, aes(fill=cluster, y=ratio, x=perturb_target)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c('#F8766D','#00BA38','#619CFF')) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dist_plot + RotatedAxis()

combine_data_select$cluster <- factor(combine_data_select$cluster,levels = c('SAN-tail','SAN-head','SAN-TZ'))
combine_data_select$perturb_target <- factor(combine_data_select$perturb_target,levels = c('CPNE5','PTPRK','L3MBTL4','PLEKHA7','MLIP','TRDN','DPF3','ATP2A2','ADCY5','TACC2',
                                                                                           "FBXO32",'NonTarget','CORIN','NPPA','TBX5','ARHGAP24','AFF3',
                                                                                           'CHRM2','EPHA4', 'CALD1','PDE4D','SFRP1','CHST11'
))

dist_plot<- ggplot(combine_data_select, aes(fill=cluster, y=ratio, x=perturb_target)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c('#00BA38','#F8766D','#619CFF')) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dist_plot + RotatedAxis()


combine_data_select$cluster <- factor(combine_data_select$cluster,levels = c('SAN-TZ','SAN-head','SAN-tail'))
combine_data_select$perturb_target <- factor(combine_data_select$perturb_target,levels = c('CHST11','ADCY5','PDE4D','SFRP1','AFF3','TBX5','EPHA4','CHRM2','NonTarget','NPPA',
                                                                                           "FBXO32",'CORIN','CALD1','DPF3','MLIP','ATP2A2','TACC2','ARHGAP24','TRDN','PTPRK',
                                                                                           'CPNE5','PLEKHA7','L3MBTL4'
))

dist_plot<- ggplot(combine_data_select, aes(fill=cluster, y=ratio, x=perturb_target)) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c('#619CFF','#F8766D','#00BA38')) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
dist_plot + RotatedAxis()


dotplot_ref <- DotPlot(subset(SAN_perturb_GBC, subset = perturb_target %in% names(select_perturb_target)),features = c("SHOX2"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot <- DotPlot(subset(SAN_perturb_GBC, subset = perturb_target %in% names(select_perturb_target)),features = c("SHOX2","TBX18","HCN4"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot_ref$data[order(dotplot_ref$data$avg.exp.scaled),]$id)


dotplot <- DotPlot(subset(SAN_perturb_GBC, subset = perturb_target %in% names(select_perturb_target)),features = c("SHOX2"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)


dotplot <- DotPlot(subset(SAN_perturb_GBC, subset = perturb_target %in% names(select_perturb_target)),features = c("NKX2-5"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)


dotplot <- DotPlot(subset(SAN_perturb_GBC, subset = perturb_target %in% names(select_perturb_target)),features = c("NPPA"),group.by = "perturb_target") + coord_flip() + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels =  dotplot$data[order(dotplot$data$avg.exp.scaled),]$id)

# get the correlation matrix of esc perturb
edata <- AverageExpression(SAN_perturb_GBC, group.by = 'perturb_target')$RNA

edata <- as.data.frame(edata)

edata_df <- edata[,colnames(edata) %in% names(select_perturb_target)]

edata_df <- edata_df[rownames(edata_df) %in% rownames(head_vs_tail_TZ[head_vs_tail_TZ$p_val_adj < .05 & head_vs_tail_TZ$avg_log2FC > 1.5,]),]

index <- match(c('L3MBTL4','CALD1','ARHGAP24','PLEKHA7','TACC2','TRDN','EPHA4',
                 'ATP2A2','CHRM2','CORIN','NPPA',"FBXO32",'SFRP1','NonTarget',
                 'TBX5','PDE4D','DPF3','AFF3','PTPRK','MLIP','CPNE5','CHST11','ADCY5'),colnames(edata_df))

edata_df <- edata_df[,index]

corr<-cor(edata_df)

corrplot(corr,method = "color",tl.cex = 0.4,tl.col = "black",is.corr = FALSE,type = "upper")


edata_df <- edata[,colnames(edata) %in% names(select_perturb_target)]

edata_df <- edata_df[rownames(edata_df) %in% rownames(tail_vs_head_TZ[tail_vs_head_TZ$p_val_adj < .05 & tail_vs_head_TZ$avg_log2FC > 1.5,]),]

index <- match(c('CPNE5','PTPRK','L3MBTL4','PLEKHA7','MLIP','TRDN','DPF3','ATP2A2','ADCY5','TACC2',
                 "FBXO32",'NonTarget','CORIN','NPPA','TBX5','ARHGAP24','AFF3',
                 'CHRM2','EPHA4', 'CALD1','PDE4D','SFRP1','CHST11'),colnames(edata_df))

edata_df <- edata_df[,index]

corr<-cor(edata_df)

corrplot(corr,method = "color",tl.cex = 0.4,tl.col = "black",is.corr = FALSE,type = "upper")



edata_df <- edata[,colnames(edata) %in% names(select_perturb_target)]

edata_df <- edata_df[rownames(edata_df) %in% rownames(TZ_vs_head_tail[TZ_vs_head_tail$p_val_adj < .05 & TZ_vs_head_tail$avg_log2FC > 1.5,]),]

index <- match(c('CHST11','ADCY5','PDE4D','SFRP1','AFF3','TBX5','EPHA4','CHRM2','NonTarget','NPPA',
                 "FBXO32",'CORIN','CALD1','DPF3','MLIP','ATP2A2','TACC2','ARHGAP24','TRDN','PTPRK',
                 'CPNE5','PLEKHA7','L3MBTL4'),colnames(edata_df))

edata_df <- edata_df[,index]

corr<-cor(edata_df)

corrplot(corr,method = "color",tl.cex = 0.4,tl.col = "black",is.corr = FALSE,type = "upper")

