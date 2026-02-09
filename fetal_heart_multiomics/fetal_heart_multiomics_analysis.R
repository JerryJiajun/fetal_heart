library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr',seqlevels(annotation))

# load the RNA and ATAC data
heart_3655_counts <- Read10X_h5(filename = "./3655/filtered_feature_bc_matrix.h5")
heart_3655_fragpath <- "./3655/atac_fragments.tsv.gz"

# create a Seurat object containing the RNA adata
heart_3655 <- CreateSeuratObject(
  counts = heart_3655_counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
heart_3655[["ATAC"]] <- CreateChromatinAssay(
  counts = heart_3655_counts$Peaks,
  sep = c(":","-"),
  fragments = heart_3655_fragpath,
  annotation = annotation
)

heart_3655


# load the RNA and ATAC data
heart_3675_counts <- Read10X_h5(filename = "./3675/filtered_feature_bc_matrix.h5")
heart_3675_fragpath <- "./3675/atac_fragments.tsv.gz"

# create a Seurat object containing the RNA adata
heart_3675 <- CreateSeuratObject(
  counts = heart_3675_counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
heart_3675[["ATAC"]] <- CreateChromatinAssay(
  counts = heart_3675_counts$Peaks,
  sep = c(":","-"),
  fragments = heart_3675_fragpath,
  annotation = annotation
)

heart_3675


# load the RNA and ATAC data
heart_7668_counts <- Read10X_h5(filename = "./7668/filtered_feature_bc_matrix.h5")
heart_7668_fragpath <- "./7668/atac_fragments.tsv.gz"

# create a Seurat object containing the RNA adata
heart_7668 <- CreateSeuratObject(
  counts = heart_7668_counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
heart_7668[["ATAC"]] <- CreateChromatinAssay(
  counts = heart_7668_counts$Peaks,
  sep = c(":","-"),
  fragments = heart_7668_fragpath,
  annotation = annotation
)

heart_7668


# load the RNA and ATAC data
heart_3688_counts <- Read10X_h5(filename = "./3688/filtered_feature_bc_matrix.h5")
heart_3688_fragpath <- "./3688/atac_fragments.tsv.gz"

# create a Seurat object containing the RNA adata
heart_3688 <- CreateSeuratObject(
  counts = heart_3688_counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
heart_3688[["ATAC"]] <- CreateChromatinAssay(
  counts = heart_3688_counts$Peaks,
  sep = c(":","-"),
  fragments = heart_3688_fragpath,
  annotation = annotation
)

heart_3688

fetal_heart <-merge(heart_3655, y = c(heart_3675,heart_7668,heart_3688),add.cell.ids = c("13.6w","14.6w","15.7w","16.5w"),project = "fetal_heart")

fetal_heart

# Quality control
DefaultAssay(fetal_heart) <- "ATAC"

fetal_heart <- NucleosomeSignal(fetal_heart)

fetal_heart <- TSSEnrichment(fetal_heart)

VlnPlot(
  object = fetal_heart,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


# filter out low quality cells
fetal_heart <- subset(
  x = fetal_heart,
  subset = nCount_ATAC < 100000 &
    nFeature_RNA < 4000 &
    nCount_ATAC > 1000 &
    nFeature_RNA > 200 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
fetal_heart

# The set of peaks identified using Cellranger often merges distinct peaks that are close together. 
# This can create a problem for certain analyses, particularly motif enrichment analysis and peak-to-gene linkage. 
# To identify a more accurate set of peaks, we can call peaks using MACS2 with the CallPeaks() function
# call peaks using MACS2
peaks <- CallPeaks(fetal_heart,macs2.path = "/Users/jiajunzhu/miniconda3/envs/ATAC/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(fetal_heart),
  features = peaks,
  cells = colnames(fetal_heart)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
fetal_heart[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  annotation = annotation
)


# Gene expression data processing
DefaultAssay(fetal_heart) <- "RNA"
fetal_heart <- SCTransform(fetal_heart)
fetal_heart <- RunPCA(fetal_heart,features = VariableFeatures(object = fetal_heart))
fetal_heart <- RunUMAP(fetal_heart,reduction = 'pca',dims = 1:18,reduction.name = 'umap.rna',reduction.key = 'rnaUMAP_')

# DNA accessibility data processing
DefaultAssay(fetal_heart) <- "ATAC"
fetal_heart <- FindTopFeatures(fetal_heart, min.cutoff = 'q0')
fetal_heart <- RunTFIDF(fetal_heart)
fetal_heart <- RunSVD(fetal_heart)
fetal_heart <- RunUMAP(fetal_heart,reduction = 'lsi',dims = 2:30,reduction.name = 'umap.atac',reduction.key = 'atacUMAP_')

# DNA accessibility data processing
DefaultAssay(fetal_heart) <- "peaks"
fetal_heart <- FindTopFeatures(fetal_heart, min.cutoff = 5)
fetal_heart <- RunTFIDF(fetal_heart)
fetal_heart <- RunSVD(fetal_heart)
fetal_heart <- RunUMAP(fetal_heart,reduction = 'lsi',dims = 2:30,reduction.name = 'umap.atac.peaks',reduction.key = 'atacpeaksUMAP_')


# calculate WNN graph
fetal_heart <- FindMultiModalNeighbors(fetal_heart, reduction.list = list("pca", "lsi"), dims.list = list(1:18, 2:30))
fetal_heart <- RunUMAP(fetal_heart, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
fetal_heart <- FindClusters(fetal_heart, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1 <- DimPlot(fetal_heart, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(fetal_heart, reduction = "umap.atac",label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(fetal_heart, reduction = "umap.atac.peaks",label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p4 <- DimPlot(fetal_heart, reduction = "wnn.umap",label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

# Linking peaks to genes
DefaultAssay(fetal_heart) <- "peaks"

# first compute the GC content for each peak
fetal_heart <- RegionStats(fetal_heart, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
fetal_heart <- LinkPeaks(
  object = fetal_heart,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("NRXN1","EML5","SYN3","PECAM1","CDH5","VWF","CD163","CD14","C1QA","CD3E","IL7R","ANK1","SPTA1","HBG1","WT1","BNC1","ACTA2","PDGFRB","MYH11","TNNT3","MYBPC1","RYR1","COL1A1","COL1A2","PDGFRA","MYH7","MYL2","NPPA","TBX5","SHOX2","HCN4","HCN1","PAX1","PAX9")
)


# Neural cell 12, 0,7,14,25,27,34
DotPlot(fetal_heart,features = c("NRXN1","NRXN3","PLP1"))
DotPlot(fetal_heart,features = c("NXPH1","EML5","SYN3"))
# Endothelial cell 4,15,21
DotPlot(fetal_heart,features = c("PECAM1","CDH5","VWF"))
# Endocardial 
DotPlot(fetal_heart,features = c("NFATC1","EMCN","ENG"))
# Myeloid 13
DotPlot(fetal_heart,features = c("CD163","CD14","C1QA"))
# Lymphoid 3,19,28,30
DotPlot(fetal_heart,features = c("CD3E","IL7R","CD40LG"))
# RBC (red blood cell) 26
DotPlot(fetal_heart,features = c("ANK1","SPTA1","HBG1"))
# Epicardial 32
DotPlot(fetal_heart,features = c("WT1","BNC1"))
# SMC (smooth muscle cell) 11,20
DotPlot(fetal_heart,features = c("ACTA2","PDGFRB","MYH11"))
# Epithelial 33
DotPlot(fetal_heart,features = c("TP63","PAX1","PAX9"))
# Ventricular CM 2,8,17,28,31
DotPlot(fetal_heart,features = c("IRX4","MYH7","MYL2"))
# Arial CM 1,22,24
DotPlot(fetal_heart,features = c("NPPA","NR2F2","MYL7"))
# Fibroblast 5,6,9,16,18,29
DotPlot(fetal_heart,features = c("COL1A1","COL3A1","DCN","LUM"))

DotPlot(fetal_heart,features = c("SHOX2","HCN1","CACNA1D","CACNA1G","CACNA2D2","GJC1")) + RotatedAxis()


# add annotations
fetal_heart <- RenameIdents(fetal_heart, '0' = 'Neuronal_1','7' = 'Neuronal_1','12' = 'Neuronal_2','14' = 'Neuronal_1','25' = 'Neuronal_1', '27' = 'Neuronal_1','34' = 'Neuronal_1')
fetal_heart <- RenameIdents(fetal_heart, '4' = 'Endothelial','15' = 'Endothelial','21' = 'Endothelial')
fetal_heart <- RenameIdents(fetal_heart, '13' = 'Myeloid')
fetal_heart <- RenameIdents(fetal_heart, '3' = 'Lymphoid','19' = 'Lymphoid','30' = 'Lymphoid')
fetal_heart <- RenameIdents(fetal_heart, '26' = 'RBC')
fetal_heart <- RenameIdents(fetal_heart, '32' = 'Epicardial')
fetal_heart <- RenameIdents(fetal_heart, '11' = 'SMC','20' = 'SMC')
fetal_heart <- RenameIdents(fetal_heart, '5' = 'Fibroblast','6' = 'Fibroblast','9' = 'Fibroblast','16' = 'Fibroblast','18' = 'Fibroblast','29' = 'Fibroblast')
fetal_heart <- RenameIdents(fetal_heart, '2' = 'VCM','8' = 'VCM','17' = 'VCM','23' = 'VCM','28' = 'VCM','31' = 'VCM')
fetal_heart <- RenameIdents(fetal_heart, '1' = 'ACM','24' = 'ACM')
fetal_heart <- RenameIdents(fetal_heart, '10' = 'SAN','22' = 'SAN')
fetal_heart <- RenameIdents(fetal_heart, '33' = 'Epithelial')

fetal_heart$celltype <- Idents(fetal_heart)

levels(fetal_heart) <- c('ACM','Endothelial','Epicardial','Epithelial','Fibroblast','Lymphoid','Myeloid','Neuronal_1','Neuronal_2','RBC','SAN','SMC','VCM')

fetal_heart$celltype <- factor(fetal_heart$celltype,levels =  c('ACM','Endothelial','Epicardial','Epithelial','Fibroblast','Lymphoid','Myeloid','Neuronal_1','Neuronal_2','RBC','SAN','SMC','VCM'))

p5 <- DimPlot(fetal_heart, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA") + labs(x = 'UMAP1',y= 'UMAP2')
p6 <- DimPlot(fetal_heart, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC") + labs(x = 'UMAP1',y= 'UMAP2')
p7 <- DimPlot(fetal_heart, reduction = "umap.atac.peaks", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC") + labs(x = 'UMAP1',y= 'UMAP2')
p8 <- DimPlot(fetal_heart, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN") + labs(x = 'UMAP1',y= 'UMAP2')

# find markers for every cell type compared to all remaining cells, report only the positive ones
fetal_heart.markers <- FindAllMarkers(fetal_heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fetal_heart.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top20 <- fetal_heart.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(fetal_heart.markers,file = "markers_of_each_cluster_RNA.csv")
write.csv(top20,file = "top20_markers_of_each_cluster_RNA.csv")

TF_list <- read.csv(file = "TF_list.csv",header = TRUE,row.names = NULL)

GPCR_list <- read.csv(file = "GPCR_list.csv",header = TRUE,row.names = NULL)

fetal_heart.markers_TF <- fetal_heart.markers[fetal_heart.markers$gene %in% TF_list$TF,]
write.csv(fetal_heart.markers_TF,file = "markers_annotated_RNA_TF.csv")
TF_SAN <- fetal_heart.markers_TF %>% 
  filter(cluster == 'SAN') %>% 
  arrange(desc(avg_log2FC))
DotPlot(fetal_heart,features = TF_SAN$gene[1:10], group.by = "celltype") + RotatedAxis()
VlnPlot(fetal_heart,features = TF_SAN$gene[1:4], group.by = "celltype",ncol = 4,pt.size = 0) + RotatedAxis()
FeaturePlot(fetal_heart,features = TF_SAN$gene[1:4],ncol = 4) + RotatedAxis()

fetal_heart.markers_GPCR <- fetal_heart.markers[fetal_heart.markers$gene %in% GPCR_list$GPCR,]
write.csv(fetal_heart.markers_GPCR,file = "markers_annotated_RNA_GPCR.csv")
DotPlot(fetal_heart,features = unique(fetal_heart.markers_GPCR$gene), group.by = "celltype") + RotatedAxis()
DotPlot(fetal_heart,features = c("CHRM2","HTR4","CELSR1"),group.by = "celltype") + RotatedAxis()
VlnPlot(fetal_heart,features = c("CHRM2","HTR4","CELSR1"),group.by = "celltype",pt.size = 0) + RotatedAxis()

meta <- fetal_heart@meta.data
meta$origin <- sub("_.*", "", rownames(meta))
fetal_heart@meta.data <- meta

p9 <- DimPlot(fetal_heart, reduction = "umap.rna", group.by = "celltype", split.by = "origin",label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p10 <- DimPlot(fetal_heart, reduction = "umap.atac", group.by = "celltype", split.by = "origin", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p11 <- DimPlot(fetal_heart, reduction = "umap.atac.peaks", group.by = "celltype", split.by = "origin", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p12 <- DimPlot(fetal_heart, reduction = "wnn.umap", group.by = "celltype", split.by = "origin", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

frags <- Fragments(fetal_heart[['ATAC']])
Fragments(fetal_heart[['ATAC']]) <- NULL
newpath1 <- "3655/atac_fragments.tsv.gz"
newpath2 <- "3675/atac_fragments.tsv.gz"
newpath3 <- "7668/atac_fragments.tsv.gz"
newpath4 <- "3688/atac_fragments.tsv.gz"
frags[[1]] <- UpdatePath(frags[[1]],new.path = newpath1)
frags[[2]] <- UpdatePath(frags[[2]],new.path = newpath2)
frags[[3]] <- UpdatePath(frags[[3]],new.path = newpath3)
frags[[4]] <- UpdatePath(frags[[4]],new.path = newpath4)
Fragments(fetal_heart[['ATAC']]) <- frags

# update the fragment files of peaks to ATAC
Fragments(fetal_heart[['peaks']]) <- Fragments(fetal_heart[['ATAC']])

CoveragePlot(fetal_heart, region = 'TBX5', features = 'TBX5', assay = 'ATAC', expression.assay = 'SCT', pextend.upstream = 10000, extend.downstream = 2000)

CoveragePlot(fetal_heart, region = 'HCN4', features = 'HCN4', assay = 'ATAC', expression.assay = 'SCT', pextend.upstream = 10000, extend.downstream = 2000)

CoveragePlot(fetal_heart, region = 'SHOX2', features = 'SHOX2', assay = 'ATAC', expression.assay = 'SCT', pextend.upstream = 10000, extend.downstream = 2000)

# find differential expressed genes of all celltypes
DefaultAssay(fetal_heart) <-'SCT'

DotPlot(fetal_heart,features = c("SHOX2","HCN4","HCN1","CACNA1G","CACNA2D2","GJC1"),group.by = "celltype") + RotatedAxis()+ labs(x='',y='')

DotPlot(fetal_heart,features = c("NPPA","CDH5","WT1","PAX1","PDGFRA","CD96","CD163","SYN3","NRXN1","HBG1","HCN1","MYH11","MYL2"),group.by = "celltype") + RotatedAxis()+ labs(x='',y='')

DotPlot(subset(fetal_heart, subset=celltype %in% c('ACM','VCM','SAN')),features = c("ERBB2","ERBB4","TEAD1","YAP1"),group.by = "celltype") + RotatedAxis() + labs(x='',y='')

DotPlot(fetal_heart,features = c("ERBB2","ERBB4","TEAD1","YAP1"),group.by = "celltype") + RotatedAxis() + labs(x='',y='')

DotPlot(fetal_heart,features = c("ERBB2","ERBB4","TEAD1","YAP1"),group.by = "celltype") + RotatedAxis() + labs(x='',y='')

DotPlot(subset(fetal_heart,subset = celltype %in% c('ACM','VCM','SAN')),features = c("HCN1","CACNA1D","CACNA1C","SCN5A","GJA1","GJA5","KCNJ2","KCNH2")) + RotatedAxis() + labs(x='',y='')

FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("NPPA")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CDH5")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("NFATC1")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("KDR")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("WT1")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("PAX1")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("PDGFRA")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("DCN")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CD96")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("LEF1")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CD163")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("SYN3")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("NRXN1")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("HBG1")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("MYH11")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("MYL2")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("HCN4")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("HCN1")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CACNA1G")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CACNA2D2")) + labs(x = 'UMAP1',y = 'UMAP2')
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("GJC1")) + labs(x = 'UMAP1',y = 'UMAP2')

VlnPlot(fetal_heart,features = c("NRXN1"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("SYN3"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("PECAM1"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("CDH5"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("ENG"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("CD163"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("CD96"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("LEF1"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("HBG1"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("WT1"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("MYH11"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("ACTA2"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("PDGFRA"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("DCN"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("MYL2"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("PAX1"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("NPPA"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("HCN1"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("GJC1"),pt.size = 0) + labs(x = '')
VlnPlot(fetal_heart,features = c("TBX5"),pt.size = 0) + labs(x = '')
FeaturePlot(fetal_heart,features = c("GJC1"))
FeaturePlot(fetal_heart,features = c("TBX5"))
FeaturePlot(fetal_heart,features = c("HCN4"))
FeaturePlot(fetal_heart,features = c("SHOX2"))

# General Autonomic Neuron Markers--PHOX2A,PHOX2B,TH,DBH
DotPlot(fetal_heart,features = c("TH","DBH","PHOX2A","PHOX2B"),group.by = "celltype") + RotatedAxis()

FeaturePlot(fetal_heart,features = c("TH"),max.cutoff = 2)
VlnPlot(fetal_heart,features = c("TH"),pt.size = 0)
FeaturePlot(fetal_heart,features = c("DBH"),max.cutoff = 3)
VlnPlot(fetal_heart,features = c("DBH"),pt.size = 0)
FeaturePlot(fetal_heart,features = c("PHOX2B"),max.cutoff = 2)
VlnPlot(fetal_heart,features = c("PHOX2B"),pt.size = 0)

# Sypathetic Neurons--TH, DBH, PRPH, NTPK1
# Parasympathetic Neurons--CHAT, PHOX2B, RET
DotPlot(fetal_heart,features = c("TH","DBH","PRPH","CHAT","PHOX2B","RET"),group.by = "celltype") + RotatedAxis()


fetal_heart.markers <- FindAllMarkers(fetal_heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fetal_heart.markers <- fetal_heart.markers %>%
  filter(p_val_adj < 0.05) %>%
  filter(!grepl("^MT-",gene))

fetal_heart.markers_top30 <- fetal_heart.markers %>% 
  group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(fetal_heart.markers_top30,file = "Top_30_markers_annotated_RNA.csv")

fetal_heart.markers_top10 <- fetal_heart.markers %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

fetal_heart.markers_unique_top10 <- fetal_heart.markers %>% 
  filter(isUnique(gene)) %>%
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

write.csv(fetal_heart.markers,file = "markers_annotated_RNA.csv")

# add the gene activity matrix to the seurat object as a new assay and normalize it
DefaultAssay(fetal_heart) <- 'peaks'
gene.activities <- GeneActivity(fetal_heart)

fetal_heart[['activity']] <- CreateAssayObject(counts = gene.activities)
fetal_heart <- NormalizeData(
  object = fetal_heart,
  assay = 'activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(fetal_heart$nCount_RNA)
)

DefaultAssay(fetal_heart) <- 'activity'

FeaturePlot(
  object = fetal_heart,
  features = c("SHOX2","HCN4","HCN1"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)

DotPlot(fetal_heart,features = c("SHOX2","HCN4","HCN1","CACNA1D","CACNA1G","CACNA2D2"),group.by = "celltype") + RotatedAxis()


# heatmap to demonstrate the top expressed genes in each cell type
DefaultAssay(fetal_heart) <-'SCT'
DoHeatmap(subset(fetal_heart,downsample = 100),features = fetal_heart.markers_top20$gene,group.by = "celltype",size = 2.5,angle = 90) + NoLegend()

DoHeatmap(subset(fetal_heart,downsample = 100),features = fetal_heart.markers_top10$gene,group.by = "celltype",size = 2.5) + NoLegend()

DoHeatmap(subset(fetal_heart,downsample = 100),features = fetal_heart.markers_unique_top10$gene,group.by = "celltype",size = 2.5) + NoLegend()

# find differential expressed genes between SAN and Ventricle, SAN and Epicardial
SAN_vs_Ventricle <- FindMarkers(fetal_heart,ident.1 = 'SAN',ident.2 = 'Ventricle')
SAN_vs_Ventricle <- SAN_vs_Ventricle %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC),desc(pct.1)) 
SAN_vs_Ventricle <- SAN_vs_Ventricle[-which(grepl("^MT-",rownames(SAN_vs_Ventricle))),]
DotPlot(fetal_heart,features = rownames(SAN_vs_Ventricle[1:20,]),group.by = "celltype") + RotatedAxis()
DotPlot(fetal_heart,features = rownames(SAN_vs_Ventricle[441:460,]),group.by = "celltype") + RotatedAxis()
write.csv(SAN_vs_Ventricle,file = "RNA_SAN_vs_Ventricle.csv")

SAN_vs_Atrial <- FindMarkers(fetal_heart,ident.1 = 'SAN',ident.2 = 'Atrial')
SAN_vs_Atrial <- SAN_vs_Atrial %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC),desc(pct.1))
DotPlot(fetal_heart,features = rownames(SAN_vs_Atrial[1:20,]),group.by = "celltype") + RotatedAxis()
DotPlot(fetal_heart,features = rownames(SAN_vs_Atrial[922:941,]),group.by = "celltype") + RotatedAxis()
write.csv(SAN_vs_Atrial,file = "RNA_SAN_vs_Atrial.csv")

SAN_vs_Epicardial <- FindMarkers(fetal_heart,ident.1 = 'SAN',ident.2 = 'Epicardial')
SAN_vs_Epicardial <- SAN_vs_Epicardial %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC),desc(pct.1))
DotPlot(fetal_heart,features = rownames(SAN_vs_Epicardial[1:20,]),group.by = "celltype") + RotatedAxis()
DotPlot(fetal_heart,features = rownames(SAN_vs_Epicardial[1048:1067,]),group.by = "celltype") + RotatedAxis()
write.csv(SAN_vs_Epicardial,file = "RNA_SAN_vs_Epicardial.csv")

DefaultAssay(fetal_heart) <-'ATAC'
fetal_heart.markers.ATAC <- FindAllMarkers(fetal_heart, test.use = 'LR',latent.vars = 'nCount_ATAC')
fetal_heart.markers.ATAC %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(fetal_heart.markers.ATAC,file = "markers_annotated_ATAC.csv")
fetal_heart.markers.ATAC.sig <- fetal_heart.markers.ATAC %>%
  filter(p_val_adj < 0.001)
write.csv(fetal_heart.markers.ATAC.sig,file = "markers_sig_annotated_ATAC.csv")

CoveragePlot(fetal_heart, region = rownames(fetal_heart.markers.ATAC)[1], assay = 'ATAC', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 20000)

CoveragePlot(fetal_heart, region = rownames(fetal_heart.markers.ATAC.sig)[13], assay = 'ATAC', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 20000)

CoveragePlot(fetal_heart, region = rownames(fetal_heart.markers.ATAC.sig)[14], assay = 'ATAC', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 15000)

CoveragePlot(fetal_heart, region = rownames(fetal_heart.markers.ATAC.sig)[15], assay = 'ATAC', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 30000)

DefaultAssay(fetal_heart) <-'peaks'

fetal_heart.markers.ATAC.peaks <- FindAllMarkers(fetal_heart, test.use = 'LR',latent.vars = 'nCount_peaks')
fetal_heart.markers.ATAC.peaks %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(fetal_heart.markers.ATAC.peaks,file = "markers_annotated_ATAC_peaks.csv")
fetal_heart.markers.ATAC.peaks.sig <- fetal_heart.markers.ATAC.peaks %>%
  filter(p_val_adj < 0.001)
write.csv(fetal_heart.markers.ATAC.peaks.sig,file = "markers_sig_annotated_ATAC_peaks.csv")

CoveragePlot(fetal_heart, region = c("NXPH1","EML5","SYN3","PECAM1","CDH5","VWF","CD163","CD14","C1QA","CD3E","IL7R","ANK1","SPTA1","HBG1","WT1","BNC1","ACTA2","PDGFRB","MYH11","TNNT3","MYBPC1","RYR1","COL1A1","COL1A2","PDGFRA","MYH7","MYL2","NPPA","TBX5","SHOX2","HCN4","HCN1","PAX1","PAX9"),annotation = FALSE,peaks = FALSE)

CoveragePlot(fetal_heart, region = c("SHOX2","HCN4","CD163","CD3E","BNC1","IRX4"),annotation = FALSE,peaks = FALSE, expression.assay = "SCT",pextend.upstream = 1000, extend.downstream = 1000, ncol = 6)

CoveragePlot(fetal_heart, region = c("SHOX2","HCN4","CACNA1G","TBX5","TBX3","HCN1"),annotation = FALSE,peaks = FALSE, expression.assay = "SCT",pextend.upstream = 1000, extend.downstream = 1000, ncol = 6)

CoveragePlot(fetal_heart, region = "MYL2",features = "MYL2",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "MYH11",features = "MYH11",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "HCN1",features = "HCN1",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "HCN4",features = "HCN4",annotation = TRUE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "HBG1",features = "HBG1",annotation = TRUE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "SYN3",features = "SYN3",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "NRXN1",features = "NRXN1",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "CD163",features = "CD163",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "CD96",features = "CD96",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "PDGFRA",features = "PDGFRA",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "PAX1",features = "PAX1",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "WT1",features = "WT1",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "CDH5",features = "CDH5",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "NPPA",features = "NPPA",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "TBX5",features = "TBX5",annotation = TRUE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)

CoveragePlot(fetal_heart, region = "TH",features = "TH",annotation = TRUE,peaks = TRUE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "PHOX2B",features = "PHOX2B",annotation = TRUE,peaks = TRUE, pextend.upstream = 1000, extend.downstream = 1000)

CoveragePlot(fetal_heart, region = "HCN4",features = "HCN4",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "TBX5",features = "TBX5",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "SHOX2",features = "SHOX2",annotation = TRUE,peaks = FALSE, pextend.upstream = 5000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = "CACNA1D",features = "CACNA1D",annotation = TRUE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "CACNA1G",features = "CACNA1G",annotation = TRUE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "CACNA2D2",features = "CACNA2D2",annotation = TRUE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "GJC1",features = "GJC1",annotation = TRUE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)

CoveragePlot(subset(fetal_heart,subset=celltype %in% c('ACM','VCM','SAN')), region = "YAP1",features = "YAP1",annotation = TRUE,peaks = TRUE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(subset(fetal_heart,subset=celltype %in% c('ACM','VCM','SAN')), region = "chr11-102120000-102160000",annotation = TRUE,peaks = TRUE)

CoveragePlot(subset(fetal_heart,subset=celltype %in% c('ACM','VCM','SAN')), region = "TEAD1",features = "YAP1",annotation = TRUE,peaks = TRUE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(subset(fetal_heart,subset=celltype %in% c('ACM','VCM','SAN')), region = "chr11-12670000-12700000",annotation = TRUE,peaks = TRUE)

# NPPA region
CoveragePlot(fetal_heart, region = "chr1-11847000-11855000",annotation = TRUE,peaks = FALSE)
# CDH5 region
CoveragePlot(fetal_heart, region = "chr16-66366000-66374000",annotation = TRUE,peaks = FALSE)
# WT1 region
CoveragePlot(fetal_heart, region = "chr11-32432000-32440000",annotation = TRUE,peaks = FALSE)
# PAX1 region
CoveragePlot(fetal_heart, region = "chr20-21705000-21713000",annotation = TRUE,peaks = FALSE)
# PDGFRA region
CoveragePlot(fetal_heart, region = "chr4-54228000-54232000",annotation = TRUE,peaks = FALSE)
# CD163 region
CoveragePlot(fetal_heart, region = "chr12-7500000-7508000",annotation = TRUE,peaks = FALSE)
# CD96
CoveragePlot(fetal_heart, region = "chr3-111290000-111300000",annotation = TRUE,peaks = FALSE)
# SYN3 region
CoveragePlot(fetal_heart, region = "chr22-33050000-33060000",annotation = TRUE,peaks = FALSE)
# NRXN1
CoveragePlot(fetal_heart, region = "chr2-51220000-51230000",annotation = TRUE,peaks = FALSE)
# MYH11 region
CoveragePlot(fetal_heart, region = "chr16-15850000-15860000",annotation = TRUE,peaks = FALSE)
# HCN1 region
CoveragePlot(fetal_heart, region = "chr5-45690000-45697000",annotation = TRUE,peaks = FALSE)
# HCN4 region
CoveragePlot(fetal_heart, region = "chr15-73340000-73380000",annotation = TRUE,peaks = FALSE)
# TBX5 region
CoveragePlot(fetal_heart, region = "chr12-114390000-114440000",annotation = TRUE,peaks = FALSE)
# GJC1 region
CoveragePlot(fetal_heart, region = "chr17-44820000-44860000",annotation = TRUE,peaks = FALSE)

# TH region
CoveragePlot(fetal_heart, region = "chr11-2164000-2174000",annotation = TRUE,peaks = FALSE)
# PHOX2B region
CoveragePlot(fetal_heart, region = "chr4-41744000-41750000",annotation = TRUE,peaks = FALSE)

CoveragePlot(fetal_heart, region = c("MYL2","MYH11","MYBPC1","HCN1","HBG1","STMN2","CD163","LEF1","PDGFRA","PAX1","WT1","CDH5","NPPA"),annotation = FALSE,peaks = FALSE)

saveRDS(fetal_heart, file = "./fetal_heart_multiomics_annotated.rds")

# find differential ATAC peaks between SAN and Ventricle, SAN and Epicardial
SAN_vs_Ventricle_ATAC <- FindMarkers(fetal_heart,ident.1 = 'SAN',ident.2 = 'Ventricle',test.use = 'LR',latent.vars = 'nCount_peaks')
SAN_vs_Ventricle_ATAC <- SAN_vs_Ventricle_ATAC %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC),desc(pct.1)) 
SAN_vs_Ventricle_ATAC_plus_genes <- ClosestFeature(object = fetal_heart,regions = rownames(SAN_vs_Ventricle_ATAC))
SAN_vs_Ventricle_ATAC <- cbind(SAN_vs_Ventricle_ATAC,SAN_vs_Ventricle_ATAC_plus_genes)
write.csv(SAN_vs_Ventricle_ATAC,file = "ATAC_SAN_vs_Ventricle.csv")
CoveragePlot(fetal_heart, region = SAN_vs_Ventricle_ATAC$closest_region[1], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = SAN_vs_Ventricle_ATAC$closest_region[8], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)

SAN_vs_Epicardial_ATAC <- FindMarkers(fetal_heart,ident.1 = 'SAN',ident.2 = 'Epicardial',test.use = 'LR',latent.vars = 'nCount_peaks')
SAN_vs_Epicardial_ATAC <- SAN_vs_Epicardial_ATAC %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC),desc(pct.1)) 
SAN_vs_Epicardial_ATAC_plus_genes <- ClosestFeature(object = fetal_heart,regions = rownames(SAN_vs_Epicardial_ATAC))
SAN_vs_Epicardial_ATAC <- cbind(SAN_vs_Epicardial_ATAC,SAN_vs_Epicardial_ATAC_plus_genes)
write.csv(SAN_vs_Epicardial_ATAC,file = "ATAC_SAN_vs_Epicardial.csv")
CoveragePlot(fetal_heart, region = SAN_vs_Epicardial_ATAC$closest_region[1], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = SAN_vs_Epicardial_ATAC$closest_region[10], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)


SAN_vs_Atrial_ATAC <- FindMarkers(fetal_heart,ident.1 = 'SAN',ident.2 = 'Atrial',test.use = 'LR',latent.vars = 'nCount_peaks')
SAN_vs_Atrial_ATAC <- SAN_vs_Atrial_ATAC %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC),desc(pct.1)) 
SAN_vs_Atrial_ATAC_plus_genes <- ClosestFeature(object = fetal_heart,regions = rownames(SAN_vs_Atrial_ATAC))
SAN_vs_Atrial_ATAC <- cbind(SAN_vs_Atrial_ATAC,SAN_vs_Atrial_ATAC_plus_genes)
write.csv(SAN_vs_Atrial_ATAC,file = "ATAC_SAN_vs_Atrial.csv")
CoveragePlot(fetal_heart, region = SAN_vs_Atrial_ATAC$closest_region[1], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = SAN_vs_Atrial_ATAC$closest_region[8], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)

# make the annotation name consistent
fetal_heart <- RenameIdents(fetal_heart, 'Epithelial' = 'Epithelial',
                            'SAN' = 'SAN',
                            'Atrial' = 'ACM',
                            'Ventricle' = 'VCM',
                            'Fibroblast' = 'Fibroblast',
                            'SMC' = 'SMC',
                            'VSM' = 'VSM',
                            'Epicardial' = 'Epicardial',
                            'RBC' = 'RBC',
                            'Lymphoid' = 'Lymphoid',
                            'Myeloid' = 'Myeloid',
                            'Endothelial' = 'Endothelial',
                            'Neuronal' = 'Neuronal')
fetal_heart$celltype <- Idents(fetal_heart)
saveRDS(fetal_heart, file = "./fetal_heart_multiomics_annotated.rds")

library(dplyr)
fetal_heart.markers.ATAC.peaks.sig %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
SAN_top_peaks <- fetal_heart.markers.ATAC.peaks.sig[fetal_heart.markers.ATAC.peaks.sig$cluster == "SAN",]
ClosestFeature(object = fetal_heart,regions = SAN_top_peaks$gene[1:10])
SAN_top_peaks_features <- ClosestFeature(object = fetal_heart,regions = SAN_top_peaks$gene)
write.csv(SAN_top_peaks_features,file = "SAN_top_peaks_closestFeatures.csv")

CoveragePlot(fetal_heart, region = SAN_top_peaks_features$closest_region[7], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = SAN_top_peaks_features$closest_region[3], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)
CoveragePlot(fetal_heart, region = SAN_top_peaks_features$closest_region[5], assay = 'peaks', expression.assay = 'SCT', pextend.upstream = 20000, extend.downstream = 5000)

# find out the differentially expressed genes between different stage in SAN cluster
meta <- fetal_heart@meta.data
meta$origin <- sub("_.*", "", rownames(meta))
fetal_heart@meta.data <- meta

celltype <- fetal_heart$celltype
origin <- fetal_heart$origin

celltype_table <- table(fetal_heart$origin,fetal_heart$celltype)
celltype_prop_table <- prop.table(celltype_table,margin = 1)

sample <- rep(rownames(celltype_prop_table),each = dim(celltype_prop_table)[2])
cell <- rep(colnames(celltype_prop_table),times = dim(celltype_prop_table)[1])
ratio <- c(t(celltype_prop_table))

combine_data <- data.frame(sample,cell,ratio)

# Stacked + percent
ggplot(combine_data, aes(fill=cell, y=ratio, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()

Idents(fetal_heart) <- paste(origin,celltype,sep = '_')
levels(fetal_heart)

SAN_16.5w_vs_13.6w <- FindMarkers(fetal_heart, ident.1 = '16.5w_SAN', ident.2 = '13.6w_SAN')
SAN_16.5w_vs_13.6w <- SAN_16.5w_vs_13.6w[SAN_16.5w_vs_13.6w$p_val_adj < 0.05,]

SAN_16.5w_vs_14.6w <- FindMarkers(fetal_heart, ident.1 = '16.5w_SAN', ident.2 = '14.6w_SAN')
SAN_16.5w_vs_14.6w <- SAN_16.5w_vs_14.6w[SAN_16.5w_vs_14.6w$p_val_adj < 0.05,]

SAN_16.5w_vs_15.7w <- FindMarkers(fetal_heart, ident.1 = '16.5w_SAN', ident.2 = '15.7w_SAN')
SAN_16.5w_vs_15.7w <- SAN_16.5w_vs_15.7w[SAN_16.5w_vs_15.7w$p_val_adj < 0.05,]

SAN_14.6w_vs_13.6w <- FindMarkers(fetal_heart, ident.1 = '14.6w_SAN', ident.2 = '13.6w_SAN')
SAN_14.6w_vs_13.6w <- SAN_14.6w_vs_13.6w[SAN_14.6w_vs_13.6w$p_val_adj < 0.05,]

SAN_15.7w_vs_14.6w <- FindMarkers(fetal_heart, ident.1 = '15.7w_SAN', ident.2 = '14.6w_SAN')
SAN_15.7w_vs_14.6w <- SAN_15.7w_vs_14.6w[SAN_15.7w_vs_14.6w$p_val_adj < 0.05,]


# subset by sample origin timepoint
fetal_heart_13.6w <- subset(fetal_heart,subset = origin == "13.6w")
saveRDS(fetal_heart_13.6w, file = "fetal_heart_13.6w.rds")

fetal_heart_14.6w <- subset(fetal_heart,subset = origin == "14.6w")
saveRDS(fetal_heart_14.6w, file = "fetal_heart_14.6w.rds")

fetal_heart_15.7w <- subset(fetal_heart,subset = origin == "15.7w")
saveRDS(fetal_heart_15.7w, file = "fetal_heart_15.7w.rds")

fetal_heart_16.5w <- subset(fetal_heart,subset = origin == "16.5w")
saveRDS(fetal_heart_16.5w, file = "fetal_heart_16.5w.rds")

# subclustering of SAN
fetal_SAN <- subset(fetal_heart,subset = celltype == 'SAN')
DefaultAssay(fetal_SAN) <- "RNA"

fetal_SAN <- FindNeighbors(fetal_SAN, dims = 1:20)
fetal_SAN <- FindClusters(fetal_SAN, graph.name = 'wsnn',resolution = 1.5)
# resolution can adjust from 0.2 to 1.5

head(Idents(fetal_SAN), 5)

#Run non-linear dimensional reduction (UMAP/tSNE)
fetal_SAN <- RunUMAP(fetal_SAN, dims = 1:20)

UMAPPlot(fetal_SAN,label=TRUE)

saveRDS(fetal_SAN, file = "./fetal_heart_SAN_Sub.rds")

fetal_SAN <- RenameIdents(fetal_SAN, '0' = 'SAN-head','3' = 'SAN-head','1' = 'SAN-tail','6' = 'SAN-tail','7' = 'SAN-tail','4' = 'SAN-TZ','2' = 'SAN-TZ','5' = 'SAN-TZ')

fetal_SAN$celltype <- Idents(fetal_SAN)

UMAPPlot(fetal_SAN,label=TRUE)

DotPlot(fetal_SAN,features = c("SHOX2","TBX18","NKX2-5","NPPA"))
VlnPlot(fetal_SAN,features = "SHOX2")
FeaturePlot(fetal_SAN,features = "SHOX2",reduction = "umap")
VlnPlot(fetal_SAN,features = "NKX2-5")
FeaturePlot(fetal_SAN,features = "NKX2-5",reduction = "umap")
VlnPlot(fetal_SAN,features = "NPPA",pt.size = 0)
FeaturePlot(fetal_SAN,features = "NPPA",reduction = "umap")

DotPlot(fetal_SAN,features = c('MEF2A','TEAD1','MGA')) + RotatedAxis()
DotPlot(fetal_SAN,features = c("MEF2A","TEAD1","MGA"))
VlnPlot(fetal_SAN,features = "MEF2A",pt.size = 0)
FeaturePlot(fetal_SAN,features = "MEF2A",reduction = "umap")
VlnPlot(fetal_SAN,features = "TEAD1",pt.size = 0)
FeaturePlot(fetal_SAN,features = "TEAD1",reduction = "umap")
VlnPlot(fetal_SAN,features = "MGA",pt.size = 0)
FeaturePlot(fetal_SAN,features = "MGA",reduction = "umap")

DotPlot(fetal_SAN,features = c('L3MBTL4','CALD1','ARHGAP24','PLEKHA7','CPNE5','PTPRK','CHST11','ADCY5','PDE4D','SFRP1')) + RotatedAxis()

DotPlot(fetal_SAN,features = c('ARHGAP24','PTPRK','CHST11','SFRP1')) + RotatedAxis()
DotPlot(fetal_SAN,features = c('ARHGAP24','PTPRK','CHST11')) + RotatedAxis()
VlnPlot(fetal_SAN,features = "ARHGAP24",pt.size = 0)
FeaturePlot(fetal_SAN,features = "ARHGAP24",reduction = "umap")
VlnPlot(fetal_SAN,features = "PTPRK",pt.size = 0)
FeaturePlot(fetal_SAN,features = "PTPRK",reduction = "umap")
VlnPlot(fetal_SAN,features = "CHST11",pt.size = 0)
FeaturePlot(fetal_SAN,features = "CHST11",reduction = "umap")
VlnPlot(fetal_SAN,features = "SFRP1",pt.size = 0)
FeaturePlot(fetal_SAN,features = "SFRP1",reduction = "umap")

