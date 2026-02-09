library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

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


# Neural cell
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("NRXN1","NRXN3","PLP1"),ncol = 3)
DotPlot(fetal_heart,features = c("NRXN1","NRXN3","PLP1"))
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("NXPH1","EML5","SYN3"),ncol = 3)
DotPlot(fetal_heart,features = c("NXPH1","EML5","SYN3"))
# Endothelial cell
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("PECAM1","CDH5","VWF"),ncol = 3)
DotPlot(fetal_heart,features = c("PECAM1","CDH5","VWF"))
# Myeloid
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CD163","CD14","C1QA"),ncol = 3)
DotPlot(fetal_heart,features = c("CD163","CD14","C1QA"))
# Lymphoid
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CD3E","IL7R","CD40LG"),ncol = 3)
DotPlot(fetal_heart,features = c("CD3E","IL7R","CD40LG"))
# RBC (red blood cell)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("ANK1","SPTA1","HBG1"),ncol = 3)
DotPlot(fetal_heart,features = c("ANK1","SPTA1","HBG1"))
# Epicardial
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("WT1","BNC1"),ncol = 3)
DotPlot(fetal_heart,features = c("WT1","BNC1"))
# Ventricular CM
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("IRX4","MYH7","MYL2"),ncol = 3)
# VSM (vascular smooth muscle cell)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("ACTA2","PDGFRB","MYH11"),ncol = 3)
DotPlot(fetal_heart,features = c("ACTA2","PDGFRB","MYH11"))
# Epithelial
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("TP63","PAX1","PAX9"),ncol = 3)
DotPlot(fetal_heart,features = c("TP63","PAX1","PAX9"))
# SMC (skeletal muscle cell)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("TNNT3","MYBPC1","RYR1"),ncol = 3)
DotPlot(fetal_heart,features = c("TNNT3","MYBPC1","RYR1"))
# Ventricular CM
DotPlot(fetal_heart,features = c("IRX4","MYH7","MYL2"))
# Arial CM
DotPlot(fetal_heart,features = c("NPPA","NR2F2","MYL7"))
# Fibroblast
DotPlot(fetal_heart,features = c("COL1A1","COL1A2","PDGFRA"))

FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("SHOX2","HCN4","HCN1","CACNA1D","CACNA1G","CACNA2D2","TBX5"),ncol = 3)
DotPlot(fetal_heart,features = c("SHOX2","HCN1","CACNA1D","CACNA1G","CACNA2D2","GJC1")) + RotatedAxis()


# find markers for every cluster compared to all remaining cells, report only the positive ones
fetal_heart.markers <- FindAllMarkers(fetal_heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
fetal_heart.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(fetal_heart.markers,file = "fetal_heart_markers_multiomics.csv")

# add annotations
fetal_heart <- RenameIdents(fetal_heart, '0' = 'Neuronal','7' = 'Neuronal','12' = 'Neuronal','14' = 'Neuronal','25' = 'Neuronal', '27' = 'Neuronal','34' = 'Neuronal')
fetal_heart <- RenameIdents(fetal_heart, '4' = 'Endothelial','15' = 'Endothelial','21' = 'Endothelial')
fetal_heart <- RenameIdents(fetal_heart, '13' = 'Myeloid')
fetal_heart <- RenameIdents(fetal_heart, '3' = 'Lymphoid','19' = 'Lymphoid','28' = 'Lymphoid','30' = 'Lymphoid')
fetal_heart <- RenameIdents(fetal_heart, '26' = 'RBC')
fetal_heart <- RenameIdents(fetal_heart, '32' = 'Epicardial')
fetal_heart <- RenameIdents(fetal_heart, '11' = 'VSM','18' = 'VSM','20' = 'VSM')
fetal_heart <- RenameIdents(fetal_heart, '29' = 'SMC')
fetal_heart <- RenameIdents(fetal_heart, '5' = 'Fibroblast','6' = 'Fibroblast','9' = 'Fibroblast','16' = 'Fibroblast')
fetal_heart <- RenameIdents(fetal_heart, '2' = 'Ventricle','8' = 'Ventricle','17' = 'Ventricle','23' = 'Ventricle','31' = 'Ventricle')
fetal_heart <- RenameIdents(fetal_heart, '1' = 'Atrial','22' = 'Atrial','24' = 'Atrial')
fetal_heart <- RenameIdents(fetal_heart, '10' = 'SAN')
fetal_heart <- RenameIdents(fetal_heart, '33' = 'Epithelial')


fetal_heart$celltype <- Idents(fetal_heart)

p5 <- DimPlot(fetal_heart, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p6 <- DimPlot(fetal_heart, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p7 <- DimPlot(fetal_heart, reduction = "umap.atac.peaks", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p8 <- DimPlot(fetal_heart, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

# find markers for every cell type compared to all remaining cells, report only the positive ones
library(dplyr)
fetal_heart.markers <- FindAllMarkers(fetal_heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fetal_heart.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

top20 <- fetal_heart.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(fetal_heart.markers,file = "markers_annotated_RNA.csv")

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
DefaultAssay(fetal_heart) <-'RNA'

levels(fetal_heart) <- c('Ventricle','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','Atrial')

levels(fetal_heart) <- c('Atrial','Endothelial','Epicardial','Epithelial','Fibroblast','Lymphoid','Myeloid','Neuronal','RBC','SAN','SMC','VSM','Ventricle')

DotPlot(fetal_heart,features = c("SHOX2","HCN4","HCN1","CACNA1D","CACNA2D2","TBX5","GJC1"),group.by = "celltype") + RotatedAxis()

DotPlot(fetal_heart,features = c("NRXN1","EML5","SYN3","PECAM1","CDH5","VWF","CD163","CD14","C1QA","CD3E","IL7R","ANK1","SPTA1","HBG1","WT1","BNC1","ACTA2","PDGFRB","MYH11","TNNT3","MYBPC1","RYR1","COL1A1","COL1A2","PDGFRA","MYH7","MYL2","NPPA","TBX5","SHOX2","HCN4","HCN1","PAX1","PAX9"),group.by = "celltype") + RotatedAxis()

DotPlot(fetal_heart,features = c("EML5","PECAM1","CD163","CD3E","HBG1","WT1","MYH11","MYBPC1","PDGFRA","MYL2","PDE3A","HCN4","PAX1"),group.by = "celltype") + RotatedAxis()

DotPlot(fetal_heart,features = c("EML5","PECAM1","CD163","CD3E","HBG1","WT1","MYH11","MYBPC1","PDGFRA","MYL2","NPPA","HCN4","PAX1"),group.by = "celltype") + RotatedAxis()

dotplot <- DotPlot(fetal_heart,features = c("MYL2","MYH11","MYBPC1","HCN1","HBG1","EML5","CD163","CD96","PDGFRA","PAX1","WT1","CDH5","NPPA"),group.by = "celltype") + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels = c('Ventricle','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','Atrial'))

dotplot <- DotPlot(fetal_heart,features = c("MYL2","MYH11","MYBPC1","HCN4","HBG1","EML5","CD163","CD96","PDGFRA","PAX1","WT1","CDH5","NPPA"),group.by = "celltype") + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels = c('Ventricle','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','Atrial'))

dotplot <- DotPlot(fetal_heart,features = c("MYL2","MYH11","MYBPC1","HCN4","HBG1","STMN2","CD163","CD96","PDGFRA","PAX1","WT1","CDH5","NPPA"),group.by = "celltype") + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels = c('Ventricle','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','Atrial'))

FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("NPPA"),max.cutoff = 20)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("MYH6"))
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("PECAM1"),max.cutoff = 10)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CDH5"),max.cutoff = 4)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("WT1"),max.cutoff = 6)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("PAX1"),max.cutoff = 3)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("PDGFRA"),max.cutoff = 10)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CD96"),max.cutoff = 10)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("LEF1"),max.cutoff = 15)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CD163"),max.cutoff = 8)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("EML5"),max.cutoff = 20)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("STMN2"),max.cutoff = 6)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("HBG1"),max.cutoff = 20)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("MYBPC1"),max.cutoff = 10)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("MYH11"),max.cutoff = 15)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("MYL2"),max.cutoff = 8)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("HCN4"),max.cutoff = 4)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("HCN1"),max.cutoff = 8)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("TBX5"),max.cutoff = 15)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CACNA1D"),max.cutoff = 12)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("CACNA2D2"),max.cutoff = 3)
FeaturePlot(fetal_heart,reduction = "umap.rna",features = c("GJC1"))

VlnPlot(fetal_heart,features = c("STMN2"),pt.size = 0)
VlnPlot(fetal_heart,features = c("EML5"),pt.size = 0)
VlnPlot(fetal_heart,features = c("PECAM1"),pt.size = 0)
VlnPlot(fetal_heart,features = c("CDH5"),pt.size = 0)
VlnPlot(fetal_heart,features = c("CD163"),pt.size = 0)
VlnPlot(fetal_heart,features = c("CD96"),pt.size = 0)
VlnPlot(fetal_heart,features = c("LEF1"),pt.size = 0)
VlnPlot(fetal_heart,features = c("HBG1"),pt.size = 0)
VlnPlot(fetal_heart,features = c("WT1"),pt.size = 0)
VlnPlot(fetal_heart,features = c("MYH11"),pt.size = 0)
VlnPlot(fetal_heart,features = c("ACTA2"),pt.size = 0)
VlnPlot(fetal_heart,features = c("MYBPC1"),pt.size = 0)
VlnPlot(fetal_heart,features = c("PDGFRA"),pt.size = 0)
VlnPlot(fetal_heart,features = c("MYL2"),pt.size = 0)
VlnPlot(fetal_heart,features = c("PAX1"),pt.size = 0)
VlnPlot(fetal_heart,features = c("NPPA"),pt.size = 0)
VlnPlot(fetal_heart,features = c("PDE3A"),pt.size = 0)
VlnPlot(fetal_heart,features = c("HCN1"),pt.size = 0)
VlnPlot(fetal_heart,features = c("CACNA1D"),pt.size = 0)
VlnPlot(fetal_heart,features = c("GJC1"),pt.size = 0)
VlnPlot(fetal_heart,features = c("TBX5"),pt.size = 0)

# General Autonomic Neuron Markers--PHOX2A,PHOX2B,TH,DBH
dotplot <- DotPlot(fetal_heart,features = c("TH","DBH","PHOX2A","PHOX2B"),group.by = "celltype") + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels = c('Ventricle','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','Atrial'))

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

dotplot <- DotPlot(fetal_heart,features = c("SHOX2","HCN4","HCN1","CACNA1D","CACNA1G","CACNA2D2","TBX5"),group.by = "celltype") + RotatedAxis()
dotplot$data$id <- factor(dotplot$data$id,levels = c('Ventricle','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','Atrial'))


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

levels(fetal_heart) <- c('Atrial','Endothelial','Epicardial','Epithelial','Fibroblast','Lymphoid','Myeloid','Neuronal','RBC','SAN','SMC','VSM','Ventricle')

fetal_heart.markers.ATAC.peaks <- FindAllMarkers(fetal_heart, test.use = 'LR',latent.vars = 'nCount_peaks')
fetal_heart.markers.ATAC.peaks %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
write.csv(fetal_heart.markers.ATAC.peaks,file = "markers_annotated_ATAC_peaks.csv")
fetal_heart.markers.ATAC.peaks.sig <- fetal_heart.markers.ATAC.peaks %>%
  filter(p_val_adj < 0.001)
write.csv(fetal_heart.markers.ATAC.peaks.sig,file = "markers_sig_annotated_ATAC_peaks.csv")

CoveragePlot(fetal_heart, region = c("NXPH1","EML5","SYN3","PECAM1","CDH5","VWF","CD163","CD14","C1QA","CD3E","IL7R","ANK1","SPTA1","HBG1","WT1","BNC1","ACTA2","PDGFRB","MYH11","TNNT3","MYBPC1","RYR1","COL1A1","COL1A2","PDGFRA","MYH7","MYL2","NPPA","TBX5","SHOX2","HCN4","HCN1","PAX1","PAX9"),annotation = FALSE,peaks = FALSE)

CoveragePlot(fetal_heart, region = c("SHOX2","HCN4","CD163","CD3E","BNC1","IRX4"),annotation = FALSE,peaks = FALSE, expression.assay = "SCT",pextend.upstream = 1000, extend.downstream = 1000, ncol = 6)

CoveragePlot(fetal_heart, region = c("SHOX2","HCN4","CACNA1G","TBX5","TBX3","HCN1"),annotation = FALSE,peaks = FALSE, expression.assay = "SCT",pextend.upstream = 1000, extend.downstream = 1000, ncol = 6)

CoveragePlot(fetal_heart, region = "MYL2",features = "MYL2",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "MYH11",features = "MYH11",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "MYBPC1",features = "MYBPC1",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "HCN1",features = "HCN1",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "HCN4",features = "HCN4",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "HBG1",features = "HBG1",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "EML5",features = "EML5",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "STMN2",features = "STMN2",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "CD163",features = "CD163",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "CD96",features = "CD96",annotation = FALSE,peaks = FALSE, pextend.upstream = 2000, extend.downstream = 200)
CoveragePlot(fetal_heart, region = "LEF1",features = "LEF1",annotation = FALSE,peaks = FALSE, pextend.upstream = 2000, extend.downstream = 200)
CoveragePlot(fetal_heart, region = "PDGFRA",features = "PDGFRA",annotation = FALSE,peaks = FALSE, pextend.upstream = 2000, extend.downstream = 200)
CoveragePlot(fetal_heart, region = "PAX1",features = "PAX1",annotation = FALSE,peaks = FALSE, pextend.upstream = 2000, extend.downstream = 200)
CoveragePlot(fetal_heart, region = "WT1",features = "WT1",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "CDH5",features = "CDH5",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 100)
CoveragePlot(fetal_heart, region = "NPPA",features = "NPPA",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)

CoveragePlot(fetal_heart, region = "TH",features = "TH",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "PHOX2B",features = "PHOX2B",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)

CoveragePlot(fetal_heart, region = "HCN4",features = "HCN4",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "TBX5",features = "TBX5",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "SHOX2",features = "SHOX2",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "CACNA1D",features = "CACNA1D",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "CACNA1G",features = "CACNA1G",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "CACNA2D2",features = "CACNA2D2",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)
CoveragePlot(fetal_heart, region = "GJC1",features = "GJC1",annotation = FALSE,peaks = FALSE, pextend.upstream = 1000, extend.downstream = 1000)


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

