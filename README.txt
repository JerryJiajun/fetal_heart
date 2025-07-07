R version 4.3.2 (2023-10-31)
Platform: x86_64-apple-darwin20 (64-bit)
Running under: macOS 15.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] dplyr_1.1.4                       BSgenome.Hsapiens.UCSC.hg38_1.4.5
 [3] BSgenome_1.70.2                   rtracklayer_1.62.0               
 [5] BiocIO_1.12.0                     Biostrings_2.70.3                
 [7] XVector_0.42.0                    EnsDb.Hsapiens.v86_2.99.0        
 [9] ensembldb_2.26.0                  AnnotationFilter_1.26.0          
[11] GenomicFeatures_1.54.4            AnnotationDbi_1.64.1             
[13] Biobase_2.62.0                    GenomicRanges_1.54.1             
[15] GenomeInfoDb_1.38.8               IRanges_2.36.0                   
[17] S4Vectors_0.40.2                  BiocGenerics_0.48.1              
[19] Signac_1.14.0                     patchwork_1.2.0                  
[21] ggplot2_3.5.1                     Seurat_5.1.0                     
[23] SeuratObject_5.0.2                sp_2.1-4                         

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22            splines_4.3.2               later_1.3.2                
  [4] bitops_1.0-8                filelock_1.0.3              tibble_3.2.1               
  [7] polyclip_1.10-7             XML_3.99-0.17               fastDummies_1.7.4          
 [10] lifecycle_1.0.4             globals_0.16.3              lattice_0.22-6             
 [13] MASS_7.3-60.0.1             magrittr_2.0.3              plotly_4.10.4              
 [16] yaml_2.3.10                 httpuv_1.6.15               sctransform_0.4.1          
 [19] spam_2.10-0                 spatstat.sparse_3.1-0       reticulate_1.39.0          
 [22] cowplot_1.1.3               pbapply_1.7-2               DBI_1.2.3                  
 [25] RColorBrewer_1.1-3          abind_1.4-8                 zlibbioc_1.48.2            
 [28] Rtsne_0.17                  purrr_1.0.2                 RCurl_1.98-1.16            
 [31] rappdirs_0.3.3              GenomeInfoDbData_1.2.11     ggrepel_0.9.6              
 [34] irlba_2.3.5.1               listenv_0.9.1               spatstat.utils_3.1-0       
 [37] goftest_1.2-3               RSpectra_0.16-2             spatstat.random_3.3-1      
 [40] fitdistrplus_1.2-1          parallelly_1.38.0           DelayedArray_0.28.0        
 [43] leiden_0.4.3.1              codetools_0.2-20            RcppRoll_0.3.1             
 [46] xml2_1.3.6                  tidyselect_1.2.1            matrixStats_1.4.1          
 [49] BiocFileCache_2.10.2        spatstat.explore_3.3-2      GenomicAlignments_1.38.2   
 [52] jsonlite_1.8.8              progressr_0.14.0            ggridges_0.5.6             
 [55] survival_3.7-0              tools_4.3.2                 progress_1.2.3             
 [58] ica_1.0-3                   Rcpp_1.0.13                 glue_1.7.0                 
 [61] SparseArray_1.2.4           gridExtra_2.3               MatrixGenerics_1.14.0      
 [64] withr_3.0.1                 fastmap_1.2.0               fansi_1.0.6                
 [67] digest_0.6.37               R6_2.5.1                    mime_0.12                  
 [70] colorspace_2.1-1            scattermore_1.2             tensor_1.5                 
 [73] spatstat.data_3.1-2         biomaRt_2.58.2              RSQLite_2.3.7              
 [76] utf8_1.2.4                  tidyr_1.3.1                 generics_0.1.3             
 [79] data.table_1.17.0           S4Arrays_1.2.1              prettyunits_1.2.0          
 [82] httr_1.4.7                  htmlwidgets_1.6.4           uwot_0.2.2                 
 [85] pkgconfig_2.0.3             gtable_0.3.5                blob_1.2.4                 
 [88] lmtest_0.9-40               htmltools_0.5.8.1           dotCall64_1.1-1            
 [91] ProtGenerics_1.34.0         scales_1.3.0                png_0.1-8                  
 [94] spatstat.univar_3.0-1       rstudioapi_0.16.0           reshape2_1.4.4             
 [97] rjson_0.2.21                nlme_3.1-166                curl_5.2.2                 
[100] zoo_1.8-12                  cachem_1.1.0                stringr_1.5.1              
[103] KernSmooth_2.23-24          parallel_4.3.2              miniUI_0.1.1.1             
[106] restfulr_0.0.15             pillar_1.9.0                grid_4.3.2                 
[109] vctrs_0.6.5                 RANN_2.6.2                  promises_1.3.0             
[112] dbplyr_2.5.0                xtable_1.8-4                cluster_2.1.6              
[115] cli_3.6.3                   compiler_4.3.2              Rsamtools_2.18.0           
[118] rlang_1.1.4                 crayon_1.5.3                future.apply_1.11.2        
[121] plyr_1.8.9                  stringi_1.8.4               viridisLite_0.4.2          
[124] deldir_2.0-4                BiocParallel_1.36.0         munsell_0.5.1              
[127] lazyeval_0.2.2              spatstat.geom_3.3-2         Matrix_1.6-5               
[130] RcppHNSW_0.6.0              hms_1.1.3                   bit64_4.0.5                
[133] future_1.34.0               KEGGREST_1.42.0             shiny_1.9.1                
[136] SummarizedExperiment_1.32.0 ROCR_1.0-11                 igraph_2.0.3               
[139] memoise_2.0.1               fastmatch_1.1-4             bit_4.0.5      



# instruction for multiomics data demonstration and visualization
put the raw data to the folder fetal_heart_multiomics
load the dependent libraries and run the R code under that folder
and below is the general plot function to demonstrate the clustering
# or you can load the saved rds file
fetal_heart <- readRDS("fetal_heart_multiomics_annotated.rds")
levels(fetal_heart) <- c('VCM','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','ACM')

# UMAP demonstrate different cell populations based on RNA expression
DimPlot(fetal_heart, reduction = "umap.rna", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("RNA") + labs(x = 'UMAP_1',y = 'UMAP_2')

# UMAP demonstrate different cell populations based on ATAC data
DimPlot(fetal_heart, reduction = "umap.atac.peaks",label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC") + labs(x = 'UMAP_1',y = 'UMAP_2')

# dotplot demonstrate the marker expression of each cell population
DefaultAssay(fetal_heart) <-'RNA'
levels(fetal_heart) <- c('VCM','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','ACM')
DotPlot(fetal_heart,features = c("MYL2","MYH11","MYBPC1","HCN4","HBG1","STMN2","CD163","CD96","PDGFRA","PAX1","WT1","CDH5","NPPA")) + RotatedAxis() + labs(x = '',y = '')

# UMAP plot showing the expression of indicated gene
DefaultAssay(fetal_heart) <-'RNA'
FeaturePlot(fetal_heart,features = c("TBX5")) + labs(x = 'UMAP_1',y = 'UMAP_2')

# Violin plot showing the expression of indicated gene
DefaultAssay(fetal_heart) <-'RNA'
VlnPlot(fetal_heart,features = c("TBX5"),pt.size = 0) + labs(x = '')

# IGV coverage plot showing the peaks at the regulatory region of  indicated gene
DefaultAssay(fetal_heart) <-'peaks'
levels(fetal_heart) <- c('VCM','VSM','SMC','SAN','RBC','Neuronal','Myeloid','Lymphoid','Fibroblast','Epithelial','Epicardial','Endothelial','ACM')
CoveragePlot(fetal_heart, region = "TBX5",features = "TBX5",annotation = TRUE,peaks = FALSE, pextend.upstream = 2000, extend.downstream = 2000)


# instruction for spatial data demonstration and visualization
put the raw data to the folder fetal_heart_spatial
load the dependent libraries and run the R code under that folder
and below is the general plot function to demonstrate the clustering
# or you can load the saved rds file

fetal_heart <- readRDS("fetal_heart_normal_spatial_annotated(006).rds")

# UMAP demonstrate different cell populations
DimPlot(fetal_heart,reduction = "umap",label = TRUE,label.size = 5) + labs(x = "UMAP_1",y = "UMAP_2")

# Spatial plot demonstrate different cell populations
SpatialDimPlot(fetal_heart,label = TRUE,label.size = 5,repel = TRUE) + theme(legend.position = "none")

# Spatial plot to show the marker expression of each cluster using Dotplot
levels(fetal_heart) <- c("VSM","SAN","Neuron","Lymphoid","CM")
DotPlot(fetal_heart,features = c("MYH11","SHOX2","PLP1",'CCL21',"ANGPT1")) + RotatedAxis() + labs(x = '',y = '')


#  Spatial plot showing the expression of indicated gene
SpatialFeaturePlot(fetal_heart,features = c("SHOX2"),alpha = c(0.1,1))

# UMAP plot showing the expression of indicated gene
FeaturePlot(fetal_heart,features = c("SHOX2")) + labs(x = "UMAP_1",y = "UMAP_2")

# Violin plot showing the expression of indicated gene
VlnPlot(fetal_heart,features = c("SHOX2"),pt.size = 0.5) + labs(x = '')

# Dot plot showing the expression of indicated gene
DotPlot(fetal_heart,features = c("SHOX2")) + RotatedAxis() + labs(x = '',y = '')


# instruction for Sinoid data demonstration and visualization
put the raw data to the folder Sinoid
load the dependent libraries and run the R code under that folder
and below is the general plot function to demonstrate the clustering
# or you can load the saved rds file
Sinoid <- readRDS("Sinoid_annotated.rds")
levels(Sinoid) <- c("SAN","Proliferating","Neuronal","Epithelial","Epicardial","Endothelial")

# UMAP demonstrate different cell populations based on RNA expression
# show in the default page
UMAPPlot(Sinoid,label=TRUE,label.size =5,repel = TRUE) + labs(x = 'UMAP_1', y = 'UMAP_2')

# dotplot demonstrate the marker expression of each cell population
# show in the default page
DotPlot(Sinoid,features = c("HCN4","TOP2A","STMN2","EPCAM","WT1","CDH5")) + RotatedAxis() + labs(x = '',y = '')

# Violin plot showing the expression of indicated gene
# show after type in the search gene 
VlnPlot(Sinoid,features = "HCN4",pt.size = 0.5) + labs(x = '')

# UMAP to demonstrate the expression of selected gene
# show after type in the search gene 
FeaturePlot(Sinoid,features = "HCN4") + labs(x = 'UMAP_1',y = "UMAP_2")

# Dotplot to demonstrate the expression of selected gene
# show after type in the search gene 
DotPlot(Sinoid,features = "HCN4") + labs(x = 'UMAP_1', y = 'UMAP_2')



# instruction for SAN-PCO data demonstration and visualization
put the raw data to the folder SAN_PCO
load the dependent libraries and run the R code under that folder
and below is the general plot function to demonstrate the clustering
# or you can load the saved rds file
# load the data
SAN_PCO <- readRDS("Combine_annotated.rds")
levels(SAN_PCO) <- c("SAN","Proliferating","Neuronal","Myofibroblasts","Mesenchymal","Epithelial","Endothelial","CM")

# UMAP demonstrate different cell populations based on RNA expression
# show in the default page
UMAPPlot(SAN_PCO,label=TRUE,label.size = 4, repel = TRUE) + labs(x = "UMAP_1",y = "UMAP_2")

# dotplot demonstrate the marker expression of each cell population
# show in the default page
DotPlot(SAN_PCO,features = c("HCN4","TOP2A","STMN2","COL1A1","ZFHX4","EPCAM","PECAM1","MYH7")) + RotatedAxis() + labs(x = "",y = "")


# Violin plot showing the expression of indicated gene
# show after type in the search gene 
VlnPlot(SAN_PCO,features = "STMN2",pt.size = 0.5) + labs(x = "")

# UMAP to demonstrate the expression of selected gene
# show after type in the search gene 
FeaturePlot(SAN_PCO,features = "STMN2") + labs(x = "UMAP_1",y = "UMAP_2")

# Dotplot to demonstrate the expression of selected gene
# show after type in the search gene 
DotPlot(SAN_PCO,features = "STMN2") + labs(x = "",y = "")

