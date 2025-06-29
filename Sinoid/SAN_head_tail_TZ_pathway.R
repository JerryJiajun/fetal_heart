# load the DEG data
dea<-read.csv("SAN_head_vs_tail_DEA.csv",row.names = 1)
# dea<-read.csv("SAN_head_vs_TZ_DEA.csv",row.names = 1)
# dea<-read.csv("SAN_tail_vs_TZ_DEA.csv",row.names = 1)

dea_sig<-subset(dea,p_val_adj < .05)
# dea_sig<-subset(dea,PValue < .05 & abs(logFC) > 1)

# load the match table which have the gene_name, gene_id,and entrez_gene_id
table_human<-read.table("table_human_index.csv",header = TRUE,sep = ",",row.names = 1)

id<-match(rownames(dea_sig),table_human$external_gene_name)

dea_sig$entrez_id<-table_human$entrezgene_id[id]

# pathway analysis using pathview
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("pathview")
# BiocManager::install("gage")
# BiocManager::install("gageData")
# KEGG pathway analysis
library(AnnotationDbi)
library(org.Hs.eg.db) # options(connectionObserver = NULL)
library(pathview)
library(gage)
library(gageData)

# The gageData package has pre-compiled databases mapping genes to KEGG pathways and GO terms for common organisms
data(kegg.sets.hs)
data(sigmet.idx.hs)
# kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
kg.hs<-kegg.gsets("hsa")
kegg.sets.hs<-kg.hs$kg.sets[kg.hs$sigmet.idx]
head(kegg.sets.hs, 3)


# Run the pathway analysis
foldchanges <- dea_sig$avg_log2FC
names(foldchanges) <- dea_sig$entrez_id
keggres <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)

# extract the top 20 upregulated pathways and downregulated pathways
up_reg<-as.data.frame(keggres$greater)
up_reg <- up_reg[order(up_reg$p.val, -up_reg$set.size, decreasing = FALSE), ]
up_reg$description<-rownames(up_reg)

down_reg<-as.data.frame(keggres$less)
down_reg<-down_reg[order(down_reg$p.val,-down_reg$set.size,decreasing = FALSE),]
down_reg$description<-rownames(down_reg)

# plot the upregulated and downregulated pathway
library(ggplot2)
ggplot(up_reg[1:15,], aes(x=description, y=set.size,fill=p.val)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  coord_flip()

ggplot(up_reg[1:10,], aes(x=p.val,y=reorder(description,-log10(p.val)),color=p.val,size=set.size)) + 
  geom_point()+
  scale_color_gradient(low = "red", high = "blue") + 
  theme_bw() + 
  xlab('pvalue')+
  ylab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggplot(down_reg[1:15,], aes(x=description, y=set.size,fill=p.val)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  coord_flip()

ggplot(down_reg[1:10,], aes(x=p.val,y=reorder(description,-log10(p.val)),color=p.val,size=set.size)) + 
  geom_point()+
  scale_color_gradient(low = "red", high = "blue") + 
  theme_bw() + 
  xlab('pvalue')+
  ylab("")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Pull out the top 10 upregulated pathways, then further process that just to get the IDs
library(dplyr)

# Get the pathways
keggrespathways <- data.frame(id=rownames(keggres$greater), keggres$greater) %>%
  tbl_df() %>%
  filter(row_number()<=10) %>%
  .$id %>%
  as.character()
keggrespathways

# Get the IDs.
keggresids <- substr(keggrespathways, start=1, stop=8)
keggresids


# the pathview() function in the pathview package makes the plots
# Define plotting function for applying later
plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# Unload dplyr since it conflicts with the next line
detach("package:dplyr", unload=T)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))



# Pull out the top 10 downregulated pathways, then further process that just to get the IDs
library(dplyr)

# Get the pathways
keggrespathways <- data.frame(id=rownames(keggres$less), keggres$less) %>%
  tbl_df() %>%
  filter(row_number()<=10) %>%
  .$id %>%
  as.character()
keggrespathways

# Get the IDs.
keggresids <- substr(keggrespathways, start=1, stop=8)
keggresids


# the pathview() function in the pathview package makes the plots
# Define plotting function for applying later
plot_pathway <- function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)

# Unload dplyr since it conflicts with the next line
detach("package:dplyr", unload=T)

# plot multiple pathways (plots saved to disk and returns a throwaway list object)
tmp <- sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))



