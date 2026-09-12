#############################################################################
# decoupleR workflow: drug-target activity inference for standard
# (gene-level) single-cell/single-nucleus RNA-seq Seurat objects
#
# Template version -- fill in the placeholders marked "EDIT ME" below with
# your own Seurat object path, assay name, and grouping variable. Everything
# else (network construction, scoring, drug-name/FDA-approval resolution,
# and visualization) works the same regardless of dataset.
#############################################################################

## ------------------------------------------------------------------------
## 1. Installation
## ------------------------------------------------------------------------
## if (!require("BiocManager")) install.packages("BiocManager")
## BiocManager::install("decoupleR")
## BiocManager::install("OmnipathR")

library(decoupleR)
library(httr2)
library(OmnipathR)
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(pheatmap)
library(ggplot2)

## OmnipathR monkey-patches httr2's internal req_perform1() to attach curl
## stats for logging, but its patched signature (req, path, handle) is
## missing req_prep/resend_count, which current httr2 (>=1.2) requires --
## every retried request then fails with "unused arguments (req_prep,
## resend_count = n)". The patch is cosmetic (stats logging only), so
## restore httr2's real implementation to fix retries/downloads.
.httr2_ns <- asNamespace("httr2")
.original_req_perform1 <- environment(get("req_perform1", .httr2_ns))$original
if (!is.null(.original_req_perform1)) {
  unlockBinding("req_perform1", .httr2_ns)
  assign("req_perform1", .original_req_perform1, .httr2_ns)
  lockBinding("req_perform1", .httr2_ns)
}
rm(.httr2_ns, .original_req_perform1)

## OmnipathR resolves every organism argument (default 9606 for human)
## through ncbi_taxid() -> get_db("organisms"), which joins Ensembl/UniProt/
## OMA species tables. OMA's file changed format upstream and now fails to
## parse into columns, breaking that join even though 9606 is already a
## valid NCBI taxon ID needing no lookup at all. Short-circuit for numeric
## input so the broken OMA join is never triggered.
.omnipathr_ns <- asNamespace("OmnipathR")
.original_ncbi_taxid <- get("ncbi_taxid", .omnipathr_ns)
.patched_ncbi_taxid <- function(name) {
  is_num <- grepl("^[0-9]+$", as.character(name))
  out <- integer(length(name))
  if (any(is_num))  out[is_num]  <- as.integer(name[is_num])
  if (any(!is_num)) out[!is_num] <- .original_ncbi_taxid(name[!is_num])
  out
}
unlockBinding("ncbi_taxid", .omnipathr_ns)
assign("ncbi_taxid", .patched_ncbi_taxid, .omnipathr_ns)
lockBinding("ncbi_taxid", .omnipathr_ns)
rm(.omnipathr_ns, .original_ncbi_taxid, .patched_ncbi_taxid)


## ------------------------------------------------------------------------
## 2. Load your Seurat object                                     [EDIT ME]
## ------------------------------------------------------------------------
seurat_obj <- readRDS("~/fetal_heart_project/Heart_organoid_scRNA/mini_heart_annotated.rds")                     # <- path to your .rds file
seurat_obj <- UpdateSeuratObject(seurat_obj)

assay_name   <- "RNA"                             # <- e.g. "RNA", "SCT"
group_by_var <- "celltype"                             # <- e.g. "celltype", "seurat_clusters"

DefaultAssay(seurat_obj) <- assay_name


## ------------------------------------------------------------------------
## 3. Build the drug -> target network at the GENE level
## ------------------------------------------------------------------------
## "DrugBank" is not a valid OmnipathR resource name (returns 0 rows), and
## the DGIdb annotation records carry no drug identity at all -- just
## gene-level "category"/"source" druggability flags. The dataset that
## actually has drug/small-molecule -> protein edges is small_molecule()
## (pools CancerDrugsDB, SIGNOR, Guide2Pharma, Cellinker). Its `source` is
## a raw PubChem CID (resolved to a name in step 6 below).
sm_interactions <- OmnipathR::small_molecule()

## The same (compound, gene) pair is sometimes curated by more than one
## source with conflicting inhibition flags; decoupleR requires unique
## edges, so collapse to one row per pair (any source flagging inhibition
## wins, since that's the more specific claim).
net <- sm_interactions %>%
  transmute(
    source = as.character(source),
    target = target_genesymbol,
    mor    = ifelse(is_inhibition, -1, 1)
  ) %>%
  group_by(source, target) %>%
  summarise(mor = if (any(mor == -1)) -1 else 1, .groups = "drop")


## ------------------------------------------------------------------------
## 4. Pseudobulk to a gene x cluster/celltype matrix and score
## ------------------------------------------------------------------------
## decoupleR densifies its input matrix internally regardless of what you
## pass in, so for anything beyond a few thousand cells, pseudobulking by
## group_by_var keeps this from blowing up memory. If your object is small
## enough for true per-cell scoring, set group_by_var to a per-cell ID
## column instead -- just watch memory (genes x cells densified fully).
mat_gene <- AggregateExpression(
  seurat_obj,
  assays        = assay_name,
  group.by      = group_by_var,
  return.seurat = FALSE
)[[assay_name]]

# AggregateExpression() returns a sparse dgCMatrix; decoupleR expects (and
# will internally coerce to) a dense numeric matrix, so densify up front.
mat_gene <- as.matrix(mat_gene)

res_gene <- decouple(
  mat_gene,
  net,
  .source = "source",
  .target = "target",
  statistics      = c("wsum", "ulm", "aucell"),
  consensus_score = TRUE,
  minsize         = 3
)


## ------------------------------------------------------------------------
## 5. Reshape gene-collapsed activity to a drug x cluster/celltype matrix
## ------------------------------------------------------------------------
## res_gene is pseudobulked per group (not per cell), so it can't be
## embedded as a Seurat assay (an assay must share the parent object's
## cells) -- work with the drug x group matrix directly instead.
drug_activity <- res_gene %>%
  filter(statistic == "ulm") %>%
  pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()


## ------------------------------------------------------------------------
## 6. Resolve PubChem CIDs to drug names
## ------------------------------------------------------------------------
## small_molecule()'s `source` is a raw PubChem CID with no name attached,
## and OmnipathR has no local CID -> name table. PubChem's PUG REST API
## resolves them live (batched, no key required), so look up only the CIDs
## that actually survived scoring rather than the whole network.
pubchem_cid_to_name <- function(cids, batch_size = 100) {
  cids <- unique(cids)
  batches <- split(cids, ceiling(seq_along(cids) / batch_size))

  purrr::map_dfr(batches, function(batch) {
    url <- sprintf(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/property/Title/JSON",
      paste(batch, collapse = ",")
    )
    resp <- tryCatch(httr2::request(url) %>% httr2::req_perform(), error = function(e) NULL)
    Sys.sleep(0.25)  # stay well under PubChem's rate limit
    if (is.null(resp) || httr2::resp_status(resp) != 200) return(tibble())

    props <- httr2::resp_body_json(resp)$PropertyTable$Properties
    tibble(
      source    = vapply(props, function(p) as.character(p$CID), character(1)),
      drug_name = vapply(props, function(p) if (is.null(p$Title)) NA_character_ else p$Title, character(1))
    )
  })
}

cid_names <- pubchem_cid_to_name(res_gene$source)

res_gene <- res_gene %>%
  left_join(cid_names, by = "source") %>%
  mutate(drug_label = coalesce(drug_name, source))


## ------------------------------------------------------------------------
## 7. Restrict to FDA-approved medications                        [OPTIONAL]
## ------------------------------------------------------------------------
## PubChem's PUG View "FDA Approved Drugs" heading is backed by the FDA's
## own Drugs@FDA database (HTTP 200 = has an approved-drug record, 404 =
## no record). Unlike the name lookup, this endpoint has no batch mode --
## one request per CID -- so results are cached to disk so re-running the
## script later doesn't repeat the check.
## Caveat: PubChem sometimes attaches this heading only to a specific salt
## form's CID rather than the free-base CID the network uses (e.g.
## doxorubicin is approved as "doxorubicin hydrochloride", a different
## CID, while the free-base CID has no heading of its own). This is a
## reasonable proxy, not a perfect one -- spot-check on PubChem directly
## if a drug you expected to see is missing.
fda_cache_path <- "fda_approval_cache.csv"           # <- EDIT ME if desired

check_fda_approved <- function(cid) {
  url <- sprintf(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON?heading=FDA%%20Approved%%20Drugs",
    cid
  )
  resp <- tryCatch(
    httr2::request(url) %>% httr2::req_error(is_error = function(resp) FALSE) %>% httr2::req_perform(),
    error = function(e) NULL
  )
  Sys.sleep(0.2)  # stay well under PubChem's rate limit
  if (is.null(resp)) return(NA)
  httr2::resp_status(resp) == 200
}

fda_cache <- if (file.exists(fda_cache_path)) {
  read.csv(fda_cache_path, stringsAsFactors = FALSE) %>% mutate(source = as.character(source))
} else {
  tibble(source = character(), fda_approved = logical())
}

cids_to_check <- setdiff(unique(res_gene$source), fda_cache$source)
if (length(cids_to_check) > 0) {
  new_checks <- tibble(
    source       = cids_to_check,
    fda_approved = vapply(cids_to_check, check_fda_approved, logical(1))
  )
  # don't cache failed (NA) lookups so they're retried on the next run
  fda_cache <- bind_rows(fda_cache, filter(new_checks, !is.na(fda_approved)))
  write.csv(fda_cache, fda_cache_path, row.names = FALSE)
}

approved_cids <- fda_cache$source[fda_cache$fda_approved]

## Comment out the next two lines to keep all compounds, not just
## FDA-approved ones.
res_gene      <- res_gene      %>% filter(source %in% approved_cids)
drug_activity <- drug_activity[rownames(drug_activity) %in% approved_cids, , drop = FALSE]


## ------------------------------------------------------------------------
## 8. Visualize
## ------------------------------------------------------------------------

## Pick compounds by activity variance across groups.
top_compounds <- names(sort(apply(drug_activity, 1, var, na.rm = TRUE), decreasing = TRUE))[1:100]

## Dot plot of drug activity across clusters/celltypes
res_gene %>%
  filter(statistic == "ulm", source %in% top_compounds) %>%
  ggplot(aes(x = condition, y = drug_label, size = abs(score), color = score)) +
  geom_point() +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = group_by_var, y = "Drug", size = "|score|", color = "Score")

## Heatmap of drug activity per cluster/celltype, labeled with drug names
drug_activity_labeled <- drug_activity
rownames(drug_activity_labeled) <- coalesce(
  cid_names$drug_name[match(rownames(drug_activity), cid_names$source)],
  rownames(drug_activity)
)

pheatmap(
  drug_activity_labeled[1:50,],
  scale = "row",cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50)
)

drug_activity_labeled_reorder <- drug_activity_labeled[,c(1,9,8,2,3,4,5,6,7)]
pheatmap(
  drug_activity_labeled_reorder[1:50,],
  scale = "row",cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50)
)

drug_activity_labeled_selected <- drug_activity_labeled[,c(1,9,8)]
pheatmap(
  drug_activity_labeled_selected[1:50,],
  scale = "row",cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50)
)


plot_top_drugs_heatmap <- function(mat, top_n = 10, scale = "row",
                                   cluster_rows = FALSE, cluster_cols = TRUE, ...) {
  
  celltypes <- colnames(mat)
  
  # 1. Top_n drug names per celltype
  top_list <- lapply(celltypes, function(ct) {
    ord <- order(mat[, ct], decreasing = TRUE)[1:top_n]
    rownames(mat)[ord]
  })
  names(top_list) <- celltypes
  
  # 2. Union of top drugs
  top_drugs <- unique(unlist(top_list))
  
  # 3. Primary celltype = first celltype (in column order) each drug is a top hit for
  primary_ct <- sapply(top_drugs, function(d) {
    celltypes[sapply(top_list, function(x) d %in% x)][1]
  })
  
  # 4. Order drugs by primary celltype
  ord_idx <- order(factor(primary_ct, levels = celltypes))
  top_drugs_ordered <- top_drugs[ord_idx]
  primary_ct_ordered <- primary_ct[ord_idx]
  
  # 5. Annotation
  annotation_row <- data.frame(
    CellType = factor(primary_ct_ordered, levels = celltypes)
  )
  rownames(annotation_row) <- top_drugs_ordered
  
  # 6. Gaps between celltype blocks — fill missing counts with 0, drop empty groups
  counts <- table(factor(primary_ct_ordered, levels = celltypes))
  counts[is.na(counts)] <- 0
  gaps <- cumsum(counts)
  gaps <- gaps[gaps > 0]                # drop leading zero-count groups
  gaps <- gaps[-length(gaps)]           # no gap after the last block
  gaps <- unique(as.integer(gaps))      # avoid duplicate gap positions
  
  # 7. Plot
  pheatmap(
    mat[top_drugs_ordered, , drop = FALSE],
    scale = scale,
    color = colorRampPalette(c("navy", "white", "firebrick"))(50),
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    annotation_row = annotation_row,
    gaps_row = gaps,
    fontsize_row = 7,
    fontsize_col = 9,
    cellwidth = 20,
    cellheight = 9,
    ...
  )
}

plot_top_drugs_heatmap(drug_activity_labeled, top_n = 30)

## Cluster/celltype-specific top drugs. res_gene is pseudobulked per group
## with one score per (drug, group) pair -- there's no per-cell replication
## for FindAllMarkers to test, so rank directly by score within each group.
top_drugs_per_group <- res_gene %>%
  filter(statistic == "ulm") %>%
  group_by(condition) %>%
  slice_max(score, n = 10) %>%
  ungroup()

print(top_drugs_per_group)


## ------------------------------------------------------------------------
## 9. Drugs specifically targeting one cluster/celltype of interest [EDIT ME]
## ------------------------------------------------------------------------
group_of_interest <- "SAN"                        # <- e.g. "AVN", or a cluster number as character

group_top_drugs <- res_gene %>%
  filter(statistic == "ulm", condition == group_of_interest) %>%
  slice_max(abs(score), n = 30)

ggplot(group_top_drugs, aes(x = reorder(drug_label, score), y = score, fill = score)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick") +
  theme_bw(base_size = 13) +
  labs(x = "Drug", y = "Activity score (ulm)", fill = "Score",
       title = paste("Top drugs targeting", group_of_interest))



## ------------------------------------------------------------------------
## 10. Get targets of selected drugs                               [EDIT ME]
## ------------------------------------------------------------------------
## Works from either drug names (resolved labels) or raw PubChem CIDs.
## Defaults to the drugs pulled out in section 9 (top hits for
## group_of_interest); swap in any character vector of drug_label/CID
## values you want instead.
selected_drugs <- c("Vorasidenib","Pomalidomide","Selpercatinib")
selected_drugs <- c("Avapritinib","Vorasidenib","Pomalidomide","Selpercatinib") 
selected_drugs <- c("Avapritinib","Chlorambucil") 
selected_drugs <- c("Ibrutinib") 
selected_drugs <- c("Dasatinib","Afatinib") 
selected_drugs <- c("Acetylcholine") 

## Map back to CIDs (drug_label may already *be* a CID if the PubChem
## lookup failed for that compound -- coalesce handles both cases).
selected_cids <- cid_names %>%
  filter(drug_name %in% selected_drugs | source %in% selected_drugs) %>%
  pull(source) %>%
  union(intersect(selected_drugs, sm_interactions$source)) %>%
  unique()

## Pull full annotation from sm_interactions (not just `net`, which was
## already collapsed to one row per compound-gene pair and dropped the
## evidence/reference info) so the source of each target call is visible.
drug_targets <- sm_interactions %>%
  filter(source %in% selected_cids) %>%
  transmute(
    source          = as.character(source),
    target_gene     = target_genesymbol,
    direction       = ifelse(is_inhibition, "inhibition",
                             ifelse(is_stimulation, "activation", "unspecified")),
    n_references    = lengths(strsplit(as.character(references), ";")),
    resource        = sources
  ) %>%
  left_join(cid_names, by = "source") %>%
  mutate(drug_label = coalesce(drug_name, source)) %>%
  distinct(drug_label, target_gene, direction, .keep_all = TRUE) %>%
  arrange(drug_label, desc(n_references))

print(drug_targets)

## ------------------------------------------------------------------------
## Plot A: Drug -> target map
## ------------------------------------------------------------------------
## Tile plot (not dot plot) so direction is legible by color alone; each
## drug is a row, each target gene a column. Faceting isn't needed since
## drug and target are both on-axis already.
ggplot(drug_targets, aes(x = target_gene, y = drug_label, fill = direction)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c(inhibition = "navy", activation = "firebrick",
                               unspecified = "grey70")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) +
  labs(x = "Target gene", y = "Drug", fill = "Effect",
       title = "Targets of selected drugs")


## ------------------------------------------------------------------------
## Plot B: Target gene expression across selected celltypes
## ------------------------------------------------------------------------
groups_of_interest <- c("SAN", "ACM", "VCM")         # <- EDIT ME: subset to specific
#    celltypes/clusters if you
#    don't want all of them,
#    e.g. c("SAN", "AVN", "CM")
target_genes <- unique(drug_targets$target_gene)
target_genes_present <- intersect(target_genes, rownames(mat_gene))
missing_targets <- setdiff(target_genes, target_genes_present)
if (length(missing_targets) > 0) {
  message("Note: ", length(missing_targets),
          " target gene(s) not found in mat_gene and will be omitted: ",
          paste(missing_targets, collapse = ", "))
}
target_expr <- mat_gene[target_genes_present, groups_of_interest, drop = FALSE]

## Drop genes with zero variance across the selected groups -- row-scaling
## (z-score) is undefined for them and would produce NA/NaN dots.
zero_var <- apply(target_expr, 1, function(x) var(x, na.rm = TRUE) == 0)
if (any(zero_var)) {
  message("Note: ", sum(zero_var),
          " target gene(s) have zero variance across selected groups and ",
          "will be omitted: ",
          paste(rownames(target_expr)[zero_var], collapse = ", "))
  target_expr <- target_expr[!zero_var, , drop = FALSE]
}

## Reshape to long format: one row per (gene, celltype), with both the raw
## pseudobulk value (for dot size) and the row-wise z-score (for color) --
## z-scoring within gene makes relative enrichment across celltypes
## comparable between genes with very different expression scales.
target_expr_long <- target_expr %>%
  as.data.frame() %>%
  rownames_to_column("target_gene") %>%
  pivot_longer(-target_gene, names_to = "celltype", values_to = "expr") %>%
  group_by(target_gene) %>%
  mutate(expr_scaled = as.numeric(scale(expr))) %>%
  ungroup()

## Order genes so the heatmap-like grouping (which celltype each gene peaks
## in) is preserved instead of alphabetical -- keeps the dot plot readable.
gene_order <- target_expr_long %>%
  group_by(target_gene) %>%
  slice_max(expr, n = 1) %>%
  arrange(factor(celltype, levels = groups_of_interest)) %>%
  pull(target_gene) %>%
  unique()

target_expr_long <- target_expr_long %>%
  mutate(target_gene = factor(target_gene, levels = rev(gene_order)),
         celltype    = factor(celltype, levels = groups_of_interest))

ggplot(target_expr_long, aes(x = celltype, y = target_gene,
                             size = expr, color = expr_scaled)) +
  geom_point() +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick") +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = "grey90")) +
  labs(x = group_by_var, y = "Target gene", size = "Pseudobulk\nexpression",
       color = "Row z-score",
       title = "Expression of drug target genes across celltypes")

## ------------------------------------------------------------------------
## Plot C: Per-drug target dot plots (pct expressed + avg expression)
## ------------------------------------------------------------------------
## Uses Seurat::DotPlot on the original per-cell object (not mat_gene),
## since pct-expressed is only meaningful at single-cell resolution.
## Subset to the celltypes of interest first so the plot only reflects them.
cells_keep <- WhichCells(seurat_obj, expression = !!sym(group_by_var) %in% groups_of_interest)
seurat_sub <- subset(seurat_obj, cells = cells_keep)

## Only keep target genes actually present as features
target_genes_present <- intersect(unique(drug_targets$target_gene),
                                  rownames(seurat_sub[[assay_name]]))

dp <- DotPlot(seurat_sub, features = target_genes_present,
              group.by = group_by_var, assay = assay_name)

## Pull the underlying data and attach which drug(s) each gene belongs to,
## so a gene targeted by multiple selected drugs appears in each facet.
dp_data <- dp$data %>%
  rename(target_gene = features.plot, celltype = id) %>%
  left_join(drug_targets %>% distinct(drug_label, target_gene),
            by = "target_gene")

ggplot(dp_data, aes(x = celltype, y = target_gene,
                    size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  scale_color_gradient2(low = "navy", mid = "white", high = "firebrick") +
  facet_wrap(~ drug_label, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = group_by_var, y = "Target gene", size = "% expressed",
       color = "Avg expression\n(scaled)")

#############################################################################
# End of workflow
#############################################################################
