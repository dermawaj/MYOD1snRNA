
#####tumorigenic MRs########################
#differential expression between MYOD1 sc vs normal skeletal muscle####
#####diff expression (cell level and pseudo-bulk) by response######
CPM_normalization <- function(dset){
  dset.log2cpm <- apply(dset,2,function(x){
    y <- 1E6*x/sum(x) + 1
    z <- log(y,2)
    return(z)
  })
  return(dset.log2cpm)
}

GES_scaled <- function(dset, ref){
  ref.mean <- apply(ref, 1, mean)
  ref.sd <- apply(ref, 1, sd)
  dset.ges <- apply(dset, 2, function(x){(x - ref.mean) / ref.sd})
  dset.ges <- dset.ges[is.finite(rowSums(dset.ges)),]
  return(dset.ges)
}

muscle.raw.obj <- readRDS(file = "/path/scRNA/RMS/tabulaSapiens_muscle.rds")
muscle.dat <- as.matrix(muscle.raw.obj@assays$RNA@counts)
muscle.cpm <- CPM_normalization(muscle.dat)
rm(muscle.dat)

mart <- readRDS("/path/scRNA/functions/mart.obj.rds")

sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
nets <- lapply(filenames, readRDS)
names(nets) <- sample
setwd("/path/scRNA/RMS")
sample.obj <- lapply(sample, function(s) readRDS(file = paste0(s, '_subset.rds')))

all_vp_RMS <- list()
for (s in seq_along(sample)) {
  MYOD1.dat <- as.matrix(sample.obj[[s]]@assays$RNA@layers$counts)
  genes <- row.names(sample.obj[[s]])
  gene_ensembl <- getBM(filters = c("hgnc_symbol"), attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        values = genes, mart = mart)
  ensembl <- gene_ensembl$ensembl_gene_id[match(genes, gene_ensembl$hgnc_symbol)]
  row.names(MYOD1.dat) <- ensembl
  MYOD1.dat <- MYOD1.dat[!is.na(rownames(MYOD1.dat)), ]
  colnames(MYOD1.dat) <- colnames(sample.obj[[s]])
  MYOD1.cpm <- CPM_normalization(MYOD1.dat) #run CPM_normalization function on raw count matrix sample by sample
  common_genes <- intersect(rownames(MYOD1.cpm), rownames(muscle.cpm))
  MYOD1.cpm <- MYOD1.cpm[match(common_genes, rownames(MYOD1.cpm)), ]
  muscle.cpm <- muscle.cpm[match(common_genes, rownames(muscle.cpm)), ]
  dat <- GES_scaled(MYOD1.cpm, muscle.cpm) #differential gene expression relative to muscle
  rm(MYOD1.dat, MYOD1.cpm)
  gc()
  colnames(dat) <- paste0(colnames(dat), "_", s)
  indices <- seq(1, ncol(dat), by = 400)
  seurat_viper_list <- list()
  for (i in 1:(length(indices) - 1)) {
    vp <- viper(dat[, indices[i]:indices[i + 1]], nets, method = 'none') #use this scaled signature in viper using methods = none 
    seurat_viper_list <- c(seurat_viper_list, list(vp))
  }
  vp <- viper(dat[, indices[length(indices)]:ncol(dat)], nets, method = 'none')
  seurat_viper_list <- c(seurat_viper_list, list(vp))
  vp_RMS <- seurat_viper_list[[1]]
  for (i in 2:length(seurat_viper_list)) {
    vp_RMS <- cbind(vp_RMS, seurat_viper_list[[i]])
  }
  all_vp_RMS[[sample[s]]] <- vp_RMS
  setwd("/path/scRNA/RMS/seurat/ARACNe")
  saveRDS(vp_RMS, paste0(sample[s], '_diffexp_muscle.viper.sample_network.rds'))
}
combined_vp_RMS <- do.call(cbind, all_vp_RMS)
genes <- row.names(combined_vp_RMS)
gene_hgnc <- getBM(filters = c("ensembl_gene_id"), attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      values = genes, mart = mart)
hgnc <- gene_hgnc$hgnc_symbol[match(genes, gene_hgnc$ensembl_gene_id)]
row.names(combined_vp_RMS) <- hgnc
combined_vp_RMS <- combined_vp_RMS[!is.na(rownames(combined_vp_RMS)), ]
setwd("/path/scRNA/RMS/seurat/ARACNe")
saveRDS(combined_vp_RMS, file = "MYOD1_muscle_diffexp.viper.sample_network.rds")
                     
##########################################
###MR matrix###
vpmat <- readRDS("MYOD1_muscle_diffexp.viper.sample_network.rds")

common_cells <- intersect(colnames(vpmat), colnames(MYOD1.integrated.vp))
vpmat <- vpmat[, common_cells]
sample_clusters <- unique(MYOD1.integrated.vp$sample_cluster)
results_list <- list()
for (cluster in sample_clusters) {
  cluster_cells <- colnames(MYOD1.integrated.vp)[MYOD1.integrated.vp$sample_cluster == cluster]
  cluster_cells <- intersect(cluster_cells, colnames(vpmat))
  if (length(cluster_cells) > 0) {
    cluster_matrix <- vpmat[, cluster_cells, drop = FALSE]
    cluster_means <- rowMeans(cluster_matrix, na.rm = TRUE)
    results_list[[cluster]] <- cluster_means
  }
}
combined_results <- do.call(cbind, results_list)
combined_matrix <- as.matrix(combined_results) 

stouffersMethod <- function(x, weights) {
  return(sum(x * weights) / sqrt(sum(weights * weights)))
}
df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv") %>% arrange(cluster)
cluster <- df$subcluster
group <- df$cluster
colnames(combined_matrix) <- gsub("_cluster_", ".cluster", colnames(combined_matrix))
combined_matrix <- combined_matrix[, match(cluster, colnames(combined_matrix))]
grouped_clusters <- split(colnames(combined_matrix), group)
results_list <- list()
for (group_name in names(grouped_clusters)) {
  subcluster_names <- grouped_clusters[[group_name]]
  subcluster_matrix <- combined_matrix[, subcluster_names, drop = FALSE]
  weights <- rep(1, ncol(subcluster_matrix))
  integrated_results <- apply(subcluster_matrix, 1, stouffersMethod, weights = weights)
  results_list[[group_name]] <- integrated_results
}
integrated_matrix <- do.call(cbind, results_list)
                     
##topMRs###
get_top_genes <- function(matrix, n = 20) {
  selected_genes <- unique(unlist(lapply(seq_len(ncol(matrix)), function(i) {
    col_data <- matrix[, i]
    names(col_data) <- rownames(matrix)
    top_genes <- names(sort(col_data, decreasing = TRUE)[seq_len(min(n, length(col_data)))])
    bottom_genes <- names(sort(col_data, decreasing = FALSE)[seq_len(min(n, length(col_data)))])
    # print(top_genes)
    return(c(top_genes,bottom_genes))
  })))
  return(selected_genes)
}
combined_matrix <- integrated_matrix
(selected_genes <- get_top_genes(combined_matrix, n = 20))

                     
############
###Oncotarget#####
# convert NES scores to pvalues ( < 10 ^-5 converted to log scale)
# OncoTarget analysis for patients, PDX, and cell lines
# Subset viper matrix for 180 regulators that have known FDA-approved or late-stage experimental drugs directly targeting them

filtered_matrix <- combined_matrix[selected_genes, ]
rownames(filtered_matrix)

OncoTarget <- read.csv("/path/oncotarget.csv")
top_dat <- filtered_matrix[rownames(filtered_matrix) %in% OncoTarget$Target,]
rownames(top_dat)
vpmat.oncotarget.pval <- apply(top_dat, 2, function(x) -log10(p.adjust(pnorm(x, lower.tail=FALSE),method='bonferroni')))
vpmat <- vpmat.oncotarget.pval >= 1 # at a fdr-corrected p-value of < 10E-1
pos <- rowSums(vpmat) >= 1 #can be 1, significant in at least one sample
vpmat <- filterRowMatrix(vpmat.oncotarget.pval,pos)
vpmat <- vpmat[order(rowMeans(vpmat), decreasing = TRUE), ]

paletteLength <- 100
myColor <- colorRampPalette(c("white", "lightcoral", "red", "darkred"))(paletteLength)
myBreaks <- c(seq(0, 2, length.out = ceiling(paletteLength * 0.2) + 1), 
              seq(2.01, 5, length.out = ceiling(paletteLength * 0.3)), 
              seq(5.01, 10, length.out = ceiling(paletteLength * 0.3)), 
              seq(10.01, max(vpmat, na.rm = TRUE), length.out = ceiling(paletteLength * 0.2)))

pdf("MYOD1_muscle_diffexp_oncotarget_integrated_heatmap.pval.pdf", height = 5, width = 3)
pheatmap(vpmat, 
          cluster_rows = F, 
          treeheight_row = 0, 
          cluster_cols = F, 
          show_rownames = TRUE, 
          show_colnames = TRUE, 
          fontsize_row = 8, 
          fontsize_col = 6,
          color = myColor,  
          breaks = myBreaks,
          show_row_dend = F,
          main = "")
dev.off()

############


##########################################
###bulk vs normal muscle sample by sample####
setwd("/path/scRNA/RMS")
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
nets <- lapply(filenames, readRDS)
names(nets) <- sample

#normal muscle (Tabula Sapiens)
muscle.obj <- readRDS(file = "/path/scRNA/RMS/tabulaSapiens_muscle.rds.rds")
muscle.obj.dat <- as.matrix(muscle.obj@assays$RNA@counts)
colnames(muscle.obj.dat) <- colnames(muscle.obj)
rownames(muscle.obj.dat) <- rownames(muscle.obj)

#Run bulk vs normal muscle diff gene expression, run viper on differentially exp genes, sample by sample. scale each gene by Mean and SD for each sample:
setwd("/path/scRNA/RMS/MYOD1_bulk")
MYOD1_bulk <- read.delim("MYOD1_COUNT.txt") 
MYOD1_bulk$SAMPLESID <- sub("\\..*", "", MYOD1_bulk$SAMPLESID)
MYOD1_bulk <- MYOD1_bulk %>% column_to_rownames("SAMPLESID")
dat <- as.matrix(MYOD1_bulk)
#normalize counts
CPM_normalization <- function(dset){
  dset.log2cpm <- apply(dset,2,function(x){
    y <- 1E6*x/sum(x) + 1
    z <- log(y,2)
    return(z)
  })
  return(dset.log2cpm)
}
dset.patient.log2cpm <- CPM_normalization(dat)
dset.muscle.log2cpm <- CPM_normalization(muscle.obj.dat)

# Get common genes between the bulk patient samples and the muscle dataset
common_genes <- intersect(rownames(dset.patient.log2cpm), rownames(dset.muscle.log2cpm))
dset.patient.log2cpm <- dset.patient.log2cpm[common_genes,]
dset.muscle.log2cpm <- dset.muscle.log2cpm[common_genes,]

# Scale function (take test and reference data set, create mean from reference 
# Generate gene expression signature (how many normalized values into a z score how far away is from normal muscle)
GES_scaled <- function(dset, ref){
  ref.mean <- apply(ref, 1, mean)
  ref.sd <- apply(ref, 1, sd)
  dset.ges <- apply(dset, 2, function(x){(x - ref.mean) / ref.sd})
  dset.ges <- dset.ges[is.finite(rowSums(dset.ges)),]
  return(dset.ges)
}
dset.ges <- GES_scaled(dset = dset.patient.log2cpm, ref = dset.muscle.log2cpm)
write.csv(dset.ges, file = "MYOD1_bulk_vs_muscle_diffexp_GES_scaled.csv", row.names = TRUE)

vpmat <- viper(dset.ges, regulon = nets, cores = 2, verbose = TRUE, method = 'none')

genes <- row.names(vpmat)
gene_hgnc <- getBM(filters= c("ensembl_gene_id"), attributes= c("ensembl_gene_id","hgnc_symbol"),
                    values=genes, mart= mart)
hgnc <- gene_hgnc$hgnc_symbol[match(genes, gene_hgnc$ensembl_gene_id)]
row.names(vpmat) <- hgnc

rows_to_remove <- grep("^RPS|^RPL", rownames(vpmat), value = TRUE)
vpmat <- vpmat[!rownames(vpmat) %in% rows_to_remove, ]

###########
###Oncotarget#####
# convert NES scores to pvalues ( < 10 ^-5 converted to log scale)
# OncoTarget analysis for patients, PDX, and cell lines
# Subset viper matrix for 180 regulators that have known FDA-approved or late-stage experimental drugs directly targeting them

get_top_genes <- function(matrix, n = 20) {
  selected_genes <- unique(unlist(lapply(seq_len(ncol(matrix)), function(i) {
    col_data <- matrix[, i]
    names(col_data) <- rownames(matrix)
    top_genes <- names(sort(col_data, decreasing = TRUE)[seq_len(min(n, length(col_data)))])
    print(top_genes)
    return(top_genes)
  })))
  return(selected_genes)
}
(selected_genes <- get_top_genes(vpmat, n = 20))
selected_genes

filtered_matrix <- vpmat[selected_genes, ]
filtered_matrix <- vpmat
rownames(filtered_matrix)

OncoTarget <- read.csv("/path/oncotarget.csv")
top_dat <- filtered_matrix[rownames(filtered_matrix) %in% OncoTarget$Target,]
rownames(top_dat)
vpmat.oncotarget.pval <- apply(top_dat, 2, function(x) -log10(p.adjust(pnorm(x, lower.tail=FALSE),method='bonferroni')))
vpmat <- vpmat.oncotarget.pval >= 1 # at a fdr-corrected p-value of < 10E-1
pos <- rowSums(vpmat) >= 1 #can be 1, significant in at least one sample
vpmat <- filterRowMatrix(vpmat.oncotarget.pval,pos)
vpmat <- vpmat[order(rowMeans(vpmat), decreasing = TRUE), ]

paletteLength <- 100
myColor <- colorRampPalette(c("white", "lightcoral", "red", "darkred"))(paletteLength)
myBreaks <- c(seq(0, 2, length.out = ceiling(paletteLength * 0.2) + 1), 
              seq(2.01, 5, length.out = ceiling(paletteLength * 0.3)), 
              seq(5.01, 10, length.out = ceiling(paletteLength * 0.3)), 
              seq(10.01, max(vpmat, na.rm = TRUE), length.out = ceiling(paletteLength * 0.2)))

pdf("MYOD1_bulk_vs_muscle_diffexp_genes_viper_oncotarget_heatmap.pval.pdf", height = 5, width = 8)
pheatmap(vpmat, 
         cluster_rows = F, 
         cluster_cols = F, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 6, 
         fontsize_col = 6,
         color = myColor,
         breaks = myBreaks,
         annotation_col = annotation,
         annotation_colors = anno.colors,
         gaps_col = gaps_col,
         main = "")
dev.off()

###########
