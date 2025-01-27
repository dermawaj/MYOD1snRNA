
#####tumorigenic MRs########################
#differential expression between MYOD1 sc vs normal skeletal muscle####
#####diff expression (cell level and pseudo-bulk) by response######
option(Seurat.object.assay.version = 'v5')
mart <- readRDS("/path/scRNA/functions/mart.obj.rds")
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.filt.integrated.v2.rds")
MYOD1.obj <- CreateAssay5Object(counts = MYOD1.obj.integrated[["SCT"]]@counts, data = MYOD1.obj.integrated[["SCT"]]@data)
##
muscle.raw.obj <- readRDS(file = "/path/scRNA/RMS/9b34eb5c-1fb2-46c5-9fe4-d90a2e5f7c4b.rds")
muscle.raw.obj <- SCTransform(muscle.raw.obj, verbose = T, conserve.memory = T)
muscle.obj <- CreateAssay5Object(counts = muscle.raw.obj[["SCT"]]@counts, data = muscle.raw.obj[["SCT"]]@data)

genes <- row.names(muscle.obj)
gene_hgnc <- getBM(filters= c("ensembl_gene_id"), attributes= c("ensembl_gene_id","hgnc_symbol"),
                   values=genes, mart= mart)
hgnc <- gene_hgnc$hgnc_symbol[match(genes, gene_hgnc$ensembl_gene_id)]
hgnc <- hgnc[hgnc != ""]
common <- intersect(row.names(MYOD1.obj), hgnc)
MYOD1.obj <- subset(MYOD1.obj, features = common)
gene_ensembl <- getBM(filters= c("hgnc_symbol"), attributes= c("ensembl_gene_id","hgnc_symbol"),
                      values=common, mart= mart)
ensembl <- gene_ensembl$ensembl_gene_id[match(common, gene_ensembl$hgnc_symbol)]
muscle.obj <- subset(muscle.obj, features = ensembl)
saveRDS(muscle.obj, file = "/path/scRNA/RMS/TabulaSapiens_muscle.rds")
muscle.obj <- readRDS(file = "/path/scRNA/RMS/TabulaSapiens_muscle.rds")
genes <- row.names(muscle.obj)
gene_hgnc <- getBM(filters= c("ensembl_gene_id"), attributes= c("ensembl_gene_id","hgnc_symbol"),
                   values=genes, mart= mart)
hgnc <- gene_hgnc$hgnc_symbol[match(genes, gene_hgnc$ensembl_gene_id)]
row.names(muscle.obj) <- hgnc
MYOD1.obj <- subset(MYOD1.obj, features = hgnc)

MYOD1.obj.dat <- as.matrix(MYOD1.obj@layers$data)
colnames(MYOD1.obj.dat) <- colnames(MYOD1.obj)
#colnames(MYOD1.obj.dat) <- sub("_.*", "", colnames(MYOD1.obj.dat))
rownames(MYOD1.obj.dat) <- rownames(MYOD1.obj)

set.seed(123)
selected_columns <- sample(ncol(MYOD1.obj.dat), ncol(muscle.obj.dat))
MYOD1.obj.dat_subset <- MYOD1.obj.dat[, selected_columns]

#by subcluster
MYOD1.obj.integrated$sample_cluster <- MYOD1.integrated.vp$sample_cluster ##tbd
library(limma)
sample_cluster <- MYOD1.integrated.vp$sample_cluster
# na_indices <- which(is.na(sample_cluster))
# if(length(na_indices) > 0) {
#   sample_cluster <- sample_cluster[-na_indices]
#   MYOD1.vp.dat <- MYOD1.vp.dat[, -na_indices]
# }
combined_matrix <- cbind(MYOD1.obj.dat, muscle.obj.dat)
group <- factor(c(as.character(sample_cluster), rep("muscle", ncol(muscle.obj.dat))))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit <- lmFit(combined_matrix, design)
subclusters <- sort(unique(MYOD1.integrated.vp@meta.data$sample_cluster))
contrast_formula <- paste(subclusters, "- muscle")
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
top_genes_list <- list()
for (contrast_name in colnames(contrast_matrix)) {
  top_genes <- topTable(fit2, coef = contrast_name, adjust = "fdr", number = Inf)
  top_genes$Comparison <- contrast_name
  top_genes_list[[contrast_name]] <- top_genes
}
combined_top_genes <- do.call(rbind, top_genes_list)
write.csv(combined_top_genes, file = "MYOD1_subcluster_muscle_diffexp_limma_gene.csv")

####viper on diff exp subclusters vs muscle####
setwd("/path/scRNA/RMS/seurat/ARACNe")
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
nets <- lapply(filenames, readRDS)
names(nets) <- sample
combined_top_genes <- read.csv("MYOD1_subcluster_muscle_diffexp_limma_gene.csv") %>% column_to_rownames("X") 
combined_top_genes$Gene <- sub(".*\\.", "", rownames(combined_top_genes))
dat <- combined_top_genes %>% dplyr::select(Comparison, Gene, t) %>%
  pivot_wider(names_from = "Comparison", values_from = t) %>%
  column_to_rownames("Gene")
genes <- rownames(dat)
mart <- readRDS("/path/scRNA/functions/mart.obj.rds")
gene_ensembl <- getBM(filters = c("hgnc_symbol"), attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      values = genes, mart = mart)
ensembl <- gene_ensembl$ensembl_gene_id[match(genes, gene_ensembl$hgnc_symbol)]
row.names(dat) <- ensembl
vp <- viper(dat, nets, method = 'none')
genes <- rownames(vp)
gene_hgnc <- getBM(filters = c("ensembl_gene_id"), attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   values = genes, mart = mart)
hgnc <- gene_hgnc$hgnc_symbol[match(genes, gene_hgnc$ensembl_gene_id)]
row.names(vp) <- hgnc
write.csv(vp, file = "MYOD1_subcluster_muscle_diffexp_viper_gene.csv")

############
#subset with OncoTarget
setwd("/path/scRNA/RMS/seurat/ARACNe")
vp <- read.csv("MYOD1_subcluster_muscle_diffexp_viper_gene.csv", row.names = 1)
colnames(vp) <- gsub("\\.\\.\\.muscle", "", colnames(vp))
colnames(vp) <- gsub("_(cluster)_", ".\\1", colnames(vp))
OncoTarget <- read.csv("/path/scRNA/functions/oncotarget.csv")
vp_ot <- vp[rownames(vp) %in% OncoTarget$Target,]
write.csv(vp_ot, file = "MYOD1_subcluster_muscle_diffexp_viper_gene_oncotarget.csv")
vp_ot <- read.csv("MYOD1_subcluster_muscle_diffexp_viper_gene_oncotarget.csv", row.names = 1)

df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv")
df <- df[!grepl("MYOD1_PDX", df$subcluster), ]
cluster <- df$subcluster
group <- df$cluster
heatmap_matrix <- vp_ot
heatmap_matrix <- heatmap_matrix[, match(cluster, colnames(heatmap_matrix))]

row_variances <- apply(heatmap_matrix, 1, var, na.rm = TRUE)
top20_indices <- order(row_variances, decreasing = TRUE)[1:20]
heatmap_matrix <- heatmap_matrix[top20_indices, ]
heatmap_matrix

annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster
anno.colors <- list(Group = c("1" = "#08519c", "2" = "#e31a1c", "3" = "#1a9850", "4" = "#6a3d9a"))
annotation_df <- data.frame(Group = group)

group_order <- annotation$Group[match(colnames(heatmap_matrix), annotation$Subcluster)]
gaps_col <- which(diff(as.numeric(factor(group_order))) != 0)
paletteLength <- 60
myColor <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(paletteLength)
myBreaks <- c(seq(min(vp_ot), 0, length.out = ceiling(paletteLength / 2) + 1), 
              seq(max(vp_ot) / paletteLength, max(vp_ot), length.out = floor(paletteLength / 2)))

pdf("MYOD1_subcluster_muscle_diffexp_viper_gene_oncotarget_heatmap.v2.top20.pdf", height = 5, width = 5)
pheatmap(heatmap_matrix, 
         cluster_rows = TRUE, 
         treeheight_row = 0, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 8, 
         fontsize_col = 6,
         color = myColor,
         breaks = myBreaks,
         annotation_col = annotation_df,
         annotation_colors = anno.colors,
         gaps_col = gaps_col,
         main = "")
dev.off()
############
###Oncotarget#####
# convert NES scores to pvalues ( < 10 ^-5 converted to log scale)
# OncoTarget analysis for patients, PDX, and cell lines
# Subset viper matrix for 180 regulators that have known FDA-approved or late-stage experimental drugs directly targeting them
setwd("/path/scRNA/RMS/seurat/ARACNe")
# vp <- read.csv("MYOD1_subcluster_muscle_diffexp_viper.csv", row.names = 1)
# OncoTarget <- read.csv("/path/scRNA/functions/oncotarget.csv")
# vp_ot <- vp[rownames(vp) %in% OncoTarget$Target,]
# write.csv(vp_ot, file = "MYOD1_subcluster_muscle_diffexp_viper_oncotarget.csv")
vp_ot <- read.csv("MYOD1_subcluster_muscle_diffexp_viper_gene_oncotarget.csv", row.names = 1)

vpmat.oncotarget.pval <- apply(vp_ot, 2, function(x) -log10(p.adjust(pnorm(x, lower.tail=FALSE),method='bonferroni')))
write.csv(vpmat.oncotarget.pval, "MYOD1_subcluster_muscle_diffexp_viper_oncotarget_log10pval.csv")

tmp10 <- vpmat.oncotarget.pval >= 1 # at a fdr-corrected p-value of < 10E-1
pos <- rowSums(tmp10) >= 1 #can be 1, significant in at least one sample
tmp <- filterRowMatrix(vpmat.oncotarget.pval,pos)

tmp <- tmp[order(rowMeans(tmp), decreasing = TRUE), ]
colnames(tmp) <- gsub("\\.\\.\\.muscle", "", colnames(tmp))
colnames(tmp) <- gsub("_(cluster)_", ".\\1", colnames(tmp))
# df <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_viper_similarity.silinfo.csv",
#                row.names = 1)
# cluster <- rownames(df)
# group <- df$cluster
df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv")
df <- df[!grepl("MYOD1_PDX", df$subcluster), ]
cluster <- df$subcluster
group <- df$cluster
tmp <- tmp[, match(cluster, colnames(tmp))]

annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster
anno.colors <- list(Group = c("1" = "#08519c", "2" = "#e31a1c", "3" = "#1a9850", "4" = "#6a3d9a"))

group_order <- annotation$Group[match(colnames(tmp), annotation$Subcluster)]
gaps_col <- which(diff(as.numeric(factor(group_order))) != 0)

paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(tmp, na.rm = TRUE), 1, length.out = ceiling(paletteLength / 2) + 1), 
              seq(1 + 1e-6, max(tmp, na.rm = TRUE), length.out = floor(paletteLength / 2)))

pdf("MYOD1_subcluster_muscle_diffexp_viper_gene_oncotarget_heatmap_Significant10Percent.v2.pdf", height = 5, width = 7)
pheatmap(tmp, 
         cluster_rows = TRUE, 
         treeheight_row = 0, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 8, 
         fontsize_col = 6,
         color = myColor,
         breaks = myBreaks,
         annotation_col = annotation,
         annotation_colors = anno.colors,
         gaps_col = gaps_col,
         main = "")
dev.off()
############

###heatmap###
rows_to_remove <- grep("^RPS|^RPL", rownames(vp), value = TRUE)
vp <- vp[!rownames(vp) %in% rows_to_remove, ]
gene_variances <- apply(vp, 1, var)
top50_genes <- names(sort(gene_variances, decreasing = TRUE))[1:50]
top50_genes <- names(sort(gene_variances, decreasing = TRUE))[1:40]
top50_vp <- vp[top50_genes, ]
#muscle genes
prefixes <- c("MYH", "MYL", "MYO", "PAX3", "PAX7", "SMOC", "CALD", "PDGFR", "ACT", "DES", "TNN")
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
muscle_genes <- grep(pattern, rownames(vp), value = TRUE)
top50_vp <- vp[muscle_genes, ]

dat <- as.matrix(top50_vp)
# Create the heatmap
dat <- dat[order(rowMeans(dat), decreasing = TRUE), ]
colnames(dat) <- gsub("\\.\\.\\.muscle", "", colnames(dat))
colnames(dat) <- gsub("_(cluster)_", ".\\1", colnames(dat))

df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv")
df <- df[!grepl("MYOD1_PDX", df$subcluster), ]
cluster <- df$subcluster
group <- df$cluster
dat <- dat[, match(cluster, colnames(dat))]

annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster
anno.colors <- list(Group = c("1" = "#08519c", "2" = "#e31a1c", "3" = "#1a9850", "4" = "#6a3d9a"))

group_order <- annotation$Group[match(colnames(dat), annotation$Subcluster)]
gaps_col <- which(diff(as.numeric(factor(group_order))) != 0)

paletteLength <- 60
myColor <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(paletteLength)
myBreaks <- c(seq(min(dat), 0, length.out = ceiling(paletteLength / 2) + 1), 
              seq(max(dat) / paletteLength, max(dat), length.out = floor(paletteLength / 2)))
pdf("MYOD1_subclusters_vp_genes_muscle_diffexp_heatmap.v4.pdf", height = 8, width = 8)
pdf("MYOD1_subclusters_vp_musclegenes_muscle_diffexp_heatmap.v3.pdf", height = 5, width = 8)
pheatmap(dat, 
         cluster_rows = TRUE, 
         treeheight_row = 0, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 8, 
         fontsize_col = 6,
         color = myColor,  
         breaks = myBreaks,
         annotation_col = annotation,
         annotation_colors = anno.colors,
         gaps_col = gaps_col,
         main = "")
dev.off()

##########################################
###bulk vs normal muscle sample by sample####
setwd("/path/scRNA/RMS")
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
nets <- lapply(filenames, readRDS)
names(nets) <- sample

#normal muscle (Tabula Sapiens)
muscle.obj <- readRDS(file = "/path/scRNA/RMS/9b34eb5c-1fb2-46c5-9fe4-d90a2e5f7c4b.rds")
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

write.csv(vpmat, file = "MYOD1_bulk_vs_muscle_diffexp_viper.csv", row.names = TRUE)

###heatmap###
vpmat <- read.csv("MYOD1_bulk_vs_muscle_diffexp_viper.csv", row.names = 1)
df <- read.csv("MYOD1_bulk.csv") %>% arrange(Type, Treatment)
df$Type <- factor(df$Type, levels = c("Biopsy", "Resection", "Local relapse", "Metastasis", "Unknown"))

vpmat <- vpmat[, match(df$ID, colnames(vpmat))]

annotation <- data.frame(Sample = df$Type, Treatment = df$Treatment)
rownames(annotation) <- df$ID
anno.colors <- list(Sample = c("Biopsy" = "#66c2a4", "Resection" = "#006d2c", "Local relapse" = "#e41a1c", "Metastasis" = "#984ea3",
                               "Unknown" = "#d9d9d9"),
                    Treatment = c("Non-treated" = "#0571b0", "Post-treatment" = "#8c510a", "Unknown" = "#d9d9d9"))

# group_order <- annotation$Sample[match(colnames(vpmat), annotation$ID)]
# gaps_col <- which(diff(as.numeric(factor(group_order))) != 0)

rows_to_remove <- grep("^RPS|^RPL", rownames(vpmat), value = TRUE)
vpmat <- vpmat[!rownames(vpmat) %in% rows_to_remove, ]
gene_variances <- apply(vpmat, 1, var)
top50_genes <- names(sort(gene_variances, decreasing = TRUE))[1:40]
top50_vp <- vpmat[top50_genes, ]
top_dat <- as.matrix(top50_vp)
# Create the heatmap
paletteLength <- 60
myColor <- colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(paletteLength)
myBreaks <- c(seq(min(top_dat), 0, length.out = ceiling(paletteLength / 2) + 1), 
              seq(max(top_dat) / paletteLength, max(top_dat), length.out = floor(paletteLength / 2)))
pdf("MYOD1_bulk_vs_muscle_diffexp_genes_viper_heatmap.v2.pdf", height = 8, width = 8)
pheatmap(top_dat, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 6, 
         fontsize_col = 6,
         color = myColor,
         breaks = myBreaks,
         annotation_col = annotation,
         annotation_colors = anno.colors,
         #  gaps_col = gaps_col,
         main = "")
dev.off()
###########
###Oncotarget#####
# convert NES scores to pvalues ( < 10 ^-5 converted to log scale)
# OncoTarget analysis for patients, PDX, and cell lines
# Subset viper matrix for 180 regulators that have known FDA-approved or late-stage experimental drugs directly targeting them
# setwd("/path/scRNA/RMS/seurat/ARACNe")
setwd("/path/scRNA/RMS/MYOD1_bulk")
vp <- read.csv("MYOD1_bulk_vs_muscle_diffexp_viper.csv", row.names = 1)
OncoTarget <- read.csv("/path/scRNA/functions/oncotarget.csv")
vp_ot <- vp[rownames(vp) %in% OncoTarget$Target,]
write.csv(vp_ot, file = "MYOD1_bulk_muscle_diffexp_viper_oncotarget.csv")

vpmat.oncotarget.pval <- apply(vp_ot, 2, function(x) -log10(p.adjust(pnorm(x, lower.tail=FALSE),method='bonferroni')))
write.csv(vpmat.oncotarget.pval, "MYOD1_bulk_muscle_diffexp_viper_oncotarget_log10pval.csv")

# PDF of a clustered heatmap including all targetable proteins that were significant in at least one sample 
# at a fdr-corrected p-value of < 10E-5
vpmat <- vpmat.oncotarget.pval >= 1 # at a fdr-corrected p-value of < 10E-1
pos <- rowSums(vpmat) >= 1 #can be 1, significant in at least one sample
vpmat <- filterRowMatrix(vpmat.oncotarget.pval,pos)
vpmat <- vpmat[order(rowMeans(vpmat), decreasing = TRUE), ]

df <- read.csv("MYOD1_bulk.csv") %>% arrange(Type, Treatment)
df$Type <- factor(df$Type, levels = c("Biopsy", "Resection", "Local relapse", "Metastasis", "Unknown"))

vpmat <- vpmat[, match(df$ID, colnames(vpmat))]

annotation <- data.frame(Sample = df$Type, Treatment = df$Treatment)
rownames(annotation) <- df$ID
anno.colors <- list(Sample = c("Biopsy" = "#66c2a4", "Resection" = "#006d2c", "Local relapse" = "#e41a1c", "Metastasis" = "#984ea3",
                               "Unknown" = "#d9d9d9"),
                    Treatment = c("Non-treated" = "#0571b0", "Post-treatment" = "#8c510a", "Unknown" = "#d9d9d9"))

paletteLength <- 100
myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
myBreaks <- c(seq(min(vpmat, na.rm=TRUE), 1, length.out = ceiling(paletteLength / 2) + 1), 
              seq(1 + 1e-6, max(vpmat, na.rm=TRUE), length.out = floor(paletteLength / 2)))

pdf("MYOD1_bulk_muscle_diffexp_viper_oncotarget_heatmap_Significant10Percent.v3.pdf", height = 4, width = 6)
pheatmap(vpmat, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         fontsize_row = 8, 
         fontsize_col = 6,
         annotation_col = annotation,
         annotation_colors = anno.colors,
         color = myColor,
         breaks = myBreaks)
dev.off()

###########
