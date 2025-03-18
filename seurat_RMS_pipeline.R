setwd("/path/scRNA")
library(scales)
library(tidyverse)
library(viper)
library(Seurat)
library(MAST)
library(igraph)
library(biomaRt)
library(uwot)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(glmGamPoi)
library(PISCES)
library(celldex)
library(cellxgene.census)
library(SingleCellExperiment)
library(scater)
library(SingleR)
library(scran)
library(scuttle)
library(infercnv)
library(ggsci)
library(HGNChelper)
library(openxlsx)
library(CytoTRACE)
library(DESeq2)
library(readxl)
library(presto)
library(CytoTRACE)
library(clusterProfiler)
library(msigdbr)
# library(fpc)
library(patchwork)
library(CellChat)
source("functions/process-utils.R")
source("functions/cluster-functions.R")
source("functions/viper-utils.R")
source("functions/misc.R")

#create integrated Seurat object
########MYOD1 integrated#########
setwd("/path/scRNA/RMS")
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
#names(sample) <- sample
obj.list <- list()
for (i in seq_along(sample)) {
  print(i)
  exp.mat <- Read10X(sample[i])
  obj <- CreateSeuratObject(counts = exp.mat, min.cells = 3, min.features = 200)
  obj <- SCTransform(obj, verbose = T, conserve.memory = T) 
  obj.list <- append(obj.list, obj)
}
setwd("/path/scRNA/RMS/seurat")
saveRDS(obj.list, file = "MYOD1.obj.list.rds")

####seurat v4 intergration#####
obj.list <- readRDS(file = "MYOD1.obj.list.rds")
features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features, verbose = T)
anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = features, verbose = T)
rm(obj.list,features)
MYOD1.obj.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = T)
rm(anchors)
saveRDS(MYOD1.obj.integrated, file = "MYOD1.obj.integrated.rds")
#find clusters
MYOD1.obj.integrated <- RunPCA(MYOD1.obj.integrated, features = VariableFeatures(object = MYOD1.obj.integrated))
pdf(paste0('MYOD1.obj.elbow.pdf'))
ElbowPlot(MYOD1.obj.integrated)
dev.off()
MYOD1.obj.integrated <- RunUMAP(MYOD1.obj.integrated, dims = 1:50, verbose = FALSE, metric="correlation")
MYOD1.obj.integrated <- FindNeighbors(MYOD1.obj.integrated, dims = 1:50, verbose = FALSE)
MYOD1.obj.integrated <- FindClusters(MYOD1.obj.integrated, resolution=seq(0.01,1,by=0.01), verbose = FALSE,algorithm=1) 
saveRDS(MYOD1.obj.integrated, file = "MYOD1.obj.integrated.rds")
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.integrated.rds")
###################################
#######SingleR######
setwd("/path/scRNA/RMS/seurat")
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.integrated.rds")
df <- as.SingleCellExperiment(MYOD1.obj.integrated, assay = "SCT")
df <- logNormCounts(df)
#SingleR with Blueprint
ref.bp.data <- BlueprintEncodeData(ensembl=FALSE)
singleR.bp.res <- SingleR(test = df, ref = ref.bp.data, assay.type.test = 1, labels = ref.bp.data$label.main)
saveRDS(singleR.bp.res, file = "MYOD1.obj.integrated.blueprint_ref_singleR.rds")
table(singleR.bp.res$labels)
write.table(table(singleR.bp.res$labels), paste0('singleR_bp.labels.txt'), row.names = F)
MYOD1.obj.integrated[["SingleR.labels"]] <- singleR.res$labels
setwd("/path/scRNA/RMS/seurat", height = 10, width = 20)
pdf(paste0('singleR_blueprint_plots.pdf'))
plotScoreHeatmap(singleR.bp.res)
dev.off()
##############################

######ARACNE networks from samples######
setwd("/path/scRNA/RMS")
project <- "RMS"
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475","MYOD1_PDX_PRMS")
sample <- "MYOD1_RMS_3_2"
sample <- "MYOD1_RMS_31_2"
sample <- "MYOD1_RMS_41B"
sample <- "MYOD1_RMS_211_2"
sample <- "MYOD1_RMS_469"
sample <- "MYOD1_RMS_475"
sample <- "MYOD1_PDX_PRMS"
ref.data <- celldex::BlueprintEncodeData(ensembl=FALSE)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#ref_group_names <- c("Endothelial cells","Macrophages","Monocytes","B-cells","CD4+ T-cells","CD8+ T-cells","NK cells","Eosinophils","DC","HSC")
#RMS
ref_group_names <- c("Endothelial cells","Macrophages","Monocytes","B-cells","CD4+ T-cells","CD8+ T-cells","NK cells","Eosinophils","Erythrocytes","HSC","Mesangial cells","DC","Melanocytes")
ref_group_names <- c("Endothelial cells","Macrophages","Monocytes","CD4+ T-cells","CD8+ T-cells","Epithelial cells","Neutrophils","Mesangial cells","Melanocytes")
ref_group_names <- c("Endothelial cells","Macrophages","Monocytes","B-cells","CD4+ T-cells","CD8+ T-cells","DC","Erythrocytes","Neutrophils")
ref_group_names <- c("Endothelial cells","Macrophages","Monocytes","B-cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Epithelial cells","HSC","Erythrocytes","Neutrophils","Keratinocytes")
ref_group_names <- c("Endothelial cells","Macrophages","Monocytes","B-cells","CD4+ T-cells","CD8+ T-cells","DC","Erythrocytes","HSC","Melanocytes")
ref_group_names <- c("Endothelial cells","Macrophages","Monocytes","B-cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Erythrocytes","HSC","Mesangial cells","Melanocytes")
ref_group_names <- c("Endothelial cells","Macrophages","Monocytes","B-cells","CD4+ T-cells","CD8+ T-cells","DC","Eosinophils","Erythrocytes","HSC","Melanocytes","Mesangial cells")
setwd("/path/scRNA/RMS")
## load data from 10x, create seurat object
exp.mat <- Read10X(sample)
pisces.obj <- CreateSeuratObject(counts = exp.mat, min.cells = 3, min.features = 200)
#counts_matrix = GetAssayData(pisces.obj, layer = "counts")
genes_to_remove <- grep("\\.[1-9]|^LINC|-", rownames(pisces.obj), value = TRUE) #remove noncoding genes
genes_to_keep <- rownames(pisces.obj)[!(rownames(pisces.obj) %in% genes_to_remove)]
pisces.obj <- subset(pisces.obj, features = genes_to_keep)
#singleR
df <- as.SingleCellExperiment(pisces.obj, assay = "RNA")
df <- logNormCounts(df)
singleR.res <- SingleR(test = df, ref = ref.data, assay.type.test = 1, labels = ref.data$label.main)
saveRDS(singleR.res, file = paste0(sample, '/singleR.res.rds'))
table(singleR.res$labels)
write.table(table(singleR.res$labels), paste0(sample, '/singleR_labels.txt'), row.names = F)
pisces.obj[["SingleR.labels"]] <- singleR.res$labels
Idents(pisces.obj) <- singleR.res$labels
pisces.obj.subset <- subset(pisces.obj, 
                            idents = ref_group_names,
                            invert = TRUE)
saveRDS(pisces.obj.subset, file = paste0(sample, '_subset.rds'))
#saveRDS(pisces.obj.subset, file = paste0(sample, '_subset_fullfeatures.rds')) #without filtering noncoding genes
pisces.obj.subset <- readRDS(file = paste0(sample, '_subset.rds'))
counts_matrix <- GetAssayData(pisces.obj.subset, layer = "counts")
 #convert to ARACNe input
genes <- row.names(counts_matrix)
gene_ensembl <- getBM(filters= c("hgnc_symbol"), attributes= c("ensembl_gene_id","hgnc_symbol"),
                    values=genes, mart= mart)
ensembl <- gene_ensembl$ensembl_gene_id[match(genes, gene_ensembl$hgnc_symbol)]
row.names(counts_matrix) <- ensembl
cpm.mat <- CPMTransform(counts_matrix)
cpm.mat <- cpm.mat[, sample(colnames(cpm.mat), min(ncol(cpm.mat), 500)) ]
setwd("/path/scRNA/RMS/seurat/ARACNe")
saveRDS(cpm.mat, file = paste0(sample, '_cpm.mat.rds'))
for(i in seq_along(sample)){
  cpm.mat <- readRDS(file = paste0(sample[i], '_cpm.mat.rds'))
  m <- as.data.frame(cpm.mat)
  m <- cbind(gene = rownames(m), m)
  colnames(m)[1] <- "gene"
  write.table(m, paste0(sample[i],'_cpm.tsv'), sep = "\t", row.names = F, col.names = T, quote = F)
  m2 <- m %>% dplyr::select(-gene)
  saveRDS(m2, file = paste0(sample[i], '_cpm.mat.rds'))
}
#run ARACNe (run_bootstraps.sh, network_final.sh) on HPC
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475","MYOD1_PDX_PRMS")
for (i in sample) {
  setwd("/path/scRNA/RMS/seurat/ARACNe")
  print(i)
  cpm.mat <- as.matrix(readRDS(paste0(i,'_cpm.mat.rds')))
  net_final <- paste0(i,"/network/network_final.txt")
  RegProcess(net_final, cpm.mat, out.dir = paste0(i,"/network/"), out.name = paste0(i,'-net-'))
}
###############################################

##########seurat v5 integration#########
#filter by integrating subset objects
setwd("/path/scRNA/RMS")
#obj.list <- readRDS(file = "MYOD1.obj.list.rds")
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
MYOD1.subset.obj.list <- list()
for (i in seq_along(sample)) {
  print(i)
  obj <- readRDS(file = paste0(sample[i], '_subset.rds'))
  obj <- SCTransform(obj, verbose = T, conserve.memory = T)
  MYOD1.subset.obj.list <- append(MYOD1.subset.obj.list, obj)
}
setwd("/path/scRNA/RMS/seurat/ARACNe")
saveRDS(MYOD1.subset.obj.list, file = "MYOD1.subset.obj.list.rds")
MYOD1.subset.obj.list <- readRDS(file = "MYOD1.subset.obj.list.rds")
features <- SelectIntegrationFeatures(object.list = MYOD1.subset.obj.list)
MYOD1.subset.obj.list <- PrepSCTIntegration(object.list = MYOD1.subset.obj.list, anchor.features = features, verbose = T)
anchors <- FindIntegrationAnchors(object.list = MYOD1.subset.obj.list, normalization.method = "SCT", anchor.features = features, verbose = T)
saveRDS(anchors, file = "MYOD1.filt.anchors.rds")
MYOD1.obj.filt.integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = T)
saveRDS(MYOD1.obj.filt.integrated, file = "MYOD1.obj.filt.integrated.rds")
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.filt.integrated.rds") #integrated from individual subset objects

MYOD1.obj.integrated <- RunPCA(MYOD1.obj.integrated, features = VariableFeatures(object = MYOD1.obj.integrated))
pdf(paste0('MYOD1.obj.filtered.elbow.v2.pdf'))
ElbowPlot(MYOD1.obj.integrated, ndims = 50, reduction = "pca") 
ElbowPlot(MYOD1.obj.integrated, ndims = 30, reduction = "pca")
dev.off()
MYOD1.obj.integrated <- RunUMAP(MYOD1.obj.integrated, dims = 1:10, verbose = T)
MYOD1.obj.integrated <- FindNeighbors(MYOD1.obj.integrated, dims = 1:10)
MYOD1.obj.integrated <- FindClusters(MYOD1.obj.integrated, resolution = seq(0.01,1,by=0.01), verbose = TRUE, algorithm = 1) 
sample_name <- c(rep(sample[[1]], ncol(MYOD1.subset.obj.list[[1]])),
                 rep(sample[[2]], ncol(MYOD1.subset.obj.list[[2]])),
                 rep(sample[[3]], ncol(MYOD1.subset.obj.list[[3]])),
                 rep(sample[[4]], ncol(MYOD1.subset.obj.list[[4]])),
                 rep(sample[[5]], ncol(MYOD1.subset.obj.list[[5]])),
                 rep(sample[[6]], ncol(MYOD1.subset.obj.list[[6]]))) 
MYOD1.obj.integrated$sample <- sample_name
saveRDS(MYOD1.obj.integrated, file = "MYOD1.obj.filt.integrated.v3.rds") #integrated from individual subset objects
#determine optimal clusters
clust=MYOD1.obj.integrated@meta.data[,which(grepl("integrated_snn_res.",colnames(MYOD1.obj.integrated@meta.data)))]
mat=MYOD1.obj.integrated@assays$SCT@scale.data
out=sil_subsample(mat,clust)
means=out[[1]]
sd=out[[2]]
x=seq(0.01,1,by=0.01)
means=means[1:100]
pdf(paste0('MYOD1.obj.filt.integrated.silhouette.v3.pdf'))
  plot(x, means, type = "n", ylab="mean silhouette score", xlab="resolution parameter")
  errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
  lines(x,means)
  best=tail(x[which(means[5:length(means)]==max(means[5:length(means)]))+4],n=1)
  legend("topright",paste("Best",best,sep = " = "))
dev.off()
MYOD1.obj.integrated$seurat_clusters=MYOD1.obj.integrated@meta.data[,which(colnames(MYOD1.obj.integrated@meta.data)==paste("integrated_snn_res.",best,sep=""))]
MYOD1.obj.integrated <- FindClusters(MYOD1.obj.integrated, resolution = 0.05)
Idents(MYOD1.obj.integrated) <- "seurat_clusters"
saveRDS(MYOD1.obj.integrated, file = "MYOD1.obj.filt.integrated.v2.rds") #integrated from individual subset objects

#create UMAP plot
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.filt.integrated.v2.rds")
# Stash cell identity classes in metadata
MYOD1.obj.integrated[["seurat_clusters"]] <- Idents(object = MYOD1.obj.integrated)
pdf(paste0('MYOD1.obj.filt.integrated.umap.v2.pdf'), height = 10, width = 20)
plot(DimPlot(MYOD1.obj.integrated, reduction = "umap",
     group.by = c("sample", "seurat_clusters"),
     label=TRUE,repel=T,label.size=5))
dev.off()
#perform differential expression
MYOD1.obj.integrated <- PrepSCTFindMarkers(MYOD1.obj.integrated, verbose = T)
markers <- FindAllMarkers(MYOD1.obj.integrated, only.pos = TRUE)
write.csv(markers, file = "MYOD1.obj.filt.integrated.allmarker.v2.csv")
markers <- read.csv("MYOD1.obj.filt.integrated.allmarker.v2.csv")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
pdf(paste0('MYOD1.obj.filt.integrated.marker.v2.pdf'))
DoHeatmap(MYOD1.obj.integrated, features = top10$gene, group.by = "ident", slot = "data", label = TRUE) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
DoHeatmap(MYOD1.obj.integrated, features = top10$gene, group.by = "ident", slot = "data", label = TRUE) + NoLegend() 
dev.off()
setwd("/path/scRNA")
source("functions/misc3.R")
setwd("/path/scRNA/RMS/seurat/ARACNe")
pdf(paste0('MYOD1.obj.filt.integrated.allmarker.geneHeatmap.v2.pdf'))
 geneHeatmap_plot(MYOD1.obj.integrated,MYOD1.obj.integrated$seurat_clusters,top10$gene,genes_by_cluster = T,n_top_genes_per_cluster = 20)
dev.off()
saveRDS(MYOD1.obj.integrated, file = "MYOD1.obj.filt.integrated.v2.rds")
###############################################
#####sample vs cluster####
metadata <- MYOD1.obj.integrated@meta.data
contingency_table <- table(metadata$seurat_clusters, metadata$sample)
contingency_df <- as.data.frame(contingency_table)
colnames(contingency_df) <- c("Cluster", "Sample", "Frequency")
pdf(paste0('MYOD1.obj.filt.integrated.cluster.sample.barplot.pdf'), height = 5, width = 10)
  ggplot(contingency_df, aes(x = factor(Sample, 
                  levels = c("MYOD1_RMS_475","MYOD1_RMS_41B","MYOD1_RMS_211_2","MYOD1_RMS_3_2","MYOD1_RMS_31_2","MYOD1_RMS_469")), 
          y = Frequency, fill = Cluster)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Distribution of Samples in Each Seurat Cluster",
        x = "Seurat Cluster",
        y = "Frequency") + 
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    theme_minimal()
  ggplot(contingency_df, aes(x = Cluster, y = Frequency, fill = Sample)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Distribution of Samples in Each Seurat Cluster",
        x = "Seurat Cluster",
        y = "Frequency") + 
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    theme_minimal()
dev.off()
###############################################
###########match gene and protein response categories#########
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.filt.integrated.v2.rds")
MYOD1.obj.integrated[["SingleR.labels"]] <- singleR.res$labels
#matching elements from MYOD1.obj.integrated$sample to MYOD1.integrated.vp
matched_elements <- unname(MYOD1.obj.integrated$sample[colnames(MYOD1.integrated.vp)])
MYOD1.integrated.vp$sample <- matched_elements
matched_elements <- unname(MYOD1.obj.integrated$sample[colnames(MYOD1.integrated.vp)])
new_matched_elements <- recode(matched_elements, "MYOD1_RMS_3_2" = "recurrence", #disease recurrent viable progressed on chemo
                                                 "MYOD1_RMS_31_2" = "progressed", #disease progression little treatment effect
                                                 "MYOD1_RMS_41B" = "stable", #stable disease (no change in size MRI) NED
                                                 "MYOD1_RMS_211_2" = "recurrence", #modest response to salvage chemo. relapsed quickly
                                                 "MYOD1_RMS_469" = "progressed", #progression despite 1-2 cycles of chemo
                                                 "MYOD1_RMS_475" = "stable") #stable disease (no change in size MRI) NED
new_matched_elements <- recode(matched_elements, "MYOD1_RMS_3_2" = "resistant", 
                                                 "MYOD1_RMS_31_2" = "resistant", 
                                                 "MYOD1_RMS_41B" = "responsive", 
                                                 "MYOD1_RMS_211_2" = "resistant", 
                                                 "MYOD1_RMS_469" = "resistant", 
                                                 "MYOD1_RMS_475" = "responsive")
MYOD1.obj.integrated$response <- new_matched_elements                           
MYOD1.integrated.vp$response <- new_matched_elements
##############################################
###CytoTRACE###
###MYOD1.obj.integrated#######
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.obj.integrated <- readRDS("MYOD1.obj.filt.integrated.v2.rds")
df <- as.SingleCellExperiment(MYOD1.obj.integrated, assay = "SCT")
counts_matrix <- as.matrix(df@assays@data@listData$counts)
MYOD1.obj.integrated.cytotrace <- CytoTRACE(counts_matrix, ncores = 8)
#CytoTRACE ouput is a list object containing:
#numeric values for CytoTRACE (values ranging from 0 (more differentiated) to 1 (less differentiated))
#ranked CytoTRACE, GCS, and gene counts
#numeric vector of the Pearson correlation between each gene and CytoTRACE
#numeric vector of the Pearson correlation between each gene and gene counts
#the IDs of filtered cells, and a normalized gene expression table
saveRDS(MYOD1.obj.integrated.cytotrace, file = "MYOD1.obj.filt.integrated.cytotrace.v2.rds")
MYOD1.obj.integrated.cytotrace <- readRDS("MYOD1.obj.filt.integrated.cytotrace.v2.rds")
MYOD1.obj.integrated <- AddMetaData(MYOD1.obj.integrated, metadata = MYOD1.obj.integrated.cytotrace$CytoTRACE, col.name = "CytoTRACE")
saveRDS(MYOD1.obj.integrated, file = "MYOD1.obj.filt.integrated.v2.rds")
pdf(paste0('MYOD1.obj.filt.integrated.cytotrace.umap.v2.pdf'), height = 10, width = 10)
  FeaturePlot(MYOD1.obj.integrated, features = c("CytoTRACE")) +
    scale_color_gradientn(colors = rev(brewer.pal(11, "RdYlBu")))
dev.off()

###############################################

####interactome#######
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.obj.integrated <- readRDS("MYOD1.obj.filt.integrated.v2.rds")
clusters <- unique(Idents(MYOD1.obj.integrated))
sample <- unique(MYOD1.obj.integrated$sample)
filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
nets=lapply(filenames,readRDS)
names(nets) <- sample
# restructure the interactome into a data frame containing each regulator-target pair
# Extract the data
# Initialize empty lists to store results
all_regulators <- list()
all_targets <- list()
all_likelihoods <- list()
all_tfmodes <- list()
all_clusters <- list()
for (i in seq_along(nets)) {
  nets.sub <- nets[[i]] 
  regulators <- unlist(lapply(names(nets.sub), function(regulator) {
    rep(regulator, times = length(nets.sub[[regulator]]$tfmode))
  }))  
  targets <- unlist(lapply(nets.sub, function(regulon) {
    names(regulon$tfmode)
  })) 
  likelihoods <- unlist(lapply(nets.sub, function(regulon) {
    regulon$likelihood
  })) 
  tfmodes <- unlist(lapply(nets.sub, function(regulon) {
    regulon$tfmode
  })) 
  clusters <- rep(names(nets)[i], times = length(regulators))
  # Append results to the lists
  all_regulators <- c(all_regulators, regulators)
  all_targets <- c(all_targets, targets)
  all_likelihoods <- c(all_likelihoods, likelihoods)
  all_tfmodes <- c(all_tfmodes, tfmodes)
  all_clusters <- c(all_clusters, clusters)
}
interactome_df <- data.frame(
  regulator = unlist(all_regulators),
  target = unlist(all_targets),
  likelihood = unlist(all_likelihoods),
  tfmode = unlist(all_tfmodes),
  cluster = unlist(all_clusters),
  stringsAsFactors = FALSE
)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
mart <- readRDS("/path/scRNA/functions/mart.obj.rds")
regulator_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = interactome_df$regulator,
  mart = mart
)
target_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = interactome_df$target,
  mart = mart
)
interactome_df$regulator <- regulator_symbols$hgnc_symbol[match(interactome_df$regulator, regulator_symbols$ensembl_gene_id)]
interactome_df$target <- target_symbols$hgnc_symbol[match(interactome_df$target, target_symbols$ensembl_gene_id)]
interactome_df <- interactome_df[!duplicated(interactome_df), ]
write.csv(interactome_df, "MYOD1.integrated.vp.network_sample.interactome.csv", row.names = FALSE)
##########################################

######run VIPER######
##load aracne networks
#run VIPER
setwd("/path/scRNA/RMS/seurat/ARACNe")
####run viper on each sample using networks from all 6 samples###
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
setwd("/path/scRNA/RMS")
sample.obj <- lapply(sample, function(s) readRDS(file = paste0(s, '_subset.rds')))
filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
nets <- lapply(filenames, readRDS)
names(nets) <- sample
all_vp_RMS <- list()
for (s in seq_along(sample)) {
  sample.obj[[s]] <- SCTransform(sample.obj[[s]], verbose = T, conserve.memory = T)
  dat <- as.matrix(sample.obj[[s]]@assays$SCT@data)
  mart <- readRDS("/data/vanderbilt/dermawaj/scRNA/functions/mart.obj.rds")
  genes <- row.names(dat)
  gene_ensembl <- getBM(filters = c("hgnc_symbol"), attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        values = genes, mart = mart)
  ensembl <- gene_ensembl$ensembl_gene_id[match(genes, gene_ensembl$hgnc_symbol)]
  row.names(dat) <- ensembl
  dat <- dat[!is.na(rownames(dat)), ]
  colnames(dat) <- paste0(colnames(dat), "_", s)
  indices <- seq(1, ncol(dat), by = 400)
  seurat_viper_list <- list()
  for (i in 1:(length(indices) - 1)) {
    vp <- viper(dat[, indices[i]:indices[i + 1]], nets, method = 'scale')
    seurat_viper_list <- c(seurat_viper_list, list(vp))
  }
  vp <- viper(dat[, indices[length(indices)]:ncol(dat)], nets, method = 'scale')
  seurat_viper_list <- c(seurat_viper_list, list(vp))
  vp_RMS <- seurat_viper_list[[1]]
  for (i in 2:length(seurat_viper_list)) {
    vp_RMS <- cbind(vp_RMS, seurat_viper_list[[i]])
  }
  all_vp_RMS[[sample[s]]] <- vp_RMS
  setwd("/data/vanderbilt/dermawaj/scRNA/RMS/seurat/ARACNe")
  saveRDS(vp_RMS, paste0(sample[s], '.viper.sample_network.rds'))
}
combined_vp_RMS <- do.call(cbind, all_vp_RMS)
setwd("/path/scRNA/RMS/seurat/ARACNe")
saveRDS(combined_vp_RMS, "MYOD1.filt.obj.filtered.viper.sample_network.indiv.sample.rds")
#############################

##########################################
#####VIPER clustering of each sample and then integrate signature######
setwd("/path/scRNA/RMS/seurat/ARACNe")
vp_RMS <- readRDS("MYOD1.filt.obj.filtered.viper.sample_network.indiv.sample.v2.rds")
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.filt.integrated.v2.rds")
common_cell <- intersect(colnames(vp_RMS), colnames(MYOD1.obj.integrated))
vp_RMS <- vp_RMS[, common_cell]
mart <- readRDS("/path/scRNA/functions/mart.obj.rds")
genes <- row.names(vp_RMS)
gene_hgnc <- getBM(filters= c("ensembl_gene_id"), attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values=genes, mart= mart)
hgnc <- gene_hgnc$hgnc_symbol[match(genes, gene_hgnc$ensembl_gene_id)]
row.names(vp_RMS) <- hgnc
vp_RMS <- vp_RMS[!is.na(rownames(vp_RMS)),]
cbcMRs <- CBCMRs(vp_RMS)
vp_RMS <- vp_RMS[!duplicated(rownames(vp_RMS)), ]
counts <- vp_RMS[cbcMRs,]
options(Seurat.object.assay.version = 'v3')
MYOD1.integrated.vp <- CreateSeuratObject(counts = counts, project = 'MYOD1')
                     
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
MYOD1.integrated.vp <- ScaleData(MYOD1.integrated.vp, assay = "RNA")
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.filt.integrated.v2.rds")
MYOD1.integrated.vp$sample <- MYOD1.obj.integrated$sample
MYOD1.integrated.vp$response <- MYOD1.obj.integrated$response
MYOD1.integrated.vp$seurat_clusters_gene <- MYOD1.obj.integrated$seurat_clusters
MYOD1.integrated.vp$CytoTRACE <- MYOD1.obj.integrated$CytoTRACE
unique_samples <- unique(MYOD1.integrated.vp$sample)
MYOD1.integrated.vp.sample.list <- list()
for (i in seq_along(unique_samples)) {
  sample_name <- unique_samples[i]
  MYOD1.integrated.vp.sample <- subset(MYOD1.integrated.vp, subset = sample == sample_name)
  MYOD1.integrated.vp.sample <- RunPCA(MYOD1.integrated.vp.sample, features = rownames(MYOD1.integrated.vp.sample))
  MYOD1.integrated.vp.sample <- RunUMAP(MYOD1.integrated.vp.sample, dims = 1:5, verbose = T)
  MYOD1.integrated.vp.sample <- FindNeighbors(MYOD1.integrated.vp.sample, dims = 1:5)
  MYOD1.integrated.vp.sample <- FindClusters(MYOD1.integrated.vp.sample, resolution = seq(0.01,1,by=0.01), verbose = TRUE, algorithm = 1) 
  clust=MYOD1.integrated.vp.sample@meta.data[,which(grepl("RNA_snn_res.",colnames(MYOD1.integrated.vp.sample@meta.data)))]
  mat=MYOD1.integrated.vp.sample@assays$RNA@scale.data
  out=sil_subsample(mat,clust)
  means=out[[1]]
  sd=out[[2]]
  x=seq(0.01,1,by=0.01)
  errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
  lines(x,means)
  best=tail(x[which(means[5:length(means)]==max(means[5:length(means)]))+4],n=1)
  Idents(MYOD1.integrated.vp.sample) <- "seurat_clusters"
  MYOD1.integrated.vp.sample.list <- append(MYOD1.integrated.vp.sample.list, MYOD1.integrated.vp.sample)   
}
saveRDS(MYOD1.integrated.vp.sample.list, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.list.rds")
MYOD1.integrated.vp.sample.list <- readRDS("MYOD1.integrated.filt.vp.network_sample.indiv.sample.list.rds")
#get all object metadata: MYOD1.integrated.vp.sample.list[[1]][[]]
#get list of metadata columns: colnames(MYOD1.integrated.vp.sample.list[[1]][[]])
##########################################
####cluster based on viperSimilarity matrix####
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.v2.rds")
MYOD1.integrated.vp.sample.list <- readRDS("MYOD1.integrated.filt.vp.network_sample.indiv.sample.list.rds")
# Define the sample list
sample_list <- unique(MYOD1.integrated.vp$sample)
names(MYOD1.integrated.vp.sample.list) <- sample_list
#sample_list <- c("MYOD1_RMS_3_2", "MYOD1_RMS_211_2", "MYOD1_RMS_31_2", "MYOD1_RMS_475", "MYOD1_RMS_41B", "MYOD1_RMS_469")
# Define the groups based on viper similarity heatmap
cell_names_list <- list()
for (i in seq_along(MYOD1.integrated.vp.sample.list)) {
  seurat_obj <- MYOD1.integrated.vp.sample.list[[i]]
  sample_name <- names(MYOD1.integrated.vp.sample.list)[i]
  for (cluster_id in levels(seurat_obj$seurat_clusters)) {
    cluster_cells <- WhichCells(seurat_obj, idents = cluster_id)
    cell_names_list[[paste0(sample_name, "_cluster_", cluster_id)]] <- cluster_cells
  }
}
MYOD1.integrated.vp$sample_cluster <- NA
for (key in names(cell_names_list)) {
  MYOD1.integrated.vp$sample_cluster[cell_names_list[[key]]] <- key
}
table(MYOD1.integrated.vp$sample_cluster)

###################################
results_list <- list()
unique_samples <- unique(MYOD1.integrated.vp$sample)
for (sample_name in seq_along(unique_samples)) {
  seurat_obj <- MYOD1.integrated.vp.sample.list[[sample_name]]
  # Get the sample name from the metadata
  sample_metadata <- unique(seurat_obj$sample)
  sample_name <- sample_metadata[1]
  # Initialize a list to store results for each cluster in the current Seurat object
  sample_results <- list()
  # Iterate over clusters 
  for (cluster_id in levels(seurat_obj$seurat_clusters)) {
    # Subset the count matrix for cells in the current cluster
    cluster_cells <- which(seurat_obj$seurat_clusters == cluster_id)
    count_matrix_cluster <- seurat_obj@assays$RNA@data[, cluster_cells]
    # count_matrix_cluster <- seurat_obj@assays$RNA@scale.data[, cluster_cells]
    # Apply the Stouffer's method function to each gene (row) in the subsetted count matrix
    # weights <- rep(1, ncol(count_matrix_cluster))  # Example: equal weights for all cells
    # cluster_results <- apply(count_matrix_cluster, 1, stouffersMethod, weights = weights)
    cluster_results <- rowMeans(count_matrix_cluster)
    # Store the results for the current cluster with appropriate column name
    sample_results[[paste0("cluster", cluster_id)]] <- cluster_results
  }
  results_list[[sample_name]] <- sample_results
}
combined_results <- do.call(cbind, unlist(results_list, recursive = FALSE))
combined_matrix <- as.matrix(combined_results) 
combined_matrix <- combined_matrix[, colSums(is.na(combined_matrix)) == 0]
write.csv(combined_matrix, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.mean.csv")

                     
#run viper on combined matrix to generate distance matrix--
#computes the similarity between VIPER signatures across samples
#Viper similarity enrichment of the top regulators in one signature in the differential protein activity and integrates the two enrichment scores)

dd <- viperSimilarity(combined_matrix, nn = 50, method = c("two.sided"))

#optimal clustering
library(fpc)
clust <- pamk(data=dd, krange=2:10)$pamobject 
clust$silinfo #silhouette score > 0.25 is considered good (can be used as weight for stouffer integration)
#generate heatmap ordered by the silhouette score
df <- as.data.frame(clust$silinfo$widths)
write.csv(df, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_viper_similarity.silinfo.csv")
cluster <- rownames(clust$silinfo$widths)
group <- df$cluster

dd <- dd[match(cluster, rownames(dd)), match(cluster, colnames(dd))]
annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster

anno.colors <- list(Group = c("differentiated" = "#08519c", "progenitor" = "#e31a1c", "intermediate" = "#6a3d9a"))

paletteLength <- 60
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)
myBreaks <- c(seq(min(dd), 0, length.out = ceiling(paletteLength/2) + 1), 
              seq(max(dd)/paletteLength, max(dd), length.out = floor(paletteLength/2)))
pdf("MYOD1.integrated.filt.vp.network_sample.individualclusters.indiv.sample.combined_matrix_viper_similarity_heatmap.top50.pamk.pdf", height = 8, width = 8)
pheatmap(dd, cluster_rows = F, cluster_cols = F, 
         show_rownames = TRUE, show_colnames = TRUE, 
         fontsize_row = 6, fontsize_col = 6,
         color = myColor,
         breaks = myBreaks,
         annotation_col = annotation,
         annotation_colors = anno.colors,
         main = "Viper Similarity Matrix") 
dev.off()
##########################################

differentiated <- c("MYOD1_RMS_475_cluster_3",
                    "MYOD1_RMS_211_2_cluster_2",
                    "MYOD1_RMS_41B_cluster_1",
                    "MYOD1_RMS_3_2_cluster_1",
                    "MYOD1_RMS_31_2_cluster_1",
                    "MYOD1_RMS_469_cluster_0")
progenitor <- c("MYOD1_RMS_41B_cluster_2",
                "MYOD1_RMS_469_cluster_1",
                "MYOD1_RMS_41B_cluster_4",
                "MYOD1_RMS_3_2_cluster_0",
                "MYOD1_RMS_31_2_cluster_0",
                "MYOD1_RMS_211_2_cluster_1",
                "MYOD1_RMS_475_cluster_1",
                "MYOD1_RMS_475_cluster_2")
intermediate <- c("MYOD1_RMS_475_cluster_0",
            "MYOD1_RMS_211_2_cluster_0",
            "MYOD1_RMS_41B_cluster_0",
            "MYOD1_RMS_469_cluster_2",
            "MYOD1_RMS_31_2_cluster_2",
            "MYOD1_RMS_41B_cluster_3")

table(MYOD1.integrated.vp$sample_cluster)

MYOD1.integrated.vp$cell_state <- NA
MYOD1.integrated.vp$cell_state[MYOD1.integrated.vp$sample_cluster %in% differentiated] <- "differentiated"
MYOD1.integrated.vp$cell_state[MYOD1.integrated.vp$sample_cluster %in% progenitor] <- "progenitor"
# MYOD1.integrated.vp$cell_state[MYOD1.integrated.vp$sample_cluster %in% group3] <- "group3"
MYOD1.integrated.vp$cell_state[MYOD1.integrated.vp$sample_cluster %in% intermediate] <- "intermediate"
table(MYOD1.integrated.vp$cell_state)

saveRDS(MYOD1.integrated.vp, "MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
##########################################
##############################
#####visualize top MRs in each group###
setwd("/data/vanderbilt/dermawaj/scRNA/RMS/seurat/ARACNe")
MYOD1.obj.integrated.cytotrace <- readRDS(file = "MYOD1.obj.filt.integrated.cytotrace.v2.rds")
cytoGenes <- MYOD1.obj.integrated.cytotrace$cytoGenes
df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv") %>% arrange(cluster)
cluster <- df$subcluster
group <- df$cluster

annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster
anno.colors <- list(Group = c("differentiated" = "#08519c", "progenitor" = "#e31a1c", "intermediate" = "#6a3d9a"))

combined_matrix <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.mean.csv", row.names = 1)

combined_matrix <- combined_matrix[, match(cluster, colnames(combined_matrix))]
group_order <- annotation$Group[match(colnames(combined_matrix), annotation$Subcluster)]
gaps_col <- which(diff(as.numeric(factor(group_order))) != 0)

combined_matrix <- combined_matrix[!grepl("RPL|RPS", rownames(combined_matrix)), ]

stouffersMethod <- function(x, weights) {
  return(sum(x * weights) / sqrt(sum(weights * weights)))
}
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
(selected_genes <- get_top_genes(integrated_matrix, n = 15))
selected_genes <- selected_genes[!selected_genes %in% c("BTG1","SDF4","CRYAB","ZNF556","DGKZ","GRIA2","CENPF")]

filtered_matrix <- integrated_matrix[selected_genes, ]
rownames(filtered_matrix)

library(ComplexHeatmap)
library(circlize)
diff_matrix <- filtered_matrix
diff_matrix <- cbind(diff_matrix, CytoTRACE = cytoGenes[rownames(diff_matrix)])
diff_matrix <- na.omit(diff_matrix)

heatmap_matrix <- diff_matrix[, -ncol(diff_matrix)]
cytotrace_data <- diff_matrix[, "CytoTRACE", drop = FALSE]
paletteLength <- 60 
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)

heatmap1 <- Heatmap(
  as.matrix(heatmap_matrix),
  name = "Activity",
  col = myColor,
  # column_split = sort(group),
  # row_split = c(rep("progenitor", 21), rep("myogenic", 21)),
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = F,
  width = unit(3, "cm"),
  row_names_gp = gpar(fontsize = 6),
  show_column_dend = F,
  column_names_gp = gpar(fontsize = 5) 
)
heatmap2 <- Heatmap(
  as.matrix(cytotrace_data),
  name = "CytoTRACE",
  col = colorRamp2(c(min(cytotrace_data), 0, max(cytotrace_data)), c("blue", "white", "red")),
  show_row_names = TRUE,
  show_column_names = TRUE,
  # row_split = c(rep("progenitor", 21), rep("myogenic", 21)),
  cluster_rows = T,
  cluster_columns = T,
  show_row_dend = F,
  width = unit(0.5, "cm")
)
combined_heatmap <- heatmap2 + heatmap1
pdf(paste0('MYOD1.integrated.filt.vp.network_sample.integrated_matrix_topMR.cytotrace.pdf'), height = 5, width = 3)
    draw(combined_heatmap)
dev.off()  
####################################################################
#####pathway enrichment######
# list2regulon1Tail
# @param x List containing the genes or proteins to be converted in a regulon object.
# @param only_names If the list is  contains only names, i.e. there are no scores
# associated to genes, set this paramter as TRUE. Defeault is FALSE
# @export
# @author Pasquale Laise
#'

list2regulon1Tail <- function(x, only_names = FALSE) {
  res <- lapply(x,function(x){
    tfmode <- rep(1, length(x))
    if(only_names==F){names(tfmode) <- names(x)}
    else{names(tfmode) <- x}
    list(tfmode=tfmode, likelihood=rep(1, length(tfmode)))
    })
    class(res) <- "regulon"
    res
}

####Lineage markers Marker list enrichment####
#sctype
db_path <- "/path/scRNA/functions/ScTypeDB_db.xlsx"
db_data <- read_excel(db_path)
cell_type_data <- db_data %>%
    filter(tissueType %in% c("Muscle","Embryo")) %>%
  dplyr::select(Cell_Type = 2, Gene_List = 3) 
gene_lists <- cell_type_data %>%
  mutate(Gene_List = strsplit(Gene_List, ",")) %>%
  deframe()

regul <- list2regulon1Tail(gene_lists,only_names = T)

combined_matrix <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.mean.csv", row.names = 1)
stouffersMethod <- function(x, weights) {
  return(sum(x * weights) / sqrt(sum(weights * weights)))
}
df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv")
cluster <- df$subcluster
group <- df$cluster
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

path_enrich <- aREA(combined_matrix, regul, minsize = 2)
path_enrich <- aREA(integrated_matrix, regul, minsize = 2) # sign_tab is the matrix of signatures (matrix of mean viper signatures for each of the subclusters)
pathways.nes <- path_enrich$nes

pathways.nes.pval <- apply(pathways.nes, 2, function(x) p.adjust(pnorm(abs(x), lower.tail=FALSE), method='fdr'))
significant <- apply(pathways.nes.pval, 1, function(x) sum(x <= 0.05))
pathways.nes.significant <- pathways.nes[which(significant > 0),]                        
pathways.nes.significant
heatmap_matrix <- pathways.nes.significant

df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv")
cluster <- df$subcluster
group <- df$cluster
heatmap_matrix <- pathways.nes.significant[, match(cluster, colnames(pathways.nes.significant))]

annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster
anno.colors <- list(Group = c("differentiated" = "#08519c", "progenitor" = "#e31a1c", "intermediate" = "#6a3d9a"))

group_order <- annotation$Group[match(colnames(heatmap_matrix), annotation$Subcluster)]
gaps_col <- which(diff(as.numeric(factor(group_order))) != 0)

paletteLength <- 60 
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)
myBreaks <- c(seq(min(heatmap_matrix), 0, length.out = ceiling(paletteLength/2) + 1), 

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
          annotation_col = annotation,
          annotation_colors = anno.colors,
          gaps_col = gaps_col,
          main = "")
##########################################
###CytoTRACE###
###MYOD1.integrated.vp#######
setwd("/data/vanderbilt/dermawaj/scRNA/RMS/seurat/ARACNe")
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.v2.rds")
MYOD1.integrated.vp.cytotrace <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.cytotrace.rds")
MYOD1.obj.integrated.cytotrace <- readRDS("MYOD1.obj.filt.integrated.cytotrace.v2.rds")
MYOD1.obj.integrated <- readRDS("MYOD1.obj.filt.integrated.v2.rds")

##cytotrace sample by sample###
MYOD1.subset.obj.list <- readRDS(file = "MYOD1.subset.obj.list.rds")
for (i in seq_along(MYOD1.subset.obj.list)) {
  seurat_obj <- MYOD1.subset.obj.list[[i]]
  seurat_obj <- SCTransform(seurat_obj, verbose = T)
  counts_matrix <- as.matrix(seurat_obj@assays$SCT@counts)
  seurat_obj_cytotrace <- CytoTRACE(counts_matrix, ncores = 8)
  cytotrace_score <- seurat_obj_cytotrace$CytoTRACE
  new_colnames <- paste0(colnames(seurat_obj), "_", i) 
  # Add cytotrace score to MYOD1.integrated.vp metadata
  MYOD1.integrated.vp@meta.data[new_colnames, "CytoTRACE.Score"] <- cytotrace_score
}

##cytotrace by subclusters###
metadata <- MYOD1.integrated.vp@meta.data
pdf(paste0('MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.viper_similarity.cytotrace.subcluster.pdf'), height = 5, width = 10)
  ggplot(metadata, aes(x = sample_cluster, y = CytoTRACE, fill = sample_cluster)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
       scale_fill_igv() +
    labs(title = "",
        x = "Subclusters",
        y = "CytoTRACE Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
dev.off()

##cytotrace by cell states###
Idents(MYOD1.integrated.vp) <- MYOD1.integrated.vp$cell_state
metadata <- MYOD1.integrated.vp@meta.data

group_colors <- c("differentiated" = "#08519c", "progenitor" = "#e31a1c", "intermediate" = "#6a3d9a")

pairwise_comparisons <- pairwise.t.test(metadata$CytoTRACE, metadata$group, p.adjust.method = "bonferroni")
comparisons1 <- list(
  c("differentiated", "intermediate"),
  c("differentiated", "progenitor")
)

##cytotrace by subclusters###
pdf(paste0('MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.viper_similarity.cytotrace.cell_state.v4.pdf'), height = 5, width = 5)
  ggplot(metadata, aes(x = cell_state, y = CytoTRACE, fill = cell_state)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
       # scale_fill_brewer(palette = "Dark2") ++
    scale_fill_manual(values = group_colors) +
    labs(title = "",
        x = "Cell States",
        y = "CytoTRACE Score") +
    stat_compare_means(comparisons = comparisons1, method = "t.test", p.adjust.method = "bonferroni", label = "p.signif") +
    theme_minimal() #+
dev.off()
              
##########################################
###MYOD1 bulk##########
#Normalize the bulk gene expression by log2CPM 
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
bulkmatrix <- CPM_normalization(dat)

#viper
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
filenames <- list.files("/path/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
nets <- lapply(filenames, readRDS)
names(nets) <- sample

vpmat <- viper(bulkmatrix, nets, method="scale") # internally scaled viper signature

mart <- readRDS("/path/mart.obj.rds")
genes <- row.names(vpmat)
gene_hgnc <- getBM(filters= c("ensembl_gene_id"), attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values=genes, mart= mart)
hgnc <- gene_hgnc$hgnc_symbol[match(genes, gene_hgnc$ensembl_gene_id)]
row.names(vpmat) <- hgnc

write.csv(vpmat, file = "MYOD1_bulk_viper.csv")

#viper similarity of bulk samples#
vpmat <- read.csv("MYOD1_bulk_viper.csv", row.names = 1)
vpsim <- viperSimilarity(vpmat, nn = 100, method = c("two.sided")) # can try all or setting nn=50 or 100 #distance matrix
vpmat.cluster <- pamk(data=vpsim, krange=2:5)$pamobject
vpmat.cluster$silinfo # see what average silhouette scores you get with optimal k or different k’s
#generate heatmap ordered by the silhouette score
df <- as.data.frame(vpmat.cluster$silinfo$widths)
write.csv(df, "MYOD1_bulk_viperSimilarity_cluster.csv")
cluster <- rownames(vpmat.cluster$silinfo$widths)
group <- df$cluster

#top 100
vpsim <- vpsim[match(cluster, rownames(vpsim)), match(cluster, colnames(vpsim))]
annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster
anno.colors <- list(Group = c("2" = "#08519c", "3" = "#e31a1c", "1" = "#1a9850"))
#cols = colorRampPalette(c("blue", "white", "red"))(100)
paletteLength <- 60
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)
myBreaks <- c(seq(min(vpsim), 0, length.out = ceiling(paletteLength/2) + 1), 
              seq(max(vpsim)/paletteLength, max(vpsim), length.out = floor(paletteLength/2)))

pheatmap(vpsim, cluster_rows = FALSE, cluster_cols = FALSE, 
         show_rownames = TRUE, show_colnames = TRUE, 
         fontsize_row = 6, fontsize_col = 6,
         color = myColor,
         breaks = myBreaks,
         main = "Viper Similarity Matrix",
         annotation_col = annotation,
         annotation_colors = anno.colors) 

#compare with integrated sc matrix###
combined_matrix <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.mean.csv", row.names = 1)
df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv") %>% arrange(cluster)
cluster <- df$subcluster
group <- df$cluster
combined_matrix <- combined_matrix[, match(cluster, colnames(combined_matrix))]
combined_matrix <- combined_matrix[!grepl("RPL|RPS", rownames(combined_matrix)), ]

stouffersMethod <- function(x, weights) {
  return(sum(x * weights) / sqrt(sum(weights * weights)))
}
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
integrated_sc_matrix <- do.call(cbind, results_list)

vpmat <- read.csv("/data/vanderbilt/dermawaj/scRNA/RMS/MYOD1_bulk/MYOD1_bulk_viper.csv", row.names = 1)
vpsim <- viperSimilarity(vpmat, nn = 100, method = c("two.sided"))
stouffersMethod <- function(x, weights) {
  return(sum(x * weights) / sqrt(sum(weights * weights)))
}
vpmat.cluster <- pamk(data=vpsim, krange=2:5)$pamobject
vpmat.cluster$silinfo # see what average silhouette scores you get with optimal k or different k’s
#generate heatmap ordered by the silhouette score
df <- as.data.frame(vpmat.cluster$silinfo$widths)
cluster <- rownames(vpmat.cluster$silinfo$widths)

group <- df$cluster
vpmat <- vpmat[, match(cluster, colnames(vpmat))]
grouped_clusters <- split(colnames(vpmat), group)
results_list <- list()
for (group_name in names(grouped_clusters)) {
  subcluster_names <- grouped_clusters[[group_name]]
  subcluster_matrix <- vpmat[, subcluster_names, drop = FALSE]
  weights <- rep(1, ncol(subcluster_matrix))
  integrated_results <- apply(subcluster_matrix, 1, stouffersMethod, weights = weights)
  results_list[[group_name]] <- integrated_results
}
integrated_bulk_matrix <- do.call(cbind, results_list)

common_genes <- intersect(rownames(integrated_sc_matrix), rownames(integrated_bulk_matrix))
integrated_sc_matrix <- integrated_sc_matrix[common_genes, ]
integrated_bulk_matrix <- integrated_bulk_matrix[common_genes, ]
colnames(integrated_bulk_matrix) <- paste0("MYOD1_bulk_", colnames(integrated_bulk_matrix))
combined_matrix <- cbind(integrated_sc_matrix, integrated_bulk_matrix)
combined_matrix <- combined_matrix[!grepl("RPL|RPS", rownames(combined_matrix)), ]

vpsim <- viperSimilarity(combined_matrix, nn = NULL, method = c("two.sided")) # can try all or setting nn=50 or 100 #distance matrix
vpmat.cluster <- pamk(data=vpsim, krange=2:5)$pamobject
vpmat.cluster$silinfo 

paletteLength <- 60
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)
myBreaks <- c(seq(min(vpsim), 0, length.out = ceiling(paletteLength/2) + 1), 
              seq(max(vpsim)/paletteLength, max(vpsim), length.out = floor(paletteLength/2)))

pheatmap(vpsim, cluster_rows = T, cluster_cols = T, 
         show_rownames = TRUE, show_colnames = TRUE, 
         fontsize_row = 6, fontsize_col = 6,
         color = myColor,
         breaks = myBreaks,
         main = "Viper Similarity Matrix") 

##############################
###Cell State Signature Enrichment###
combined_matrix <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.mean.csv", row.names = 1)
df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv") %>% arrange(cluster)
cluster <- df$subcluster
group <- df$cluster
combined_matrix <- combined_matrix[, match(cluster, colnames(combined_matrix))]
combined_matrix <- combined_matrix[!grepl("RPL|RPS", rownames(combined_matrix)), ]

stouffersMethod <- function(x, weights) {
  return(sum(x * weights) / sqrt(sum(weights * weights)))
}
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
integrated_sc_matrix <- do.call(cbind, results_list)

#create regulon from top 50 and bottom 50 MRs from sc viper matrix integrated by cell state
get_top_genes <- function(vector, n = 50) {
  top_genes <- names(sort(vector, decreasing = TRUE)[seq_len(min(n, length(vector)))])
  # bottom_genes <- names(sort(vector, decreasing = FALSE)[seq_len(min(n, length(vector)))])
  return(top_genes)
}
gene_lists <- list()
cell_states <- c("differentiated", "intermediate", "progenitor")
for (state in cell_states) {
  state_matrix <- integrated_sc_matrix[, grepl(state, colnames(integrated_sc_matrix))]
  gene_lists[[state]] <- get_top_genes(state_matrix, n = 50)
}

list2regulon1Tail <- function(x, only_names = FALSE) {
  res <- lapply(x,function(x){
    tfmode <- rep(1, length(x))
    if(only_names==F){names(tfmode) <- names(x)}
    else{names(tfmode) <- x}
    list(tfmode=tfmode, likelihood=rep(1, length(tfmode)))
    })
    class(res) <- "regulon"
    res
}

regul <- list2regulon1Tail(gene_lists, only_names = T)

integrated_matrix <- read.csv("MYOD1_bulk.integrated.viper.csv", row.names = 1)
vpmat <- read.csv("/path/RMS/MYOD1_bulk/MYOD1_bulk_viper.csv", row.names = 1)
common_genes <- intersect(rownames(integrated_sc_matrix), rownames(vpmat))
integrated_matrix <- integrated_matrix[common_genes, ]
vpmat <- vpmat[common_genes, ]

path_enrich <- aREA(vpmat, regul, minsize = 2) # sign_tab is the matrix of signatures (matrix of mean viper signatures for each of the subclusters)
pathways.nes <- path_enrich$nes
heatmap_matrix <- pathways.nes

df2 <- read.csv("MYOD1_bulk_viperSimilarity_cluster.csv")
setwd("/path/RMS/MYOD1_bulk")
df <- read.csv("MYOD1_bulk.csv") %>% 
       left_join(df2, by = c("ID" = "X")) %>%
       arrange(cluster)
df$Type <- factor(df$Type, levels = c("Biopsy", "Resection", "Local relapse", "Metastasis", "Unknown"))
cluster <- df$ID
group <- df$cluster

heatmap_matrix <- heatmap_matrix[, match(df$ID, colnames(heatmap_matrix))]

annotation <- data.frame(Sample = df$Type, Treatment = df$Treatment, Group = df$cluster, Subcluster = df$ID)
rownames(annotation) <- df$ID
anno.colors <- list(Sample = c("Biopsy" = "#66c2a4", "Resection" = "#006d2c", "Local relapse" = "#e41a1c", "Metastasis" = "#984ea3",
                               "Unknown" = "#d9d9d9"),
                    Treatment = c("Non-treated" = "#0571b0", "Post-treatment" = "#8c510a", "Unknown" = "#d9d9d9"),
                    Group = c("1" = "#1a9850", "2" = "#08519c", "3" = "#e31a1c"))

group_order <- annotation$Group[match(colnames(heatmap_matrix), annotation$Subcluster)]
gaps_col <- which(diff(as.numeric(factor(group_order))) != 0)

paletteLength <- 60 
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)
myBreaks <- c(seq(min(heatmap_matrix), 0, length.out = ceiling(paletteLength/2) + 1), 
              seq(max(heatmap_matrix)/paletteLength, max(heatmap_matrix), length.out = floor(paletteLength/2))) 

pheatmap(heatmap_matrix, 
          cluster_rows = F, 
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
############################################################################
####cell cycle scoring######
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.filt.integrated.v2.rds")
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.v2.rds")

common_cells <- intersect(colnames(MYOD1.obj.integrated), colnames(MYOD1.integrated.vp))
MYOD1.obj.integrated <- subset(MYOD1.obj.integrated, cells = common_cells)
MYOD1.integrated.vp <- subset(MYOD1.integrated.vp, cells = common_cells)
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# assign cell-cycle scoring
MYOD1.obj.integrated <- CellCycleScoring(object = MYOD1.obj.integrated, 
                        s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

MYOD1.subset.obj.list <- readRDS(file = "MYOD1.subset.obj.list.rds")
for (i in seq_along(MYOD1.subset.obj.list)) {
  seurat_obj <- MYOD1.subset.obj.list[[i]]
  seurat_obj <- SCTransform(seurat_obj, verbose = T)
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
  s_score <- seurat_obj@meta.data$S.Score
  g2m_score <- seurat_obj@meta.data$G2M.Score
  new_colnames <- paste0(colnames(seurat_obj), "_", i) 
  # Add S.Score and G2M.Score to MYOD1.obj.integrated metadata
  MYOD1.obj.integrated@meta.data[new_colnames, "S.Score"] <- s_score
  MYOD1.obj.integrated@meta.data[new_colnames, "G2M.Score"] <- g2m_score
}
saveRDS(MYOD1.obj.integrated, file = "MYOD1.obj.filt.integrated.v2.rds")
              
MYOD1.obj.integrated$sample_cluster <- MYOD1.integrated.vp$sample_cluster
metadata <- MYOD1.obj.integrated@meta.data

pdf(paste0('MYOD1.obj.filt.integrated.clusters.sample.cellcycle.pdf'), height = 10, width = 15)
  ggplot(metadata, aes(x = sample_cluster, y = S.Score, fill = sample_cluster)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    scale_fill_igv() +
    labs(title = "",
        x = "Cluster",
        y = "S Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggplot(metadata, aes(x = sample_cluster, y = G2M.Score, fill = sample_cluster)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    scale_fill_igv() +
    labs(title = "",
        x = "Cluster",
        y = "G2M Score") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

##########################################
