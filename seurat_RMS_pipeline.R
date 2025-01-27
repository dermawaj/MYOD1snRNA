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
MYOD1.obj.integrated <- readRDS("MYOD1.obj.filt.integrated.v2.rds")
#copy net-pruned.rds into "/path/scRNA/RMS/seurat/ARACNe/network"
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
#by sample
filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
#by cluster
#filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/clusters/network_clusters", pattern="*.rds", full.names=TRUE)
nets=lapply(filenames,readRDS)
names(nets) <- sample
#names(nets) <- clusters
###
#dat=RMS.integrated.meta@assays$integrated@scale.data
df <- as.SingleCellExperiment(MYOD1.obj.integrated, assay = "SCT")
dat <- as.matrix(df@assays@data@listData$counts)
#dat=MYOD1.obj.integrated@assays$SCT@data
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(dat)
gene_ensembl <- getBM(filters= c("hgnc_symbol"), attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values=genes, mart= mart)
ensembl <- gene_ensembl$ensembl_gene_id[match(genes, gene_ensembl$hgnc_symbol)]
row.names(dat) <- ensembl
dat <- dat[!is.na(rownames(dat)),]
###
indices=seq(1,ncol(dat),by=400)
seurat_viper_list=list()
for(i in 1:(length(indices)-1)){
  #vp=viper(dat[,indices[i]:indices[i+1]], nets, method = 'none') #set to scale
  vp=viper(dat[,indices[i]:indices[i+1]], nets, method = 'scale') #set to scale
  seurat_viper_list=c(seurat_viper_list,list(vp))
}
#vp=viper(dat[,indices[i+1]:ncol(dat)],nets,method="none")
vp=viper(dat[,indices[i+1]:ncol(dat)],nets,method="scale")
seurat_viper_list=c(seurat_viper_list,list(vp))
rm(dat)
vp_RMS=seurat_viper_list[[1]]
for(i in 2:length(seurat_viper_list)){
  vp_RMS=cbind(vp_RMS,seurat_viper_list[[i]])
}
setwd("/path/scRNA/RMS/seurat/ARACNe")
saveRDS(vp_RMS, "MYOD1.filt.obj.filtered.viper.sample_network.rds")
#########################
####run viper on each sample using networks from all 6 samples###
sample <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
setwd("/path/scRNA/RMS")
sample.obj <- lapply(sample, function(s) readRDS(file = paste0(s, '_subset.rds')))
filenames <- list.files("/path/scRNA/RMS/seurat/ARACNe/network_samples", pattern="*.rds", full.names=TRUE)
nets <- lapply(filenames, readRDS)
names(nets) <- sample
all_vp_RMS <- list()
for (s in seq_along(sample)) {
  df <- as.SingleCellExperiment(sample.obj[[s]], assay = "RNA")
  dat <- as.matrix(df@assays@data@listData$counts)
  mart <- readRDS("/path/scRNA/functions/mart.obj.rds")
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
  setwd("/path/scRNA/RMS/seurat/ARACNe")
  saveRDS(vp_RMS, paste0(sample[s], '.viper.sample_network.rds'))
}
combined_vp_RMS <- do.call(cbind, all_vp_RMS)
setwd("/path/scRNA/RMS/seurat/ARACNe")
saveRDS(combined_vp_RMS, "MYOD1.filt.obj.filtered.viper.sample_network.indiv.sample.rds")
#############################

##########################################
#####VIPER clustering of each sample and then integrate signature######
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
#MYOD1.integrated.vp@assays$RNA@scale.data=as.matrix(MYOD1.integrated.vp@assays$RNA@data)
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
  #plot(x, means, type = "n", ylab="mean silhouette score", xlab="resolution parameter")
  errbar(x,means,means+sd,means-sd,ylab="mean silhouette score",xlab="resolution parameter")
  lines(x,means)
  best=tail(x[which(means[5:length(means)]==max(means[5:length(means)]))+4],n=1)
  #MYOD1.integrated.vp.sample$seurat_clusters=MYOD1.integrated.vp.sample@meta.data[,which(colnames(MYOD1.integrated.vp.sample@meta.data)==paste("RNA_snn_res.",best,sep=""))]
  MYOD1.integrated.vp.sample <- FindClusters(MYOD1.integrated.vp.sample, resolution = best)
  #MYOD1.integrated.vp.sample <- FindClusters(MYOD1.integrated.vp.sample, resolution = 0.03)
  Idents(MYOD1.integrated.vp.sample) <- "seurat_clusters"
  MYOD1.integrated.vp.sample.list <- append(MYOD1.integrated.vp.sample.list, MYOD1.integrated.vp.sample)   
}
saveRDS(MYOD1.integrated.vp.sample.list, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.list.rds")
MYOD1.integrated.vp.sample.list <- readRDS("MYOD1.integrated.filt.vp.network_sample.indiv.sample.list.rds")
#get all object metadata: MYOD1.integrated.vp.sample.list[[1]][[]]
#get list of metadata columns: colnames(MYOD1.integrated.vp.sample.list[[1]][[]])
##########################################
####cluster based on viperSimilarity matrix####
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
MYOD1.integrated.vp.sample.list <- readRDS("MYOD1.integrated.filt.vp.network_sample.indiv.sample.list.rds")
# Define the sample list
sample_list <- unique(MYOD1.integrated.vp$sample)
#sample_list <- unique_samples[!is.na(unique_samples)]
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

group1 <- c("MYOD1_RMS_211_2_cluster_2",
            "MYOD1_RMS_31_2_cluster_0",
            "MYOD1_RMS_3_2_cluster_0",
            "MYOD1_RMS_41B_cluster_3",
            "MYOD1_RMS_469_cluster_0",
            "MYOD1_RMS_475_cluster_2")
group2 <- c("MYOD1_RMS_3_2_cluster_1",
            "MYOD1_RMS_41B_cluster_0",
            "MYOD1_RMS_469_cluster_1",
            "MYOD1_RMS_475_cluster_3",
            "MYOD1_RMS_31_2_cluster_1")
group3 <- c("MYOD1_RMS_211_2_cluster_0",
            "MYOD1_RMS_211_2_cluster_1",
            "MYOD1_RMS_475_cluster_1",
            "MYOD1_RMS_475_cluster_0")
group4 <- c("MYOD1_RMS_31_2_cluster_2",
            "MYOD1_RMS_41B_cluster_2",
            "MYOD1_RMS_469_cluster_2",
            "MYOD1_RMS_475_cluster_4",
            "MYOD1_RMS_41B_cluster_1")

MYOD1.integrated.vp$group <- NA
MYOD1.integrated.vp$group[MYOD1.integrated.vp$sample_cluster %in% group1] <- "group1"
MYOD1.integrated.vp$group[MYOD1.integrated.vp$sample_cluster %in% group2] <- "group2"
MYOD1.integrated.vp$group[MYOD1.integrated.vp$sample_cluster %in% group3] <- "group3"
MYOD1.integrated.vp$group[MYOD1.integrated.vp$sample_cluster %in% group4] <- "group4"
table(MYOD1.integrated.vp$group)
saveRDS(MYOD1.integrated.vp, "MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
##########################################
#dark2_palette <- brewer.pal(n = 8, name = "Dark2")
#igv_palette <- pal_igv("default")(20)
pdf(paste0('MYOD1.integrated.vp.sample_network.indiv.sample.viper_similarity.pca.umap.v2.pdf'), height = 10, width = 10)
  plot(DimPlot(MYOD1.integrated.vp, reduction = "umap",
      group.by = c("group"),
      label=TRUE,repel=T,label.size=5),
      cols = dark2_palette) 
  plot(DimPlot(MYOD1.integrated.vp, reduction = "umap",
      group.by = c("response"),
      label=TRUE,repel=T,label.size=5),
      cols = dark2_palette) 
  plot(DimPlot(MYOD1.integrated.vp, reduction = "umap",
      group.by = c("sample"),
      label=TRUE,repel=T,label.size=5)) 
  plot(DimPlot(MYOD1.integrated.vp, reduction = "umap",
      group.by = c("sample_cluster"),
      label=TRUE,repel=T,label.size=5)) +
      guides(color = guide_legend(override.aes = list(size=4), ncol=1))
  plot(DimPlot(MYOD1.integrated.vp, reduction = "pca",
      group.by = c("group"),
      label=TRUE,repel=T,label.size=5),
      cols = dark2_palette) 
  plot(DimPlot(MYOD1.integrated.vp, reduction = "pca",
      group.by = c("response"),
      label=TRUE,repel=T,label.size=5),
      cols = dark2_palette) 
  plot(DimPlot(MYOD1.integrated.vp, reduction = "pca",
      group.by = c("sample"),
      label=TRUE,repel=T,label.size=5)) 
  plot(DimPlot(MYOD1.integrated.vp, reduction = "pca",
      group.by = c("sample_cluster"),
      label=TRUE,repel=T,label.size=5)) +
      guides(color = guide_legend(override.aes = list(size=4), ncol=1))
dev.off()
write.csv(rownames(MYOD1.integrated.vp), "MYOD1.integrated.vp.network_sample.indiv.sample.viper_similarity.genes.csv")
##############################
###gene expression clusters vs viper groups###
groups <- MYOD1.integrated.vp$group
clusters <- MYOD1.obj.integrated$seurat_clusters
df <- data.frame(
  Cell = colnames(MYOD1.obj.integrated),
  Cluster = clusters,
  Group = groups
)
df_percent <- df %>%
  group_by(Cluster, Group) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Cluster) %>%
  mutate(percentage = Count / sum(Count) * 100)
group_colors <- c("group1" = "#08519c", "group2" = "#e31a1c", "group3" = "#1a9850", "group4" = "#6a3d9a")

pdf("MYOD1_geneclusters_v_vipergroup.barplot.pdf", width = 5, height = 5)
ggplot(df_percent, aes(x = Cluster, y = percentage, fill = Group)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  scale_fill_manual(values = group_colors) +
  labs(title = "",
       x = "Seurat Clusters",
       y = "Percentage",
       fill = "Viper Group") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
dev.off()

pdf("MYOD1_cluster_frequency_by_group.v2.pdf", width = 3, height = 6)
ggplot(df, aes(x = Group, fill = Cluster)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Proportion of Seurat Clusters within Each Group",
       x = "Group",
       y = "Proportion") +
  scale_fill_brewer(palette = "Set1")
dev.off()
##############################
#Stouffer method
stouffersMethod <- function(x, weights) {
  return(sum(x * weights) / sqrt(sum(weights * weights)))
}
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
    weights <- rep(1, ncol(count_matrix_cluster))  # Example: equal weights for all cells
    cluster_results <- apply(count_matrix_cluster, 1, stouffersMethod, weights = weights)
    # Store the results for the current cluster with appropriate column name
    sample_results[[paste0("cluster", cluster_id)]] <- cluster_results
  }
  results_list[[sample_name]] <- sample_results
}
combined_results <- do.call(cbind, unlist(results_list, recursive = FALSE))
combined_matrix <- as.matrix(combined_results) 
combined_matrix <- combined_matrix[, colSums(is.na(combined_matrix)) == 0]
write.csv(combined_matrix, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.v2.csv")
combined_matrix <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.v2.csv", row.names = 1)

#run viper on combined matrix to generate distance matrix--
#computes the similarity between VIPER signatures across samples
#Viper similarity enrichment of the top regulators in one signature in the differential protein activity and integrates the two enrichment scores)
dd <- viperSimilarity(combined_matrix, nn = NULL, method = c("two.sided"))
write.csv(dd, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_viper_similarity.v2.csv")
dd <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_viper_similarity.v2.csv", row.names = 1)
#optimal clustering
library(fpc)
clust <- pamk(data=dd, krange=2:10)$pamobject 
clust$silinfo #silhouette score > 0.25 is considered good (can be used as weight for stouffer integration)
#generate heatmap ordered by the silhouette score
df <- as.data.frame(clust$silinfo$widths)
write.csv(df, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_viper_similarity.silinfo.csv")
cluster <- rownames(clust$silinfo$widths)
group <- df$cluster

df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv")
df <- df[!grepl("MYOD1_PDX", df$subcluster), ]
cluster <- df$subcluster
group <- df$cluster
#top 100
dd <- viperSimilarity(combined_matrix, nn = 100, method = c("two.sided"))
write.csv(dd, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_viper_similarity.top100.csv")
dd <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_viper_similarity.top100.csv", row.names = 1)
dd <- dd[match(cluster, rownames(dd)), match(cluster, colnames(dd))]
annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster

anno.colors <- list(Group = c("1" = "#08519c", "2" = "#e31a1c", "3" = "#1a9850", "4" = "#6a3d9a"))

paletteLength <- 60
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)
myBreaks <- c(seq(min(dd), 0, length.out = ceiling(paletteLength/2) + 1), 
              seq(max(dd)/paletteLength, max(dd), length.out = floor(paletteLength/2)))
pdf("MYOD1.integrated.filt.vp.network_sample.individualclusters.indiv.sample.combined_matrix_viper_similarity_heatmap.top100.pamk.v2.pdf", height = 8, width = 8)
pheatmap(dd, cluster_rows = FALSE, cluster_cols = FALSE, 
         show_rownames = TRUE, show_colnames = TRUE, 
         fontsize_row = 6, fontsize_col = 6,
         color = myColor,
         breaks = myBreaks,
         main = "Viper Similarity Matrix",
         annotation_col = annotation,
         annotation_colors = anno.colors) 
dev.off()

setwd("/path/scRNA/RMS/seurat/ARACNe")
# df <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_viper_similarity.silinfo.csv",
#                row.names = 1)
df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv") %>% arrange(cluster)
df <- df[!grepl("MYOD1_PDX", df$subcluster), ]
cluster <- df$subcluster
group <- df$cluster

annotation <- data.frame(Subcluster = cluster, Group = group)
rownames(annotation) <- annotation$Subcluster
anno.colors <- list(Group = c("1" = "#08519c", "2" = "#e31a1c", "3" = "#1a9850", "4" = "#6a3d9a"))

combined_matrix <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.v2.csv", row.names = 1)

combined_matrix <- combined_matrix[, match(cluster, colnames(combined_matrix))]
group_order <- annotation$Group[match(colnames(combined_matrix), annotation$Subcluster)]
gaps_col <- which(diff(as.numeric(factor(group_order))) != 0)

combined_matrix <- combined_matrix[!grepl("RPL|RPS", rownames(combined_matrix)), ]
gene_variances <- apply(combined_matrix, 1, var)
top_highest_variance_genes <- names(sort(gene_variances, decreasing = TRUE)[1:50])
filtered_matrix <- combined_matrix[top_highest_variance_genes, ]

#Top 20 and bottom 20 MRs
row_means <- rowMeans(combined_matrix, na.rm = TRUE)
top_20_rows <- names(sort(row_means, decreasing = TRUE)[1:20])
bottom_20_rows <- names(sort(row_means, decreasing = FALSE)[1:20])
top_20_rows <- names(sort(row_means, decreasing = TRUE)[1:10])
bottom_20_rows <- names(sort(row_means, decreasing = FALSE)[1:10])
selected_rows <- c(top_20_rows, bottom_20_rows)
filtered_matrix <- combined_matrix[selected_rows, ]

#muscle genes
db_path <- "/path/scRNA/functions/ScTypeDB_db.xlsx"
db_data <- read_excel(db_path)
cell_type_data <- db_data %>%
  # filter(tissueType %in% c("Muscle","Embryo")) %>%
    filter(tissueType %in% c("Muscle")) %>%
  dplyr::select(Cell_Type = 2, Gene_List = 3) %>%
  mutate(Gene_List = strsplit(Gene_List, ","))
cells <- cell_type_data %>% unnest(Gene_List)
muscle_genes <- cells %>% 
  filter(Cell_Type %in% c("Mesenchymal stem cell",
                          "Muscle satellite cell",
                          "Myoblasts",
                          "Myocytes",
                          "Smooth muscle cells")) %>% pull(Gene_List)
muscle_genes <- unique(muscle_genes)

prefixes <- c("MYH", "MYL", "MYO", "PAX3", "PAX7", "SMOC", "CALD", "PDGFR", "ACT", "DES", "TNN")
pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
muscle_genes <- grep(pattern, rownames(combined_matrix), value = TRUE)
muscle_genes <- muscle_genes[!muscle_genes %in% c("MYO10","MYO6")]

filtered_matrix <- combined_matrix[muscle_genes, ]
filtered_matrix <- filtered_matrix[complete.cases(filtered_matrix), ]

paletteLength <- 60
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)
myBreaks <- c(seq(min(filtered_matrix), 0, length.out = ceiling(paletteLength/2) + 1), 
              seq(max(filtered_matrix)/paletteLength, max(filtered_matrix), length.out = floor(paletteLength/2)))

setwd("/path/scRNA/RMS/seurat/ARACNe")
pdf("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_muscle_genes_heatmap.v4.selected.pdf", height = 5, width = 10)
pdf("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_topMR_heatmap.v3.pdf", height = 5, width = 10)
pdf("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix_heatmap.v2.pdf", height = 10, width = 10)
pheatmap(filtered_matrix, cluster_rows = F, cluster_cols = FALSE, 
         treeheight_row = 0, 
         show_rownames = TRUE, show_colnames = TRUE, 
         fontsize_row = 6, fontsize_col = 6,
         color = myColor, breaks = myBreaks,
         annotation_col = annotation,
         annotation_colors = anno.colors,
         gaps_col = gaps_col,
         main = "")
dev.off()

###use silhouette widths as weights for stouffer integration for each subcluster
weights <- df$sil_width
stouffersMethod <- function(x, weights) {
  return(sum(x * weights) / sqrt(sum(weights * weights)))
}
results_list <- list()
unique_samples <- unique(MYOD1.integrated.vp$sample)
#unique_samples <- unique_samples[!is.na(unique_samples)]
#unique_samples <- c(unique_samples[unique_samples != "MYOD1_RMS_211_2"],"MYOD1_RMS_211_2")
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
    weights <- weights  # use weights from pamk silhouette widths
    cluster_results <- apply(count_matrix_cluster, 1, stouffersMethod, weights = weights)
    # Store the results for the current cluster with appropriate column name
    sample_results[[paste0("cluster", cluster_id)]] <- cluster_results
  }
  results_list[[sample_name]] <- sample_results
}
combined_results <- do.call(cbind, unlist(results_list, recursive = FALSE))
combined_matrix <- as.matrix(combined_results) 
combined_matrix <- combined_matrix[, colSums(is.na(combined_matrix)) == 0]
write.csv(combined_matrix, file = "MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.silweighted.csv")

##############################
####cell cycle scoring######
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.obj.integrated <- readRDS(file = "MYOD1.obj.filt.integrated.v2.rds")
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
# assign cell-cycle scoring
MYOD1.obj.integrated <- CellCycleScoring(object = MYOD1.obj.integrated, 
                        s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
S.Score <- c("DTL", "HELLS", "ATAD2", "POLA1", "BRIP1")
G2M.Score <- c("HMGB2", "NUSAP1", "TPX2" , "TOP2A" , "NUF2", "MKI67", "CENPF", "SMC4", "CBX5",
               "BUB1",  "KIF11", "KIF20B", "ECT2", "ANLN",  "CKAP5", "CENPE", "GAS2L3")
metadata <- MYOD1.obj.integrated@meta.data
group_colors <- c("group1" = "#08519c", "group2" = "#e31a1c", "group3" = "#1a9850", "group4" = "#6a3d9a")

comparisons1 <- list(
  c("0", "1"),
  c("0", "2"),
  c("0", "3")
)
comparisons2 <- list(
  c("group1", "group2"),
  c("group1", "group3"),
  c("group1", "group4")
)

# Idents(MYOD1.obj.integrated) <- MYOD1.obj.integrated$seurat_clusters
pdf(paste0('MYOD1.obj.filt.integrated.clusters.vipergroups.cellcycle.pdf'), height = 5, width = 5)
  # VlnPlot(MYOD1.obj.integrated, features = c("S.Score", "G2M.Score"), pt.size = 0.000001, alpha = 0) + NoLegend()
  # RidgePlot(MYOD1.obj.integrated, features = c("S.Score", "G2M.Score"), 
  #           layer = "SCT", slot = "scale.data", ncol = 2)
  # RidgePlot(MYOD1.obj.integrated, features = c(S.Score, G2M.Score), 
  #           layer = "SCT", slot = "scale.data", ncol = 2)
  ggplot(metadata, aes(x = seurat_clusters, y = S.Score, fill = seurat_clusters)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "",
        x = "Cluster",
        y = "S Score") +
    theme_minimal() +
    stat_compare_means(comparisons = comparisons1, method = "t.test", p.adjust.method = "bonferroni", label = "p.signif")
  ggplot(metadata, aes(x = seurat_clusters, y = G2M.Score, fill = seurat_clusters)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "",
        x = "Cluster",
        y = "G2M Score") +
    theme_minimal() +
    stat_compare_means(comparisons = comparisons1, method = "t.test", p.adjust.method = "bonferroni", label = "p.signif")
  ggplot(metadata, aes(x = group, y = S.Score, fill = group)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    labs(title = "",
        x = "Viper Groups",
        y = "S Score") +
    scale_fill_manual(values = group_colors) +
    theme_minimal() +
    stat_compare_means(comparisons = comparisons2, method = "t.test", p.adjust.method = "bonferroni", label = "p.signif")
  ggplot(metadata, aes(x = group, y = G2M.Score, fill = group)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
    labs(title = "",
        x = "Viper Groups",
        y = "G2M Score") +
    scale_fill_manual(values = group_colors) +
    theme_minimal() +
    stat_compare_means(comparisons = comparisons2, method = "t.test", p.adjust.method = "bonferroni", label = "p.signif")
dev.off()

Idents(MYOD1.obj.integrated) <- MYOD1.obj.integrated$response
MYOD1.obj.integrated<- subset(MYOD1.obj.integrated, idents = c("responsive","resistant"))
pdf(paste0('MYOD1.obj.filt.integrated.response.cellcycle.pdf'), height = 20, width = 10)
  VlnPlot(MYOD1.obj.integrated, features = c("S.Score", "G2M.Score"), pt.size = 0.000001, alpha = 0) + NoLegend()
  RidgePlot(MYOD1.obj.integrated, features = c("S.Score", "G2M.Score"), 
            layer = "SCT", slot = "scale.data", ncol = 2)
  RidgePlot(MYOD1.obj.integrated, features = c(S.Score, G2M.Score), 
            layer = "SCT", slot = "scale.data", ncol = 2)
dev.off()
##########################################
#response vs group
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
metadata <- MYOD1.integrated.vp@meta.data
#subsample bootstrapping comparing group1 vs group2 in resistant vs responsive cells
bootstrap_results <- list()
#rm(.Random.seed)
set.seed(12345)
MYOD1.integrated.vp.resistant <- subset(MYOD1.integrated.vp, subset = response == "resistant")
MYOD1.integrated.vp.responsive <- subset(MYOD1.integrated.vp, subset = response == "responsive")
for (i in 1:1000) {
  resistant_cells_id <- sample(colnames(MYOD1.integrated.vp.resistant), 500, replace = TRUE)
  responsive_cells_id <- sample(colnames(MYOD1.integrated.vp.responsive), 500, replace = TRUE)
  subsampled_df <- metadata[c(resistant_cells_id, responsive_cells_id),]
  group_percentages <- subsampled_df %>%
    group_by(group, response) %>%
    summarise(Count = n(), .groups = "drop") %>%
    # group_by(group) %>%
    mutate(percentage = Count / sum(Count)) %>%
    dplyr::select(group, response, percentage)
  bootstrap_results[[i]] <- group_percentages
}
bootstrap_df <- bind_rows(bootstrap_results, .id = "iteration")
write.csv(bootstrap_df, "MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.viper_similarity.group.response.bootstrap.v3.csv", row.names = FALSE)
bootstrap_df %>%
  group_by(response, group) %>%
  summarise(Median = median(percentage), Mean = mean(percentage), SD = sd(percentage),
            Min = min(percentage), Max = max(percentage), .groups = "drop")

split_data <- split(bootstrap_df, bootstrap_df$group)
results_list <- list()
for (group_name in names(split_data)) {
  group_data <- split_data[[group_name]]
  model <- lm(percentage ~ response, data = group_data)
  anova_result <- anova(model)
  p_value <- anova_result["response", "Pr(>F)"]
  results_list[[group_name]] <- data.frame(
    group = group_name,
    p_value = p_value
  )
}
results_df <- do.call(rbind, results_list)
results_df$adjusted_p_value <- p.adjust(results_df$p_value, method = "bonferroni")
print(results_df)

annotation_df <- results_df %>%
  mutate(label = paste0("p.adj = ", format(adjusted_p_value, digits = 3)))

group_colors <- c("group1" = "#08519c", "group2" = "#e31a1c", "group3" = "#1a9850", "group4" = "#6a3d9a")

pdf("MYOD1.integrated.vp.sample_network.indiv.sample.viper_similarity.group.response.bootstrap.v3.pdf", width = 5, height = 5)
  ggplot(bootstrap_df, aes(x = response, y = percentage, fill = group, color = group)) +
    # geom_boxplot(size = 0.5) + 
    #geom_violin() +
    geom_jitter(position = position_jitterdodge(), alpha = 0.5) +
    labs(title = "Distribution of Cells Belonging to Group1 vs Group2 in Subsampled Resistant vs Responsive Cells",
        x = "Group",
        y = "Percentage") +
    theme_minimal() +
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~ group, scales = "free") +
    # stat_compare_means(method = "wilcox.test", label = "p.format")+
    geom_text(data = annotation_df, 
      aes(x = 1.5, y = max(bootstrap_df$percentage) * 0.05, label = label), 
      inherit.aes = FALSE, vjust = -0.5, size = 3)
dev.off()

#contingency table
contingency_table1 <- table(metadata$group, metadata$response)
p_value1 <- format.pval(chisq.test(contingency_table1)$p.value, digits = 3, scientific = T)
contingency_df1 <- as.data.frame(contingency_table1)
colnames(contingency_df1) <- c("Group", "Response", "Frequency")
contingency_table2 <- table(metadata$group, metadata$sample)
chisq.test(contingency_table2)
p_value2 <- format.pval(chisq.test(contingency_table2)$p.value, digits = 3, scientific = T)
contingency_df2 <- as.data.frame(contingency_table2)
colnames(contingency_df2) <- c("Group", "Sample", "Frequency")

pdf(paste0('MYOD1.integrated.vp.sample_network.indiv.sample.viper_similarity.group.response.barplot.pdf'), height = 3, width = 5)
  ggplot(contingency_df1, aes(x = Response, 
          y = Frequency, fill = Group)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Distribution of Viper Similarity Group by Response Category",
        x = "Response",
        y = "Frequency") + 
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("Chi-square p-value:", 
            format(p_value1, digits = 3)), hjust = 1.1, vjust = 1.1)

  ggplot(contingency_df1, aes(x = Group, 
          y = Frequency, fill = Response)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Distribution of Response by Viper Similarity Group",
        x = "Viper Similarity Group",
        y = "Frequency") + 
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("Chi-square p-value:", 
            format(p_value1, digits = 3)), hjust = 1.1, vjust = 1.1)
  
  ggplot(contingency_df2, aes(x = Sample, 
          y = Frequency, fill = Group)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Distribution of Viper Similarity Group by Sample",
        x = "Sample",
        y = "Frequency") + 
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("Chi-square p-value:", 
            format(p_value2, digits = 3)), hjust = 1.1, vjust = 1.1)
  ggplot(contingency_df2, aes(x = Group, 
          y = Frequency, fill = Sample)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(title = "Distribution of Sample by Viper Similarity Group",
        x = "Viper Similarity Group",
        y = "Frequency") + 
    scale_fill_brewer(palette = "Set1") +
    coord_flip() +
    theme_minimal() +
  annotate("text", x = Inf, y = Inf, label = paste("Chi-square p-value:", 
              format(p_value2, digits = 3)), hjust = 1.1, vjust = 1.1)
dev.off()

###########################################


##############################################
###CytoTRACE###
###MYOD1.integrated.vp#######
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
MYOD1.obj.integrated.cytotrace <- readRDS("MYOD1.obj.filt.integrated.cytotrace.v2.rds")
MYOD1.obj.integrated <- readRDS("MYOD1.obj.filt.integrated.v2.rds")
#plot
Idents(MYOD1.integrated.vp) <- MYOD1.integrated.vp$group
metadata <- MYOD1.integrated.vp@meta.data

group_colors <- c("group1" = "#08519c", "group2" = "#e31a1c", "group3" = "#1a9850", "group4" = "#6a3d9a")

pairwise_comparisons <- pairwise.t.test(metadata$CytoTRACE, metadata$group, p.adjust.method = "bonferroni")
comparisons1 <- list(
  c("group1", "group2"),
  c("group1", "group3"),
  c("group1", "group4")
)
pdf(paste0('MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.viper_similarity.cytotrace.group.pdf'), height = 5, width = 5)
  ggplot(metadata, aes(x = group, y = CytoTRACE, fill = group)) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5) +
       # scale_fill_brewer(palette = "Dark2") ++
    scale_fill_manual(values = group_colors) +
    labs(title = "Distribution of CytoTRACE Scores by VIPER Group",
        x = "VIPER Group",
        y = "CytoTRACE Score") +
    stat_compare_means(comparisons = comparisons1, method = "t.test", p.adjust.method = "bonferroni", label = "p.signif") +
    theme_minimal() #+
dev.off()
##############################################
######
combined_matrix <- read.csv("MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.combined_matrix.v2.csv", row.names = 1)
combined_matrix <- read.csv("/path/scRNA/RMS/MYOD1_PDX_PRMS/seurat/MYOD1.integrated.filt.vp.network_sample.indiv.sample.individualclusters.PDX.combined_matrix.csv", row.names = 1)

setwd("/path/scRNA/RMS/seurat/ARACNe")
markers.vp <- read.csv("MYOD1.integrated.filt.vp.network_sample.groups.allmarker.csv")
top10 <- markers.vp %>% group_by(cluster) %>% 
         filter(avg_log2FC > 1) %>% top_n(n = 25, wt = avg_log2FC) %>% 
         ungroup() %>% distinct(gene, .keep_all = TRUE)
library(ComplexHeatmap)
library(circlize)
diff_matrix <- combined_matrix[top10$gene, ]
diff_matrix$CytoTRACE <- cytoGenes[rownames(diff_matrix)]
diff_matrix <- na.omit(diff_matrix)
top10 <- top10 %>% filter(gene %in% rownames(diff_matrix))
df <- read.csv("MYOD1.integrated.filt.vp.subclusters.group.csv") %>% arrange(cluster)
df <- df[!grepl("MYOD1_PDX", df$subcluster), ]
cluster <- df$subcluster
group <- df$cluster
heatmap_matrix <- diff_matrix[, match(cluster, colnames(diff_matrix))]
cytotrace_data <- diff_matrix[, "CytoTRACE", drop = FALSE]
paletteLength <- 60 
myColor <- colorRampPalette((rev(brewer.pal(n = 8, name = "RdBu"))))(paletteLength)
heatmap1 <- Heatmap(
  as.matrix(heatmap_matrix),
  name = "Expression",
  col = myColor,
  column_split = group,
  row_split = top10$cluster,
  show_row_names = TRUE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  width = unit(8, "cm"),
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 4) 
)
heatmap2 <- Heatmap(
  as.matrix(cytotrace_data),
  name = "CytoTRACE",
  col = colorRamp2(c(min(cytotrace_data), 0, max(cytotrace_data)), c("blue", "white", "red")),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_split = top10$cluster,
  cluster_rows = T,
  cluster_columns = FALSE,
  width = unit(0.5, "cm")
)
combined_heatmap <- heatmap2 + heatmap1
##############################
