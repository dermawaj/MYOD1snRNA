



#############################
####Cases swimmers plot####
cases_swimmers_plot <- read.csv("cases_swimmers_plot.csv")

cases_swimmers_plot %>% 
  ggplot(aes(y = id, x = month)) +
  geom_line(data = cases_swimmers_plot %>% filter(treatment != ""),
            aes(x = month, col = factor(treatment)), linewidth = 1.8) +
  #  geom_bar(stat = "identity", aes(fill = factor(treatment)), width = 0.5) + 
  #  scale_color_brewer(palette = "Dark2") +
  geom_point(data = cases_swimmers_plot %>% filter(status!=""),
             aes(shape = fct_relevel(status, "Diagnosis","PR","CR","POD","Surgery")),
             size = 3, col = "black") +
  geom_segment(data = cases_swimmers_plot %>% filter(continued == 1),
               aes(y = id, yend = id, x = month - 2, xend = month + 1),
               pch = 15, size = 0.8, arrow = arrow(type = "closed", length = unit(0.1, "in"))) +
  geom_segment(data = cases_swimmers_plot %>% filter(continued == 2),
               aes(y = id, yend = id, x = month - 2, xend = month + 1),
               pch = 15, size = 0.8, arrow = arrow(angle = 90, type = "closed", length = unit(0.1, "in"))) +
  scale_color_igv(name = "Treatment",
                  labels = c("VAC","VDC","VBC","VTC","GTX","IT","IE","AIM","other","doxil")) +
  scale_shape(name = "Status", solid = T) +
  scale_x_continuous((expand = c(0,0)))  +
  labs(x = "Months since diagnosis",
       y = "Patient",
       title = "Treatment timeline") +
  theme_bw() +
  theme(title = element_text(angle = 0, vjust=.5,
                             size=12, face="bold"),
        axis.title.y = element_text(angle = 90, vjust=.5,
                                    size=12, face="bold"),
        axis.title.x = element_text(size=12, face="bold"),
        axis.text.y = element_text(size=10,
                                   hjust=1),
        axis.ticks.y = element_blank()) +
  theme(#legend.position = c(0.8, 0.3),
    legend.title = element_text(colour="black",
                                size=13,
                                face=4),
    legend.text = element_text(colour="black",
                               size=10),
    legend.background = element_rect(linewidth = 0.5,
                                     linetype="solid",
                                     colour ="gray30")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()
  )
##############################
############################
###Cell Chat####
setwd("/path/scRNA/RMS/seurat/ARACNe")
MYOD1.integrated.vp <- readRDS("MYOD1.integrated.filt.vp.network.sample_network.indiv.sample.rds")
MYOD1.obj.integrated <- readRDS("MYOD1.obj.filt.integrated.v2.rds")

# data.input <- MYOD1.integrated.vp[["RNA"]]$data
# labels <- MYOD1.integrated.vp$sample_cluster
# samples <- MYOD1.integrated.vp$sample

# data.input <- MYOD1.obj.integrated[["SCT"]]@data
# MYOD1.obj.integrated$sample_cluster <- MYOD1.integrated.vp$sample_cluster
# labels <- MYOD1.obj.integrated$sample_cluster
# samples <- MYOD1.obj.integrated$sample

# meta <- data.frame(labels = labels, rownames = names(labels), samples = factor(samples))

setwd("/path/scRNA/RMS")
sample_list <- c("MYOD1_RMS_3_2", "MYOD1_RMS_31_2", "MYOD1_RMS_41B", "MYOD1_RMS_211_2", "MYOD1_RMS_469", "MYOD1_RMS_475")
# MYOD1.subset.obj.list <- readRDS(file = "MYOD1.subset.obj.list.rds")
position <- which(sample_list == "MYOD1_RMS_41B")
position <- which(sample_list == "MYOD1_RMS_469")
position <- which(sample_list == "MYOD1_RMS_475")
position <- which(sample_list == "MYOD1_RMS_31_2")
position <- which(sample_list == "MYOD1_RMS_3_2")
position <- which(sample_list == "MYOD1_RMS_211_2")

MYOD1.obj <- readRDS(paste0(sample_list[[position]], "_subset.rds"))
MYOD1.obj <- NormalizeData(MYOD1.obj, verbose = T)
colnames(MYOD1.obj) <- paste0(colnames(MYOD1.obj), "_", position)
MYOD1.vp <- subset(MYOD1.integrated.vp, sample == sample_list[[position]])
common_cells <- intersect(colnames(MYOD1.obj), colnames(MYOD1.vp))
MYOD1.obj <- subset(MYOD1.obj, cells = common_cells)
MYOD1.obj$sample_cluster <- MYOD1.vp$sample_cluster
MYOD1.obj$samples <- MYOD1.vp$sample

#create cell chat object
# cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

Idents(MYOD1.obj) <- "sample_cluster"
#MYOD1_RMS_41B
levels(MYOD1.obj) <- c("MYOD1_RMS_41B_cluster_2", "MYOD1_RMS_41B_cluster_0", "MYOD1_RMS_41B_cluster_1", "MYOD1_RMS_41B_cluster_3")
new_levels <- c("group4", "group2", "groupx", "group1")
#MYOD1_RMS_469
levels(MYOD1.obj) <- c("MYOD1_RMS_469_cluster_0", "MYOD1_RMS_469_cluster_2", "MYOD1_RMS_469_cluster_1")
new_levels <- c("group1", "group4", "group2")
#MYOD1_RMS_475
levels(MYOD1.obj) <- c("MYOD1_RMS_475_cluster_2", "MYOD1_RMS_475_cluster_0", "MYOD1_RMS_475_cluster_1", "MYOD1_RMS_475_cluster_3",
                       "MYOD1_RMS_475_cluster_4")
new_levels <- c("group1", "groupx", "group3", "group2", "group4")
#MYOD1_RMS_31_2
levels(MYOD1.obj) <- c("MYOD1_RMS_31_2_cluster_2", "MYOD1_RMS_31_2_cluster_0", "MYOD1_RMS_31_2_cluster_1")
new_levels <- c("group4", "group1", "group2")
#MYOD1_RMS_3_2
levels(MYOD1.obj) <- c("MYOD1_RMS_3_2_cluster_0", "MYOD1_RMS_3_2_cluster_1")
new_levels <- c("group1", "group2")
#MYOD1_RMS_211_2
levels(MYOD1.obj) <- c("MYOD1_RMS_211_2_cluster_1", "MYOD1_RMS_211_2_cluster_0", "MYOD1_RMS_211_2_cluster_2")
new_levels <- c("group2", "group3", "group1")
#MYOD1_PDX
levels(MYOD1.obj) <- c("0","1")
new_levels <- c("group1", "group2")

names(new_levels) <- levels(MYOD1.obj)
MYOD1.obj <- RenameIdents(MYOD1.obj, new_levels)
levels(MYOD1.obj)
MYOD1.obj <- subset(MYOD1.obj, idents = c("group1", "group2", "group4"))
Idents(MYOD1.obj) <- factor(Idents(MYOD1.obj), levels = c("group1", "group2", "group4"))
levels(MYOD1.obj)

data.input <- MYOD1.obj[["RNA"]]$data
labels <- MYOD1.obj$sample_cluster
meta <- data.frame(labels = labels, rownames = names(labels), samples = factor(sample_list[[position]]))


setwd("/path/scRNA/RMS/MYOD1_PDX_PRMS/seurat")
MYOD1.obj <- readRDS("MYOD1_PDX_PRMS_subset.rds")
MYOD1.obj <- NormalizeData(MYOD1.obj, verbose = T)
MYOD1.vp <- readRDS("MYOD1_PDX_PRMS.viper.sample_network.rds")
common_cells <- intersect(colnames(MYOD1.obj), colnames(MYOD1.vp))
MYOD1.obj <- subset(MYOD1.obj, cells = common_cells)
MYOD1.vp <- subset(MYOD1.vp, cells = common_cells)
MYOD1.obj$sample_cluster <- MYOD1.vp$seurat_clusters

# MYOD1.obj <- subset(MYOD1.obj, idents = c("MYOD1_RMS_469_cluster_0", "MYOD1_RMS_469_cluster_2", "MYOD1_RMS_469_cluster_1"))
cellchat <- createCellChat(object = MYOD1.obj, group.by = "ident", assay = "RNA")

# Set the ligand-receptor interaction database
# Show the structure of the database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)
# use all CellChatDB for cell-cell communication analysis
# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# Compute the communication probability
# ptm = Sys.time()
options(future.globals.maxSize = 8000 * 1024^2)

# project gene expression data onto PPI 
#(Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` 
#in order to use the projected data)
# cellchat <- smoothData(cellchat, adj = PPI.human)
# cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)
# saveRDS(cellchat, file = "MYOD1.integrated.vp_cellchat_PPI.rds")

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Extract the inferred cellular communication network as a data frame
# df.net <- subsetCommunication(cellchat)
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = paste0(sample_list[[position]], "_cellchat.rds"))

saveRDS(cellchat, file = paste0("MYOD1_PDX_PRMS", "_cellchat.rds"))

cellchat <- readRDS(paste0(sample_list[[position]], "_cellchat.rds"))
#Visualization of cell-cell communication network
pathways.show.all <- cellchat@netP$pathways
pathways.show <- c("IGF")
netVisual_aggregate(cellchat, signaling = pathways.show, 
                    layout = "circle", color.use = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.IGF <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.IGF[1,]
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

#Visualization of cell–cell communication mediated by multiple L–R or signaling pathways
#Visualize the inferred significant interactions using a Bubble plot
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2), remove.isolate = FALSE)
#Show all the L–R mediated interactions sending from ‘Inflam.FIB’ defined by ‘sources.use’
netVisual_chord_gene(cellchat, sources.use = 3, targets.use = c(1:2), lab.cex = 0.5,legend.pos.y = 30)
#Show all the signaling pathways mediated interactions by setting ‘slot.name’ as ‘netP’
netVisual_chord_gene(cellchat, sources.use = 3, targets.use = c(1:2), slot.name = "netP", legend.pos.x = 10)

#Visualize signaling gene expression using CellChat built-in function
plotGeneExpression(cellchat, signaling = "IGF", enriched.only = TRUE, type = "violin")

#Identify the signaling roles and major contributing signaling events of cell groups
#Compute the network centrality scores of the inferred cell–cell communication network
pathways.show <- c("JAM")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
###############################################