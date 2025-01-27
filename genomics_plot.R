################################################
##################################
##COMPLEX HEATMAP
library(ComplexHeatmap)
MYOD1_oncoprint <- read.csv("path/MYOD1_oncoprint.csv")
MYOD1_oncoprint <- MYOD1_oncoprint %>% 
  arrange(Treatment, IGF2, PIK3CA, PTEN, NRAS, MDM2, CDK4, CDKN2A, GATA3) %>% 
  select(SAMPLE_ID,	Sex,	Age_at_Diagnosis,	Primary_Tumor_Site,	Tumor_Size_cm,	Treatment,	Sample,
         MYOD1, IGF2, PIK3CA, PTEN, AKT, PIK3R3, PIK3C2G, NRAS, NF1, NF2, BCOR, FGFR4, 
         MDM2, CDK4, 
         MGA, CDKN2A, CHEK2, GATA3, SMARCB1, ARID1A, NOTCH4, 
         chr1q,	chr5,	chr10,	chr11p,	chr13q,	chr16q,	chr19p,	chr20) %>% 
  filter(SAMPLE_ID != "P-0086519-T01-IM7")

MYOD1_oncoprint <- MYOD1_oncoprint %>% 
  filter(Treatment == "Non-treated")
MYOD1_oncoprint <- MYOD1_oncoprint %>% 
  filter(Treatment == "Post-treatment")

MYOD1_ha1 <- HeatmapAnnotation(
  Age = anno_barplot(as.numeric(MYOD1_oncoprint[[3]]),
                     ylim = c(0, 62), height = unit(1, "cm"),
                     bar_width = 1, gp = gpar(fill = 3, col = "white"),
                     axis_param = list(side = "right",
                                       at = c(0, 20, 40, 60),
                                       labels = c("0", "20", "40","60"))),
  Sex = MYOD1_oncoprint[[2]],
  Site = MYOD1_oncoprint[[4]],
  col = list(
    Sex = c("FEMALE" = "#e41a1c", "MALE" = "#377eb8"),
    Site = c("Head and neck" = "#66a61e", "Upper Extremity" = "#6a51a3", "Lower Extremity" = "#7570b3", 
             "Thorax" = "#fec44f", "Abdomen" = "#2b8cbe", "Trunk" = "#d9f0a3",
             "Pelvis" = "#e7298a")
  ),
  show_annotation_name = T,
  annotation_name_rot = 0,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Sex = list(title = "Sex",
               at = c("FEMALE", "MALE"),
               labels = c("Female", "Male")
    ),
    Site = list(title = "Primary Site",
                at = c("Head and neck", "Upper Extremity", "Lower Extremity", 
                       "Thorax", "Abdomen", "Trunk",
                       "Pelvis"),
                labels = c("Head and neck", "Upper Extremity", "Lower Extremity", 
                           "Thorax", "Abdomen", "Trunk",
                           "Pelvis")
    )
  )
)

MYOD1_ha2 <- HeatmapAnnotation(
  Size = anno_barplot(as.numeric(MYOD1_oncoprint[[5]]),
                      ylim = c(0, 30), height = unit(1, "cm"),
                      bar_width = 1, gp = gpar(fill = 1),
                      axis_param = list(side = "right",
                                        at = c(0, 10, 20, 30),
                                        labels = c("0", "10", "20","30"))),
  Treatment = MYOD1_oncoprint[[6]],
  Sample = MYOD1_oncoprint[[7]],
  col = list(
    Treatment = c("Non-treated" = "#0571b0", "Post-treatment" = "#8c510a"),
    Sample = c("Biopsy" = "#66c2a4", "Resection" = "#006d2c", 
               "Local recurrence" = "#e41a1c", 
               "Metastasis" = "#984ea3")
  ),
  show_annotation_name = T,
  annotation_name_rot = 0,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Treatment = list(title = "Treatment Status",
                     at = c("Non-treated", "Post-treatment"),
                     labels = c("Untreated", "Post-treatment")
    ),
    Sample = list(title = "Sample Type",
                  at = c("Biopsy", "Resection", "Local recurrence", 
                         "Metastasis"),
                  labels = c("Biopsy", "Resection", "Local recurrence", 
                             "Metastasis")
    )
  )
)

col_mut <- c("nonsynonymous_SNV"="#005a32",
             "nonframeshift_deletion"="#bf812d",
             "stopgain_SNV"="#252525",
             "stoploss_SNV"="#737373",
             "frameshift_insertion" ="#de2d26",
             "frameshift_deletion"="#67001f",
             # "upstream"="#542788",
             "splicing"="#542788",
             "INTRAGENIC"="#2166ac", 
             "Amplification"="#b2182b",
             "Deletion"="#2166ac"#,
             # "structural_variant" = "#88419d"
)

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = "white"))
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Amplification"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.9, 
              gp = gpar(fill = col_mut["Deletion"], col = NA))
  },
  nonframeshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.3, "mm"), h*0.5, 
              gp = gpar(fill = col_mut["nonframeshift_deletion"], col = NA))
  },
  frameshift_deletion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_deletion"], col = NA))
  },
  frameshift_insertion = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.5, 
              gp = gpar(fill = col_mut["frameshift_insertion"], col = NA))
  },
  stopgain_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stopgain_SNV"], col = NA))
  },
  stoploss_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8,  h*0.5, 
              gp = gpar(fill = col_mut["stoploss_SNV"], col = NA))
  },
  nonsynonymous_SNV = function(x, y, w, h) {
    grid.rect(x, y, w*0.8, h*0.3, 
              gp = gpar(fill = col_mut["nonsynonymous_SNV"], col = NA))
  },
  structural_variant = function(x, y, w, h) {
    grid.polygon(
      unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
      unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
      gp = gpar(fill = col_mut["structural_variant"], col = NA))
  },
  splicing = function(x, y, w, h) {
    grid.segments(x - w*0.3, y - h*0.3, x + w*0.3, y + h*0.3,
                  gp = gpar(lwd = 1))
  }
)

heatmap_legend_mut <- list(title = "Alteration Type",
                           at = c("nonsynonymous_SNV",
                                  "nonframeshift_deletion",
                                  "stopgain_SNV",
                                  #"stoploss_SNV",
                                  "frameshift_insertion",
                                  "frameshift_deletion",
                                  "splicing",
                                  "structural_variant",
                                  "Amplification",
                                  "Deletion"),
                           labels = c("Missense Mutation",
                                      "In-frame deletion",
                                      "Nonsense Mutation",
                                      #"Stop Loss Mutation",
                                      "Frameshift Insertion",
                                      "Frameshift Deletion ",
                                      "Splicing",
                                      "Structural Variant",
                                      "Amplification",
                                      "Deletion"),
                           title_position = "topleft",
                           ncol = 3)

MYOD1_mut_mat <- t(MYOD1_oncoprint[8:ncol(MYOD1_oncoprint)])
colnames(MYOD1_mut_mat) <- MYOD1_oncoprint[,1]
MYOD1_mut_mat[is.na(MYOD1_mut_mat)] <- ""

MYOD1_mut_mat <- oncoPrint(MYOD1_mut_mat,
                           column_title = "MYOD1", cluster_column_slices = F, 
                           column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                           remove_empty_rows = F, remove_empty_columns = F,
                           column_order = 1:ncol(MYOD1_mut_mat),
                           alter_fun = alter_fun,
                           col = col_mut, 
                           show_heatmap_legend = T,
                           heatmap_legend_param = heatmap_legend_mut,
                           height = unit(0.4 * nrow(MYOD1_mut_mat), "cm"),
                           width = unit(0.2 * ncol(MYOD1_mut_mat), "cm"),
                           bottom_annotation = MYOD1_ha2,
                           top_annotation = MYOD1_ha1,
                           # top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(show_fraction = FALSE,
                           #                                                                  height = unit(1.0, "cm"),
                           #                                                                  ylim = c(0, 10))),
                           right_annotation = NULL,
                           # right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot(show_fraction = TRUE,
                           #                                                                width = unit(1, "cm"),
                           #                                                                ylim = c(0,1))),
                           show_pct = TRUE, show_row_names = T,
                           row_names_side = "right", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                           #row_title = "Alterations",
                           row_title_gp = gpar(fontsize = 10),
                           row_order = 1:nrow(MYOD1_mut_mat),
                           row_title_rot = 0, row_title_side = "left"
)

ht_list <- MYOD1_mut_mat
ht_list <- MYOD1_mut_mat1 + MYOD1_mut_mat2
ComplexHeatmap::draw(ht_list, heatmap_legend_side = "bottom", merge_legend = FALSE)
#################################################
############average window######
average_in_window = function(window, gr, v, method = "weighted", empty_v = NA) {
  
  if(missing(v)) v = rep(1, length(gr))
  if(is.null(v)) v = rep(1, length(gr))
  if(is.atomic(v) && is.vector(v)) v = cbind(v)
  
  v = as.matrix(v)
  if(is.character(v) && ncol(v) > 1) {
    stop("`v` can only be a character vector.")
  }
  
  if(length(empty_v) == 1) {
    empty_v = rep(empty_v, ncol(v))
  }
  
  u = matrix(rep(empty_v, each = length(window)), nrow = length(window), ncol = ncol(v))
  
  mtch = as.matrix(findOverlaps(window, gr))
  intersect = pintersect(window[mtch[,1]], gr[mtch[,2]])
  w = width(intersect)
  v = v[mtch[,2], , drop = FALSE]
  n = nrow(v)
  
  ind_list = split(seq_len(n), mtch[, 1])
  window_index = as.numeric(names(ind_list))
  window_w = width(window)
  
  if(is.character(v)) {
    for(i in seq_along(ind_list)) {
      ind = ind_list[[i]]
      if(is.function(method)) {
        u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
      } else {
        tb = tapply(w[ind], v[ind], sum)
        u[window_index[i], ] = names(tb[which.max(tb)])
      }
    }
  } else {
    if(method == "w0") {
      gr2 = reduce(gr, min.gapwidth = 0)
      mtch2 = as.matrix(findOverlaps(window, gr2))
      intersect2 = pintersect(window[mtch2[, 1]], gr2[mtch2[, 2]])
      
      width_intersect = tapply(width(intersect2), mtch2[, 1], sum)
      ind = unique(mtch2[, 1])
      width_setdiff = width(window[ind]) - width_intersect
      
      w2 = width(window[ind])
      
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        x = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
        u[window_index[i], ] = (x*width_intersect[i] + empty_v*width_setdiff[i])/w2[i]
      }
      
    } else if(method == "absolute") {
      for(i in seq_along(ind_list)) {
        u[window_index[i], ] = colMeans(v[ind_list[[i]], , drop = FALSE])
      }
      
    } else if(method == "weighted") {
      for(i in seq_along(ind_list)) {
        ind = ind_list[[i]]
        u[window_index[i], ] = colSums(v[ind, , drop = FALSE]*w[ind])/sum(w[ind])
      }
    } else {
      if(is.function(method)) {
        for(i in seq_along(ind_list)) {
          ind = ind_list[[i]]
          u[window_index[i], ] = method(v[ind], w[ind], window_w[i])
        }
      } else {
        stop("wrong method.")
      }
    }
  }
  
  return(u)
}
########################################################################################
library(CNTools)
library(factoextra)
library(GenomicRanges)
library(EnrichedHeatmap)
library(circlize)
library(ComplexHeatmap)

MYOD1_oncoprint <- read.csv("path/MYOD1_oncoprint.csv")
NULL
####MYOD1 copy number clustering####
#download copy number .seg profile from cbioportal for all MYOD1 samples
MYOD1_copynumber <- read.delim("path/MYOD1_copynumber.seg")

MYOD1_cn_seg <- CNTools::CNSeg(segList = MYOD1_copynumber, chromosome = "chrom", 
                               end = "loc.end", start = "loc.start", segMean = "seg.mean",
                               id = "ID")
MYOD1_cn_reduced_seg <- CNTools::getRS(MYOD1_cn_seg, XY = FALSE, what = "mean", imput = FALSE, by = "region")
MYOD1_cn <- CNTools::rs(MYOD1_cn_reduced_seg)
MYOD1_cn <- MYOD1_cn %>% mutate(chrom = as.numeric(chrom),start=as.numeric(start),end=as.numeric(end)) %>% 
  mutate(seglength = as.numeric(end) - as.numeric(start)) %>% filter(seglength > 1) %>% select(-seglength)

chr_gr <- GRanges(seqnames = MYOD1_cn[,1], ranges = IRanges(MYOD1_cn[,2] + 1, MYOD1_cn[,3]))
chr_window <- EnrichedHeatmap::makeWindows(chr_gr, w = 1e6)
MYOD1_cn_mat <- matrix(as.numeric(unlist(MYOD1_cn)),nrow=nrow(MYOD1_cn))
num_mat <- average_in_window(chr_window, chr_gr, MYOD1_cn_mat[, -(1:3)])
MYOD1_cn_df <- cbind(as.data.frame(chr_window),as.data.frame(num_mat))
MYOD1_cn_df <- MYOD1_cn_df %>% select(-c(4:7))
colnames(MYOD1_cn_df) <- colnames(MYOD1_cn)

MYOD1_cn_df_wide <- MYOD1_cn_df %>% 
  pivot_longer(cols = starts_with("P-"),
               names_to = "SAMPLE_ID",
               values_to = "seg_mean") %>% 
  pivot_wider(names_from = c(chrom,start,end),
              names_sep = ".",
              values_from = seg_mean) %>% 
  filter(SAMPLE_ID != "P-0086519-T01-IM7") %>% 
  column_to_rownames("SAMPLE_ID")

MYOD1_cn_mat <- matrix(as.numeric(unlist(MYOD1_cn_df_wide)),nrow=nrow(MYOD1_cn_df_wide))
row.names(MYOD1_cn_mat) <- row.names(MYOD1_cn_df_wide)
colnames(MYOD1_cn_mat) <- colnames(MYOD1_cn_df_wide)

chrom_anno <- factor(MYOD1_cn_df$chrom)
levels(chrom_anno) <- 1:22

col_fun = colorRamp2(c(-1,0,1), 
                     c("#377EB8", "#f7f7f7", "#E41A1C"))

chr_level <- 1:22

MYOD1_anno <- data.frame(SAMPLE_ID = rownames(MYOD1_cn_mat)) %>% 
  left_join(MYOD1_oncoprint %>% select(SAMPLE_ID, Treatment, Sample, cluster))
MYOD1_cluster <- MYOD1_anno$cluster
MYOD1_cluster <- MYOD1_anno$Treatment

MYOD1_ha <- rowAnnotation(
  Treatment = MYOD1_anno[[2]],
  Sample = MYOD1_anno[[3]],
  col = list(
    Treatment = c("Non-treated" = "#0571b0", "Post-treatment" = "#8c510a"),
    Sample = c("Biopsy" = "#66c2a4", "Resection" = "#006d2c", 
               "Local recurrence" = "#e41a1c", 
               "Metastasis" = "#984ea3")
  ),
  show_annotation_name = T,
  annotation_name_rot = 90,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Treatment = list(title = "Treatment Status",
                     at = c("Non-treated", "Post-treatment"),
                     labels = c("Untreated", "Post-treatment")
    ),
    Sample = list(title = "Sample Type",
                  at = c("Biopsy", "Resection", "Local recurrence", 
                         "Metastasis"),
                  labels = c("Biopsy", "Resection", "Local recurrence", 
                             "Metastasis")
    )
  )
)

set.seed(123)
MYOD1_cn_mat_ht <-
  Heatmap(MYOD1_cn_mat, col = col_fun, 
          cluster_columns = FALSE, show_column_names = FALSE, column_split = chrom_anno, cluster_column_slices = FALSE, 
          cluster_row = TRUE, clustering_distance_rows = "euclidean", row_gap = unit(0.5, "mm"), cluster_row_slices = FALSE, 
          #row_title_side = "left", border = TRUE, row_title_gp = gpar(col = "black", fontsize = 8),
          row_split = MYOD1_cluster,
          show_row_names = F, 
          row_names_side = "left", row_names_gp = gpar(fontsize = 10, fontface = "italic"),
          row_title_rot = 0,
          column_title = (ifelse(1:22 %% 2 == 0, paste0("\n", chr_level), paste0(chr_level, "\n"))),
          column_title_gp = gpar(fill = c(rep(c("azure4","azure3"), 11)), col = "white", fontsize = 10), 
          column_gap = unit(0, "mm"), column_title_side = "top", 
          #row_title = ifelse(1:22 %% 2 ==0, paste0("    ",levels(chrom_anno)), paste0(levels(chrom_anno),"    ")), 
          height = unit(0.4*nrow(MYOD1_cn_df_wide), "cm"), width = unit(20, "cm"),
          right_annotation = MYOD1_ha,
          #bottom_annotation = MYOD1_cn_ha, top_annotation = ha,
          heatmap_legend_param = list(title = "Copy Number",legend_direction = "horizontal", legend_width = unit(3,"cm")))

ht <- MYOD1_cn_mat_ht
ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
########################################################