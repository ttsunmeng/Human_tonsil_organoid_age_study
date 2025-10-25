# install.packages('remotes')
# remotes::install_version("Seurat", version = "3.2.3")
# library(Seurat)

library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(ggrepel)
library("readxl")
library('stringr')
setwd("/Volumes/GoogleDrive/My\ Drive/Stanford/RNA-seq/data_analysis/tonsil_LAIV_Rhapsody_120a-s/")

table1 <- 'tonsil_LAIV_120a_s_background_list.xlsx'
all_marker_table <- read_excel(table1, sheet = "Sheet1")

Bcell.subset <- readRDS("tonsil_LAIV_120a_s_cleaned_log1p_merge_scale_Bcell_cca_integrated.rds")
all.markers <- rownames(Bcell.subset)
condition_list <- unique(Bcell.subset$condition)
age_group_donorID_list <- list()
age_group_donorID_list[[1]] <- c("02yrs F IMD030","02yrs M IMD082","02yrs M IMD085","03yrs F IMD035")
age_group_donorID_list[[2]] <- c("07yrs F IMD150","07yrs M IMD170","09yrs M IMD107")
age_group_donorID_list[[3]] <- c("27yrs F VIP031","33yrs M VIP024","39yrs M VIP015")
Bcell.subset$age_group[Bcell.subset$age_group == "06-09yrs"] <- "07-09yrs"
Bcell.subset$gender <- substring(Bcell.subset$donor_ID, first = 7, last = 7)
Bcell.subset$day_group <- 'day0'
Bcell.subset$day_group[Bcell.subset$days %in% c('day02')] <- 'day02'
Bcell.subset$day_group[Bcell.subset$days %in% c('day04')] <- 'day04'
Bcell.subset$day_group[Bcell.subset$days %in% c('day06','day07')] <- 'day06-07'
Bcell.subset$day_group[Bcell.subset$days %in% c('day08')] <- 'day08'
Bcell.subset$day_group[Bcell.subset$days %in% c('day10')] <- 'day10'
Bcell.subset$day_group[Bcell.subset$days %in% c('day12','day14')] <- 'day12-14'
Bcell.subset$stimulation <- Bcell.subset$treatment
age_group_list <- sort(unique(Bcell.subset$age_group))
all.genes <- all.markers[!grepl('ab-',all.markers)]
all.Ab <- all.markers[grepl('ab-',all.markers)]
batch_marker_list <- list()
batch_marker_list[[1]] <- readRDS('tonsil_LAIV_120a.rds')
batch_marker_list[[2]] <- readRDS('tonsil_LAIV_120b.rds')
batch_marker_list[[3]] <- readRDS('tonsil_LAIV_120c.rds')
batch_marker_list[[4]] <- readRDS('tonsil_LAIV_120j.rds')
overlapped_marker_list <- Reduce(intersect, batch_marker_list)
batch1_donorID_list <- unique(Bcell.subset$donor_ID[Bcell.subset$batch == 'batch1'])
batch2_donorID_list <- unique(Bcell.subset$donor_ID[Bcell.subset$batch == 'batch2'])
batch3_donorID_list <- unique(Bcell.subset$donor_ID[Bcell.subset$batch == 'batch3'])
batch4_donorID_list <- unique(Bcell.subset$donor_ID[Bcell.subset$batch == 'batch4'])
all_donorID_list <- unique(Bcell.subset$donor_ID)
day0_condition_list <- condition_list[grepl('day0$',condition_list)]
nonday0_condition_list <- condition_list[!grepl('day0$',condition_list)]

############## UMAP clustering ######################################################
Bcell.subset[["lognorm"]] <- CreateAssayObject(counts = as.matrix(Bcell.subset@assays[["raw"]]@counts))
DefaultAssay(Bcell.subset) <- "lognorm"
Bcell.subset <- NormalizeData(Bcell.subset)

DefaultAssay(Bcell.subset) <- "integrated"
s.genes <- cc.genes$s.genes
s.genes <- s.genes[s.genes %in% all.markers]
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- g2m.genes[g2m.genes %in% all.markers]
Bcell.subset <- CellCycleScoring(Bcell.subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Bcell.subset <- ScaleData(Bcell.subset,vars.to.regress = c("S.Score", "G2M.Score"),features = all.markers,)#,scale.max = 6) # 3.5: ;3: 1.3% of data 

scaled_integrated_assay <- CreateAssayObject(counts = Bcell.subset@assays[["integrated"]]@scale.data)
Bcell.subset[["integrated_scale"]] <- scaled_integrated_assay
Bcell.subset <- CellCycleScoring(Bcell.subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Bcell.subset <- ScaleData(Bcell.subset,vars.to.regress = c("S.Score", "G2M.Score"),features = all.markers)
Bcell.subset@assays[["integrated_scale"]]@scale.data <- as.matrix(Bcell.subset@assays[["integrated_scale"]]@data)
DefaultAssay(Bcell.subset) <- "integrated_scale"
rm(scaled_integrated_assay)

Bcell.subset <- RunPCA(Bcell.subset, features = all.markers, npcs = 50)
# ElbowPlot(object = Bcell.subset,ndims = 50) + theme(axis.text = element_text(size = 20))
# dev.print(pdf, 'tonsil_LAIV_120a_scleaned_log1p_cca_scale_CycleRegressOut_PC_elbow.pdf',width = 700, height = 4.34)
n_pca_selected <- 30
Bcell.subset <- RunUMAP(Bcell.subset, reduction = "pca", dims = 1:n_pca_selected)

DefaultAssay(Bcell.subset) <- 'raw'
Bcell.subset[['lognorm']] <- NULL
all.markers <- rownames(Bcell.subset)
all.genes <- all.markers[!grepl('ab-',all.markers)]
all.Ab <- all.markers[grepl('ab-',all.markers)]
Bcell.subset[["RNA"]] <- CreateAssayObject(counts = as.matrix(Bcell.subset@assays[["raw"]]@counts)[all.genes,])
Bcell.subset[["ADT"]] <- CreateAssayObject(counts = as.matrix(Bcell.subset@assays[["raw"]]@counts)[all.Ab,])
DefaultAssay(Bcell.subset) <- 'RNA'

VlnPlot(Bcell.subset, feature = "nCount_raw",group.by = 'batch',log = F,pt.size = 0.1) + geom_hline(yintercept = c(8*10^4))
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_raw_nCount_batch.pdf',width = 8, height = 6)

VlnPlot(Bcell.subset, feature = "nFeature_raw",group.by = 'batch',pt.size = 0.1) + geom_hline(yintercept = c(30,200))
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_raw_nFeature_batch.pdf',width = 5, height = 4)

FeatureScatter(Bcell.subset, feature1 = "nCount_raw", feature2 = "nFeature_raw") + geom_hline(yintercept = c(30,200))# + geom_vline(xintercept = 1e5)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_raw_nCount_nFeature.pdf',width = 5, height = 4)

Bcell.subset <- subset(Bcell.subset, nFeature_raw > 60 & nFeature_RNA > 30 & nFeature_raw < 240 & nCount_raw < 8*10^4)

FeatureScatter(Bcell.subset, feature1 = "nCount_raw", feature2 = "nFeature_raw") 
FeatureScatter(Bcell.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 

DefaultAssay(Bcell.subset) <- 'RNA'
Bcell.subset <- NormalizeData(Bcell.subset)
DefaultAssay(Bcell.subset) <- 'ADT'
Bcell.subset <- NormalizeData(Bcell.subset, normalization.method = 'CLR', margin = 2)
Bcell.subset[['normalized']] <- CreateAssayObject(rbind(as.matrix(Bcell.subset@assays[["ADT"]]@data),as.matrix(Bcell.subset@assays[["RNA"]]@data)))
DefaultAssay(Bcell.subset) <- 'normalized'

Bcell.subset <- subset(Bcell.subset, clusters != 13)

Bcell_LAIV.subset <- subset(Bcell.subset, stimulation != 'LAIV-')
Bcell_LAIV.subset <- subset(Bcell_LAIV.subset,days != 'day02')
Bcell_LAIV.subset <- subset(Bcell_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))

DimPlot(Bcell_LAIV.subset, label = TRUE, reduction = "umap",group.by = 'donor_ID')
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_donorID.pdf',width = 5, height = 4.34)

DimPlot(Bcell_LAIV.subset, label = TRUE, reduction = "umap",group.by = 'batch')
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_batch.pdf',width = 5, height = 4.34)
DimPlot(Bcell_LAIV.subset, reduction = "umap", split.by = "stimulation", ncol = 3)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_stim_split.pdf',width = 10, height = 8.8)
DimPlot(Bcell_LAIV.subset, reduction = "umap", group.by = "stimulation")
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_stim.pdf',width = 5, height = 4.34)

DimPlot(Bcell_LAIV.subset, label = TRUE, reduction = "umap",group.by = 'days')
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_days.pdf',width = 5, height = 4.34)
DimPlot(Bcell_LAIV.subset, reduction = "umap", split.by = "days", ncol = 2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_days.pdf',width = 10, height = 8.8)
DimPlot(Bcell_LAIV.subset, reduction = "umap", split.by = "day_group", ncol = 2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_day_group.pdf',width = 10, height = 8.8)
DimPlot(Bcell_LAIV.subset, reduction = "umap", group.by = "age_group") +
  scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99'))
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_age_group.pdf',width = 7, height = 5)

KO_cells <- 1:dim(Bcell_LAIV.subset)[2]
downsampled_KO_cells <- sample(KO_cells, dim(Bcell_LAIV.subset)[2]/10)
temp <- Bcell_LAIV.subset[,downsampled_KO_cells]
# DimPlot(WT_KO_integrated_downsampled, reduction = "umap", split.by ="orig.ident", ncol=2)
DimPlot(temp, reduction = "umap", pt.size = 0.03, group.by = "cluster_annotated",split.by = "age_group", ncol = 3) + 
  theme(text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold")) 
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_age_group_split_cluster_annotated.pdf',width = 7, height = 2.5)

Bcell.subset <- FindNeighbors(Bcell.subset, reduction = "pca", dims = 1:n_pca_selected)
Bcell.subset <- FindClusters(Bcell.subset, resolution = 0.4)
Bcell.subset$subclusters <- Bcell.subset$seurat_clusters
Bcell.subset$clusters <- Bcell.subset$seurat_clusters
DimPlot(Bcell_LAIV.subset, label = TRUE, reduction = "umap",group.by = 'clusters',repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_clusters.pdf',width = 4.1, height = 3.6)

Bcell_LAIV.subset <- subset(Bcell.subset,stimulation != 'LAIV-')
DimPlot(Bcell_LAIV.subset, label = TRUE, reduction = "umap",group.by = 'subclusters',repel = TRUE)


FeaturePlot(Bcell.subset, features = c("XBP1","PRDM1"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_cleaned_umap_XBP1_PRDM1.pdf',width = 1200, height = 4.34)

Bcell.subset$if_PB <- 'nonPB Bcells'
# Bcell.subset$if_PB <- 0
PB_index <- c(13,15,23,30,32) 
Bcell.subset$if_PB[Bcell.subset$seurat_clusters %in% PB_index] <- 'PB'
Bcell.subset$differentiation <- Bcell.subset$if_PB
DimPlot(Bcell.subset, label = TRUE, reduction = "umap",repel = TRUE,group.by = 'differentiation') + ggtitle('tonsil orgnaoid LAIV B cells') + scale_fill_manual(values=c("#999999", "#E69F00"))


Bcell.markers_subcluster <- FindAllMarkers(Bcell.subset)
Bcell.markers_subcluster <- Bcell.markers_subcluster %>% filter(p_val_adj <= 0.05)
p <- match(Bcell.markers_subcluster$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
Bcell.markers_subcluster$ENTREZID <- temp$ENTREZID
write.xlsx(Bcell.markers_subcluster, "tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_subcluster_gene.xlsx")
for (cluster_name in sort(unique(Bcell.subset$subclusters))[31:length(unique(Bcell.subset$subclusters))]){
  temp <- Bcell.markers_subcluster[Bcell.markers_subcluster$cluster == cluster_name,]
  temp <- temp %>% arrange(desc(avg_log2FC))
  write.xlsx(temp, "tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_subcluster_gene copy.xlsx",sheetName = paste('cluster',cluster_name,sep = ''),append = T)
}
Bcell.markers_subcluster.top10 <- Bcell.markers_subcluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Bcell.markers_subcluster.top10, "tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_subcluster_gene_top10.csv")

Bcell.markers_cluster <- FindAllMarkers(Bcell.subset)
Bcell.markers_cluster <- Bcell.markers_cluster %>% filter(p_val_adj <= 0.05)
p <- match(Bcell.markers_cluster$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
Bcell.markers_cluster$ENTREZID <- temp$ENTREZID
write.xlsx(Bcell.markers_cluster, "tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_cluster_gene.xlsx")
for (cluster_name in sort(unique(Bcell.subset$clusters))[10:15]){
  temp <- Bcell.markers_cluster[Bcell.markers_cluster$cluster == cluster_name,]
  temp <- temp %>% arrange(desc(avg_log2FC))
  write.xlsx(temp, "tonsil_LAIV_120a_s_Bcell_cluster_gene.xlsx",sheetName = paste('cluster',cluster_name,sep = ''),append = T)
}
Bcell.markers_cluster.top10 <- Bcell.markers_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Bcell.markers_cluster.top10, "tonsil_LAIV_120a_s_Bcell_cluster_gene_top10.csv")

nonPB.subset <- subset(Bcell.subset,differentiation != 'PB')
nonPB.markers_subcluster <- FindAllMarkers(nonPB.subset)
nonPB.markers_subcluster <- nonPB.markers_subcluster %>% filter(p_val_adj <= 0.05)
p <- match(nonPB.markers_subcluster$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
nonPB.markers_subcluster$ENTREZID <- temp$ENTREZID
write.xlsx(nonPB.markers_subcluster, "tonsil_LAIV_120a_s_nonPB_log1p_cca_scale_CycleRegressOut_subcluster_gene.xlsx")
for (cluster_name in sort(unique(nonPB.markers_subcluster$cluster))[26:length(unique(nonPB.markers_subcluster$cluster))]){
  temp <- nonPB.markers_subcluster[nonPB.markers_subcluster$cluster == cluster_name,]
  temp <- temp %>% arrange(desc(avg_log2FC))
  write.xlsx(temp, "tonsil_LAIV_120a_s_nonPB_log1p_cca_scale_CycleRegressOut_subcluster_gene.xlsx",sheetName = paste('cluster',cluster_name,sep = ''),append = T)
}
nonPB.markers_subcluster.top10 <- nonPB.markers_subcluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(nonPB.markers_subcluster.top10, "tonsil_LAIV_120a_s_nonPB_log1p_cca_scale_CycleRegressOut_cluster_gene_top10.csv")

nonPB.markers_cluster <- FindAllMarkers(nonPB.subset)
nonPB.markers_cluster <- nonPB.markers_cluster %>% filter(p_val_adj <= 0.05)
p <- match(nonPB.markers_cluster$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
nonPB.markers_cluster$ENTREZID <- temp$ENTREZID
write.xlsx(nonPB.markers_cluster, "tonsil_LAIV_120a_s_nonPB_log1p_cca_scale_CycleRegressOut_cluster_gene.xlsx")
for (cluster_name in sort(unique(nonPB.markers_cluster$cluster))[26:length(unique(nonPB.markers_cluster$cluster))]){
  temp <- nonPB.markers_cluster[nonPB.markers_cluster$cluster == cluster_name,]
  temp <- temp %>% arrange(desc(avg_log2FC))
  write.xlsx(temp, "tonsil_LAIV_120a_s_nonPB_log1p_cca_scale_CycleRegressOut_cluster_gene.xlsx",sheetName = paste('cluster',cluster_name,sep = ''),append = T)
}
nonPB.markers_cluster.top10 <- nonPB.markers_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(nonPB.markers_cluster.top10, "tonsil_LAIV_120a_s_nonPB_cluster_gene_top10.csv")

Idents(Bcell.subset) <- 'CD27CD38_cluster'
temp <- FindAllMarkers(Bcell.subset)
for (cluster_name in sort(unique(Bcell.subset$CD27CD38_cluster))){
  temp1 <- temp[temp$cluster == cluster_name,]
  temp1 <- temp1 %>% arrange(desc(avg_log2FC))
  write.xlsx(temp1, "tonsil_LAIV_120a_s_Bcell_CD27CD38_cluster_gene.xlsx",sheetName = cluster_name,append = T)
}
write.xlsx(as.data.frame(temp %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)), "tonsil_LAIV_120a_s_Bcell_CD27CD38_cluster_gene_top10.xlsx")

memory_B_marker <- c('ab-CD27','CD27','FCRL4','ab-CD11c')
GC_marker <- c('ab-CD38','CD38','ab-CXCR5','CXCR5','BCL6','AICDA','MME','MYC','BACH2','MKI67','BCL2','BCL2A1','ab-CD83','CD83','ab-CD184','CXCR4')
naive_marker <- c('TCL1A')
PB_marker <- c('XBP1','PRDM1','IRF4','MZB1')
features <- list("memory" = memory_B_marker, "GC" = GC_marker, 'naive' = naive_marker,"PB" = PB_marker)
Bcell.subset$CD27CD38_cluster <- factor(Bcell.subset$CD27CD38_cluster,
                                        levels = c('PB',"BCL6highB","BCL6+B","CD27+CD38+B","naive B",'memory B'))
DotPlot(Bcell.subset, group.by = 'CD27CD38_cluster',dot.scale = 5,features=features) + RotatedAxis() +
  theme(axis.text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_CD27CD38_cluster_dotplot.pdf',sep = ''),width = 8, height = 3)

Bcell.subset$clusters <- factor(Bcell.subset$clusters,
                                        levels = c(7,8,11,9,1,14,5,3,12,0,6,10,13,2,4))
DotPlot(Bcell.subset, group.by = 'clusters',dot.scale = 5,features=features) + RotatedAxis() +
  theme(axis.text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_clusters_dotplot.pdf',sep = ''),width = 8, height = 4)

Bcell.subset$if_act <- 'others'
Bcell.subset$if_act[Bcell.subset$clusters %in% c(5,6,9)] <- 'act.'

DimPlot(Bcell.subset, label = T, reduction = "umap",group.by = 'if_act',repel = TRUE)

temp <- FindMarkers(Bcell.subset,ident.1 = 'act.',ident.2 = 'others',group.by = 'if_act')
temp <- temp %>% arrange(p_val_adj,desc(avg_log2FC))
############ B cell annotation #####################################################
grep('IGH', all.markers, value=TRUE)

# day0.subset <- subset(Bcell.subset, days == 'day0')
# day0_GC.subset <- subset(day0.subset,clusters == 11)
gene_name <- 'IGHM-secreted'
# View(Bcell.markers_subcluster %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))
c((Bcell.markers_cluster %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))['cluster'])
# View(nonPB.markers_subcluster %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))
c((nonPB.markers_cluster %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))['cluster'])
FeaturePlot(Bcell.subset, features = c(gene_name),ncol = 1, min.cutoff = 0)#, max.cutoff = 4)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_umap_',gene_name,'.pdf',sep = ''),width = 5, height = 4.3)
VlnPlot(Bcell.subset, features = c(gene_name),pt.size = 0.1, ncol = 1,split.by = 'age_group',log = F,assay = 'integrated_scale',group.by = 'major_clusters2')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_violin_log1p_age_group_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

VlnPlot(Bcell_LAIV.subset, features = c('IgM_sum'),pt.size = 0.1, ncol = 1,split.by = 'age_group',log = F, group.by = 'major_clusters2')

temp <- subset(Bcell_LAIV.subset, days %in% c('day0','day04','day06','day07','day10','day12','day14'))
temp$day_group <- 'day0'
temp$day_group[temp$days %in% c('day04')] <- 'day04'
temp$day_group[temp$days %in% c('day06','day07')] <- 'day06-07'
temp$day_group[temp$days %in% c('day10')] <- 'day10'
temp$day_group[temp$days %in% c('day12','day14')] <- 'day12-14'
# temp_PB <- subset(temp, if_PB == 'PB')
# temp_preGC <- subset(temp, major_clusters2 == 'preGC/GC B')
gene_name <- 'HLA-A'
VlnPlot(temp, features = c(gene_name),pt.size = 0.1, ncol = 1,split.by = 'age_group',log = F,assay = 'normalized',group.by = 'day_group')
dev.print(pdf, paste('tonsil_LAIV_120a_s_LAIV_Bcell_violin_normalized_age_group_',gene_name,'.pdf',sep = ''),width = 15, height = 6)


batch_marker_list <- list()
batch_marker_list[[1]] <- readRDS('tonsil_LAIV_120a.rds')
batch_marker_list[[2]] <- readRDS('tonsil_LAIV_120b.rds')
batch_marker_list[[3]] <- readRDS('tonsil_LAIV_120c.rds')
batch_marker_list[[4]] <- readRDS('tonsil_LAIV_120j.rds')
minimal_marker_list <- Reduce(intersect,batch_marker_list[2:4])
feature_list <- c('STAT1','SLC2A3','CLEC2D','S100A4','ALDOC','MX1','IFI6','LAP3','BNIP3L','TNFSF10','OAS1')
min_feature_list <- feature_list[feature_list %in% minimal_marker_list]
Bcell_LAIV.subset <- AddModuleScore(
  Bcell_LAIV.subset,
  features = list(min_feature_list),
  name = 'MX1_module',
  assay = 'normalized',
  slot = 'data',
  ctrl = 25
)
gene_name <- 'MX1_module1'
VlnPlot(temp, features = c(gene_name),pt.size = 0.1, ncol = 1,split.by = 'age_group',log = F,assay = 'normalized',group.by = 'day_group')
dev.print(pdf, paste('tonsil_LAIV_120a_s_LAIV_Bcell_violin_normalized_age_group_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

VlnPlot(Bcell.subset, features = c(gene_name),pt.size = 0, ncol = 1,group.by = 'donor_ID',log = F,assay = 'integrated_scale') + geom_hline(yintercept = 0)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_violin_integrated_scale_donorID_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

VlnPlot(Bcell.subset, features = c(gene_name),pt.size = 0.1, ncol = 1,group.by = 'batch',log = F,assay = 'lognorm') + geom_hline(yintercept = 6.5)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_violin_lognorm_batch_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

# BCL6: 1 11
# ab-CXCR5: 14 12 4  1 
# ab-CD27: 4  2  10 13 7  11
# ab-CD11c: 4  13 2 
# ab-CD38: 1  11 12 4  13
# ab-CD83: 13 9  14 1 
# ab-CD184: 5  3  12 2  13 0 

# ab-IgG: 10 13 11
# ab-IgD: 14 1  6  12 3  13 0 
# ab-IgM: 13 14 1  12 3 

# ab-CD69: 14 4  1  2  9  12
# FOSB: 10 6  9  5  7 
# JUNB: 9  10 6  8  5 
# IL4R: 11 1  9  0 
# TCL1A: 1  11 14 0 
Bcell.subset$cluster_annotated <- as.character(Bcell.subset$clusters)
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '0'] <- 'naive'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '1'] <- 'GC LZ'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '2'] <- 'memory'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '3'] <- 'cytokine-rich'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '4'] <- 'GC-like'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '5'] <- 'CXCR4+act.'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '6'] <- 'CD23+act.'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '7'] <- 'IgG+PB'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '8'] <- 'IgA+/IgM+PB'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '9'] <- 'CD83+folli.'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '10'] <- 'switched memory'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '11'] <- 'switched GC'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '12'] <- 'preGC/folli.'
Bcell.subset$cluster_annotated[Bcell.subset$cluster_annotated == '14'] <- 'prolif. CD83+folli.'

PB.subset <- subset(Bcell.subset,if_PB == 'PB')
FeatureScatter(PB.subset, feature1 = "IGHG1-secreted", feature2 = "IGHG4-secreted", group.by = 'cluster_annotated') 
FeatureScatter(PB.subset, feature1 = "IgM_sum_normalized", feature2 = "IGHA1-secreted", group.by = 'cluster_annotated') 

PB.subset$IgM_sum_normalized <- as.matrix(PB.subset@assays[["normalized"]]@counts)['IGHM-secreted',] + as.matrix(PB.subset@assays[["normalized"]]@counts)['IGHM-membrane',]

PB_secreted <- data.frame(t(as.matrix(Bcell.subset@assays$raw@data)),check.names = F)[,rownames(Bcell.subset)[grepl('-secreted',rownames(Bcell.subset))]]
PB_secreted$max <- colnames(PB_secreted)[apply(PB_secreted,1,which.max)]
Bcell.subset$subisotype <- PB_secreted$max
Bcell.subset$subisotype[!Bcell.subset$clusters %in% c(7,8)] <- 'B cells'

PB_secreted$isotype <- PB_secreted$max
PB_secreted$isotype[grepl('IGHG',PB_secreted$max)] <- 'IgG+PB'
PB_secreted$isotype[grepl('IGHM',PB_secreted$max)] <- 'IgM+PB'
PB_secreted$isotype[grepl('IGHA',PB_secreted$max)] <- 'IgA+PB'
PB_secreted$isotype[grepl('IGHE',PB_secreted$max)] <- 'IgE+PB'
Bcell.subset$isotype <- PB_secreted$isotype
Bcell.subset$isotype[!Bcell.subset$clusters %in% c(7,8)] <- 'B cells'


ordered_cluster_annotated_list <- c('naive','cytokine-rich','CXCR4+act.','CD23+act.','CD83+folli.','prolif. CD83+folli.','preGC/folli.',
                                    'switched GC','GC LZ','GC-like','memory','switched memory','IgA+/IgM+PB','IgG+PB')
Bcell_LAIV.subset$cluster_annotated <- factor(Bcell_LAIV.subset$cluster_annotated,
                                              levels = ordered_cluster_annotated_list)
DimPlot(Bcell_LAIV.subset, label = T, reduction = "umap",group.by = 'cluster_annotated',repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_cluster_annotated.pdf',width = 7, height = 5)

p <- DimPlot(Bcell_LAIV.subset, group.by = 'cluster_annotated', label = TRUE, reduction = "umap")
pbuild <- ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # Pull the data used for the plot
pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
cols_group <- unique(pdata$colour) # Get a vector of unique colors
names(cols_group) <- ordered_cluster_annotated_list


Bcell.subset$major_clusters <- as.character(Bcell.subset$clusters)
Bcell.subset$major_clusters[Bcell.subset$major_clusters %in% c('0','3')] <- 'Naive B'
Bcell.subset$major_clusters[Bcell.subset$major_clusters == '1'] <- 'GC LZ'
Bcell.subset$major_clusters[Bcell.subset$major_clusters %in% c('11')] <- 'switched GC'
Bcell.subset$major_clusters[Bcell.subset$major_clusters %in% c('4')] <- 'GC-like B'
Bcell.subset$major_clusters[Bcell.subset$major_clusters %in% c('2','10')] <- 'Memory B'
Bcell.subset$major_clusters[Bcell.subset$major_clusters %in% c('5','6')] <- 'Act. B'
Bcell.subset$major_clusters[Bcell.subset$major_clusters %in% c('7','8')] <- 'PB'
Bcell.subset$major_clusters[Bcell.subset$major_clusters %in% c('9','12','14')] <- 'preGC/folli. B'

DimPlot(Bcell.subset, label = T, reduction = "umap",group.by = 'major_clusters',repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_major_clusters.pdf',width = 7, height = 5)

Bcell.subset$major_clusters2 <- as.character(Bcell.subset$clusters)
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('0','3','5','6')] <- 'Naive/Act. B'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('11')] <- 'switched GC'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('4')] <- 'GC-like B'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('2','10')] <- 'Memory B'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('7','8')] <- 'PB'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('1','9','12','14')] <- 'preGC/GC B'

DimPlot(Bcell.subset, label = T, reduction = "umap",group.by = 'major_clusters2',repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_major_clusters2.pdf',width = 7, height = 5)


Bcell.subset$major_clusters3 <- as.character(Bcell.subset$clusters)
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('0','3','5','6')] <- 'Naive/Act. B'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('4')] <- 'GC-like B'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('2','10')] <- 'Memory B'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('1','9','12','14','11')] <- 'preGC/GC B'
Bcell.subset$major_clusters2[Bcell.subset$major_clusters2 %in% c('7','8')] <- 'PB'

DimPlot(Bcell.subset, label = T, reduction = "umap",group.by = 'major_clusters2',repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_LAIV_umap_major_clusters2.pdf',width = 7, height = 5)

Bcell_LAIV.subset$cluster_annotated <- factor(Bcell_LAIV.subset$cluster_annotated,
                                              levels = rev(ordered_cluster_annotated_list))

# memory_B_marker <- c('ab-CD27','CD27','FCRL4','ab-CD11c')
# GC_marker <- c('ab-CD38','ab-CXCR5','CXCR5','BCL6','AICDA','MME','MYC','BACH2','MKI67','BCL2','BCL2A1','ab-CD83','CD83','ab-CD184','CXCR4')
# naive_marker <- c('TCL1A')
# PB_marker <- c('XBP1','PRDM1','IRF4','MZB1')
# features <- list("memory" = memory_B_marker, "GC" = GC_marker, 'naive' = naive_marker,"PB" = PB_marker)
DefaultAssay(Bcell_LAIV.subset) <- 'integrated_scale'
features <- c('IL4R','TCL1A','ab-IgD','CSF3','IL2','IL21','ab-CD184','FOSB','JUNB','ab-CD23','ab-CD83','ab-CXCR5','MKI67','ab-IgG','BCL6','ab-CD38','ab-CD27','IGHD-membrane','IGHM-membrane',
              'XBP1','PRDM1','IGHA1-secreted','IGHM-secreted','IGHG1-secreted','IGHG2-secreted')

DotPlot(Bcell_LAIV.subset, group.by = 'cluster_annotated',dot.scale = 5,features=features) + RotatedAxis() +
  theme(axis.text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) +
  ylab('B-cell subsets')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_cluster_annotated_dotplot.pdf',sep = ''),width = 8, height = 4)

lineage19_table <- read.csv('tonsil_LAIV_120a_s_Bcell_log1p_cca_CycleRegressOut_cluster_lineage19_cluster_annotated.csv')
p <- match(colnames(Bcell.subset),lineage19_table$X)
lineage19_table <- lineage19_table[p,]
Bcell.subset$lineage19_annotation <- lineage19_table$lineage19_annotation

BCR_index <- data.frame(project_ID = Bcell.subset$project_ID, 
                        cell_index = gsub('_.*','',colnames(Bcell.subset)),
                        donorID = Bcell.subset$donor_ID, 
                        condition = Bcell.subset$condition, 
                        cluster_annotated = Bcell.subset$cluster_annotated,
                        major_clusters2 = Bcell.subset$major_clusters2,
                        major_clusters = Bcell.subset$major_clusters)
write.csv(BCR_index,'tonsil_LAIV_120a_s_Bcell_log1p_cca_CycleRegressOut_3annotations_012324.csv')

Bcell_LAIV.subset <- subset(Bcell.subset,stimulation != 'LAIV-')
Bcell_LAIV.subset[["raw"]] <- NULL
Bcell_LAIV.subset[["integrated"]] <- NULL
Bcell_LAIV.subset[["RNA"]] <- NULL
Bcell_LAIV.subset[["ADT"]] <- NULL
Bcell_LAIV.subset <- subset(Bcell_LAIV.subset,days != 'day02')
Bcell_LAIV.subset <- subset(Bcell_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
saveRDS(Bcell_LAIV.subset,'tonsil_LAIV_120a_s_lognorm_nonPB_LAIV.rds')

Bcell.subset[["raw"]] <- NULL
Bcell.subset[["integrated"]] <- NULL
Bcell.subset[["RNA"]] <- NULL
Bcell.subset[["ADT"]] <- NULL

saveRDS(Bcell.subset,'tonsil_LAIV_120a_s_lognorm_Bcell_LAIVns.rds')

temp <- subset(Bcell.subset,major_clusters2 == 'PB')
saveRDS(temp,'tonsil_LAIV_120a_s_LAIVns_major_cluster2_PB.rds')

temp <- subset(Bcell.subset,major_clusters2 == 'Memory B')
saveRDS(temp,'tonsil_LAIV_120a_s_LAIVns_major_cluster2_MemoryB.rds')

temp <- subset(Bcell.subset,major_clusters2 == 'Naive/Act. B')
saveRDS(temp,'tonsil_LAIV_120a_s_LAIVns_major_cluster2_NaiveOrActB.rds')

temp <- subset(Bcell.subset,major_clusters2 == 'GC-like B')
saveRDS(temp,'tonsil_LAIV_120a_s_LAIV_major_cluster2_GClikeB.rds')

temp <- subset(Bcell.subset,major_clusters2 == 'preGC/GC B')
saveRDS(temp,'tonsil_LAIV_120a_s_LAIVns_major_cluster2_preGCOrGCB.rds')

temp <- subset(Bcell.subset,if_PB != 'PB')
saveRDS(temp,'tonsil_LAIV_120a_s_nonPB_B_LAIVns.rds')

temp <- subset(Bcell_LAIV.subset,BCL6_pos == 'BCL6+B')
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6pos_nonPB_B_LAIV.rds')

temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (clusters %in% c(2,4,10)))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_memory_B_LAIV.rds')

temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (clusters %in% c(5,6)))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_activated_nonmemory_B_LAIV.rds')

temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (lineage19_annotation %in% c('unswitched act.','switched act. memory')))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_activated_B_LAIV.rds')

temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (differentiation != 'memory B'))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_nonmemory_B_LAIV.rds')

temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (lineage19_annotation == 'naive'))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_lineage19_naive_B_LAIV.rds')

temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (cluster_annotated == 'naive'))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_cluster_annotated_naive_B_LAIV.rds')

#CD27_CD38_clusters
temp <- subset(Bcell_LAIV.subset,AbCD27AbCD38_pos == 'AbCD27+AbCD38+B')
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_AbCD27_AbCD38_clusters_AbCD27+AbCD38+nonPB_B_LAIV.rds')
temp <- subset(Bcell_LAIV.subset,AbCD27AbCD38_pos == 'AbCD27+AbCD38-B')
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_AbCD27_AbCD38_clusters_AbCD27+AbCD38-nonPB_B_LAIV.rds')
temp <- subset(Bcell_LAIV.subset,AbCD27AbCD38_pos == 'AbCD27-AbCD38+B')
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_AbCD27_AbCD38_clusters_AbCD27-AbCD38+nonPB_B_LAIV.rds')
temp <- subset(Bcell_LAIV.subset,AbCD27AbCD38_pos == 'AbCD27-AbCD38-B')
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_AbCD27_AbCD38_clusters_AbCD27-AbCD38-nonPB_B_LAIV.rds')

#BCL6-CD27_CD38_clusters
temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (AbCD27AbCD38_pos == 'AbCD27+AbCD38+B'))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_AbCD27_AbCD38_clusters_AbCD27+AbCD38+nonPB_B_LAIV.rds')
temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (AbCD27AbCD38_pos == 'AbCD27+AbCD38-B'))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_AbCD27_AbCD38_clusters_AbCD27+AbCD38-nonPB_B_LAIV.rds')
temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (AbCD27AbCD38_pos == 'AbCD27-AbCD38+B'))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_AbCD27_AbCD38_clusters_AbCD27-AbCD38+nonPB_B_LAIV.rds')
temp <- subset(Bcell_LAIV.subset,(BCL6_pos == 'BCL6-B') & (AbCD27AbCD38_pos == 'AbCD27-AbCD38-B'))
# DimPlot(temp)
saveRDS(temp,'tonsil_LAIV_120a_s_BCL6neg_AbCD27_AbCD38_clusters_AbCD27-AbCD38-nonPB_B_LAIV.rds')

########### Bcell gating ###########################################
Bcell.subset$BCL6_pos <- 'BCL6-B'
Bcell.subset$BCL6_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['BCL6']) > 0)] <- 'BCL6+B'
Bcell.subset$BCL6_pos[Bcell.subset$if_PB == 'PB'] <- 'PB'

Bcell.subset$IGHG2_pos <- 'Bcell'
Bcell.subset$IGHG2_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['IGHG2-secreted']) > 5) & 
                              (!(Bcell.subset$BCL6_differentiation %in% c('PB','memory B')))] <- 'IGHG2+B'
Bcell.subset$IGHG2_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['IGHG2-secreted']) > 5) & 
                              (Bcell.subset$BCL6_differentiation %in% c('memory B'))] <- 'IGHG2+memory B'

Bcell.subset$MKI67_pos <- 'MKI67-B'
Bcell.subset$MKI67_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['MKI67']) > 0)] <- 'MKI67+B'

Bcell.subset$BCL6_gate <- 'BCL6-non-memory B'
Bcell.subset$BCL6_gate[Bcell.subset$differentiation == 'memory B'] <- 'BCL6-memory B'
Bcell.subset$BCL6_gate[(as.matrix(Bcell.subset@assays[["log1p"]]['BCL6']) > 0)] <- 'BCL6+B'
Bcell.subset$BCL6_gate[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AICDA_pos <- 'AICDA-B'
Bcell.subset$AICDA_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['AICDA']) > 0)] <- 'AICDA+B'
Bcell.subset$AICDA_pos[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AbCXCR5_pos <- 'AbCXCR5-B'
Bcell.subset$AbCXCR5_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0)] <- 'AbCXCR5+B'
Bcell.subset$AbCXCR5_pos[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AbCD21_pos <- 'AbCD21-B'
Bcell.subset$AbCD21_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD21']) > 0)] <- 'AbCD21+B'
Bcell.subset$AbCD21_pos[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AbCD83_pos <- 'AbCD83-B'
Bcell.subset$AbCD83_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD83']) > 0)] <- 'AbCD83+B'
Bcell.subset$AbCD83_pos[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AbCXCR4_pos <- 'AbCXCR4-B'
Bcell.subset$AbCXCR4_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD184']) > 0)] <- 'AbCXCR4+B'
Bcell.subset$AbCXCR4_pos[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AbCXCR5AbCD83_pos <- 'Bcell'
Bcell.subset$AbCXCR5AbCD83_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD83']) > 0)] <- 'AbCXCR5+AbCD83+B'
Bcell.subset$AbCXCR5AbCD83_pos[Bcell.subset$differentiation == 'PB'] <- 'Bcell'

Bcell.subset$AbCXCR5AbCXCR4_pos <- 'Bcell'
Bcell.subset$AbCXCR5AbCXCR4_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & 
                                       (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD184']) > 0)] <- 'AbCXCR5+AbCXCR4+B'
Bcell.subset$AbCXCR5AbCXCR4_pos[Bcell.subset$differentiation == 'PB'] <- 'Bcell'

Bcell.subset$AbCD27AbCD38_pos <- 'Bcell'
Bcell.subset$AbCD27AbCD38_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD27']) > 0) & 
                                     (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD38']) > 0)] <- 'GC B'
Bcell.subset$AbCD27AbCD38_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD27']) > 0) & 
                                (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD38']) < 0)] <- 'Memory B'
Bcell.subset$AbCD27AbCD38_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD27']) < 0) & 
                                (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD38']) > 0)] <- 'PreGC B'
Bcell.subset$AbCD27AbCD38_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD27']) < 0) & 
                                (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD38']) < 0)] <- 'Naive B'
Bcell.subset$AbCD27AbCD38_pos[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AbIgDAbIgM_pos <- 'Bcell'
Bcell.subset$AbIgDAbIgM_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-IgD']) > 0) & 
                                (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-IgM']) > 0)] <- 'AbIgD+AbIgM+B'
Bcell.subset$AbIgDAbIgM_pos[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AbCD24AbCD38_pos <- 'Bcell'
Bcell.subset$AbCD24AbCD38_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD38']) > 0) & 
                              (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD24']) > 0)] <- 'AbCD24AbCD38+B'
Bcell.subset$AbCD24AbCD38_pos[Bcell.subset$differentiation == 'PB'] <- 'Bcell'

Bcell.subset$COX_PFKP_pos <- 'COX-PFKP-B'
Bcell.subset$COX_PFKP_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['MT-CO1']) > 4) & 
                            (as.matrix(Bcell.subset@assays[["log1p"]]['PFKP']) <= 0)] <- 'COX+PFKP-B'
Bcell.subset$COX_PFKP_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['MT-CO1']) <= 4) & 
                            (as.matrix(Bcell.subset@assays[["log1p"]]['PFKP']) > 0)] <- 'COX-PFKP+B'
Bcell.subset$COX_PFKP_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['MT-CO1']) > 4) & 
                            (as.matrix(Bcell.subset@assays[["log1p"]]['PFKP']) > 0)] <- 'COX+PFKP+B'
Bcell.subset$COX_PFKP_pos[Bcell.subset$differentiation == 'PB'] <- 'PB'

Bcell.subset$AbCD11cTBX21_pos <- 'Bcell'
Bcell.subset$AbCD11cTBX21_pos[(as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-CD11c']) > 0) & 
                                 (as.matrix(Bcell.subset@assays[["log1p"]]['TBX21']) > 0)] <- 'AbCD11c+TBX21+B'
Bcell.subset$AbCD11cTBX21_pos[Bcell.subset$differentiation %in% c('PB')] <- 'PB'

Bcell.subset$GZMH_pos <- 'GZMH-B'
Bcell.subset$GZMH_pos[(as.matrix(Bcell.subset@assays[["log1p"]]['GZMH']) > 0)] <- 'GZMH+B'

Bcell.subset$BCL6_differentiation <- 'BCL6-B'
Bcell.subset$BCL6_differentiation[Bcell.subset$BCL6_pos == 'BCL6+B'] <- 'BCL6+B'
Bcell.subset$BCL6_differentiation[(Bcell.subset$BCL6_pos == 'BCL6-B') & (Bcell.subset$differentiation != 'memory B')] <- 'BCL6-B'
Bcell.subset$BCL6_differentiation[(Bcell.subset$BCL6_pos == 'BCL6-B') & (Bcell.subset$differentiation == 'memory B')] <- 'memory B'
Bcell.subset$BCL6_differentiation[(Bcell.subset$differentiation == 'PB')] <- 'PB'

Bcell.subset$AbCXCR5_differentiation <- 'AbCXCR5-B'
Bcell.subset$AbCXCR5_differentiation[Bcell.subset$AbCXCR5_pos == 'AbCXCR5+B'] <- 'AbCXCR5+B'
Bcell.subset$AbCXCR5_differentiation[(Bcell.subset$AbCXCR5_pos == 'AbCXCR5-B') & (Bcell.subset$differentiation != 'memory B')] <- 'AbCXCR5-B'
Bcell.subset$AbCXCR5_differentiation[(Bcell.subset$AbCXCR5_pos == 'AbCXCR5-B') & (Bcell.subset$differentiation == 'memory B')] <- 'memory B'
Bcell.subset$AbCXCR5_differentiation[(Bcell.subset$differentiation == 'PB')] <- 'PB'

if_NA <- is.na(Bcell.subset@assays[["RNA"]]@counts['FOXO1',])
Bcell.subset <- subset(Bcell.subset, if_NA != 1)
Bcell.subset$FOXO1_MYC_gating <- 'B cells'
Bcell.subset$FOXO1_MYC_gating[(as.matrix(Bcell.subset@assays[["RNA"]]@counts['FOXO1',]) > 0) & (as.matrix(Bcell.subset@assays[["RNA"]]@counts['MYC',]) > 0)] <- 'FOXO1+MYC+B'
Bcell.subset$FOXO1_MYC_gating[(Bcell.subset$differentiation == 'PB')] <- 'PB'

temp <- readRDS('tonsil_LAIV_120a_s_BCR_isotype.rds')
Bcell.subset$BCR_isotype <- 'unknown'
Bcell.subset$BCR_isotype <- temp$c_call[match(Bcell.subset$ProjectID_CellID,temp$ProjectID_CellID)]
Bcell.subset$BCR_isotype[is.na(Bcell.subset$BCR_isotype)] <- 'unknown'

Bcell.subset$naive <- 'Bcell'
Bcell.subset$naive[((as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-IgD']) > 0) & 
                      (as.matrix(Bcell.subset@assays[["integrated_scale"]]['ab-IgM']) < 0) &
                      (as.matrix(Bcell.subset@assays[["log1p"]]['IGHM-membrane']) > 0) &
                      (as.matrix(Bcell.subset@assays[["log1p"]]['IGHD-membrane']) > 0) &
                      (as.matrix(Bcell.subset@assays[["log1p"]]['TCL1A']) > 0) &
                      (as.matrix(Bcell.subset@assays[["log1p"]]['IL4R']) > 0)) ] <- 'naive B'
Bcell.subset$naive[Bcell.subset$differentiation == 'PB'] <- 'Bcell'


DimPlot(Bcell.subset, label = T, reduction = "umap",group.by = 'FOXO1_MYC_gating')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_umap_AbCD27AbCD38_pos.pdf',sep = ''),width = 6, height = 4)

temp <- subset(Bcell.subset, (if_PB != 'PB') & (stimulation != 'LAIV-'))
Idents(temp) <- temp$cluster_annotated

DimPlot(temp,reduction = "umap", split.by = "AbCD27AbCD38_pos", ncol = 2)
dev.print(pdf, 'tonsil_LAIV_120a_s_nonPB_B_log1p_cca_scale_CycleRegressOut_LAIV_umap_AbCD27AbCD38_pos_split.pdf',width = 10, height = 8.8)

Bcell_log1p_dataframe <- data.frame(t(as.matrix(Bcell.subset@assays[["log1p"]]@data)),check.names = F)
Bcell_integrated_scale_dataframe <- data.frame(t(as.matrix(Bcell.subset@assays[["integrated_scale"]]@data)),check.names = F)

Bcell_integrated_scale_dataframe$AbCD27AbCD38_pos <- Bcell.subset$AbCD27AbCD38_pos
ggplot(Bcell_integrated_scale_dataframe, aes(x = !!sym(c('ab-CD38')), y = !!sym(c('ab-CD27')))) + ggtitle('Bcell cells') +#!!sym(c('TBX21'))
  geom_point(alpha = 0.1,size = 0.5,aes(color = !!sym(c('AbCD27AbCD38_pos')))) +
  # geom_point(alpha = 0.1,size = 0.5) +
  # scale_colour_manual(values = c("naive Bcell" = "#F8766D", "Bcell" = '#00BFC4'))+
  # scale_colour_manual(values = c("Bcell5RA+CD83+ nonPB B" = "gold", "CD184-CD83- nonPB B" = "blue","CD184+CD83- nonPB B" = 'cyan',"CD184-CD83+ nonPB B" = 'magenta')) +
  geom_density_2d(color = 'black')# +
# ylim(0,8) + xlim(0,8)
# geom_hline(aes(yintercept=3.25,color = 'red')) + guides(color = 'none') +
# geom_vline(aes(xintercept=4.1,color = 'red')) + guides(color = 'none')
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_AbCD27AbCD38_pos_contour_integrated_scale.pdf',width = 7, height = 4.6)

DefaultAssay(Bcell.subset) <- 'lognorm'
FeatureScatter(Bcell.subset, feature1 = "MT-CO1", feature2 = "PFKP") 
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_scatter_abBcell5RA_abBcell5RO.pdf',width = 700, height = 4.34)

######### add gene module and visualize ###########################
library("xlsx")
gene_module_list <- read.xlsx('../GO_terms/GO0001666_response_to_hypoxia.xlsx', sheetName = "Sheet1",header = F)
gene_module_list <- gene_module_list[gene_module_list$X5 == 'cellular response to hypoxia',]
gene_module_list <- gene_module_list[toupper(gene_module_list$X1) %in% toupper(all.markers),]
gene_module_list <- toupper(unique(gene_module_list$X1))
# gene_module_list <- read.table('../GO_terms/GO0061621_canonical_glycolysis.txt',sep = "\t",header = T, quote = "#",check.names = FALSE)
# gene_module_list <- gene_module_list[toupper(gene_module_list$Symbol) %in% toupper(all.markers),]
# gene_module_list <- toupper(unique(gene_module_list$Symbol))

GO_term <- 'HIF1a_signaling_pathway'
Bcell.subset <- AddModuleScore(Bcell.subset,features = list(gene_module_list),name = GO_term,ctrl = 25)
FeaturePlot(Bcell.subset, features = paste(GO_term,"1",sep = ''), label = F, repel = TRUE,min.cutoff = 0) 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_umap_',GO_term,'.pdf',sep = ''),width = 5, height = 4)

# FindMarkers(TB_shrink.seurat, group.by = "Group",ident.1 = 'RSTR',ident.2 = 'LTBI',features = c("early_differentiation1"))
wilcox_test <- wilcox.test(Bcell.subset$early_differentiation1[TB_shrink.seurat$Group == 'RSTR'],
                           Bcell.subset$early_differentiation1[TB_shrink.seurat$Group == 'LTBI'])




##### cytokine analysis ###########################################
# cytokine_global_table <- read_excel('../Global landscape of cytokines Supplementary Table S1.xlsx', sheet = "Cytokines")
# cytokine_global_list <- unique(cytokine_global_table$`HGNC symbol`)
# cytokine_global_list <- c(cytokine_global_list,'GZMB','PRF1','MIF')
# cytokine_Rhapsody_list <- all.markers[all.markers %in% cytokine_global_list]
# 
# for (cytokine_name in cytokine_Rhapsody_list){
#   graphics.off()
#   plot <- FeaturePlot(Bcell.subset, features = c(cytokine_name),pt.size = 0.2, ncol = 1, sort.cell = TRUE,min.cutoff = 0)#, max.cutoff = 15)
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_umap_',cytokine_name,'.pdf',sep = ''),width = 5, height = 4.3)
#   
# }

##### all Bcells: cell cluster fraction visualization #####################
flow_count_table <- read.xlsx('tonsil_LAIV_120a_s_flow_cell_count.xlsx', sheetName = "Sheet1")
flow_count_table <- flow_count_table[,!colnames(flow_count_table) == '...1']
flow_count_table$condition <- gsub('day0 LAIV-','day0',flow_count_table$condition)
flow_count_table$B_cells[grepl('IMD170',flow_count_table$donor_ID)] <- 0

Bcell.database <- data.frame(cbind(Bcell.subset$donor_ID,as.character(Bcell.subset$subisotype),Bcell.subset$condition,Bcell.subset$days,Bcell.subset$stimulation,Bcell.subset$age_group,Bcell.subset$day_group,Bcell.subset$age))
colnames(Bcell.database) <- c('donor_ID','subisotype','condition','days','stimulation','age_group','day_group','age')
# sapply(Bcell.database,class)
Bcell_cluster_condition_donor_count <- dplyr::count(Bcell.database, donor_ID, subisotype,condition, days,stimulation,age_group,day_group,age)
Bcell_cluster_condition_donor_count <- Bcell_cluster_condition_donor_count %>% group_by(donor_ID, condition, days,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
Bcell_cluster_condition_donor_count$days <- as.numeric(substring(Bcell_cluster_condition_donor_count$days,first = 4,last = 5))
Bcell_cluster_condition_donor_count$percentage <- Bcell_cluster_condition_donor_count$n/Bcell_cluster_condition_donor_count$total*100
Bcell_cluster_condition_donor_count$age <- as.numeric(Bcell_cluster_condition_donor_count$age)
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(Bcell.subset$donor_ID)){
  temp_condition_list <- sort(unique(Bcell.subset$condition[Bcell.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$donor_ID == donor_name) & (Bcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Bcell.subset$subisotype)){
      if_row <- ((Bcell_cluster_condition_donor_count$donor_ID == donor_name) & (Bcell_cluster_condition_donor_count$subisotype == cluster_name) & (Bcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount$subisotype <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        Bcell_cluster_condition_donor_count[nrow(Bcell_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(Bcell.database)

# write.csv(Bcell_cluster_condition_donor_count, "tonsil_Rhapsody_Bcell_cluster_percentage.csv")
# two figures:
# 1. only responders
# 2. all donors but only at day0
Bcell_cluster_condition_donor_count_day0 <- Bcell_cluster_condition_donor_count[Bcell_cluster_condition_donor_count$days == '0',]
Bcell_cluster_condition_donor_count_double_day0 <- Bcell_cluster_condition_donor_count
temp <- Bcell_cluster_condition_donor_count_double_day0[Bcell_cluster_condition_donor_count_double_day0$days == '0',]
Bcell_cluster_condition_donor_count_double_day0$stimulation[Bcell_cluster_condition_donor_count_double_day0$stimulation == 'day0'] <- 'LAIV+'
temp$stimulation[temp$stimulation == 'day0'] <- 'LAIV-'
Bcell_cluster_condition_donor_count_double_day0 <- rbind(Bcell_cluster_condition_donor_count_double_day0,temp)
Bcell_cluster_condition_donor_count_double_day0 <- Bcell_cluster_condition_donor_count_double_day0[Bcell_cluster_condition_donor_count_double_day0$donor_ID != "02yrs M IMD085",]
p <- match(Bcell_cluster_condition_donor_count_double_day0$condition, flow_count_table$condition)
temp <- flow_count_table[p,]
Bcell_cluster_condition_donor_count_double_day0$Bcell_count <- temp$B_cells
Bcell_cluster_condition_donor_count_LAIV <- Bcell_cluster_condition_donor_count_double_day0[Bcell_cluster_condition_donor_count_double_day0$stimulation == 'LAIV+',]
annotate_cluster_list <- unique(Bcell.subset$subisotype)#'AbCD278+B'#
annotate_cluster_list <- annotate_cluster_list[grepl('\\+',annotate_cluster_list)]# | grepl('CD278\\+PD1\\-CD4',annotate_cluster_list) | grepl('CD278\\-PD1\\+CD4',annotate_cluster_list)]

for (cluster_name in annotate_cluster_list)
{
  print(cluster_name)
  temp_Bcell_cluster <- Bcell_cluster_condition_donor_count_double_day0[Bcell_cluster_condition_donor_count_double_day0$subisotype == cluster_name,]
  plot <- ggplot(temp_Bcell_cluster,aes(x=days, y=percentage,color = stimulation, shape = donor_ID)) +  ggtitle(cluster_name) +
    # facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
    facet_wrap(~ age_group + donor_ID, scales = "fixed") +
    geom_point(size = 1)+ geom_line(aes(group = interaction(donor_ID,stimulation))) + 
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('Bcell %') +
    scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
    # scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) +
    scale_shape_manual(values=1:length(unique(Bcell_cluster_condition_donor_count_double_day0$donor_ID))) +
    # geom_rect(aes(xmin=4, xmax=10, ymin=0, ymax=Inf),alpha = 0.1,fill = 'gray',color = 'gray') +
    theme_bw()
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_log1p_scale_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 7, height = 5)
  temp_Bcell_cluster$cluster_count <- temp_Bcell_cluster$Bcell_count*temp_Bcell_cluster$percentage/100
  plot <- ggplot(temp_Bcell_cluster,aes(x=days, y=cluster_count,color = stimulation, shape = donor_ID)) +  ggtitle(cluster_name) +
    # facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
    facet_wrap(~ age_group + donor_ID, scales = "fixed") +
    geom_point(size = 1)+ geom_line(aes(group = interaction(donor_ID,stimulation))) + 
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('#cells') +
    scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
    # scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) +
    scale_shape_manual(values=1:length(unique(Bcell_cluster_condition_donor_count_double_day0$donor_ID))) +
    # geom_rect(aes(xmin=4, xmax=10, ymin=0, ymax=Inf),alpha = 0.1,fill = 'gray',color = 'gray') +
    theme_bw()
    # annotate('rect', xmin=4, xmax=6, ymin=0, ymax=1e6, alpha=.5, fill='gray')
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_log1p_scale_',cluster_name,'_condition_count2.pdf',sep = ''),width = 8, height = 5)
}

######## cell fraction comparison visualization ###########################
age_group_list <- sort(unique(Bcell.subset$age_group))
LAIV_day_group_dict <- list()
LAIV_day_group_dict[['day0']] <- c(0)
LAIV_day_group_dict[['day04']] <- c(4)
LAIV_day_group_dict[['day06-07']] <- c(6,7)
LAIV_day_group_dict[['day07-08']] <- c(7,8)
LAIV_day_group_dict[['day10']] <- c(10)
LAIV_day_group_dict[['day12-14']] <- c(12,14)

age_pair_list <- data.frame(age_group1 = c('02-03yrs','02-03yrs','07-09yrs'),
                            age_group2 = c('07-09yrs','27-39yrs','27-39yrs'))

LAIV_day_group_list <- names(LAIV_day_group_dict)
cluster_ttest_age_group <- expand.grid(cluster = annotate_cluster_list, day_group = LAIV_day_group_list, age_group1 = c('02-03yrs','07-09yrs'),age_group2 = c('07-09yrs','27-39yrs'))
cluster_ttest_age_group$x <- NaN
cluster_ttest_age_group$y <- NaN
cluster_ttest_age_group$pval <- NaN
cluster_ttest_age_group$age_group1 <- as.character(cluster_ttest_age_group$age_group1)
cluster_ttest_age_group$age_group2 <- as.character(cluster_ttest_age_group$age_group2)

cluster_ttest_age_group <- cluster_ttest_age_group[as.character(cluster_ttest_age_group$age_group1) != as.character(cluster_ttest_age_group$age_group2),]
cluster_correlate_age <- expand.grid(cluster = annotate_cluster_list, day_group = LAIV_day_group_list)
cluster_correlate_age$r <- NaN
cluster_correlate_age$pval <- NaN

library(dplyr)
for (cluster_name in annotate_cluster_list) {
  print(cluster_name)
  graphics.off()
  temp <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$subisotype == cluster_name) & 
                                                (Bcell_cluster_condition_donor_count$stimulation != 'LAIV-'),]
  temp <- temp[temp$donor_ID != '02yrs M IMD085',]
  if ((!all(temp$percentage == 0)) & (!all(temp$percentage == 100))) {
    # One should do chi-square test!!
    for (day_group_name in LAIV_day_group_list) {
      temp_day0 <- temp[temp$days %in% LAIV_day_group_dict[[day_group_name]],]
      # mean
      if (!all(temp_day0$percentage == 0)) {
        for (age_group_index in c(1:3)) {
          ttest_result <- t.test(temp_day0$percentage[temp_day0$age_group %in% age_pair_list$age_group1[age_group_index]],
                                 temp_day0$percentage[temp_day0$age_group %in% age_pair_list$age_group2[age_group_index]])
          row_select <- (cluster_ttest_age_group$cluster == cluster_name) & (cluster_ttest_age_group$age_group1 == age_pair_list$age_group1[age_group_index]) & 
            (cluster_ttest_age_group$age_group2 == age_pair_list$age_group2[age_group_index]) & (cluster_ttest_age_group$day_group == day_group_name)
          cluster_ttest_age_group$x[row_select] <- ttest_result[["estimate"]][["mean of x"]]
          cluster_ttest_age_group$y[row_select] <- ttest_result[["estimate"]][["mean of y"]]
          cluster_ttest_age_group$pval[row_select] <- ttest_result$p.value
        }
        row_select <- (cluster_correlate_age$cluster == cluster_name) & (cluster_correlate_age$day_group == day_group_name)
        cor_test <- cor.test(temp_day0$percentage, temp_day0$age)
        cluster_correlate_age$r[row_select] <- cor_test[["estimate"]][["cor"]]
        cluster_correlate_age$pval[row_select] <- cor_test$p.value  
      }
    }
  }
}

# cluster_ttest_age_group$pval_adj <- p.adjust(cluster_ttest_age_group$pval,method = 'fdr')
# cluster_correlate_age$pval_adj <- p.adjust(cluster_correlate_age$pval,method = 'fdr')
cluster_ttest_age_group <- cluster_ttest_age_group %>% arrange(pval)
cluster_ttest_age_group$FC <- cluster_ttest_age_group$y/cluster_ttest_age_group$x
cluster_ttest_age_group$log2FC <- log2(cluster_ttest_age_group$FC)
cluster_correlate_age <- cluster_correlate_age %>% arrange(pval)
write.xlsx(cluster_ttest_age_group,'tonsil_LAIV_120a_s_Bcell_subisotype_age_group_ttest.xlsx')
write.xlsx(cluster_correlate_age,'tonsil_LAIV_120a_s_Bcell_subisotype_age_corr.xlsx')

####### visualization ####################################
cluster_name <- 'AbCD11c+TBX21+B'
day_group_name <- 'day10'
age_group_name <- '02-03yrs'
marker_name <- 'ab-CD11c'
x_lab <- 1.7
plot_corr <- 0
graphics.off()
temp_donorID_batch_list <- unlist(batch_donorID_list[which(sapply(batch_marker_list, function(x) (marker_name %in% x)))])

temp <- Bcell_cluster_condition_donor_count[Bcell_cluster_condition_donor_count$AbCD11cTBX21_pos == cluster_name,]
rownames(temp) <- temp$condition
if (day_group_name == 'day0') {
  temp_day0 <- temp[temp$days == '0',]
  temp_day0 <- temp_day0[temp_day0$donor_ID != '02yrs M IMD085',] %>% arrange(donor_ID)

} else {
  # paired diff comparison
  condition_paired_cluster_table <- condition_paired_table
  condition_paired_cluster_table$LAIV_cluster <- unlist(temp[condition_paired_table$LAIV,][,'percentage'])
  condition_paired_cluster_table$ns_cluster <- unlist(temp[condition_paired_table$ns,][,'percentage'])
  condition_paired_cluster_table$diff_cluster <- condition_paired_cluster_table$LAIV_cluster - condition_paired_cluster_table$ns_cluster
  
  ggplot(condition_paired_cluster_table,aes(x=days, y=diff_cluster,shape = donor_ID)) +  
    ggtitle(cluster_name) +
    geom_hline(yintercept = 0,linetype=2) + 
    # facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
    facet_wrap(~ age_group + donor_ID, scales = "fixed") +
    geom_point(size = 1)+ geom_line() + 
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('diff in % B cells\nbetween LAIV+ and LAIV-') +
    # scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) +
    scale_shape_manual(values=1:length(unique(condition_paired_cluster_table$donor_ID))) +
    # geom_rect(aes(xmin=4, xmax=10, ymin=0, ymax=Inf),alpha = 0.1,fill = 'gray',color = 'gray') +
    theme_bw()
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_AbCD11cTBX21_pos_diff_',cluster_name,'.pdf',sep = ''),width = 7, height = 5)
  
  temp_day_nonpaired <- condition_paired_cluster_table[condition_paired_cluster_table$donor_ID %in% temp_donorID_batch_list,]
  temp_day_nonpaired <- temp_day_nonpaired[(temp_day_nonpaired$days %in% day_group_dict[[day_group_name]]),]
  if (plot_corr == 1) {
    cor_test <- cor.test(temp_day_nonpaired$diff_cluster,temp_day_nonpaired$age)
    plot <- ggplot(temp_day_nonpaired, aes(y=diff_cluster,x=age)) +
    geom_point() + RotatedAxis() +
    ylab('diff in % B cells\nbetween LAIV+ and LAIV-') +
    geom_smooth(method='lm', formula= y~x, color = 'red') +
    ggtitle(paste(cluster_name, day_group_name)) +
    theme_bw() +
    annotate("text", x=(temp_day_nonpaired$age[which.max(temp_day_nonpaired$diff_cluster)] + 25)%%40, 
             y=mean(sort(temp_day_nonpaired$diff_cluster)[(length(temp_day_nonpaired$diff_cluster) - 1):length(temp_day_nonpaired$diff_cluster)]),
             label = paste("r =",sprintf(cor_test$estimate, fmt = '%#.2f')),size = 6)
    print(plot)
    dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_AbCD11cTBX21_pos_age_corr_',cluster_name,'_',day_group_name,'.pdf',sep = ''),width = 3, height = 3)
  } else {
    ttest_result <- t.test(temp_day_nonpaired$diff_cluster[temp_day_nonpaired$age_group == age_group_name],
                           temp_day_nonpaired$diff_cluster[temp_day_nonpaired$age_group != age_group_name])
    plot <- ggplot(temp_day_nonpaired, aes(y=diff_cluster,x=age_group)) +
      geom_point() + RotatedAxis() +
      geom_hline(yintercept = 0,linetype=2) + 
      ylab('diff in % B cells\nbetween LAIV+ and LAIV-') +
      # geom_smooth(method='lm', formula= y~x, color = 'red') +
      ggtitle(paste(cluster_name, day_group_name)) +
      theme_bw() +
      annotate("text", x = x_lab,
               y=max(temp_day_nonpaired$diff_cluster) + 0.5*abs(max(temp_day_nonpaired$diff_cluster)),
               label = paste("p =",sprintf(ttest_result$p.value,fmt = '%#.3f')),size = 5)
    print(plot)
    dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_AbCD11cTBX21_pos_age_group_',cluster_name,'_',day_group_name,'_',age_group_name,'.pdf',sep = ''),width = 3, height = 3)
  }
}      
 


####### subcluster dynamic ranking ###############################
library(ggridges)
library(reshape2)
Bcell_subcluster_mean_dynamics <- data.frame()
for (cluster_name in unique(Bcell_cluster_condition_donor_count_double_day0$subclusters)) {
  temp <- Bcell_cluster_condition_donor_count_double_day0[Bcell_cluster_condition_donor_count_double_day0$subclusters == cluster_name,]
  graphics.off()
  temp_model <- loess(percentage ~ days, data=temp)
  temp_xseq <- c(0:14)
  temp_pred <- predict(temp_model, newdata = data.frame(days = temp_xseq), se=TRUE)
  temp_y = temp_pred$fit
  temp_ci <- temp_pred$se.fit * qt(0.95 / 2 + .5, temp_pred$df)
  temp_ymin = temp_y - temp_ci
  temp_ymax = temp_y + temp_ci
  temp_loess.DF <- data.frame(days = temp_xseq, percentage = temp_y, cluster = cluster_name, ymin = temp_ymin, ymax = temp_ymax, se = temp_pred$se.fit)
  temp_loess.DF$percentage[temp_loess.DF$percentage < 0] <- 0
  temp_mean <- data.frame(t(temp_loess.DF$percentage))
  colnames(temp_mean) <- temp_xseq
  rownames(temp_mean) <- cluster_name
  Bcell_subcluster_mean_dynamics <- rbind(Bcell_subcluster_mean_dynamics,temp_mean)
#   plot <- ggplot(temp,aes(x=days,y=percentage)) +  ggtitle(cluster_name) +
#     geom_point(aes(x=days,y=percentage, color = stimulation, shape = donor_ID), size = 5)+
#     geom_line(aes(x=days,y=percentage, color = stimulation, shape = donor_ID,group = interaction(donor_ID,stimulation))) +
#     theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
#     ylab('Bcell %') +
#     scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
#     scale_shape_manual(values=1:length(unique(Bcell_cluster_condition_donor_count_double_day0$donor_ID))) +
#     geom_smooth(aes_auto(temp_loess.DF), data=temp_loess.DF, stat="identity")
#   print(plot)
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_smooth_pattern_subclusters_',cluster_name,'_2.pdf',sep = ''),width = 7, height = 4)
# # this figure cannot eliminate points negative since it is doing the fitting automatically
#   plot <- ggplot(temp,aes(x=days,y=percentage)) +  ggtitle(cluster_name) +
#     geom_point(aes(x=days,y=percentage, color = stimulation, shape = donor_ID), size = 5)+
#     geom_line(aes(x=days,y=percentage, color = stimulation, shape = donor_ID,group = interaction(donor_ID,stimulation))) +
#     theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
#     ylab('Bcell %') +
#     scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
#     scale_shape_manual(values=1:length(unique(Bcell_cluster_condition_donor_count_double_day0$donor_ID)))+
#     stat_smooth(method = "loess", formula = y ~ x, size = 1)
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_smooth_pattern_subclusters_',cluster_name,'.pdf',sep = ''),width = 7, height = 4)

}

Bcell_subcluster_mean_dynamics_norm <- Bcell_subcluster_mean_dynamics/apply(Bcell_subcluster_mean_dynamics, 1, max)

# find the peak
temp <- data.frame(which(Bcell_subcluster_mean_dynamics_norm == 1,arr.ind=TRUE))
temp$subcluster <- rownames(temp)
temp$row <- NULL
temp$day_peak <- temp$col
temp <- temp %>% arrange(day_peak)


Bcell_subcluster_mean_dynamics_norm$subcluster <- rownames(Bcell_subcluster_mean_dynamics_norm)
temp2 <- melt(Bcell_subcluster_mean_dynamics_norm,id = 'subcluster', value.name = 'smoothed_percentage_norm') 
temp2$subcluster <- factor(temp2$subcluster,levels = temp$subcluster)
temp2$days <- as.numeric(temp2$variable) - 1
ggplot(temp2, aes(x=days, y=subcluster, height = smoothed_percentage_norm, group = subcluster, fill=subcluster)) + geom_density_ridges(stat = "identity", scale = 3)
# ggplot(temp2, aes(x=days,y=smoothed_percentage_norm,color = subcluster)) +  geom_line() 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_smooth_pattern_subclusters_norm_ordered.pdf',sep = ''),width = 7, height = 4)


Bcell_subcluster_stimulation_mean_dynamics <- data.frame()
for (cluster_name in unique(Bcell_cluster_condition_donor_count_double_day0$subclusters)) {
  temp <- Bcell_cluster_condition_donor_count_double_day0[(Bcell_cluster_condition_donor_count_double_day0$subclusters == cluster_name),]
  graphics.off()
  temp1 <- Bcell_cluster_condition_donor_count_double_day0[(Bcell_cluster_condition_donor_count_double_day0$subclusters == cluster_name) & (Bcell_cluster_condition_donor_count_double_day0$stimulation == 'LAIV+'),]
  temp2 <- Bcell_cluster_condition_donor_count_double_day0[(Bcell_cluster_condition_donor_count_double_day0$subclusters == cluster_name) & (Bcell_cluster_condition_donor_count_double_day0$stimulation == 'LAIV-'),]
  temp_model <- loess(percentage ~ days, data=temp1)
  temp_xseq <- c(0:14)
  temp_pred <- predict(temp_model, newdata = data.frame(days = temp_xseq), se=TRUE)
  temp_y = temp_pred$fit
  temp_ci <- temp_pred$se.fit * qt(0.95 / 2 + .5, temp_pred$df)
  temp_ymin = temp_y - temp_ci
  temp_ymax = temp_y + temp_ci
  temp_loess.DF <- data.frame(days = temp_xseq, percentage = temp_y, cluster = cluster_name, ymin = temp_ymin, ymax = temp_ymax, se = temp_pred$se.fit)
  temp_loess.DF$percentage[temp_loess.DF$percentage < 0] <- 0
  temp1 <- data.frame(t(temp_loess.DF$percentage))
  colnames(temp1) <- temp_xseq
  rownames(temp1) <- paste(cluster_name,'LAIV')

  temp_model <- loess(percentage ~ days, data=temp2)
  temp_xseq <- c(0:14)
  temp_pred <- predict(temp_model, newdata = data.frame(days = temp_xseq), se=TRUE)
  temp_y = temp_pred$fit
  temp_ci <- temp_pred$se.fit * qt(0.95 / 2 + .5, temp_pred$df)
  temp_ymin = temp_y - temp_ci
  temp_ymax = temp_y + temp_ci
  temp_loess.DF <- data.frame(days = temp_xseq, percentage = temp_y, cluster = cluster_name, ymin = temp_ymin, ymax = temp_ymax, se = temp_pred$se.fit)
  temp_loess.DF$percentage[temp_loess.DF$percentage < 0] <- 0
  temp2 <- data.frame(t(temp_loess.DF$percentage))
  colnames(temp2) <- temp_xseq
  rownames(temp2) <- paste(cluster_name,'ns')
  if (max(temp1) > max(temp2)){
    Bcell_subcluster_stimulation_mean_dynamics <- rbind(Bcell_subcluster_stimulation_mean_dynamics,temp1)
  } else {
    Bcell_subcluster_stimulation_mean_dynamics <- rbind(Bcell_subcluster_stimulation_mean_dynamics,temp2)
  }
  
  # # this figure cannot eliminate points negative since it is doing the fitting automatically
  #   plot <- ggplot(temp,aes(x=days,y=percentage, color = stimulation)) +  ggtitle(cluster_name) +
  #     geom_point(aes(x=days,y=percentage, shape = donor_ID), size = 5)+
  #     geom_line(aes(x=days,y=percentage, shape = donor_ID,group = interaction(donor_ID,stimulation))) +
  #     theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  #     ylab('Bcell %') +
  #     scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "grey45")) +
  #     scale_shape_manual(values=1:length(unique(Bcell_cluster_condition_donor_count_double_day0$donor_ID)))+
  #     stat_smooth(method = "loess", formula = y ~ x, size = 1) 
  #   print(plot)
  #   dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_smooth_pattern_subclusters_',cluster_name,'_LAIV_ns.pdf',sep = ''),width = 7, height = 4)

}

Bcell_subcluster_stimulation_mean_dynamics_norm <- Bcell_subcluster_stimulation_mean_dynamics/apply(Bcell_subcluster_stimulation_mean_dynamics, 1, max)
# sum(Bcell_subcluster_stimulation_mean_dynamics_norm > 0.95)
# sum(Bcell_subcluster_stimulation_mean_dynamics_norm == 1)
# temp <- data.frame(melt(Bcell_subcluster_stimulation_mean_dynamics_norm))
# temp <- temp %>% arrange(value)
# temp$variable <- c(1:dim(temp)[1])
# ggplot(temp,aes(x = variable,y = value)) + geom_point() + geom_hline(yintercept = 0.95) # no clear cutoff


# find the peak, check whether LAIV % or ns% dominates, then choose the peak from the dominant
temp <- data.frame(which(Bcell_subcluster_stimulation_mean_dynamics_norm == 1,arr.ind=TRUE))
temp$subcluster_stimulation <- rownames(temp)
temp$subcluster <- gsub(' .*','',temp$subcluster_stimulation)
temp$stimulation <- gsub('.* ','',temp$subcluster_stimulation)
temp$row <- NULL
temp$day_peak <- temp$col
temp <- temp %>% arrange(day_peak)
temp2 <- dcast(temp,   subcluster ~ stimulation, value.var = 'day_peak')
temp2$trend <- apply(temp2[,c('LAIV','ns')], 1, min)
write.xlsx(temp,'tonsil_LAIV_120a_s_Bcell_smooth_pattern_subclusters_LAIV_ns_dominant.xlsx')


Bcell_subcluster_mean_cumsum <- data.frame(t(apply(Bcell_subcluster_mean_dynamics, 1, cumsum))) # smoothed cumulative sum
Bcell_subcluster_mean_cumsum_norm <- Bcell_subcluster_mean_cumsum/Bcell_subcluster_mean_cumsum$X14
Bcell_subcluster_mean_cumsum_norm$subcluster <- rownames(Bcell_subcluster_mean_cumsum_norm)
temp <- melt(Bcell_subcluster_mean_cumsum_norm,id = 'subcluster', value.name = 'smoothed_percentage_cumsum') 
temp$days <- as.numeric(substring(temp$variable,first = 2, last = 4))

ggplot(temp, aes(x=days,y=smoothed_percentage_cumsum,color = subcluster)) +  geom_line() 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_smooth_pattern_subclusters_cumsum_norm.pdf',sep = ''),width = 7, height = 4)







library(hrbrthemes)
library(GGally)
library(viridis)


Bcell_cluster_condition_donor_count_LAIV

Bcell_cluster_condition_donor_count_ns <- Bcell_cluster_condition_donor_count_double_day0[Bcell_cluster_condition_donor_count_double_day0$stimulation == 'LAIV-',]


###### GC subset #########################
GC.subset <- subset(Bcell.subset,differentiation == 'GC B')
saveRDS(GC.subset,'tonsil_LAIV_120a_s_Bcell_GC.rds')
GC.subset <- RunPCA(GC.subset, features = all.markers, npcs = 50)
GC.subset <- RunUMAP(GC.subset, reduction = "pca", dims = 1:n_pca_selected)
GC.subset <- FindNeighbors(GC.subset, reduction = "pca", dims = 1:n_pca_selected)
GC.subset <- FindClusters(GC.subset, resolution = 2)

DimPlot(GC.subset, label = TRUE, reduction = "umap",repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_GC_umap_subcluster.pdf',width = 700, height = 4.34)

BCL6pos.subset <- subset(Bcell_LAIV.subset,BCL6_pos == 'BCL6+B')
saveRDS(BCL6pos.subset,'tonsil_LAIV_120a_s_Bcell_BCL6pos_nonPB.rds')
rm(BCL6pos.subset)

BCL6pos.subset <- subset(Bcell.subset,AbCXCR5_pos == 'AbCXCR5+B')
saveRDS(BCL6pos.subset,'tonsil_LAIV_120a_s_Bcell_AbCXCR5pos_nonPB.rds')
rm(BCL6pos.subset)

naive.subset <- subset(Bcell.subset,differentiation == 'naive B')
saveRDS(naive.subset,'tonsil_LAIV_120a_s_Bcell_naive.rds')
rm(naive.subset)

BCL6neg_memory.subset <- subset(Bcell.subset,(BCL6_pos == 'BCL6-B') & (differentiation == 'memory B'))
saveRDS(BCL6neg_memory.subset,'tonsil_LAIV_120a_s_Bcell_BCL6neg_memory_B.rds')

BCL6neg_memory.subset <- subset(Bcell.subset,CD27CD38_cluster == 'CD27+CD38+B')
saveRDS(BCL6neg_memory.subset,'tonsil_LAIV_120a_s_Bcell_CD27CD38_cluster_CD27posCD38pos_B.rds')

BCL6neg_memory.subset <- subset(Bcell.subset,CD27CD38_cluster == 'memory B')
saveRDS(BCL6neg_memory.subset,'tonsil_LAIV_120a_s_Bcell_CD27CD38_cluster_memory_B.rds')

BCL6neg_memory.subset <- subset(Bcell.subset,CD27CD38_cluster == 'BCL6+B')
saveRDS(BCL6neg_memory.subset,'tonsil_LAIV_120a_s_Bcell_CD27CD38_cluster_BCL6pos_B.rds')

BCL6neg.subset <- subset(Bcell.subset,AbCXCR5_pos == 'AbCXCR5-B')
saveRDS(BCL6neg.subset,'tonsil_LAIV_120a_s_Bcell_AbCXCR5neg_nonmemory_nonPB_B.rds')
rm(BCL6neg.subset)

######## nonmemory nonPB B cells ####################
nonmemory_nonPB.subset <- subset(Bcell.subset, (differentiation != 'memory B') & (differentiation != 'PB'))
DefaultAssay(nonmemory_nonPB.subset) <- "integrated_scale"
nonmemory_nonPB.subset <- ScaleData(nonmemory_nonPB.subset,features = all.markers)

nonmemory_nonPB.subset <- RunPCA(nonmemory_nonPB.subset, features = all.markers, npcs = 50)
# ElbowPlot(object = Bcell.subset,ndims = 50) + theme(axis.text = element_text(size = 20))
# dev.print(pdf, 'tonsil_LAIV_120a_scleaned_log1p_cca_scale_CycleRegressOut_PC_elbow.pdf',width = 700, height = 4.34)
n_pca_selected <- 30
nonmemory_nonPB.subset <- RunUMAP(nonmemory_nonPB.subset, reduction = "pca", dims = 1:n_pca_selected)
nonmemory_nonPB.subset <- FindNeighbors(nonmemory_nonPB.subset, reduction = "pca", dims = 1:n_pca_selected)
nonmemory_nonPB.subset <- FindClusters(nonmemory_nonPB.subset, resolution = 0.5)
# nonmemory_nonPB.subset$subclusters <- nonmemory_nonPB.subset$seurat_clusters
nonmemory_nonPB.subset$clusters <- nonmemory_nonPB.subset$seurat_clusters
DimPlot(nonmemory_nonPB.subset, label = TRUE, reduction = "umap",group.by = 'clusters',repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_nonmemory_nonPB_umap_cluster.pdf',width = 570, height = 400)

FeaturePlot(nonmemory_nonPB.subset, features = c("XBP1","PRDM1"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_cleaned_umap_XBP1_PRDM1.pdf',width = 1200, height = 4.34)

nonmemory_nonPB.markers_cluster <- FindAllMarkers(nonmemory_nonPB.subset)
nonmemory_nonPB.markers_cluster <- nonmemory_nonPB.markers_cluster %>% filter(p_val_adj <= 0.05)
p <- match(nonmemory_nonPB.markers_cluster$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
nonmemory_nonPB.markers_cluster$ENTREZID <- temp$ENTREZID
write.xlsx(nonmemory_nonPB.markers_cluster, "tonsil_LAIV_120a_s_Bcell_nonmemory_nonPB_cluster_gene.xlsx")
for (cluster_name in sort(unique(nonmemory_nonPB.subset$clusters))){
  temp <- nonmemory_nonPB.markers_cluster[nonmemory_nonPB.markers_cluster$cluster == cluster_name,]
  temp <- temp %>% arrange(desc(avg_log2FC))
  write.xlsx(temp, "tonsil_LAIV_120a_s_Bcell_nonmemory_nonPB_cluster_gene.xlsx",sheetName = paste('cluster',cluster_name,sep = ''),append = T)
}
Bcell.markers_cluster.top10 <- Bcell.markers_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Bcell.markers_cluster.top10, "tonsil_LAIV_120a_s_Bcell_cluster_gene_top10.csv")


temp <- FindMarkers(Bcell.subset,group.by = 'GZMH_pos',ident.1 = 'GZMH+B')
temp$gene <- rownames(temp)
temp$significant <- 'no'
temp$significant[(abs(temp$avg_log2FC) >= 0.25) & (temp$p_val_adj <= 0.05)] <- 'yes'
library("org.Hs.eg.db")
hs <- org.Hs.eg.db
geneID <- AnnotationDbi::select(hs, 
                                keys = temp$gene,
                                columns = c("ENTREZID", "SYMBOL"),
                                keytype = "SYMBOL")
# There is a repeat of gene Entrez, pick the later one
# if using the ENTREZID does not work because there is NA value in this column.
p <- match(temp$gene,geneID$SYMBOL)
sum((p[2:length(p)] - p[1:length(p) - 1]) != 1)
# which((p[2:length(p)] - p[1:length(p) - 1]) != 1)
# View(geneID)
# # # TEC & MEMO1
# geneID <- geneID[((rownames(geneID) != 25385) & (rownames(geneID) != 22838)),]
temp$geneID <- geneID$ENTREZID

temp$labels <- rownames(temp)
temp$labels[temp$significant == 'no'] <- ''
temp$pct_diff <- temp$pct.1 - temp$pct.2
temp <- temp %>% arrange(desc(significant),desc(avg_log2FC))
write.xlsx(temp,'tonsil_LAIV_120a_s_Bcell_GZMHpos_vs_neg_DE.xlsx')

graphics.off()
ggplot(temp, aes(x=avg_log2FC, y=-log10(p_val_adj),color = significant,legend = significant)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_colour_manual(values = c("yes" = "red", "no" = "darkgrey")) +
  ggtitle('DE GZMH+ vs GZMH- Bcells') + #PP1') +
  xlab("log2 Fold Change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 20)) +
  theme_bw() +
  geom_text_repel(aes(x=avg_log2FC, y=-log10(p_val_adj),label = labels), color = "black", size = 3) +
  geom_vline(xintercept=c(-0.25,0.25), linetype="dotted") +
  geom_hline(yintercept=c(-log10(0.05)), linetype="dotted")#+
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_GZMHpos_vs_neg_DE.pdf',width = 5.5, height = 5)

# nonmemory_nonPB.markers_subcluster <- nonmemory_nonPB.markers_subcluster %>% filter(p_val_adj <= 0.05)
# p <- match(nonmemory_nonPB.markers_subcluster$gene, all_marker_table$marker)
# temp <- all_marker_table[p,]
# nonmemory_nonPB.markers_subcluster$ENTREZID <- temp$ENTREZID
# write.xlsx(nonmemory_nonPB.markers_subcluster, "tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_subcluster_gene.xlsx")
# for (cluster_name in sort(unique(nonmemory_nonPB.subset$subclusters))[31:length(unique(nonmemory_nonPB.subset$subclusters))]){
#   temp <- nonmemory_nonPB.markers_subcluster[nonmemory_nonPB.markers_subcluster$cluster == cluster_name,]
#   temp <- temp %>% arrange(desc(avg_log2FC))
#   write.xlsx(temp, "tonsil_LAIV_120a_s_Bcell_nonmemory_nonPB_log1p_cca_scale_CycleRegressOut_subcluster_gene.xlsx",sheetName = paste('cluster',cluster_name,sep = ''),append = T)
# }
# nonmemory_nonPB.markers_subcluster.top10 <- nonmemory_nonPB.markers_subcluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# write.csv(nonmemory_nonPB.markers_subcluster.top10, "tonsil_LAIV_120a_s_Bcell_log1p_cca_scale_CycleRegressOut_subcluster_gene_top10.csv")

