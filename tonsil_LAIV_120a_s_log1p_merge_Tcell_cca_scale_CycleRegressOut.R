# install.packages('remotes')
# remotes::install_version("Seurat", version = "3.2.3")
# library(Seurat)

library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(ggrepel)
setwd("/Volumes/GoogleDrive/My\ Drive/Stanford/RNA-seq/data_analysis/tonsil_LAIV_Rhapsody_120a-s/")

##### data integration ###########################################################
# read from Nov2022 T-cell data
DefaultAssay(Tcell.subset) <- 'raw'
all.markers <- rownames(Tcell.subset)
all.genes <- all.markers[!grepl('ab-',all.markers)]
all.Ab <- all.markers[grepl('ab-',all.markers)]
Tcell.subset[["RNA"]] <- CreateAssayObject(counts = as.matrix(Tcell.subset@assays[["raw"]]@counts)[all.genes,])
Tcell.subset[["ADT"]] <- CreateAssayObject(counts = as.matrix(Tcell.subset@assays[["raw"]]@counts)[all.Ab,])
DefaultAssay(Tcell.subset) <- 'RNA'

VlnPlot(Tcell.subset, feature = "nCount_raw",group.by = 'batch',log = F,pt.size = 0.1) + geom_hline(yintercept = c(8*10^4))
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_raw_nCount_batch.pdf',width = 8, height = 6)

VlnPlot(Tcell.subset, feature = "nFeature_raw",group.by = 'batch',pt.size = 0.1) + geom_hline(yintercept = c(60,240)) # mean ~ 120
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_raw_nFeature_batch.pdf',width = 5, height = 4)

FeatureScatter(Tcell.subset, feature1 = "nCount_raw", feature2 = "nFeature_raw") + geom_hline(yintercept = c(60,240))# + geom_vline(xintercept = 1e5)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_raw_nCount_nFeature.pdf',width = 5, height = 4)

VlnPlot(Tcell.subset, feature = "nFeature_RNA",group.by = 'batch',pt.size = 0.1) + geom_hline(yintercept = c(30,200)) # mean ~ 120
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_RNA_nFeature_batch.pdf',width = 5, height = 4)

FeatureScatter(Tcell.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = c(30,200))# + geom_vline(xintercept = 1e5)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_RNA_nCount_nFeature.pdf',width = 5, height = 4)

Tcell.subset <- subset(Tcell.subset,nFeature_RNA > 30 & nFeature_RNA < 200 & nFeature_raw > 60 & nFeature_raw < 240 & nCount_raw < 8*10^4)

DefaultAssay(Tcell.subset) <- "integrated"
s.genes <- cc.genes$s.genes
s.genes <- s.genes[s.genes %in% all.markers]
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- g2m.genes[g2m.genes %in% all.markers]
Tcell.subset <- CellCycleScoring(Tcell.subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Tcell.subset <- ScaleData(Tcell.subset,vars.to.regress = c("S.Score", "G2M.Score"),features = all.markers,)#,scale.max = 6) # 3.5: ;3: 1.3% of data 

scaled_integrated_assay <- CreateAssayObject(counts = Tcell.subset@assays[["integrated"]]@scale.data)
Tcell.subset[["integrated_scale"]] <- scaled_integrated_assay
Tcell.subset <- CellCycleScoring(Tcell.subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Tcell.subset <- ScaleData(Tcell.subset,vars.to.regress = c("S.Score", "G2M.Score"),features = all.markers)
Tcell.subset@assays[["integrated_scale"]]@scale.data <- as.matrix(Tcell.subset@assays[["integrated_scale"]]@data)
DefaultAssay(Tcell.subset) <- "integrated_scale"
rm(scaled_integrated_assay)

Tcell.subset <- RunPCA(Tcell.subset, features = all.markers, npcs = 50)
# ElbowPlot(object = Tcell.subset,ndims = 50) + theme(axis.text = element_text(size = 20))
# dev.print(pdf, 'tonsil_LAIV_120a_scleaned_log1p_cca_scale_CycleRegressOut_PC_elbow.pdf',width = 7, height = 4.34)
n_pca_selected <- 30
Tcell.subset <- RunUMAP(Tcell.subset, reduction = "pca", dims = 1:n_pca_selected)
DimPlot(Tcell.subset, label = TRUE, reduction = "umap",group.by = 'donor_ID')
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_donorID.pdf',width = 5, height = 4.34)
DimPlot(Tcell.subset, reduction = "umap", split.by = "donor_ID", ncol = 3)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_donorID_split.pdf',width = 18, height = 12)

DimPlot(Tcell.subset, label = TRUE, reduction = "umap",group.by = 'batch')
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_batch.pdf',width = 5, height = 4.34)
DimPlot(Tcell.subset, reduction = "umap", split.by = "batch", ncol = 2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_batch_split.pdf',width = 10, height = 8.8)

DimPlot(Tcell.subset, label = TRUE, reduction = "umap",group.by = 'days')
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_days.pdf',width = 5, height = 4.34)
DimPlot(Tcell.subset, reduction = "umap", split.by = "day_group", ncol = 2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_day_group.pdf',width = 10, height = 8.8)
DimPlot(Tcell.subset, reduction = "umap", group.by = "age_group")
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_age_group.pdf',width = 5, height = 4.34)
DimPlot(Tcell.subset, reduction = "umap", split.by = "age_group", ncol = 2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_age_group_split.pdf',width = 10, height = 8.8)

Tcell.subset <- FindNeighbors(Tcell.subset, reduction = "pca", dims = 1:n_pca_selected)
Tcell.subset <- FindClusters(Tcell.subset, resolution = 0.5)
Tcell.subset$clusters <- Tcell.subset$seurat_clusters
# Tcell.subset$seurat_clusters_clusters <- Tcell.subset$seurat_clusters
DimPlot(Tcell.subset, label = TRUE, reduction = "umap",repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_log1p_cca_scale_CycleRegressOut_umap_clusters.pdf',width = 6, height = 4)

Tcell.markers_cluster <- FindAllMarkers(Tcell.subset)
Tcell.markers_cluster <- Tcell.markers_cluster %>% filter(p_val_adj <= 0.05)
Tcell.markers_cluster.top10 <- Tcell.markers_cluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Tcell.markers_cluster.top10, "tonsil_LAIV_120a_s_Tcell_cluster_gene_top10.csv")

write.xlsx('', "tonsil_LAIV_120a_s_Tcell_cluster_gene.xlsx")
for (cluster_name in c(8,10,13,14)){
  temp <- Tcell.markers_cluster[Tcell.markers_cluster$cluster == cluster_name,]
  temp <- temp %>% arrange(desc(avg_log2FC))
  write.xlsx(temp, "tonsil_LAIV_120a_s_Tcell_cluster_gene.xlsx",sheetName = paste('cluster',cluster_name,sep = ''),append = T)
}

temp <- FindMarkers(Tcell.subset,ident.1 = 13,ident.2 = 14)
Tcell.markers_subcluster <- FindAllMarkers(Tcell.subset)
Tcell.markers_subcluster <- Tcell.markers_subcluster %>% filter(p_val_adj <= 0.05)
Tcell.markers_subcluster.top10 <- Tcell.markers_subcluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Tcell.markers_subcluster, "tonsil_LAIV_120a_s_Tcell_subcluster_gene.csv")
write.csv(Tcell.markers_subcluster.top10, "tonsil_LAIV_120a_s_Tcell_subcluster_gene_top10.csv")

FeaturePlot(Tcell.subset, features = c("CD8A", "CD8B"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_CD8A_CD8B.pdf',width = 12, height = 4.34)
FeaturePlot(Tcell.subset, features = c("ab-TCR-gamma-delta", "ab-TCR-alpha-beta"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_ab-TCRgd_TCRab.pdf',width = 12, height = 4.34)
FeaturePlot(Tcell.subset, features = c("ab-CD8", "ab-CD4"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_ab-CD8_CD4.pdf',width = 12, height = 4.34)
FeaturePlot(Tcell.subset, features = c("CD4", "TRDC"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_CD4_TRDC.pdf',width = 12, height = 4.34)

FeaturePlot(Tcell.subset, features = c("CD3", "TRDC"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_CD4_TRDC.pdf',width = 12, height = 4.34)

VlnPlot(Tcell.subset, features = c("ab-CD4"),pt.size = 0.2,group.by = 'seurat_clusters',sort = 'increasing',log = F) + NoLegend()# + geom_hline(yintercept = -0)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_violin_abCD4.pdf',width = 12, height = 4.34)

Tcell.subset$meta_cluster <- 'CD4'
Tcell.subset$meta_cluster[Tcell.subset$seurat_clusters %in% c(28,29,18,14,24,25,31,34)] <- 'CD8 and gdT' # half BCL6+
Tcell.subset$meta_cluster[Tcell.subset$seurat_clusters %in% c(30)] <- 'B-cell' # half BCL6+
Tcell.subset$meta_cluster[Tcell.subset$seurat_clusters %in% c(33)] <- 'noise' # half BCL6+

DimPlot(Tcell.subset, label = TRUE, reduction = "umap",group.by = 'meta_cluster')

Tcell.subset <- subset(Tcell.subset,meta_cluster != 'B-cell')
Tcell.subset <- subset(Tcell.subset,meta_cluster != 'noise')

DefaultAssay(Tcell.subset) <- 'RNA'
Tcell.subset <- NormalizeData(Tcell.subset)
DefaultAssay(Tcell.subset) <- 'ADT'
Tcell.subset <- NormalizeData(Tcell.subset, normalization.method = 'CLR', margin = 2)
Tcell.subset[['normalized']] <- CreateAssayObject(rbind(as.matrix(Tcell.subset@assays[["ADT"]]@data),as.matrix(Tcell.subset@assays[["RNA"]]@data)))
DefaultAssay(Tcell.subset) <- "integrated_scale"


temp <- subset(Tcell_LAIV.subset, days %in% c('day0','day04','day06','day07','day10','day12','day14'))
temp$day_group <- 'day0'
temp$day_group[temp$days %in% c('day04')] <- 'day04'
temp$day_group[temp$days %in% c('day06','day07')] <- 'day06-07'
temp$day_group[temp$days %in% c('day10')] <- 'day10'
temp$day_group[temp$days %in% c('day12','day14')] <- 'day12-14'

gene_name <- 'SELL'
VlnPlot(temp, features = c(gene_name),pt.size = 0.1, ncol = 1,split.by = 'age_group',log = F,assay = 'normalized',group.by = 'day_group')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_violin_normalized_age_group_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

############ CD8 and gdT reclustering #####################################################
CD8_gdT.subset <- subset(Tcell.subset,meta_cluster == 'CD8 and gdT')
CD8_gdT.subset <- RunPCA(CD8_gdT.subset, features = all.markers, npcs = 50)
n_pca_selected <- 30
CD8_gdT.subset <- RunUMAP(CD8_gdT.subset, reduction = "pca", dims = 1:n_pca_selected)
CD8_gdT.subset <- FindNeighbors(CD8_gdT.subset, reduction = "pca", dims = 1:n_pca_selected)
CD8_gdT.subset <- FindClusters(CD8_gdT.subset, resolution = 2)
DimPlot(CD8_gdT.subset, label = TRUE, reduction = "umap",repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_subcluster.pdf',width = 7, height = 4.34)
DimPlot(CD8_gdT.subset, label = TRUE, reduction = "umap",group.by = 'batch')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_batch.pdf',width = 7, height = 4.34)

# lineage_marker <- c('ab-CD3','ab-CD20','ab-CD19','ab-CD4','ab-CD8','ab-CD56','ab-TCR-gamma-delta','ab-TCR-alpha-beta')
# lineage_gene <- c('CD3D','CD3E','CD3G','MS4A1','CD4','CD8A','CD8B','FCGR3A','TRAC','TRBC2','TRDC')
# Tcell_marker <- c('ab-CD27','ab-CD38','ab-HLA-DR','ab-CD45RA','ab-CD45RO','ab-CCR7','ab-CD62L','ab-CD28','ab-Tim3','ab-KIR-NKAT2','ab-CD158e1')
# Tcell_gene <- c('CD27','CD38','HLA-DRA','CD69','CCR7','CD44','SELL','CD28','MKI67','RORC','TBX21','EOMES','GNLY','GZMK','GZMA','GZMB','NKG7','KLRC3','KLRC4','IFNG')
# features <- list('lineage_marker' = lineage_marker, 'lineage_gene'= lineage_gene,"Tcell marker" = Tcell_marker, "Tcell gene" = Tcell_gene)
# CD8_gdT.subset$clusters <- as.character(CD8_gdT.subset$seurat_clusters)
# # CD8_gdT.subset$clusters <- factor(CD8_gdT.subset$clusters,levels = c(gdT_index,CD8_index,CD4_index))
# DotPlot(object = CD8_gdT.subset, group.by = 'clusters',dot.scale = 10,features=features) + RotatedAxis() + ggtitle('Tcell clusters') +
#   theme(axis.text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold"))
# dev.print(pdf, paste('tonsil_LAIV_120a_s_CD8_gdT_subcluster_dotplot.pdf',sep = ''),width = 2200, height = 800)

FeaturePlot(CD8_gdT.subset, features = c("CD8A", "CD8B"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_CD8A_CD8B.pdf',width = 12, height = 4.34)
FeaturePlot(CD8_gdT.subset, features = c("ab-TCR-gamma-delta", "ab-TCR-alpha-beta"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_ab-TCRgd_TCRab.pdf',width = 12, height = 4.34)
FeaturePlot(CD8_gdT.subset, features = c("ab-CD8", "ab-CD4"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_ab-CD8_CD4.pdf',width = 12, height = 4.34)
FeaturePlot(CD8_gdT.subset, features = c("TRAC", "TRDC"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_TRAC_TRDC.pdf',width = 12, height = 4.34)
FeaturePlot(CD8_gdT.subset, features = c("TRBC1", "TRBC2"),pt.size = 0.2, ncol = 2, sort.cell = TRUE,min.cutoff = 0,max.cutoff = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_TRBC1_TRBC2.pdf',width = 12, height = 4.34)

# CD8_gdT.markers_subcluster <- FindAllMarkers(CD8_gdT.subset)
# CD8_gdT.markers_subcluster <- CD8_gdT.markers_subcluster %>% filter(p_val_adj <= 0.05)
# CD8_gdT.markers_subcluster.top10 <- CD8_gdT.markers_subcluster %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# write.csv(CD8_gdT.markers_subcluster, "tonsil_LAIV_120a_s_CD8_gdT_subcluster_gene.csv")
# write.csv(CD8_gdT.markers_subcluster.top10, "tonsil_LAIV_120a_s_CD8_gdT_subcluster_gene_top10.csv")
# DoHeatmap(CD8_gdT.subset, features = c(CD8_gdT.markers_subcluster.top10$gene)) +
#   theme(text = element_text(size=20))
# dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_gene_heatmap_subclusters.pdf',width = 18, height = 25)
# 
# gene_name <- 'TRDC'#'ab-TCR-gamma-delta'
# c((CD8_gdT.markers_subcluster %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))['cluster'])
# 
# FeaturePlot(Bcell.subset, features = c(gene_name),pt.size = 0.2, ncol = 1, order = TRUE,min.cutoff = 0)#, max.cutoff = 4)
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_umap_',gene_name,'.pdf',sep = ''),width = 5, height = 4.3)
# 
# CD8_gdT.subset$meta_cluster <- 'CD8'
# CD8_gdT.subset$meta_cluster[CD8_gdT.subset$seurat_clusters %in% c(5,0,25,6,20,8,21,14)] <- 'gdT'
# DimPlot(CD8_gdT.subset, label = TRUE, reduction = "umap",group.by = 'meta_cluster')
# dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_meta_cluster.pdf',width = 7, height = 4.34)


CD8_gdT.subset$CD8_gdT_gate <- 'CD8'
CD8_gdT.subset$CD8_gdT_gate[as.matrix(CD8_gdT.subset@assays[["log1p"]]['TRDC']) > 0] <- 'gdT'
DimPlot(CD8_gdT.subset, label = TRUE, reduction = "umap",group.by = 'CD8_gdT_gate')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_CD8_gdT_gate.pdf',width = 7, height = 4.34)

# FeatureScatter(CD8_gdT.subset, feature1 = "TRDC", feature2 = "CD8A",group.by = 'CD8_gdT_gate') + ggtitle('CD8 and gdT cells')
# dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_scatter_TRDC_CD8A.pdf',width = 7, height = 4.34)
# 
# CD8_gdT.subset$TCRgd_TCRab_gate <- 'CD8'
# CD8_gdT.subset$TCRgd_TCRab_gate[as.matrix(CD8_gdT.subset@assays[["log1p"]]['ab-TCR-gamma-delta']) > as.matrix(CD8_gdT.subset@assays[["log1p"]]['ab-TCR-alpha-beta'])] <- 'gdT'
# DimPlot(CD8_gdT.subset, label = TRUE, reduction = "umap",group.by = 'TCRgd_TCRab_gate')
# dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_gdT_umap_TCRgd_TCRab_gate.pdf',width = 7, height = 4.34)
# 
# FeatureScatter(CD8_gdT.subset, feature1 = "ab-TCR-gamma-delta", feature2 = "ab-TCR-alpha-beta",group.by = 'CD8_gdT_gate')
# FeatureScatter(CD8_gdT.subset, feature1 = "TRAC", feature2 = "TRDC",group.by = 'CD8_gdT_gate')
# FeatureScatter(CD8_gdT.subset, feature1 = "TRBC2", feature2 = "TRDC",group.by = 'CD8_gdT_gate')

CD8_gdT.subset$meta_cluster <- CD8_gdT.subset$CD8_gdT_gate
Tcell.subset$meta_cluster[Tcell.subset$meta_cluster != 'CD4'] <- CD8_gdT.subset$meta_cluster
# Tcell.subset$TCR_meta_cluster <- 'CD4'
# Tcell.subset$TCR_meta_cluster[Tcell.subset$meta_cluster != 'CD4'] <- CD8_gdT.subset$CD8_gdT_gate
DimPlot(Tcell.subset, label = TRUE, reduction = "umap",group.by = 'meta_cluster')
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_meta_clusters.pdf',width = 7, height = 4.34)

Tcell_LAIV.subset <- subset(Tcell.subset, stimulation != 'LAIV-')
Tcell_LAIV.subset <- subset(Tcell_LAIV.subset,days != 'day02')
Tcell_LAIV.subset <- subset(Tcell_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))

saveRDS(CD8_gdT.subset,'tonsil_LAIV_120a_s_CD8_gdT_cycleRegress.rds')
saveRDS(Tcell.subset,'tonsil_LAIV_120a_s_Tcell.rds')

############ CD4 reclustering #####################################################
CD4.subset <- subset(Tcell.subset,meta_cluster == 'CD4')
CD4.subset <- RunPCA(CD4.subset, features = all.markers, npcs = 50)
CD4.subset <- RunUMAP(CD4.subset, reduction = "pca", dims = 1:n_pca_selected)
CD4.subset <- FindNeighbors(CD4.subset, reduction = "pca", dims = 1:n_pca_selected)
CD4.subset <- FindClusters(CD4.subset, resolution = 0.5)
CD4.subset$cluster <- CD4.subset$seurat_clusters
CD4.subset <- subset(CD4.subset, cluster != 12)

# Tcell.subset$CD4_subclusters  <- Tcell.subset$meta_cluster
# Tcell.subset$CD4_subclusters[Tcell.subset$meta_cluster == 'CD4'] <- CD4.subset$seurat_clusters
# Tcell.subset <- subset(Tcell.subset, CD4_subclusters != '26') #if factor is 25, then changed to character is 26
# 
# saveRDS(CD4.subset,'tonsil_LAIV_120a_s_CD4_unfiltered.rds')
# CD4.subset <- readRDS('tonsil_LAIV_120a_s_CD4_unfiltered.rds')
# 
# CD4.subset <- subset(CD4.subset, seurat_clusters != '26')
# CD4.subset <- readRDS('tonsil_LAIV_120a_s_CD4.rds')
# 
# # saveRDS(,'tonsil_LAIV_120a_s_CD4_unfiltered.rds')

DimPlot(CD4.subset, label = TRUE, reduction = "umap",repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD4_umap_cluster.pdf',width = 7, height = 4.34)
DimPlot(CD4.subset, label = TRUE, reduction = "umap",group.by = 'project_ID')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD4_umap_projectID.pdf',width = 7, height = 4.34)
DimPlot(CD4.subset, label = TRUE, reduction = "umap",group.by = 'donor_ID')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD4_umap_donorID.pdf',width = 7, height = 4.34)
DimPlot(CD4.subset, label = TRUE, reduction = "umap",group.by = 'days')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD4_umap_days.pdf',width = 7, height = 4.34)

# CD4.cluster2 <- FindMarkers(CD4.subset, ident.1 = '14')
CD4.cluster.markers <- FindAllMarkers(CD4.subset)
CD4.cluster.markers <- CD4.cluster.markers %>% filter(p_val_adj <= 0.05)
p <- match(CD4.cluster.markers$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
CD4.cluster.markers$ENTREZID <- temp$ENTREZID
CD4.cluster.markers.top10 <- CD4.cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# DoHeatmap(CD4.subset, features = c(CD4.cluster.markers.top10$gene)) +
#   theme(text = element_text(size=20))
# dev.print(pdf, 'tonsil_LAIV_120a_s_CD4_cluster_gene_heatmap.pdf',width = 12, height = 25)
write.xlsx('', "tonsil_LAIV_120a_s_CD4_cluster_gene.xlsx")
for (cluster_name in sort(unique(CD4.cluster.markers$cluster))){
  temp <- CD4.cluster.markers[CD4.cluster.markers$cluster == cluster_name,]
  temp <- temp %>% arrange(desc(avg_log2FC))
  write.xlsx(temp, "tonsil_LAIV_120a_s_CD4_cluster_gene.xlsx",sheetName = paste('cluster',cluster_name,sep = ''),append = T)
}
write.csv(CD4.cluster.markers.top10, "tonsil_LAIV_120a_s_CD4_cluster_gene_top10.csv")

table1 <- 'tonsil_LAIV_120a_s_Tfh_cluster_nonday0_subcluster_labeled_gene.xlsx'
table1_DE_data <- read_excel(table1, sheet = "Sheet1")
# table1_DE_data <- table1_DE_data[,!colnames(table1_DE_data) == '...1']
p <- match(table1_DE_data$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
table1_DE_data$ENTREZID <- temp$ENTREZID
write.xlsx(table1_DE_data,sheetName = 'Sheet2',table1,append = T)


##### CD4 annotation #####################
grep('IL15R', all.markers, value=TRUE)

DefaultAssay(CD4.subset) <- 'integrated_scale'
gene_name <- 'IL2RA'
View(CD4.cluster.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))
c((CD4.cluster.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))['cluster'])
FeaturePlot(CD4.subset, features = c(gene_name),pt.size = 0.2, ncol = 1, order = TRUE, min.cutoff = 0)#, max.cutoff = 15)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_umap_',gene_name,'.pdf',sep = ''),width = 5, height = 4.3)

VlnPlot(CD4.subset, features = c(gene_name),pt.size = 0.1, ncol = 1,group.by = 'donor_ID',log = F,assay = 'log1p') + geom_hline(yintercept = 6.5)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_violin_log1p_donorID_',gene_name,'.pdf',sep = ''),width = 15, height = 6)
VlnPlot(CD4.subset, features = c(gene_name),pt.size = 0.1, ncol = 1,group.by = 'donor_ID',log = F,assay = 'integrated_scale') + geom_hline(yintercept = 0)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_violin_integrated_scale_donorID_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

temp <- subset(CD4_LAIV.subset, days %in% c('day0','day04','day06','day07','day10','day12','day14'))
temp$day_group <- 'day0'
temp$day_group[temp$days %in% c('day04')] <- 'day04'
temp$day_group[temp$days %in% c('day06','day07')] <- 'day06-07'
temp$day_group[temp$days %in% c('day10')] <- 'day10'
temp$day_group[temp$days %in% c('day12','day14')] <- 'day12-14'

temp_Th1_17 <- subset(temp, major_cluster == 'Th1*')
temp_Tfh <- subset(temp, major_cluster == 'Tfh')
gene_name <- 'LTB'
VlnPlot(temp_Th1_17, features = c(gene_name),pt.size = 0.1, ncol = 1,split.by = 'age_group',log = F,assay = 'normalized',group.by = 'day_group')
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_Th1_Th17_violin_normalized_age_group_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

VlnPlot(CD4.subset, features = c(gene_name),pt.size = 0.1, ncol = 1,sort = 'increasing',log = F,assay = 'integrated_scale') + geom_hline(yintercept = -0.2)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_violin_integrated_scale_',gene_name,'.pdf',sep = ''),width = 15, height = 10)
VlnPlot(CD4.subset, features = c(gene_name),pt.size = 0.1, ncol = 1,sort = 'increasing',log = F,assay = 'log1p') + geom_hline(yintercept = 2.7)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_violin_log1p_',gene_name,'.pdf',sep = ''),width = 15, height = 10)
VlnPlot(CD4_responder.subset, features = gene_name,pt.size = 0.1,split.by = 'donor_ID',group.by = 'days_treatment',assay = 'log1p',slot = 'data')
# VlnPlot(CD4.subset, features = gene_name,pt.size = 0.2,split.by = 'treatment',group.by = 'cluster_manual_label',assay = 'RNA',slot = 'data',sort = 'increasing')
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_responder_violin_',gene_name,'_log1p.pdf',sep = ''),width = 15, height = 7)
VlnPlot(CD4_day0.subset, features = gene_name,pt.size = 0.1,split.by = 'donor_ID',group.by = 'LAIV_response',assay = 'log1p',slot = 'data')
# VlnPlot(CD4.subset, features = gene_name,pt.size = 0.2,split.by = 'treatment',group.by = 'cluster_manual_label',assay = 'RNA',slot = 'data',sort = 'increasing')
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_day0_violin_',gene_name,'_log1p.pdf',sep = ''),width = 10, height = 7)
# ab-CD45RA: 10 11 12
# ab-CD45RO: everything but CD45RA positive 
# ab-CCR7: 12 10 5  11 9 

# RORC: 0 3 (8)
# TBX21: 0 3
# IFNG: 0 3
# TNF: 8
# IL2: 8
# IL17F: 3
# IL17A: none
# ab-CD161: 0 3 12
# FOXP3: 9
# MZB1: 3 9
# ab-CD62L: 
# ab-CD27: 
# ab-CD183: 12 9  0  5 
# CXCR3: 0 9
# CCR6: 0
# CD40LG: 5 8

# BCL6: 8 1 5 4
# CXCR5: 5 4 1
# ab-CXCR5: 4 1 6 5 2

# ab-CD279: 5 8 1 4 3
# PDCD1: 
# ab-CD278 (ICOS): 3  9  12 5  8 
# ICOS: 
# IL21: 5 8
# IL21R: 
# STAT4: 
# AICDA: 1
# CXCL13: 
# CXCR4: 
# ab-CD184: 


# SELPLG: 
# GATA3: 


# cytolytic markers (better to use the gating!)
# GZMA: 0  13 4  9 
# PRF1: 0  13 23 15
# GNLY: 0  13 23 4 
# NKG7: 0  13

# MKI67:12
# GAPDH: 12 24 29 22 17 13 20 4  32
# ab-CD154: 2  26 19 30 31
# ab-CD69: 2  26 4  20 13 5  7  23 17 30
# CD69: 28 30 8  16
# CD40LG: 4  32 12 2  24
# ab-CD56: 

# FOXP3: 15 24
# IL10: 24 5  15 2 
# ab-CD25: 15 23 2  26 0 
# IL2RA: 15 23 2  18

# TGFB1: 0  13
# GATA3: 15 18 0  17 23 13



# ab-CD38: 2  0  23 17 9  26 13 18 21
# CD38: 18 23 2  17 28 0 
# ab-HLA-DR: 23 2  0  15 13 26 9  24 27
# HLA-DRA: 23

# EGR3: 16 11
# EOMES: 19 29 25 23 26
# B3GAT1(CD57): 35 32 31 34
# PDCD1: 20 24 11 7  2 (align with TFH)
# ab-CD279(PD1): 35 1  31 27 3  2  33 14 10
# ICOS: 12 10 17 27 1  2  35 4 
# ab-Tim3: 9  29 15 6  14 24 7  18 23 16
# CTLA4: 15 16 8  7  6  18 11
# ARG1 (checkpoint): 22 15 29 11 4  0 not much
# TNFRSF4 (OX40): 12 0  10 11 17 3  14 2 
# TNFRSF8 (CD30): 7  0  15
# TNFSF8: 11 0  16 7 
# SLAMF1: 13 12 36 8  4  17 11 22 0 

# MX1: 17 21 0  18 16 2  23 7 
# IFIT3: 
# STAT1: 

# VEGFA: 7  11 19 5  6  18 2 
# ALDOC: 
# TIGIT 
# CD200: 
# EGR1: 
# EGR3: 
# ADGRG1 (GPR56):

# TNFSF13: 0
# TNFSF13B: 0  13 18 4  17 9 


CD4.subset$cluster_labeled <- 'CD4'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(0)] <- 'cytotoxic Th1/17'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(1)] <- 'AID+Tfh'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(2)] <- 'VEGFA+Tfh-like'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(3)] <- 'IL17+Th1/17'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(4)] <- 'Tfh'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(5)] <- 'early act. CD4' 
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(6)] <- 'Tfh-like'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(7)] <- 'Tcm'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(8)] <- 'IL2+Tfh'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(9)] <- 'Treg'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(10)] <- 'naive CD4'
CD4.subset$cluster_labeled[CD4.subset$cluster %in% c(11)] <- 'early act. naive CD4'
ordered_cluster_annotated_list <- c('Treg','cytotoxic Th1/17','IL17+Th1/17','IL2+Tfh','AID+Tfh','Tfh','VEGFA+Tfh-like',
                                    'Tfh-like','Tcm','naive CD4','early act. naive CD4','early act. CD4')
CD4.subset$cluster_labeled <- factor(CD4.subset$cluster_labeled,
                                         levels = ordered_cluster_annotated_list)

CD4_LAIV.subset <- subset(CD4.subset, stimulation != 'LAIV-')
CD4_LAIV.subset <- subset(CD4_LAIV.subset,days != 'day02')
CD4_LAIV.subset <- subset(CD4_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
CD4_LAIV.subset <- subset(CD4_LAIV.subset,donor_ID != "02yrs M IMD085")


DimPlot(CD4_LAIV.subset, label = T, reduction = "umap",group.by = 'cluster_labeled', repel = T)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_LAIV_umap_cluster_labeled.pdf',sep = ''),width = 7, height = 5)

CD4_LAIV.subset$cluster_labeled <- factor(CD4_LAIV.subset$cluster_labeled,
                                              levels = rev(ordered_cluster_annotated_list))

# memory_marker <- c('ab-CD45RA','ab-CD45RO','ab-CCR7','ab-CD62L')
# Tfh_marker <- c('ab-CXCR5','BCL6','ab-CD279','ab-CD278')
# Th17_marker <- c('ab-CD161','RORC')
# Treg_marker <- c('FOXP3','ab-CD25')
# features <- list("memory" = memory_marker, "Tfh" = Tfh_marker, 'Th17' = Th17_marker,"Treg" = Treg_marker)
DefaultAssay(CD4_LAIV.subset) <- 'integrated_scale'
features <- c('FOXP3','ab-CD25','GNLY','GZMA','RORC','ab-CD161','TBX21','IFNG','IL17F','IL2','TNF',
              'BCL6','ab-CXCR5','AICDA','VEGFA','ab-CCR7','ab-CD62L','ab-CD45RA','ab-CD45RO')

DotPlot(CD4_LAIV.subset, group.by = 'cluster_labeled',dot.scale = 5,features=features) + RotatedAxis() +
  theme(axis.text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) +
  ylab('CD4 subsets')
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_LAIV_cluster_labeled_dotplot.pdf',sep = ''),width = 8, height = 3.5)

Tcell.subset$cluster_labeled <- Tcell.subset$meta_cluster
p <- match(colnames(Tcell.subset)[Tcell.subset$cluster_labeled == 'CD4'],colnames(CD4.subset))
Tcell.subset$cluster_labeled[Tcell.subset$cluster_labeled == 'CD4'] <- as.character(CD4.subset$cluster_labeled[p])

# Tcell.subset$cluster_labeled[Tcell.subset$subclusters %in% c(34,28)] <- 'naive CD8'
# Tcell.subset$cluster_labeled[(Tcell.subset$subclusters %in% c(31,18,24)) & (Tcell.subset$meta_cluster == 'CD8')] <- 'RORC-CD8'
# Tcell.subset$cluster_labeled[(Tcell.subset$subclusters %in% c(31,18,24)) & (Tcell.subset$meta_cluster == 'gdT')] <- 'RORC-gdT'
# Tcell.subset$cluster_labeled[(Tcell.subset$subclusters %in% c(14,25,29)) & (Tcell.subset$meta_cluster == 'CD8')] <- 'RORC+CD8'
# Tcell.subset$cluster_labeled[(Tcell.subset$subclusters %in% c(14,25,29)) & (Tcell.subset$meta_cluster == 'gdT')] <- 'RORC+gdT'

# Tcell.subset$cluster_labeled[is.na(Tcell.subset$cluster_labeled)] <-'noise'
# Tcell.subset <- subset(Tcell.subset,cluster_labeled != 'noise')
DimPlot(Tcell.subset, label = T, reduction = "umap",group.by = 'cluster_labeled', repel = T)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_LAIV_umap_cluster_labeled.pdf',sep = ''),width = 7, height = 5)

ordered_cluster_annotated_list <- c('Treg','cytotoxic Th1/17','IL17+Th1/17','IL2+Tfh','AID+Tfh','Tfh','VEGFA+Tfh-like',
                                    'Tfh-like','Tcm','naive CD4','early act. naive CD4','early act. CD4','CD8','gdT')
Tcell.subset$cluster_labeled <- factor(Tcell.subset$cluster_labeled)

Tcell_LAIV.subset <- subset(Tcell.subset, stimulation != 'LAIV-')
Tcell_LAIV.subset <- subset(Tcell_LAIV.subset,days != 'day02')
Tcell_LAIV.subset <- subset(Tcell_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
Tcell_LAIV.subset <- subset(Tcell_LAIV.subset,donor_ID != "02yrs M IMD085")

p <- DimPlot(Tcell.subset, label = T, reduction = "umap",group.by = 'cluster_labeled', repel = T)
pbuild <- ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # Pull the data used for the plot
pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
cols_group <- unique(pdata$colour) # Get a vector of unique colors
names(cols_group) <- levels(pbuild[["data"]][[2]][["label"]])

Tcell_LAIV.subset$cluster_labeled <- factor(Tcell_LAIV.subset$cluster_labeled,
                                            levels = ordered_cluster_annotated_list)

cols_group <- cols_group[ordered_cluster_annotated_list]
DimPlot(Tcell_LAIV.subset, cols = cols_group, label = T, reduction = "umap",group.by = 'cluster_labeled', repel = T)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_umap_cluster_labeled.pdf',sep = ''),width = 7, height = 5)

Tcell_LAIV.subset$cluster_labeled <- factor(Tcell_LAIV.subset$cluster_labeled,
                                          levels = rev(ordered_cluster_annotated_list))

DefaultAssay(Tcell_LAIV.subset) <- 'integrated_scale'
features <- c('FOXP3','ab-CD25','GNLY','GZMA','RORC','ab-CD161','TBX21','IFNG','IL17F','IL2','TNF',
              'BCL6','ab-CXCR5','AICDA','VEGFA','ab-CCR7','ab-CD62L','ab-CD45RA','ab-CD45RO','CD8A','CD8B','TRDC')

DotPlot(Tcell_LAIV.subset, group.by = 'cluster_labeled',dot.scale = 5,features=features) + RotatedAxis() +
  theme(axis.text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) +
  ylab('T-cell subsets')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_cluster_labeled_dotplot.pdf',sep = ''),width = 8, height = 3.5)

KO_cells <- 1:dim(Tcell_LAIV.subset)[2]
downsampled_KO_cells <- sample(KO_cells, dim(Tcell_LAIV.subset)[2]/5)
temp <- Tcell_LAIV.subset[,downsampled_KO_cells]
# DimPlot(WT_KO_integrated_downsampled, reduction = "umap", split.by ="orig.ident", ncol=2)
DimPlot(temp, cols = cols_group, reduction = "umap", pt.size = 0.03, group.by = "cluster_labeled",split.by = "age_group", ncol = 3) + 
  theme(text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold")) 
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_LAIV_umap_age_group_split_cluster_annotated2.pdf',width = 7, height = 2.5)

CD4.subset$major_cluster <- 'CD4'
CD4.subset$major_cluster[CD4.subset$cluster %in% c(0,3)] <- 'Th1*'
CD4.subset$major_cluster[CD4.subset$cluster %in% c(1,4,8)] <- 'Tfh'
CD4.subset$major_cluster[CD4.subset$cluster %in% c(2,6)] <- 'Tfh-like'
CD4.subset$major_cluster[CD4.subset$cluster %in% c(5,11)] <- 'early activated CD4' 
CD4.subset$major_cluster[CD4.subset$cluster %in% c(7,10)] <- 'naive/Tcm'
CD4.subset$major_cluster[CD4.subset$cluster %in% c(9)] <- 'Treg'
DimPlot(CD4.subset, label = T, reduction = "umap",group.by = 'major_cluster', repel = T)
DimPlot(CD4_LAIV.subset, label = T, reduction = "umap",group.by = 'major_cluster', repel = T)

cytokine_global_table <- read.xlsx('../Global landscape of cytokines Supplementary Table S1.xlsx', sheetName = "Cytokines")
cytokine_global_list <- unique(cytokine_global_table$HGNC.symbol)
cytokine_global_list <- c(cytokine_global_list,'GZMB','PRF1','MIF','FASLG','VEGFA','TGFA')
seurat_cytokine_list <- all.markers[all.markers %in% cytokine_global_list]

temp <- CD4.cluster.markers[CD4.cluster.markers$avg_log2FC > 0,]
temp <- temp[temp$gene %in% seurat_cytokine_list,]
write.csv(temp,'tonsil_LAIV_120a_s_CD4_cluster_cytokine.csv')

temp_shrink_CD4 <- CD4.subset
temp_shrink_CD4[["raw"]] <- NULL
temp_shrink_CD4[["integrated"]] <- NULL
temp_shrink_CD4[["RNA"]] <- NULL
temp_shrink_CD4[["ADT"]] <- NULL
saveRDS(temp_shrink_CD4,'tonsil_LAIV_120a_s_CD4.rds')

temp_shrink_CD4 <- subset(temp_shrink_CD4, stimulation != 'LAIV-')
saveRDS(temp_shrink_CD4,'tonsil_LAIV_120a_s_lognorm_CD4_LAIV.rds')

temp <- subset(temp_shrink_CD4,major_cluster == 'Treg')
saveRDS(temp,'tonsil_LAIV_120a_s_CD4_Treg.rds')
temp <- subset(temp_shrink_CD4,major_cluster == 'Th1*')
saveRDS(temp,'tonsil_LAIV_120a_s_CD4_Th1*.rds')
temp <- subset(temp_shrink_CD4,major_cluster == 'Tfh')
saveRDS(temp,'tonsil_LAIV_120a_s_CD4_Tfh.rds')
temp <- subset(temp_shrink_CD4,major_cluster == 'Tfh-like')
saveRDS(temp,'tonsil_LAIV_120a_s_CD4_Tfh_like.rds')
temp <- subset(temp_shrink_CD4,major_cluster == 'naive/Tcm')
saveRDS(temp,'tonsil_LAIV_120a_s_CD4_NaiveOrTcm.rds')
rm(temp)
rm(temp_shrink_CD4)

CD8.subset <- subset(Tcell.subset,meta_cluster == 'CD8')
saveRDS(CD8.subset,'tonsil_LAIV_120a_s_CD8.rds')

CD4_LAIV.subset <- subset(CD4.subset, stimulation != 'LAIV-')
CD4_LAIV.subset <- subset(CD4_LAIV.subset,days != 'day02')
CD4_LAIV.subset <- subset(CD4_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
CD4_LAIV.subset <- subset(CD4_LAIV.subset,donor_ID != "02yrs M IMD085")

temp <- subset(CD4_LAIV.subset, days %in% c('day0','day04','day06','day07','day10','day12','day14'))
temp$day_group <- 'day0'
temp$day_group[temp$days %in% c('day04')] <- 'day04'
temp$day_group[temp$days %in% c('day06','day07')] <- 'day06-07'
temp$day_group[temp$days %in% c('day10')] <- 'day10'
temp$day_group[temp$days %in% c('day12','day14')] <- 'day12-14'
gene_name <- 'GATA3'
VlnPlot(temp, features = c(gene_name),pt.size = 0.1, ncol = 1,split.by = 'age_group',log = F,assay = 'normalized',group.by = 'day_group')
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_violin_normalized_age_group_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

temp <- subset(Tcell.subset,meta_cluster == 'gdT')
saveRDS(temp,'tonsil_LAIV_120a_s_gdT.rds')

Tcell.subset$differentiation  <- Tcell.subset$meta_cluster
Tcell.subset$differentiation[Tcell.subset$meta_cluster == 'CD4'] <- CD4.subset$differentiation2
DimPlot(Tcell.subset, label = TRUE, repel = T, reduction = "umap",group.by = 'differentiation')
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_differentiation2.pdf',width = 600, height = 4.34)

saveRDS(Tcell.subset,'tonsil_LAIV_120a_s_Tcell.rds')
saveRDS(Tcell_LAIV.subset,'tonsil_LAIV_120a_s_Tcell_LAIV.rds')

CD4_nonday0.subset <- subset(CD4.subset,days != 'day0')
temp <- FindMarkers(CD4_nonday0.subset,group.by = 'differentiation',ident.1 = 'CD4',ident.2 = NULL)
CD4_nonday0.subset <- SetIdent(CD4_nonday0.subset, value = 'differentiation')
temp <- FindAllMarkers(CD4_nonday0.subset)
write.xlsx(temp,'tonsil_LAIV_120a_s_CD4_differentiation.xlsx')
DimPlot(CD4.subset, label = FALSE, reduction = "umap",split.by = 'differentiation', ncol = 4)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD4_umap_differentiation_split.pdf',width = 12, height = 6)

Tfh.subset <- subset(CD4.subset, BCL6_pos == 'BCL6+CD4')
CXCR5posCD4.subset <- subset(CD4.subset, AbCXCR5AbPD1_BCL6_pos != 'CXCR5-BCL6-CD4')


temp <- subset(Tcell.subset,batch != 'batch4')
saveRDS(temp,'tonsil_LAIV_120a_s_Tcell_120abc.rds')


Tcell.subset$AbCD38_pos <- 'CD38-Tcell'
Tcell.subset$AbCD38_pos[(as.matrix(Tcell.subset@assays[["integrated_scale"]]['ab-CD38']) > 0)] <- 'CD38+Tcell'
Tcell.subset$AbHLADR_pos <- 'HLA-DR-Tcell'
Tcell.subset$AbHLADR_pos[(as.matrix(Tcell.subset@assays[["integrated_scale"]]['ab-HLA-DR']) > 0)] <- 'HLA-DR+Tcell'
Tcell.subset$AbCD38_AbHLADR_pos <- paste(gsub('Tcell','',Tcell.subset$AbCD38_pos), Tcell.subset$AbHLADR_pos, sep = '')

Tcell.subset$CD38_pos <- 'CD38-Tcell'
Tcell.subset$CD38_pos[(as.matrix(Tcell.subset@assays[["log1p"]]['CD38']) > 0)] <- 'CD38+Tcell'
Tcell.subset$HLADR_pos <- 'HLA-DR-Tcell'
Tcell.subset$HLADR_pos[(as.matrix(Tcell.subset@assays[["log1p"]]['HLA-DRA']) > 0)] <- 'HLA-DR+Tcell'
Tcell.subset$CD38_HLADR_pos <- paste(gsub('Tcell','',Tcell.subset$CD38_pos), Tcell.subset$HLADR_pos, sep = '')

TCR_index <- data.frame(project = Tcell.subset$project_ID, cell_index = gsub('_.*','',colnames(Tcell.subset)), sample = Tcell.subset$donor_ID, sample_condition = Tcell.subset$condition, 
                        Tcell_subset = Tcell.subset$differentiation,CD38_HLADR_pos = Tcell.subset$CD38_HLADR_pos, AbCD38_AbHLADR_pos = Tcell.subset$AbCD38_AbHLADR_pos)
write.csv(TCR_index,'tonsil_LAIV_120a_s_Tcell_project_cell_index_differentiation_CD38_HLADR.csv')
CD4.subset$cell_index <- gsub('_.*','',colnames(CD4.subset))
CD4.subset$project_cellindex <- paste(CD4.subset$project_ID,CD4.subset$cell_index)
saveRDS(CD4.subset,'tonsil_LAIV_120a_s_CD4_orginial_cellindex.rds')

for (day_group_name in unique(Tcell.subset$day_group)){
  graphics.off()
  plot <- DimPlot(subset(Tcell.subset,day_group == day_group_name), split.by = 'age_group',reduction = "umap",group.by = 'stimulation',pt.size = 0.1) + ggtitle(day_group_name) +
    scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_umap_',day_group_name,'.pdf',sep = ''),width = 10, height = 4)
}

for (day_group_name in unique(CD4.subset$day_group)){
  graphics.off()
  plot <- DimPlot(subset(CD4.subset,day_group == day_group_name), split.by = 'age_group',reduction = "umap",group.by = 'stimulation',pt.size = 0.1) + ggtitle(day_group_name) +
    scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_umap_',day_group_name,'.pdf',sep = ''),width = 10, height = 4)
}

CD4_nonday0.subset <- subset(CD4.subset, days != 'day0')
CD4_memory.subset <- subset(CD4.subset, AbCD45RA_AbCD45RO_pos == 'memory CD4')
Tfh_cluster.subset <- subset(CD4.subset, cluster_labeled == 'Tfh')
Tfh_cluster_nonday0.subset <- subset(Tfh_cluster.subset, days != 'day0')
# Tfh_cluster_nonday0.subset$temp <- 0
# Tfh_cluster_nonday0.subset$temp[Tfh_cluster_nonday0.subset$subcluster %in% c(10,11)] <- 1
Idents(Tfh_cluster_nonday0.subset) <- Tfh_cluster_nonday0.subset$subcluster
# CD4_nonday0.markers_differentiation <- FindMarkers(Tfh_cluster_nonday0.subset, group.by = "temp",ident.1 = 1,ident.2 = NULL,logfc.threshold = 0.25,min.pct = 0)
# CD4_nonday0.markers_differentiation$gene <- rownames(CD4_nonday0.markers_differentiation)
CD4_nonday0.markers_differentiation <- FindAllMarkers(Tfh_cluster_nonday0.subset)
CD4_nonday0.markers_differentiation <- CD4_nonday0.markers_differentiation %>% filter(p_val_adj <= 0.05)
p <- match(CD4_nonday0.markers_differentiation$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
CD4_nonday0.markers_differentiation$ENTREZID <- temp$ENTREZID
write.xlsx(CD4_nonday0.markers_differentiation, "tonsil_LAIV_120a_s_CD4_nonday0_BCL6+Tfh_cluster_gene.xlsx")

table1 <- 'CD4/tonsil_LAIV_120a_s_CD4_subcluster_gene.xlsx'
table1_DE_data <- read_excel(table1, sheet = "cluster10")
table1_DE_data <- table1_DE_data[,!colnames(table1_DE_data) == '...1']
p <- match(table1_DE_data$gene, all_marker_table$marker)
temp <- all_marker_table[p,]
table1_DE_data$ENTREZID <- temp$ENTREZID
write.xlsx(table1_DE_data,sheetName = 'cluster10_copy',table1,append = T)

############# Tcell TCR mapping ##########################################
temp <- read.csv('tonsil_LAIV_120a_s_LW_Tcell_project_cell_index_clonal_expansion.csv')
p <- match(CD4.subset$project_cellindex,temp$project_cellindex)
temp <- temp[p,]
CD4.subset$clonal_expansion <- temp$clonal_expanded
DimPlot(CD4.subset, label = FALSE, reduction = "umap",group.by = 'clonal_expansion', pt.size = 0.1)#order = T,
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_umap_clonal_expansion_2.pdf',sep = ''),width = 6, height = 4)

########### CD4 gating ###########################################
CD4_dataframe_lintegrated_scale <- data.frame(scale(t(as.matrix(CD4.subset@assays[["integrated"]]@data))),check.names = F)
# CD4.subset$AbCD45RA_AbCD45RO_pos <- 'CD45RA+CD45RO- CD4'
# CD4.subset$AbCD45RA_AbCD45RO_pos[(as.matrix(CD4.subset@assays[["log1p"]]['ab-CD45RO']) > (2/3*as.matrix(CD4.subset@assays[["log1p"]]['ab-CD45RA']) + 1))] <- 'CD45RA-CD45RO+ CD4'
# CD4.subset$AbCD45RA_AbCD45RO_pos[((as.matrix(CD4.subset@assays[["log1p"]]['ab-CD45RO'])) > 3.2) & ((as.matrix(CD4.subset@assays[["log1p"]]['ab-CD45RA'])) > 3.4)] <- 'CD45RA+CD45RO+ CD4'
CD4.subset$AbCD45RA_AbCD45RO_pos <- 'naive CD4'
CD4.subset$AbCD45RA_AbCD45RO_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD45RO']) > (2/3*as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD45RA']) - 2.5))] <- 'memory CD4'
CD4.subset$AbCD45RA_AbCD45RO_pos[((as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD45RO'])) > -1) & ((as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD45RA'])) > 1.8)] <- 'stem memory CD4'

CD4.subset$AbCCR7_pos <- 'AbCCR7-CD4'
CD4.subset$AbCCR7_pos[as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CCR7']) > 0] <- 'AbCCR7+CD4'

CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos <- CD4.subset$AbCD45RA_AbCD45RO_pos
CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CCR7']) > 0) & (CD4.subset$AbCD45RA_AbCD45RO_pos == 'memory CD4')] <- 'CM'
CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CCR7']) <= 0) & (CD4.subset$AbCD45RA_AbCD45RO_pos == 'memory CD4')] <- 'EM'
CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CCR7']) <= 0) & (CD4.subset$AbCD45RA_AbCD45RO_pos == 'naive CD4')] <- 'TEMRA'

CD4.subset$AbCD62L_pos <- 'AbCD62L-CD4'
CD4.subset$AbCD62L_pos[as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD62L']) > 0] <- 'AbCD62L+CD4'

CD4.subset$AbCD62L_AbCD45RA_AbCD45RO_pos <- CD4.subset$AbCD45RA_AbCD45RO_pos
CD4.subset$AbCD62L_AbCD45RA_AbCD45RO_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD62L']) > 0) & (CD4.subset$AbCD45RA_AbCD45RO_pos == 'memory CD4')] <- 'CM'
CD4.subset$AbCD62L_AbCD45RA_AbCD45RO_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD62L']) <= 0) & (CD4.subset$AbCD45RA_AbCD45RO_pos == 'memory CD4')] <- 'EM'
CD4.subset$AbCD62L_AbCD45RA_AbCD45RO_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD62L']) <= 0) & (CD4.subset$AbCD45RA_AbCD45RO_pos == 'naive CD4')] <- 'TEMRA'

# naive CD4 follows CCR7 instead of CD62L, in particular for day0 cluster 28. 
CD4.subset$cluster_memory <- CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos#'EM'
# CD4.subset$cluster_memory[(CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos == 'CM') & (CD4.subset$AbCD62L_AbCD45RA_AbCD45RO_pos == 'CM')] <- 'CM'
# CD4.subset$cluster_memory[CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos == 'naive CD4'] <- 'naive CD4'
# CD4.subset$cluster_memory[CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos == 'stem memory CD4'] <- 'stem memory CD4'
# CD4.subset$cluster_memory[CD4.subset$AbCCR7_AbCD45RA_AbCD45RO_pos == 'TEMRA'] <- 'TEMRA'

# CD4.subset$FOXP3_pos <- 'CD4'
# CD4.subset$FOXP3_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['FOXP3']) > 0.5)] <- 'Treg'
CD4.subset$FOXP3_pos <- 'FOXP3-CD4'
CD4.subset$FOXP3_pos[(as.matrix(CD4.subset@assays[["log1p"]]['FOXP3']) > 0)] <- 'FOXP3+CD4'

CD4.subset$AbCD25_FOXP3_pos <- 'CD4'
CD4.subset$AbCD25_FOXP3_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['FOXP3']) > 0.5) & (as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD25']) > 0)] <- 'Treg'

CD4.subset$AbCD161_pos <- 'CD4'
CD4.subset$AbCD161_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD161']) > 0.25)] <- 'MAIT'

CD4.subset$BCL6_pos <- 'BCL6-CD4'
# CD4.subset$BCL6_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['BCL6']) > 0)] <- 'BCL6+CD4'# not much difference
CD4.subset$BCL6_pos[(as.matrix(CD4.subset@assays[["log1p"]]['BCL6']) > 0)] <- 'BCL6+CD4'

CD4.subset$BCL6_AbCD161_pos <- 'CD4'
CD4.subset$BCL6_AbCD161_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD161']) > 0.25) & (as.matrix(CD4.subset@assays[["log1p"]]['BCL6']) <= 0)] <- 'CD161+CD4'
CD4.subset$BCL6_AbCD161_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD161']) <= 0.25) & (as.matrix(CD4.subset@assays[["log1p"]]['BCL6']) > 0)] <- 'Tfh'
CD4.subset$BCL6_AbCD161_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD161']) > 0.25) & (as.matrix(CD4.subset@assays[["log1p"]]['BCL6']) > 0)] <- 'BCL6+CD161+CD4'

CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_pos <- CD4.subset$BCL6_AbCD161_pos
CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_pos[(CD4.subset$BCL6_AbCD161_pos == 'CD161+CD4') & (CD4.subset$AbCXCR5AbPD1_pos == 'AbCXCR5+CD4')] <- 'CXCR5+CD161+CD4'
CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_pos[(CD4.subset$BCL6_AbCD161_pos == 'CD4') & (CD4.subset$AbCXCR5AbPD1_pos == 'AbCXCR5+CD4')] <- 'BCL6-CXCR5+CD4'

CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_FOXP3_pos <- CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_pos
CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_FOXP3_pos[(CD4.subset$FOXP3_pos == 'FOXP3+CD4') & (CD4.subset$BCL6_pos == 'BCL6+CD4')] <- 'Tfh Treg'
CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_FOXP3_pos[(CD4.subset$FOXP3_pos == 'FOXP3+CD4') & (CD4.subset$BCL6_pos != 'BCL6+CD4')] <- 'Treg'

# CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_FOXP3_pos[(CD4.subset$FOXP3_pos == 'FOXP3+CD4') <- gsub('CD4','Treg',CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_FOXP3_pos[CD4.subset$FOXP3_pos == 'FOXP3+CD4'])
# CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_FOXP3_pos[(CD4.subset$FOXP3_pos == 'FOXP3+CD4') & (CD4.subset$BCL6_AbCXCR5AbPD1_AbCD161_FOXP3_pos == 'Tfh')] <- 'Tfh Treg'

# CD4.subset$FOXP3_BCL6_AbCD161_pos <- paste(gsub('CD4','',CD4.subset$FOXP3_pos),CD4.subset$BCL6_AbCD161_pos,sep = '')

CD4.subset$AbCXCR5_pos <- 'AbCXCR5-CD4'
CD4.subset$AbCXCR5_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0)] <- 'AbCXCR5+CD4'

CD4.subset$AbCXCR5_AbCD161_pos <- 'CD4'
CD4.subset$AbCXCR5_AbCD161_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD161']) > 0.25) & (as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CXCR5']) <= 0.5)] <- 'CD161+CD4'
CD4.subset$AbCXCR5_AbCD161_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD161']) <= 0.25) & (as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0.5)] <- 'CXCR5+CD4'
CD4.subset$AbCXCR5_AbCD161_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD161']) > 0.25) & (as.matrix(CD4.subset@assays[["integrated_scale"]]['ab=CXCR5']) > 0.5)] <- 'CXCR5+CD161+CD4'

CD4.subset$AbCXCR5AbPD1_pos <- 'AbCXCR5-CD4'
CD4.subset$AbCXCR5AbPD1_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & (as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD279']) > 0)] <- 'AbCXCR5+CD4'

CD4.subset$AbCXCR5AbPD1_BCL6_pos <- paste(substring(CD4.subset$AbCXCR5AbPD1_pos,first = 3, last = 8), CD4.subset$BCL6_pos, sep = '')
CD4.subset$differentiation_BCL6_pos <- paste(substring(CD4.subset$BCL6_pos,first = 1, last = 5), CD4.subset$differentiation, sep = '')
CD4.subset$differentiation2_AbCXCR5AbPD1_BCL6_pos <- paste(gsub('CD4','',CD4.subset$AbCXCR5AbPD1_BCL6_pos), CD4.subset$differentiation, sep = '')

# functional 
CD4.subset$BAFF_pos <- 'BAFF-CD4'
CD4.subset$BAFF_pos[(as.matrix(CD4.subset@assays[["log1p"]]['TNFSF13B']) > 0)] <- 'BAFF+CD4' 
# CD4.subset$BAFF_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['TNFSF13B']) > 0.7)] <- 'BAFF+CD4'

CD4.subset$differentiation3_BAFF_pos <- paste(substring(CD4.subset$BAFF_pos,first = 1, last = 5), CD4.subset$differentiation3, sep = '')

CD4.subset$VEGFA_pos <- 'VEGFA-CD4'
CD4.subset$VEGFA_pos[(as.matrix(CD4.subset@assays[["log1p"]]['VEGFA']) > 0)] <- 'VEGFA+CD4'
# CD4.subset$VEGFA_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['VEGFA']) > 0.7)] <- 'VEGFA+CD4'

CD4.subset$differentiation3_VEGFA_pos <- paste(substring(CD4.subset$VEGFA_pos,first = 1, last = 6), CD4.subset$differentiation3, sep = '')

CD4.subset$CXCR3_pos <- 'CXCR3-CD4'
CD4.subset$CXCR3_pos[(as.matrix(CD4.subset@assays[["log1p"]]['CXCR3']) > 0)] <- 'CXCR3+CD4'# not much difference
CD4.subset$CCR6_pos <- 'CCR6-CD4'
CD4.subset$CCR6_pos[(as.matrix(CD4.subset@assays[["log1p"]]['CCR6']) > 0)] <- 'CCR6+CD4'# not much difference
CD4.subset$CXCR3_CCR6_pos <- paste(substring(CD4.subset$CXCR3_pos,first = 1, last = 6),CD4.subset$CCR6_pos,sep = '')

CD4.subset$differentiation3_CXCR3_CCR6_pos <- paste(substring(CD4.subset$CXCR3_CCR6_pos,first = 1, last = 11), CD4.subset$differentiation3, sep = '')

CD4.subset$IL21_pos <- 'IL21-CD4'
CD4.subset$IL21_pos[(as.matrix(CD4.subset@assays[["log1p"]]['IL21']) > 0)] <- 'IL21+CD4'
CD4.subset$differentiation3_IL21_pos <- paste(substring(CD4.subset$IL21_pos,first = 1, last = 5), CD4.subset$differentiation3, sep = '')

CD4.subset$IL12_pos <- 'IL12-CD4'
CD4.subset$IL12_pos[(as.matrix(CD4.subset@assays[["log1p"]]['IL12A']) > 0)] <- 'IL12+CD4'
CD4.subset$differentiation3_IL12_pos <- paste(substring(CD4.subset$IL12_pos,first = 1, last = 5), CD4.subset$differentiation3, sep = '')

CD4.subset$IL4_pos <- 'IL4-CD4'
CD4.subset$IL4_pos[(as.matrix(CD4.subset@assays[["log1p"]]['IL4']) > 0)] <- 'IL4+CD4'

CD4.subset$IL2_pos <- 'IL2-CD4'
CD4.subset$IL2_pos[(as.matrix(CD4.subset@assays[["log1p"]]['IL2']) > 0)] <- 'IL2+CD4'
CD4.subset$differentiation3_IL2_pos <- paste(substring(CD4.subset$IL2_pos,first = 1, last = 4), CD4.subset$differentiation3, sep = '')

CD4.subset$TBET_pos <- 'TBET-CD4'
CD4.subset$TBET_pos[(as.matrix(CD4.subset@assays[["log1p"]]['TBX21']) > 0)] <- 'TBET+CD4'
CD4.subset$GATA3_pos <- 'GATA3-CD4'
CD4.subset$GATA3_pos[(as.matrix(CD4.subset@assays[["log1p"]]['GATA3']) > 0)] <- 'GATA3+CD4'
CD4.subset$RORC_pos <- 'RORC-CD4'
CD4.subset$RORC_pos[(as.matrix(CD4.subset@assays[["log1p"]]['RORC']) > 0)] <- 'RORC+CD4'

CD4.subset$TBET_RORC_pos <- paste(gsub('CD4','',CD4.subset$TBET_pos), CD4.subset$RORC_pos, sep = '')
CD4.subset$GATA3_RORC_pos <- paste(gsub('CD4','',CD4.subset$GATA3_pos), CD4.subset$RORC_pos, sep = '')

CD4.subset$CXCL13_pos <- 'CD4'
CD4.subset$CXCL13_pos[(as.matrix(CD4.subset@assays[["log1p"]]['CXCL13']) > 0)] <- 'CXCL13+CD4'

CD4.subset$ICOS_pos <- 'ICOS-CD4'
CD4.subset$ICOS_pos[(as.matrix(CD4.subset@assays[["log1p"]]['ICOS']) > 0)] <- 'ICOS+CD4'

CD4.subset$AbPD1_pos <- 'PD1-CD4'
CD4.subset$AbPD1_pos[(as.matrix(CD4.subset@assays[["integrated_scale"]]['ab-CD279']) > 0)] <- 'PD1+CD4'

CD4.subset$AbCXCR5_AbPD1_pos <- paste(gsub('Ab','',gsub('CD4','',CD4.subset$AbCXCR5_pos)), CD4.subset$AbPD1_pos, sep = '')

CD4.subset$AbCD38_pos <- 'AbCD38-CD4'
CD4.subset$AbCD38_pos[(CD4_dataframe_lintegrated_scale$`ab-CD38` > 0)] <- 'AbCD38+CD4'
CD4.subset$AbHLADR_pos <- 'AbHLA-DR-CD4'
CD4.subset$AbHLADR_pos[(CD4_dataframe_lintegrated_scale$`ab-HLA-DR` > 0)] <- 'AbHLA-DR+CD4'
CD4.subset$AbCD38_AbHLADR_pos <- paste(gsub('CD4','',CD4.subset$AbCD38_pos), CD4.subset$AbHLADR_pos, sep = '')
CD4.subset$AbCD38_AbHLADR_pos[CD4.subset$AbCD38_AbHLADR_pos != 'AbCD38+AbHLA-DR+CD4'] <- 'AbCD38-/AbHLA-DR-CD4'

CD4.subset$CD38_pos <- 'CD38-CD4'
CD4.subset$CD38_pos[(as.matrix(CD4.subset@assays[["log1p"]]['CD38']) > 0)] <- 'CD38+CD4'
CD4.subset$HLADRA_pos <- 'HLA-DR-CD4'
CD4.subset$HLADRA_pos[(as.matrix(CD4.subset@assays[["log1p"]]['HLA-DRA']) > 0)] <- 'HLA-DR+CD4'
CD4.subset$CD38_HLADRA_pos <- paste(gsub('CD4','',CD4.subset$CD38_pos), CD4.subset$HLADRA_pos, sep = '')
CD4.subset$CD38_HLADRA_pos[CD4.subset$CD38_HLADRA_pos != 'CD38+HLA-DR+CD4'] <- 'CD38-/HLA-DR-CD4'
CD4.subset$AbCD38_HLADR_pos <- paste(gsub('CD4','',CD4.subset$AbCD38_pos), CD4.subset$HLADRA_pos, sep = '')
CD4.subset$AbCD38_HLADR_pos[CD4.subset$AbCD38_HLADR_pos != 'AbCD38+HLA-DR+CD4'] <- 'AbCD38-/HLA-DR-CD4'

CD4.subset$CD69_pos <- 'CD69-CD4'
CD4.subset$CD69_pos[(as.matrix(CD4.subset@assays[["log1p"]]['CD69']) > 0)] <- 'CD69+CD4'
CD4.subset$CD40L_pos <- 'CD40L-CD4'
CD4.subset$CD40L_pos[(as.matrix(CD4.subset@assays[["log1p"]]['CD40LG']) > 0)] <- 'CD40L+CD4'
CD4.subset$CD69_CD40L_pos <- paste(gsub('CD4','',CD4.subset$CD69_pos), CD4.subset$CD40L_pos, sep = '')
CD4.subset$CD69_CD40L_pos[CD4.subset$CD69_CD40L_pos != 'CD69+CD40L+CD4'] <- 'CD69-/CD40L-CD4'

CD4.subset$AbCD69_pos <- 'AbCD69-CD4'
CD4.subset$AbCD69_pos[(CD4_dataframe_lintegrated_scale['ab-CD69'] > 0)] <- 'AbCD69+CD4'
CD4.subset$AbCD40L_pos <- 'AbCD154-CD4'
CD4.subset$AbCD40L_pos[(CD4_dataframe_lintegrated_scale['ab-CD154'] > 0)] <- 'AbCD154+CD4'
CD4.subset$AbCD69_AbCD40L_pos <- paste(gsub('CD4','',CD4.subset$AbCD69_pos), CD4.subset$AbCD40L_pos, sep = '')
CD4.subset$AbCD69_AbCD40L_pos[CD4.subset$AbCD69_AbCD40L_pos != 'AbCD69+AbCD154+CD4'] <- 'AbCD69-/AbCD154-CD4'

CD4.subset$differentiation3_AbCD38_AbHLADR_pos <- paste(gsub('CD4','',CD4.subset$AbCD38_AbHLADR_pos), CD4.subset$differentiation3, sep = '')

DimPlot(CD4.subset, label = FALSE, reduction = "umap",group.by = 'AbCD38_AbHLADR_pos', order = T,pt.size = 0.1)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_umap_AbCD38_AbHLADR_pos.pdf',sep = ''),width = 600, height = 400)

CD4_log1p_dataframe <- data.frame(t(as.matrix(CD4.subset@assays[["log1p"]]@data)),check.names = F)
CD4_integrated_scale_dataframe <- data.frame(t(as.matrix(CD4.subset@assays[["integrated_scale"]]@data)),check.names = F)

CD4_integrated_scale_dataframe$AbCD45RA_AbCD45RO_pos <- CD4.subset$AbCD45RA_AbCD45RO_pos
ggplot(CD4_integrated_scale_dataframe, aes(x = !!sym(c('ab-CD45RA')), y = !!sym(c('ab-CD45RO')))) + ggtitle('CD4 cells') +#!!sym(c('TBX21'))
  geom_point(alpha = 0.1,size = 0.5,aes(color = !!sym(c('AbCD45RA_AbCD45RO_pos')))) +
  # geom_point(alpha = 0.1,size = 0.5) +
  # scale_colour_manual(values = c("naive CD4" = "#F8766D", "CD4" = '#00BFC4'))+
  # scale_colour_manual(values = c("CD45RA+CD83+ nonPB B" = "gold", "CD184-CD83- nonPB B" = "blue","CD184+CD83- nonPB B" = 'cyan',"CD184-CD83+ nonPB B" = 'magenta')) +
  geom_density_2d(color = 'black')# +
# ylim(0,8) + xlim(0,8)
# geom_hline(aes(yintercept=3.25,color = 'red')) + guides(color = 'none') +
# geom_vline(aes(xintercept=4.1,color = 'red')) + guides(color = 'none')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD4_AbCD45RA_AbCD45RO_pos_contour_integrated_scale.pdf',width = 7, height = 4.6)

DefaultAssay(CD4.subset) <- 'integrated_scale'
temp <- subset(CD4.subset,cluster_labeled == 'naive CD4')
FeatureScatter(temp, feature1 = "ab-CD45RA", feature2 = "ab-CCR7") 
dev.print(pdf, 'tonsil_LAIV_120a_s_CD4_scatter_abCD45RA_abCD45RO.pdf',width = 7, height = 4.34)

### cytokine analysis ###########################################
cytokine_global_table <- read_excel('../Global landscape of cytokines Supplementary Table S1.xlsx', sheet = "Cytokines")
cytokine_global_list <- unique(cytokine_global_table$`HGNC symbol`)
cytokine_global_list <- c(cytokine_global_list,'GZMB','PRF1','MIF')
cytokine_Rhapsody_list <- all.markers[all.markers %in% cytokine_global_list]

for (cytokine_name in cytokine_Rhapsody_list){
  graphics.off()
  plot <- FeaturePlot(CD4.subset, features = c(cytokine_name),pt.size = 0.2, ncol = 1, sort.cell = TRUE,min.cutoff = 0)#, max.cutoff = 15)
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_umap_',cytokine_name,'.pdf',sep = ''),width = 5, height = 4.3)
  
}



### all CD4s: cell cluster fraction visualization #####################
selected_day_group_list <- c('day04-06','day07-08','day10','day12-14')
subset_name <- 'differentiation'
CD4.database <- data.frame(cbind(CD4.subset$donor_ID,as.character(CD4.subset@meta.data[[subset_name]]),CD4.subset$condition,CD4.subset$days,CD4.subset$stimulation,CD4.subset$age_group,CD4.subset$day_group,CD4.subset$age))
colnames(CD4.database) <- c('donor_ID',subset_name,'condition','day','stimulation','age_group','day_group','age')
# sapply(CD4.database,class)
CD4_cluster_condition_donor_count <- count(CD4.database, donor_ID, !!sym(subset_name),condition, day,stimulation,age_group,day_group,age)
CD4_cluster_condition_donor_count <- CD4_cluster_condition_donor_count %>% group_by(donor_ID, condition, day,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
CD4_cluster_condition_donor_count$day <- as.numeric(substring(CD4_cluster_condition_donor_count$day,first = 4,last = 5))
CD4_cluster_condition_donor_count$percentage <- CD4_cluster_condition_donor_count$n/CD4_cluster_condition_donor_count$total*100
CD4_cluster_condition_donor_count$age <- as.numeric(CD4_cluster_condition_donor_count$age)
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(CD4.subset$donor_ID)){
  temp_condition_list <- sort(unique(CD4.subset$condition[CD4.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- CD4_cluster_condition_donor_count[(CD4_cluster_condition_donor_count$donor_ID == donor_name) & (CD4_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(CD4.subset@meta.data[[subset_name]])){
      if_row <- ((CD4_cluster_condition_donor_count$donor_ID == donor_name) & (CD4_cluster_condition_donor_count[,subset_name] == cluster_name) & (CD4_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount[,subset_name] <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        CD4_cluster_condition_donor_count[nrow(CD4_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(CD4.database)

# statistics & visualization
CD4_cluster_condition_donor_count_day0 <- CD4_cluster_condition_donor_count[CD4_cluster_condition_donor_count$day == '0',]
CD4_cluster_condition_donor_count_double_day0 <- CD4_cluster_condition_donor_count
temp <- CD4_cluster_condition_donor_count_double_day0[CD4_cluster_condition_donor_count_double_day0$day == '0',]
CD4_cluster_condition_donor_count_double_day0$stimulation[CD4_cluster_condition_donor_count_double_day0$stimulation == 'day0'] <- 'LAIV+'
temp$stimulation[temp$stimulation == 'day0'] <- 'LAIV-'
CD4_cluster_condition_donor_count_double_day0 <- rbind(CD4_cluster_condition_donor_count_double_day0,temp)
CD4_cluster_condition_donor_count_double_day0 <- CD4_cluster_condition_donor_count_double_day0[CD4_cluster_condition_donor_count_double_day0$donor_ID != "02yrs M IMD085",]
p <- match(CD4_cluster_condition_donor_count_double_day0$condition, flow_count_table$condition)
temp <- flow_count_table[p,]
CD4_cluster_condition_donor_count_double_day0$CD4_count <- temp$CD4_count
annotate_cluster_list <- unique(CD4.subset@meta.data[[subset_name]])
# annotate_cluster_list <- annotate_cluster_list[grepl('\\+',annotate_cluster_list)]# | grepl('CXCR5\\+PD1\\-CD4',annotate_cluster_list) | grepl('CXCR5\\-PD1\\+CD4',annotate_cluster_list)]
day_ttest_age_group <- data.frame()
day_correlate_age <- data.frame()
for (cluster_name in annotate_cluster_list)
{
  print(cluster_name)
  temp1 <- CD4_cluster_condition_donor_count[CD4_cluster_condition_donor_count[,subset_name] == cluster_name,]
  temp1$age <- as.numeric(substring(temp1$donor_ID,first = 1, last = 2))
  temp1$day_group <- 'day02'
  # temp1$day_group[temp1$day == 4] <- 'day04'
  temp1$day_group[(temp1$day >= 4) & (temp1$day <= 6)] <- 'day04-06'
  temp1$day_group[(temp1$day >= 7) & (temp1$day <= 8)] <- 'day07-08'
  temp1$day_group[(temp1$day >= 10) & (temp1$day <= 10)] <- 'day10'
  temp1$day_group[temp1$day >= 12] <- 'day12-14'
  # temp1$day_group2 <- temp1$day_group
  # temp1$day_group2[temp1$day_group2 %in% c('day04','day05-06')] <- 'day04-06'
  temp1$age_group <- 2
  temp1$age_group[temp1$age < 5] <- 1
  temp1$age_group[temp1$age > 18] <- 3
  for (day_group_name in day_group_list){
    temp2 <- temp1[temp1$day %in% day_group_dict[[day_group_name]],]
    if (day_group_name != 'day0') {
      temp_LAIV <- temp2[temp2$stimulation == 'LAIV+',]
      temp_ns <- temp2[temp2$stimulation == 'LAIV-',]
      corr_result <- cor.test(temp_LAIV$percentage, temp_LAIV$age)
      temp <- data.frame(subset = cluster_name, day_group = day_group_name, stimulation = 'LAIV+', corr = corr_result$estimate, p_val = corr_result$p.value)
      day_correlate_age <- rbind(day_correlate_age,temp)
      corr_result <- cor.test(temp_ns$percentage, temp_ns$age)
      temp <- data.frame(subset = cluster_name, day_group = day_group_name, stimulation = 'LAIV-', corr = corr_result$estimate, p_val = corr_result$p.value)
      day_correlate_age <- rbind(day_correlate_age,temp)
      for (age_group_index in c(1:3)) {
        if ((sum(temp_LAIV$age_group == age_group_index) >= 2) & (sum(temp_LAIV$age_group != age_group_index) >= 2)) {
          ttest_result <- t.test(temp_LAIV$percentage[temp_LAIV$age_group == age_group_index],
                                 temp_LAIV$percentage[temp_LAIV$age_group != age_group_index])
          temp <- data.frame(subset = cluster_name,
                             day_group = day_group_name,
                             age_group = age_group_list[age_group_index], stimulation = 'LAIV+',
                             diff = ttest_result[["estimate"]][["mean of x"]] - ttest_result[["estimate"]][["mean of y"]], 
                             p_val = ttest_result$p.value)
          day_ttest_age_group <- rbind(day_ttest_age_group,temp)
        }
        if ((sum(temp_ns$age_group == age_group_index) >= 2) & (sum(temp_ns$age_group != age_group_index) >= 2)) {
          ttest_result <- t.test(temp_ns$percentage[temp_ns$age_group == age_group_index],
                                 temp_ns$percentage[temp_ns$age_group != age_group_index])
          temp <- data.frame(subset = cluster_name,
                             day_group = day_group_name,
                             age_group = age_group_list[age_group_index], stimulation = 'LAIV-',
                             diff = ttest_result[["estimate"]][["mean of x"]] - ttest_result[["estimate"]][["mean of y"]], 
                             p_val = ttest_result$p.value)
          day_ttest_age_group <- rbind(day_ttest_age_group,temp)
        }
      }
    } else {
      corr_result <- cor.test(temp2$percentage, temp2$age)
      temp <- data.frame(subset = cluster_name, day_group = day_group_name, stimulation = 'day0', corr = corr_result$estimate, p_val = corr_result$p.value)
      day_correlate_age <- rbind(day_correlate_age,temp)
      for (age_group_index in c(1:3)) {
        if ((sum(temp2$age_group == age_group_index) >= 2) & (sum(temp2$age_group != age_group_index) >= 2)) {
          ttest_result <- t.test(temp2$percentage[temp2$age_group == age_group_index],
                                 temp2$percentage[temp2$age_group != age_group_index])
          temp <- data.frame(subset = cluster_name,
                             day_group = day_group_name, 
                             age_group = age_group_list[age_group_index], stimulation = 'day0',
                             diff = ttest_result[["estimate"]][["mean of x"]] - ttest_result[["estimate"]][["mean of y"]], 
                             p_val = ttest_result$p.value)
          day_ttest_age_group <- rbind(day_ttest_age_group,temp)
        }
      }
    }
  }
  graphics.off()
  cols_age_group <- c('#99CCFF','#3399FF','#003366')
  names(cols_age_group) <- c('02-04yrs','07-09yrs','27-39yrs')
  shape_age_group <- c(25,4,3)
  names(shape_age_group) <- c('02-04yrs','07-09yrs','27-39yrs')
  temp <- CD4_cluster_condition_donor_count[(CD4_cluster_condition_donor_count$day == '0') & (CD4_cluster_condition_donor_count[,subset_name] == cluster_name),]
  temp$age <- as.numeric(substring(temp$donor_ID,first = 1, last = 2))
  temp$age_group <- '07-09yrs'
  temp$age_group[temp$age < 5] <- '02-04yrs'
  temp$age_group[temp$age > 18] <- '27-39yrs'
  plot <- ggplot(temp,aes(x=age, y=percentage)) +
    ggtitle(paste('uncultured',cluster_name)) +
    geom_point(size = 1,aes(fill = age_group, shape = age_group)) + 
    theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) + RotatedAxis() +
    ylab(paste('CD4 (%)')) +
    geom_smooth(method='lm', formula= y~x,color = 'red') +
    annotate("text", x=(temp$age[which.max(temp$percentage)] + 25)%%40, y=mean(sort(temp$percentage)[(length(temp$percentage) - 1): length(temp$percentage)]), label = paste("r =",sprintf(day_correlate_age$corr[(day_correlate_age$day_group == 'day0') & (day_correlate_age$subset == cluster_name)], fmt = '%#.2f')),size = 3) +
    scale_fill_manual(values = cols_age_group) + 
    scale_shape_manual(values = shape_age_group)
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_',cluster_name,'_condition_ratio_day0.pdf',sep = ''),width = 3, height = 2)
  temp <- CD4_cluster_condition_donor_count_double_day0[CD4_cluster_condition_donor_count_double_day0[,subset_name] == cluster_name,]
  plot <- ggplot(temp,aes(x=day, y=percentage,color = stimulation)) +  ggtitle(cluster_name) +
    facet_wrap( ~ donor_ID, scales = "fixed",ncol = 3) +
    # facet_grid(age_group~donor_ID, scales = "free_x",ncol = 3) +
    geom_point(size = 1)+ geom_line() + 
    theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) + RotatedAxis() +
    ylab('CD4 %') +
    scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_',cluster_name,'_condition_ratio2.pdf',sep = ''),width = 5, height = 4)
  plot <- ggplot(temp,aes(x=day, y=percentage,color = stimulation, shape = donor_ID)) +  ggtitle(cluster_name) +
    facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
    # facet_grid(age_group~donor_ID, scales = "free_x") +
    geom_point(size = 5)+ geom_line(aes(group = interaction(donor_ID,stimulation))) + 
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('CD4 %') +
    scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
    # scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) +
    scale_shape_manual(values=1:length(unique(CD4_cluster_condition_donor_count_double_day0$donor_ID)))
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 15, height = 4)
  temp$cluster_count <- temp$CD4_count*temp$percentage/100
  plot <- ggplot(temp[temp$donor_ID != "07yrs M IMD170",],aes(x=day, y=cluster_count,color = stimulation, shape = donor_ID)) +  ggtitle(cluster_name) +
    facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
    # facet_grid(age_group~donor_ID, scales = "free_x") +
    geom_point(size = 5)+ geom_line(aes(group = interaction(donor_ID,stimulation))) + 
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('#cells') +
    scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
    # scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) +
    scale_shape_manual(values=c(1:5,7:length(unique(CD4_cluster_condition_donor_count_double_day0$donor_ID))))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_',cluster_name,'_condition_count.pdf',sep = ''),width = 15, height = 4)
  
}

day_ttest_age_group <- day_ttest_age_group %>% arrange(p_val)
write.xlsx(day_ttest_age_group,paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_day_age_ttest_group.xlsx',sep = ''))
day_correlate_age <- day_correlate_age %>% arrange(desc(abs(corr)))
write.xlsx(day_correlate_age,paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_day_age_corr.xlsx',sep = ''))


######## all cell diff visualization ###########################
# LAIV_condition_list <- unique(ALL_Rhapsody_TCR_clonotype_table$sample_condition)
# condition_paired_table <- data.frame(LAIV = LAIV_condition_list[grepl('\\+',LAIV_condition_list)])
# condition_paired_table$ns <- gsub('\\+','\\-',condition_paired_table$LAIV)
# condition_paired_table <- condition_paired_table[(condition_paired_table$LAIV %in% LAIV_condition_list) & (condition_paired_table$ns %in% LAIV_condition_list),]
# condition_paired_table$donor_ID <- substring(condition_paired_table$LAIV,first = 1, last = 14) 
# condition_paired_table$age <- as.numeric(as.numeric(substring(condition_paired_table$LAIV,first = 1, last = 2)))
# condition_paired_table$day <- as.numeric(gsub('day','',gsub('.*(day.*) .*','\\1',condition_paired_table$LAIV)))
# condition_paired_table$age_group <- '07-09yrs'
# condition_paired_table$age_group[condition_paired_table$age < 5] <- '02-04yrs'
# condition_paired_table$age_group[condition_paired_table$age > 18] <- '27-39yrs'
# condition_paired_table$donorID_day <- substring(condition_paired_table$LAIV, first = 1, last = 20)

# One should do chi-square test!!
nonday0_diff_ttest_age_group <- data.frame()
nonday0_diff_correlate_age <- data.frame()

for (cluster_name in annotate_cluster_list) {
  print(cluster_name)
  graphics.off()
  temp <- CD4_cluster_condition_donor_count[CD4_cluster_condition_donor_count[,subset_name] == cluster_name,]
  rownames(temp) <- temp$condition
  if ((!all(temp$percentage == 0)) & (!all(temp$percentage == 100))) {
    condition_paired_cluster_table <- condition_paired_table
    condition_paired_cluster_table$LAIV_cluster <- unlist(temp[condition_paired_table$LAIV,][,'percentage'])
    condition_paired_cluster_table$ns_cluster <- unlist(temp[condition_paired_table$ns,][,'percentage'])
    # condition_paired_cluster_table$FC_cluster <- condition_paired_cluster_table$LAIV_cluster/condition_paired_cluster_table$ns_cluster
    # condition_paired_cluster_table$FC_cluster[condition_paired_cluster_table$ns_cluster == 0] <- NaN
    # condition_paired_cluster_table$log2FC_cluster <- log2(condition_paired_cluster_table$FC_cluster)
    # condition_paired_cluster_table$log2FC_cluster[condition_paired_cluster_table$FC_cluster == 0] <- NaN
    condition_paired_cluster_table$diff_cluster <- condition_paired_cluster_table$LAIV_cluster - condition_paired_cluster_table$ns_cluster
    condition_paired_cluster_table$day_group <- 'day02'
    # condition_paired_cluster_table$day_group[condition_paired_cluster_table$day == 4] <- 'day04'
    condition_paired_cluster_table$day_group[(condition_paired_cluster_table$day >= 4) & (condition_paired_cluster_table$day <= 6)] <- 'day04-06'
    condition_paired_cluster_table$day_group[(condition_paired_cluster_table$day >= 7) & (condition_paired_cluster_table$day <= 8)] <- 'day07-08'
    condition_paired_cluster_table$day_group[(condition_paired_cluster_table$day >= 10) & (condition_paired_cluster_table$day <= 10)] <- 'day10'
    condition_paired_cluster_table$day_group[condition_paired_cluster_table$day >= 12] <- 'day12-14'
    # condition_paired_cluster_table$day_group2 <- condition_paired_cluster_table$day_group
    # condition_paired_cluster_table$day_group2[condition_paired_cluster_table$day_group2 %in% c('day04','day05-06')] <- 'day04-06'
    
    temp1 <- condition_paired_cluster_table
    temp1$age_group2 <- 2
    temp1$age_group2[temp1$age < 5] <- 1
    temp1$age_group2[temp1$age > 18] <- 3
    for (day_group_name in day_group_list[2:18]) {
      temp2 <- temp1[temp1$day %in% day_group_dict[[day_group_name]],]
      corr_result <- cor.test(temp2$diff_cluster, temp2$age)
      temp <- data.frame(subset = cluster_name,
                         day_group = day_group_name,
                         corr = corr_result$estimate, p_val = corr_result$p.value)
      nonday0_diff_correlate_age <- rbind(nonday0_diff_correlate_age,temp)
      for (age_group_index in c(1:3)) {
        if ((sum(temp2$age_group2 == age_group_index) >= 2) & (sum(temp2$age_group2 != age_group_index) >= 2)) {
          ttest_result <- t.test(temp2$diff_cluster[temp2$age_group2 == age_group_index],
                                 temp2$diff_cluster[temp2$age_group2 != age_group_index])
          temp <- data.frame(subset = cluster_name,
                             day_group = day_group_name,
                             age_group = age_group_list[age_group_index], 
                             diff = ttest_result[["estimate"]][["mean of x"]] - ttest_result[["estimate"]][["mean of y"]], 
                             p_val = ttest_result$p.value)
          nonday0_diff_ttest_age_group <- rbind(nonday0_diff_ttest_age_group,temp)
        }  
      }
    }
    
    cols_age_group <- c('#99CCFF','#3399FF','#003366')
    names(cols_age_group) <- c('02-04yrs','07-09yrs','27-39yrs')
    shape_age_group <- c(25,4,3)
    names(shape_age_group) <- c('02-04yrs','07-09yrs','27-39yrs')
    
    plot <- ggplot(condition_paired_cluster_table,aes(x=day, y=diff_cluster,shape = age_group,color = age_group,group = donor_ID)) +
      geom_hline(yintercept = 0,color = 'red',linetype = 'dotted') +
      ggtitle(paste(cluster_name,'diff in LAIV and ns')) +
      geom_point(size = 1)+ geom_line() +
      theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) + RotatedAxis() +
      ylab(paste('diff in CD4 (%)')) +
      scale_color_manual(values = cols_age_group) +
      scale_shape_manual(values = shape_age_group)
    print(plot)
    dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_nonday0_diff_trend_',cluster_name,'.pdf',sep = ''),width = 4, height = 2.5)
    
    plot <- ggplot(condition_paired_cluster_table,aes(x=day, y=diff_cluster,shape = age_group,fill = age_group, group = donor_ID)) +
      geom_hline(yintercept = 0,color = 'red',linetype = 'dotted') +
      ggtitle(paste(cluster_name,'diff in LAIV and ns')) +
      geom_point(size = 2) + 
      theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) + RotatedAxis() +
      ylab(paste('diff in CD4 (%)')) +
      scale_fill_manual(values = cols_age_group) + 
      scale_shape_manual(values = shape_age_group)
    print(plot)
    dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_nonday0_diff_',cluster_name,'.pdf',sep = ''),width = 4, height = 2.5)
    

    temp <- condition_paired_cluster_table[condition_paired_cluster_table$day_group %in% selected_day_group_list,]
    temp1 <- (   temp
                 %>% group_by( day_group )
                 %>% dplyr::slice( which.max( diff_cluster ) )
                 %>% as.data.frame
    )
    temp2 <- data.frame(diff_cluster = temp1$diff_cluster,
                           age = (temp1$age + 23)%%40,
                           label = sprintf(nonday0_diff_correlate_age$corr[(nonday0_diff_correlate_age$subset == cluster_name) & (nonday0_diff_correlate_age$day_group %in% selected_day_group_list)], fmt = '%#.2f'),
                           day_group = nonday0_diff_correlate_age$day[(nonday0_diff_correlate_age$subset == cluster_name) & (nonday0_diff_correlate_age$day_group %in% selected_day_group_list)])
    plot <- ggplot(temp, aes(y=diff_cluster,x=age)) +
      facet_wrap(~day_group, scales = "free_x",ncol = 4) +
      geom_point() + RotatedAxis() +
      geom_smooth(method='lm', formula= y~x, color = 'red') +
      ylab('diff of CD4 %\nbetween LAIV+ and LAIV-') +
      ggtitle(cluster_name)+
      geom_text(data = temp2, aes(label = label),size = 4)
    print(plot)
    dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_nonday0_diff_day_group_',cluster_name,'.pdf',sep = ''),width = 5, height = 3)
    
  }
}
nonday0_diff_ttest_age_group <- nonday0_diff_ttest_age_group %>% arrange(p_val)
write.xlsx(nonday0_diff_ttest_age_group,paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_nonday0_diff_day_age_ttest_group.xlsx',sep = ''))
nonday0_diff_correlate_age <- nonday0_diff_correlate_age %>% arrange(desc(abs(corr)))
write.xlsx(nonday0_diff_correlate_age,paste('tonsil_LAIV_120a_s_CD4_',subset_name,'_nonday0_diff_day_age_corr.xlsx',sep = ''))


####### subcluster dynamic ranking ###############################
library(ggridges)
library(reshape2)
library(mixtools)
CD4_subcluster_stimulation_mean_dynamics <- data.frame()
for (cluster_name in unique(CD4_cluster_condition_donor_count_double_day0$differentiation2)) {
  temp <- CD4_cluster_condition_donor_count_double_day0[(CD4_cluster_condition_donor_count_double_day0$differentiation2 == cluster_name),]
  graphics.off()
  temp1 <- CD4_cluster_condition_donor_count_double_day0[(CD4_cluster_condition_donor_count_double_day0$differentiation2 == cluster_name) & (CD4_cluster_condition_donor_count_double_day0$stimulation == 'LAIV+'),]
  temp2 <- CD4_cluster_condition_donor_count_double_day0[(CD4_cluster_condition_donor_count_double_day0$differentiation2 == cluster_name) & (CD4_cluster_condition_donor_count_double_day0$stimulation == 'LAIV-'),]
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
    CD4_subcluster_stimulation_mean_dynamics <- rbind(CD4_subcluster_stimulation_mean_dynamics,temp1)
  } else {
    CD4_subcluster_stimulation_mean_dynamics <- rbind(CD4_subcluster_stimulation_mean_dynamics,temp2)
  }
  
  # this figure cannot eliminate points negative since it is doing the fitting automatically
    # plot <- ggplot(temp,aes(x=days,y=percentage, color = stimulation)) +  ggtitle(cluster_name) +
    #   geom_point(aes(x=days,y=percentage, shape = donor_ID), size = 5)+
    #   geom_line(aes(x=days,y=percentage, shape = donor_ID,group = interaction(donor_ID,stimulation))) +
    #   theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    #   ylab('CD4 %') +
    #   scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "grey45")) +
    #   scale_shape_manual(values=1:length(unique(CD4_cluster_condition_donor_count_double_day0$donor_ID)))+
    #   stat_smooth(method = "loess", formula = y ~ x, size = 1)
    # print(plot)
    # dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_smooth_pattern_differentiation2_',cluster_name,'_LAIV_ns.pdf',sep = ''),width = 7, height = 4)
  
}
CD4_subcluster_stimulation_mean_dynamics_norm <- CD4_subcluster_stimulation_mean_dynamics - apply(CD4_subcluster_stimulation_mean_dynamics, 1, min)
CD4_subcluster_stimulation_mean_dynamics_norm <- CD4_subcluster_stimulation_mean_dynamics_norm/apply(CD4_subcluster_stimulation_mean_dynamics_norm, 1, max)

# # if peak has more than one.
# sum(CD4_subcluster_stimulation_mean_dynamics_norm > 0.95)
# sum(CD4_subcluster_stimulation_mean_dynamics_norm == 1)
# temp <- temp %>% arrange(value)
# temp$variable <- c(1:dim(temp)[1])
# ggplot(temp,aes(x = variable,y = value)) + geom_point() + geom_hline(yintercept = 0.95) # no clear cutoff

temp <- data.frame(CD4_subcluster_stimulation_mean_dynamics_norm)
temp <- temp[temp$X0 == 1,]
temp$cluster <- rownames(temp)
temp <- melt(temp,id = 'cluster', value.name = 'norm_pct')
temp$days <- as.numeric(substring(temp$variable,first = 2, last = 4))
ggplot(temp, aes(x=days,y=norm_pct,color = cluster)) +  geom_line()
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_smooth_pattern_differentiation2_norm_day0.pdf',sep = ''),width = 7, height = 4)

# CD4_subcluster_mean_cumsum <- data.frame(t(apply(temp, 1, cumsum))) # smoothed cumulative sum
# CD4_subcluster_mean_cumsum_norm <- CD4_subcluster_mean_cumsum/CD4_subcluster_mean_cumsum$X14
# CD4_subcluster_mean_cumsum_norm$cluster <- rownames(CD4_subcluster_mean_cumsum_norm)
# temp <- melt(CD4_subcluster_mean_cumsum_norm,id = 'cluster', value.name = 'smoothed_percentage_cumsum')
# temp$days <- as.numeric(substring(temp$variable,first = 2, last = 4))
# ggplot(temp, aes(x=days,y=smoothed_percentage_cumsum,color = cluster)) +  geom_line()
# dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_smooth_pattern_differentiation2_cumsum_norm_day0.pdf',sep = ''),width = 7, height = 4)

temp <- data.frame(which(CD4_subcluster_stimulation_mean_dynamics_norm == 1,arr.ind=TRUE))
temp$day_peak <- c(0:15)[temp$col]
temp$cluster_stimulation <- rownames(temp)
temp$lambda1 <- 0
temp$mu1 <- 0
temp$lambda2 <- 0
temp$mu2 <- 0
temp$sigma1 <- 0
temp$sigma2 <- 0
graphics.off()
for (cluster_name in rownames(CD4_subcluster_stimulation_mean_dynamics_norm)){
  print(cluster_name)
  temp1 <- data.frame(norm_pct = as.numeric(CD4_subcluster_stimulation_mean_dynamics_norm[cluster_name,]),day = 0:14)
  plot(as.numeric(temp1$norm_pct))
  if (cluster_name %in% c("EM CD161+CD4 LAIV")){
    fit <- nls(norm_pct ~ (lambda1 * exp(-(day-mu1)**2/(2 * sigma1**2)) +
                             lambda2 * exp(-(day-mu2)**2/(2 * sigma2**2))), 
               data=temp1, algorithm="port",
               start=list(lambda1=0.5, mu1=8, sigma1=2,lambda2=2, mu2=10, sigma2=10),
               low = c(0,-40,0,0,-40,0),upper =c(4,40,40,4,40,40)) # cannot work it out using only positive lambda!
              # start=list(lambda1=1, mu1=1, sigma1=1,lambda2=1, mu2=13, sigma2=1), 
              # low = c(-2,-15,0,-2,-15,0),upper =c(2,30,20,2,30,20))
    
  } else if (cluster_name %in% c("CM non-Tfh CD4 LAIV")) {
    fit <- nls(norm_pct ~ (lambda1 * exp(-(day-mu1)**2/(2 * sigma1**2)) +
                             lambda2 * exp(-(day-mu2)**2/(2 * sigma2**2))), 
               data=temp1, algorithm="port",
               start=list(lambda1=0.5, mu1=4, sigma1=5,lambda2=0.5, mu2=10, sigma2=5),
               low = c(0,-40,0,0,-40,0),upper =c(4,40,40,4,40,40)) # cannot work it out using only positive lambda!
    # start=list(lambda1=1, mu1=1, sigma1=1,lambda2=1, mu2=13, sigma2=1), 
    # low = c(-2,-15,0,-2,-15,0),upper =c(2,30,20,2,30,20))
  } else if (cluster_name %in% c("EM non-Tfh CD4 ns")) {
    fit <- nls(norm_pct ~ (lambda1 * exp(-(day-mu1)**2/(2 * sigma1**2)) +
                             lambda2 * exp(-(day-mu2)**2/(2 * sigma2**2))), 
               data=temp1, algorithm="port",
               start=list(lambda1=1, mu1=0, sigma1=2.5,lambda2=1, mu2=14, sigma2=5),
               low = c(0,-10,0,0,14,0),upper =c(4,0,20,4,20,20)) # cannot work it out using only positive lambda!
    # start=list(lambda1=1, mu1=1, sigma1=1,lambda2=1, mu2=13, sigma2=1), 
    # low = c(-2,-15,0,-2,-15,0),upper =c(2,30,20,2,30,20))
  } else {
    fit <- nls(norm_pct ~ (lambda1 * exp(-(day-mu1)**2/(2 * sigma1**2)) +
                             lambda2 * exp(-(day-mu2)**2/(2 * sigma2**2))), 
               data=temp1, algorithm="port",
               start=list(lambda1=0.5, mu1=4, sigma1=2,lambda2=0.5, mu2=10, sigma2=2),
               low = c(0,-40,0,0,-40,0),upper =c(4,40,40,4,40,40)) # cannot work it out using only positive lambda!
    # start=list(lambda1=1, mu1=1, sigma1=1,lambda2=1, mu2=13, sigma2=1), 
    # low = c(-2,-15,0,-2,-15,0),upper =c(2,30,20,2,30,20))
  }
  fit
  temp1 <- data.frame(coef(fit))
  temp1$variable <- gsub('[0-9]','',rownames(temp1))
  temp1$order <- gsub('.*([0-9])','\\1',rownames(temp1))
  temp1 <- dcast( temp1 , order ~ variable,value.var = 'coef.fit.' )
  temp1$order <- NULL
  temp1 <- temp1 %>% plyr::arrange(sigma)
  temp$lambda1[temp$cluster == cluster_name] <- temp1$lambda[1]
  temp$lambda2[temp$cluster == cluster_name] <- temp1$lambda[2]
  temp$sigma1[temp$cluster == cluster_name] <- temp1$sigma[1]
  temp$sigma2[temp$cluster == cluster_name] <- temp1$sigma[2]
  temp$mu1[temp$cluster == cluster_name] <- temp1$mu[1]
  temp$mu2[temp$cluster == cluster_name] <- temp1$mu[2]
  temp$modified_mu1[temp$cluster == cluster_name] <- temp1$mu[1]
  temp$modified_mu2[temp$cluster == cluster_name] <- temp1$mu[2]
  if (temp1$lambda[1] < 0) {
    if (temp$mu1[temp$cluster == cluster_name] < temp$mu2[temp$cluster == cluster_name]) {
      temp$modified_mu1[temp$cluster == cluster_name] <- temp$modified_mu1[temp$cluster == cluster_name] + temp1$sigma[1]
    } else {
      temp$modified_mu1[temp$cluster == cluster_name] <- temp$modified_mu1[temp$cluster == cluster_name] - temp1$sigma[1]
    }
  }
  if (temp1$lambda[2] < 0) {
    if (temp$modified_mu2[temp$cluster == cluster_name] < temp$modified_mu1[temp$cluster == cluster_name]) {
      temp$modified_mu2[temp$cluster == cluster_name] <- temp$modified_mu2[temp$cluster == cluster_name] + temp1$sigma[2]
    } else {
      temp$modified_mu2[temp$cluster == cluster_name] <- temp$modified_mu2[temp$cluster == cluster_name] - temp1$sigma[2]
    }
  }
  if (temp$modified_mu1[temp$cluster == cluster_name] < 0) {
    temp$modified_mu1[temp$cluster == cluster_name] <- 0
  }
  if (temp$modified_mu1[temp$cluster == cluster_name] > 14) {
    temp$modified_mu1[temp$cluster == cluster_name] <- 14
  }
  if (temp$modified_mu2[temp$cluster == cluster_name] < 0) {
    temp$modified_mu2[temp$cluster == cluster_name] <- 0
  }
  if (temp$modified_mu2[temp$cluster == cluster_name] > 14) {
    temp$modified_mu2[temp$cluster == cluster_name] <- 14
  }
  if (abs(temp$modified_mu1[temp$cluster == cluster_name] - temp$day_peak[temp$cluster == cluster_name]) < 1) {
    temp1 <- c(temp$lambda1[temp$cluster == cluster_name],temp$modified_mu1[temp$cluster == cluster_name],temp$sigma1[temp$cluster == cluster_name])
    temp$lambda1[temp$cluster == cluster_name] <- temp$lambda2[temp$cluster == cluster_name]
    temp$modified_mu1[temp$cluster == cluster_name] <- temp$modified_mu2[temp$cluster == cluster_name]
    temp$sigma1[temp$cluster == cluster_name] <- temp$sigma2[temp$cluster == cluster_name]
    temp$lambda2[temp$cluster == cluster_name] <- temp1[1]
    temp$modified_mu2[temp$cluster == cluster_name] <- temp1[2]
    temp$sigma2[temp$cluster == cluster_name] <- temp1[3]
  }
  if (temp$modified_mu1[temp$cluster == cluster_name] < temp$day_peak[temp$cluster == cluster_name]) {
    temp$change[temp$cluster == cluster_name] <- temp$modified_mu1[temp$cluster == cluster_name] + temp$sigma1[temp$cluster == cluster_name]
  } else {
    temp$change[temp$cluster == cluster_name] <- temp$modified_mu1[temp$cluster == cluster_name] - temp$sigma1[temp$cluster == cluster_name]
  }
  
}

temp$cluster <- gsub(' ns','',gsub(' LAIV','',temp$cluster_stimulation))
temp$stimulation <- gsub('.* ','',temp$cluster_stimulation)
temp$row <- NULL

temp <- temp %>% arrange(day_peak,modified_mu1)
# temp2 <- dcast(temp,   cluster ~ stimulation, value.var = 'day_peak')
# temp2$trend <- apply(temp2[,c('LAIV','ns')], 1, min)
write.xlsx(temp,'tonsil_LAIV_120a_s_CD4_smooth_pattern_differentiation2_LAIV_ns_dominant.xlsx')

########## cell subset geom_flow visualization ####################
ordered_differentiation2_subset <- c("naive non-Tfh CD4","TEMRA","stem memory CD4","CM non-Tfh CD4","CM CXCR5+BCL6-CD4","EM non-Tfh CD4","Treg",
                                     "CM Tfh","EM Tfh","naive Tfh","EM CXCR5+BCL6-CD4",
                                     "CM CD161+CD4","EM CD161+CD4")
CD4_cluster_condition_donor_count_double_day0$differentiation2 <- factor(CD4_cluster_condition_donor_count_double_day0$differentiation2,
                                                                         levels = ordered_differentiation2_subset)
temp <- merge(x = CD4_cluster_condition_donor_count_double_day0, y = data.frame(differentiation2 = ordered_differentiation2_subset, order = c(1:nb.cols)), by = "differentiation2", all.x=TRUE)
temp <- temp %>% arrange(order)

library(RColorBrewer)
# subset_colors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
nb.cols <- length(unique(CD4_cluster_condition_donor_count_double_day0$differentiation2))
subset_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
names(subset_colors) <- ordered_differentiation2_subset
for (donor_name in unique(CD4_cluster_condition_donor_count_double_day0$donor_ID)){
  graphics.off()
  print(donor_name)
  plot <- ggplot(temp[temp$donor_ID == donor_name,], aes(x = days, stratum = stimulation, alluvium = differentiation2, y = percentage, label = stimulation,fill = differentiation2)) +#, 
    geom_alluvium(aes(fill = differentiation2), width = 1) +
    facet_grid(stimulation ~ .,scales = "fixed") +
    geom_stratum(width = 1/12, fill = "black") +
    scale_fill_manual(values = subset_colors) +
    ggtitle(donor_name) + ylab('CD4 %') +
    theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_',donor_name,'_differentiation2.pdf',sep = ''),width = 6, height = 3.5)
}
  
############ CD8 #####################################################
CD8.subset <- subset(Tcell.subset,meta_cluster == 'CD8')
CD8.subset <- RunPCA(CD8.subset, features = all.markers, npcs = 50)
CD8.subset <- RunUMAP(CD8.subset, reduction = "pca", dims = 1:n_pca_selected)
CD8.subset <- FindNeighbors(CD8.subset, reduction = "pca", dims = 1:n_pca_selected)
CD8.subset <- FindClusters(CD8.subset, resolution = 2)
DimPlot(CD8.subset, label = TRUE, reduction = "umap",repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_umap_subcluster.pdf',width = 7, height = 4.34)
DimPlot(CD8.subset, label = TRUE, reduction = "umap",group.by = 'project_ID')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_umap_projectID.pdf',width = 7, height = 4.34)
DimPlot(CD8.subset, label = TRUE, reduction = "umap",group.by = 'donor_ID')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_umap_donorID.pdf',width = 7, height = 4.34)

DimPlot(CD8.subset, label = TRUE, reduction = "umap",group.by = 'days')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_umap_days.pdf',width = 7, height = 4.34)
CD8.subcluster.markers <- FindAllMarkers(CD8.subset)
CD8.subcluster.markers <- CD8.subcluster.markers %>% filter(p_val_adj <= 0.05)
CD8.subcluster.markers.top10 <- CD8.subcluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(CD8.subset, features = c(CD8.subcluster.markers.top10$gene)) +
  theme(text = element_text(size=20))
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_subcluster_gene_heatmap.pdf',width = 12, height = 25)
write.csv(CD8.subcluster.markers, "tonsil_LAIV_120a_s_CD8_subcluster_gene.csv")
write.csv(CD8.subcluster.markers.top10, "tonsil_LAIV_120a_s_CD8_subcluster_gene_top10.csv")

saveRDS(CD8.subset,'tonsil_LAIV_120a_s_CD8.rds')
CD8.subset$AbCD45RA_AbCD45RO_pos <- 'CD45RA+CD45RO- CD8'
CD8.subset$AbCD45RA_AbCD45RO_pos[as.matrix(CD8.subset@assays[["log1p"]]['ab-CD45RO']) > 2/3*(as.matrix(CD8.subset@assays[["log1p"]]['ab-CD45RA']) )] <- 'CD45RA-CD45RO+ CD8'
# CD8.subset$AbCD45RA_AbCD45RO_pos[as.matrix(CD8.subset@assays[["log1p"]]['ab-CD45RO']) > (as.matrix(CD8.subset@assays[["log1p"]]['ab-CD45RA']) - 2)] <- 'CD45RA-CD45RO+ CD8'
CD8.subset$AbCD45RA_AbCD45RO_pos[((as.matrix(CD8.subset@assays[["log1p"]]['ab-CD45RO'])) > 3) & ((as.matrix(CD8.subset@assays[["log1p"]]['ab-CD45RA'])) > 4)] <- 'CD45RA+CD45RO+ CD8'

CD8.subset$AbCCR7_AbCD45RA_AbCD45RO_pos <- paste('ab-CCR7-',CD8.subset$AbCD45RA_AbCD45RO_pos,sep = '')
CD8.subset$AbCCR7_AbCD45RA_AbCD45RO_pos[((as.matrix(CD8.subset@assays[["log1p"]]['ab-CCR7'])) > 1.3)] <- paste('ab-CCR7+',CD8.subset$AbCD45RA_AbCD45RO_pos[((as.matrix(CD8.subset@assays[["log1p"]]['ab-CCR7'])) > 1.3)],sep = '')

CD8.subset$TCF7_AbCD45RA_AbCD45RO_pos <- paste('TCF7-',CD8.subset$AbCD45RA_AbCD45RO_pos,sep = '')
CD8.subset$TCF7_AbCD45RA_AbCD45RO_pos[((as.matrix(CD8.subset@assays[["log1p"]]['TCF7'])) > 0.5)] <- paste('TCF7+',CD8.subset$AbCD45RA_AbCD45RO_pos[((as.matrix(CD8.subset@assays[["log1p"]]['TCF7'])) > 0.5)],sep = '')

CD8.subset$BCL6_pos <- 'BCL6-CD8'
CD8.subset$BCL6_pos[((as.matrix(CD8.subset@assays[["log1p"]]['BCL6'])) > 0.5)] <- 'BCL6+CD8'

CD8.subset$CXCR5_pos <- 'CXCR5-CD8'
CD8.subset$CXCR5_pos[((as.matrix(CD8.subset@assays[["log1p"]]['CXCR5'])) > 0.5)] <- 'CXCR5+CD8'

CD8.subset$AbCXCR5_pos <- 'ab-CXCR5-CD8'
CD8.subset$AbCXCR5_pos[((as.matrix(CD8.subset@assays[["log1p"]]['ab-CXCR5'])) > 2.7)] <- 'ab-CXCR5+CD8'

CD8.subset$AbCD161_pos <- 'abCD161-CD8'
CD8.subset$AbCD161_pos[(as.matrix(CD8.subset@assays[["log1p"]]['ab-CD161']) > 2.5)] <- 'abCD161+CD8'

CD8.subset$AbCD161_level <- 'abCD161-CD8'
CD8.subset$AbCD161_level[(as.matrix(CD8.subset@assays[["log1p"]]['ab-CD161']) > 2.5)] <- 'abCD161mid CD8'
CD8.subset$AbCD161_level[(as.matrix(CD8.subset@assays[["log1p"]]['ab-CD161']) > 4)] <- 'abCD161++CD8'

CD8.subset$IFNG_pos <- 'IFNG-CD8'
CD8.subset$IFNG_pos[((as.matrix(CD8.subset@assays[["log1p"]]['IFNG'])) > 0.5)] <- 'IFNG+CD8'

DimPlot(CD8.subset, label = FALSE, reduction = "umap",group.by = 'AbCD161_pos')#,order = c("CD27+BCL6+ B","CD27-BCL6+ B","nonPB B")) + 
# scale_colour_manual(values = c("CD27+BCL6+ B" = "red", "nonPB B" = 'grey',"CD27-BCL6+ B" = 'cyan'))
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD8_cleaned_umap_AbCD161_pos.pdf',sep = ''),width = 600, height = 400)

# DimPlot(CD8.subset, label = FALSE, reduction = "umap",group.by = 'CD184_CD83_pos')+#,order = c("CD83+ B","CXCR4+ B")) + 
#   # scale_colour_manual(values = c("CXCR4+ B" = "#F8766D", "CXCR4- B" = '#00BFC4'))
#   scale_colour_manual(values = c("CD184+CD83+ nonPB B" = "gold", "CD184-CD83- nonPB B" = "blue","CD184+CD83- nonPB B" = 'cyan',"CD184-CD83+ nonPB B" = 'magenta'))
# dev.print(pdf, paste('tonsil_LAIV_120a_s_CD8_cleaned_umap_CD184_CD83_pos.pdf',sep = ''),width = 600, height = 400)

CD8_log1p_dataframe <- data.frame(t(as.matrix(CD8.subset@assays[["log1p"]]@data)),check.names = F)
CD8_log1p_scale_dataframe <- data.frame(t(as.matrix(CD8.subset@assays[["log1p_scale"]]@data)),check.names = F)

CD8_log1p_dataframe$AbCXCR5_pos <- CD8.subset$AbCXCR5_pos
ggplot(CD8_log1p_dataframe, aes(x = !!sym(c('ab-CD161')), y = !!sym(c('ab-CXCR5')))) + ggtitle('CD8 cells') +
  geom_point(alpha = 0.1,size = 0.5,aes(color = !!sym(c('AbCXCR5_pos')))) +
  # geom_point(alpha = 0.1,size = 0.5) +
  # scale_colour_manual(values = c("naive CD8" = "#F8766D", "CD8" = '#00BFC4'))+
  # scale_colour_manual(values = c("CD45RA+CD83+ nonPB B" = "gold", "CD184-CD83- nonPB B" = "blue","CD184+CD83- nonPB B" = 'cyan',"CD184-CD83+ nonPB B" = 'magenta')) +
  # ylim(0,8) + xlim(0,8) +
  geom_density_2d(color = 'black')
# geom_hline(aes(yintercept=3.25,color = 'red')) + guides(color = 'none') +
# geom_vline(aes(xintercept=4.1,color = 'red')) + guides(color = 'none')
dev.print(pdf, 'tonsil_LAIV_120a_s_CD8_AbCD161_AbCXCR5_contour_log1p.pdf',width = 5, height = 4.6)

CD8.subset <- subset(Tcell.subset,meta_cluster == 'CD8')
saveRDS(CD8.subset,'tonsil_LAIV_120a_s_CD8.rds')

CD8_LAIV.subset <- subset(CD8.subset, stimulation != 'LAIV-')
CD8_LAIV.subset <- subset(CD8_LAIV.subset,days != 'day02')
CD8_LAIV.subset <- subset(CD8_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
CD8_LAIV.subset <- subset(CD8_LAIV.subset,donor_ID != "02yrs M IMD085")

temp <- subset(CD8_LAIV.subset, days %in% c('day0','day04','day06','day07','day10','day12','day14'))
temp$day_group <- 'day0'
temp$day_group[temp$days %in% c('day04')] <- 'day04'
temp$day_group[temp$days %in% c('day06','day07')] <- 'day06-07'
temp$day_group[temp$days %in% c('day10')] <- 'day10'
temp$day_group[temp$days %in% c('day12','day14')] <- 'day12-14'
gene_name <- 'KLRC3'
VlnPlot(temp, features = c(gene_name),pt.size = 0.1, ncol = 1,split.by = 'age_group',log = F,assay = 'normalized',group.by = 'day_group')
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD8_violin_normalized_age_group_',gene_name,'.pdf',sep = ''),width = 15, height = 6)

gene_name <- 'GNLY'
View(CD8.subcluster.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))
c((CD8.subcluster.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))['cluster'])
# VlnPlot(CD8.subset, features = c(gene_name),pt.size = 0.2, ncol = 1,sort = 'increasing',log = F)# + geom_hline(yintercept = -0.7)
# dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_violin_ab-CD38.pdf',width = 12, height = 4.34)
FeaturePlot(CD8.subset, features = c(gene_name),pt.size = 0.2, ncol = 1, sort.cell = TRUE,min.cutoff = 0, max.cutoff = 3)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD8_umap_',gene_name,'.pdf',sep = ''),width = 5, height = 4.3)
VlnPlot(Tcell.subset, features = gene_name,pt.size = 0.2,group.by = 'donor_ID',assay = 'integrated_scale',slot = 'data',sort = 'descreasing')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_violin_subcluster_labeled_donors_nCount.pdf',sep = ''),width = 10, height =5)
VlnPlot(Tcell.subset, features = gene_name,pt.size = 0.2,split.by = 'donor_ID',group.by = 'subcluster_manual_label',assay = 'integrated_scale',slot = 'data',sort = 'descreasing')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_violin_subcluster_labeled_donors_nCount.pdf',sep = ''),width = 10, height =5)
VlnPlot(Tcell.subset, features = gene_name,pt.size = 0.2,group.by = 'cluster_manual_label',assay = 'integrated_scale',slot = 'data',sort = 'descreasing')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_violin_cluster_labeled_nCount.pdf',sep = ''),width = 10, height =5)
# cluster_subset <- subset(Tcell.subset,donor_ID == '02yrs M IMD170')
VlnPlot(Tcell.subset, features = gene_name,pt.size = 0.2,split.by = 'treatment',group.by = 'donor_ID',assay = 'integrated_scale',slot = 'data',sort = 'increasing')
# VlnPlot(Tcell.subset, features = gene_name,pt.size = 0.2,split.by = 'treatment',group.by = 'cluster_manual_label',assay = 'RNA',slot = 'data',sort = 'increasing')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_violin_',gene_name,'_cluster_labeled_log1p.pdf',sep = ''),width = 10, height =5)

# ab-CD45RA: 26 6  18 25 9  21 27 24 2  8  13
# ab-CD45RO: 19 1  20 0  22 5  12 14 17
# ab-CCR7: 18 4  26 6  21 24 9  25 27
# ab-CD62L: 18 19 26 6  21
# ab-CD127: 18 24 26 6  9  13 8 
# CCR7: 18 21 26 24 6  25 27 9 
# SELL: 23 12 24 19 11 16 6  22 34 17
# TCF7: 21 24 26 13 18 25 6  8  9  2  27
# IL7R: 18 24 3  6  13 8  26
# LEF1: 26 27 25 21 24 18 9  6  8  2  13

# IKZF2: 19 2
# ab-CXCR5: 20 7  5  4  1  14 11
# CXCR5: 14 5  20 1  7  4  11
# BCL6: 24 21
# ab-CD161: 0  22 3  15 19 16 17 12
# KLRB1: 0  22 19 3  15 12 17 16
# RORC: 22 19 0  3  12

# MKI67:7  29
# GAPDH: 7  29 18 11 15
# CD80LG: 16 11 7 
# ADGRE1: 21 7  3  15 11 4 
# RORC: 7 9
# ab-CD161: 9  10 6  5  23 25 29
# IL17A: 9  7  10 27 6  12 5 
# IL17F: 18 7  9  11 27
# SELPLG: 
# GATA3: 15 29 7  8 
# TBX21: 29 15 16 25 7  26 23 9 
# IFNG: 7  16 9  29 19 8 

# FOXP3: 15 22
# IL10: 15 11 6  16
# ab-CD25: 15 18 7  22 6  29 9 
# IL2RA: 12 16 25 2 

# CXCR5: 0  14 1  2  11
# ab-CXCR5: 14 1  20 24 0  12 3  13
# BCL6: 11 16 0  14 2  21
# CD69: 28 30 8  16
# DUSP1: 8  30 28
# ITGAE: 7  29 23 18 25 21 0 

# ab-CD38: 6  9  18 7  29 21 5  15 10 14 12
# ab-HLA-DR: 14 18 7  9  6  15
# CD38: 25 12
# HLA-DRA: 12 25 2 

# TGFBI: 1  0  13
# TNF: 16 7  11
# EGR3: 16 11
# GZMB: 29 25 23 19 9  26 7 
# GZMK: 19 29 25 26 23
# EOMES: 19 29 25 23 26
# B3GAT1(CD57): 11 19 8  0  16 20 2 (align with the TFH)
# IL21: 16 7  18 11
# PDCD1: 20 24 11 7  2 (align with TFH)
# ab-CD279(PD1): 0  11 6  1  16 2  7  everything except 4,28,3,27,22
# ab-Tim3: 9  29 15 6  14 24 7  18 23 16
# CTLA4: 15 16 8  7  6  18 11
# ARG1 (checkpoint): 22 15 29 11 4  0 not much
# TNFRSF4 (OX40): 11 20 22
# TNFRSF8 (CD30): 7  0  15
# TNFSF8: 11 0  16 7 
# CXCL13: 0  11 18 7  16 2 

# MX1: 7  18 0  29 21 11
# IFIT3: 7  14 29 9  21
# TNFSF13B: 
# GZMA: 
# GZMK: 
# GNLY: 
# NKG7: 
# VEGFA: 2  29 11 15 16
# ALDOC: 2  18 11 14 7  0  15
# TIGIT 
# CD200: 
# CXCL13: 
# EGR1: 
# EGR3: 8
# ADGRG1 (GPR56): 2  6  14 4
### CD8 cell cluster fraction visualization #####################
CD8.database <- data.frame(cbind(CD8.subset$donor_ID,as.character(CD8.subset$IFNG_pos),CD8.subset$condition,CD8.subset$days,CD8.subset$treatment,CD8.subset$responder_group,CD8.subset$age_group,CD8.subset$condition_group,CD8.subset$LAIV_response))
colnames(CD8.database) <- c('donor_ID','IFNG_pos','condition','days','treatment','responder_group','age_group','condition_group','LAIV_response')
# sapply(CD8.database,class)
CD8_cluster_condition_donor_count <- count(CD8.database, donor_ID, IFNG_pos,condition, days,treatment,responder_group,age_group,condition_group,LAIV_response)
CD8_cluster_condition_donor_count <- CD8_cluster_condition_donor_count %>% group_by(donor_ID, condition,days,treatment,responder_group,age_group,condition_group,LAIV_response) %>% mutate(total = sum(n))
CD8_cluster_condition_donor_count$percentage <- CD8_cluster_condition_donor_count$n/CD8_cluster_condition_donor_count$total*100
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(Tcell.subset$donor_ID)){
  condition_list <- sort(unique(Tcell.subset$condition[Tcell.subset$donor_ID == donor_name]))
  for (condition_name in condition_list){
    temp_Bcount <- CD8_cluster_condition_donor_count[(CD8_cluster_condition_donor_count$donor_ID == donor_name) & (CD8_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(CD8.subset$IFNG_pos)){
      if_row <- ((CD8_cluster_condition_donor_count$donor_ID == donor_name) & (CD8_cluster_condition_donor_count$IFNG_pos == cluster_name) & (CD8_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Bcount <- CD8_cluster_condition_donor_count[CD8_cluster_condition_donor_count$donor_ID == donor_name,][1,]
        temp_Bcount$IFNG_pos <- cluster_name
        temp_Bcount$n <- 0
        temp_Bcount$percentage <- 0
        temp_Bcount$condition <- condition_name
        temp_Bcount$days <- Tcell.subset$days[Tcell.subset$condition == condition_name][1]
        temp_Bcount$treatment <- Tcell.subset$treatment[Tcell.subset$condition == condition_name][1]
        CD8_cluster_condition_donor_count[nrow(CD8_cluster_condition_donor_count) + 1,] <- temp_Bcount
      }
    }
  }
}
# CD8_cluster_condition_donor_count$responder_group <- factor(CD8_cluster_condition_donor_count$responder_group,levels = c("non-LAIV-responder","LAIV-responder F","LAIV-responder M","LAIV-responder"))

# write.csv(CD8_cluster_condition_donor_count, "tonsil_Rhapsody_CD8_cluster_percentage.csv")
# two figures:
# 1. only responders
# 2. all donors but only at day0
TB_cluster_condition_donor_count_day0 <- CD8_cluster_condition_donor_count[CD8_cluster_condition_donor_count$days == 'day0',]
TB_cluster_condition_donor_count_double_day0 <- CD8_cluster_condition_donor_count
temp <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$days == 'day0',]
TB_cluster_condition_donor_count_double_day0$treatment[TB_cluster_condition_donor_count_double_day0$treatment == 'day0'] <- 'LAIV+'
temp$treatment[temp$treatment == 'day0'] <- 'LAIV-'
TB_cluster_condition_donor_count_double_day0 <- rbind(TB_cluster_condition_donor_count_double_day0,temp)
TB_cluster_condition_donor_count_double_day0 <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$condition_group == 'all timepoints',]
for (cluster_name in unique(CD8.subset$IFNG_pos))
{
  print(cluster_name)
  subdata_allTP <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$IFNG_pos == cluster_name,]
  plot <- ggplot(subdata_allTP,aes(x=days, y=percentage,color = treatment, shape = donor_ID)) +  ggtitle(cluster_name) +
    facet_wrap(~ responder_group + donor_ID, scales = "free_x",ncol = 3) +
    # facet_grid(gender~age_group, scales = "fixed") +
    geom_point(size = 5)+ geom_line(aes(group = interaction(donor_ID,treatment))) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    ylab('% in CD8 T cells') +
    # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) + 
    scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) + 
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_double_day0$donor_ID)))
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD8_log1p_scale_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 15, height = 9)
  subdata_LAIV <- subdata_allTP[subdata_allTP$treatment == 'LAIV+',] %>% arrange(condition)
  subdata_NS <- subdata_allTP[subdata_allTP$treatment == 'LAIV-',] %>% arrange(condition)
  subdata_LAIV$fold_change <- subdata_LAIV$percentage/subdata_NS$percentage
  graphics.off()
  plot <- ggplot(subdata_LAIV,aes(x=days, y=fold_change,shape = donor_ID)) +  ggtitle(cluster_name) +
    facet_wrap(~ responder_group + donor_ID, scales = "free_x",ncol = 3) +
    # facet_grid(gender~age_group, scales = "fixed") +
    geom_point(size = 5)+ geom_line(aes(group = interaction(donor_ID,treatment))) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    ylab('Fold-change of\n% in CD8 T cells in LAIV+ vs LAIV-') + geom_hline(aes(yintercept=1,color = 'red')) + guides(color = 'none') +
    # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) + 
    # scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) + 
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_double_day0$donor_ID)))
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD8_log1p_scale_',cluster_name,'_condition_FC.pdf',sep = ''),width = 15, height = 9)
  graphics.off()
  ttest_result <- t.test(TB_cluster_condition_donor_count_day0$percentage[(TB_cluster_condition_donor_count_day0$LAIV_response =='non-responder') & (TB_cluster_condition_donor_count_day0$IFNG_pos == cluster_name)],TB_cluster_condition_donor_count_day0$percentage[(TB_cluster_condition_donor_count_day0$LAIV_response =='responder') & (TB_cluster_condition_donor_count_day0$IFNG_pos == cluster_name)])
  plot <- ggplot(TB_cluster_condition_donor_count_day0[TB_cluster_condition_donor_count_day0$IFNG_pos == cluster_name,],aes(x=LAIV_response, y=percentage,color = LAIV_response, shape = donor_ID)) +  ggtitle(paste('day0 pval:',as.integer(10*ttest_result$p.value)/10,'\n',cluster_name)) +
    facet_grid(~age_group, scales = "free_x") +
    geom_point(size = 5) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    scale_colour_manual(values = c("non-responder" = "#00BFC4", "responder" = "#F8766D")) +
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_day0$donor_ID))) +
    ylab('% in CD8 T cells')
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD8_log1p_scale_day0_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 11, height = 7)
  
}

########### CD8 subset visualization in T cells ####################################
Tcell.subset$CD8_subset <- 'nonCD8 T cells'
Tcell.subset$CD8_subset[Tcell.subset$meta_cluster == 'CD8'] <- CD8.subset$IFNG_pos
Tcell.database <- data.frame(cbind(Tcell.subset$donor_ID,as.character(Tcell.subset$CD8_subset),Tcell.subset$condition,Tcell.subset$days,Tcell.subset$treatment,Tcell.subset$responder_group,Tcell.subset$age_group,Tcell.subset$condition_group,Tcell.subset$LAIV_response))
colnames(Tcell.database) <- c('donor_ID','CD8_subset','condition','days','treatment','responder_group','age_group','condition_group','LAIV_response')
# sapply(Tcell.database,class)
Tcell_cluster_condition_donor_count <- count(Tcell.database, donor_ID, CD8_subset,condition, days,treatment,responder_group,age_group,condition_group,LAIV_response)
Tcell_cluster_condition_donor_count <- Tcell_cluster_condition_donor_count %>% group_by(donor_ID, condition,days,treatment,responder_group,age_group,condition_group,LAIV_response) %>% mutate(total = sum(n))
Tcell_cluster_condition_donor_count$percentage <- Tcell_cluster_condition_donor_count$n/Tcell_cluster_condition_donor_count$total*100
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(Tcell.subset$donor_ID)){
  condition_list <- sort(unique(Tcell.subset$condition[Tcell.subset$donor_ID == donor_name]))
  for (condition_name in condition_list){
    temp_Bcount <- Tcell_cluster_condition_donor_count[(Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Tcell.subset$CD8_subset)){
      if_row <- ((Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$CD8_subset == cluster_name) & (Tcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Bcount <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$donor_ID == donor_name,][1,]
        temp_Bcount$CD8_subset <- cluster_name
        temp_Bcount$n <- 0
        temp_Bcount$percentage <- 0
        temp_Bcount$condition <- condition_name
        temp_Bcount$days <- Tcell.subset$days[Tcell.subset$condition == condition_name][1]
        temp_Bcount$treatment <- Tcell.subset$treatment[Tcell.subset$condition == condition_name][1]
        Tcell_cluster_condition_donor_count[nrow(Tcell_cluster_condition_donor_count) + 1,] <- temp_Bcount
      }
    }
  }
}
# Tcell_cluster_condition_donor_count$responder_group <- factor(Tcell_cluster_condition_donor_count$responder_group,levels = c("non-LAIV-responder","LAIV-responder F","LAIV-responder M","LAIV-responder"))

# write.csv(Tcell_cluster_condition_donor_count, "tonsil_Rhapsody_Tcell_cluster_percentage.csv")
# two figures:
# 1. only responders
# 2. all donors but only at day0
TB_cluster_condition_donor_count_day0 <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$days == 'day0',]
TB_cluster_condition_donor_count_double_day0 <- Tcell_cluster_condition_donor_count
temp <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$days == 'day0',]
TB_cluster_condition_donor_count_double_day0$treatment[TB_cluster_condition_donor_count_double_day0$treatment == 'day0'] <- 'LAIV+'
temp$treatment[temp$treatment == 'day0'] <- 'LAIV-'
TB_cluster_condition_donor_count_double_day0 <- rbind(TB_cluster_condition_donor_count_double_day0,temp)
TB_cluster_condition_donor_count_double_day0 <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$condition_group == 'all timepoints',]
for (cluster_name in unique(Tcell.subset$CD8_subset))
{
  print(cluster_name)
  subdata_allTP <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$CD8_subset == cluster_name,]
  plot <- ggplot(subdata_allTP,aes(x=days, y=percentage,color = treatment, shape = donor_ID)) +  ggtitle(cluster_name) +
    facet_wrap(~ responder_group + donor_ID, scales = "free_x",ncol = 3) +
    # facet_grid(gender~age_group, scales = "fixed") +
    geom_point(size = 5)+ geom_line(aes(group = interaction(donor_ID,treatment))) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    ylab('% in T cells') +
    # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) + 
    scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) + 
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_double_day0$donor_ID)))
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_log1p_scale_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 15, height = 9)
  subdata_LAIV <- subdata_allTP[subdata_allTP$treatment == 'LAIV+',] %>% arrange(condition)
  subdata_NS <- subdata_allTP[subdata_allTP$treatment == 'LAIV-',] %>% arrange(condition)
  subdata_LAIV$fold_change <- subdata_LAIV$percentage/subdata_NS$percentage
  graphics.off()
  plot <- ggplot(subdata_LAIV,aes(x=days, y=fold_change,shape = donor_ID)) +  ggtitle(cluster_name) +
    facet_wrap(~ responder_group + donor_ID, scales = "free_x",ncol = 3) +
    # facet_grid(gender~age_group, scales = "fixed") +
    geom_point(size = 5)+ geom_line(aes(group = interaction(donor_ID,treatment))) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    ylab('Fold-change of\n% in T cells in LAIV+ vs LAIV-') + geom_hline(aes(yintercept=1,color = 'red')) + guides(color = 'none') +
    # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) + 
    # scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) + 
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_double_day0$donor_ID)))
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_log1p_scale_',cluster_name,'_condition_FC.pdf',sep = ''),width = 15, height = 9)
  graphics.off()
  ttest_result <- t.test(TB_cluster_condition_donor_count_day0$percentage[(TB_cluster_condition_donor_count_day0$LAIV_response =='non-responder') & (TB_cluster_condition_donor_count_day0$CD8_subset == cluster_name)],TB_cluster_condition_donor_count_day0$percentage[(TB_cluster_condition_donor_count_day0$LAIV_response =='responder') & (TB_cluster_condition_donor_count_day0$CD8_subset == cluster_name)])
  plot <- ggplot(TB_cluster_condition_donor_count_day0[TB_cluster_condition_donor_count_day0$CD8_subset == cluster_name,],aes(x=LAIV_response, y=percentage,color = LAIV_response, shape = donor_ID)) +  ggtitle(paste('day0 pval:',as.integer(10*ttest_result$p.value)/10,'\n',cluster_name)) +
    facet_grid(~age_group, scales = "free_x") +
    geom_point(size = 5) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    scale_colour_manual(values = c("non-responder" = "#00BFC4", "responder" = "#F8766D")) +
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_day0$donor_ID))) +
    ylab('% in T cells')
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_log1p_scale_day0_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 11, height = 7)
  
}

############ gdT #####################################################
gdT.subset <- subset(Tcell.subset,meta_cluster == 'gdT')
saveRDS(gdT.subset,'tonsil_LAIV_120a_s_gdT.rds')
gdT.subset <- RunPCA(gdT.subset, features = all.markers, npcs = 50)
gdT.subset <- RunUMAP(gdT.subset, reduction = "pca", dims = 1:n_pca_selected)
gdT.subset <- FindNeighbors(gdT.subset, reduction = "pca", dims = 1:n_pca_selected)
gdT.subset <- FindClusters(gdT.subset, resolution = 2)
# gdT.subset$subclusters <- gdT.subset$seurat_clusters
# gdT.subset$clusters <- gdT.subset$seurat_clusters
DimPlot(gdT.subset, label = TRUE, reduction = "umap",repel = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_gdT_umap_cluster.pdf',width = 7, height = 4.34)
DimPlot(gdT.subset, label = TRUE, reduction = "umap",group.by = 'project_ID')
dev.print(pdf, 'tonsil_LAIV_120a_s_gdT_umap_projectID.pdf',width = 7, height = 4.34)
DimPlot(gdT.subset, label = TRUE, reduction = "umap",group.by = 'donor_ID')
dev.print(pdf, 'tonsil_LAIV_120a_s_gdT_umap_donorID.pdf',width = 7, height = 4.34)

DimPlot(gdT.subset, label = TRUE, reduction = "umap",group.by = 'days')
dev.print(pdf, 'tonsil_LAIV_120a_s_gdT_umap_days.pdf',width = 7, height = 4.34)
gdT.cluster.markers <- FindAllMarkers(gdT.subset)
gdT.cluster.markers <- gdT.cluster.markers %>% filter(p_val_adj <= 0.05)
gdT.cluster.markers.top10 <- gdT.cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(gdT.subset, features = c(gdT.cluster.markers.top10$gene)) +
  theme(text = element_text(size=20))
dev.print(pdf, 'tonsil_LAIV_120a_s_gdT_cluster_gene_heatmap.pdf',width = 12, height = 25)
write.csv(gdT.cluster.markers, "tonsil_LAIV_120a_s_gdT_cluster_gene.csv")
write.csv(gdT.cluster.markers.top10, "tonsil_LAIV_120a_s_gdT_cluster_gene_top10.csv")

gdT.subset$AbCD45RA_AbCD45RO_pos <- 'CD45RA+CD45RO- gdT'
gdT.subset$AbCD45RA_AbCD45RO_pos[as.matrix(gdT.subset@assays[["log1p"]]['ab-CD45RO']) > 2/3*(as.matrix(gdT.subset@assays[["log1p"]]['ab-CD45RA']))] <- 'CD45RA-CD45RO+ gdT'
# gdT.subset$AbCD45RA_AbCD45RO_pos[as.matrix(gdT.subset@assays[["log1p"]]['ab-CD45RO']) > (as.matrix(gdT.subset@assays[["log1p"]]['ab-CD45RA']) - 2)] <- 'CD45RA-CD45RO+ gdT'
gdT.subset$AbCD45RA_AbCD45RO_pos[((as.matrix(gdT.subset@assays[["log1p"]]['ab-CD45RO'])) > 3) & ((as.matrix(gdT.subset@assays[["log1p"]]['ab-CD45RA'])) > 4)] <- 'CD45RA+CD45RO+ gdT'

gdT.subset$AbCCR7_AbCD45RA_AbCD45RO_pos <- paste('ab-CCR7-',gdT.subset$AbCD45RA_AbCD45RO_pos,sep = '')
gdT.subset$AbCCR7_AbCD45RA_AbCD45RO_pos[((as.matrix(gdT.subset@assays[["log1p"]]['ab-CCR7'])) > 1.3)] <- paste('ab-CCR7+',gdT.subset$AbCD45RA_AbCD45RO_pos[((as.matrix(gdT.subset@assays[["log1p"]]['ab-CCR7'])) > 1.3)],sep = '')

gdT.subset$TCF7_AbCD45RA_AbCD45RO_pos <- paste('TCF7-',gdT.subset$AbCD45RA_AbCD45RO_pos,sep = '')
gdT.subset$TCF7_AbCD45RA_AbCD45RO_pos[((as.matrix(gdT.subset@assays[["log1p"]]['TCF7'])) > 0.5)] <- paste('TCF7+',gdT.subset$AbCD45RA_AbCD45RO_pos[((as.matrix(gdT.subset@assays[["log1p"]]['TCF7'])) > 0.5)],sep = '')

gdT.subset$GZMH_pos <- 'non-cytotoxic gdT'
gdT.subset$GZMH_pos[((as.matrix(gdT.subset@assays[["log1p"]]['GZMH'])) > 0.5)] <- 'GZMH+ gdT'

gdT.subset$GZMB_PRF1_pos <- 'non-cytotoxic gdT'
gdT.subset$GZMB_PRF1_pos[((as.matrix(gdT.subset@assays[["log1p"]]['GZMB'])) > 0.5) & ((as.matrix(gdT.subset@assays[["log1p"]]['PRF1'])) > 0.5)] <- 'GZMB+PRF1+ gdT'

gdT.subset$MKI67_GZMB_PRF1_pos <- paste('MKI67-',gdT.subset$GZMB_PRF1_pos,sep = '')
gdT.subset$MKI67_GZMB_PRF1_pos[((as.matrix(gdT.subset@assays[["log1p"]]['MKI67'])) > 0.5)] <- paste('MKI67+',gdT.subset$GZMB_PRF1_pos[((as.matrix(gdT.subset@assays[["log1p"]]['MKI67'])) > 0.5)],sep = '')

gdT.subset$TYMS_GZMB_PRF1_pos <- paste('TYMS-',gdT.subset$GZMB_PRF1_pos,sep = '')
gdT.subset$TYMS_GZMB_PRF1_pos[((as.matrix(gdT.subset@assays[["log1p"]]['TYMS'])) > 0.5)] <- paste('TYMS+',gdT.subset$GZMB_PRF1_pos[((as.matrix(gdT.subset@assays[["log1p"]]['TYMS'])) > 0.5)],sep = '')

gdT.subset$ZBTB16_TYMS_GZMB_PRF1_pos <- paste('ZBTB16-',gdT.subset$TYMS_GZMB_PRF1_pos,sep = '')
gdT.subset$ZBTB16_TYMS_GZMB_PRF1_pos[((as.matrix(gdT.subset@assays[["log1p"]]['ZBTB16'])) > 0.5)] <- paste('ZBTB16+',gdT.subset$TYMS_GZMB_PRF1_pos[((as.matrix(gdT.subset@assays[["log1p"]]['ZBTB16'])) > 0.5)],sep = '')

gdT.subset$CCL4L1_pos <- 'CCL4L1- gdT'
gdT.subset$CCL4L1_pos[((as.matrix(gdT.subset@assays[["log1p"]]['CCL4L1'])) > 0.5)] <- 'CCL4L1+ gdT'

gdT.subset$RORC_pos <- 'RORC- gdT'
gdT.subset$RORC_pos[((as.matrix(gdT.subset@assays[["log1p"]]['RORC'])) > 0.5)] <- 'RORC+ gdT'


DimPlot(gdT.subset, label = FALSE, reduction = "umap",group.by = 'RORC_pos')#,order = c("CD27+BCL6+ B","CD27-BCL6+ B","nonPB B")) + 
# scale_colour_manual(values = c("CD27+BCL6+ B" = "red", "nonPB B" = 'grey',"CD27-BCL6+ B" = 'cyan'))
dev.print(pdf, paste('tonsil_LAIV_120a_s_gdT_cleaned_umap_RORC_pos.pdf',sep = ''),width = 600, height = 400)

# DimPlot(gdT.subset, label = FALSE, reduction = "umap",group.by = 'CD184_gdT3_pos')+#,order = c("gdT3+ B","CXCR4+ B")) + 
#   # scale_colour_manual(values = c("CXCR4+ B" = "#F8766D", "CXCR4- B" = '#00BFC4'))
#   scale_colour_manual(values = c("CD184+gdT3+ nonPB B" = "gold", "CD184-gdT3- nonPB B" = "blue","CD184+gdT3- nonPB B" = 'cyan',"CD184-gdT3+ nonPB B" = 'magenta'))
# dev.print(pdf, paste('tonsil_LAIV_120a_s_gdT_cleaned_umap_CD184_gdT3_pos.pdf',sep = ''),width = 600, height = 400)

gdT_log1p_dataframe <- data.frame(t(as.matrix(gdT.subset@assays[["log1p"]]@data)),check.names = F)
gdT_log1p_scale_dataframe <- data.frame(t(as.matrix(gdT.subset@assays[["log1p_scale"]]@data)),check.names = F)

gdT_log1p_scale_dataframe$AbCD45RA_AbCD45RO_pos <- gdT.subset$AbCD45RA_AbCD45RO_pos
ggplot(gdT_log1p_dataframe, aes(x = !!sym(c('ab-CD25')), y = !!sym(c('ab-CXCR5')))) + ggtitle('gdT cells') +
  geom_point(alpha = 0.1,size = 0.5)+#,aes(color = !!sym(c('AbCD45RA_AbCD45RO_pos')))) +
  # geom_point(alpha = 0.1,size = 0.5) +
  # scale_colour_manual(values = c("naive gdT" = "#F8766D", "gdT" = '#00BFC4'))+
  # scale_colour_manual(values = c("CD45RA+gdT3+ nonPB B" = "gold", "CD184-gdT3- nonPB B" = "blue","CD184+gdT3- nonPB B" = 'cyan',"CD184-gdT3+ nonPB B" = 'magenta')) +
  # ylim(0,8) + xlim(0,8) +
  geom_density_2d(color = 'black')
# geom_hline(aes(yintercept=3.25,color = 'red')) + guides(color = 'none') +
# geom_vline(aes(xintercept=4.1,color = 'red')) + guides(color = 'none')
dev.print(pdf, 'tonsil_LAIV_120a_s_gdT_AbCD25_AbCXCR5_contour_log1p_scale.pdf',width = 5, height = 4.6)

grep('CLTC', all.markers, value=TRUE)
gene_name <- 'GNLY'
View(gdT.subcluster.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))
c((gdT.subcluster.markers %>% filter((avg_log2FC > 0) & (gene == gene_name)) %>% arrange(desc(avg_log2FC)))['cluster'])
FeaturePlot(gdT.subset, features = c(gene_name),pt.size = 0.2, ncol = 1, sort.cell = TRUE,min.cutoff = 0, max.cutoff = 3)
dev.print(pdf, paste('tonsil_LAIV_120a_s_gdT_umap_',gene_name,'.pdf',sep = ''),width = 5, height = 4.3)

# ab-CD45RA: 4  8  11 17 15 14
# ab-CD45RO: 12 10 0  9  6  13 3  2 
# ab-CCR7: 8  12 (none visually)
# ab-CD62L: 12 15 0  7 
# ab-CD127: all expressing
# CCR7: 15 4  14
# SELL: 12
# TCF7: 15 17 8  11 14 4 
# IL7R: 18 (most TCF7- expressing)
# LEF1: 15 8  17 4  11


# MKI67:7  29
# GAPDH: 7  29 18 11 15
# CD80LG: 16 11 7 
# ADGRE1: 21 7  3  15 11 4 
# RORC: 7 9
# ab-CD161: 9  10 6  5  23 25 29
# IL17A: 9  7  10 27 6  12 5 
# IL17F: 18 7  9  11 27
# SELPLG: 
# GATA3: 15 29 7  8 
# TBX21: 29 15 16 25 7  26 23 9 
# IFNG: 7  16 9  29 19 8 

# FOXP3: 15 22
# IL10: 15 11 6  16
# ab-CD25: 15 18 7  22 6  29 9 
# IL2RA: 12 16 25 2 

# CXCR5: 0  14 1  2  11
# ab-CXCR5: 14 1  20 24 0  12 3  13
# BCL6: 11 16 0  14 2  21
# CD69: 28 30 8  16
# DUSP1: 8  30 28
# ITGAE: 7  29 23 18 25 21 0 

# ab-CD38: 6  9  18 7  29 21 5  15 10 14 12
# ab-HLA-DR: 14 18 7  9  6  15
# CD38: 25 12
# HLA-DRA: 12 25 2 

# TGFBI: 1  0  13
# TNF: 16 7  11
# EGR3: 16 11
# GZMB: 29 25 23 19 9  26 7 
# GZMK: 19 29 25 26 23
# EOMES: 19 29 25 23 26
# B3GAT1(CD57): 11 19 8  0  16 20 2 (align with the TFH)
# IL21: 16 7  18 11
# PDCD1: 20 24 11 7  2 (align with TFH)
# ab-CD279(PD1): 0  11 6  1  16 2  7  everything except 4,28,3,27,22
# ab-Tim3: 9  29 15 6  14 24 7  18 23 16
# CTLA4: 15 16 8  7  6  18 11
# ARG1 (checkpoint): 22 15 29 11 4  0 not much
# TNFRSF4 (OX40): 11 20 22
# TNFRSF8 (CD30): 7  0  15
# TNFSF8: 11 0  16 7 
# CXCL13: 0  11 18 7  16 2 

# MX1: 7  18 0  29 21 11
# IFIT3: 7  14 29 9  21
# TNFSF13B: 
# GZMA: 
# GZMK: 
# GNLY: 
# NKG7: 
# VEGFA: 2  29 11 15 16
# ALDOC: 2  18 11 14 7  0  15
# TIGIT 
# CD200: 
# CXCL13: 
# EGR1: 
# EGR3: 8
# ADGRG1 (GPR56): 2  6  14 4

### gdT cell cluster fraction visualization #####################
gdT.database <- data.frame(cbind(gdT.subset$donor_ID,as.character(gdT.subset$seurat_clusters),gdT.subset$condition,gdT.subset$days,gdT.subset$treatment,gdT.subset$responder_group,gdT.subset$age_group,gdT.subset$condition_group,gdT.subset$LAIV_response))
colnames(gdT.database) <- c('donor_ID','seurat_clusters','condition','days','treatment','responder_group','age_group','condition_group','LAIV_response')
# sapply(gdT.database,class)
gdT_cluster_condition_donor_count <- count(gdT.database, donor_ID, seurat_clusters,condition, days,treatment,responder_group,age_group,condition_group,LAIV_response)
gdT_cluster_condition_donor_count <- gdT_cluster_condition_donor_count %>% group_by(donor_ID, condition,days,treatment,responder_group,age_group,condition_group,LAIV_response) %>% mutate(total = sum(n))
gdT_cluster_condition_donor_count$percentage <- gdT_cluster_condition_donor_count$n/gdT_cluster_condition_donor_count$total*100
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(Tcell.subset$donor_ID)){
  condition_list <- sort(unique(Tcell.subset$condition[Tcell.subset$donor_ID == donor_name]))
  for (condition_name in condition_list){
    temp_Bcount <- gdT_cluster_condition_donor_count[(gdT_cluster_condition_donor_count$donor_ID == donor_name) & (gdT_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(gdT.subset$seurat_clusters)){
      if_row <- ((gdT_cluster_condition_donor_count$donor_ID == donor_name) & (gdT_cluster_condition_donor_count$seurat_clusters == cluster_name) & (gdT_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Bcount <- gdT_cluster_condition_donor_count[gdT_cluster_condition_donor_count$donor_ID == donor_name,][1,]
        temp_Bcount$seurat_clusters <- cluster_name
        temp_Bcount$n <- 0
        temp_Bcount$percentage <- 0
        temp_Bcount$condition <- condition_name
        temp_Bcount$days <- Tcell.subset$days[Tcell.subset$condition == condition_name][1]
        temp_Bcount$treatment <- Tcell.subset$treatment[Tcell.subset$condition == condition_name][1]
        gdT_cluster_condition_donor_count[nrow(gdT_cluster_condition_donor_count) + 1,] <- temp_Bcount
      }
    }
  }
}
# gdT_cluster_condition_donor_count$responder_group <- factor(gdT_cluster_condition_donor_count$responder_group,levels = c("non-LAIV-responder","LAIV-responder F","LAIV-responder M","LAIV-responder"))

# write.csv(gdT_cluster_condition_donor_count, "tonsil_Rhapsody_gdT_cluster_percentage.csv")
# two figures:
# 1. only responders
# 2. all donors but only at day0
TB_cluster_condition_donor_count_day0 <- gdT_cluster_condition_donor_count[gdT_cluster_condition_donor_count$days == 'day0',]
TB_cluster_condition_donor_count_double_day0 <- gdT_cluster_condition_donor_count
temp <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$days == 'day0',]
TB_cluster_condition_donor_count_double_day0$treatment[TB_cluster_condition_donor_count_double_day0$treatment == 'day0'] <- 'LAIV+'
temp$treatment[temp$treatment == 'day0'] <- 'LAIV-'
TB_cluster_condition_donor_count_double_day0 <- rbind(TB_cluster_condition_donor_count_double_day0,temp)
TB_cluster_condition_donor_count_double_day0 <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$condition_group == 'all timepoints',]
for (cluster_name in unique(gdT.subset$seurat_clusters))
{
  print(cluster_name)
  subdata_allTP <- TB_cluster_condition_donor_count_double_day0[TB_cluster_condition_donor_count_double_day0$seurat_clusters == cluster_name,]
  plot <- ggplot(subdata_allTP,aes(x=days, y=percentage,color = treatment, shape = donor_ID)) +  ggtitle(cluster_name) +
    facet_wrap(~ responder_group + donor_ID, scales = "free_x",ncol = 3) +
    # facet_grid(gender~age_group, scales = "fixed") +
    geom_point(size = 5)+ geom_line(aes(group = interaction(donor_ID,treatment))) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    ylab('% in gdT T cells') +
    # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) + 
    scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) + 
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_double_day0$donor_ID)))
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_gdT_log1p_scale_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 15, height = 9)
  subdata_LAIV <- subdata_allTP[subdata_allTP$treatment == 'LAIV+',] %>% arrange(condition)
  subdata_NS <- subdata_allTP[subdata_allTP$treatment == 'LAIV-',] %>% arrange(condition)
  subdata_LAIV$fold_change <- subdata_LAIV$percentage/subdata_NS$percentage
  graphics.off()
  plot <- ggplot(subdata_LAIV,aes(x=days, y=fold_change,shape = donor_ID)) +  ggtitle(cluster_name) +
    facet_wrap(~ responder_group + donor_ID, scales = "free_x",ncol = 3) +
    # facet_grid(gender~age_group, scales = "fixed") +
    geom_point(size = 5)+ geom_line(aes(group = interaction(donor_ID,treatment))) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    ylab('Fold-change of\n% in gdT T cells in LAIV+ vs LAIV-') + geom_hline(aes(yintercept=1,color = 'red')) + guides(color = 'none') +
    # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) + 
    # scale_colour_manual(values = c("LAIV+" = "#F8766D", "LAIV-" = "#00BFC4")) + 
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_double_day0$donor_ID)))
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_gdT_log1p_scale_',cluster_name,'_condition_FC.pdf',sep = ''),width = 15, height = 9)
  graphics.off()
  ttest_result <- t.test(TB_cluster_condition_donor_count_day0$percentage[(TB_cluster_condition_donor_count_day0$LAIV_response =='non-responder') & (TB_cluster_condition_donor_count_day0$seurat_clusters == cluster_name)],TB_cluster_condition_donor_count_day0$percentage[(TB_cluster_condition_donor_count_day0$LAIV_response =='responder') & (TB_cluster_condition_donor_count_day0$seurat_clusters == cluster_name)])
  plot <- ggplot(TB_cluster_condition_donor_count_day0[TB_cluster_condition_donor_count_day0$seurat_clusters == cluster_name,],aes(x=LAIV_response, y=percentage,color = LAIV_response, shape = donor_ID)) +  ggtitle(paste('day0 pval:',as.integer(10*ttest_result$p.value)/10,'\n',cluster_name)) +
    facet_grid(~age_group, scales = "free_x") +
    geom_point(size = 5) + theme(text = element_text(size = 28),plot.title = element_text(size = 30, face = "bold")) + RotatedAxis()+
    scale_colour_manual(values = c("non-responder" = "#00BFC4", "responder" = "#F8766D")) +
    scale_shape_manual(values=1:length(unique(TB_cluster_condition_donor_count_day0$donor_ID))) +
    ylab('% in gdT T cells')
  # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_gdT_log1p_scale_day0_',cluster_name,'_condition_ratio.pdf',sep = ''),width = 11, height = 7)
  
}


naive <- c('ab-CD45RO','ab-CD45RA','ab-CCR7','ab-CD62L','CCR7','SELL','LEF1','TCF7','PASK','IL7R','CD27','ab-CD27','ab-CD38','ab-HLA-DR','HLA-DRA','CD38','CD44','CD28','CD40LG')
cytotoxic <- c('NKG7','GNLY','GZMK','GZMA','GZMB','KLRC3','KLRC4','IFNG','CCL5','CCL3','CCL4','IL2RB')
NK_like <- c('FCGR3A','ADGRG1','ADGRG3','ab-CD56','NCAM1','B3GAT1')
inhibitory <- c('PDCD1','ab-CD279','CTLA4','ab-Tim3','HAVCR2','TNFRSF4','ICOS','LAG3','TNFRSF18','GAPDH','ab-CD25','IL2RA','FOXP3')
# inhibitory <- c('ab-KIR-NKAT2','ab-CD158e1','ab-Tim3','HAVCR2','PDCD1','ab-CD279','CTLA4','TNFRSF4','ICOS','LAG3','TNFRSF18','ab-CD25','IL2RA','FOXP3')
# LAMP3_set <- c('LAMP3','TNFSF13B','CD274','STAT1', 'LGALS9', 'OAS1', 'OASL', 'GIMAP5', 'ITK', 'HLA-A', 'LAP3')
IFN_stim <- c('MX1','IFI6','IFITM3','IFI44L','XAF1','IFIT3','ISG15','IFIT1','STAT1','LAP3','OAS1','TNFSF10')
# IFN_stim <- c('MX1','IFI6','IFITM3','IFI44L','XAF1','IFIT3','ISG15','IFIT1','LAMP3','TNFSF13B','RSAD2','STAT1', 'LGALS9', 'OAS1', 'OASL', 'FAS','HLA-A', 'LAP3','GIMAP5', 'ITK','TNFSF10')
RORC_set <- c('RORC','IL17A','IL17F','ab-CD161','KLRB1','CXCR6')#,'CCR2','CCR5','KLRB1','LGALS3','ab-CD39','ENTPD1','DPP4','NCR3','CCR6','SLAMF1','CD58','ITGA4','LGALS1','S100A4','CD40LG','ARL4C','IL32')
activation <- c('ITGAE','CD69','FOS','DUSP1','FOSB','JUN','JUNB')
VEGF_set <- c('VEGFA','ALDOC', 'BNIP3L','SLC2A3')
CD200_set <- c('EGR3','TNF','ERG1','CXCL13','CD200','TNFSF8','BTLA','TIGIT')
# tisssue_resident <- c('ITGA1','ITGA1')
# functional <- c('MKI67','ab-CD71','TBX21','GATA3','EOMES','CXCR5','ab-CXCR5','BCL6','CXCL13','ADGRE1','CXCR3','FASLG','TXK','KLRG1','ICAM1','NRP1','TFRC','IRF4')
functional <- c('MKI67','ab-CD71','TBX21','GATA3','EOMES','CXCR5','ab-CXCR5','BCL6','ADGRE1','CXCR3','FASLG','TXK','KLRG1','ICAM1','NRP1','TFRC','IRF4','TGFBI')
IL_cytokines <- c('IL21','IL23A','IL26','IL10','IL6','IL12RB2')#IL32
# gene_features <- list('lineage_marker' = lineage_marker,'lineage_gene' = lineage_gene,'naive / memory' = naive, 'cytotoxic'= cytotoxic, 'NK_like' = NK_like,"activation / exhaustion / inhibitory" = inhibitory,'IFN_stim' = IFN_stim,'transient activation' = activation,'functional' = functional,'IL_cytokines' = IL_cytokines)
gene_features <- list('naive / memory' = naive, 'cytotoxic'= cytotoxic, 'NK_like' = NK_like,"activation / exhaustion / inhibitory" = inhibitory,'EGR3 set' = CD200_set,'IFN_stim' = IFN_stim,'transient activation' = activation,'RORC_set' = RORC_set,'VEGF_set' = VEGF_set,'functional' = functional,'IL_cytokines' = IL_cytokines)
# gdT.subset$subclusters <- as.character(Tcell.subset$seurat_clusters)
# gdT.subset$subclusters <- factor(gdT.subset$subclusters,levels = c(25,29,23,19,26,17,30,8,28,4,3,22,27,12,20,13,1,2,0,21,11,16,7,18,15,6,9,5,14,10,24))
DotPlot(gdT.subset,features = gene_features, dot.scale = 5,group.by = 'clusters') + RotatedAxis() + ggtitle('Tonsil Organoid 120a_s gdT cells') +
  theme(axis.text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_gdT_cluster_dotplot.pdf',sep = ''),width = 25, height = 600)

FeaturePlot(Tcell.subset, features = c("CXCR5", "CXCR5.1"),pt.size = 1, ncol = 2, sort.cell = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_CXCR5.pdf',width = 12, height = 4.34)

FeaturePlot(Tcell.subset, features = c("BCL6"),pt.size = 1, ncol = 1, sort.cell = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_BCL6.pdf',width = 600, height = 4.34)

FeaturePlot(Tcell.subset, features = c("CTLA4", "ICOS"),pt.size = 1, ncol = 2, sort.cell = TRUE)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_CTLA4_ICOS.pdf',width = 12, height = 4.34)

## There is RORC+ Tcell cells, but they are not standing out as a cluster due to a lack of resolution
Tcell.subset$cluster_manual_label <- 'T cells'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(4)] <- 'naive CD4'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(25)] <- 'gdT'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(17)] <- 'naive CD8'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(30)] <- 'transient activated CD8'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(8,28)] <- 'transient activated CD4'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(29)] <- 'proliferating cytotoxic CD8'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(7)] <- 'proliferating activated CM' # half BCL6+
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(23,19,26)] <- 'cytotoxic CD8'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(5,9,10,27)] <- 'IL17+EM or Th17'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(15,22)] <- 'Treg'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(18,6)] <- 'activated EM'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(2,11,16,0,21,14)] <- 'TFH'
Tcell.subset$cluster_manual_label[Tcell.subset$seurat_clusters %in% c(1,20,24,12,3,13)] <- 'CXCR5+EM'
DimPlot(Tcell.subset, label = TRUE, reduction = "umap",repel = TRUE,group.by = 'cluster_manual_label') + ggtitle('tonsil orgnaoid T cells') + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_cluster_labeled.pdf',width = 7, height = 4.34)

ordered_Tcell_list <- c('cytotoxic gdT','proliferating cytotoxic CD8','cytotoxic KIR+CD8','cytotoxic CD8','transient activated CD8','naive CD8','naive CD4','non-TFH-activated EM','CXCR5+EM','TFH','transient activated CD4','activated EM','proliferating activated CM','activated Treg','Treg')
Tcell.subset$cluster_manual_label <- factor(Tcell.subset$cluster_manual_label,levels = ordered_Tcell_list)

lineage_marker <- c('ab-CD3','ab-CD20','ab-CD19','ab-CD4','ab-CD8','ab-TCR-gamma-delta','ab-TCR-alpha-beta')
lineage_gene <- c('CD3D','CD3E','CD3G','MS4A1','CD4','CD8A','CD8B','TRAC','TRBC2','TRDC')
naive <- c('ab-CD45RO','ab-CD45RA','ab-CCR7','CCR7','ab-CD62L','SELL','LEF1','TCF7','PASK','IL7R','CD27','ab-CD27','CD38','ab-CD38','HLA-DRA','ab-HLA-DR','CD44','CD28','CD40LG')
cytotoxic <- c('NKG7','GNLY','GZMK','GZMA','GZMB','KLRC3','KLRC4','IFNG','TNF','CCL5','CCL3','CCL4','IL2RB')
NK_like <- c('FCGR3A','ADGRG1','ADGRG3','ab-CD56','NCAM1','B3GAT1')
inhibitory <- c('PDCD1','ab-CD279','CTLA4','ab-Tim3','HAVCR2','TNFRSF4','ICOS','LAG3','TNFRSF18','GAPDH','ab-CD25','IL2RA','FOXP3')
IFN_stim <- c('MX1','IFI6','IFITM3','IFI44L','XAF1','IFIT3','ISG15','IFIT1','STAT1','LAP3','OAS1','TNFSF10')
RORC_set <- c('RORC','IL17A','IL17F','ab-CD161','KLRB1','CXCR6')#,'CCR2','CCR5','KLRB1','LGALS3','ab-CD39','ENTPD1','DPP4','NCR3','CCR6','SLAMF1','CD58','ITGA4','LGALS1','S100A4','CD40LG','ARL4C','IL32')
activation_edited <- c('ITGAE','CD69','EGR3','CXCL13','CD200','ERG1','TNFSF8','BTLA','TIGIT','FOS','DUSP1','FOSB','JUN','JUNB')
VEGF_set <- c('VEGFA','ALDOC', 'BNIP3L','SLC2A3')
# tisssue_resident <- c('ITGA1','ITGA1')
# functional <- c('MKI67','ab-CD71','TBX21','GATA3','EOMES','CXCR5','ab-CXCR5','BCL6','CXCL13','ADGRE1','CXCR3','FASLG','TXK','KLRG1','ICAM1','NRP1','TFRC','IRF4')
functional <- c('MKI67','ab-CD71','TBX21','GATA3','EOMES','CXCR5','ab-CXCR5','BCL6','ADGRE1','CXCR3','FASLG','TXK','KLRG1','ICAM1','NRP1','TFRC','IRF4','TGFBI')
IL_cytokines <- c('IL21','IL23A','IL26','IL10','IL6','IL12RB2')#IL32
# gene_features <- list('lineage_marker' = lineage_marker,'lineage_gene' = lineage_gene,'naive / memory' = naive, 'cytotoxic'= cytotoxic, 'NK_like' = NK_like,"activation / exhaustion / inhibitory" = inhibitory,'IFN_stim' = IFN_stim,'transient activation' = activation,'functional' = functional,'IL_cytokines' = IL_cytokines)
gene_features_edited <- list('naive / memory' = naive, 'cytotoxic'= cytotoxic, 'NK_like' = NK_like,"activation / exhaustion / inhibitory" = inhibitory,'IFN_stim' = IFN_stim,'transient activation' = activation_edited,'RORC_set' = RORC_set,'VEGF_set' = VEGF_set,'functional' = functional,'IL_cytokines' = IL_cytokines)

DotPlot(Tcell.subset,features = gene_features_edited, dot.scale = 8,group.by = 'subcluster_manual_label') + RotatedAxis() + ggtitle('Tonsil Organoid treatment T cells') +
  theme(axis.text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_subcluster_dotplot_labeled.pdf',sep = ''),width = 3000, height = 800)

Tcell.subset$subcluster_manual_label <- 'T cells'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(4)] <- 'naive CD4'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(25)] <- 'gdT'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(17)] <- 'naive CD8'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(14)] <- 'activated TFH SCM'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(30)] <- 'transient activated CD8'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(8,28)] <- 'transient activated CD4'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(29)] <- 'proliferating cytotoxic CD8'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(7)] <- 'proliferating activated CM' # half BCL6+
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(23,19)] <- 'cytotoxic CD8'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(26)] <- 'KIR+cytotoxic CD8'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(5)] <- 'Th17'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(9)] <- 'activated Th17'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(10,27)] <- 'IL17+CD4'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(15)] <- 'activated Treg'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(22)] <- 'Treg'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(18,6)] <- 'activated EM'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(11,16)] <- 'EGR3+TFH'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(2)] <- 'VEGFA+TFH'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(0,21)] <- 'TFH'
Tcell.subset$subcluster_manual_label[Tcell.subset$seurat_clusters %in% c(1,20,24,12,3,13)] <- 'CXCR5+EM'

DimPlot(Tcell.subset, label = TRUE, reduction = "umap",repel = TRUE,group.by = 'subcluster_manual_label') + ggtitle('tonsil orgnaoid T cells') + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_subcluster_labeled.pdf',width = 7, height = 4.34)

# Tcell.subset$subcluster_manual_label <- factor(Tcell.subset$subcluster_manual_label,levels = c('cytotoxic gdT', 'cytotoxic CD8','naive CD8','naive CD4','EM','CXCR5+EM','TFH','transient activated TFH','activated EM','proliferating activated CM','Th17','Treg'))
DotPlot(Tcell.subset,features = gene_features, dot.scale = 5,group.by = 'subcluster_manual_label') + RotatedAxis() + ggtitle('Tonsil Organoid treatment T cells') +
  theme(axis.text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_subcluster_dotplot_labeled.pdf',sep = ''),width = 3000, height = 600)

# 
# Tcell.subset$treatment <- 'day0'
# Tcell.subset$treatment[grepl('I53',Tcell.subset$condition)] <- 'I53'
# Tcell.subset$treatment[grepl('ns',Tcell.subset$condition)] <- 'ns'

DimPlot(Tcell.subset, label = TRUE, reduction = "umap",repel = TRUE,group.by = 'cluster_manual_label') + ggtitle('tonsil orgnaoid T cells') + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_umap_cluster_labeled.pdf',width = 7, height = 4.34)


Idents(Tcell.subset) <- Tcell.subset$cluster_manual_label
Tcell_cluster.markers <- FindAllMarkers(Tcell.subset,base = exp(1))
Tcell_cluster.markers <- Tcell_cluster.markers %>% filter(p_val_adj <= 0.05)
Tcell_cluster.markers.top10 <- Tcell_cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Tcell.subset, features = c(Tcell_cluster.markers.top10$gene)) +
  theme(text = element_text(size=20))
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_cluster_labeled_gene_heatmap.pdf',width = 12, height = 25)
write.csv(Tcell_cluster.markers, "tonsil_LAIV_120a_s_Tcell_cluster_labeled_gene.csv")
write.csv(Tcell_cluster.markers.top10, "tonsil_LAIV_120a_s_Tcell_cluster_labeld_gene_top10.csv")

Tcell.subset$cluster_treatment <- paste(Tcell.subset$cluster_manual_label,Tcell.subset$treatment)
Tcell.subset$donor_condition <- paste(Tcell.subset$donor_ID,Tcell.subset$condition)
Tcell.subset$meta_cluster_treatment <- paste(Tcell.subset$meta_cluster,Tcell.subset$treatment)
Tcell.subset$meta_cluster_condition <- paste(Tcell.subset$meta_cluster,Tcell.subset$condition)
Tcell.subset$cluster_donor <- paste(as.character(Tcell.subset$cluster_manual_label),Tcell.subset$donor_ID)

DotPlot(object = Tcell.subset,group.by = 'cluster_donor',cols = c( "blue","red"),dot.scale = 10,features=gene_features) + RotatedAxis() + ggtitle('tonsil response age comparison Tcell clusters') + theme(axis.text = element_text(size = 20),plot.title = element_text(size = 30, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_cluster_dotplot_labeled_donor.pdf',sep = ''),width = 3200, height = 12)


VlnPlot(Tcell.subset, features = c("VEGFA"),pt.size = 0.2,group.by = 'subVEGFA_pos',sort = 'increasing',log = F) + NoLegend()# + geom_hline(yintercept = -0)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_violin_VEGFA.pdf',width = 12, height = 4.34)
Tcell.subset$TBX21_pos <- 'TBX21- Tcells'
Tcell.subset$TBX21_pos[(as.matrix(Tcell.subset@assays[["RNA"]]['TBX21']) > 0)] <- 'TBX21+ Tcells'
# Tcell.subset$TBX21_pos[(as.matrix(Tcell.subset@assays[["RNA"]]['HLA-DRA']) > -0.7) & (as.matrix(Tcell.subset@assays[["RNA"]]['ab-CD38']) > 0)] <- 'activated Tcells'
DimPlot(Tcell.subset, label = FALSE, reduction = "umap",group.by = 'TBX21_pos')

### day0 pie chart 
colnames(Tcell_subcluster_condition_donor_count) <- c('donor','Tcell_type','condition','n','total','cell.fraction','stim')
ordered_Tcell_list <- unique(Tcell.subset$meta_cluster)#c('gdT','naive CD8','RSAD2+ IFN-stim CD8','transient activated CD8','effector memory CD8','RORC+ effector memory CD8','activated effector memory CD8','proliferating activated effector memory CD8','NK-like CD8','PD1++TFH','Treg','proliferating actviated effector memory CD4','actviated Th17','actviated effector memory CD4','VEGF+TFH','CXCR5+CD4','central memory TFH','RSAD2+ IFN-induced central memory CD4','Th2','central memroy CD4','transient activated CD4','naive CD4')#
cell_colors <- c('gray85', 'red', 'orange', 'blue', 'limegreen', 'skyblue', '#88fcd1', '#ee00a4', 'purple', 'black', 'pink', 'gold', 'firebrick', 'green', 'slateblue','yellow','cyan','wheat4','violetred2','turquoise4','tomato3','thistle3','brown2')
cell_colors <- cell_colors[1:length(ordered_Tcell_list)]
condition_list <- unique(Tcell.subset$condition)
Tcell_subcluster_condition_donor_count <- Tcell_subcluster_condition_donor_count %>%
  arrange(desc(Tcell_type)) %>%
  mutate(cell.fraction_label = replace(Tcell_type,cell.fraction == 0, ''))
Tcell_subcluster_condition_donor_count <- Tcell_subcluster_condition_donor_count %>%
  arrange(desc(Tcell_type)) %>%
  mutate(label.ypos = cumsum(cell.fraction)- 0.5*cell.fraction )

for (condition in condition_list)
{
  print(condition)
  plot <- ggplot(Tcell_subcluster_condition_donor_count[Tcell_subcluster_condition_donor_count$condition == condition,], aes(x = "", y = cell.fraction, fill = Tcell_type)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    ggtitle(paste('tonsil organoid T cells ',condition)) +
    coord_polar("y")+
    facet_grid(~donor) + 
    guides(fill=guide_legend(ncol=1)) + 
    # geom_text_repel(aes(x = 1.4, y = label.ypos, label = cell.fraction_label), color = "black", size = 5, nudge_x = 0.3, segment.size = 0.7, show.legend = FALSE)+
    scale_fill_manual(values = cell_colors) +
    theme_void()# +
  # theme(plot.title = element_text(size = 20))# + NoLegend()
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_',condition,'_frequency_condition_lineage.pdf',sep = ''),width = 7, height = 5)
}

Tcell.subset$MX1_pos <- 'MX1-'
Tcell.subset$MX1_pos[(as.matrix(Tcell.subset@assays[["integrated_scale"]]['MX1']) > 0) & (Tcell.subset$TB_division == 'Tcells')] <- 'MX1+ T cells'
DimPlot(Tcell.subset, label = FALSE, reduction = "umap",group.by = 'MX1_pos')

# single volcano plot:
logFC_cutoff <- log(1.5)
p_cutoff <- 0.05
cluster_subset_markers_subset <- FindMarkers(Tcell.subset, group.by = "MX1_pos",ident.1 = 'MX1+ T cells',ident.2 = NULL,logfc.threshold = 0)
cluster_subset_markers_subset$gene <- rownames(cluster_subset_markers_subset)
cluster_subset_markers_subset$significant <- 'no'
cluster_subset_markers_subset$significant[(cluster_subset_markers_subset$p_val_adj <= p_cutoff) & (abs(cluster_subset_markers_subset$avg_logFC) >= logFC_cutoff)] <- 'yes'
cluster_subset_markers_subset$labels <- cluster_subset_markers_subset$gene
cluster_subset_markers_subset$labels[cluster_subset_markers_subset$significant == 'no'] <- ''
write.xlsx(cluster_subset_markers_subset,paste('tonsil_LAIV_120a_s_Tcell_MX1+_vs_rest.xlsx'))
graphics.off()
cluster_subset_markers_subset$log_p_val_adj <- -log10(cluster_subset_markers_subset$p_val_adj)
plot<-ggplot(cluster_subset_markers_subset, aes(x=avg_logFC, y=log_p_val_adj,color = significant,legend = significant)) +
  geom_point() +
  scale_colour_manual(values = c("yes" = "red", "no" = "darkgrey")) +
  ggtitle(paste('DE of MX1+ vs rest T cells')) +
  xlab("avg_logFC") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 20)) +
  geom_text_repel(aes(x=avg_logFC, y=log_p_val_adj,label = labels), color = "black", size = 5) +
  geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_MX1+_vs_rest.pdf',sep = ''),width = 600, height = 4.34*1.2)



##### DE analysis ###################################
# Tcell.subset$donor_ID <- Tcell.subset$project_ID
# Tcell.subset$donor_ID[Tcell.subset$project_ID == '120e'] <- '03yrs F IMD120'
# Tcell.subset$donor_ID[(Tcell.subset$project_ID == '120f') | (Tcell.subset$project_ID == '120g')] <- '08yrs M IMD164'
# Tcell.subset$donor_ID[(Tcell.subset$project_ID == '120h') | (Tcell.subset$project_ID == '120i')] <- '35yrs M BBENT001'
# 
# Tcell.subset$condition <- gsub('I53','I53+',Tcell.subset$condition)
# Tcell.subset$condition <- gsub('ns','I53-',Tcell.subset$condition)
# Tcell.subset$treatment <- gsub('I53','I53+',Tcell.subset$treatment)
# Tcell.subset$treatment <- gsub('ns','I53-',Tcell.subset$treatment)
# Tcell.subset$days <- substring(Tcell.subset$condition,first = 1, last = 5)
CD4.subset <- subset(Tcell.subset,meta_cluster == 'CD4')
gdT.subset <- subset(Tcell.subset,meta_cluster == 'gdT')
CD8.subset <- subset(Tcell.subset,meta_cluster == 'CD8')
# volcano plot
library(ggrepel)
logFC_cutoff <- 0.25
p_cutoff <- 0.05
#Seurat findmarker volcano
DefaultAssay(CD4.subset) <- "log1p"


library(ggrepel)
logFC_cutoff <- 0.25
p_cutoff <- 0.05
#Seurat findmarker volcano
DefaultAssay(CD4.subset) <- "log1p"

for (row_index in 1:dim(condition_nonpaired_table)[1]){
  donor_name <- condition_nonpaired_table$donor_ID[row_index]
  day_name <- condition_nonpaired_table$day_group[row_index]
  CD4_LAIV_ns.markers <- FindMarkers(CD4.subset, group.by = "condition",ident.1 = condition_nonpaired_table$LAIV[row_index],ident.2 = condition_nonpaired_table$ns[row_index],logfc.threshold = 0,min.pct = 0)
  CD4_LAIV_ns.markers$FoldChange <- exp(CD4_LAIV_ns.markers$avg_log2FC)
  CD4_LAIV_ns.markers$pct_diff <- CD4_LAIV_ns.markers$pct.1 - CD4_LAIV_ns.markers$pct.2
  CD4_LAIV_ns.markers$pct_FC <- CD4_LAIV_ns.markers$pct.1/CD4_LAIV_ns.markers$pct.2
  CD4_LAIV_ns.markers$significant <- 'n.s.'
  CD4_LAIV_ns.markers$significant[(abs(CD4_LAIV_ns.markers$avg_log2FC) <= logFC_cutoff) & (CD4_LAIV_ns.markers$p_val_adj <= p_cutoff)] <- 'significant'
  CD4_LAIV_ns.markers$significant[(abs(CD4_LAIV_ns.markers$avg_log2FC) > logFC_cutoff) & (CD4_LAIV_ns.markers$p_val_adj <= p_cutoff)] <- 'significant\nand distinct'
  CD4_LAIV_ns.markers$gene <- rownames(CD4_LAIV_ns.markers)
  CD4_LAIV_ns.markers$sign <- 'upregulated'
  CD4_LAIV_ns.markers$sign[CD4_LAIV_ns.markers$avg_log2FC < 0] <- 'downregulated'
  CD4_LAIV_ns.markers$labels <- rownames(CD4_LAIV_ns.markers)
  CD4_LAIV_ns.markers$labels[CD4_LAIV_ns.markers$significant == 'no'] <- ''
  CD4_LAIV_ns.markers <- CD4_LAIV_ns.markers %>% arrange(desc(significant),desc(sign),p_val_adj)
  p <- match(CD4_LAIV_ns.markers$gene, all_marker_table$marker)
  temp <- all_marker_table[p,]
  CD4_LAIV_ns.markers$ENTREZID <- temp$ENTREZID
  write.xlsx(CD4_LAIV_ns.markers,paste('tonsil_LAIV_120a_s_CD4_LAIV_vs_ns_',donor_name,'_',day_name,'_log1p.xlsx',sep = ''))        
  graphics.off()
  plot<-ggplot(CD4_LAIV_ns.markers, aes(x=avg_log2FC, y=-log10(p_val_adj),color = significant,legend = significant)) +
    geom_point() +
    scale_colour_manual(values = c("significant\nand distinct" = "red", "significant" = 'orange',"n.s." = "darkgrey")) +
    ggtitle(paste('DE of CD4 with LAIV+ vs LAIV-\n in',donor_name,day_name)) +
    xlab("log2 of Fold Change") + ylab("-log10 adjusted p-value") +
    theme(text = element_text(size = 20)) +
    geom_text_repel(aes(x=avg_log2FC, y=-log10(p_val_adj),label = labels), color = "black", size = 5) +
    geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff), linetype="dotted") +
    geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_LAIV_vs_ns_',donor_name,'_',day_name,'_log1p.pdf',sep = ''),width = 600, height = 4.34*1.2)
}

# DefaultAssay(CD4.subset) <- "log1p"
# for (day_name in unique(CD4.subset$days)[3:9]){
#   print(day_name)
#   if (sum(CD4.subset$days == day_name) > 6){
#     selected_subset  <- subset(CD4.subset, days == day_name)
#     for (donor_name in unique(selected_subset$donor_ID))
#     {
#       if (sum(selected_subset$donor_ID == donor_name) > 6){
#         cluster_subset <- subset(selected_subset, donor_ID == donor_name)
#         print(donor_name)
#         if ((sum(cluster_subset$treatment == 'LAIV+') > 3) & (sum(cluster_subset$treatment == 'LAIV-') > 3)){
#           CD4_LAIV_ns.markers <- FindMarkers(cluster_subset, group.by = "treatment",ident.1 = 'LAIV+',ident.2 = 'LAIV-',logfc.threshold = 0,min.pct = 0)
#           CD4_LAIV_ns.markers$FoldChange <- exp(CD4_LAIV_ns.markers$avg_log2FC)
#           CD4_LAIV_ns.markers$pct_diff <- CD4_LAIV_ns.markers$pct.1 - CD4_LAIV_ns.markers$pct.2
#           CD4_LAIV_ns.markers$pct_FC <- CD4_LAIV_ns.markers$pct.1/CD4_LAIV_ns.markers$pct.2
#           CD4_LAIV_ns.markers$significant <- 'no'
#           CD4_LAIV_ns.markers$significant[(abs(CD4_LAIV_ns.markers$avg_log2FC) > logFC_cutoff) & (CD4_LAIV_ns.markers$p_val_adj <= p_cutoff)] <- 'yes'
#           # CD4_LAIV_ns.markers$significant[(((CD4_LAIV_ns.markers$avg_log2FC >= logFC_cutoff) & (CD4_LAIV_ns.markers$pct_diff >= 0.05)) | ((CD4_LAIV_ns.markers$avg_log2FC <= -logFC_cutoff) & (CD4_LAIV_ns.markers$pct_diff <= -0.05))) & (CD4_LAIV_ns.markers$p_val_adj <= p_cutoff)] <- 'yes'
#           CD4_LAIV_ns.markers$gene <- rownames(CD4_LAIV_ns.markers)
#           CD4_LAIV_ns.markers$sign <- 'upregulated'
#           CD4_LAIV_ns.markers$sign[CD4_LAIV_ns.markers$avg_log2FC < 0] <- 'downregulated'
#           CD4_LAIV_ns.markers$labels <- rownames(CD4_LAIV_ns.markers)
#           CD4_LAIV_ns.markers$labels[CD4_LAIV_ns.markers$significant == 'no'] <- ''
#           CD4_LAIV_ns.markers <- CD4_LAIV_ns.markers %>% arrange(desc(significant),desc(sign),p_val_adj)
#           p <- match(CD4_LAIV_ns.markers$gene, all_marker_table$marker)
#           temp <- all_marker_table[p,]
#           CD4_LAIV_ns.markers$ENTREZID <- temp$ENTREZID
#           write.xlsx(CD4_LAIV_ns.markers,paste('tonsil_LAIV_120a_s_CD4_LAIV_vs_ns_',donor_name,'_',day_name,'_log1p.xlsx',sep = ''))
#            graphics.off()
#           plot<-ggplot(CD4_LAIV_ns.markers, aes(x=avg_log2FC, y=-log10(p_val_adj),color = significant,legend = significant)) +
#             geom_point() +
#             scale_colour_manual(values = c("yes" = "red", "no" = "darkgrey")) +
#             ggtitle(paste('DE of CD4 with LAIV+ vs LAIV-\n in',donor_name,day_name)) +
#             xlab("log2 of Fold Change") + ylab("-log10 adjusted p-value") +
#             theme(text = element_text(size = 20)) +
#             geom_text_repel(aes(x=avg_log2FC, y=-log10(p_val_adj),label = labels), color = "black", size = 5) +
#             geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff), linetype="dotted") +
#             geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")
#           print(plot)
#           dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_LAIV_vs_ns_',donor_name,'_',day_name,'_log1p.pdf',sep = ''),width = 600, height = 4.34*1.2)
#         }
#       }
#     }
#   }
# }
# 


DefaultAssay(CD4.subset) <- "log1p"

# MAIT.subset <- subset(CD4.subset,differentiation == 'CD161+CD4')
Tfh.subset <- subset(CD4.subset,BCL6_AbCD161_pos == 'Tfh')

CD4_LAIV_ns.markers <- FindMarkers(Tfh.subset, group.by = "clonal_expansion",ident.1 = 'yes',ident.2 = 'no',logfc.threshold = 0,min.pct = 0)
CD4_LAIV_ns.markers$FoldChange <- exp(CD4_LAIV_ns.markers$avg_log2FC)
CD4_LAIV_ns.markers$pct_diff <- CD4_LAIV_ns.markers$pct.1 - CD4_LAIV_ns.markers$pct.2
CD4_LAIV_ns.markers$pct_FC <- CD4_LAIV_ns.markers$pct.1/CD4_LAIV_ns.markers$pct.2
CD4_LAIV_ns.markers$significant <- 'no'
CD4_LAIV_ns.markers$significant[(abs(CD4_LAIV_ns.markers$avg_log2FC) > logFC_cutoff) & (CD4_LAIV_ns.markers$p_val_adj <= p_cutoff)] <- 'yes'
CD4_LAIV_ns.markers$gene <- rownames(CD4_LAIV_ns.markers)
CD4_LAIV_ns.markers$sign <- 'upregulated'
CD4_LAIV_ns.markers$sign[CD4_LAIV_ns.markers$avg_log2FC < 0] <- 'downregulated'
CD4_LAIV_ns.markers$labels <- rownames(CD4_LAIV_ns.markers)
CD4_LAIV_ns.markers$labels[CD4_LAIV_ns.markers$significant == 'no'] <- ''
CD4_LAIV_ns.markers <- CD4_LAIV_ns.markers %>% arrange(desc(significant),desc(sign),p_val_adj)
write.xlsx(CD4_LAIV_ns.markers,paste('tonsil_LAIV_120a_s_CD4_BCL6_AbCD161_pos_Tfh_clonal_expansion_DE.xlsx',sep = ''))
plot<-ggplot(CD4_LAIV_ns.markers, aes(x=avg_log2FC, y=-log10(p_val_adj),color = significant,legend = significant)) +
  geom_point() +
  scale_colour_manual(values = c("yes" = "red", "no" = "darkgrey")) +
  ggtitle(paste('DE of clonal expanded CD161-Tfh vs singlets')) +
  xlab("log2 of Fold Change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 20)) +
  geom_text_repel(aes(x=avg_log2FC, y=-log10(p_val_adj),label = labels), color = "black", size = 5) +
  geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_BCL6_AbCD161_pos_Tfh_clonal_expansion_DE_log1p.pdf',sep = ''),width = 600, height = 4.34)
# dev.print(pdf, paste('tonsil_LAIV_120a_s_CD4_CD4_clonal_expansion_DE_log1p.pdf',sep = ''),width = 5, height = 4.34*1.5)
