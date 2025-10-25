library(ComplexHeatmap)
library(dendextend)
library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(ggrepel)
library("readxl")
library('stringr')

batch_marker_list <- list()
batch_marker_list[[1]] <- readRDS('tonsil_LAIV_120a.rds')
batch_marker_list[[2]] <- readRDS('tonsil_LAIV_120b.rds')
batch_marker_list[[3]] <- readRDS('tonsil_LAIV_120c.rds')
batch_marker_list[[4]] <- readRDS('tonsil_LAIV_120j.rds')
minimal_marker_list <- Reduce(intersect,batch_marker_list[2:4])

Tcell_LAIV.subset <- subset(Tcell_LAIV.subset,donor_ID != '02yrs M IMD085')

Tcell_dataframe_lognorm <- data.frame(t(as.matrix(Tcell_LAIV.subset@assays[["normalized"]]@data)),check.names = F)
Tcell_dataframe_lognorm <- Tcell_dataframe_lognorm[,colSums(Tcell_dataframe_lognorm) != 0]
cor_mat <- cor(Tcell_dataframe_lognorm,use = 'complete.obs')#[selected_table_day0_day04_day07$age_group == 'Adults',]

HM2 <- Heatmap(cor_mat, cluster_columns = T, cluster_rows = T,
               name = "mean\nnormalized\nexpression", show_column_names = T,
               show_row_names = F,#row_km = 20, column_km = 20,
               column_title_side = "top",row_title_side = "left",row_names_side = 'left')
HM2 <- draw(HM2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_LAIV_gene_module_cor.pdf',width = 80, height = 80)

r.dend <- row_dend(HM2)  #Extract row dendrogram
cor_cluster_table <- cutree(r.dend, h = 1.5)
sort(table(cor_cluster_table))
# library(factoextra)
# library(cluster)
# fviz_nbclust(cor_mat, k.max = 40,FUN = hcut, method = "wss")
# dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_gene_module_cor_hcluster_elbow.pdf',width = 10, height = 5)
# fviz_nbclust(cor_mat, k.max = 40,FUN = hcut, method = "silhouette")
# dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_gene_module_corhcluster_silhouette.pdf',width = 10, height = 5)
# gap_stat <- clusGap(cor_mat, FUN = hcut, nstart = 25, K.max = 40, B = 50)
# fviz_gap_stat(gap_stat)
# dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_gene_module_corhcluster_GasStat.pdf',width = 10, height = 5)

# k_cluster_table <- data.frame(k = c(3:200),n_count5 = NA, n_MX1 = NA,n_VEGFA = NA,n_AURKB = NA,n_FOSB = NA,n_BCL6 = NA)
# for (k_val in k_cluster_table$k) {
k_val <- 43
cor_cluster_table <- data.frame(cutree(r.dend, k = k_val))
colnames(cor_cluster_table) <- 'cluster'
cor_cluster_table$marker <- rownames(cor_cluster_table)
cor_cluster_table$k <- k_val
cluster_count_table <- data.frame(rev(sort(table(cor_cluster_table$cluster))))
colnames(cluster_count_table) <- c('cluster','Freq')
k_cluster_table$n_count5[k_cluster_table$k == k_val] <- sum(cluster_count_table$Freq>=5)

cor_cluster_table <- merge(cluster_count_table,cor_cluster_table,by = 'cluster')
# k_cluster_table$n_MX1[k_cluster_table$k == k_val] <- cor_cluster_table$Freq[cor_cluster_table$marker == 'MX1']
# k_cluster_table$n_VEGFA[k_cluster_table$k == k_val] <- cor_cluster_table$Freq[cor_cluster_table$marker == 'VEGFA']
# k_cluster_table$n_AURKB[k_cluster_table$k == k_val] <- cor_cluster_table$Freq[cor_cluster_table$marker == 'AURKB']
# k_cluster_table$n_FOSB[k_cluster_table$k == k_val] <- cor_cluster_table$Freq[cor_cluster_table$marker == 'FOSB']
# k_cluster_table$n_BCL6[k_cluster_table$k == k_val] <- cor_cluster_table$Freq[cor_cluster_table$marker == 'BCL6']
# 
# }

# ggplot(k_cluster_table,aes(x = k,y = n_sig01)) +
#   geom_line()
# write.csv(k_cluster_table,'tonsil_LAIV_120a_s_Tcell_LAIV_marker_module_hcluster_cor_k_list.csv')
cor_cluster_table <- cor_cluster_table %>% arrange(Freq,cluster)
write.csv(cor_cluster_table,paste('tonsil_LAIV_120a_s_Tcell_LAIV_marker_module_hcluster_cor_',k_val,'.csv',sep = ''))
cluster_list <- sort(as.character(unique(cor_cluster_table$cluster[(cor_cluster_table$Freq >= 5) & (cor_cluster_table$Freq <= 40)])))
##### gene module visaluzation ############################
DefaultAssay(Tcell_LAIV.subset) <- 'normalized'
Tcell_LAIV.subset <- ScaleData(Tcell_LAIV.subset)
# cluster_list <- setdiff(cluster_list,c('35','5'))
cor_cluster_table$name <- ''
cor_cluster_table$name[cor_cluster_table$cluster == 43] <- 'IRF7Or9_module'
cor_cluster_table$name[cor_cluster_table$cluster == 40] <- 'CXCR5_module'
cor_cluster_table$name[cor_cluster_table$cluster == 37] <- 'IL32_module'
cor_cluster_table$name[cor_cluster_table$cluster == 22] <- 'MT_module'
cor_cluster_table$name[cor_cluster_table$cluster == 11] <- 'IL10_CD28_module'
cor_cluster_table$name[cor_cluster_table$cluster == 14] <- 'CD8_module'
cor_cluster_table$name[cor_cluster_table$cluster == 42] <- 'MX1_module'
cor_cluster_table$name[cor_cluster_table$cluster == 41] <- 'IFIT_module'
cor_cluster_table$name[cor_cluster_table$cluster == 30] <- 'BCL6_module'
cor_cluster_table$name[cor_cluster_table$cluster == 21] <- 'glucose2_module'
cor_cluster_table$name[cor_cluster_table$cluster == 31] <- 'activation_module'
cor_cluster_table$name[cor_cluster_table$cluster == 38] <- 'NFAT_module'
cor_cluster_table$name[cor_cluster_table$cluster == 15] <- 'KIR_module'
cor_cluster_table$name[cor_cluster_table$cluster == 5] <- 'CD161_module'
cor_cluster_table$name[cor_cluster_table$cluster == 26] <- 'glucose_module'
cor_cluster_table$name[cor_cluster_table$cluster == 16] <- 'RORC_exhaustion_module'
cor_cluster_table$name[cor_cluster_table$cluster == 28] <- 'prolif_module'
cor_cluster_table$name[cor_cluster_table$cluster == 24] <- 'STAT3_module'
cor_cluster_table$name[cor_cluster_table$cluster == 27] <- 'TBET_module'
cor_cluster_table$name[cor_cluster_table$cluster == 33] <- 'gdT_module'
cor_cluster_table$name[cor_cluster_table$cluster == 29] <- 'TLR9_module'

# saveRDS(cor_cluster_table,'Tcell_LAIV_module_cor_cluster_table.rds')
cluster_ttest_table <- expand.grid(cluster = cluster_list,days = c(0,4,7,10,14))
cluster_ttest_table <- merge(cluster_ttest_table,unique(cor_cluster_table[,c('cluster','Freq','name')]), by = 'cluster', all.x = T)
cluster_ttest_table$p <- NA
cluster_ttest_table$diff <- NA
cluster_ttest_table$FC <- NA
module_list <- unique(cor_cluster_table$name[cor_cluster_table$name != ''])
min_module_list <- c()
for (module_name in module_list) {
  graphics.off()
  feature_list <- cor_cluster_table$marker[cor_cluster_table$name == module_name]
  min_feature_list <- feature_list[feature_list %in% minimal_marker_list]
  Tcell_LAIV.subset <- AddModuleScore(
    Tcell_LAIV.subset,
    features = list(feature_list),
    name = module_name,
    assay = 'normalized',
    slot = 'slot.data',
    ctrl = 25
  )
  if (length(min_feature_list) > 1){
    min_module_list <- c(min_module_list,module_name)
    Tcell_LAIV.subset <- AddModuleScore(
      Tcell_LAIV.subset,
      features = list(min_feature_list),
      name = paste('min_',module_name,sep = ''),
      assay = 'normalized',
      slot = 'slot.data',
      ctrl = 25
    )
    temp <- data.frame(score = Tcell_LAIV.subset@meta.data[[paste('min_',module_name,'1',sep = '')]],
                       days = as.numeric(gsub('day','',Tcell_LAIV.subset$days)),
                       age_group = Tcell_LAIV.subset$age_group)
    temp <- temp[sample(1:dim(Tcell_LAIV.subset)[2],10000),]
    temp <- temp[temp$days != 8,]
    temp$day_group <- temp$days
    temp$day_group[temp$days %in% c(6,7)] <- 7
    temp$day_group[temp$days %in% c(12,14)] <- 14
    # for (day_name in c(0,4,7,10,14)){
    #   temp1 <- temp[temp$day_group == day_name,]
    #   wilcoxon_test <- wilcox.test(temp1$score[temp1$age_group == '27-39yrs'],
    #                   temp1$score[temp1$age_group == '02-03yrs'])
    #   cluster_ttest_table$p[(cluster_ttest_table$cluster == cluster_name) &
    #                                (cluster_ttest_table$days == day_name)] <- wilcoxon_test$p.value
    #   cluster_ttest_table$diff[(cluster_ttest_table$cluster == cluster_name) &
    #                           (cluster_ttest_table$days == day_name)] <-
    #     mean(temp1$score[temp1$age_group == '27-39yrs']) -
    #     mean(temp1$score[temp1$age_group == '02-03yrs'])
    #   cluster_ttest_table$FC[(cluster_ttest_table$cluster == cluster_name) &
    #                              (cluster_ttest_table$days == day_name)] <-
    #     cluster_ttest_table$diff[(cluster_ttest_table$cluster == cluster_name) &
    #                                (cluster_ttest_table$days == day_name)]/sd(temp1$score)
    # }
    # plot violin and trendline
    plot <- ggplot(temp, aes(x=day_group, y=score,color = age_group, linetype = age_group)) +
      geom_violin(aes(group = cut_width(day_group, 1)), scale = "width") +
      # geom_jitter(size = 0.1) +
      geom_smooth(method = "loess",se = F, linewidth=0.5) +
      scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
      scale_linetype_manual(values = c('02-03yrs' = 'dotted','07-09yrs'= 'longdash','27-39yrs' = 'solid')) +
      ggtitle(module_name) + theme_bw() +
      ylab('Module Score') +
      theme(text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold"))
    print(plot)
    dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_gene_module_k43_',module_name,'.pdf',sep = ''),width = 2.7, height = 1.5)
    # plot lines and trendlines
    # plot <- ggplot(temp, aes(x=day_group, y=score,color = age_group, linetype = age_group)) +
    #   # geom_violin(aes(group = cut_width(day_group, 1,color = age_group, linetype = age_group)), scale = "width") +
    #   # geom_jitter(size = 0.1) +
    #   geom_smooth(method = "loess",se = F, linewidth=0.5) 
    # for (feature_name in feature_list){
    #   Tcell_LAIV.subset <- AddModuleScore(
    #     Tcell_LAIV.subset,
    #     features = list(c(feature_name)),
    #     name = 'temp',
    #     assay = 'normalized',
    #     slot = 'slot.data',
    #     ctrl = 25
    #   )
    #   temp1 <- data.frame(index = colnames(Tcell_LAIV.subset),
    #                       temp = Tcell_LAIV.subset@meta.data[['temp1']])
    #   p <- match(rownames(temp),temp1$index)
    #   temp$gene <- temp1$temp[p]
    #   plot <- plot + geom_line(data = temp,stat = "smooth", method = 'loess', aes(x = day_group,y = gene),se = F, linewidth=0.3, alpha = 0.5)
    #   # plot <- plot + geom_smooth(data = temp,aes(x = day_group,y = gene),method = "loess",se = F, linewidth=0.1)#,color = 'black'
    # }
    # plot <- plot + 
    #   scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
    #   scale_linetype_manual(values = c('02-03yrs' = 'dotted','07-09yrs'= 'longdash','27-39yrs' = 'solid')) +
    #   ggtitle(module_name) + theme_bw() +
    #   ylab('Module Score') +
    #   theme(text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold"))
    # print(plot)
    # dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_gene_module_',module_name,'_lines.pdf',sep = ''),width = 2.7, height = 1.5)
    # 
  }  


  
}
# cluster_ttest_table$p_adj <- p.adjust(cluster_ttest_table$p,method = 'fdr')
# cluster_ttest_table <- cluster_ttest_table %>% arrange(p_adj)
# write.xlsx(cluster_ttest_table,'tonsil_LAIV_120a_s_Tcell_LAIV_gene_module_wilcoxon_test.xlsx')
# p_cutoff <- 0.05
# cluster_ttest_table$label <- paste('day',cluster_ttest_table$days,'_',cluster_ttest_table$name,sep = '')
# cluster_ttest_table$significant <- 'no'
# cluster_ttest_table$significant[(cluster_ttest_table$p_adj <= p_cutoff)] <- 'yes'
# ggplot(cluster_ttest_table, aes(x=FC, y=-log10(p_adj),color = significant)) +
#   geom_point(size = 1,alpha = 0.5) +
#   scale_colour_manual(values = c("yes" = "red","no" = "darkgrey")) +
#   ggtitle(paste('DEG of Tcells not presented at baseline\nhigher in adults higher in toddlers\n')) +
#   xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
#   theme(text = element_text(size = 10)) +
#   theme_bw() +
#   geom_text_repel(aes(x=FC, y=-log10(p_adj),label = label), color = "black", size = 3) +
#   # geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
#   geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
# # dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n300_toddler_vs_adult_big.pdf',sep = ''),width = 15, height = 15)
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_gene_module_volcano.pdf',sep = ''),width = 6.5, height = 5)
min_module_list <- colnames(Tcell_LAIV.subset@meta.data)
min_module_list <- min_module_list[grepl('module',min_module_list) & grepl('^min_',min_module_list)]
# plot UMAP
for (module_name in min_module_list){
  # Tcell_LAIV.subset@meta.data[,module_name] <- NULL
  graphics.off()
  KO_cells <- 1:dim(Tcell_LAIV.subset)[2]
  downsampled_KO_cells <- sample(KO_cells, dim(Tcell_LAIV.subset)[2]/3)
  temp <- Tcell_LAIV.subset[,downsampled_KO_cells]
  # DimPlot(WT_KO_integrated_downsampled, reduction = "umap", split.by ="orig.ident", ncol=2)
  plot <- FeaturePlot(temp, features = module_name,pt.size = 0.01, ncol = 1, sort.cell = TRUE,min.cutoff = 0)+#,max.cutoff = 4)+
    theme(text = element_text(size = 7),axis.text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_k43_umap_',module_name,'.pdf',sep = ''),width = 3, height = 2.3)
  
}

module_name <- 'CD8_module'
temp <- data.frame(score = Tcell_LAIV.subset@meta.data[[paste('min_',module_name,'1',sep = '')]],
                   days = as.numeric(gsub('day','',Tcell_LAIV.subset$days)),
                   age_group = Tcell_LAIV.subset$age_group)
temp <- temp[sample(1:dim(Tcell_LAIV.subset)[2],10000),]
temp <- temp[temp$days != 8,]
temp$day_group <- temp$days
temp$day_group[temp$days %in% c(6,7)] <- 7
temp$day_group[temp$days %in% c(12,14)] <- 14
# for (day_name in c(0,4,7,10,14)){
#   temp1 <- temp[temp$day_group == day_name,]
#   wilcoxon_test <- wilcox.test(temp1$score[temp1$age_group == '27-39yrs'],
#                   temp1$score[temp1$age_group == '02-03yrs'])
#   cluster_ttest_table$p[(cluster_ttest_table$cluster == cluster_name) &
#                                (cluster_ttest_table$days == day_name)] <- wilcoxon_test$p.value
#   cluster_ttest_table$diff[(cluster_ttest_table$cluster == cluster_name) &
#                           (cluster_ttest_table$days == day_name)] <-
#     mean(temp1$score[temp1$age_group == '27-39yrs']) -
#     mean(temp1$score[temp1$age_group == '02-03yrs'])
#   cluster_ttest_table$FC[(cluster_ttest_table$cluster == cluster_name) &
#                              (cluster_ttest_table$days == day_name)] <-
#     cluster_ttest_table$diff[(cluster_ttest_table$cluster == cluster_name) &
#                                (cluster_ttest_table$days == day_name)]/sd(temp1$score)
# }
# plot violin and trendline
plot <- ggplot(temp, aes(x=day_group, y=score,color = age_group, linetype = age_group)) +
  geom_violin(aes(group = cut_width(day_group, 1)), scale = "width") +
  # geom_jitter(size = 0.1) +
  scale_y_log10() +
  geom_smooth(method = "loess",se = F, linewidth=0.5) +
  scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
  scale_linetype_manual(values = c('02-03yrs' = 'dotted','07-09yrs'= 'longdash','27-39yrs' = 'solid')) +
  ggtitle(module_name) + theme_bw() +
  ylab('Module Score') +
  theme(text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold"))
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_gene_module_',module_name,'.pdf',sep = ''),width = 2.7, height = 1.5)
##### gene and module visualization #########################
module_name <- 'activation_module'
day_name <- 7
temp <- data.frame(score = Tcell_LAIV.subset@meta.data[[paste('min_',module_name,'1',sep = '')]],
                   days = as.numeric(gsub('day','',Tcell_LAIV.subset$days)),
                   age_group = Tcell_LAIV.subset$age_group)
temp$age_group[temp$age_group == '02-03yrs'] <- 'Toddlers'
temp$age_group[temp$age_group == '07-09yrs'] <- 'Older Children'
temp$age_group[temp$age_group == '27-39yrs'] <- 'Adults'
temp$age_group <- factor(temp$age_group,levels = c('Toddlers','Older Children','Adults'))
temp <- temp[temp$days != 8,]
temp$day_group <- temp$days
temp$day_group[temp$days %in% c(6,7)] <- 7
temp$day_group[temp$days %in% c(12,14)] <- 14
temp <- temp[temp$day_group == day_name,]

label_table <- data.frame(x = 2)
label_table$y <- min(temp$score) + 1.2*(max(temp$score) - min(temp$score))
label_table$label <- 'p_adj = 0.000'
# label_table$sig <- 'yes'#ifelse(label_table$p < 0.05,'yes','no')

ggplot(temp, aes(x=age_group, y=score,)) + 
  geom_violin(aes(color = age_group)) +
  geom_boxplot(aes(color = age_group),outlier.size = 0.1,width=0.1) + 
  # geom_jitter(size = 0.1) +
  scale_colour_manual(values = c('Toddlers' = '#F16C23','Older Children'= '#1b7c3d','Adults' = '#2b6a99')) +
  ggtitle(paste(module_name,'on day',day_name)) + theme_bw() +
  ylab('Module Score') + RotatedAxis() + 
  geom_text(label_table,mapping = aes(x = 2,y = y,label = label,fill = NULL),color = 'red',
            size = 5/14*7) +
  ylim(min(temp$score),min(temp$score) + 1.3*(max(temp$score) - min(temp$score))) +
  theme(text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_gene_module_min_',module_name,'_day',day_name,'.pdf',sep = ''),width = 2.4, height = 1.7)


#### partial correlation network #####################
library('bootnet')
library("readxl")
library("xlsx")
library('qgraph')
library('ppcor')
library(huge)
library('NetworkComparisonTest')
library(igraph)
library(RCy3)
library(reshape2)
temp <- Tcell_LAIV.subset@meta.data
temp <- data.frame(temp[,grepl('^min_',colnames(temp))],
                   days = as.numeric(gsub('day','',Tcell_LAIV.subset$days)),
                   age_group = Tcell_LAIV.subset$age_group,
                   condition = Tcell_LAIV.subset$condition)
colnames(temp) <- gsub('min_','',gsub('1$','',colnames(temp)))
temp <- temp[temp$days != 8,]
temp$day_group <- temp$days
temp$day_group[temp$days %in% c(6,7)] <- 7
temp$day_group[temp$days %in% c(12,14)] <- 14
# saveRDS(temp,'Tcell_LAIV_gene_module_table.rds')
day_name <- 7
temp1 <- temp[temp$day_group == day_name,]
temp2 <- temp1[temp1$age_group == '27-39yrs',]
my_data <- data.frame(temp2[,grepl('module',colnames(temp2))])#dim(selected_table_day0_day04_day07)[2]]
my_data_npn = huge.npn(my_data) # normalization?
results_Adults <- estimateNetwork(my_data_npn,default = 'EBICglasso',corMethod = 'cor_auto', tuning = 0.5)

temp2 <- temp1[temp1$age_group == '02-03yrs',]
my_data <- data.frame(temp2[,grepl('module',colnames(temp2))])#dim(selected_table_day0_day04_day07)[2]]
my_data_npn = huge.npn(my_data) # normalization?
results_Toddlers <- estimateNetwork(my_data_npn,default = 'EBICglasso',corMethod = 'cor_auto', tuning = 0.5)

temp2 <- temp1[temp1$age_group == '07-09yrs',]
my_data <- data.frame(temp2[,grepl('module',colnames(temp2))])#dim(selected_table_day0_day04_day07)[2]]
my_data_npn = huge.npn(my_data) # normalization?
results_Older_Children <- estimateNetwork(my_data_npn,default = 'EBICglasso',corMethod = 'cor_auto', tuning = 0.5)

temp <- Tcell_LAIV_module_wilcoxon_youth_vs_adult[Tcell_LAIV_module_wilcoxon_youth_vs_adult$days == day_name,]
# graphics.off()

network <- graph_from_adjacency_matrix(results_Adults[["graph"]], weighted=T, mode="undirected", diag=F)
V(network)$degree <- degree(network)
# p <- match(names(V(network)),temp$module)
# V(network)$FC <- temp$FC_std_mean[p]
V(network)$name <- gsub('_module','',V(network)$name)
createNetworkFromIgraph(network)

# boot1 <- bootnet(results, nCores = 8, nBoots = 1000, type = 'nonparametric')
# boot2 <- bootnet(results, nCores = 8, nBoots = 1000, type = 'case')


NCT_comparison <- NCT(results_Toddlers, results_Adults, it=5000, test.edges=T, edges = "all",progressbar=TRUE)
NCT_output <- NCT_comparison[["einv.pvals"]]
colnames(NCT_output)[3] <- 'p'
NCT_output$fdr <- p.adjust(NCT_output$p,method = 'fdr')
NCT_output <- NCT_output %>% arrange(p)

# NCT_output <- read.csv('Tcell_LAIV_gene_module_network_diff_day0_age_NCT.csv')
temp1 <- data.frame(results_Adults$graph)
temp1$Var1 <- rownames(temp1)
temp1 <- melt(temp1,id.vars = 'Var1', variable.name = 'Var2',value.name = 'adults')
temp2 <- data.frame(results_Toddlers$graph)
temp2$Var1 <- rownames(temp2)
temp2 <- melt(temp2,id.vars = 'Var1', variable.name = 'Var2',value.name = 'toddlers')
temp <- merge(temp1,temp2,by = c('Var1','Var2'),all.x = T)
temp$diff <- temp$adults - temp$toddlers
temp$FC <- temp$diff/abs(temp$toddlers)

NCT_output <- merge(NCT_output,temp,by = c('Var1','Var2'),all.x = T)
NCT_output$sig <- 'no'
NCT_output$sig[NCT_output$fdr <= 0.1] <- 'yes'
NCT_output <- NCT_output %>% arrange(desc(sig),desc(FC),desc(diff))
write.xlsx(NCT_output,'Tcell_LAIV_gene_module_network_diff_day07_age_NCT.xlsx')

NCT_output_day07 <- read.xlsx('Tcell_LAIV_gene_module_network_diff_day07_age_NCT.xlsx',sheetName = 'Sheet1')
NCT_output_day0 <- NCT_output
NCT_output_day07$name <- paste(NCT_output_day07$Var1,NCT_output_day07$Var2,sep = '_')
NCT_output_day0$name <- paste(NCT_output_day0$Var1,NCT_output_day0$Var2,sep = '_')

temp1 <- NCT_output_day07$name[(NCT_output_day07$sig == 'yes') & (NCT_output_day07$diff > 0)]
temp2 <- NCT_output_day0$name[(NCT_output_day0$sig == 'yes') & (NCT_output_day0$diff > 0)]


##### 
DEGs <- as.character(unique(Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$marker[Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$day_group != 'day0']))
sum(DEGs %in% cor_cluster_table$marker[cor_cluster_table$Freq < 5])
sum(DEGs %in% cor_cluster_table$marker[cor_cluster_table$Freq > 31])
sum(DEGs %in% cor_cluster_table$marker[(cor_cluster_table$Freq <= 31) & (cor_cluster_table$Freq >= 5) & (cor_cluster_table$name != '')])