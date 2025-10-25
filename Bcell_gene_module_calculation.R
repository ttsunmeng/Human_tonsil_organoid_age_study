library(ComplexHeatmap)
library(dendextend)
library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(ggrepel)
library("readxl")
library('stringr')
library(reshape2)

batch_marker_list <- list()
batch_marker_list[[1]] <- readRDS('tonsil_LAIV_120a.rds')
batch_marker_list[[2]] <- readRDS('tonsil_LAIV_120b.rds')
batch_marker_list[[3]] <- readRDS('tonsil_LAIV_120c.rds')
batch_marker_list[[4]] <- readRDS('tonsil_LAIV_120j.rds')
minimal_marker_list <- Reduce(intersect,batch_marker_list[2:4])

Bcell_LAIV.subset <- subset(Bcell_LAIV.subset,donor_ID != '02yrs M IMD085')

Bcell_dataframe_lognorm <- data.frame(t(as.matrix(Bcell_LAIV.subset@assays[["normalized"]]@data)),check.names = F)
Bcell_dataframe_lognorm <- Bcell_dataframe_lognorm[,colSums(Bcell_dataframe_lognorm) != 0]
cor_mat <- cor(Bcell_dataframe_lognorm[,colnames(Bcell_dataframe_lognorm) %in% minimal_marker_list],use = 'complete.obs')#[selected_table_day0_day04_day07$age_group == 'Adults',]

HM2 <- Heatmap(cor_mat, cluster_columns = T, cluster_rows = T,
               name = "mean\nnormalized\nexpression", show_column_names = T,
               show_row_names = F,#row_km = 20, column_km = 20,
               column_title_side = "top",row_title_side = "left",row_names_side = 'left')
HM2 <- draw(HM2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_LAIV_gene_module_cor.pdf',width = 80, height = 80)

r.dend <- row_dend(HM2)  #Extract row dendrogram
cor_cluster_table <- cutree(r.dend, h = 1.5)
sort(table(cor_cluster_table))
# library(factoextra)
# library(cluster)
# fviz_nbclust(cor_mat, k.max = 40,FUN = hcut, method = "wss")
# dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_gene_module_cor_hcluster_elbow.pdf',width = 10, height = 5)
# fviz_nbclust(cor_mat, k.max = 40,FUN = hcut, method = "silhouette")
# dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_gene_module_corhcluster_silhouette.pdf',width = 10, height = 5)
# gap_stat <- clusGap(cor_mat, FUN = hcut, nstart = 25, K.max = 40, B = 50)
# fviz_gap_stat(gap_stat)
# dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_gene_module_corhcluster_GasStat.pdf',width = 10, height = 5)

# k_cluster_table <- data.frame(k = c(3:200),n_count5 = NA, n_MX1 = NA,n_VEGFA = NA,n_AURKB = NA,n_FOSB = NA,n_BCL6 = NA)
# for (k_val in k_cluster_table$k) {
k_val <- 51
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
# write.csv(k_cluster_table,'tonsil_LAIV_120a_s_Bcell_LAIV_marker_module_hcluster_cor_k_list.csv')
cor_cluster_table <- cor_cluster_table %>% arrange(Freq,cluster)
write.csv(cor_cluster_table,paste('tonsil_LAIV_120a_s_Bcell_LAIV_min_marker_module_hcluster_cor_',k_val,'.csv',sep = ''))
cluster_list <- sort(as.character(unique(cor_cluster_table$cluster[cor_cluster_table$Freq >= 5])))
##### gene module visaluzation ############################
DefaultAssay(Bcell_LAIV.subset) <- 'normalized'
Bcell_LAIV.subset <- ScaleData(Bcell_LAIV.subset)
cluster_list <- setdiff(cluster_list,c('1','2','6','22','25'))
cor_cluster_table$name <- ''
cor_cluster_table$name[cor_cluster_table$cluster == 41] <- 'glucose_module'
cor_cluster_table$name[cor_cluster_table$cluster == 39] <- 'naive_module'
cor_cluster_table$name[cor_cluster_table$cluster == 38] <- 'activation_module'
cor_cluster_table$name[cor_cluster_table$cluster == 47] <- 'glucose_module2'
cor_cluster_table$name[cor_cluster_table$cluster == 45] <- 'ISG_module'
cor_cluster_table$name[cor_cluster_table$cluster == 43] <- 'PB_module'
cor_cluster_table$name[cor_cluster_table$cluster == 20] <- 'IRF_module'
cor_cluster_table$name[cor_cluster_table$cluster == 2] <- 'CD16_module'
cor_cluster_table$name[cor_cluster_table$cluster == 31] <- 'activation_module2'
cor_cluster_table$name[cor_cluster_table$cluster == 32] <- 'presenting_module'
cor_cluster_table$name[cor_cluster_table$cluster == 30] <- 'BCL6_module'
cor_cluster_table$name[cor_cluster_table$cluster == 37] <- 'LZ_module'
cor_cluster_table$name[cor_cluster_table$cluster == 33] <- 'switched_module'
cor_cluster_table$name[cor_cluster_table$cluster == 23] <- 'switched_module2'
cor_cluster_table$name[cor_cluster_table$cluster == 24] <- 'MX1_module'
cor_cluster_table$name[cor_cluster_table$cluster == 44] <- 'IFN_module'
cor_cluster_table$name[cor_cluster_table$cluster == 46] <- 'NFkB_module'
cor_cluster_table$name[cor_cluster_table$cluster == 29] <- 'prolif_module'
cor_cluster_table$name[cor_cluster_table$cluster == 27] <- 'NFAT_module'
cor_cluster_table$name[cor_cluster_table$cluster == 21] <- 'TBET_module'
cor_cluster_table$name[cor_cluster_table$cluster == 28] <- 'MYC_module'
cor_cluster_table$name[cor_cluster_table$cluster == 26] <- 'CXCR5_module'
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
  # Bcell_LAIV.subset <- AddModuleScore(
  #   Bcell_LAIV.subset,
  #   features = list(feature_list),
  #   name = module_name,
  #   assay = 'normalized',
  #   slot = 'slot.data',
  #   ctrl = 25
  # )
  if (length(min_feature_list) > 1){
    min_module_list <- c(min_module_list,module_name)
  #   Bcell_LAIV.subset <- AddModuleScore(
  #     Bcell_LAIV.subset,
  #     features = list(min_feature_list),
  #     name = paste('min_',module_name,sep = ''),
  #     assay = 'normalized',
  #     slot = 'slot.data',
  #     ctrl = 25
  #   )
  }
 #  temp <- data.frame(score = Bcell_LAIV.subset@meta.data[[paste('min_',module_name,'1',sep = '')]],
 #                     days = as.numeric(gsub('day','',Bcell_LAIV.subset$days)),
 #                     age_group = Bcell_LAIV.subset$age_group)
 #  temp <- temp[sample(1:dim(Bcell_LAIV.subset)[2],10000),]
 #  temp <- temp[temp$days != 8,]
 #  temp$day_group <- temp$days
 #  temp$day_group[temp$days %in% c(6,7)] <- 7
 #  temp$day_group[temp$days %in% c(12,14)] <- 14
 #  # for (day_name in c(0,4,7,10,14)){
 #  #   temp1 <- temp[temp$day_group == day_name,]
 #  #   wilcoxon_test <- wilcox.test(temp1$score[temp1$age_group == '27-39yrs'],
 #  #                   temp1$score[temp1$age_group == '02-03yrs'])
 #  #   cluster_ttest_table$p[(cluster_ttest_table$cluster == cluster_name) & 
 #  #                                (cluster_ttest_table$days == day_name)] <- wilcoxon_test$p.value
 #  #   cluster_ttest_table$diff[(cluster_ttest_table$cluster == cluster_name) & 
 #  #                           (cluster_ttest_table$days == day_name)] <- 
 #  #     mean(temp1$score[temp1$age_group == '27-39yrs']) - 
 #  #     mean(temp1$score[temp1$age_group == '02-03yrs'])
 #  #   cluster_ttest_table$FC[(cluster_ttest_table$cluster == cluster_name) & 
 #  #                              (cluster_ttest_table$days == day_name)] <- 
 #  #     cluster_ttest_table$diff[(cluster_ttest_table$cluster == cluster_name) & 
 #  #                                (cluster_ttest_table$days == day_name)]/sd(temp1$score)
 #  # }
 # # plot violin and trendline
 #  plot <- ggplot(temp, aes(x=day_group, y=score,color = age_group, linetype = age_group)) +
 #    geom_violin(aes(group = cut_width(day_group, 1)), scale = "width") +
 #    # geom_jitter(size = 0.1) +
 #    geom_smooth(method = "loess",se = F, linewidth=0.5) +
 #    scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
 #    scale_linetype_manual(values = c('02-03yrs' = 'dotted','07-09yrs'= 'longdash','27-39yrs' = 'solid')) +
 #    ggtitle(module_name) + theme_bw() +
 #    ylab('Module Score') +
 #    theme(text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold"))
 #  print(plot)
 #  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_gene_module_',module_name,'.pdf',sep = ''),width = 2.7, height = 1.5)
  # plot lines and trendlines
  # plot <- ggplot(temp, aes(x=day_group, y=score,color = age_group, linetype = age_group)) +
  #   # geom_violin(aes(group = cut_width(day_group, 1,color = age_group, linetype = age_group)), scale = "width") +
  #   # geom_jitter(size = 0.1) +
  #   geom_smooth(method = "loess",se = F, linewidth=0.5) 
  # for (feature_name in feature_list){
  #   Bcell_LAIV.subset <- AddModuleScore(
  #     Bcell_LAIV.subset,
  #     features = list(c(feature_name)),
  #     name = 'temp',
  #     assay = 'normalized',
  #     slot = 'slot.data',
  #     ctrl = 25
  #   )
  #   temp1 <- data.frame(index = colnames(Bcell_LAIV.subset),
  #                       temp = Bcell_LAIV.subset@meta.data[['temp1']])
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
  # dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_gene_module_',module_name,'_lines.pdf',sep = ''),width = 2.7, height = 1.5)
  # 
}
cluster_ttest_table$p_adj <- p.adjust(cluster_ttest_table$p,method = 'fdr')
cluster_ttest_table <- cluster_ttest_table %>% arrange(p_adj)
write.xlsx(cluster_ttest_table,'tonsil_LAIV_120a_s_Bcell_LAIV_gene_module_wilcoxon_test.xlsx')
p_cutoff <- 0.05
cluster_ttest_table$label <- paste('day',cluster_ttest_table$days,'_',cluster_ttest_table$name,sep = '')
cluster_ttest_table$significant <- 'no'
cluster_ttest_table$significant[(cluster_ttest_table$p_adj <= p_cutoff)] <- 'yes'
ggplot(cluster_ttest_table, aes(x=FC, y=-log10(p_adj),color = significant)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_colour_manual(values = c("yes" = "red","no" = "darkgrey")) +
  ggtitle(paste('DEG of Bcells not presented at baseline\nhigher in adults higher in toddlers\n')) +
  xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 10)) +
  theme_bw() +
  geom_text_repel(aes(x=FC, y=-log10(p_adj),label = label), color = "black", size = 3) +
  # geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_wilcoxon_bootstrap_n300_toddler_vs_adult_big.pdf',sep = ''),width = 15, height = 15)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_gene_module_volcano.pdf',sep = ''),width = 6.5, height = 5)

min_module_list <- colnames(Bcell_LAIV.subset@meta.data)
min_module_list <- min_module_list[grepl('module',min_module_list) & grepl('^min_',min_module_list)]
min_module_list <- gsub('1$','',gsub('^min_','',min_module_list))
for (module_name in min_module_list) {
  graphics.off()
   temp <- data.frame(score = Bcell_LAIV.subset@meta.data[[paste('min_',module_name,'1',sep = '')]],
                      days = as.numeric(gsub('day','',Bcell_LAIV.subset$days)),
                      age_group = Bcell_LAIV.subset$age_group)
   temp <- temp[sample(1:dim(Bcell_LAIV.subset)[2],10000),]
   temp <- temp[temp$days != 8,]
   temp$day_group <- temp$days
   temp$day_group[temp$days %in% c(6,7)] <- 7
   temp$day_group[temp$days %in% c(12,14)] <- 14
   temp <- temp %>% group_by(age_group,day_group) %>% mutate(mean = mean(score))

  # plot violin and trendline
   plot <- ggplot(temp, aes(x=day_group, y=score,color = age_group, linetype = age_group,shape = age_group,fill = age_group)) +
     # geom_point(aes(x = day_group, y = mean),size = 2) + # 21 is filled circle
     geom_smooth(method = "loess",se = F, linewidth=0.5) +
     scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
     # scale_fill_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
     scale_linetype_manual(values = c('02-03yrs' = 'dotted','07-09yrs'= 'longdash','27-39yrs' = 'solid')) +
     # scale_shape_manual(values = c('02-03yrs' = 25,'07-09yrs'= 23,'27-39yrs' = 22)) +
     ggtitle(module_name) + theme_bw() +
     ylab('Module Score') +
     theme(text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold"))
   print(plot)
   dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_gene_module_lines_',module_name,'.pdf',sep = ''),width = 2.7, height = 1.5)
}
##### gene and module visualization #########################
Bcell_dataframe_lognorm_scale <- Bcell_LAIV.subset@assays[["normalized"]]@scale.data
Bcell_dataframe_lognorm_scale <- data.frame(t(Bcell_dataframe_lognorm_scale))
cluster_name <- 41
graphics.off()
feature_list <- cor_cluster_table$marker[cor_cluster_table$cluster == cluster_name]
module_name <- cor_cluster_table$name[cor_cluster_table$cluster == cluster_name][1]
temp <- data.frame(score = Bcell_LAIV.subset@meta.data[[paste(module_name,'1',sep = '')]],
                   days = as.numeric(gsub('day','',Bcell_LAIV.subset$days)),
                   age_group = Bcell_LAIV.subset$age_group)
# temp1 <- Bcell_dataframe_lognorm_scale[,feature_list]
# temp <- cbind(temp,temp1)
temp1 <- sample(1:dim(Bcell_LAIV.subset)[2],10000)
temp <- temp[temp1,]
temp <- temp[temp$days != 8,]
temp$day_group <- temp$days
temp$day_group[temp$days %in% c(6,7)] <- 7
temp$day_group[temp$days %in% c(12,14)] <- 14

##### heatmap visualization #########################
temp <- cor_cluster_table[cor_cluster_table$cluster %in% cluster_list,]
temp <- Bcell_dataframe_lognorm[,(colnames(Bcell_dataframe_lognorm) %in% cor_cluster_table$marker[cor_cluster_table$cluster %in% cluster_list])]
cor_mat <- cor(temp,use = 'complete.obs')#[selected_table_day0_day04_day07$age_group == 'Adults',]

HM2 <- Heatmap(cor_mat, cluster_columns = T, cluster_rows = T,
               name = "mean\nnormalized\nexpression", show_column_names = F,
               show_row_names = F,#row_km = 20, column_km = 20,
               column_title_side = "top",row_title_side = "left",row_names_side = 'left')
HM2 <- draw(HM2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_LAIV_gene_module_cor_selected.pdf',width = 4, height = 4)

#### module volcano visualization #########################
module_name <- 'MX1_module'
day_name <- 4
temp <- data.frame(score = Bcell_LAIV.subset@meta.data[[paste(module_name,'1',sep = '')]],
                   days = as.numeric(gsub('day','',Bcell_LAIV.subset$days)),
                   age_group = Bcell_LAIV.subset$age_group)
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
label_table$label <- paste('p_adj =',sprintf(cluster_ttest_table$p_adj[(cluster_ttest_table$name == module_name) & 
                                                                       (cluster_ttest_table$days == day_name)], fmt = '%#.3f'))
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
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_gene_module_',module_name,'_day',day_name,'.pdf',sep = ''),width = 2.4, height = 1.7)

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
temp <- data.frame(project_ID = Bcell_LAIV.subset$project_ID, 
                   cell_index = gsub('_.*','',colnames(Bcell_LAIV.subset)),
                   donorID = Bcell_LAIV.subset$donor_ID, 
                   condition = Bcell_LAIV.subset$condition, 
                   Bcell_LAIV.subset@meta.data[,paste('min_',min_module_list,'1',sep = '')],
                   days = as.numeric(gsub('day','',Bcell_LAIV.subset$days)),
                   age_group = Bcell_LAIV.subset$age_group)
colnames(temp) <- gsub('min_','',gsub('1$','',colnames(temp)))
write.csv(temp,'tonsil_LAIV_120a_s_Bcell_log1p_cca_CycleRegressOut_min_module.csv')

temp <- temp[temp$days != 8,]
temp$day_group <- temp$days
temp$day_group[temp$days %in% c(6,7)] <- 7
temp$day_group[temp$days %in% c(12,14)] <- 14
# saveRDS(temp,'Bcell_LAIV_gene_module_table.rds')
day_name <- 0
temp1 <- temp[temp$day_group == day_name,]
temp2 <- temp1[temp1$age_group == '27-39yrs',]
my_data <- data.frame(temp2[,grepl('module',colnames(temp2))])#dim(selected_table_day0_day04_day07)[2]]
my_data_npn = huge.npn(my_data) # normalization?
results_Adults <- estimateNetwork(my_data_npn,default = 'EBICglasso',corMethod = 'cor_auto', tuning = 0.5)

temp2 <- temp1[temp1$age_group == '02-03yrs',]
my_data <- data.frame(temp2[,grepl('module',colnames(temp2))])#dim(selected_table_day0_day04_day07)[2]]
my_data_npn = huge.npn(my_data) # normalization?
results_Toddlers <- estimateNetwork(my_data_npn,default = 'EBICglasso',corMethod = 'cor_auto', tuning = 0.5)

temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult[Bcell_LAIV_module_wilcoxon_youth_vs_adult$days == day_name,]
graphics.off()

network <- graph_from_adjacency_matrix(results[["graph"]], weighted=T, mode="undirected", diag=F)
V(network)$degree <- degree(network)
p <- match(V(network)$name,temp$module)
V(network)$FC <- temp$FC_mean[p]
V(network)$name <- gsub('_module','',V(network)$name)
createNetworkFromIgraph(network)

# boot1 <- bootnet(results, nCores = 8, nBoots = 1000, type = 'nonparametric')
# boot2 <- bootnet(results, nCores = 8, nBoots = 1000, type = 'case')
min_module_list <- colnames(Bcell_LAIV.subset@meta.data)
min_module_list <- min_module_list[grepl('module',min_module_list) & grepl('^min_',min_module_list)]
# plot UMAP
for (module_name in min_module_list){
  graphics.off()
  plot <- FeaturePlot(Bcell_LAIV.subset, features = paste('min_',module_name,'1',sep = ''),pt.size = 0.01, ncol = 1, sort.cell = TRUE,min.cutoff = 0)+#,max.cutoff = 4)+
    theme(text = element_text(size = 7),axis.text = element_text(size = 7),plot.title = element_text(size = 7, face = "bold"))
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_umap_min_',module_name,'.pdf',sep = ''),width = 3, height = 2.3)
}  

NCT_comparison <- NCT(results_Toddlers, results_Adults, it=5000, test.edges=T, edges = "all",progressbar=TRUE)
NCT_output <- NCT_comparison[["einv.pvals"]]
colnames(NCT_output)[3] <- 'p'
NCT_output$fdr <- p.adjust(NCT_output$p,method = 'fdr')
NCT_output <- NCT_output %>% arrange(p)

NCT_output <- read.csv('Bcell_LAIV_gene_module_network_diff_day0_age_NCT.csv')
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
write.xlsx(NCT_output,'Bcell_LAIV_gene_module_network_diff_day0_age_NCT.xlsx')

NCT_output_day07 <- read.xlsx('Bcell_LAIV_gene_module_network_diff_day07_age_NCT.xlsx',sheetName = 'Sheet1')
NCT_output_day0 <- NCT_output
NCT_output_day07$name <- paste(NCT_output_day07$Var1,NCT_output_day07$Var2,sep = '_')
NCT_output_day0$name <- paste(NCT_output_day0$Var1,NCT_output_day0$Var2,sep = '_')

temp1 <- NCT_output_day07$name[(NCT_output_day07$sig == 'yes') & (NCT_output_day07$diff > 0)]
temp2 <- NCT_output_day0$name[(NCT_output_day0$sig == 'yes') & (NCT_output_day0$diff > 0)]


##### 
DEGs <- as.character(unique(Bcell_LAIV_wilcoxon_youth_vs_adult_downregulation$marker[Bcell_LAIV_wilcoxon_youth_vs_adult_downregulation$day_group != 'day0']))
sum(DEGs %in% cor_cluster_table$marker[cor_cluster_table$Freq < 5])
sum(DEGs %in% cor_cluster_table$marker[cor_cluster_table$Freq > 31])
sum(DEGs %in% cor_cluster_table$marker[(cor_cluster_table$Freq <= 31) & (cor_cluster_table$Freq >= 5) & (cor_cluster_table$name != '')])


##### visualize dynamics distribution ##########
temp <- data.frame(Bcell_LAIV.subset@meta.data[,paste('min_',min_module_list,'1',sep = '')],
                   days = as.numeric(gsub('day','',Bcell_LAIV.subset$days)),
                   age_group = Bcell_LAIV.subset$age_group)
colnames(temp) <- gsub('min_','',gsub('1$','',colnames(temp)))
temp <- temp[temp$days != 8,]
temp$day_group <- temp$days
temp$day_group[temp$days %in% c(6,7)] <- 7
temp$day_group[temp$days %in% c(12,14)] <- 14
# scale
scale_temp <- temp
# scale_temp[,grepl('_module',colnames(scale_temp))] <- scale(scale_temp[,grepl('_module',colnames(scale_temp))])

# calculating average by day groups and age groups
scale_temp <- melt(scale_temp, id.vars=c("age_group", "day_group",'days'),variable.name="module",value.name="score")
scale_temp$module <- gsub('_module','',scale_temp$module)
average_temp <- 
  scale_temp %>% group_by(day_group,age_group,module) %>% summarize(mean = mean(score, na.rm=TRUE))
# add up
average_temp$padded_mean <- average_temp$mean - min(average_temp$mean)
average_temp <- 
  average_temp %>% group_by(day_group,age_group) %>% mutate(total = sum(padded_mean))
average_temp$mean_percentage <- average_temp$padded_mean/average_temp$total*100
# visualize
# library(purrr) 
# library(ggalluvial)

#line up using day07 adult data for colors
temp <- average_temp[(average_temp$day_group == 7) & (average_temp$age_group == '27-39yrs'),]
temp <- temp %>% arrange(mean)
average_temp$module <- factor(average_temp$module,levels = temp$module)
# ggplot(average_temp, aes(x = day_group, stratum = age_group, alluvium = module, y = padded_mean, label = module,fill = module)) +#, 
#   geom_alluvium(aes(fill = module), width = 1,alpha = 1,color = 'black',size = 0.1) +
#   facet_grid(age_group ~ .,scales = "fixed") +
#   # geom_stratum(width = 1/12, fill = "black",size = 0.5) +
#   # scale_fill_manual(values = cols_group) +
#   ggtitle('') + ylab('#cells') + #ylim(0,300) + xlim(0,14)
#   theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
# # dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_age_group_cluster_labeled_count.pdf',sep = ''),width = 5, height = 5)
# 
# ggplot(average_temp, aes(x = day_group, stratum = age_group, alluvium = module, y = mean_percentage, label = module,fill = module)) +#, 
#   geom_alluvium(aes(fill = module), width = 0.5,alpha = 1,color = 'black',size = 0.1) +
#   facet_grid(age_group ~ .,scales = "fixed") +
#   # geom_stratum(width = 1/12, fill = "black") +
#   # scale_fill_manual(values = cols_group) +
#   ggtitle('') + ylab('% in T cells') + #ylim(0,300) + xlim(0,14)
#   theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
# # dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_age_group_cluster_labeled_fraction.pdf',sep = ''),width = 5, height = 5)

ggplot(average_temp, aes(x = day_group, y = mean_percentage, color = module)) +
  geom_smooth(method = "loess",se = F, linewidth=1) +
  # geom_line() +
  facet_grid(. ~ age_group,scales = "fixed") +
  # geom_stratum(width = 1/12, fill = "black") +
  # scale_fill_manual(values = cols_group) +
  ggtitle('B-cells') + ylab('avg. module score') + 
  theme_bw() +
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_module_dynamics.pdf',sep = ''),width = 10, height = 2.5)
