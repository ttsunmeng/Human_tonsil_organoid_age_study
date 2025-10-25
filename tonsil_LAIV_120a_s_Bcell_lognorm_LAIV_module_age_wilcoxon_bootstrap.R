library(dplyr)
library(Seurat)
library(gridExtra)
library(grid)
library(gtable)
library(dplyr)
library(xlsx)
library(Seurat)
library(ggplot2)
library(ggrepel)
log2FC_cutoff <- 0.25
p_cutoff <- 0.05

Bcell_LAIV_module_wilcoxon <- readRDS('Bcell_gene_module/tonsil_LAIV_120a_s_Bcell_LAIV_lognorm_scale_module_age_wilcoxon_output_parallel_n300.rds')
############# examine output overlaps #######################################
Bcell_LAIV_module_wilcoxon <- Bcell_LAIV_module_wilcoxon[!is.na(Bcell_LAIV_module_wilcoxon$pval),]
Bcell_LAIV_module_wilcoxon <- Bcell_LAIV_module_wilcoxon[Bcell_LAIV_module_wilcoxon$n_min_feature > 1,]
Bcell_LAIV_module_wilcoxon$diff <- Bcell_LAIV_module_wilcoxon$min_mean2 - Bcell_LAIV_module_wilcoxon$min_mean1
Bcell_LAIV_module_wilcoxon$FC_min <- Bcell_LAIV_module_wilcoxon$diff/abs(Bcell_LAIV_module_wilcoxon$min_mean1)#
Bcell_LAIV_module_wilcoxon$FC_std <- Bcell_LAIV_module_wilcoxon$diff/Bcell_LAIV_module_wilcoxon$min_std
# Bcell_LAIV_module_wilcoxon$log2FC <- sign(Bcell_LAIV_module_wilcoxon$FC)*log2(abs(Bcell_LAIV_module_wilcoxon$FC))

Bcell_LAIV_module_wilcoxon <- Bcell_LAIV_module_wilcoxon %>% group_by(age_group1,age_group2) %>% 
  mutate(pval_adj = p.adjust(min_pval,method = 'fdr'))

Bcell_LAIV_module_wilcoxon <- Bcell_LAIV_module_wilcoxon %>% group_by(day_group,module,age_group1,age_group2) %>% 
  mutate(pval_adj_mean = exp(mean(log(pval_adj),na.rm = T)),diff_mean = mean(diff,na.rm = T), FC_min_mean = mean(FC_min,na.rm = T),FC_std_mean = mean(FC_std,na.rm = T))

Bcell_LAIV_module_wilcoxon$effect <- '+'
Bcell_LAIV_module_wilcoxon$effect[Bcell_LAIV_module_wilcoxon$diff_mean < 0] <- '-'
Bcell_LAIV_module_wilcoxon$significant <- 'no'
Bcell_LAIV_module_wilcoxon$significant[(Bcell_LAIV_module_wilcoxon$pval_adj <= p_cutoff)] <- 'yes'

Bcell_LAIV_module_wilcoxon <- Bcell_LAIV_module_wilcoxon %>% group_by(day_group,module,age_group1,age_group2,effect) %>% 
  mutate(n_significant = sum(significant == 'yes'))
n_bootstrap_final <- max(Bcell_LAIV_module_wilcoxon$bootstrap_index)
Bcell_LAIV_module_wilcoxon <- Bcell_LAIV_module_wilcoxon %>% group_by(day_group,module,age_group1,age_group2,effect) %>% 
  mutate(significant = ifelse((n_significant > floor(n_bootstrap_final*0.95)), "yes", "no"))#
Bcell_LAIV_module_wilcoxon$label <- as.character(Bcell_LAIV_module_wilcoxon$module)
Bcell_LAIV_module_wilcoxon$label[Bcell_LAIV_module_wilcoxon$significant == 'no'] <- ''
Bcell_LAIV_module_wilcoxon <- Bcell_LAIV_module_wilcoxon %>% group_by(day_group,module,age_group1,age_group2) %>% 
  mutate(effect = ifelse(diff_mean > 0, "+", "-"))#

Bcell_LAIV_module_wilcoxon$days <- Bcell_LAIV_module_wilcoxon$day_group
Bcell_LAIV_module_wilcoxon$day_group <- paste('day',sprintf('%#02s',Bcell_LAIV_module_wilcoxon$day_group),sep = '')
Bcell_LAIV_module_wilcoxon$day_group[Bcell_LAIV_module_wilcoxon$day_group == 'day00'] <- 'day0'
LAIV_day_group_list <- as.character(unique(Bcell_LAIV_module_wilcoxon$day_group))

Bcell_LAIV_module_wilcoxon_youth_vs_adult <- Bcell_LAIV_module_wilcoxon[(Bcell_LAIV_module_wilcoxon$age_group1 == '02-03yrs') & (Bcell_LAIV_module_wilcoxon$age_group2 == '27-39yrs'),]
Bcell_LAIV_module_wilcoxon_youth_vs_adult <- unique(Bcell_LAIV_module_wilcoxon_youth_vs_adult[,c('day_group','days','module','FC_std_mean','FC_min_mean','pval_adj_mean','n_min_feature','significant','effect')])
Bcell_LAIV_module_wilcoxon_youth_vs_adult <- Bcell_LAIV_module_wilcoxon_youth_vs_adult %>% arrange(desc(significant),desc(FC_std_mean)) 
Bcell_LAIV_module_wilcoxon_youth_vs_adult <- data.frame(Bcell_LAIV_module_wilcoxon_youth_vs_adult)
Bcell_LAIV_module_wilcoxon_youth_vs_adult$sig_effect <- 0
Bcell_LAIV_module_wilcoxon_youth_vs_adult$sig_effect[(Bcell_LAIV_module_wilcoxon_youth_vs_adult$effect == '+') & 
                                                       (Bcell_LAIV_module_wilcoxon_youth_vs_adult$significant == 'yes')] <- 1
Bcell_LAIV_module_wilcoxon_youth_vs_adult$sig_effect[(Bcell_LAIV_module_wilcoxon_youth_vs_adult$effect == '-') & 
                                                       (Bcell_LAIV_module_wilcoxon_youth_vs_adult$significant == 'yes')] <- -1                                                      
write.xlsx(Bcell_LAIV_module_wilcoxon_youth_vs_adult,
           'tonsil_LAIV_120a_s_Bcells_LAIV_lognorm_scale_module_wilcoxon_bootstrap_n300_toddler_vs_adult.xlsx')
# for (day_name in LAIV_day_group_list){
#   temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult[Bcell_LAIV_module_wilcoxon_youth_vs_adult$day_group == day_name,]
#   write.xlsx(temp,sheetName = day_name,
#              'tonsil_LAIV_120a_s_Bcells_LAIV_lognorm_scale_module_wilcoxon_bootstrap_n300_toddler_vs_adult.xlsx',append = T)
#   
# }

temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult
temp$significant[(abs(temp$FC_mean) < 0.50) & (temp$significant == 'yes')] <- 'no'
temp$significant[(abs(temp$FC_mean) >= 0.50) & (temp$significant == 'yes')] <- 'yes'

temp$label <- as.character(temp$module)
temp$label[temp$significant == 'no'] <- ''
temp$label[temp$label != ''] <- paste(temp$day_group[temp$label != ''],temp$label[temp$label != ''],sep = '_')

ggplot(temp, aes(x=FC_mean, y=-log10(pval_adj_mean),color = significant)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_colour_manual(values = c("yes" = "red","no" = "darkgrey")) +
  ggtitle(paste('DEG of Bcells not presented at baseline\nhigher in adults higher in toddlers\n')) +
  xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 10)) +
  theme_bw() +
  geom_text_repel(aes(x=FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
  geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n300_toddler_vs_adult_big.pdf',sep = ''),width = 15, height = 15)
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n300_toddler_vs_adult.pdf',sep = ''),width = 6.5, height = 5)


temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult
temp$log2FC_mean <- -temp$log2FC_mean
temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'no'
temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'only at baseline'
# need to be same direction
temp$significant[(temp$marker %in% temp$marker[(temp$significant == 'only at baseline') & (temp$day_group != 'day0') & (temp$log2FC_mean > 0)]) & temp$log2FC_mean > 0] <- 'persist'
temp$significant[(temp$marker %in% temp$marker[(temp$significant == 'only at baseline') & (temp$day_group != 'day0') & (temp$log2FC_mean < 0)]) & temp$log2FC_mean < 0] <- 'persist'
temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'no'
temp$label <- as.character(temp$label)
# temp$color <- 'black'
# temp$color[temp$significant == 'only at baseline'] <- 'grey'
temp$label[temp$significant == 'no'] <- ''
# temp$label[temp$marker %in% temp$marker[(temp$significant == 'yes') & (temp$day_group == 'day0')]] <- ''
temp <- temp[!temp$label %in% c('ab-CD3','ab-Bcell','ab-TCR-gamma-delta','ab-CD8','ab-TCR-alpha-beta'),] 
temp <- temp[temp$day_group == 'day0',]
ggplot(temp, aes(x=log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_colour_manual(values = c("persist" = "red","only at baseline" = "orange","no" = "darkgrey")) +
  ggtitle(paste('DEG of Bcells not presented at baseline\nhigher in adults higher in toddlers\n')) +
  xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 10)) +
  theme_bw() +
  geom_text_repel(aes(x=log2FC_mean, y=-log10(pval_adj_mean),label = label), color = 'black', size = 3) +
  geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n300_toddler_vs_adult_day0.pdf',sep = ''),width = 6.5, height = 5)


temp <- read.xlsx('./Bcell_lognorm_CCA/tonsil_LAIV_120a_s_Bcell_BCL6pos_lognorm_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'Sheet1')
temp <- temp[temp$log2FC_mean < -1,]
####### 
# select all cytokines
cytokine_global_table <- read.xlsx('../Global landscape of cytokines Supplementary Table S1.xlsx', sheetName = "Cytokines")
cytokine_global_list <- unique(cytokine_global_table$HGNC.symbol)
cytokine_global_list <- c(cytokine_global_list,'GZMB','PRF1','MIF','FASLG','VEGFA','TGFA')
seurat_cytokine_list <- all.markers[all.markers %in% cytokine_global_list]

temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult[!Bcell_LAIV_module_wilcoxon_youth_vs_adult$day_group %in% c('day07-08'),]
temp$log2FC_mean <- -temp$log2FC_mean
temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'no'
temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'yes'
# as.character(temp$day_group[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')])
temp <- temp[temp$marker %in% seurat_cytokine_list,]
temp$label <- as.character(temp$label)
temp$label[temp$significant == 'no'] <- ''
temp$label[temp$label != ''] <- paste(temp$day_group[temp$label != ''],temp$label[temp$label != ''],sep = '_')

ggplot(temp, aes(x=log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_colour_manual(values = c("yes" = "red","no" = "darkgrey")) +
  ggtitle(paste('DEG of Bcells not presented at baseline\nhigher in adults higher in toddlers\n')) +
  xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 10)) +
  theme_bw() +
  geom_text_repel(aes(x=log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
  geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n300_toddler_vs_adult_big.pdf',sep = ''),width = 15, height = 15)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n300_toddler_vs_adult_cytokine_only.pdf',sep = ''),width = 10, height = 10)

all_Bcell_LAIV_module_wilcoxon_youth_vs_adult <- Bcell_LAIV_module_wilcoxon_youth_vs_adult
NaiveOrActB_LAIV_wilcoxon_youth_vs_adult <- Bcell_LAIV_module_wilcoxon_youth_vs_adult
PB_LAIV_wilcoxon_youth_vs_adult <- Bcell_LAIV_module_wilcoxon_youth_vs_adult
GClikeB_LAIV_wilcoxon_youth_vs_adult <- Bcell_LAIV_module_wilcoxon_youth_vs_adult
MemoryB_LAIV_wilcoxon_youth_vs_adult <- Bcell_LAIV_module_wilcoxon_youth_vs_adult
preGCOrGCB_LAIV_wilcoxon_youth_vs_adult <- Bcell_LAIV_module_wilcoxon_youth_vs_adult

gene_name <- 'CXCL10'
gene_name %in% NaiveOrActB_LAIV_wilcoxon_youth_vs_adult$marker[(NaiveOrActB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                                 (NaiveOrActB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                                 (NaiveOrActB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')]
gene_name %in% preGCOrGCB_LAIV_wilcoxon_youth_vs_adult$marker[(preGCOrGCB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                                (preGCOrGCB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                                (preGCOrGCB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')]
gene_name %in% GClikeB_LAIV_wilcoxon_youth_vs_adult$marker[(GClikeB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                             (GClikeB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                             (GClikeB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')]
gene_name %in% MemoryB_LAIV_wilcoxon_youth_vs_adult$marker[(MemoryB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                             (MemoryB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                             (MemoryB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')]
gene_name %in% PB_LAIV_wilcoxon_youth_vs_adult$marker[(PB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                        (PB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                        (PB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')]

temp <- unique(as.character(c(NaiveOrActB_LAIV_wilcoxon_youth_vs_adult$marker[(NaiveOrActB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                            (NaiveOrActB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                            (NaiveOrActB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')],
          preGCOrGCB_LAIV_wilcoxon_youth_vs_adult$marker[(preGCOrGCB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                           (preGCOrGCB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                           (preGCOrGCB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')],
          GClikeB_LAIV_wilcoxon_youth_vs_adult$marker[(GClikeB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                        (GClikeB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                        (GClikeB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')],
          MemoryB_LAIV_wilcoxon_youth_vs_adult$marker[(MemoryB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                        (MemoryB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                        (MemoryB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')],
          PB_LAIV_wilcoxon_youth_vs_adult$marker[(PB_LAIV_wilcoxon_youth_vs_adult$day_group != 'day0') & 
                                                   (PB_LAIV_wilcoxon_youth_vs_adult$log2FC_mean < 0) & 
                                                   (PB_LAIV_wilcoxon_youth_vs_adult$significant == 'yes')])))
########### plotting heatmap selecting the highest logFC #########################
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays[(Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays$effect == '+'),] 
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays[(Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays$effect == '-'),] 

Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation %>% group_by(marker) %>% 
  summarize(log2FCmax = -max(log2FC_mean),day_group = day_group[which.max(log2FC_mean)],if_significant = significant[which.max(log2FC_mean)])
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation %>% group_by(marker) %>% 
  summarize(log2FCmax = max(-log2FC_mean),day_group = day_group[which.max(-log2FC_mean)],if_significant = significant[which.max(-log2FC_mean)])

Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation_sig <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays[(Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays$significant == 'yes') & 
                                               (Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays$effect == '+'),] 
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation_sig <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays[(Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays$significant == 'yes') & 
                                               (Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays$effect == '-'),] 
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation_sig <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation_sig %>% group_by(marker) %>% 
  summarize(log2FCmax = -max(log2FC_mean),day_group = day_group[which.max(log2FC_mean)])
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation_sig <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation_sig %>% group_by(marker) %>% 
  summarize(log2FCmax = max(-log2FC_mean),day_group = day_group[which.max(-log2FC_mean)])

Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation$celltype <- 'BCL6-nonmemory B'
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation$celltype <- 'BCL6-nonmemory B'
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation_sig$celltype <- 'BCL6-nonmemory B'
Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation_sig$celltype <- 'BCL6-nonmemory B'

Bcell_subset_DEG_downregulation <- rbind(Bcell_subset_DEG_downregulation,Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation)
Bcell_subset_DEG_upregulation <- rbind(Bcell_subset_DEG_upregulation,Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation)
Bcell_subset_DEG_downregulation_sig <- rbind(Bcell_subset_DEG_downregulation_sig,Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_upregulation_sig)
Bcell_subset_DEG_upregulation_sig <- rbind(Bcell_subset_DEG_upregulation_sig,Bcell_LAIV_module_wilcoxon_youth_vs_adult_alldays_downregulation_sig)

saveRDS(Bcell_subset_DEG_downregulation,'Bcell_lognorm_DEG_age_downregulation.rds')
saveRDS(Bcell_subset_DEG_upregulation,'Bcell_lognorm_DEG_age_upregulation.rds')
saveRDS(Bcell_subset_DEG_downregulation_sig,'Bcell_lognorm_DEG_age_downregulation_sig.rds')
saveRDS(Bcell_subset_DEG_upregulation_sig,'Bcell_lognorm_DEG_age_upregulation_sig.rds')

Bcell_subset_DEG_noday0_downregulation_sig <- Bcell_subset_DEG_downregulation_sig[Bcell_subset_DEG_downregulation_sig$day_group != 'day0',]
Bcell_subset_DEG_noday0_downregulation_sig <- Bcell_subset_DEG_noday0_downregulation_sig %>% group_by(marker) %>%
  mutate(count = n())
Bcell_subset_DEG_noday0_upregulation_sig <- Bcell_subset_DEG_upregulation_sig[Bcell_subset_DEG_upregulation_sig$day_group != 'day0',]
Bcell_subset_DEG_noday0_upregulation_sig <- Bcell_subset_DEG_noday0_upregulation_sig %>% group_by(marker) %>%
  mutate(count = n()) 


length(unique(Bcell_subset_DEG_upregulation_sig$marker[Bcell_subset_DEG_upregulation_sig$count == 3]))
########## heatmap of same markers across different population ################
library(reshape2)
library(ComplexHeatmap)
temp <- Bcell_subset_DEG_noday0_upregulation_sig[Bcell_subset_DEG_noday0_upregulation_sig$count == 3,]
temp$log2FCmax <- as.numeric(temp$log2FCmax)

temp <- dcast(temp, marker ~ celltype, value.var="log2FCmax")
rownames(temp) <- temp$marker
temp$marker <- NULL
temp <- temp %>% arrange(desc(rowMeans(temp[,c('BCL6-nonmemory B','BCL6+B')])-temp$`memory B`))

# temp <- temp %>% arrange(desc(rowMeans(temp)))
temp <- temp[rowSums(temp > 0.25) == 3,]

Heatmap(as.matrix(temp), cluster_columns = F, cluster_rows = F, 
                name = "log2FC", show_column_names = T,#
                show_row_names = T,column_title = 'Bcell', #row_km = 3, #column_km = 20, 
                column_title_side = "top",column_names_side = "top",row_title_side = "left",row_names_side = 'left')
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_celltype_logFCmax_upregulation.pdf',sep = ''),width = 5, height = 6)



############ ClusterProfiler ############################
# library("org.Hs.eg.db")
# hs <- org.Hs.eg.db
library(enrichplot)
library(clusterProfiler)
for (day_group_name in LAIV_day_group_list){
  print(day_group_name)
  temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult_upregulation[Bcell_LAIV_module_wilcoxon_youth_vs_adult_upregulation$day_group == day_group_name,]
  
  ego <- enrichGO(temp$ENTREZID[!is.na(temp$ENTREZID)], OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.05, 
                  universe = all_marker_table$ENTREZID, ont="BP", readable=TRUE)
  
  if (dim(ego)[1] != 0) {
    print(paste(day_group_name,'upregulation'))
    write.xlsx(ego@result,paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n316_toddler_vs_adult_GO_ClusterProfiler.xlsx',sep = ''),
               sheetName = paste(day_group_name,'upregulation'),append = T)
  }
  temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult_downregulation[Bcell_LAIV_module_wilcoxon_youth_vs_adult_downregulation$day_group == day_group_name,]
  
  ego <- enrichGO(temp$ENTREZID[!is.na(temp$ENTREZID)], OrgDb = "org.Hs.eg.db", pvalueCutoff = 0.05, 
                  universe = all_marker_table$ENTREZID, ont="BP", readable=TRUE)
  
  if (dim(ego)[1] != 0) {
    print(paste(day_group_name,'downregulation'))
    write.xlsx(ego@result,paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n316_toddler_vs_adult_GO_ClusterProfiler.xlsx',sep = ''),
               sheetName = paste(day_group_name,'downregulation'),append = T)
    graphics.off()
    plot <- dotplot(ego) + ggtitle(paste('GO of Bcell LAIV vs ns\ndownregulated on',day_name)) +
      theme(axis.text = element_text(size = 20),plot.title = element_text(size = 25, face = "bold"))
    print(plot)
    dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n316_toddler_vs_adult_GO_ClusterProfiler_downregulation_',day_name,'.pdf',sep = ''),width = 8, height = 10)
  }
}
   
# network visualization
bp <- pairwise_termsim(ego)
plot <- emapplot(bp,showCategory = 200)
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n316_toddler_vs_adult_GO_ClusterProfiler_upregulation_day07-08.pdf',sep = ''),
          width = 10, height = 10)
# simplify
bp <- simplify(ego, by="p.adjust", select_fun=min)
write.xlsx(bp@result,paste('tonsil_LAIV_120a_s_Bcell_DE_NS_vs_day0_',donor_name,'_',day_name,'_log1p_paired_GO_ClusterProfiler_downregulation_simplify.xlsx',sep = ''))
plot <- dotplot(bp, showCategory=20) + ggtitle(paste('GO of Bcell LAIV vs ns\ndownregulated on',day_name)) +
  theme(axis.text = element_text(size = 20),plot.title = element_text(size = 25, face = "bold"))
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_DE_NS_vs_day0_',donor_name,'_',day_name,'_log1p_paired_GO_ClusterProfiler_downregulation_simplify.pdf',sep = ''),width = 8, height = 10)
bp <- pairwise_termsim(bp)
plot <- emapplot(bp)
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_DE_NS_vs_day0_',donor_name,'_',day_name,'_log1p_paired_GO_ClusterProfiler_downregulation_network_simplify.pdf',sep = ''),
          width = 7, height = 7)
  
  
table1 <- 'tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_toddler_vs_adult_downregulation_day04-06_DAVID_UP_KW_BIOLOGICAL_PROCESS.txt'
table1_DE_data <- read.table(table1, sep = '\t',header = T)
day_group_name <- 'day04-06'
temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult_downregulation %>% dplyr::filter(day_group == day_group_name)
temp <- temp[,c('marker','ENTREZID')]
temp$ENTREZID <- as.numeric(temp$ENTREZID)
rownames(temp) <- temp$ENTREZID
temp <- temp %>% arrange(desc(ENTREZID))
# for (geneID_index in rownames(RSTR_signficant_genes)) {
#   table1_DE_data$Genes <- gsub(paste(', ',geneID_index,',',sep = ''),paste(', ',RSTR_signficant_genes[geneID_index,][,c('gene')],',',sep = ''),table1_DE_data$Genes)
#   table1_DE_data$Genes <- gsub(paste('^',geneID_index,',',sep = ''),paste(RSTR_signficant_genes[geneID_index,][,c('gene')],',',sep = ''),table1_DE_data$Genes)
#   table1_DE_data$Genes <- gsub(paste(', ',geneID_index,'$',sep = ''),paste(', ',RSTR_signficant_genes[geneID_index,][,c('gene')],sep = ''),table1_DE_data$Genes)
# }
# write.xlsx(table1_DE_data,sheetName = 'Sheet2',table1,append = T)

table1_DE_data$Term <- gsub('.*\\~','',table1_DE_data$Term)
table1_DE_data <- table1_DE_data %>% dplyr::filter(FDR <= 0.05)
# BP_direct_select <- c('translation','mitochondrial electron transport, NADH to ubiquinone',
#                       'cell-cell adhesion','movement of cell or suBcellular component','Wnt signaling pathway, planar cell polarity pathway')
# 
# BP_direct_select <- c('NIK/NF-kappaB signaling','T cell receptor signaling pathway','MAPK cascade',
#                       'tumor necrosis factor-mediated signaling pathway')
# table1_DE_data <- table1_DE_data[table1_DE_data$Term %in% BP_direct_select,]
library(forcats)
library(scales)

table1_DE_data %>%
  mutate(Term = fct_reorder(Term, Count)) %>%
  ggplot( aes(x=Term, y=Count, fill = FDR)) +
  geom_bar(stat="identity", colour="black") +
  coord_flip() +
  # ggtitle('metabolic and basic activities') +
  ggtitle(paste('Bcell expression diff by LAIV on',day_group_name,'\nlower in toddlers')) +
  xlab(expression(paste('UniprotKB Keywords\nBiological Processes'))) + ylab('#genes') + 
  theme_bw() +
  scale_fill_gradientn(colors = c('yellow','#984EA3'),limits = c(0,0.05)) +
  scale_x_discrete(labels = wrap_format(35)) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10),
        plot.title = element_text(size = 12, face = "bold"))
# theme(axis.text = element_text(size = 20),plot.title = element_text(size = 25, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_toddler_vs_adult_downregulation_day04-06_UP_KW_BP.pdf',sep = ''),width = 4.5, height = 2.5)



############# two age group comparison ##################
age_group1 <- '02-03yrs'
age_group2 <- '07-09yrs'
temp <- Bcell_LAIV_module_wilcoxon[(Bcell_LAIV_module_wilcoxon$age_group1 == age_group1),]
temp <- temp[(temp$age_group2 == age_group2),]

temp_upregulation <- temp[(temp$significant == 'yes') & (temp$effect == '+'),]
temp_upregulation <-
  unique(temp_upregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
temp_upregulation$days <- 
  substring(temp_upregulation$day_group,first = 4,last = 30)
temp_upregulation <- temp_upregulation %>% arrange(days,desc(log2FC_mean)) 
temp_upregulation <- temp_upregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(temp_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_sBcell_lognorm_LAIV_bootstrap_toddlers_vs_tweens.xlsx')

temp_downregulation <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
temp_downregulation <-
  unique(temp_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
temp_downregulation$days <- 
  substring(temp_downregulation$day_group,first = 4,last = 30)
temp_downregulation <- temp_downregulation %>% arrange(days,log2FC_mean) 
temp_downregulation <- temp_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(temp_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_sBcell_lognorm_LAIV_bootstrap_toddlers_vs_tweens.xlsx',append = T)

# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp1 <- unique(temp[(temp$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp1, aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of Bcell expression change by LAIV\nbetween toddlers vs tweens on',day_group_name)) +
#     xlab("log2 fold-change of expression diff") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-FC_cutoff,FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_toddlers_vs_tweens_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }

age_group1 <- '07-09yrs'
age_group2 <- '27-39yrs'
temp <- Bcell_LAIV_module_wilcoxon[(Bcell_LAIV_module_wilcoxon$age_group1 == age_group1),]
temp <- temp[(temp$age_group2 == age_group2),]

temp_upregulation <- temp[(temp$significant == 'yes') & (temp$effect == '+'),]
temp_upregulation <-
  unique(temp_upregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
temp_upregulation$days <- 
  substring(temp_upregulation$day_group,first = 4,last = 30)
temp_upregulation <- temp_upregulation %>% arrange(days,desc(log2FC_mean)) 
temp_upregulation <- temp_upregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(temp_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_sBcell_lognorm_LAIV_bootstrap_tweens_vs_adults.xlsx')

temp_downregulation <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
temp_downregulation <-
  unique(temp_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
temp_downregulation$days <- 
  substring(temp_downregulation$day_group,first = 4,last = 30)
temp_downregulation <- temp_downregulation %>% arrange(days,log2FC_mean) 
temp_downregulation <- temp_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(temp_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_sBcell_lognorm_LAIV_bootstrap_tweens_vs_adults.xlsx',append = T)

# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp1 <- unique(temp[(temp$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp1, aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of Bcell expression change by LAIV\nbetween tweens vs adults on',day_group_name)) +
#     xlab("log2 fold-change of expression diff") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-FC_cutoff,FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_tweens_vs_adults_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }
############# stringent in toddler vs both ########################
Bcell_LAIV_module_wilcoxon_youth_stringent <- Bcell_LAIV_module_wilcoxon_youth %>% group_by(day_group,marker) %>% 
  mutate(significant = ifelse(all(significant == 'yes') & (prod(log2FC_mean) > 0),'yes','no'))
        
        
# temp <- unique(Bcell_LAIV_module_wilcoxon_youth_stringent[,c('day_group','marker','age_group1','age_group2','significant','significant2','log2FC_mean','effect')])
# temp <- temp %>% arrange(day_group,marker)

Bcell_LAIV_module_wilcoxon_youth_stringent <- Bcell_LAIV_module_wilcoxon_youth_stringent %>% group_by(day_group,marker) %>% 
  mutate(pval_adj_mean = exp(mean(log(unique(pval_adj_mean)))),log2FC_mean = mean(unique(log2FC_mean)),n2 = mean(unique(n2)))
Bcell_LAIV_module_wilcoxon_youth_stringent <- Bcell_LAIV_module_wilcoxon_youth_stringent %>% group_by(day_group,marker) %>% 
  mutate(effect = ifelse(log2FC_mean > 0, "+", "-"))
Bcell_LAIV_module_wilcoxon_youth_stringent$label[Bcell_LAIV_module_wilcoxon_youth_stringent$significant == 'no'] <- ''

Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <- 
  Bcell_LAIV_module_wilcoxon_youth_stringent[(Bcell_LAIV_module_wilcoxon_youth_stringent$significant == 'yes') & 
                                       (Bcell_LAIV_module_wilcoxon_youth_stringent$effect == '+'),]
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <-
  unique(Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','n1','n2','ENTREZID')])
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation$days <- 
  substring(Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation$day_group,first = 4,last = 30)
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <- Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation %>% arrange(days,desc(log2FC_mean)) 
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <- Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_sBcell_lognorm_LAIV_bootstrap_stringent_toddler_vs_both.xlsx')

Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <- 
  Bcell_LAIV_module_wilcoxon_youth_stringent[(Bcell_LAIV_module_wilcoxon_youth_stringent$significant == 'yes') & 
                                       (Bcell_LAIV_module_wilcoxon_youth_stringent$effect == '-'),]
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <-
  unique(Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','n1','n2','ENTREZID')])
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation$days <- 
  substring(Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation$day_group,first = 4,last = 30)
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <- Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation %>% arrange(days,log2FC_mean) 
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <- Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_sBcell_lognorm_LAIV_bootstrap_stringent_toddler_vs_both.xlsx',append = T)

# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp <- unique(Bcell_LAIV_module_wilcoxon_youth_stringent[(Bcell_LAIV_module_wilcoxon_youth_stringent$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$log2_log2FC_mean <- sign(temp$log2FC_mean)*log2(abs(temp$log2FC_mean) + 1)
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp, aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of Bcell expression change by LAIV\nbetween toddlers and olders on',day_group_name)) +
#     xlab("log2 fold-change of expression diff") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-FC_cutoff,FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_stringent_toddler_vs_both_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }

############# stringent in adult vs both ########################
Bcell_LAIV_module_wilcoxon_adult <- Bcell_LAIV_module_wilcoxon[(Bcell_LAIV_module_wilcoxon$age_group2 == '27-39yrs'),]
Bcell_LAIV_module_wilcoxon_adult_stringent <- Bcell_LAIV_module_wilcoxon_adult %>% group_by(day_group,marker) %>% 
  mutate(significant = ifelse(all(significant == 'yes') & (prod(unique(log2FC_mean)) > 0),'yes','no'))
Bcell_LAIV_module_wilcoxon_adult_stringent <- Bcell_LAIV_module_wilcoxon_adult_stringent %>% group_by(day_group,marker) %>% 
  mutate(pval_adj_mean = exp(mean(log(unique(pval_adj_mean)))),log2FC_mean = mean(unique(log2FC_mean)),n1 = mean(unique(n1)))
Bcell_LAIV_module_wilcoxon_adult_stringent <- Bcell_LAIV_module_wilcoxon_adult_stringent %>% group_by(day_group,marker) %>% 
  mutate(effect = ifelse(log2FC_mean > 0, "+", "-"))
Bcell_LAIV_module_wilcoxon_adult_stringent$label[Bcell_LAIV_module_wilcoxon_adult_stringent$significant == 'no'] <- ''

Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation <- 
  Bcell_LAIV_module_wilcoxon_adult_stringent[(Bcell_LAIV_module_wilcoxon_adult_stringent$significant == 'yes') & 
                                        (Bcell_LAIV_module_wilcoxon_adult_stringent$effect == '+'),]
Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation <-
  unique(Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','n1','n2','ENTREZID')])
Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation$days <- 
  substring(Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation$day_group,first = 4,last = 30)
Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation <- Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation %>% arrange(days,desc(log2FC_mean)) 
Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation <- Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(Bcell_LAIV_module_wilcoxon_adult_stringent_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_sBcell_lognorm_LAIV_bootstrap_stringent_adult_vs_both.xlsx')

Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation <- 
  Bcell_LAIV_module_wilcoxon_adult_stringent[(Bcell_LAIV_module_wilcoxon_adult_stringent$significant == 'yes') & 
                                        (Bcell_LAIV_module_wilcoxon_adult_stringent$effect == '-'),]
Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation <-
  unique(Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','n1','n2','ENTREZID')])
Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation$days <- 
  substring(Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation$day_group,first = 4,last = 30)
Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation <- Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation %>% arrange(days,log2FC_mean) 
Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation <- Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(Bcell_LAIV_module_wilcoxon_adult_stringent_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_sBcell_lognorm_LAIV_bootstrap_stringent_adult_vs_both.xlsx',append = T)

# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp <- unique(Bcell_LAIV_module_wilcoxon_adult_stringent[(Bcell_LAIV_module_wilcoxon_adult_stringent$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$log2_log2FC_mean <- sign(temp$log2FC_mean)*log2(abs(temp$log2FC_mean) + 1)
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp, aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of Bcell expression change by LAIV\nbetween adults and olders on',day_group_name)) +
#     xlab("log2 fold-change of expression diff") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-FC_cutoff,FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_stringent_adult_vs_both_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }
# 
############## venn diagram #########################
# compare BCL6+, BCL6-memory and BCL6- nonmemory nonPB 
temp1 <- read.xlsx('./Bcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Bcell_BCL6neg_memory_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'downregulation')
temp2 <- read.xlsx('./Bcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Bcell_BCL6pos_nonPB_n238_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'downregulation')
temp3 <- read.xlsx('./Bcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Bcell_BCL6pos_lognorm_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'downregulation')
temp <- read.xlsx('./Bcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n316_toddler_vs_adult.xlsx',sheetName = 'downregulation')
temp11 <- unique(temp1$marker[temp1$day_group != 'day0'])
Bcell_GC_downregulated_marker_list <- unique(temp2$marker[temp2$day_group != 'day0'])
temp33 <- unique(temp3$marker[temp3$day_group != 'day0'])
Bcell_downregulated_marker_list <- unique(temp$marker[temp$day_group != 'day0'])
# all_Bcell_sig_genes <- Bcell_downregulated_marker_list
# all_Bcell_T17_sig_genes <- temp11
# all_Bcell_Tfh_sig_genes <- temp22
# all_Bcell_non_Tfh_sig_genes <- temp33
# all_Bcell_Tfh_like_sig_genes <- temp44
overlapped_downregulation_marker_list <- Bcell_downregulated_marker_list[Bcell_downregulated_marker_list %in% unique(temp11[(temp11 %in% Bcell_GC_downregulated_marker_list) & (temp11 %in% temp33)])]

# length(unique(temp11[(temp11 %in% temp22) & (temp11 %in% temp33)]))
# length(unique(temp11[(temp11 %in% temp22) & (!temp11 %in% temp33)]))
# length(unique(temp11[(!temp11 %in% temp22) & (temp11 %in% temp33)]))
# length(unique(temp11[(!temp11 %in% temp22) & (!temp11 %in% temp33)]))
# length(unique(temp22[(!temp22 %in% temp11) & (temp22 %in% temp33)]))
# length(unique(temp22[(!temp22 %in% temp11) & (!temp22 %in% temp33)]))
# length(unique(temp33[(!temp33 %in% temp11) & (!temp33 %in% temp22)]))

temp1 <- read.xlsx('./Bcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Bcell_BCL6neg_memory_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'upregulation')
temp2 <- read.xlsx('./Bcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Bcell_BCL6pos_nonPB_n238_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'upregulation')
temp3 <- read.xlsx('./Bcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Bcell_BCL6pos_lognorm_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'upregulation')
temp <- read.xlsx('./Bcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_bootstrap_n316_toddler_vs_adult.xlsx',sheetName = 'upregulation')
temp11 <- unique(temp1$marker[temp1$day_group != 'day0'])
Bcell_GC_upregulated_marker_list <- unique(temp2$marker[temp2$day_group != 'day0'])
temp33 <- unique(temp3$marker[temp3$day_group != 'day0'])
Bcell_upregulated_marker_list <- unique(temp$marker[temp$day_group != 'day0'])
# all_Bcell_sig_genes <- Bcell_upregulated_marker_list
# all_Bcell_T17_sig_genes <- temp11
# all_Bcell_Tfh_sig_genes <- temp22
# all_Bcell_non_Tfh_sig_genes <- temp33
# all_Bcell_Tfh_like_sig_genes <- temp44
overlapped_upregulation_marker_list <- Bcell_upregulated_marker_list[Bcell_upregulated_marker_list %in% unique(temp11[(temp11 %in% temp22) & (temp11 %in% temp33)])]

marker_name <- 'IGHA1-secreted'
marker_name %in% temp11
marker_name %in% temp22
marker_name %in% temp33

temp_overlap <- unique(temp11[(temp11 %in% temp22) & (temp11 %in% temp33)])

cytokine_global_table <- read.xlsx('../Global landscape of cytokines Supplementary Table S1.xlsx', sheetName = "Cytokines")
cytokine_global_list <- unique(cytokine_global_table$HGNC.symbol)
cytokine_global_list <- c(cytokine_global_list,'GZMB','PRF1','MIF')
all_Bcell_sig_genes <- unique(all_Bcell_sig_genes)
all_Bcell_Tfh_sig_genes <- unique(all_Bcell_sig_genes)
seurat_cytokine_list <- cytokine_global_list[cytokine_global_list %in% all_Bcell_sig_genes]
seurat_cytokine_list <- all.markers[all.markers %in% cytokine_global_list]
############## visualization #################################
library(reshape2)
Bcell_GC_overlapped_upregulation_downregulation_marker_list <- Bcell_GC_downregulated_marker_list[Bcell_GC_downregulated_marker_list %in% Bcell_GC_upregulated_marker_list]
# temp <- temp[!(temp$marker %in% c('ab-TCR-gamma-delta','ab-TCR-alpha-beta','ab-CD3','ab-Bcell_GC',"ab-CD8")),]
# temp <- temp[!(temp$marker %in% c('CD38',"CXCR3","IGHG1-membrane","IL2RA","ENTPD1",'GPI',"IL2RB","PFKP","TPI1","LDHA","PGK1","ENO1")),]

Bcell_GC_LAIV_wilcoxon_youth_vs_adult_upregulation_cytokine <- 
  Bcell_GC_LAIV_wilcoxon_youth_vs_adult_upregulation[Bcell_GC_LAIV_wilcoxon_youth_vs_adult_upregulation$marker %in% Bcell_GC_overlapped_upregulation_downregulation_marker_list,]
Bcell_GC_LAIV_wilcoxon_youth_vs_adult_downregulation_cytokine <- 
  Bcell_GC_LAIV_wilcoxon_youth_vs_adult_downregulation[Bcell_GC_LAIV_wilcoxon_youth_vs_adult_downregulation$marker %in% Bcell_GC_overlapped_upregulation_downregulation_marker_list,]
# combine the upregulation and downregulation 
temp <- Bcell_GC_LAIV_wilcoxon_youth_vs_adult_upregulation_cytokine[,c('marker','days','log2FC_mean')]
temp$day_group <- paste('day',temp$days,sep = '')
temp$effect <- 'downregulation'
temp2 <- data.frame(table(temp$marker,temp$days))
colnames(temp2) <- c('marker','days','Freq')
temp2$Freq[temp2$Freq > 0] <- 1
temp2$effect <- 'up'

temp4 <- Bcell_GC_LAIV_wilcoxon_youth_vs_adult_downregulation_cytokine[,c('marker','days','log2FC_mean')]
temp4$day_group <- paste('day',temp4$days,sep = '')
temp4$effect <- 'upregulation'
temp3 <- data.frame(table(temp4$marker,temp4$days))
colnames(temp3) <- c('marker','days','Freq')
temp3$Freq[temp3$Freq > 0] <- 1
temp3$effect <- 'down'
temp2 <- rbind(temp2,temp3)

temp1 <- dcast(temp2, marker ~ days + effect, value.var="Freq")
temp1[is.na(temp1)] <- 0
temp1 <- data.frame(temp1)
temp1 <- temp1 %>% arrange(X0_down,X04_down,X12.14_down,X0_up,X06.07_up,X10_up,X12.14_up)
# temp1 <- temp1 %>% arrange(X0_down,X04_down,X06.07_down,X10_down,X12.14_down,X0_up,X04_up,X06.07_up,X10_up,X12.14_up)

# Hierarchical clustering using Complete Linkage
# hc1 <- hclust(dist(temp1, method = "euclidean"), method = "complete")
temp_ordered_marker <- temp1$marker #as.character(temp1$marker[rev(hc1$order)])

temp2 <- rbind(temp,temp4)
temp2$marker <- factor(temp2$marker, levels = temp_ordered_marker)
temp2 <- temp2[temp2$day_group != 'day07-08',]
# # temp2 <- temp2[!(temp2$marker %in% c('ab-TCR-gamma-delta','ab-TCR-alpha-beta','ab-CD3','ab-Bcell_GC','CD3E')),]
# ggplot(temp2, aes(x=days, y=marker,color = log2FC_mean)) + geom_point() + #(1 + )
#   RotatedAxis() + 
#   facet_wrap(~effect) +
#   # scale_colour_manual(values = c("upregulation" = "red", "downregulation" = "blue")) +
#   scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   ggtitle('significant DEGs in Tfh\nbetween toddlers and adults\nin LAIV subtracted by ns') #+
#   # theme(legend.position="none")
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_GC_LAIV_wilcoxon_toddler_vs_adult_days_cytokines.pdf',sep = ''),width = 7, height = 18)

temp2$log2FC <- abs(temp2$log2FC_mean)
ggplot(temp2, aes(x=day_group, y=marker,color = effect)) + geom_point(aes(size = log2FC)) + #(1 + )
  RotatedAxis() + 
  facet_wrap(~effect) +
  scale_colour_manual(values = c("upregulation" = "red", "downregulation" = "blue")) +
  # scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(breaks = 1:4) +
  ggtitle('significant DEGs in Bcell_GC\nin adults vs toddlers') #+
# theme(legend.position="none")
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_GC_LAIV_wilcoxon_toddler_vs_adult_days_overlap.pdf',sep = ''),width = 5, height = 6)

##### visualization one direction #################################
Bcell_LAIV_module_wilcoxon_youth_vs_adult_downregulation_shrink <- 
  Bcell_LAIV_module_wilcoxon_youth_vs_adult_downregulation[Bcell_LAIV_module_wilcoxon_youth_vs_adult_downregulation$marker %in% 
                                                    overlapped_downregulation_marker_list,]
temp <- Bcell_LAIV_module_wilcoxon_youth_vs_adult_downregulation_shrink[,c('marker','days','log2FC_mean')]
temp <- temp[temp$days != '07-08',]
temp$day_group <- paste('day',temp$days,sep = '')
temp$day_group[temp$day_group == 'day12-14'] <- 'day14'
temp$effect <- 'upregulation'
temp2 <- data.frame(table(temp$marker,temp$day_group))
colnames(temp2) <- c('marker','day_group','Freq')
temp2$Freq[temp2$Freq > 0] <- 1
temp2$effect <- '+'

temp1 <- dcast(temp2, marker ~ day_group, value.var="Freq")
temp11 <- data.frame(temp1)
temp11 <- temp11 %>% arrange(day0,day04,day06.07,day10,day14)
# Hierarchical clustering using Complete Linkage

temp_ordered_marker <- temp11$marker
temp <- temp[!(temp$marker %in% c('ab-TCR-gamma-delta','ab-TCR-alpha-beta','ab-CD3','ab-CD4',"ab-CD8",'TRAC','TRBC2')),]
temp <- temp[!(temp$marker %in% c('CD38',"CXCR3","IGHG1-membrane","IL2RA","ENTPD1",'GPI',"IL2RB","PFKP","TPI1","LDHA","PGK1","ENO1")),]

temp$marker <- factor(temp$marker, levels = temp_ordered_marker)
temp$log2FC_mean <- abs(temp$log2FC_mean)
ggplot(temp, aes(x=day_group, y=marker,color = log2FC_mean)) + geom_point(aes(size = log2FC_mean)) + #aes(size = log2FC_mean)
  RotatedAxis() + 
  scale_size(range = c(0,4)) +
  scale_color_gradient(low="lightpink", high="#400000",limits=c(0,4)) +
  ggtitle('Upregulated DEGs in Bcell\nin adults than younger children')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_bootstrap_toddler_vs_adult_days_downregulation_overlapped.pdf',sep = ''),
          width = 6, height = 11)



######
temp_table1 <- read_excel('tonsil_LAIV_120a_s_Bcell_age_wilcoxon_DAVID_upregulation.xlsx')

View(temp_table2[duplicated(temp_table2$ENTREZID),])
rownames(temp_table2) <- temp_table2$ENTREZID
temp_table2 <- temp_table2 %>% arrange(desc(ENTREZID))
for (geneID_index in rownames(temp_table2)) {
  temp_table1$Genes <- gsub(paste(', ',geneID_index,',',sep = ''),paste(', ',temp_table2[geneID_index,][,c('marker')],',',sep = ''),temp_table1$Genes)
  temp_table1$Genes <- gsub(paste('^',geneID_index,',',sep = ''),paste(temp_table2[geneID_index,][,c('marker')],',',sep = ''),temp_table1$Genes)
  temp_table1$Genes <- gsub(paste(', ',geneID_index,'$',sep = ''),paste(', ',temp_table2[geneID_index,][,c('marker')],sep = ''),temp_table1$Genes)
}
write.xlsx(temp_table1,sheetName = 'Sheet2',temp_table1_name,append = T)

temp_table1$Term <- gsub('.*\\~','',temp_table1$Term)




####### 
temp <- Bcell_LAIV_module_wilcoxon[(Bcell_LAIV_module_wilcoxon$pval_adj <= 0.05) & 
                                              (Bcell_LAIV_module_wilcoxon$age_group2 == '27-39yrs') & 
                                              (!(is.na(Bcell_LAIV_module_wilcoxon$pval_adj))),]
temp <- unique(temp[,c('marker','effect','n_significant')])
temp <- temp %>% arrange(effect,desc(n_significant))

############## stringent comparison ##############
# stringent
Bcell_LAIV_module_wilcoxon_youth_stringent <- Bcell_LAIV_module_wilcoxon_youth
Bcell_LAIV_module_wilcoxon_stringent <- Bcell_LAIV_module_wilcoxon_stringent %>% group_by(day_group,marker,bootstrap_index) %>% 
  mutate(significant = ifelse(all((pval_adj <= 0.05) & (prod(diff) > 0)), 'yes','no'))
Bcell_LAIV_module_wilcoxon_stringent <- Bcell_LAIV_module_wilcoxon_stringent %>% group_by(day_group,marker,effect) %>% 
  mutate(n_significant = sum(significant == 'yes')/2)
Bcell_LAIV_module_wilcoxon_stringent <- Bcell_LAIV_module_wilcoxon_stringent %>% group_by(day_group,marker,effect) %>% 
  mutate(significant = ifelse(n_significant > floor(n_bootstrap_final*0.95), "yes", "no"))
Bcell_LAIV_module_wilcoxon_stringent <- Bcell_LAIV_module_wilcoxon_stringent %>% group_by(day_group,marker) %>% 
  mutate(diff_mean = mean(diff),p_mean = exp(mean(log(pval_adj))))
Bcell_LAIV_module_wilcoxon_stringent$label <- Bcell_LAIV_module_wilcoxon_stringent$marker
Bcell_LAIV_module_wilcoxon_stringent$label[Bcell_LAIV_module_wilcoxon_stringent$significant == 'no'] <- ''

Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <- Bcell_LAIV_module_wilcoxon_stringent[(Bcell_LAIV_module_wilcoxon_stringent$significant == 'yes') & 
                                                                                              (Bcell_LAIV_module_wilcoxon_stringent$effect == '+'),]
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <- unique(Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation[,c('day_group','marker','diff_mean','p_mean','n1','n2','ENTREZID')])
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <- unique(Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation %>% group_by(day_group,marker) %>% mutate(n2 = mean(n2)))
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation$days <- substring(Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation$day_group,first = 4,last = 30)
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <- Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation %>% arrange(days,p_mean) 
Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation <- Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation %>% group_by(day_group) %>% mutate(n = length(unique(marker)))
write.xlsx(data.frame(Bcell_LAIV_module_wilcoxon_youth_stringent_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_diff_stringent_toddler_vs_both_older.xlsx')

Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <- Bcell_LAIV_module_wilcoxon_stringent[(Bcell_LAIV_module_wilcoxon_stringent$significant == 'yes') & 
                                                                                                (Bcell_LAIV_module_wilcoxon_stringent$effect == '-'),]
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <- unique(Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation[,c('day_group','marker','diff_mean','p_mean','n1','n2','ENTREZID')])
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <- unique(Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation %>% group_by(day_group,marker) %>% mutate(n2 = mean(n2)))
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation$days <- substring(Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation$day_group,first = 4,last = 30)
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <- Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation %>% arrange(days,p_mean) 
Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation <- Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation %>% group_by(day_group) %>% mutate(n = length(unique(marker)))
write.xlsx(data.frame(Bcell_LAIV_module_wilcoxon_youth_stringent_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_diff_stringent_toddler_vs_both_older.xlsx',append = T)

for (day_group_name in LAIV_day_group_list[2:13]) {
  print(day_group_name)
  graphics.off()
  temp <- unique(Bcell_LAIV_module_wilcoxon_stringent[(Bcell_LAIV_module_wilcoxon_stringent$day_group == day_group_name),][,c('marker','diff_mean','p_mean','effect','significant','label')])
  plot<-ggplot(temp, aes(x=diff_mean, y=-log10(p_mean),color = significant)) +
    geom_point(size = 1,alpha = 0.5) +
    scale_colour_manual(values = c("yes" = "red", "no" = "darkgrey")) +
    ggtitle(paste('DE of Bcell expression change by LAIV\nbetween toddlers and non-toddlers on',day_group_name)) +
    xlab("Bcell gene expression change by LAIV") + ylab("-log10 adjusted p-value") +
    theme(text = element_text(size = 20)) +
    theme_bw() +
    geom_text_repel(aes(x=diff_mean, y=-log10(p_mean),label = label), color = "black", size = 3) +
    # geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff), linetype="dotted") +
    geom_hline(yintercept=c(-log10(0.05)), linetype="dotted")#+
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_module_wilcoxon_diff_stringent_toddler_vs_both_older_',day_group_name,'.pdf',sep = ''),width = 5.5, height = 5)
  
}
