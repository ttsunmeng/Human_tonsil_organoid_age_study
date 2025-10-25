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

Tcell_LAIV_wilcoxon <- readRDS('Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_LAIV_lognorm_age_wilcoxon_output_parallel_n300.rds')
############# examine output overlaps #######################################
Tcell_LAIV_wilcoxon <- Tcell_LAIV_wilcoxon[!is.na(Tcell_LAIV_wilcoxon$pval),]
Tcell_LAIV_wilcoxon$log2FC <- log2(Tcell_LAIV_wilcoxon$mean1/Tcell_LAIV_wilcoxon$mean2)
Tcell_LAIV_wilcoxon$log2FC[Tcell_LAIV_wilcoxon$log2FC %in% c(Inf,-Inf,NaN)] <-
  log2(Tcell_LAIV_wilcoxon$mean1[Tcell_LAIV_wilcoxon$log2FC %in% c(Inf,-Inf,NaN)] + 1)
all_marker_table <- read.xlsx('tonsil_LAIV_120a_s_background_list.xlsx', sheetName = "Sheet1")
p <- match(Tcell_LAIV_wilcoxon$marker, all_marker_table$marker)
temp <- all_marker_table[p,]
Tcell_LAIV_wilcoxon$ENTREZID <- temp$ENTREZID
Tcell_LAIV_wilcoxon <- Tcell_LAIV_wilcoxon %>% group_by(day_group,age_group1,age_group2) %>% 
  mutate(pval_adj = p.adjust(pval,method = 'fdr'))

Tcell_LAIV_wilcoxon <- Tcell_LAIV_wilcoxon %>% group_by(day_group,marker,age_group1,age_group2) %>% 
  mutate(pval_adj_mean = exp(mean(log(pval_adj))),log2FC_mean = mean(log2FC),pct1_mean = mean(pct1),pct2_mean = mean(pct2))

Tcell_LAIV_wilcoxon$effect <- '+'
Tcell_LAIV_wilcoxon$effect[Tcell_LAIV_wilcoxon$log2FC < 0] <- '-'
Tcell_LAIV_wilcoxon$significant <- 'no'
Tcell_LAIV_wilcoxon$significant[(Tcell_LAIV_wilcoxon$pval_adj <= p_cutoff)] <- 'yes'

Tcell_LAIV_wilcoxon <- Tcell_LAIV_wilcoxon %>% group_by(day_group,marker,age_group1,age_group2,effect) %>% 
  mutate(n_significant = sum(significant == 'yes'))
n_bootstrap_final <- max(Tcell_LAIV_wilcoxon$bootstrap_index)
Tcell_LAIV_wilcoxon <- Tcell_LAIV_wilcoxon %>% group_by(day_group,marker,age_group1,age_group2,effect) %>% 
  mutate(significant = ifelse((n_significant > floor(n_bootstrap_final*0.95)), "yes", "no"))#
Tcell_LAIV_wilcoxon$label <- Tcell_LAIV_wilcoxon$marker
Tcell_LAIV_wilcoxon$label[Tcell_LAIV_wilcoxon$significant == 'no'] <- ''
Tcell_LAIV_wilcoxon <- Tcell_LAIV_wilcoxon %>% group_by(day_group,marker,age_group1,age_group2) %>% 
  mutate(effect = ifelse(log2FC_mean > 0, "+", "-"))#

Tcell_LAIV_wilcoxon <- Tcell_LAIV_wilcoxon[!Tcell_LAIV_wilcoxon$day_group %in% c("day14",'day07'),]
LAIV_day_group_list <- as.character(unique(Tcell_LAIV_wilcoxon$day_group))
# LAIV_day_group_list_youth <- LAIV_day_group_list#LAIV_day_group_list[!LAIV_day_group_list %in% c("day14",'day07')]

Tcell_LAIV_wilcoxon_youth <- Tcell_LAIV_wilcoxon[(Tcell_LAIV_wilcoxon$age_group1 == '02-03yrs'),]
Tcell_LAIV_wilcoxon_youth_vs_adult <- Tcell_LAIV_wilcoxon_youth[(Tcell_LAIV_wilcoxon_youth$age_group2 == '27-39yrs'),]
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays <- Tcell_LAIV_wilcoxon_youth_vs_adult
Tcell_LAIV_wilcoxon_youth_vs_adult <- Tcell_LAIV_wilcoxon_youth_vs_adult[Tcell_LAIV_wilcoxon_youth_vs_adult$day_group %in% LAIV_day_group_list,]

# temp <- Tcell_LAIV_wilcoxon_youth_vs_adult[Tcell_LAIV_wilcoxon_youth_vs_adult$marker == 'CCR10',]

Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult[(Tcell_LAIV_wilcoxon_youth_vs_adult$significant == 'yes') & 
                                            (Tcell_LAIV_wilcoxon_youth_vs_adult$effect == '+'),]
Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation <-
  unique(Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation$days <- 
  substring(Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation$day_group,first = 4,last = 30)
Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation <- Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation %>% arrange(days,desc(log2FC_mean)) 
Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation <- Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.5))
# Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation <- 
#   Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation[abs(Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation$log2FC_mean) > 0.25,]
# 
write.xlsx(data.frame(Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_s_Tcells_LAIV_lognorm_wilcoxon_bootstrap_n300_toddler_vs_adult.xlsx')
# unique(Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation[,c('day_group','n_significant','n_distinct')])


Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult[(Tcell_LAIV_wilcoxon_youth_vs_adult$significant == 'yes') & 
                                            (Tcell_LAIV_wilcoxon_youth_vs_adult$effect == '-'),]
Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation <-
  unique(Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$days <- 
  substring(Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$day_group,first = 4,last = 30)
Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation <- Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation %>% arrange(days,log2FC_mean) 
Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation <- Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.5))
# Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation <- 
#   Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation[abs(Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$log2FC_mean) > 0.5,]
write.xlsx(data.frame(Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_s_Tcells_LAIV_lognorm_wilcoxon_bootstrap_n300_toddler_vs_adult.xlsx',append = T)

# temp1 <- Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation[Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$day_group != 'day0',]
# temp1 <- temp1[temp1$log2FC_mean < -1,]
# temp2 <- Tcell_memory_LAIV_wilcoxon_youth_vs_adult_downregulation$marker[Tcell_memory_LAIV_wilcoxon_youth_vs_adult_downregulation$day_group == 'day0']
# # 
# # day_group_name <- 'day14'
# # temp1 <- Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$marker[Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$day_group == day_group_name]
# # temp2 <- Tcell_memory_LAIV_wilcoxon_youth_vs_adult_downregulation$marker[Tcell_memory_LAIV_wilcoxon_youth_vs_adult_downregulation$day_group == day_group_name]
# # temp1[!temp1 %in% temp2]
# # # temp2[!temp2 %in% temp1]
# 
# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp <- unique(Tcell_LAIV_wilcoxon_youth_vs_adult[(Tcell_LAIV_wilcoxon_youth_vs_adult$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$log2FC_mean <- -temp$log2FC_mean
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp, aes(x=log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of GC Tcell expression with LAIV\nbetween toddlers and adults on',day_group_name,'\nhigher in adults higher in toddlers\n')) +
#     xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_bootstrap_toddler_vs_adult_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }

# visualization of all day groups together
temp <- Tcell_LAIV_wilcoxon_youth_vs_adult[!Tcell_LAIV_wilcoxon_youth_vs_adult$day_group %in% c('day07-08'),]
temp$log2FC_mean <- -temp$log2FC_mean
temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','effect','significant','label')])

temp$sig_sign <- 'not sig.'
temp$sig_sign[(temp$log2FC_mean > 0.25) & (temp$significant == 'yes')] <- 'Higher in Adults'
temp$sig_sign[(temp$log2FC_mean <= 0.25) & (temp$significant == 'yes')] <- 'Lower in Adults'
DEG_count <- data.frame(dplyr::count(temp,day_group,sig_sign))
DEG_count <- temp %>% group_by(day_group,sig_sign) %>%
  mutate(count = n())
DEG_count <- unique(DEG_count[,c('day_group','sig_sign','count')])
DEG_count <- DEG_count[DEG_count$sig_sign != 'not sig.',]
DEG_count$day_group <- factor(DEG_count$day_group, levels = c('day0','day04','day06-07','day10','day12-14'))
ggplot(DEG_count, aes(fill=sig_sign, y=count, x=day_group)) + 
  geom_bar(position="stack", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = c("Higher in Adults" = "red",'Lower in Adults' = 'blue')) +
  ylab('#DEGs of Tcells') + RotatedAxis()
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_lognorm_bootstrap_n300_toddler_vs_adult_numDEGs.pdf',sep = ''),width = 3.5, height = 3)

temp <- Tcell_LAIV_wilcoxon_youth_vs_adult[!Tcell_LAIV_wilcoxon_youth_vs_adult$day_group %in% c('day07-08'),]
temp$log2FC_mean <- -temp$log2FC_mean
temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'no'
temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'yes'
  # as.character(temp$day_group[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')])
temp <- temp[!(temp$marker %in% temp$marker[(temp$significant == 'yes') & (temp$day_group == 'day0')]),]
temp <- temp[!(temp$marker %in% intersect(temp$marker[(temp$sig_sign == 'Higher in Adults')],temp$marker[(temp$sig_sign == 'Lower in Adults')])),]

temp <- temp[temp$day_group != 'day0',]

temp$label <- as.character(temp$label)
temp$label[temp$significant == 'no'] <- ''
# temp$label[temp$marker %in% temp$marker[(temp$significant == 'yes') & (temp$day_group == 'day0')]] <- ''
temp$label[temp$label %in% c('ab-CD3','ab-Tcell','ab-TCR-gamma-delta','ab-CD8','ab-TCR-alpha-beta')] <- ''
temp$label[temp$label != ''] <- paste(temp$day_group[temp$label != ''],temp$label[temp$label != ''],sep = '_')
# temp$label[!(temp$label %in% c('day06-07_GPI','day10_CD27','day10_GPI','day04_GPI','day10_ISG15','day06-07_JCHAIN','day12-14-GPI','day06-07_CD27','day04_ITGAL','day10_JCHAIN','day10_IGHG1_secreted','day06-07_IGHG2_secreted','day06-07_XBP1',
#                                'day10_CCR2','day12-14_CCR2','day04_MKI67','day04_STAT3','day10_PRDM1','day10_SELPLG','day06-07_SELPLG','day10_LGALS3',
#                                # 'day06-07_CCR10','day10_CCR10','day12-14_CCR10','day04_CXCR3','day06-07_CXCR3','day10_CXCR3','day12-14_CXCR3','day06-07_MKI67','day04_IL2RA','day04_IL2RB','day06-07_IL2RA','day06-07_IL2RB','day10_IL2RA','day10_IL2RB','day12-14_IL2RA','day12-14_IL2RB',
#                                # 'day04_VEGFA','day12-14_VEGFA','day04_TBX21','day04_ITGAX','day10_ITGAX','day12-14_ITGAX','day04_ICAM1','day06-07_ICAM1','day10_ICAM1','day12-14_ICAM1','day04_ITGB7',
#                                'day06-07_CD1A','day06-07_ab-IgD','day12-14_ab-IgD','day10_ab-IgD','day04_ab-IgD','day10_CD83','day12-14_IL4R','day06-07_KLRK1','day06-07_CD83','day10_IL4R','day12-14_IFNGR1'))] <- ''

ggplot(temp, aes(x=log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_colour_manual(values = c("yes" = "red","no" = "darkgrey")) +
  ggtitle(paste('DEG of Tcells not presented at baseline\nhigher in adults higher in toddlers\n')) +
  xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 10)) +
  theme_bw() +
  geom_text_repel(aes(x=log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
  geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n300_toddler_vs_adult_big.pdf',sep = ''),width = 15, height = 15)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n300_toddler_vs_adult_2.pdf',sep = ''),width = 6.5, height = 5)


temp <- Tcell_LAIV_wilcoxon_youth_vs_adult
temp$log2FC_mean <- -temp$log2FC_mean
temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'no'
temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'only at baseline'
# need to be same direction
temp$significant[(temp$marker %in% temp$marker[(temp$significant == 'only at baseline') & (temp$day_group != 'day0') & (temp$log2FC_mean > 0)]) & 
                   (temp$log2FC_mean > 0) & (temp$significant == 'only at baseline')] <- 'persist'
temp$significant[(temp$marker %in% temp$marker[(temp$significant == 'only at baseline') & (temp$day_group != 'day0') & (temp$log2FC_mean < 0)]) & 
                   (temp$log2FC_mean < 0) & (temp$significant == 'only at baseline')] <- 'persist'
temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'no'
temp$label <- as.character(temp$label)
# temp$color <- 'black'
# temp$color[temp$significant == 'only at baseline'] <- 'grey'
temp$label[temp$significant == 'no'] <- ''
# temp$label[temp$marker %in% temp$marker[(temp$significant == 'yes') & (temp$day_group == 'day0')]] <- ''
temp <- temp[!temp$label %in% c('ab-CD3','ab-Tcell','ab-TCR-gamma-delta','ab-CD8','ab-TCR-alpha-beta'),] 
temp <- temp[temp$day_group == 'day0',]
ggplot(temp, aes(x=log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
  geom_point(size = 1,alpha = 0.5) +
  scale_colour_manual(values = c("persist" = "red","only at baseline" = "orange","no" = "darkgrey")) +
  ggtitle(paste('DEG of Tcells not presented at baseline\nhigher in adults higher in toddlers\n')) +
  xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 10)) +
  theme_bw() +
  geom_text_repel(aes(x=log2FC_mean, y=-log10(pval_adj_mean),label = label), color = 'black', size = 3) +
  geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n300_toddler_vs_adult_day0.pdf',sep = ''),width = 6.5, height = 5)


temp <- read.xlsx('./Tcell_lognorm_CCA/tonsil_LAIV_120a_s_Tcell_BCL6pos_lognorm_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'Sheet1')
temp <- temp[temp$log2FC_mean < -1,]

###### venn diagram visualization ###############################
library('xlsx')
Tcell_lognorm_DEG_upreg <- read.xlsx('Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcells_LAIV_lognorm_wilcoxon_bootstrap_n300_toddler_vs_adult.xlsx',sheetName = 'downregulation')
Tcell_lognorm_DEG_downreg <- read.xlsx('Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcells_LAIV_lognorm_wilcoxon_bootstrap_n300_toddler_vs_adult.xlsx',sheetName = 'upregulation')
Tcell_lognorm_DEG_upreg_day0 <- Tcell_lognorm_DEG_upreg[Tcell_lognorm_DEG_upreg$day_group == 'day0',]
Tcell_lognorm_DEG_upreg_nonday0 <- Tcell_lognorm_DEG_upreg[Tcell_lognorm_DEG_upreg$day_group != 'day0',]
Tcell_lognorm_DEG_downreg_day0 <- Tcell_lognorm_DEG_downreg[Tcell_lognorm_DEG_downreg$day_group == 'day0',]
Tcell_lognorm_DEG_downreg_nonday0 <- Tcell_lognorm_DEG_downreg[Tcell_lognorm_DEG_downreg$day_group != 'day0',]

# devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)

x <- list(
  up_day0 = unique(Tcell_lognorm_DEG_upreg_day0$marker), 
  up_later = unique(Tcell_lognorm_DEG_upreg_nonday0$marker), 
  down_later = unique(Tcell_lognorm_DEG_downreg_nonday0$marker),
  down_day0 = unique(Tcell_lognorm_DEG_downreg_day0$marker)
)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_lognorm_age_DEG_venn_up&down.pdf',width = 5, height = 5)

# devtools::install_github("gaospecial/ggVennDiagram")
# https://cran.r-project.org/web/packages/ggVennDiagram/readme/README.html
# There are other types of visualization
library(ggVennDiagram)
x <- list(
  up_day0 = unique(Tcell_lognorm_DEG_upreg_day0$marker), 
  up_day4 = unique(Tcell_lognorm_DEG_upreg_nonday0$marker[Tcell_lognorm_DEG_upreg_nonday0$day_group == 'day04']), 
  up_day7 = unique(Tcell_lognorm_DEG_upreg_nonday0$marker[Tcell_lognorm_DEG_upreg_nonday0$day_group == 'day06-07']), 
  up_day10 = unique(Tcell_lognorm_DEG_upreg_nonday0$marker[Tcell_lognorm_DEG_upreg_nonday0$day_group == 'day10']),
  up_day14 = unique(Tcell_lognorm_DEG_upreg_nonday0$marker[Tcell_lognorm_DEG_upreg_nonday0$day_group == 'day12-14'])
)
ggVennDiagram(
  x, 
  fill_color = c("#868686FF","#0073C2FF", "#EFC000FF", "#CD534CFF",'#BAE4B3'),
  stroke_size = 0.5
)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_lognorm_age_DEG_upregulation_venn_days.pdf',width = 5, height = 5)

ggVennDiagram(x,force_upset = TRUE,order.set.by = "none",order.intersect.by = "none")
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_lognorm_age_DEG_upregulation_bars_intersect_days.pdf',width = 5, height = 4)

x <- list(
  down_day0 = unique(Tcell_lognorm_DEG_downreg_day0$marker), 
  down_day4 = unique(Tcell_lognorm_DEG_downreg_nonday0$marker[Tcell_lognorm_DEG_downreg_nonday0$day_group == 'day04']), 
  down_day7 = unique(Tcell_lognorm_DEG_downreg_nonday0$marker[Tcell_lognorm_DEG_downreg_nonday0$day_group == 'day06-07']), 
  down_day10 = unique(Tcell_lognorm_DEG_downreg_nonday0$marker[Tcell_lognorm_DEG_downreg_nonday0$day_group == 'day10']),
  down_day14 = unique(Tcell_lognorm_DEG_downreg_nonday0$marker[Tcell_lognorm_DEG_downreg_nonday0$day_group == 'day12-14'])
)
ggVennDiagram(
  x, 
  fill_color = c("#868686FF","#0073C2FF", "#EFC000FF", "#CD534CFF",'#BAE4B3'),
  stroke_size = 0.5
)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_lognorm_age_DEG_downregulation_venn_days.pdf',width = 5, height = 5)

ggVennDiagram(x,force_upset = TRUE,order.set.by = "none",order.intersect.by = "none")
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_lognorm_age_DEG_downregulation_bars_intersect_days.pdf',width = 5, height = 4)

####### 
# select all cytokines
cytokine_global_table <- read.xlsx('../Global landscape of cytokines Supplementary Table S1.xlsx', sheetName = "Cytokines")
cytokine_global_list <- unique(cytokine_global_table$HGNC.symbol)
cytokine_global_list <- c(cytokine_global_list,'GZMB','PRF1','MIF','FASLG','VEGFA','TGFA')
seurat_cytokine_list <- all.markers[all.markers %in% cytokine_global_list]

temp <- Tcell_LAIV_wilcoxon_youth_vs_adult[!Tcell_LAIV_wilcoxon_youth_vs_adult$day_group %in% c('day07-08'),]
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
  ggtitle(paste('DEG of Tcells not presented at baseline\nhigher in adults higher in toddlers\n')) +
  xlab("log2 fold-change") + ylab("-log10 adjusted p-value") +
  theme(text = element_text(size = 10)) +
  theme_bw() +
  geom_text_repel(aes(x=log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
  geom_vline(xintercept=c(-log2FC_cutoff,log2FC_cutoff), linetype="dotted") +
  geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n300_toddler_vs_adult_big.pdf',sep = ''),width = 15, height = 15)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n300_toddler_vs_adult_cytokine_only.pdf',sep = ''),width = 10, height = 10)

all_Tcell_LAIV_wilcoxon_youth_vs_adult <- Tcell_LAIV_wilcoxon_youth_vs_adult
NaiveOrActB_LAIV_wilcoxon_youth_vs_adult <- Tcell_LAIV_wilcoxon_youth_vs_adult
PB_LAIV_wilcoxon_youth_vs_adult <- Tcell_LAIV_wilcoxon_youth_vs_adult
GClikeB_LAIV_wilcoxon_youth_vs_adult <- Tcell_LAIV_wilcoxon_youth_vs_adult
MemoryB_LAIV_wilcoxon_youth_vs_adult <- Tcell_LAIV_wilcoxon_youth_vs_adult
preGCOrGCB_LAIV_wilcoxon_youth_vs_adult <- Tcell_LAIV_wilcoxon_youth_vs_adult

gene_name <- 'IL15'
gene_name_list <- c('TNFSF10','IL15','CXCL10','CCL5','TNFSF13','VEGFA','CXCL16','IL32','GZMB')
FC_table <- expand.grid(gene = gene_name_list,subset = c('Naive_activated','preGC_Folli_GC','GC_like','memory','PB'))
FC_table$log2FC <- NA
library(reshape2)
FC_table <- dcast(FC_table, gene ~ subset, value.var="log2FC")
rownames(FC_table) <- FC_table$gene
for (gene_name in gene_name_list) {
  temp <- NaiveOrActB_LAIV_wilcoxon_youth_vs_adult
  # temp <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
  # if (length(temp) == 0) {
  #   FC_table[FC_table$gene == gene_name,][,'Naive_activated'] <- 0
  # } else {
    temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
    temp <- temp[(temp$day_group != 'day0') & (temp$marker == gene_name),]
    temp$log2FC_mean <- -temp$log2FC_mean
    temp <- temp %>% arrange(desc(log2FC_mean))
    FC_table[FC_table$gene == gene_name,][,'Naive_activated'] <- mean(temp$log2FC_mean)
  # }
  
  temp <- preGCOrGCB_LAIV_wilcoxon_youth_vs_adult
  # temp <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
  # if (length(temp) == 0) {
  #   FC_table[FC_table$gene == gene_name,][,'preGC_Folli_GC'] <- 0
  # } else {
    temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
    temp <- temp[(temp$day_group != 'day0') & (temp$marker == gene_name),]
    temp$log2FC_mean <- -temp$log2FC_mean
    temp <- temp %>% arrange(desc(log2FC_mean))
    FC_table[FC_table$gene == gene_name,][,'preGC_Folli_GC'] <- mean(temp$log2FC_mean)
  # }

  temp <- GClikeB_LAIV_wilcoxon_youth_vs_adult
  # temp <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
  # if (length(temp) == 0) {
  #   FC_table[FC_table$gene == gene_name,][,'GC_like'] <- 0
  # } else {
    temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
    temp <- temp[(temp$day_group != 'day0') & (temp$marker == gene_name),]
    temp$log2FC_mean <- -temp$log2FC_mean
    temp <- temp %>% arrange(desc(log2FC_mean))
    FC_table[FC_table$gene == gene_name,][,'GC_like'] <- mean(temp$log2FC_mean)
  # }

  temp <- MemoryB_LAIV_wilcoxon_youth_vs_adult
  # temp <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
  # if (length(temp) == 0) {
  #   FC_table[FC_table$gene == gene_name,][,'memory'] <- 0
  # } else {
    temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
    temp <- temp[(temp$day_group != 'day0') & (temp$marker == gene_name),]
    temp$log2FC_mean <- -temp$log2FC_mean
    temp <- temp %>% arrange(desc(log2FC_mean))
    FC_table[FC_table$gene == gene_name,][,'memory'] <- mean(temp$log2FC_mean)
  # }

  temp <- PB_LAIV_wilcoxon_youth_vs_adult
  # temp <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
  # if (length(temp) == 0) {
  #   FC_table[FC_table$gene == gene_name,][,'PB'] <- 0
  # } else {
    temp <- unique(temp[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
    temp <- temp[(temp$day_group != 'day0') & (temp$marker == gene_name),]
    temp$log2FC_mean <- -temp$log2FC_mean
    temp <- temp %>% arrange(desc(log2FC_mean))
    FC_table[FC_table$gene == gene_name,][,'PB'] <- mean(temp$log2FC_mean)
  # }

}
FC_table$gene <- NULL
FC_table <- mutate_all(FC_table, function(x) as.numeric(x))
# temp <- melt(FC_table)
# ggplot(temp, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile() +
#   labs(title = "Correlation Heatmap",
#        x = "Variable 1",
#        y = "Variable 2")
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
col_fun(seq(-3, 3))
HM2 <- Heatmap(FC_table, cluster_columns = F, cluster_rows = T,
               name = "mean\nnormalized\nexpression", show_column_names = T,col = col_fun,
               show_row_names = T,#row_km = 20, column_km = 20,
               column_title_side = "top",row_title_side = "left",row_names_side = 'left')
HM2 <- draw(HM2)
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_LAIV_cytokine_DEG_log2FC.pdf',width = 6, height = 4)


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
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_alldays[(Tcell_LAIV_wilcoxon_youth_vs_adult_alldays$effect == '+'),] 
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_alldays[(Tcell_LAIV_wilcoxon_youth_vs_adult_alldays$effect == '-'),] 

Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation %>% group_by(marker) %>% 
  summarize(log2FCmax = -max(log2FC_mean),day_group = day_group[which.max(log2FC_mean)],if_significant = significant[which.max(log2FC_mean)])
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation %>% group_by(marker) %>% 
  summarize(log2FCmax = max(-log2FC_mean),day_group = day_group[which.max(-log2FC_mean)],if_significant = significant[which.max(-log2FC_mean)])

Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation_sig <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_alldays[(Tcell_LAIV_wilcoxon_youth_vs_adult_alldays$significant == 'yes') & 
                                               (Tcell_LAIV_wilcoxon_youth_vs_adult_alldays$effect == '+'),] 
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation_sig <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_alldays[(Tcell_LAIV_wilcoxon_youth_vs_adult_alldays$significant == 'yes') & 
                                               (Tcell_LAIV_wilcoxon_youth_vs_adult_alldays$effect == '-'),] 
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation_sig <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation_sig %>% group_by(marker) %>% 
  summarize(log2FCmax = -max(log2FC_mean),day_group = day_group[which.max(log2FC_mean)])
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation_sig <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation_sig %>% group_by(marker) %>% 
  summarize(log2FCmax = max(-log2FC_mean),day_group = day_group[which.max(-log2FC_mean)])

Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation$celltype <- 'BCL6-nonmemory B'
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation$celltype <- 'BCL6-nonmemory B'
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation_sig$celltype <- 'BCL6-nonmemory B'
Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation_sig$celltype <- 'BCL6-nonmemory B'

Tcell_subset_DEG_downregulation <- rbind(Tcell_subset_DEG_downregulation,Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation)
Tcell_subset_DEG_upregulation <- rbind(Tcell_subset_DEG_upregulation,Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation)
Tcell_subset_DEG_downregulation_sig <- rbind(Tcell_subset_DEG_downregulation_sig,Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_upregulation_sig)
Tcell_subset_DEG_upregulation_sig <- rbind(Tcell_subset_DEG_upregulation_sig,Tcell_LAIV_wilcoxon_youth_vs_adult_alldays_downregulation_sig)

saveRDS(Tcell_subset_DEG_downregulation,'Tcell_lognorm_DEG_age_downregulation.rds')
saveRDS(Tcell_subset_DEG_upregulation,'Tcell_lognorm_DEG_age_upregulation.rds')
saveRDS(Tcell_subset_DEG_downregulation_sig,'Tcell_lognorm_DEG_age_downregulation_sig.rds')
saveRDS(Tcell_subset_DEG_upregulation_sig,'Tcell_lognorm_DEG_age_upregulation_sig.rds')

Tcell_subset_DEG_noday0_downregulation_sig <- Tcell_subset_DEG_downregulation_sig[Tcell_subset_DEG_downregulation_sig$day_group != 'day0',]
Tcell_subset_DEG_noday0_downregulation_sig <- Tcell_subset_DEG_noday0_downregulation_sig %>% group_by(marker) %>%
  mutate(count = n())
Tcell_subset_DEG_noday0_upregulation_sig <- Tcell_subset_DEG_upregulation_sig[Tcell_subset_DEG_upregulation_sig$day_group != 'day0',]
Tcell_subset_DEG_noday0_upregulation_sig <- Tcell_subset_DEG_noday0_upregulation_sig %>% group_by(marker) %>%
  mutate(count = n()) 


length(unique(Tcell_subset_DEG_upregulation_sig$marker[Tcell_subset_DEG_upregulation_sig$count == 3]))
########## heatmap of same markers across different population ################
library(reshape2)
library(ComplexHeatmap)
temp <- Tcell_subset_DEG_noday0_upregulation_sig[Tcell_subset_DEG_noday0_upregulation_sig$count == 3,]
temp$log2FCmax <- as.numeric(temp$log2FCmax)

temp <- dcast(temp, marker ~ celltype, value.var="log2FCmax")
rownames(temp) <- temp$marker
temp$marker <- NULL
temp <- temp %>% arrange(desc(rowMeans(temp[,c('BCL6-nonmemory B','BCL6+B')])-temp$`memory B`))

# temp <- temp %>% arrange(desc(rowMeans(temp)))
temp <- temp[rowSums(temp > 0.25) == 3,]

Heatmap(as.matrix(temp), cluster_columns = F, cluster_rows = F, 
                name = "log2FC", show_column_names = T,#
                show_row_names = T,column_title = 'Tcell', #row_km = 3, #column_km = 20, 
                column_title_side = "top",column_names_side = "top",row_title_side = "left",row_names_side = 'left')
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_celltype_logFCmax_upregulation.pdf',sep = ''),width = 5, height = 6)



############ ClusterProfiler ############################
# DAVID website 

# library("org.Hs.eg.db")
# hs <- org.Hs.eg.db
library(enrichplot)
library(clusterProfiler)
for (day_group_name in LAIV_day_group_list){
  print(day_group_name)
  temp <- Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation[Tcell_LAIV_wilcoxon_youth_vs_adult_upregulation$day_group == day_group_name,]
  
  ego <- enrichGO(temp$ENTREZID[!is.na(temp$ENTREZID)], OrgDb = "org.Hs.eg.db",
                  universe = all_marker_table$ENTREZID, ont="BP", readable=TRUE)
  
  if (dim(ego)[1] != 0) {
    print(paste(day_group_name,'upregulation'))
    write.xlsx(ego@result,paste('tonsil_LAIV_120a_s_lognorm_Tcell_preGCOrGCB_LAIV_wilcoxon_bootstrap_nmean_toddler_vs_adult_GO_ClusterProfiler.xlsx',sep = ''),
               sheetName = paste(day_group_name,'upregulation'),append = T)
  }
  temp <- Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation[Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$day_group == day_group_name,]
  
  ego <- enrichGO(temp$ENTREZID[!is.na(temp$ENTREZID)], OrgDb = "org.Hs.eg.db", 
                  universe = all_marker_table$ENTREZID, ont="BP", readable=TRUE)
  
  if (dim(ego)[1] != 0) {
    print(paste(day_group_name,'downregulation'))
    write.xlsx(ego@result,paste('tonsil_LAIV_120a_s_lognorm_Tcell_preGCOrGCB_LAIV_wilcoxon_bootstrap_nmean_toddler_vs_adult_GO_ClusterProfiler.xlsx',sep = ''),
               sheetName = paste(day_group_name,'downregulation'),append = T)
    graphics.off()
    plot <- dotplot(ego) + ggtitle(paste('GO of Tcell LAIV vs ns\ndownregulated on',day_group_name)) +
      theme(axis.text = element_text(size = 20),plot.title = element_text(size = 25, face = "bold"))
    print(plot)
    dev.print(pdf, paste('tonsil_LAIV_120a_s_lognorm_Tcell_LAIV_preGCOrGCB_wilcoxon_bootstrap_nmean_toddler_vs_adult_GO_ClusterProfiler_downregulation_',day_group_name,'.pdf',sep = ''),width = 8, height = 10)
  }
}
   
# network visualization
bp <- pairwise_termsim(ego)
plot <- emapplot(bp,showCategory = 200)
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n316_toddler_vs_adult_GO_ClusterProfiler_upregulation_day07-08.pdf',sep = ''),
          width = 10, height = 10)
# simplify
bp <- simplify(ego, by="p.adjust", select_fun=min)
write.xlsx(bp@result,paste('tonsil_LAIV_120a_s_Tcell_DE_NS_vs_day0_',donor_name,'_',day_name,'_log1p_paired_GO_ClusterProfiler_downregulation_simplify.xlsx',sep = ''))
plot <- dotplot(bp, showCategory=20) + ggtitle(paste('GO of Tcell LAIV vs ns\ndownregulated on',day_name)) +
  theme(axis.text = element_text(size = 20),plot.title = element_text(size = 25, face = "bold"))
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_DE_NS_vs_day0_',donor_name,'_',day_name,'_log1p_paired_GO_ClusterProfiler_downregulation_simplify.pdf',sep = ''),width = 8, height = 10)
bp <- pairwise_termsim(bp)
plot <- emapplot(bp)
print(plot)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_DE_NS_vs_day0_',donor_name,'_',day_name,'_log1p_paired_GO_ClusterProfiler_downregulation_network_simplify.pdf',sep = ''),
          width = 7, height = 7)
  
  
table1 <- 'tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_toddler_vs_adult_downregulation_day04-06_DAVID_UP_KW_BIOLOGICAL_PROCESS.txt'
table1_DE_data <- read.table(table1, sep = '\t',header = T)
day_group_name <- 'day04-06'
temp <- Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation %>% dplyr::filter(day_group == day_group_name)
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
#                       'cell-cell adhesion','movement of cell or suTcellular component','Wnt signaling pathway, planar cell polarity pathway')
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
  ggtitle(paste('Tcell expression diff by LAIV on',day_group_name,'\nlower in toddlers')) +
  xlab(expression(paste('UniprotKB Keywords\nBiological Processes'))) + ylab('#genes') + 
  theme_bw() +
  scale_fill_gradientn(colors = c('yellow','#984EA3'),limits = c(0,0.05)) +
  scale_x_discrete(labels = wrap_format(35)) +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10),
        plot.title = element_text(size = 12, face = "bold"))
# theme(axis.text = element_text(size = 20),plot.title = element_text(size = 25, face = "bold"))
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_toddler_vs_adult_downregulation_day04-06_UP_KW_BP.pdf',sep = ''),width = 4.5, height = 2.5)

#### GSEA ########################
library(msigdbr)
# msigdbr_species()
# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures
m_t2g <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene)

CD4_RSTR_LTBI.markers <- CD4_RSTR_LTBI.markers %>% arrange(desc(avg_log2FC))

CD4_RSTR_LTBI.markers_filtered <- CD4_RSTR_LTBI.markers[!is.na(CD4_RSTR_LTBI.markers$geneID),]
geneList <- CD4_RSTR_LTBI.markers_filtered$avg_log2FC
names(geneList) <- CD4_RSTR_LTBI.markers_filtered$geneID
em2 <- GSEA(geneList, TERM2GENE = m_t2g)
temp_output <- em2@result
rownames(temp_output) <- NULL
temp_output$Description <- NULL
all_gene_table <- CD4_RSTR_LTBI.markers
all_gene_table <- all_gene_table[!is.na(all_gene_table$geneID),][,c('gene','geneID')]
all_gene_table$geneID <- as.numeric(all_gene_table$geneID)
rownames(all_gene_table) <- all_gene_table$geneID
all_gene_table <- all_gene_table %>% arrange(desc(geneID))
for (geneID_index in rownames(all_gene_table)) {
  temp_output$core_enrichment <- gsub(paste('\\/',geneID_index,'\\/',sep = ''),paste('\\/',all_gene_table[geneID_index,][,c('gene')],'\\/',sep = ''),temp_output$core_enrichment)
  temp_output$core_enrichment <- gsub(paste('^',geneID_index,'\\/',sep = ''),paste(all_gene_table[geneID_index,][,c('gene')],'\\/',sep = ''),temp_output$core_enrichment)
  temp_output$core_enrichment <- gsub(paste('\\/',geneID_index,'$',sep = ''),paste('\\/',all_gene_table[geneID_index,][,c('gene')],sep = ''),temp_output$core_enrichment)
}
temp_output$core_enrichment <- gsub('\\/',', ',temp_output$core_enrichment)
temp_output <- temp_output %>% arrange(desc(NES))
write.xlsx(temp_output,'TB_RSTR_PP1_CD4_lognorm_cleaned_DE_RSTR_vs_LTBI_GSEA_msigdbr_C5.xlsx')


############# two age group comparison ##################
age_group1 <- '02-03yrs'
age_group2 <- '07-09yrs'
temp <- Tcell_LAIV_wilcoxon[(Tcell_LAIV_wilcoxon$age_group1 == age_group1),]
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
           'tonsil_LAIV_120a_sTcell_lognorm_LAIV_bootstrap_toddlers_vs_tweens.xlsx')

temp_downregulation <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
temp_downregulation <-
  unique(temp_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
temp_downregulation$days <- 
  substring(temp_downregulation$day_group,first = 4,last = 30)
temp_downregulation <- temp_downregulation %>% arrange(days,log2FC_mean) 
temp_downregulation <- temp_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(temp_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_sTcell_lognorm_LAIV_bootstrap_toddlers_vs_tweens.xlsx',append = T)

# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp1 <- unique(temp[(temp$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp1, aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of Tcell expression change by LAIV\nbetween toddlers vs tweens on',day_group_name)) +
#     xlab("log2 fold-change of expression diff") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-FC_cutoff,FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_toddlers_vs_tweens_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }

age_group1 <- '07-09yrs'
age_group2 <- '27-39yrs'
temp <- Tcell_LAIV_wilcoxon[(Tcell_LAIV_wilcoxon$age_group1 == age_group1),]
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
           'tonsil_LAIV_120a_sTcell_lognorm_LAIV_bootstrap_tweens_vs_adults.xlsx')

temp_downregulation <- temp[(temp$significant == 'yes') & (temp$effect == '-'),]
temp_downregulation <-
  unique(temp_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','pct1_mean','pct2_mean','n1','n2','ENTREZID')])
temp_downregulation$days <- 
  substring(temp_downregulation$day_group,first = 4,last = 30)
temp_downregulation <- temp_downregulation %>% arrange(days,log2FC_mean) 
temp_downregulation <- temp_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(temp_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_sTcell_lognorm_LAIV_bootstrap_tweens_vs_adults.xlsx',append = T)

# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp1 <- unique(temp[(temp$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp1, aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of Tcell expression change by LAIV\nbetween tweens vs adults on',day_group_name)) +
#     xlab("log2 fold-change of expression diff") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-FC_cutoff,FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_tweens_vs_adults_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }
############# stringent in toddler vs both ########################
Tcell_LAIV_wilcoxon_youth_stringent <- Tcell_LAIV_wilcoxon_youth %>% group_by(day_group,marker) %>% 
  mutate(significant = ifelse(all(significant == 'yes') & (prod(log2FC_mean) > 0),'yes','no'))
        
        
# temp <- unique(Tcell_LAIV_wilcoxon_youth_stringent[,c('day_group','marker','age_group1','age_group2','significant','significant2','log2FC_mean','effect')])
# temp <- temp %>% arrange(day_group,marker)

Tcell_LAIV_wilcoxon_youth_stringent <- Tcell_LAIV_wilcoxon_youth_stringent %>% group_by(day_group,marker) %>% 
  mutate(pval_adj_mean = exp(mean(log(unique(pval_adj_mean)))),log2FC_mean = mean(unique(log2FC_mean)),n2 = mean(unique(n2)))
Tcell_LAIV_wilcoxon_youth_stringent <- Tcell_LAIV_wilcoxon_youth_stringent %>% group_by(day_group,marker) %>% 
  mutate(effect = ifelse(log2FC_mean > 0, "+", "-"))
Tcell_LAIV_wilcoxon_youth_stringent$label[Tcell_LAIV_wilcoxon_youth_stringent$significant == 'no'] <- ''

Tcell_LAIV_wilcoxon_youth_stringent_upregulation <- 
  Tcell_LAIV_wilcoxon_youth_stringent[(Tcell_LAIV_wilcoxon_youth_stringent$significant == 'yes') & 
                                       (Tcell_LAIV_wilcoxon_youth_stringent$effect == '+'),]
Tcell_LAIV_wilcoxon_youth_stringent_upregulation <-
  unique(Tcell_LAIV_wilcoxon_youth_stringent_upregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','n1','n2','ENTREZID')])
Tcell_LAIV_wilcoxon_youth_stringent_upregulation$days <- 
  substring(Tcell_LAIV_wilcoxon_youth_stringent_upregulation$day_group,first = 4,last = 30)
Tcell_LAIV_wilcoxon_youth_stringent_upregulation <- Tcell_LAIV_wilcoxon_youth_stringent_upregulation %>% arrange(days,desc(log2FC_mean)) 
Tcell_LAIV_wilcoxon_youth_stringent_upregulation <- Tcell_LAIV_wilcoxon_youth_stringent_upregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(Tcell_LAIV_wilcoxon_youth_stringent_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_sTcell_lognorm_LAIV_bootstrap_stringent_toddler_vs_both.xlsx')

Tcell_LAIV_wilcoxon_youth_stringent_downregulation <- 
  Tcell_LAIV_wilcoxon_youth_stringent[(Tcell_LAIV_wilcoxon_youth_stringent$significant == 'yes') & 
                                       (Tcell_LAIV_wilcoxon_youth_stringent$effect == '-'),]
Tcell_LAIV_wilcoxon_youth_stringent_downregulation <-
  unique(Tcell_LAIV_wilcoxon_youth_stringent_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','n1','n2','ENTREZID')])
Tcell_LAIV_wilcoxon_youth_stringent_downregulation$days <- 
  substring(Tcell_LAIV_wilcoxon_youth_stringent_downregulation$day_group,first = 4,last = 30)
Tcell_LAIV_wilcoxon_youth_stringent_downregulation <- Tcell_LAIV_wilcoxon_youth_stringent_downregulation %>% arrange(days,log2FC_mean) 
Tcell_LAIV_wilcoxon_youth_stringent_downregulation <- Tcell_LAIV_wilcoxon_youth_stringent_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(Tcell_LAIV_wilcoxon_youth_stringent_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_sTcell_lognorm_LAIV_bootstrap_stringent_toddler_vs_both.xlsx',append = T)

# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp <- unique(Tcell_LAIV_wilcoxon_youth_stringent[(Tcell_LAIV_wilcoxon_youth_stringent$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$log2_log2FC_mean <- sign(temp$log2FC_mean)*log2(abs(temp$log2FC_mean) + 1)
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp, aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of Tcell expression change by LAIV\nbetween toddlers and olders on',day_group_name)) +
#     xlab("log2 fold-change of expression diff") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-FC_cutoff,FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_stringent_toddler_vs_both_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }

############# stringent in adult vs both ########################
Tcell_LAIV_wilcoxon_adult <- Tcell_LAIV_wilcoxon[(Tcell_LAIV_wilcoxon$age_group2 == '27-39yrs'),]
Tcell_LAIV_wilcoxon_adult_stringent <- Tcell_LAIV_wilcoxon_adult %>% group_by(day_group,marker) %>% 
  mutate(significant = ifelse(all(significant == 'yes') & (prod(unique(log2FC_mean)) > 0),'yes','no'))
Tcell_LAIV_wilcoxon_adult_stringent <- Tcell_LAIV_wilcoxon_adult_stringent %>% group_by(day_group,marker) %>% 
  mutate(pval_adj_mean = exp(mean(log(unique(pval_adj_mean)))),log2FC_mean = mean(unique(log2FC_mean)),n1 = mean(unique(n1)))
Tcell_LAIV_wilcoxon_adult_stringent <- Tcell_LAIV_wilcoxon_adult_stringent %>% group_by(day_group,marker) %>% 
  mutate(effect = ifelse(log2FC_mean > 0, "+", "-"))
Tcell_LAIV_wilcoxon_adult_stringent$label[Tcell_LAIV_wilcoxon_adult_stringent$significant == 'no'] <- ''

Tcell_LAIV_wilcoxon_adult_stringent_upregulation <- 
  Tcell_LAIV_wilcoxon_adult_stringent[(Tcell_LAIV_wilcoxon_adult_stringent$significant == 'yes') & 
                                        (Tcell_LAIV_wilcoxon_adult_stringent$effect == '+'),]
Tcell_LAIV_wilcoxon_adult_stringent_upregulation <-
  unique(Tcell_LAIV_wilcoxon_adult_stringent_upregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','n1','n2','ENTREZID')])
Tcell_LAIV_wilcoxon_adult_stringent_upregulation$days <- 
  substring(Tcell_LAIV_wilcoxon_adult_stringent_upregulation$day_group,first = 4,last = 30)
Tcell_LAIV_wilcoxon_adult_stringent_upregulation <- Tcell_LAIV_wilcoxon_adult_stringent_upregulation %>% arrange(days,desc(log2FC_mean)) 
Tcell_LAIV_wilcoxon_adult_stringent_upregulation <- Tcell_LAIV_wilcoxon_adult_stringent_upregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(Tcell_LAIV_wilcoxon_adult_stringent_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_sTcell_lognorm_LAIV_bootstrap_stringent_adult_vs_both.xlsx')

Tcell_LAIV_wilcoxon_adult_stringent_downregulation <- 
  Tcell_LAIV_wilcoxon_adult_stringent[(Tcell_LAIV_wilcoxon_adult_stringent$significant == 'yes') & 
                                        (Tcell_LAIV_wilcoxon_adult_stringent$effect == '-'),]
Tcell_LAIV_wilcoxon_adult_stringent_downregulation <-
  unique(Tcell_LAIV_wilcoxon_adult_stringent_downregulation[,c('day_group','marker','log2FC_mean','pval_adj_mean','n1','n2','ENTREZID')])
Tcell_LAIV_wilcoxon_adult_stringent_downregulation$days <- 
  substring(Tcell_LAIV_wilcoxon_adult_stringent_downregulation$day_group,first = 4,last = 30)
Tcell_LAIV_wilcoxon_adult_stringent_downregulation <- Tcell_LAIV_wilcoxon_adult_stringent_downregulation %>% arrange(days,log2FC_mean) 
Tcell_LAIV_wilcoxon_adult_stringent_downregulation <- Tcell_LAIV_wilcoxon_adult_stringent_downregulation %>% group_by(day_group) %>% 
  mutate(n_significant = length(unique(marker)),n_distinct = sum(abs(log2FC_mean) >= 0.25))
write.xlsx(data.frame(Tcell_LAIV_wilcoxon_adult_stringent_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_sTcell_lognorm_LAIV_bootstrap_stringent_adult_vs_both.xlsx',append = T)

# for (day_group_name in LAIV_day_group_list) {
#   print(day_group_name)
#   graphics.off()
#   temp <- unique(Tcell_LAIV_wilcoxon_adult_stringent[(Tcell_LAIV_wilcoxon_adult_stringent$day_group == day_group_name),][,c('marker','log2FC_mean','pval_adj_mean','effect','significant','label')])
#   temp$log2_log2FC_mean <- sign(temp$log2FC_mean)*log2(abs(temp$log2FC_mean) + 1)
#   temp$significant[(abs(temp$log2FC_mean) >= 0.25) & (temp$significant == 'yes')] <- 'significant and distinct'
#   temp$significant[(abs(temp$log2FC_mean) < 0.25) & (temp$significant == 'yes')] <- 'significant'
#   plot<-ggplot(temp, aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),color = significant)) +
#     geom_point(size = 1,alpha = 0.5) +
#     scale_colour_manual(values = c("significant and distinct" = "red",'significant' = 'blue', "no" = "darkgrey")) +
#     ggtitle(paste('DEG of Tcell expression change by LAIV\nbetween adults and olders on',day_group_name)) +
#     xlab("log2 fold-change of expression diff") + ylab("-log10 adjusted p-value") +
#     theme(text = element_text(size = 20)) +
#     theme_bw() +
#     geom_text_repel(aes(x=log2_log2FC_mean, y=-log10(pval_adj_mean),label = label), color = "black", size = 3) +
#     geom_vline(xintercept=c(-FC_cutoff,FC_cutoff), linetype="dotted") +
#     geom_hline(yintercept=c(-log10(p_cutoff)), linetype="dotted")#+
#   print(plot)
#   dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_stringent_adult_vs_both_',day_group_name,'.pdf',sep = ''),width = 6.5, height = 5)
# }
# 
############## venn diagram #########################
# compare BCL6+, BCL6-memory and BCL6- nonmemory nonPB 
temp1 <- read.xlsx('./Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_BCL6neg_memory_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'downregulation')
temp2 <- read.xlsx('./Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_BCL6pos_nonPB_n238_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'downregulation')
temp3 <- read.xlsx('./Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_BCL6pos_lognorm_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'downregulation')
temp <- read.xlsx('./Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n316_toddler_vs_adult.xlsx',sheetName = 'downregulation')
temp11 <- unique(temp1$marker[temp1$day_group != 'day0'])
Tcell_GC_downregulated_marker_list <- unique(temp2$marker[temp2$day_group != 'day0'])
temp33 <- unique(temp3$marker[temp3$day_group != 'day0'])
Tcell_downregulated_marker_list <- unique(temp$marker[temp$day_group != 'day0'])
# all_Tcell_sig_genes <- Tcell_downregulated_marker_list
# all_Tcell_T17_sig_genes <- temp11
# all_Tcell_Tfh_sig_genes <- temp22
# all_Tcell_non_Tfh_sig_genes <- temp33
# all_Tcell_Tfh_like_sig_genes <- temp44
overlapped_downregulation_marker_list <- Tcell_downregulated_marker_list[Tcell_downregulated_marker_list %in% unique(temp11[(temp11 %in% Tcell_GC_downregulated_marker_list) & (temp11 %in% temp33)])]

# length(unique(temp11[(temp11 %in% temp22) & (temp11 %in% temp33)]))
# length(unique(temp11[(temp11 %in% temp22) & (!temp11 %in% temp33)]))
# length(unique(temp11[(!temp11 %in% temp22) & (temp11 %in% temp33)]))
# length(unique(temp11[(!temp11 %in% temp22) & (!temp11 %in% temp33)]))
# length(unique(temp22[(!temp22 %in% temp11) & (temp22 %in% temp33)]))
# length(unique(temp22[(!temp22 %in% temp11) & (!temp22 %in% temp33)]))
# length(unique(temp33[(!temp33 %in% temp11) & (!temp33 %in% temp22)]))

temp1 <- read.xlsx('./Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_BCL6neg_memory_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'upregulation')
temp2 <- read.xlsx('./Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_BCL6pos_nonPB_n238_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'upregulation')
temp3 <- read.xlsx('./Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_BCL6pos_lognorm_LAIV_bootstrap_toddler_vs_adult.xlsx',sheetName = 'upregulation')
temp <- read.xlsx('./Tcell_age_LAIV_wilcoxon_bootstrap/tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_bootstrap_n316_toddler_vs_adult.xlsx',sheetName = 'upregulation')
temp11 <- unique(temp1$marker[temp1$day_group != 'day0'])
Tcell_GC_upregulated_marker_list <- unique(temp2$marker[temp2$day_group != 'day0'])
temp33 <- unique(temp3$marker[temp3$day_group != 'day0'])
Tcell_upregulated_marker_list <- unique(temp$marker[temp$day_group != 'day0'])
# all_Tcell_sig_genes <- Tcell_upregulated_marker_list
# all_Tcell_T17_sig_genes <- temp11
# all_Tcell_Tfh_sig_genes <- temp22
# all_Tcell_non_Tfh_sig_genes <- temp33
# all_Tcell_Tfh_like_sig_genes <- temp44
overlapped_upregulation_marker_list <- Tcell_upregulated_marker_list[Tcell_upregulated_marker_list %in% unique(temp11[(temp11 %in% temp22) & (temp11 %in% temp33)])]

marker_name <- 'IGHA1-secreted'
marker_name %in% temp11
marker_name %in% temp22
marker_name %in% temp33

temp_overlap <- unique(temp11[(temp11 %in% temp22) & (temp11 %in% temp33)])

cytokine_global_table <- read.xlsx('../Global landscape of cytokines Supplementary Table S1.xlsx', sheetName = "Cytokines")
cytokine_global_list <- unique(cytokine_global_table$HGNC.symbol)
cytokine_global_list <- c(cytokine_global_list,'GZMB','PRF1','MIF')
all_Tcell_sig_genes <- unique(all_Tcell_sig_genes)
all_Tcell_Tfh_sig_genes <- unique(all_Tcell_sig_genes)
seurat_cytokine_list <- cytokine_global_list[cytokine_global_list %in% all_Tcell_sig_genes]
seurat_cytokine_list <- all.markers[all.markers %in% cytokine_global_list]
############## visualization #################################
library(reshape2)
Tcell_GC_overlapped_upregulation_downregulation_marker_list <- Tcell_GC_downregulated_marker_list[Tcell_GC_downregulated_marker_list %in% Tcell_GC_upregulated_marker_list]
# temp <- temp[!(temp$marker %in% c('ab-TCR-gamma-delta','ab-TCR-alpha-beta','ab-CD3','ab-Tcell_GC',"ab-CD8")),]
# temp <- temp[!(temp$marker %in% c('CD38',"CXCR3","IGHG1-membrane","IL2RA","ENTPD1",'GPI',"IL2RB","PFKP","TPI1","LDHA","PGK1","ENO1")),]

Tcell_GC_LAIV_wilcoxon_youth_vs_adult_upregulation_cytokine <- 
  Tcell_GC_LAIV_wilcoxon_youth_vs_adult_upregulation[Tcell_GC_LAIV_wilcoxon_youth_vs_adult_upregulation$marker %in% Tcell_GC_overlapped_upregulation_downregulation_marker_list,]
Tcell_GC_LAIV_wilcoxon_youth_vs_adult_downregulation_cytokine <- 
  Tcell_GC_LAIV_wilcoxon_youth_vs_adult_downregulation[Tcell_GC_LAIV_wilcoxon_youth_vs_adult_downregulation$marker %in% Tcell_GC_overlapped_upregulation_downregulation_marker_list,]
# combine the upregulation and downregulation 
temp <- Tcell_GC_LAIV_wilcoxon_youth_vs_adult_upregulation_cytokine[,c('marker','days','log2FC_mean')]
temp$day_group <- paste('day',temp$days,sep = '')
temp$effect <- 'downregulation'
temp2 <- data.frame(table(temp$marker,temp$days))
colnames(temp2) <- c('marker','days','Freq')
temp2$Freq[temp2$Freq > 0] <- 1
temp2$effect <- 'up'

temp4 <- Tcell_GC_LAIV_wilcoxon_youth_vs_adult_downregulation_cytokine[,c('marker','days','log2FC_mean')]
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
# # temp2 <- temp2[!(temp2$marker %in% c('ab-TCR-gamma-delta','ab-TCR-alpha-beta','ab-CD3','ab-Tcell_GC','CD3E')),]
# ggplot(temp2, aes(x=days, y=marker,color = log2FC_mean)) + geom_point() + #(1 + )
#   RotatedAxis() + 
#   facet_wrap(~effect) +
#   # scale_colour_manual(values = c("upregulation" = "red", "downregulation" = "blue")) +
#   scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
#   ggtitle('significant DEGs in Tfh\nbetween toddlers and adults\nin LAIV subtracted by ns') #+
#   # theme(legend.position="none")
# dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_GC_LAIV_wilcoxon_toddler_vs_adult_days_cytokines.pdf',sep = ''),width = 7, height = 18)

temp2$log2FC <- abs(temp2$log2FC_mean)
ggplot(temp2, aes(x=day_group, y=marker,color = effect)) + geom_point(aes(size = log2FC)) + #(1 + )
  RotatedAxis() + 
  facet_wrap(~effect) +
  scale_colour_manual(values = c("upregulation" = "red", "downregulation" = "blue")) +
  # scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(breaks = 1:4) +
  ggtitle('significant DEGs in Tcell_GC\nin adults vs toddlers') #+
# theme(legend.position="none")
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_GC_LAIV_wilcoxon_toddler_vs_adult_days_overlap.pdf',sep = ''),width = 5, height = 6)

##### visualization one direction #################################
Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation_shrink <- 
  Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation[Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation$marker %in% 
                                                    overlapped_downregulation_marker_list,]
temp <- Tcell_LAIV_wilcoxon_youth_vs_adult_downregulation_shrink[,c('marker','days','log2FC_mean')]
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
  ggtitle('Upregulated DEGs in Tcell\nin adults than younger children')
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_bootstrap_toddler_vs_adult_days_downregulation_overlapped.pdf',sep = ''),
          width = 6, height = 11)



######
temp_table1 <- read_excel('tonsil_LAIV_120a_s_Tcell_age_wilcoxon_DAVID_upregulation.xlsx')

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
temp <- Tcell_LAIV_wilcoxon[(Tcell_LAIV_wilcoxon$pval_adj <= 0.05) & 
                                              (Tcell_LAIV_wilcoxon$age_group2 == '27-39yrs') & 
                                              (!(is.na(Tcell_LAIV_wilcoxon$pval_adj))),]
temp <- unique(temp[,c('marker','effect','n_significant')])
temp <- temp %>% arrange(effect,desc(n_significant))

############## stringent comparison ##############
# stringent
Tcell_LAIV_wilcoxon_youth_stringent <- Tcell_LAIV_wilcoxon_youth
Tcell_LAIV_wilcoxon_stringent <- Tcell_LAIV_wilcoxon_stringent %>% group_by(day_group,marker,bootstrap_index) %>% 
  mutate(significant = ifelse(all((pval_adj <= 0.05) & (prod(diff) > 0)), 'yes','no'))
Tcell_LAIV_wilcoxon_stringent <- Tcell_LAIV_wilcoxon_stringent %>% group_by(day_group,marker,effect) %>% 
  mutate(n_significant = sum(significant == 'yes')/2)
Tcell_LAIV_wilcoxon_stringent <- Tcell_LAIV_wilcoxon_stringent %>% group_by(day_group,marker,effect) %>% 
  mutate(significant = ifelse(n_significant > floor(n_bootstrap_final*0.95), "yes", "no"))
Tcell_LAIV_wilcoxon_stringent <- Tcell_LAIV_wilcoxon_stringent %>% group_by(day_group,marker) %>% 
  mutate(diff_mean = mean(diff),p_mean = exp(mean(log(pval_adj))))
Tcell_LAIV_wilcoxon_stringent$label <- Tcell_LAIV_wilcoxon_stringent$marker
Tcell_LAIV_wilcoxon_stringent$label[Tcell_LAIV_wilcoxon_stringent$significant == 'no'] <- ''

Tcell_LAIV_wilcoxon_youth_stringent_upregulation <- Tcell_LAIV_wilcoxon_stringent[(Tcell_LAIV_wilcoxon_stringent$significant == 'yes') & 
                                                                                              (Tcell_LAIV_wilcoxon_stringent$effect == '+'),]
Tcell_LAIV_wilcoxon_youth_stringent_upregulation <- unique(Tcell_LAIV_wilcoxon_youth_stringent_upregulation[,c('day_group','marker','diff_mean','p_mean','n1','n2','ENTREZID')])
Tcell_LAIV_wilcoxon_youth_stringent_upregulation <- unique(Tcell_LAIV_wilcoxon_youth_stringent_upregulation %>% group_by(day_group,marker) %>% mutate(n2 = mean(n2)))
Tcell_LAIV_wilcoxon_youth_stringent_upregulation$days <- substring(Tcell_LAIV_wilcoxon_youth_stringent_upregulation$day_group,first = 4,last = 30)
Tcell_LAIV_wilcoxon_youth_stringent_upregulation <- Tcell_LAIV_wilcoxon_youth_stringent_upregulation %>% arrange(days,p_mean) 
Tcell_LAIV_wilcoxon_youth_stringent_upregulation <- Tcell_LAIV_wilcoxon_youth_stringent_upregulation %>% group_by(day_group) %>% mutate(n = length(unique(marker)))
write.xlsx(data.frame(Tcell_LAIV_wilcoxon_youth_stringent_upregulation),sheetName = 'upregulation',
           'tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_diff_stringent_toddler_vs_both_older.xlsx')

Tcell_LAIV_wilcoxon_youth_stringent_downregulation <- Tcell_LAIV_wilcoxon_stringent[(Tcell_LAIV_wilcoxon_stringent$significant == 'yes') & 
                                                                                                (Tcell_LAIV_wilcoxon_stringent$effect == '-'),]
Tcell_LAIV_wilcoxon_youth_stringent_downregulation <- unique(Tcell_LAIV_wilcoxon_youth_stringent_downregulation[,c('day_group','marker','diff_mean','p_mean','n1','n2','ENTREZID')])
Tcell_LAIV_wilcoxon_youth_stringent_downregulation <- unique(Tcell_LAIV_wilcoxon_youth_stringent_downregulation %>% group_by(day_group,marker) %>% mutate(n2 = mean(n2)))
Tcell_LAIV_wilcoxon_youth_stringent_downregulation$days <- substring(Tcell_LAIV_wilcoxon_youth_stringent_downregulation$day_group,first = 4,last = 30)
Tcell_LAIV_wilcoxon_youth_stringent_downregulation <- Tcell_LAIV_wilcoxon_youth_stringent_downregulation %>% arrange(days,p_mean) 
Tcell_LAIV_wilcoxon_youth_stringent_downregulation <- Tcell_LAIV_wilcoxon_youth_stringent_downregulation %>% group_by(day_group) %>% mutate(n = length(unique(marker)))
write.xlsx(data.frame(Tcell_LAIV_wilcoxon_youth_stringent_downregulation),sheetName = 'downregulation',
           'tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_diff_stringent_toddler_vs_both_older.xlsx',append = T)

for (day_group_name in LAIV_day_group_list[2:13]) {
  print(day_group_name)
  graphics.off()
  temp <- unique(Tcell_LAIV_wilcoxon_stringent[(Tcell_LAIV_wilcoxon_stringent$day_group == day_group_name),][,c('marker','diff_mean','p_mean','effect','significant','label')])
  plot<-ggplot(temp, aes(x=diff_mean, y=-log10(p_mean),color = significant)) +
    geom_point(size = 1,alpha = 0.5) +
    scale_colour_manual(values = c("yes" = "red", "no" = "darkgrey")) +
    ggtitle(paste('DE of Tcell expression change by LAIV\nbetween toddlers and non-toddlers on',day_group_name)) +
    xlab("Tcell gene expression change by LAIV") + ylab("-log10 adjusted p-value") +
    theme(text = element_text(size = 20)) +
    theme_bw() +
    geom_text_repel(aes(x=diff_mean, y=-log10(p_mean),label = label), color = "black", size = 3) +
    # geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff), linetype="dotted") +
    geom_hline(yintercept=c(-log10(0.05)), linetype="dotted")#+
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_wilcoxon_diff_stringent_toddler_vs_both_older_',day_group_name,'.pdf',sep = ''),width = 5.5, height = 5)
  
}
