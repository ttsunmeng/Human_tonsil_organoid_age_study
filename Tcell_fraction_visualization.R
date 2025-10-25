library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(ggrepel)
library('stringr')
# setwd("/Volumes/GoogleDrive/My\ Drive/Stanford/RNA-seq/data_analysis/tonsil_LAIV_Rhapsody_120a-s/")

Tcell_LAIV.subset <- readRDS('tonsil_LAIV_120a_s_Tcell_lognorm_cca_LAIV.rds')
Tcell_LAIV.subset <- subset(Tcell.subset,stimulation != 'LAIV-')
Tcell_LAIV.subset <- subset(Tcell_LAIV.subset,days != 'day02')
Tcell_LAIV.subset <- subset(Tcell_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
         
# unique(Tcell_LAIV.subset$condition[Tcell_LAIV.subset$days == 'day12'])               
saveRDS(Tcell_LAIV.subset,'tonsil_LAIV_120a_s_Tcell_lognorm_cca_LAIV.rds')

LAIV_day_group_dict <- list()
LAIV_day_group_dict[['day0']] <- c(0)
LAIV_day_group_dict[['day04']] <- c(4)
LAIV_day_group_dict[['day10']] <- c(10)
LAIV_day_group_dict[['day06-07']] <- c(6,7)
LAIV_day_group_dict[['day12-14']] <- c(12,14)
LAIV_day_group_list <- names(LAIV_day_group_dict)

age_group_list <- sort(unique(Tcell_LAIV.subset$age_group))

batch_donorID_list <- list()
batch_donorID_list[[1]] <- unique(Tcell_LAIV.subset$donor_ID[Tcell_LAIV.subset$batch == 'batch1'])
batch_donorID_list[[2]] <- unique(Tcell_LAIV.subset$donor_ID[Tcell_LAIV.subset$batch == 'batch2'])
batch_donorID_list[[3]] <- unique(Tcell_LAIV.subset$donor_ID[Tcell_LAIV.subset$batch == 'batch3'])
batch_donorID_list[[4]] <- unique(Tcell_LAIV.subset$donor_ID[Tcell_LAIV.subset$batch == 'batch4'])

batch_marker_list <- list()
batch_marker_list[[1]] <- readRDS('tonsil_LAIV_120a.rds')
batch_marker_list[[2]] <- readRDS('tonsil_LAIV_120b.rds')
batch_marker_list[[3]] <- readRDS('tonsil_LAIV_120c.rds')
batch_marker_list[[4]] <- readRDS('tonsil_LAIV_120j.rds')
########### Tcell gating ###########################################
Tcell_LAIV.subset$BCL6_pos <- 'BCL6-B'
Tcell_LAIV.subset$BCL6_pos[(as.matrix(Tcell_LAIV.subset@assays[["log1p"]]['BCL6']) > 0)] <- 'BCL6+B'
Tcell_LAIV.subset$BCL6_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Tcell_LAIV.subset$AICDA_pos <- 'AICDA-B'
Tcell_LAIV.subset$AICDA_pos[(as.matrix(Tcell_LAIV.subset@assays[["log1p"]]['AICDA']) > 0)] <- 'AICDA+B'
Tcell_LAIV.subset$AICDA_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Tcell_LAIV.subset$MZB1_pos <- 'MZB1-B'
Tcell_LAIV.subset$MZB1_pos[(as.matrix(Tcell_LAIV.subset@assays[["log1p"]]['MZB1']) > 0)] <- 'MZB1+B'
Tcell_LAIV.subset$MZB1_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Tcell_LAIV.subset$AbCXCR5_pos <- 'AbCXCR5-B'
Tcell_LAIV.subset$AbCXCR5_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0)] <- 'AbCXCR5+B'
Tcell_LAIV.subset$AbCXCR5_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'PB'
Tcell_LAIV.subset$AbCXCR5_pos[Tcell_LAIV.subset$differentiation == 'memory B'] <- 'PB'

Tcell_LAIV.subset$AbCD83_pos <- 'AbCD83-B'
Tcell_LAIV.subset$AbCD83_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD83']) > 0)] <- 'AbCD83+B'
Tcell_LAIV.subset$AbCD83_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Tcell_LAIV.subset$AbCXCR4_pos <- 'AbCXCR4-B'
Tcell_LAIV.subset$AbCXCR4_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD184']) > 0)] <- 'AbCXCR4+B'
Tcell_LAIV.subset$AbCXCR4_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Tcell_LAIV.subset$AbCXCR5AbCD83_pos <- 'Tcell'
Tcell_LAIV.subset$AbCXCR5AbCD83_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD83']) > 0)] <- 'AbCXCR5+AbCD83+B'
Tcell_LAIV.subset$AbCXCR5AbCD83_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'Tcell'

Tcell_LAIV.subset$AbCXCR5AbCXCR4_pos <- 'Tcell'
Tcell_LAIV.subset$AbCXCR5AbCXCR4_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & 
                                  (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD184']) > 0)] <- 'AbCXCR5+AbCXCR4+B'
Tcell_LAIV.subset$AbCXCR5AbCXCR4_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'Tcell'

Tcell_LAIV.subset$AbCD27AbCD38_pos <- 'Tcell'
Tcell_LAIV.subset$AbCD27AbCD38_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD27']) > 0) & 
                                       (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD38']) > 0)] <- 'AbCD27+AbCD38+B'
Tcell_LAIV.subset$AbCD27AbCD38_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'Tcell'

Tcell_LAIV.subset$AbCD11cTBET_pos <- 'Tcell'
Tcell_LAIV.subset$AbCD11cTBET_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD11c']) > 0) & 
                                 (as.matrix(Tcell_LAIV.subset@assays[["log1p"]]['TBX21']) > 0)] <- 'AbCD11c+TBX21+B'
Tcell_LAIV.subset$AbCD11cTBET_pos[Tcell_LAIV.subset$differentiation %in% c('PB')] <- 'PB'

Tcell_LAIV.subset$AbCD24AbCD38_pos <- 'Tcell'
Tcell_LAIV.subset$AbCD24AbCD38_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD38']) > 0) & 
                                (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD24']) > 0)] <- 'AbCD24AbCD38+B'
Tcell_LAIV.subset$AbCD24AbCD38_pos[Tcell_LAIV.subset$differentiation == 'PB'] <- 'Tcell'


Tcell_LAIV.subset$BCL6_differentiation <- 'BCL6-B'
Tcell_LAIV.subset$BCL6_differentiation[Tcell_LAIV.subset$BCL6_pos == 'BCL6+B'] <- 'BCL6+B'
Tcell_LAIV.subset$BCL6_differentiation[(Tcell_LAIV.subset$BCL6_pos == 'BCL6-B') & (Tcell_LAIV.subset$differentiation != 'memory B')] <- 'BCL6-B'
Tcell_LAIV.subset$BCL6_differentiation[(Tcell_LAIV.subset$BCL6_pos == 'BCL6-B') & (Tcell_LAIV.subset$differentiation == 'memory B')] <- 'memory B'
Tcell_LAIV.subset$BCL6_differentiation[(Tcell_LAIV.subset$differentiation == 'PB')] <- 'PB'

Tcell_LAIV.subset$AbCXCR5_differentiation <- 'AbCXCR5-B'
Tcell_LAIV.subset$AbCXCR5_differentiation[Tcell_LAIV.subset$AbCXCR5_pos == 'AbCXCR5+B'] <- 'AbCXCR5+B'
Tcell_LAIV.subset$AbCXCR5_differentiation[(Tcell_LAIV.subset$AbCXCR5_pos == 'AbCXCR5-B') & (Tcell_LAIV.subset$differentiation != 'memory B')] <- 'AbCXCR5-B'
Tcell_LAIV.subset$AbCXCR5_differentiation[(Tcell_LAIV.subset$AbCXCR5_pos == 'AbCXCR5-B') & (Tcell_LAIV.subset$differentiation == 'memory B')] <- 'memory B'
Tcell_LAIV.subset$AbCXCR5_differentiation[(Tcell_LAIV.subset$differentiation == 'PB')] <- 'PB'

Tcell_LAIV.subset$AbCD11cTBET_BCL6_differentiation <- 'Tcells'
Tcell_LAIV.subset$AbCD11cTBET_BCL6_differentiation[Tcell_LAIV.subset$AbCD11cTBET_pos == 'AbCD11c+TBX21+B'] <- 
  paste('AbCD11c+TBX21+',Tcell_LAIV.subset$BCL6_differentiation[Tcell_LAIV.subset$AbCD11cTBET_pos == 'AbCD11c+TBX21+B'],sep = '')

Tcell_LAIV.subset$AbCD11c_AbCXCR5_differentiation <- 'Tcells'
Tcell_LAIV.subset$AbCD11c_AbCXCR5_differentiation[(Tcell_LAIV.subset$AbCXCR5_differentiation == 'AbCXCR5-B') & 
                                                    (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD11c']) > 0) & 
                                                    (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD24']) > 0) & 
                                                    (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD38']) > 0)] <- 
  'aNAV'
Tcell_LAIV.subset$AbCD11c_AbCXCR5_differentiation[(Tcell_LAIV.subset$differentiation == 'PB')] <- 'PB'

Tcell_LAIV.subset$DN2 <- 'Tcells'
Tcell_LAIV.subset$DN2[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) < 0) & 
                        (as.matrix(Tcell_LAIV.subset@assays[["log1p"]]['BCL6']) <= 0) &
                        (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD21']) < 0) & 
                        (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-IgD']) < 0) & 
                        (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD27']) < 0)] <- 
  'DN2'
Tcell_LAIV.subset$DN2[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & 
                        (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD21']) > 0) & 
                        (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-IgD']) < 0) & 
                        (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD27']) < 0)] <- 
  'DN1'
Tcell_LAIV.subset$DN2[(Tcell_LAIV.subset$differentiation == 'PB')] <- 'PB'

Tcell_LAIV.subset$MKI67_BCL6_differentiation <- 
  paste(ifelse(as.matrix(Tcell_LAIV.subset@assays[["log1p"]]['MKI67']) > 0,'MKI67+','MKI67-'),Tcell_LAIV.subset$BCL6_differentiation, sep = '')

Tcell_LAIV.subset$MZB <- 'Tcells'
Tcell_LAIV.subset$MZB[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-IgD']) > 0) & 
                        (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-IgM']) > 0) & 
                        (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD27']) > 0)] <- 
  'MZB'
Tcell_LAIV.subset$MZB[(Tcell_LAIV.subset$differentiation == 'PB')] <- 'Tcells'

Tcell_LAIV.subset$AbCD11cAbCD21_pos <- 'Tcell'
Tcell_LAIV.subset$AbCD11cAbCD21_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD11c']) > 0) & 
                                    (as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD21']) < 0)] <- 'AbCD11c+AbCD21-B'
Tcell_LAIV.subset$AbCD11cAbCD21_pos[Tcell_LAIV.subset$differentiation %in% c('PB')] <- 'Tcell'

PB_secreted <- data.frame(t(as.matrix(Tcell_LAIV.subset@assays$raw@data)),check.names = F)[,rownames(Tcell_LAIV.subset)[grepl('-secreted',rownames(Tcell_LAIV.subset))]]
PB_secreted$max <- colnames(PB_secreted)[apply(PB_secreted,1,which.max)]
PB_secreted$isotype <- PB_secreted$max
PB_secreted$isotype[grepl('IGHG',PB_secreted$max)] <- 'IgG+PB'
PB_secreted$isotype[grepl('IGHM',PB_secreted$max)] <- 'IgM+PB'
PB_secreted$isotype[grepl('IGHA',PB_secreted$max)] <- 'IgA+PB'
PB_secreted$isotype[grepl('IGHE',PB_secreted$max)] <- 'IgE+PB'
Tcell_LAIV.subset$isotype <- PB_secreted$isotype
Tcell_LAIV.subset$isotype_secreted <- PB_secreted$max
Tcell_LAIV.subset$isotype[Tcell_LAIV.subset$differentiation != 'PB'] <- 'B cells'
Tcell_LAIV.subset$isotype_secreted[Tcell_LAIV.subset$differentiation != 'PB'] <- 'B cells'

PB_secreted <- data.frame(t(as.matrix(Tcell_LAIV.subset@assays$raw@data)),check.names = F)[,rownames(Tcell_LAIV.subset)[grepl('-secreted',rownames(Tcell_LAIV.subset)) | grepl('-membrane',rownames(Tcell_LAIV.subset))]]
PB_secreted$max <- colnames(PB_secreted)[apply(PB_secreted,1,which.max)]
Tcell_LAIV.subset$isotype_all <- PB_secreted$max

Tcell_LAIV.subset$isotype_both <- Tcell_LAIV.subset$isotype_all
Tcell_LAIV.subset$isotype_both <- gsub('-.*','+',Tcell_LAIV.subset$isotype_both)
Tcell_LAIV.subset$isotype_both[Tcell_LAIV.subset$differentiation == 'PB'] <- paste(Tcell_LAIV.subset$isotype_both[Tcell_LAIV.subset$differentiation == 'PB'],'PB',sep = '')
Tcell_LAIV.subset$isotype_both[Tcell_LAIV.subset$differentiation != 'PB'] <- paste(Tcell_LAIV.subset$isotype_both[Tcell_LAIV.subset$differentiation != 'PB'],'B',sep = '')

# Tcell_LAIV.subset$isotype_main_both <- Tcell_LAIV.subset$isotype_all
# Tcell_LAIV.subset$isotype_main_both <- gsub('-.*','+',Tcell_LAIV.subset$isotype_both)
# Tcell_LAIV.subset$isotype_main_both[Tcell_LAIV.subset$differentiation == 'PB'] <- paste(Tcell_LAIV.subset$isotype_both[Tcell_LAIV.subset$differentiation == 'PB'],'PB',sep = '')
# Tcell_LAIV.subset$isotype_main_both[Tcell_LAIV.subset$differentiation != 'PB'] <- paste(Tcell_LAIV.subset$isotype_both[Tcell_LAIV.subset$differentiation != 'PB'],'B',sep = '')

Tcell_LAIV.subset$isotype_cluster <- Tcell_LAIV.subset$isotype_all
Tcell_LAIV.subset$isotype_cluster <- gsub('-.*','max_',Tcell_LAIV.subset$isotype_cluster)
Tcell_LAIV.subset$isotype_cluster <- paste(Tcell_LAIV.subset$isotype_cluster,Tcell_LAIV.subset$CD27CD38_cluster,sep = '')

Tcell_LAIV.subset$AbCXCR3_pos <- 'Tcell'
Tcell_LAIV.subset$AbCXCR3_pos[(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD183']) > 0) &
                                      Tcell_LAIV.subset$BCL6_differentiation == 'BCL6-B'] <- 'AbCXCR3+BCL6-B'
Tcell_LAIV.subset$AbCXCR3_pos[Tcell_LAIV.subset$differentiation %in% c('PB')] <- 'Tcell'

Tcell_LAIV.subset$IGHG2_pos <- 'Tcell'
Tcell_LAIV.subset$IGHG2_pos[(as.matrix(Tcell_LAIV.subset@assays[["log1p"]]['IGHG2-secreted']) > 5) & 
                              (!(Tcell_LAIV.subset$BCL6_differentiation %in% c('PB','memory B')))] <- 'IGHG2+B'
Tcell_LAIV.subset$IGHG2_pos[(as.matrix(Tcell_LAIV.subset@assays[["log1p"]]['IGHG2-secreted']) > 5) & 
                              (Tcell_LAIV.subset$BCL6_differentiation %in% c('memory B'))] <- 'IGHG2+memory B'


DimPlot(Tcell_LAIV.subset, label = T, reduction = "umap",group.by = 'isotype', order = F,pt.size = 0.1)
dev.print(png, paste('tonsil_LAIV_120a_s_Tcell_umap_AbCD11cTBET_pos.png',sep = ''),width = 600, height = 400)

Tcell_log1p_dataframe <- data.frame(t(as.matrix(Tcell_LAIV.subset@assays[["log1p"]]@data)),check.names = F)
Tcell_integrated_scale_dataframe <- data.frame(t(as.matrix(Tcell_LAIV.subset@assays[["integrated_scale"]]@data)),check.names = F)

Tcell_integrated_scale_dataframe$AbTcell5RA_AbTcell5RO_pos <- Tcell_LAIV.subset$AbTcell5RA_AbTcell5RO_pos
ggplot(Tcell_integrated_scale_dataframe, aes(x = !!sym(c('ab-CXCR5')), y = !!sym(c('ab-CD279')))) + ggtitle('Tcell cells') +#!!sym(c('TBX21'))
  geom_point(alpha = 0.1,size = 0.5,aes(color = !!sym(c('AbTcell5RA_AbTcell5RO_pos')))) +
  # geom_point(alpha = 0.1,size = 0.5) +
  # scale_colour_manual(values = c("naive Tcell" = "#F8766D", "Tcell" = '#00BFC4'))+
  # scale_colour_manual(values = c("Tcell5RA+CD83+ nonPB B" = "gold", "CD184-CD83- nonPB B" = "blue","CD184+CD83- nonPB B" = 'cyan',"CD184-CD83+ nonPB B" = 'magenta')) +
  geom_density_2d(color = 'black')# +
# ylim(0,8) + xlim(0,8)
# geom_hline(aes(yintercept=3.25,color = 'red')) + guides(color = 'none') +
# geom_vline(aes(xintercept=4.1,color = 'red')) + guides(color = 'none')
dev.print(pdf, 'tonsil_LAIV_120a_s_Tcell_AbCXCR5_AbCD279_contour_integrated_scale.pdf',width = 7, height = 4.6)

DefaultAssay(Tcell_LAIV.subset) <- 'integrated_scale'
FeatureScatter(Tcell_LAIV.subset, feature1 = "ab-Tcell5RA", feature2 = "ab-Tcell5RO") 
dev.print(png, 'tonsil_LAIV_120a_s_Tcell_scatter_abTcell5RA_abTcell5RO.png',width = 700, height = 434)

##### all Tcells: cell cluster fraction visualization #####################
flow_count_table <- read.xlsx('tonsil_LAIV_120a_s_flow_cell_count.xlsx', sheetName = "Sheet1")
flow_count_table <- flow_count_table[,!colnames(flow_count_table) == '...1']
flow_count_table$condition <- gsub('day0 LAIV-','day0',flow_count_table$condition)
flow_count_table$B_cells[grepl('IMD170',flow_count_table$donor_ID)] <- 0

Tcell_LAIV.subset$cluster_labeled2 <- Tcell_LAIV.subset$cluster_labeled
Tcell_LAIV.subset$cluster_labeled2[Tcell_LAIV.subset$cluster_labeled2 == 'early act. naive CD4'] <- 'naive CD4'

Tcell.database <- data.frame(cbind(Tcell_LAIV.subset$donor_ID,as.character(Tcell_LAIV.subset$cluster_labeled2),Tcell_LAIV.subset$condition,Tcell_LAIV.subset$days,Tcell_LAIV.subset$stimulation,Tcell_LAIV.subset$age_group,Tcell_LAIV.subset$day_group,Tcell_LAIV.subset$age))
colnames(Tcell.database) <- c('donor_ID','clusters','condition','days','stimulation','age_group','day_group','age')
# sapply(Tcell.database,class)
Tcell_cluster_condition_donor_count <- dplyr::count(Tcell.database, donor_ID, clusters,condition, days,stimulation,age_group,day_group,age)
Tcell_cluster_condition_donor_count <- Tcell_cluster_condition_donor_count %>% group_by(donor_ID, condition, days,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
Tcell_cluster_condition_donor_count$days <- as.numeric(substring(Tcell_cluster_condition_donor_count$days,first = 4,last = 5))
Tcell_cluster_condition_donor_count$percentage <- Tcell_cluster_condition_donor_count$n/Tcell_cluster_condition_donor_count$total*100
Tcell_cluster_condition_donor_count$age <- as.numeric(Tcell_cluster_condition_donor_count$age)
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(Tcell_LAIV.subset$donor_ID)){
  temp_condition_list <- sort(unique(Tcell_LAIV.subset$condition[Tcell_LAIV.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- Tcell_cluster_condition_donor_count[(Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Tcell_LAIV.subset$cluster_labeled2)){
      if_row <- ((Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$clusters == cluster_name) & (Tcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount$clusters <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        Tcell_cluster_condition_donor_count[nrow(Tcell_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(Tcell.database)

p <- match(Tcell_cluster_condition_donor_count$condition, flow_count_table$condition)
temp <- flow_count_table[p,]
Tcell_cluster_condition_donor_count$Tcell_count <- temp$B_cells
Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '02yrs M IMD085 day0'] <- 
  Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '02yrs M IMD085 day0']*6/10
Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '03yrs F IMD035 day0'] <- 
  Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '03yrs F IMD035 day0']*6/10
Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '07yrs F IMD150 day0'] <- 
  Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '07yrs F IMD150 day0']*6/12
Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '09yrs M IMD107 day0'] <- 
  Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '09yrs M IMD107 day0']*6/12
Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '27yrs F VIP031 day0'] <- 
  Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '27yrs F VIP031 day0']*6/10
Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '39yrs M VIP015 day0'] <- 
  Tcell_cluster_condition_donor_count$Tcell_count[Tcell_cluster_condition_donor_count$condition == '39yrs M VIP015 day0']*6/12
Tcell_cluster_condition_donor_count$cluster_count <- Tcell_cluster_condition_donor_count$Tcell_count*Tcell_cluster_condition_donor_count$percentage/100
# write.csv(Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$clusters == 'CD4',],'tonsil_LAIV_120a_s_CD4_flow_cell_count.csv')
# saveRDS(Tcell_cluster_condition_donor_count,'tonsil_LAIV_120a_s_Tcell_LAIV_log1p_cca_CycleRegressOut_cluster_labeled2_fraction.rds')
annotate_cluster_list <- unique(Tcell_LAIV.subset$cluster_labeled2)#'AbCD278+B'#
# annotate_cluster_list <- annotate_cluster_list[grepl('MKI67\\+',annotate_cluster_list)]# | grepl('CD278\\+PD1\\-CD4',annotate_cluster_list) | grepl('CD278\\-PD1\\+CD4',annotate_cluster_list)]
annotate_cluster_list <- setdiff(annotate_cluster_list,c('PB','B cells'))
for (cluster_name in annotate_cluster_list)
{
  graphics.off()
  print(cluster_name)
  temp_Tcell_cluster <- Tcell_cluster_condition_donor_count[(Tcell_cluster_condition_donor_count$clusters == cluster_name),]
  # plot <- ggplot(temp_Tcell_cluster,aes(x=days, y=percentage,color = age_group)) +  ggtitle(cluster_name) +
  #   facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
  #   geom_point(size = 1)+ geom_line(aes(group = donor_ID)) + 
  #   scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
  #   theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  #   ylab('Tcell %') +
  #   # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  #   theme_bw()
  # print(plot)
  # dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_lognorm_LAIV_age_group_cluster_labeled2_',cluster_name,'_condition_ratio2.pdf',sep = ''),width = 5, height = 4)
  plot <- ggplot(temp_Tcell_cluster,aes(x=days, y=percentage,color = age_group,fill = age_group,shape = age_group, linetype = age_group)) +  ggtitle(cluster_name) +
    geom_point(size = 1.5) + #,aes(shape = donor_ID)
    geom_smooth(method = "loess",se = F, linewidth=0.5) +
    ylim(0,max(temp_Tcell_cluster$percentage)) +
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('% in T cells') +
    scale_shape_manual(values = c('02-03yrs' = 25,'07-09yrs'= 23,'27-39yrs' = 22)) +
    scale_linetype_manual(values = c('02-03yrs' = 'dotted','07-09yrs'= 'longdash','27-39yrs' = 'solid')) +
    scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
    # geom_rect(aes(xmin=4, xmax=10, ymin=0, ymax=Inf),alpha = 0.1,fill = 'gray',color = 'gray') +
    theme_bw()
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_age_group_cluster_labeled2_',gsub('\\/','or',cluster_name),'_condition_ratio.pdf',sep = ''),width = 2.6, height = 1.5)
  temp_Tcell_cluster <- temp_Tcell_cluster[temp_Tcell_cluster$donor_ID != '07yrs M IMD170',]
  # plot <- ggplot(temp_Tcell_cluster,aes(x=days, y=cluster_count,color = age_group)) +  ggtitle(cluster_name) +
  #   facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
  #   geom_point(size = 1)+ geom_line(aes(group = donor_ID)) + 
  #   scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
  #   theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  #   ylab('#Tcells') +
  #   # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  #   theme_bw()
  # print(plot)
  # dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_age_group_cluster_labeled2_',cluster_name,'_condition_count2.pdf',sep = ''),width = 5, height = 4)
  plot <- ggplot(temp_Tcell_cluster,aes(x=days, y=cluster_count,color = age_group,fill = age_group,shape = age_group)) +  ggtitle(cluster_name) +
    geom_point(size = 2) + #,aes(shape = donor_ID)
    geom_smooth(method = "loess",se = F) +
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('#cells') +
    ylim(0,max(temp_Tcell_cluster$cluster_count)) +
    scale_shape_manual(values = c('02-03yrs' = 25,'07-09yrs'= 23,'27-39yrs' = 22)) +
    scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
    # geom_rect(aes(xmin=4, xmax=10, ymin=0, ymax=Inf),alpha = 0.1,fill = 'gray',color = 'gray') +
    theme_bw()
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_age_group_cluster_labeled2_',gsub('\\/','or',cluster_name),'_condition_count.pdf',sep = ''),width = 3, height = 2)
}
#### split by age ################ 
# mean of each cell count
library(purrr) 
library(ggalluvial)
temp <- Tcell_cluster_condition_donor_count
temp <- temp[temp$day_group %in% c('day0','day04','day06-07','day10','day12-14'),]
temp <- temp %>% group_by(age_group,clusters,day_group) %>% mutate(avg_cluster_count = mean(cluster_count),avg_fraction = mean(percentage))
temp <- unique(temp[,c('age_group','clusters','day_group','avg_cluster_count','avg_fraction')])
temp$days <- 0
temp$days[temp$day_group == 'day04'] <- 4
temp$days[temp$day_group == 'day06-07'] <- 7
temp$days[temp$day_group == 'day10'] <- 10
temp$days[temp$day_group == 'day12-14'] <- 14
# Idents(TB_shrink.seurat_update) <- TB_shrink.seurat_update@meta.data$sample_tag
temp$clusters <- factor(temp$clusters,levels = ordered_cluster_annotated_list)
ggplot(temp, aes(x = days, stratum = age_group, alluvium = clusters, y = avg_cluster_count, label = clusters,fill = clusters)) +#, 
  geom_alluvium(aes(fill = clusters), width = 1,alpha = 1,color = 'black',size = 0.1) +
  facet_grid(age_group ~ .,scales = "fixed") +
  # geom_stratum(width = 1/12, fill = "black",size = 0.5) +
  scale_fill_manual(values = cols_group) +
  ggtitle('') + ylab('#cells') + #ylim(0,300) + xlim(0,14)
theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_age_group_cluster_labeled_count.pdf',sep = ''),width = 5, height = 5)

ggplot(temp, aes(x = days, stratum = age_group, alluvium = clusters, y = avg_fraction, label = clusters,fill = clusters)) +#, 
  geom_alluvium(aes(fill = clusters), width = 0.5,alpha = 1,color = 'black',size = 0.1) +
  facet_grid(age_group ~ .,scales = "fixed") +
  # geom_stratum(width = 1/12, fill = "black") +
  scale_fill_manual(values = cols_group) +
  ggtitle('') + ylab('% in T cells') + #ylim(0,300) + xlim(0,14)
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_age_group_cluster_labeled_fraction.pdf',sep = ''),width = 5, height = 5)

######## fraction: compare each age_group pairs ###########################
Tcell_LAIV_cluster_ttest_age_group_paired <- expand.grid(clusters = annotate_cluster_list, age_group1 = age_group_list[1], age_group2 = age_group_list[3],
                                                  day_group = LAIV_day_group_list)
Tcell_LAIV_cluster_ttest_age_group_paired <- transform(Tcell_LAIV_cluster_ttest_age_group_paired, age_group1 = as.character(age_group1),
                                                age_group2 = as.character(age_group2))
Tcell_LAIV_cluster_ttest_age_group_paired <- Tcell_LAIV_cluster_ttest_age_group_paired[Tcell_LAIV_cluster_ttest_age_group_paired$age_group1 != Tcell_LAIV_cluster_ttest_age_group_paired$age_group2,]
Tcell_LAIV_cluster_ttest_age_group_paired$diff <- NaN
Tcell_LAIV_cluster_ttest_age_group_paired$pval <- NaN
Tcell_LAIV_cluster_ttest_age_group_paired$n1 <- NaN
Tcell_LAIV_cluster_ttest_age_group_paired$n2 <- NaN
Tcell_LAIV_cluster_ttest_age_group_paired$mean1 <- NaN
Tcell_LAIV_cluster_ttest_age_group_paired$mean2 <- NaN

age_pair_list <- data.frame(age_group1 = c('02-03yrs','02-03yrs','07-09yrs'),
                            age_group2 = c('07-09yrs','27-39yrs','27-39yrs'))

library(dplyr)
for (cluster_name in annotate_cluster_list) {
  print(cluster_name)
  graphics.off()
  temp <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$clusters == cluster_name,]
  rownames(temp) <- temp$condition
  # if ((!all(temp$percentage == 0)) & (!all(temp$percentage == 100))) {
  for (day_group_name in LAIV_day_group_list) {
    temp1 <- temp[temp$days %in% LAIV_day_group_dict[[day_group_name]],]
    if (!all(temp1$percentage == 0)) {
      for (age_group_index in c(2)) {
        temp_age_group1 <- age_pair_list$age_group1[age_group_index]
        temp_age_group2 <- age_pair_list$age_group2[age_group_index]
        ttest_result <- t.test(temp1$percentage[temp1$age_group == temp_age_group1],
                               temp1$percentage[temp1$age_group == temp_age_group2])
        row_select <- (Tcell_LAIV_cluster_ttest_age_group_paired$clusters == cluster_name) &
          (Tcell_LAIV_cluster_ttest_age_group_paired$age_group1 == temp_age_group1) &
          (Tcell_LAIV_cluster_ttest_age_group_paired$age_group2 == temp_age_group2) &
          (Tcell_LAIV_cluster_ttest_age_group_paired$day_group == day_group_name)
        Tcell_LAIV_cluster_ttest_age_group_paired$n1[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == temp_age_group1]))
        Tcell_LAIV_cluster_ttest_age_group_paired$n2[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == temp_age_group2]))
        Tcell_LAIV_cluster_ttest_age_group_paired$mean1[row_select] <- ttest_result[["estimate"]][["mean of x"]]
        Tcell_LAIV_cluster_ttest_age_group_paired$mean2[row_select] <- ttest_result[["estimate"]][["mean of y"]]
        Tcell_LAIV_cluster_ttest_age_group_paired$pval[row_select] <- ttest_result$p.value
        Tcell_LAIV_cluster_ttest_age_group_paired$sd[row_select] <- sd(temp1$percentage[(temp1$age_group == temp_age_group1) |
                                                                                          (temp1$age_group == temp_age_group2)])
      }
    }

    # }

  }
}

Tcell_LAIV_cluster_ttest_age_group_paired <- Tcell_LAIV_cluster_ttest_age_group_paired %>% arrange(pval)
Tcell_LAIV_cluster_ttest_age_group_paired$diff <- Tcell_LAIV_cluster_ttest_age_group_paired$mean1 - Tcell_LAIV_cluster_ttest_age_group_paired$mean2

Tcell_LAIV_cluster_ttest_age_group_paired$fdr <- p.adjust(Tcell_LAIV_cluster_ttest_age_group_paired$pval,method = 'fdr')
Tcell_LAIV_cluster_ttest_age_group_paired <- Tcell_LAIV_cluster_ttest_age_group_paired %>% arrange(fdr)
Tcell_LAIV_cluster_ttest_age_group_paired$FC <- Tcell_LAIV_cluster_ttest_age_group_paired$diff/Tcell_LAIV_cluster_ttest_age_group_paired$sd
Tcell_LAIV_cluster_ttest_age_group_paired <- Tcell_LAIV_cluster_ttest_age_group_paired[Tcell_LAIV_cluster_ttest_age_group_paired$cluster != 'PB',]
# write.xlsx(Tcell_LAIV_cluster_ttest_age_group_paired,'tonsil_LAIV_120a_s_Tcell_log1p_cca_LAIV_age_cluster_labeled_age_group_ttest_paired.xlsx')


###### prepare for the correlation network ############################
library(reshape2)
temp <- Tcell_cluster_condition_donor_count[!((Tcell_cluster_condition_donor_count$donor_ID == '02yrs F IMD030') & (Tcell_cluster_condition_donor_count$days == 12)),]
temp <- temp[!((temp$donor_ID == '33yrs M VIP024') & (temp$days == 12)),]
temp <- temp[temp$donor_ID != '02yrs M IMD085',]
temp <- dcast(temp[,c('clusters','day_group','donor_ID','percentage')], clusters + day_group ~ donor_ID, value.var="percentage")
temp <- temp[temp$day_group %in% LAIV_day_group_list,]

temp_final <- merge(temp,Tcell_LAIV_cluster_ttest_age_group_paired,by = c('clusters','day_group'))
temp_final <- temp_final[(!is.na(temp_final$n1)) & (!is.na(temp_final$n2)) & (!is.na(temp_final$pval)),]
temp_final <- temp_final[rowSums(temp_final[,3:11] == 0) < 9,]
temp_final$cluster <- as.character(temp_final$cluster)
temp_final$fdr <- p.adjust(temp_final$pval,method = 'fdr')
temp_final <- temp_final %>% arrange(fdr)
temp_final$sig <- 'no'
temp_final$sig[temp_final$pval <= 0.05] <- 'yes'
rownames(temp_final) <- paste(temp_final$day_group,temp_final$clusters,sep = '_')

write.xlsx(temp_final,'tonsil_LAIV_120a_s_Tcell_log1p_cca_CycleRegressOut_cluster_labeled_LAIV_fraction_cleaned_CytoScape.xlsx')


library(igraph)
library(RCy3)
cor_mat <- cor(t(temp_final[,3:11]))

cor_vec <- c(cor_mat)
cor_vec <- cor_vec[cor_vec != 1]
graphics.off()
hist(cor_vec, breaks = 10)
cutoff <- 0.9
abline(v=cutoff, col = 'red')# cutoff around 
abline(v=-cutoff, col = 'red')# cutoff around 

adj_mat <- cor_mat
adj_mat[(adj_mat < cutoff) & (adj_mat > -cutoff)] <- 0

# color by the FC of age, size by the connectivity
graphics.off()

network <- graph_from_adjacency_matrix(adj_mat, weighted=T, mode="undirected", diag=F)
V(network)$degree <- degree(network)

createNetworkFromIgraph(network)

######## fraction: all cell FC diff visualization ###########################
Tcell_LAIV_cluster_ttest_age_group <- expand.grid(cluster = annotate_cluster_list, age_group = age_group_list, day_group = LAIV_day_group_list)
Tcell_LAIV_cluster_ttest_age_group$diff <- NaN
Tcell_LAIV_cluster_ttest_age_group$pval <- NaN
Tcell_LAIV_cluster_ttest_age_group$n1 <- NaN
Tcell_LAIV_cluster_ttest_age_group$n2 <- NaN
Tcell_LAIV_cluster_ttest_age_group$mean1 <- NaN
Tcell_LAIV_cluster_ttest_age_group$mean2 <- NaN
Tcell_LAIV_cluster_correlate_age <- expand.grid(cluster = annotate_cluster_list, day_group = LAIV_day_group_list)
Tcell_LAIV_cluster_correlate_age$r <- NaN
Tcell_LAIV_cluster_correlate_age$pval <- NaN
Tcell_LAIV_cluster_correlate_age$n <- NaN

library(dplyr)
for (cluster_name in annotate_cluster_list) {
  print(cluster_name)
  graphics.off()
  temp <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$clusters == cluster_name,]
  rownames(temp) <- temp$condition
  # if ((!all(temp$percentage == 0)) & (!all(temp$percentage == 100))) {
  for (day_group_name in LAIV_day_group_list) {
    temp1 <- temp[temp$days %in% LAIV_day_group_dict[[day_group_name]],]
    if (!all(temp1$percentage == 0)) {
      for (age_group_index in c(1:3)) {
        age_group_name <- age_group_list[age_group_index]
        ttest_result <- t.test(temp1$percentage[temp1$age_group == age_group_name],
                               temp1$percentage[temp1$age_group != age_group_name])
        row_select <- (Tcell_LAIV_cluster_ttest_age_group$cluster == cluster_name) & 
          (Tcell_LAIV_cluster_ttest_age_group$age_group == age_group_list[age_group_index]) &
          (Tcell_LAIV_cluster_ttest_age_group$day_group == day_group_name)
        Tcell_LAIV_cluster_ttest_age_group$n1[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == age_group_name]))
        Tcell_LAIV_cluster_ttest_age_group$n2[row_select] <- length(unique(temp1$donor_ID[temp1$age_group != age_group_name]))
        Tcell_LAIV_cluster_ttest_age_group$mean1[row_select] <- ttest_result[["estimate"]][["mean of x"]]
        Tcell_LAIV_cluster_ttest_age_group$mean2[row_select] <- ttest_result[["estimate"]][["mean of y"]]
        Tcell_LAIV_cluster_ttest_age_group$pval[row_select] <- ttest_result$p.value
      }
      row_select <- (Tcell_LAIV_cluster_correlate_age$cluster == cluster_name) &
        (Tcell_LAIV_cluster_correlate_age$day_group == day_group_name)
      cor_test <- cor.test(temp1$percentage, temp1$age)
      Tcell_LAIV_cluster_correlate_age$n[row_select] <- length(unique(temp$donor_ID))
      Tcell_LAIV_cluster_correlate_age$r[row_select] <- cor_test[["estimate"]][["cor"]]
      Tcell_LAIV_cluster_correlate_age$pval[row_select] <- cor_test$p.value  
    }  
    
    # }
    
  }
}

Tcell_LAIV_cluster_ttest_age_group <- Tcell_LAIV_cluster_ttest_age_group %>% arrange(pval)
Tcell_LAIV_cluster_correlate_age <- Tcell_LAIV_cluster_correlate_age %>% arrange(pval)
Tcell_LAIV_cluster_ttest_age_group$diff <- Tcell_LAIV_cluster_ttest_age_group$mean1 - Tcell_LAIV_cluster_ttest_age_group$mean2

# Tcell_LAIV_cluster_ttest_age_group$pval_adj <- p.adjust(Tcell_LAIV_cluster_ttest_age_group$pval,method = 'fdr')
# Tcell_LAIV_cluster_correlate_age$pval_adj <- p.adjust(Tcell_LAIV_cluster_correlate_age$pval,method = 'fdr')
Tcell_LAIV_cluster_ttest_age_group <- Tcell_LAIV_cluster_ttest_age_group %>% arrange(cluster,pval)
Tcell_LAIV_cluster_correlate_age <- Tcell_LAIV_cluster_correlate_age %>% arrange(cluster,desc(r))

Tcell_LAIV_cluster_ttest_age_group <- Tcell_LAIV_cluster_ttest_age_group[Tcell_LAIV_cluster_ttest_age_group$age_group != '07-09yrs',]
Tcell_LAIV_cluster_correlate_age <- Tcell_LAIV_cluster_correlate_age[Tcell_LAIV_cluster_correlate_age$cluster != 'PB',]

write.xlsx(Tcell_LAIV_cluster_ttest_age_group,'tonsil_LAIV_120a_s_Tcell_lognorm_LAIV_age_cluster_labeled_age_group_ttest.xlsx')
write.xlsx(Tcell_LAIV_cluster_correlate_age,'tonsil_LAIV_120a_s_Tcell_lognorm_LAIV_age_cluster_labeled_age_corr.xlsx')


####### fraction: visualization ####################################
cluster_name <- 'IgG+PB'
day_group_list <- c(0,4,10)
day_group_name <- 'day0'
age_group_name <- '27-39yrs'
marker_name <- 'CXCR5'
x_lab <- 2.2
plot_corr <- 0
graphics.off()
temp_donorID_batch_list <- unlist(batch_donorID_list[which(sapply(batch_marker_list, function(x) (marker_name %in% x)))])
for (day_group_name in setdiff(LAIV_day_group_list,c('day07-08','day07'))) {
  Tcell_cluster_condition_donor_count$day_group[Tcell_cluster_condition_donor_count$days %in% LAIV_day_group_dict[[day_group_name]]] <-
    day_group_name
}

for (cluster_name in cluster_name){#annotate_cluster_list){
  print(cluster_name)
  temp <- Tcell_cluster_condition_donor_count[(Tcell_cluster_condition_donor_count$clusters == cluster_name) &
                                                # (Tcell_cluster_condition_donor_count$days %in% LAIV_day_group_dict[[day_group_name]]) &
                                                (Tcell_cluster_condition_donor_count$day_group %in% setdiff(LAIV_day_group_list,c('day07-08','day07'))) &
                                                Tcell_cluster_condition_donor_count$donor_ID %in% temp_donorID_batch_list,]
  
  plot <- ggplot(temp, aes(y = percentage, x=age_group)) +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter() +
    facet_wrap( ~ day_group, scales = "fixed",ncol = 5) +
    geom_bar(position = 'dodge', stat = 'summary', fun = mean,fill = 'white',color = 'black',width = 0.5) +
    geom_jitter(width = 0.2,size = 1) +
    # geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height=0.4, dodge.width=0.9)) +
    geom_hline(yintercept = 0,linetype=2) + 
    ylab('% B cells') +
    ggtitle(cluster_name) +
    theme_bw() +
    ylim(0,max(temp$percentage) + 0.5*abs(max(temp$percentage))) +
    RotatedAxis()
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_cluster_labeled_age_group_',cluster_name,'.pdf',sep = ''),width = 6, height = 3)
  
}

ggplot(temp, aes(y=percentage,x=age)) +
  facet_wrap( ~ day_group, scales = "fixed",ncol = 5) +
  geom_point() + RotatedAxis() +
  ylab('% B cells') +
  geom_smooth(method='lm', formula= y~x, color = 'red') +
  ggtitle(cluster_name) +
  theme_bw()#  +
  # annotate("text", x=(temp$age[which.max(temp$percentage)] + 25)%%40, 
  #          y=mean(sort(temp$percentage)[(length(temp$percentage) - 1):length(temp$percentage)]),
  #          label = paste("r =",sprintf(cor_test$estimate, fmt = '%#.2f')),size = 6)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_cluster_labeled_age_corr_',cluster_name,'.pdf',sep = ''),width = 6, height = 3)

# ttest_result <- t.test(temp$percentage[temp$age_group == age_group_name],
#                        temp$percentage[temp$age_group != age_group_name])

cor_test <- cor.test(temp$percentage,temp$age)
ggplot(temp, aes(y=percentage,x=age)) +
  geom_point() + RotatedAxis() +
  ylab('diff in % B cells\nbetween LAIV+ and LAIV-') +
  geom_smooth(method='lm', formula= y~x, color = 'red') +
  ggtitle(paste(cluster_name, day_group_name)) +
  theme_bw() +
  annotate("text", x=(temp$age[which.max(temp$percentage)] + 25)%%40, 
           y=mean(sort(temp$percentage)[(length(temp$percentage) - 1):length(temp$percentage)]),
           label = paste("r =",sprintf(cor_test$estimate, fmt = '%#.2f')),size = 6)

dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_cluster_labeled_age_corr_',cluster_name,'_',day_group_name,'.pdf',sep = ''),width = 3, height = 3)


######## count: compare each age_group pairs ###########################
p <- match(Tcell_cluster_condition_donor_count$condition, flow_count_table$condition)
temp <- flow_count_table[p,]
Tcell_cluster_condition_donor_count$Tcell_count <- temp$B_cells
Tcell_cluster_condition_donor_count$cluster_count <- Tcell_cluster_condition_donor_count$percentage*Tcell_cluster_condition_donor_count$Tcell_count/100
Tcell_cluster_condition_donor_count_cleaned <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$donor_ID != '07yrs M IMD170',]
Tcell_LAIV_cluster_ttest_age_group_paired <- expand.grid(cluster = annotate_cluster_list, age_group1 = age_group_list[1:2], age_group2 = age_group_list[2:3],
                                                         day_group = LAIV_day_group_list)
Tcell_LAIV_cluster_ttest_age_group_paired <- transform(Tcell_LAIV_cluster_ttest_age_group_paired, age_group1 = as.character(age_group1),
                                                       age_group2 = as.character(age_group2))
Tcell_LAIV_cluster_ttest_age_group_paired <- Tcell_LAIV_cluster_ttest_age_group_paired[Tcell_LAIV_cluster_ttest_age_group_paired$age_group1 != Tcell_LAIV_cluster_ttest_age_group_paired$age_group2,]
# Tcell_LAIV_cluster_ttest_age_group_paired$diff <- NaN
# Tcell_LAIV_cluster_ttest_age_group_paired$pval <- NaN
# Tcell_LAIV_cluster_ttest_age_group_paired$n1 <- NaN
# Tcell_LAIV_cluster_ttest_age_group_paired$n2 <- NaN
# Tcell_LAIV_cluster_ttest_age_group_paired$mean1 <- NaN
# Tcell_LAIV_cluster_ttest_age_group_paired$mean2 <- NaN

age_pair_list <- data.frame(age_group1 = c('02-03yrs','02-03yrs','07-09yrs'),
                            age_group2 = c('07-09yrs','27-39yrs','27-39yrs'))
# 
# library(dplyr)
# for (cluster_name in annotate_cluster_list) {
#   print(cluster_name)
#   graphics.off()
#   temp <- Tcell_cluster_condition_donor_count_cleaned[Tcell_cluster_condition_donor_count_cleaned$clusters == cluster_name,]
#   rownames(temp) <- temp$condition
#   # if ((!all(temp$cluster_count == 0)) & (!all(temp$cluster_count == 100))) {
#   for (day_group_name in LAIV_day_group_list) {
#     temp1 <- temp[temp$days %in% LAIV_day_group_dict[[day_group_name]],]
#     if (!all(temp1$cluster_count == 0)) {
#       for (age_group_index in c(1:3)) {
#         temp_age_group1 <- age_pair_list$age_group1[age_group_index]
#         temp_age_group2 <- age_pair_list$age_group2[age_group_index]
#         ttest_result <- t.test(temp1$cluster_count[temp1$age_group == temp_age_group1],
#                                temp1$cluster_count[temp1$age_group == temp_age_group2])
#         row_select <- (Tcell_LAIV_cluster_ttest_age_group_paired$cluster == cluster_name) & 
#           (Tcell_LAIV_cluster_ttest_age_group_paired$age_group1 == temp_age_group1) &
#           (Tcell_LAIV_cluster_ttest_age_group_paired$age_group2 == temp_age_group2) &
#           (Tcell_LAIV_cluster_ttest_age_group_paired$day_group == day_group_name)
#         Tcell_LAIV_cluster_ttest_age_group_paired$n1[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == temp_age_group1]))
#         Tcell_LAIV_cluster_ttest_age_group_paired$n2[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == temp_age_group2]))
#         Tcell_LAIV_cluster_ttest_age_group_paired$mean1[row_select] <- ttest_result[["estimate"]][["mean of x"]]
#         Tcell_LAIV_cluster_ttest_age_group_paired$mean2[row_select] <- ttest_result[["estimate"]][["mean of y"]]
#         Tcell_LAIV_cluster_ttest_age_group_paired$pval[row_select] <- ttest_result$p.value
#       }
#     }  
#     
#     # }
#     
#   }
# }
# 
# Tcell_LAIV_cluster_ttest_age_group_paired <- Tcell_LAIV_cluster_ttest_age_group_paired %>% arrange(pval)
# Tcell_LAIV_cluster_ttest_age_group_paired$diff <- Tcell_LAIV_cluster_ttest_age_group_paired$mean1 - Tcell_LAIV_cluster_ttest_age_group_paired$mean2
# 
# # Tcell_LAIV_cluster_ttest_age_group_paired$pval_adj <- p.adjust(Tcell_LAIV_cluster_ttest_age_group_paired$pval,method = 'fdr')
# Tcell_LAIV_cluster_ttest_age_group_paired <- Tcell_LAIV_cluster_ttest_age_group_paired %>% arrange(cluster,pval)
# 
# Tcell_LAIV_cluster_ttest_age_group_paired <- Tcell_LAIV_cluster_ttest_age_group_paired[Tcell_LAIV_cluster_ttest_age_group_paired$cluster != 'PB',]
# write.xlsx(Tcell_LAIV_cluster_ttest_age_group_paired,'tonsil_LAIV_120a_s_Tcell_LAIV_age_cluster_labeled_age_group_ttest_count_paired.xlsx')
# 

######## count: all cell FC diff visualization ###########################
Tcell_LAIV_cluster_ttest_age_group <- expand.grid(cluster = annotate_cluster_list, age_group = age_group_list, day_group = LAIV_day_group_list)
Tcell_LAIV_cluster_ttest_age_group$diff <- NaN
Tcell_LAIV_cluster_ttest_age_group$pval <- NaN
Tcell_LAIV_cluster_ttest_age_group$n1 <- NaN
Tcell_LAIV_cluster_ttest_age_group$n2 <- NaN
Tcell_LAIV_cluster_ttest_age_group$mean1 <- NaN
Tcell_LAIV_cluster_ttest_age_group$mean2 <- NaN
Tcell_LAIV_cluster_correlate_age <- expand.grid(cluster = annotate_cluster_list, day_group = LAIV_day_group_list)
Tcell_LAIV_cluster_correlate_age$r <- NaN
Tcell_LAIV_cluster_correlate_age$pval <- NaN
Tcell_LAIV_cluster_correlate_age$n <- NaN

library(dplyr)
for (cluster_name in annotate_cluster_list) {
  print(cluster_name)
  graphics.off()
  temp <- Tcell_cluster_condition_donor_count_cleaned[Tcell_cluster_condition_donor_count_cleaned$clusters == cluster_name,]
  rownames(temp) <- temp$condition
  # if ((!all(temp$cluster_count == 0)) & (!all(temp$cluster_count == 100))) {
  for (day_group_name in LAIV_day_group_list) {
    temp1 <- temp[temp$days %in% LAIV_day_group_dict[[day_group_name]],]
    if (!all(temp1$cluster_count == 0)) {
      for (age_group_index in c(1:3)) {
        age_group_name <- age_group_list[age_group_index]
        ttest_result <- t.test(temp1$cluster_count[temp1$age_group == age_group_name],
                               temp1$cluster_count[temp1$age_group != age_group_name])
        row_select <- (Tcell_LAIV_cluster_ttest_age_group$cluster == cluster_name) & 
          (Tcell_LAIV_cluster_ttest_age_group$age_group == age_group_list[age_group_index]) &
          (Tcell_LAIV_cluster_ttest_age_group$day_group == day_group_name)
        Tcell_LAIV_cluster_ttest_age_group$n1[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == age_group_name]))
        Tcell_LAIV_cluster_ttest_age_group$n2[row_select] <- length(unique(temp1$donor_ID[temp1$age_group != age_group_name]))
        Tcell_LAIV_cluster_ttest_age_group$mean1[row_select] <- ttest_result[["estimate"]][["mean of x"]]
        Tcell_LAIV_cluster_ttest_age_group$mean2[row_select] <- ttest_result[["estimate"]][["mean of y"]]
        Tcell_LAIV_cluster_ttest_age_group$pval[row_select] <- ttest_result$p.value
      }
      row_select <- (Tcell_LAIV_cluster_correlate_age$cluster == cluster_name) &
        (Tcell_LAIV_cluster_correlate_age$day_group == day_group_name)
      cor_test <- cor.test(temp1$cluster_count, temp1$age)
      Tcell_LAIV_cluster_correlate_age$n[row_select] <- length(unique(temp$donor_ID))
      Tcell_LAIV_cluster_correlate_age$r[row_select] <- cor_test[["estimate"]][["cor"]]
      Tcell_LAIV_cluster_correlate_age$pval[row_select] <- cor_test$p.value  
    }  
    
    # }
    
  }
}

Tcell_LAIV_cluster_ttest_age_group <- Tcell_LAIV_cluster_ttest_age_group %>% arrange(pval)
Tcell_LAIV_cluster_correlate_age <- Tcell_LAIV_cluster_correlate_age %>% arrange(pval)
Tcell_LAIV_cluster_ttest_age_group$diff <- Tcell_LAIV_cluster_ttest_age_group$mean1 - Tcell_LAIV_cluster_ttest_age_group$mean2

# Tcell_LAIV_cluster_ttest_age_group$pval_adj <- p.adjust(Tcell_LAIV_cluster_ttest_age_group$pval,method = 'fdr')
# Tcell_LAIV_cluster_correlate_age$pval_adj <- p.adjust(Tcell_LAIV_cluster_correlate_age$pval,method = 'fdr')
Tcell_LAIV_cluster_ttest_age_group <- Tcell_LAIV_cluster_ttest_age_group %>% arrange(cluster,pval)
Tcell_LAIV_cluster_correlate_age <- Tcell_LAIV_cluster_correlate_age %>% arrange(cluster,desc(r))

Tcell_LAIV_cluster_ttest_age_group <- Tcell_LAIV_cluster_ttest_age_group[Tcell_LAIV_cluster_ttest_age_group$cluster != 'PB',]
Tcell_LAIV_cluster_correlate_age <- Tcell_LAIV_cluster_correlate_age[Tcell_LAIV_cluster_correlate_age$cluster != 'PB',]

Tcell_LAIV_cluster_ttest_age_group <- Tcell_LAIV_cluster_ttest_age_group[Tcell_LAIV_cluster_ttest_age_group$age_group != '07-09yrs',]
write.xlsx(Tcell_LAIV_cluster_ttest_age_group,'tonsil_LAIV_120a_s_Tcell_LAIV_age_cluster_labeled_age_group_ttest_count.xlsx')
write.xlsx(Tcell_LAIV_cluster_correlate_age,'tonsil_LAIV_120a_s_Tcell_LAIV_age_cluster_labeled_age_corr_count.xlsx')
####### count: visualization ####################################
cluster_name <- 'memory B'
day_group_list <- c(0,4,10)
day_group_name <- 'day0'
age_group_name <- '27-39yrs'
marker_name <- 'CXCR5'
x_lab <- 2.2
plot_corr <- 0
graphics.off()
temp_donorID_batch_list <- unlist(batch_donorID_list[which(sapply(batch_marker_list, function(x) (marker_name %in% x)))])
for (day_group_name in setdiff(LAIV_day_group_list,c('day07-08','day07'))) {
  Tcell_cluster_condition_donor_count_cleaned$day_group[Tcell_cluster_condition_donor_count_cleaned$days %in% LAIV_day_group_dict[[day_group_name]]] <-
    day_group_name
}
for (cluster_name in cluster_name){#annotate_cluster_list){
  print(cluster_name)
  temp <- Tcell_cluster_condition_donor_count_cleaned[(Tcell_cluster_condition_donor_count_cleaned$clusters == cluster_name) &
                                                # (Tcell_cluster_condition_donor_count_cleaned$days %in% LAIV_day_group_dict[[day_group_name]]) &
                                                (Tcell_cluster_condition_donor_count_cleaned$day_group %in% setdiff(LAIV_day_group_list,c('day07-08','day07'))) &
                                                Tcell_cluster_condition_donor_count_cleaned$donor_ID %in% temp_donorID_batch_list,]
  
  plot <- ggplot(temp, aes(y = cluster_count, x=age_group)) +
    # geom_boxplot(outlier.shape = NA) +
    # geom_jitter() +
    facet_wrap( ~ day_group, scales = "fixed",ncol = 5) +
    geom_bar(position = 'dodge', stat = 'summary', fun = mean,fill = 'white',color = 'black',width = 0.5) +
    geom_jitter(width = 0.2,size = 1) +
    # geom_point(position = position_jitterdodge(jitter.width = 0.5, jitter.height=0.4, dodge.width=0.9)) +
    geom_hline(yintercept = 0,linetype=2) + 
    ylab('# B cells') +
    ggtitle(cluster_name) +
    theme_bw() +
    ylim(0,max(temp$cluster_count) + 0.5*abs(max(temp$cluster_count))) +
    RotatedAxis()
  # annotate("text", x = x_lab,
  #          y=max(temp$cluster_count) + 0.5*abs(max(temp$cluster_count)),
  #          label = paste("p =",sprintf(ttest_result$p.value,fmt = '%#.3f')),size = 5)
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_cluster_labeled_age_group_count_',cluster_name,'.pdf',sep = ''),width = 6, height = 3)
}

ggplot(temp, aes(y=cluster_count,x=age)) +
  facet_wrap( ~ day_group, scales = "fixed",ncol = 5) +
  geom_point() + RotatedAxis() +
  ylab('% B cells') +
  geom_smooth(method='lm', formula= y~x, color = 'red') +
  ggtitle(cluster_name) +
  theme_bw()#  +
# annotate("text", x=(temp$age[which.max(temp$cluster_count)] + 25)%%40, 
#          y=mean(sort(temp$cluster_count)[(length(temp$cluster_count) - 1):length(temp$cluster_count)]),
#          label = paste("r =",sprintf(cor_test$estimate, fmt = '%#.2f')),size = 6)
dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_cluster_labeled_age_corr_',cluster_name,'.pdf',sep = ''),width = 6, height = 3)

# ttest_result <- t.test(temp$cluster_count[temp$age_group == age_group_name],
#                        temp$cluster_count[temp$age_group != age_group_name])

cor_test <- cor.test(temp$cluster_count,temp$age)
ggplot(temp, aes(y=cluster_count,x=age)) +
  geom_point() + RotatedAxis() +
  ylab('diff in % B cells\nbetween LAIV+ and LAIV-') +
  geom_smooth(method='lm', formula= y~x, color = 'red') +
  ggtitle(paste(cluster_name, day_group_name)) +
  theme_bw() +
  annotate("text", x=(temp$age[which.max(temp$cluster_count)] + 25)%%40, 
           y=mean(sort(temp$cluster_count)[(length(temp$cluster_count) - 1):length(temp$cluster_count)]),
           label = paste("r =",sprintf(cor_test$estimate, fmt = '%#.2f')),size = 6)

dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_cluster_labeled_age_corr_',cluster_name,'_',day_group_name,'.pdf',sep = ''),width = 3, height = 3)


##### within PB cells #####################
Tcell.database <- data.frame(cbind(Tcell_LAIV.subset$donor_ID,as.character(Tcell_LAIV.subset$cluster_labeled),Tcell_LAIV.subset$condition,Tcell_LAIV.subset$days,Tcell_LAIV.subset$stimulation,Tcell_LAIV.subset$age_group,Tcell_LAIV.subset$day_group,Tcell_LAIV.subset$age))
colnames(Tcell.database) <- c('donor_ID','clusters','condition','days','stimulation','age_group','day_group','age')
# sapply(Tcell.database,class)
Tcell_cluster_condition_donor_count <- dplyr::count(Tcell.database, donor_ID, clusters,condition, days,stimulation,age_group,day_group,age)
Tcell_cluster_condition_donor_count <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$clusters != 'B cells',]

Tcell_cluster_condition_donor_count <- Tcell_cluster_condition_donor_count %>% group_by(donor_ID, condition, days,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
Tcell_cluster_condition_donor_count$days <- as.numeric(substring(Tcell_cluster_condition_donor_count$days,first = 4,last = 5))
Tcell_cluster_condition_donor_count$percentage <- Tcell_cluster_condition_donor_count$n/Tcell_cluster_condition_donor_count$total*100
Tcell_cluster_condition_donor_count$age <- as.numeric(Tcell_cluster_condition_donor_count$age)

for (donor_name in unique(Tcell_LAIV.subset$donor_ID)){
  temp_condition_list <- sort(unique(Tcell_LAIV.subset$condition[Tcell_LAIV.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- Tcell_cluster_condition_donor_count[(Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Tcell_cluster_condition_donor_count$clusters)){
      if_row <- ((Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$clusters == cluster_name) & (Tcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount$clusters <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        Tcell_cluster_condition_donor_count[nrow(Tcell_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(Tcell.database)

annotate_cluster_list <- unique(Tcell_LAIV.subset$cluster_labeled)
# temp <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$days != 8,]
# for (day_group_name in setdiff(LAIV_day_group_list,'day07-08')) {
#   temp$day_group[temp$days %in% LAIV_day_group_dict[[day_group_name]]] <- day_group_name
# }
# # temp <- temp[temp$donor_ID != '02yrs M IMD085',]
# ggplot(temp,aes(x=donor_ID, y=percentage,fill = clusters,color = 'black')) +  
#   facet_wrap( ~ day_group, scales = "free_x",ncol = 3) +
#   geom_bar(position="stack", stat="identity") + 
#   theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
#   ylab('Tcell %') +
#   # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
#   theme_bw() + RotatedAxis() 

temp <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$donor_ID != '02yrs M IMD085',]
temp <- temp[temp$clusters != 'IGHG3-secreted',]
ggplot(temp,aes(x=days, y=percentage,color = clusters)) +
  facet_wrap( ~ donor_ID, scales = "free_x",ncol = 3) +
  geom_point(size = 1)+ geom_line(aes(group = interaction(donor_ID,clusters))) + 
  theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  ylab('PB %') +
  # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  theme_bw()

dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_PB_LAIV_age_group_cluster_labeled_donorID_ratio.pdf',sep = ''),width = 6, height = 5)

##### within PB cells #####################
Tcell_ns.subset <- subset(Tcell.subset,stimulation != 'LAIV+')
  
Tcell.database <- data.frame(cbind(Tcell_ns.subset$donor_ID,as.character(Tcell_ns.subset$cluster_labeled),Tcell_ns.subset$condition,Tcell_ns.subset$days,Tcell_ns.subset$stimulation,Tcell_ns.subset$age_group,Tcell_ns.subset$day_group,Tcell_ns.subset$age))
colnames(Tcell.database) <- c('donor_ID','clusters','condition','days','stimulation','age_group','day_group','age')
# sapply(Tcell.database,class)
Tcell_cluster_condition_donor_count <- dplyr::count(Tcell.database, donor_ID, clusters,condition, days,stimulation,age_group,day_group,age)
Tcell_cluster_condition_donor_count <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$clusters != 'B cells',]

Tcell_cluster_condition_donor_count <- Tcell_cluster_condition_donor_count %>% group_by(donor_ID, condition, days,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
Tcell_cluster_condition_donor_count$days <- as.numeric(substring(Tcell_cluster_condition_donor_count$days,first = 4,last = 5))
Tcell_cluster_condition_donor_count$percentage <- Tcell_cluster_condition_donor_count$n/Tcell_cluster_condition_donor_count$total*100
Tcell_cluster_condition_donor_count$age <- as.numeric(Tcell_cluster_condition_donor_count$age)

for (donor_name in unique(Tcell_ns.subset$donor_ID)){
  temp_condition_list <- sort(unique(Tcell_ns.subset$condition[Tcell_ns.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- Tcell_cluster_condition_donor_count[(Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Tcell_cluster_condition_donor_count$clusters)){
      if_row <- ((Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$clusters == cluster_name) & (Tcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount$clusters <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        Tcell_cluster_condition_donor_count[nrow(Tcell_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(Tcell.database)

annotate_cluster_list <- unique(Tcell_ns.subset$cluster_labeled)
# temp <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$days != 8,]
# for (day_group_name in setdiff(LAIV_day_group_list,'day07-08')) {
#   temp$day_group[temp$days %in% LAIV_day_group_dict[[day_group_name]]] <- day_group_name
# }
# # temp <- temp[temp$donor_ID != '02yrs M IMD085',]
# ggplot(temp,aes(x=donor_ID, y=percentage,fill = clusters,color = 'black')) +  
#   facet_wrap( ~ day_group, scales = "free_x",ncol = 3) +
#   geom_bar(position="stack", stat="identity") + 
#   theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
#   ylab('Tcell %') +
#   # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
#   theme_bw() + RotatedAxis() 

temp <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$donor_ID != '02yrs M IMD085',]
temp <- temp[temp$clusters != 'IGHG3-secreted',]
ggplot(temp,aes(x=days, y=percentage,color = clusters)) +
  facet_wrap( ~ donor_ID, scales = "free_x",ncol = 3) +
  geom_point(size = 1)+ geom_line(aes(group = interaction(donor_ID,clusters))) + 
  theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  ylab('PB %') +
  # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  theme_bw()

dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_PB_ns_age_group_cluster_labeled_donorID_ratio.pdf',sep = ''),width = 6, height = 5)

########## cell number #####################
flow_count_table <- read.xlsx('tonsil_LAIV_120a_s_flow_cell_count.xlsx', sheetName = "Sheet1")
flow_count_table <- flow_count_table[,!colnames(flow_count_table) == '...1']
flow_count_table$condition <- gsub('day0 LAIV-','day0',flow_count_table$condition)
flow_count_table$B_cells[grepl('IMD170',flow_count_table$donor_ID)] <- 0

Tcell.database <- data.frame(cbind(Tcell_LAIV.subset$donor_ID,as.character(Tcell_LAIV.subset$cluster_labeled),Tcell_LAIV.subset$condition,Tcell_LAIV.subset$days,Tcell_LAIV.subset$stimulation,Tcell_LAIV.subset$age_group,Tcell_LAIV.subset$day_group,Tcell_LAIV.subset$age))
colnames(Tcell.database) <- c('donor_ID','clusters','condition','days','stimulation','age_group','day_group','age')
# sapply(Tcell.database,class)
Tcell.database <- Tcell.database[Tcell.database$clusters == 'memory B',]
Tcell_cluster_condition_donor_count <- dplyr::count(Tcell.database, donor_ID, clusters,condition, days,stimulation,age_group,day_group,age)
Tcell_cluster_condition_donor_count <- Tcell_cluster_condition_donor_count %>% group_by(donor_ID, condition, days,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
Tcell_cluster_condition_donor_count$days <- as.numeric(substring(Tcell_cluster_condition_donor_count$days,first = 4,last = 5))
Tcell_cluster_condition_donor_count$percentage <- Tcell_cluster_condition_donor_count$n/Tcell_cluster_condition_donor_count$total*100
Tcell_cluster_condition_donor_count$age <- as.numeric(Tcell_cluster_condition_donor_count$age)
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(Tcell_LAIV.subset$donor_ID)){
  temp_condition_list <- sort(unique(Tcell_LAIV.subset$condition[Tcell_LAIV.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- Tcell_cluster_condition_donor_count[(Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Tcell_LAIV.subset$cluster_labeled)){
      if_row <- ((Tcell_cluster_condition_donor_count$donor_ID == donor_name) & (Tcell_cluster_condition_donor_count$clusters == cluster_name) & (Tcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount$clusters <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        Tcell_cluster_condition_donor_count[nrow(Tcell_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(Tcell.database)

Tcell_cluster_condition_donor_count_cleaned <- Tcell_cluster_condition_donor_count 
p <- match(Tcell_cluster_condition_donor_count$condition, flow_count_table$condition)
temp <- flow_count_table[p,]
Tcell_cluster_condition_donor_count$Tcell_count <- temp$B_cells
Tcell_cluster_condition_donor_count$cluster_count <- Tcell_cluster_condition_donor_count$percentage*Tcell_cluster_condition_donor_count$Tcell_count/100
temp <- Tcell_cluster_condition_donor_count[Tcell_cluster_condition_donor_count$donor_ID != '02yrs M IMD085',]
temp <- temp[temp$clusters != 'IGHG3-secreted',]
temp <- temp[temp$clusters != 'B cells',]
ggplot(temp,aes(x=days, y=cluster_count,color = clusters)) +
  facet_wrap( ~ donor_ID, scales = "free_x",ncol = 3) +
  geom_point(size = 1)+ geom_line(aes(group = interaction(donor_ID,clusters))) + 
  theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  ylab('#Tcells') +
  # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  theme_bw()

dev.print(pdf, paste('tonsil_LAIV_120a_s_Tcell_LAIV_age_group_cluster_labeled_donorID_count.pdf',sep = ''),width = 6, height = 5)
