library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(ggrepel)
library('stringr')
# setwd("/Volumes/GoogleDrive/My\ Drive/Stanford/RNA-seq/data_analysis/tonsil_LAIV_Rhapsody_120a-s/")

Bcell_LAIV.subset <- readRDS('tonsil_LAIV_120a_s_Bcell_lognorm_cca_LAIV.rds')
Bcell_LAIV.subset <- subset(Bcell.subset,stimulation != 'LAIV-')
Bcell_LAIV.subset <- subset(Bcell_LAIV.subset,days != 'day02')
Bcell_LAIV.subset <- subset(Bcell_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
         
# unique(Bcell_LAIV.subset$condition[Bcell_LAIV.subset$days == 'day12'])               
saveRDS(Bcell_LAIV.subset,'tonsil_LAIV_120a_s_Bcell_lognorm_cca_LAIV.rds')

LAIV_day_group_dict <- list()
LAIV_day_group_dict[['day0']] <- c(0)
LAIV_day_group_dict[['day04']] <- c(4)
LAIV_day_group_dict[['day10']] <- c(10)
LAIV_day_group_dict[['day06-07']] <- c(6,7)
LAIV_day_group_dict[['day12-14']] <- c(12,14)
LAIV_day_group_list <- names(LAIV_day_group_dict)

age_group_list <- sort(unique(Bcell_LAIV.subset$age_group))

batch_donorID_list <- list()
batch_donorID_list[[1]] <- unique(Bcell_LAIV.subset$donor_ID[Bcell_LAIV.subset$batch == 'batch1'])
batch_donorID_list[[2]] <- unique(Bcell_LAIV.subset$donor_ID[Bcell_LAIV.subset$batch == 'batch2'])
batch_donorID_list[[3]] <- unique(Bcell_LAIV.subset$donor_ID[Bcell_LAIV.subset$batch == 'batch3'])
batch_donorID_list[[4]] <- unique(Bcell_LAIV.subset$donor_ID[Bcell_LAIV.subset$batch == 'batch4'])

batch_marker_list <- list()
batch_marker_list[[1]] <- readRDS('tonsil_LAIV_120a.rds')
batch_marker_list[[2]] <- readRDS('tonsil_LAIV_120b.rds')
batch_marker_list[[3]] <- readRDS('tonsil_LAIV_120c.rds')
batch_marker_list[[4]] <- readRDS('tonsil_LAIV_120j.rds')
########### Bcell gating ###########################################
Bcell_LAIV.subset$BCL6_pos <- 'Bcells'
Bcell_LAIV.subset$BCL6_pos[(as.matrix(Bcell_LAIV.subset@assays[["log1p"]]@data['BCL6',]) > 0)] <- 'BCL6+B'
Bcell_LAIV.subset$BCL6_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'Bcells'

Bcell_LAIV.subset$AICDA_pos <- 'AICDA-B'
Bcell_LAIV.subset$AICDA_pos[(as.matrix(Bcell_LAIV.subset@assays[["log1p"]]['AICDA']) > 0)] <- 'AICDA+B'
Bcell_LAIV.subset$AICDA_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Bcell_LAIV.subset$MZB1_pos <- 'MZB1-B'
Bcell_LAIV.subset$MZB1_pos[(as.matrix(Bcell_LAIV.subset@assays[["log1p"]]['MZB1']) > 0)] <- 'MZB1+B'
Bcell_LAIV.subset$MZB1_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Bcell_LAIV.subset$AbCXCR5_pos <- 'AbCXCR5-B'
Bcell_LAIV.subset$AbCXCR5_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0)] <- 'AbCXCR5+B'
Bcell_LAIV.subset$AbCXCR5_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'PB'
Bcell_LAIV.subset$AbCXCR5_pos[Bcell_LAIV.subset$differentiation == 'memory B'] <- 'PB'

Bcell_LAIV.subset$AbCD83_pos <- 'AbCD83-B'
Bcell_LAIV.subset$AbCD83_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD83']) > 0)] <- 'AbCD83+B'
Bcell_LAIV.subset$AbCD83_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Bcell_LAIV.subset$AbCXCR4_pos <- 'AbCXCR4-B'
Bcell_LAIV.subset$AbCXCR4_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD184']) > 0)] <- 'AbCXCR4+B'
Bcell_LAIV.subset$AbCXCR4_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'PB'

Bcell_LAIV.subset$AbCXCR5AbCD83_pos <- 'Bcell'
Bcell_LAIV.subset$AbCXCR5AbCD83_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD83']) > 0)] <- 'AbCXCR5+AbCD83+B'
Bcell_LAIV.subset$AbCXCR5AbCD83_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'Bcell'

Bcell_LAIV.subset$AbCXCR5AbCXCR4_pos <- 'Bcell'
Bcell_LAIV.subset$AbCXCR5AbCXCR4_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & 
                                  (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD184']) > 0)] <- 'AbCXCR5+AbCXCR4+B'
Bcell_LAIV.subset$AbCXCR5AbCXCR4_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'Bcell'

Bcell_LAIV.subset$AbCD27AbCD38_pos <- 'Bcell'
Bcell_LAIV.subset$AbCD27AbCD38_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD27']) > 0) & 
                                       (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD38']) > 0)] <- 'AbCD27+AbCD38+B'
Bcell_LAIV.subset$AbCD27AbCD38_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'Bcell'

Bcell_LAIV.subset$AbCD11cTBET_pos <- 'B cells'
Bcell_LAIV.subset$AbCD11cTBET_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data)['ab-CD11c',] > 0) & 
                                 (as.matrix(Bcell_LAIV.subset@assays[["log1p"]]@data)['TBX21',] > 0)] <- 'AbCD11c+TBX21+B'
Bcell_LAIV.subset$AbCD11cTBET_pos[Bcell_LAIV.subset$differentiation %in% c('PB')] <- 'PB'

Bcell_LAIV.subset$TBET_pos <- 'B cells'
Bcell_LAIV.subset$TBET_pos[(as.matrix(Bcell_LAIV.subset@assays[["log1p"]]@data)['TBX21',] > 0)] <- 'TBX21+B'
Bcell_LAIV.subset$TBET_pos[Bcell_LAIV.subset$differentiation %in% c('PB')] <- 'PB'

Bcell_LAIV.subset$AbCD24AbCD38_pos <- 'Bcell'
Bcell_LAIV.subset$AbCD24AbCD38_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD38']) > 0) & 
                                (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD24']) > 0)] <- 'AbCD24AbCD38+B'
Bcell_LAIV.subset$AbCD24AbCD38_pos[Bcell_LAIV.subset$differentiation == 'PB'] <- 'Bcell'


Bcell_LAIV.subset$BCL6_differentiation <- 'BCL6-B'
Bcell_LAIV.subset$BCL6_differentiation[Bcell_LAIV.subset$BCL6_pos == 'BCL6+B'] <- 'BCL6+B'
Bcell_LAIV.subset$BCL6_differentiation[(Bcell_LAIV.subset$BCL6_pos == 'BCL6-B') & (Bcell_LAIV.subset$differentiation != 'memory B')] <- 'BCL6-B'
Bcell_LAIV.subset$BCL6_differentiation[(Bcell_LAIV.subset$BCL6_pos == 'BCL6-B') & (Bcell_LAIV.subset$differentiation == 'memory B')] <- 'memory B'
Bcell_LAIV.subset$BCL6_differentiation[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'PB'

Bcell_LAIV.subset$AbCXCR5_differentiation <- 'AbCXCR5-B'
Bcell_LAIV.subset$AbCXCR5_differentiation[Bcell_LAIV.subset$AbCXCR5_pos == 'AbCXCR5+B'] <- 'AbCXCR5+B'
Bcell_LAIV.subset$AbCXCR5_differentiation[(Bcell_LAIV.subset$AbCXCR5_pos == 'AbCXCR5-B') & (Bcell_LAIV.subset$differentiation != 'memory B')] <- 'AbCXCR5-B'
Bcell_LAIV.subset$AbCXCR5_differentiation[(Bcell_LAIV.subset$AbCXCR5_pos == 'AbCXCR5-B') & (Bcell_LAIV.subset$differentiation == 'memory B')] <- 'memory B'
Bcell_LAIV.subset$AbCXCR5_differentiation[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'PB'


Bcell_LAIV.subset$AbCD11c_AbCXCR5_differentiation <- 'Bcells'
Bcell_LAIV.subset$AbCD11c_AbCXCR5_differentiation[(Bcell_LAIV.subset$AbCXCR5_differentiation == 'AbCXCR5-B') & 
                                                    (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD11c']) > 0) & 
                                                    (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD24']) > 0) & 
                                                    (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD38']) > 0)] <- 
  'aNAV'
Bcell_LAIV.subset$AbCD11c_AbCXCR5_differentiation[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'PB'

Bcell_LAIV.subset$DN2 <- 'Bcells'
Bcell_LAIV.subset$DN2[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data['ab-CD21',]) < 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data['ab-CD27',]) < 0)] <- 'DN2'

Bcell_LAIV.subset$DN2[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'Bcells'

Bcell_LAIV.subset$DN2 <- 'Bcells'
Bcell_LAIV.subset$DN2[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) < 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["log1p"]]['BCL6']) <= 0) &
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD21']) < 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-IgD']) < 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD27']) < 0)] <- 
  'DN2'
Bcell_LAIV.subset$DN2[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CXCR5']) > 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD21']) > 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-IgD']) < 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD27']) < 0)] <- 
  'DN1'
Bcell_LAIV.subset$DN2[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'PB'

Bcell_LAIV.subset$MKI67_BCL6_differentiation <- 
  paste(ifelse(as.matrix(Bcell_LAIV.subset@assays[["log1p"]]['MKI67']) > 0,'MKI67+','MKI67-'),Bcell_LAIV.subset$BCL6_differentiation, sep = '')


Bcell_LAIV.subset$MZB_abgate <- 'B cells'
Bcell_LAIV.subset$MZB_abgate[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data)['ab-IgD',] > 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data)['ab-IgM',] > 0) & 
                        (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data)['ab-CD27',] > 0)] <- 
  'MZB'
Bcell_LAIV.subset$MZB_abgate[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'B cells'

Bcell_LAIV.subset$MZB_abgate2 <- 'B cells'
Bcell_LAIV.subset$MZB_abgate2[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data)['ab-CD20',] > 0) & 
                               (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data)['ab-IgM',] > 0) & 
                               (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data)['ab-CD21',] > 0)] <- 
  'MZB'
Bcell_LAIV.subset$MZB_abgate2[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'B cells'

Bcell_LAIV.subset$IGHM_sum_lognorm <- Bcell_LAIV.subset@assays[["normalized"]]@data['IGHM-secreted'] + Bcell_LAIV.subset@assays[["normalized"]]@data['IGHM-membrane']
Bcell_LAIV.subset$IGHG1_sum_lognorm <- Bcell_LAIV.subset@assays[["normalized"]]@data['IGHG1-secreted'] + Bcell_LAIV.subset@assays[["normalized"]]@data['IGHG1-membrane']
Bcell_LAIV.subset$MZB_genegate <- 'B cells'
Bcell_LAIV.subset$MZB_genegate[(as.matrix(Bcell_LAIV.subset@assays[["normalized"]]@data)['IGHD-membrane',] > 0) & 
                                 (as.matrix(Bcell_LAIV.subset@assays[["normalized"]]@data)['IGHM-membrane',] > 0) & 
                               (as.matrix(Bcell_LAIV.subset@assays[["normalized"]]@data)['CD27',] > 0)] <- 
  'MZB'
Bcell_LAIV.subset$MZB_genegate[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'B cells'

Bcell_LAIV.subset$CD1C_pos <- 'B cells'
Bcell_LAIV.subset$CD1C_pos[(as.matrix(Bcell_LAIV.subset@assays[["normalized"]]@data)['CD1C',] > 0)] <- 
  'CD1C+B'
Bcell_LAIV.subset$CD1C_pos[(Bcell_LAIV.subset$differentiation == 'PB')] <- 'B cells'


Bcell_LAIV.subset$AbCD11cAbCD21_pos <- 'Bcell'
Bcell_LAIV.subset$AbCD11cAbCD21_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD11c']) > 0) & 
                                    (as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD21']) < 0)] <- 'AbCD11c+AbCD21-B'
Bcell_LAIV.subset$AbCD11cAbCD21_pos[Bcell_LAIV.subset$differentiation %in% c('PB')] <- 'B cells'

PB_secreted <- data.frame(t(as.matrix(Bcell_LAIV.subset@assays$raw@data)),check.names = F)[,rownames(Bcell_LAIV.subset)[grepl('-secreted',rownames(Bcell_LAIV.subset))]]
PB_secreted$max <- colnames(PB_secreted)[apply(PB_secreted,1,which.max)]
PB_secreted$isotype <- PB_secreted$max
PB_secreted$isotype[grepl('IGHG',PB_secreted$max)] <- 'IgG+PB'
PB_secreted$isotype[grepl('IGHM',PB_secreted$max)] <- 'IgM+PB'
PB_secreted$isotype[grepl('IGHA',PB_secreted$max)] <- 'IgA+PB'
PB_secreted$isotype[grepl('IGHE',PB_secreted$max)] <- 'IgE+PB'
Bcell_LAIV.subset$isotype <- PB_secreted$isotype
Bcell_LAIV.subset$isotype_secreted <- PB_secreted$max
Bcell_LAIV.subset$isotype[Bcell_LAIV.subset$differentiation != 'PB'] <- 'B cells'
Bcell_LAIV.subset$isotype_secreted[Bcell_LAIV.subset$differentiation != 'PB'] <- 'B cells'

Bcell_LAIV.subset$IgM_sum <- as.matrix(Bcell_LAIV.subset@assays[["normalized"]]@data)['IGHM-membrane',] + as.matrix(Bcell_LAIV.subset@assays[["normalized"]]@data)['IGHM-secreted',]


PB_secreted <- data.frame(t(as.matrix(Bcell_LAIV.subset@assays$raw@data)),check.names = F)[,rownames(Bcell_LAIV.subset)[grepl('-secreted',rownames(Bcell_LAIV.subset))]]
PB_secreted$max <- colnames(PB_secreted)[apply(PB_secreted,1,which.max)]
PB_secreted$isotype <- PB_secreted$max
PB_secreted$isotype[grepl('IGHG',PB_secreted$max)] <- 'IgG+PB'
PB_secreted$isotype[grepl('IGHM',PB_secreted$max)] <- 'IgM+PB'
PB_secreted$isotype[grepl('IGHA',PB_secreted$max)] <- 'IgA+PB'
PB_secreted$isotype[grepl('IGHE',PB_secreted$max)] <- 'IgE+PB'
Bcell_LAIV.subset$isotype <- PB_secreted$isotype
Bcell_LAIV.subset$isotype_secreted <- PB_secreted$max
Bcell_LAIV.subset$isotype[Bcell_LAIV.subset$differentiation != 'PB'] <- 'B cells'
Bcell_LAIV.subset$isotype_secreted[Bcell_LAIV.subset$differentiation != 'PB'] <- 'B cells'


PB_secreted <- data.frame(t(as.matrix(Bcell_LAIV.subset@assays$raw@data)),check.names = F)[,rownames(Bcell_LAIV.subset)[grepl('-secreted',rownames(Bcell_LAIV.subset)) | grepl('-membrane',rownames(Bcell_LAIV.subset))]]
PB_secreted$max <- colnames(PB_secreted)[apply(PB_secreted,1,which.max)]
Bcell_LAIV.subset$isotype_all <- PB_secreted$max

Bcell_LAIV.subset$isotype_both <- Bcell_LAIV.subset$isotype_all
Bcell_LAIV.subset$isotype_both <- gsub('-.*','+',Bcell_LAIV.subset$isotype_both)
Bcell_LAIV.subset$isotype_both[Bcell_LAIV.subset$differentiation == 'PB'] <- paste(Bcell_LAIV.subset$isotype_both[Bcell_LAIV.subset$differentiation == 'PB'],'PB',sep = '')
Bcell_LAIV.subset$isotype_both[Bcell_LAIV.subset$differentiation != 'PB'] <- paste(Bcell_LAIV.subset$isotype_both[Bcell_LAIV.subset$differentiation != 'PB'],'B',sep = '')

# Bcell_LAIV.subset$isotype_main_both <- Bcell_LAIV.subset$isotype_all
# Bcell_LAIV.subset$isotype_main_both <- gsub('-.*','+',Bcell_LAIV.subset$isotype_both)
# Bcell_LAIV.subset$isotype_main_both[Bcell_LAIV.subset$differentiation == 'PB'] <- paste(Bcell_LAIV.subset$isotype_both[Bcell_LAIV.subset$differentiation == 'PB'],'PB',sep = '')
# Bcell_LAIV.subset$isotype_main_both[Bcell_LAIV.subset$differentiation != 'PB'] <- paste(Bcell_LAIV.subset$isotype_both[Bcell_LAIV.subset$differentiation != 'PB'],'B',sep = '')

Bcell_LAIV.subset$isotype_cluster <- Bcell_LAIV.subset$isotype_all
Bcell_LAIV.subset$isotype_cluster <- gsub('-.*','max_',Bcell_LAIV.subset$isotype_cluster)
Bcell_LAIV.subset$isotype_cluster <- paste(Bcell_LAIV.subset$isotype_cluster,Bcell_LAIV.subset$CD27CD38_cluster,sep = '')

Bcell_LAIV.subset$AbCXCR3_pos <- 'Bcell'
Bcell_LAIV.subset$AbCXCR3_pos[(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]['ab-CD183']) > 0) &
                                      Bcell_LAIV.subset$BCL6_differentiation == 'BCL6-B'] <- 'AbCXCR3+BCL6-B'
Bcell_LAIV.subset$AbCXCR3_pos[Bcell_LAIV.subset$differentiation %in% c('PB')] <- 'Bcell'

Bcell_LAIV.subset$IGHG2_pos <- 'Bcell'
Bcell_LAIV.subset$IGHG2_pos[(as.matrix(Bcell_LAIV.subset@assays[["log1p"]]['IGHG2-secreted']) > 5) & 
                              (!(Bcell_LAIV.subset$BCL6_differentiation %in% c('PB','memory B')))] <- 'IGHG2+B'
Bcell_LAIV.subset$IGHG2_pos[(as.matrix(Bcell_LAIV.subset@assays[["log1p"]]['IGHG2-secreted']) > 5) & 
                              (Bcell_LAIV.subset$BCL6_differentiation %in% c('memory B'))] <- 'IGHG2+memory B'


DimPlot(Bcell_LAIV.subset, label = T, reduction = "umap",group.by = 'isotype', order = F,pt.size = 0.1)
dev.print(png, paste('tonsil_LAIV_120a_s_Bcell_umap_AbCD11cTBET_pos.png',sep = ''),width = 600, height = 400)

Bcell_log1p_dataframe <- data.frame(t(as.matrix(Bcell_LAIV.subset@assays[["log1p"]]@data)),check.names = F)
Bcell_integrated_scale_dataframe <- data.frame(t(as.matrix(Bcell_LAIV.subset@assays[["integrated_scale"]]@data)),check.names = F)

Bcell_integrated_scale_dataframe$AbBcell5RA_AbBcell5RO_pos <- Bcell_LAIV.subset$AbBcell5RA_AbBcell5RO_pos
ggplot(Bcell_integrated_scale_dataframe, aes(x = !!sym(c('ab-CXCR5')), y = !!sym(c('ab-CD279')))) + ggtitle('Bcell cells') +#!!sym(c('TBX21'))
  geom_point(alpha = 0.1,size = 0.5,aes(color = !!sym(c('AbBcell5RA_AbBcell5RO_pos')))) +
  # geom_point(alpha = 0.1,size = 0.5) +
  # scale_colour_manual(values = c("naive Bcell" = "#F8766D", "Bcell" = '#00BFC4'))+
  # scale_colour_manual(values = c("Bcell5RA+CD83+ nonPB B" = "gold", "CD184-CD83- nonPB B" = "blue","CD184+CD83- nonPB B" = 'cyan',"CD184-CD83+ nonPB B" = 'magenta')) +
  geom_density_2d(color = 'black')# +
# ylim(0,8) + xlim(0,8)
# geom_hline(aes(yintercept=3.25,color = 'red')) + guides(color = 'none') +
# geom_vline(aes(xintercept=4.1,color = 'red')) + guides(color = 'none')
dev.print(pdf, 'tonsil_LAIV_120a_s_Bcell_AbCXCR5_AbCD279_contour_integrated_scale.pdf',width = 7, height = 4.6)

DefaultAssay(Bcell_LAIV.subset) <- 'integrated_scale'
FeatureScatter(Bcell_LAIV.subset, feature1 = "ab-Bcell5RA", feature2 = "ab-Bcell5RO") 
dev.print(png, 'tonsil_LAIV_120a_s_Bcell_scatter_abBcell5RA_abBcell5RO.png',width = 700, height = 434)

##### all Bcells: cell cluster fraction visualization #####################
flow_count_table <- read.xlsx('tonsil_LAIV_120a_s_flow_cell_count.xlsx', sheetName = "Sheet1")
flow_count_table <- flow_count_table[,!colnames(flow_count_table) == '...1']
flow_count_table$condition <- gsub('day0 LAIV-','day0',flow_count_table$condition)
flow_count_table$B_cells[grepl('IMD170',flow_count_table$donor_ID)] <- 0

Bcell.database <- data.frame(cbind(Bcell_LAIV.subset$donor_ID,as.character(Bcell_LAIV.subset$BCL6_pos),Bcell_LAIV.subset$condition,Bcell_LAIV.subset$days,Bcell_LAIV.subset$stimulation,Bcell_LAIV.subset$age_group,Bcell_LAIV.subset$day_group,Bcell_LAIV.subset$age))
colnames(Bcell.database) <- c('donor_ID','clusters','condition','days','stimulation','age_group','day_group','age')
# sapply(Bcell.database,class)
Bcell_cluster_condition_donor_count <- dplyr::count(Bcell.database, donor_ID, clusters,condition, days,stimulation,age_group,day_group,age)
Bcell_cluster_condition_donor_count <- Bcell_cluster_condition_donor_count %>% group_by(donor_ID, condition, days,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
Bcell_cluster_condition_donor_count$days <- as.numeric(substring(Bcell_cluster_condition_donor_count$days,first = 4,last = 5))
Bcell_cluster_condition_donor_count$percentage <- Bcell_cluster_condition_donor_count$n/Bcell_cluster_condition_donor_count$total*100
Bcell_cluster_condition_donor_count$age <- as.numeric(Bcell_cluster_condition_donor_count$age)
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(Bcell_LAIV.subset$donor_ID)){
  temp_condition_list <- sort(unique(Bcell_LAIV.subset$condition[Bcell_LAIV.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$donor_ID == donor_name) & (Bcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Bcell_LAIV.subset$BCL6_pos)){
      if_row <- ((Bcell_cluster_condition_donor_count$donor_ID == donor_name) & (Bcell_cluster_condition_donor_count$clusters == cluster_name) & (Bcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount$clusters <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        Bcell_cluster_condition_donor_count[nrow(Bcell_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(Bcell.database)

p <- match(Bcell_cluster_condition_donor_count$condition, flow_count_table$condition)
temp <- flow_count_table[p,]
Bcell_cluster_condition_donor_count$Bcell_count <- temp$B_cells
Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '02yrs M IMD085 day0'] <- 
  Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '02yrs M IMD085 day0']*6/10
Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '03yrs F IMD035 day0'] <- 
  Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '03yrs F IMD035 day0']*6/10
Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '07yrs F IMD150 day0'] <- 
  Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '07yrs F IMD150 day0']*6/12
Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '09yrs M IMD107 day0'] <- 
  Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '09yrs M IMD107 day0']*6/12
Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '27yrs F VIP031 day0'] <- 
  Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '27yrs F VIP031 day0']*6/10
Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '39yrs M VIP015 day0'] <- 
  Bcell_cluster_condition_donor_count$Bcell_count[Bcell_cluster_condition_donor_count$condition == '39yrs M VIP015 day0']*6/12
Bcell_cluster_condition_donor_count$cluster_count <- Bcell_cluster_condition_donor_count$Bcell_count*Bcell_cluster_condition_donor_count$percentage/100
# write.csv(Bcell_cluster_condition_donor_count,'tonsil_LAIV_120a_s_Bcell_log1p_cca_CycleRegressOut_BCL6_pos_LAIV_fraction.csv')
# saveRDS(Bcell_cluster_condition_donor_count,'tonsil_LAIV_120a_s_Bcell_LAIV_log1p_cca_CycleRegressOut_BCL6_pos_fraction.rds')
annotate_cluster_list <- unique(Bcell_LAIV.subset$BCL6_pos)#'AbCD278+B'#
# annotate_cluster_list <- annotate_cluster_list[grepl('MKI67\\+',annotate_cluster_list)]# | grepl('CD278\\+PD1\\-CD4',annotate_cluster_list) | grepl('CD278\\-PD1\\+CD4',annotate_cluster_list)]
annotate_cluster_list <- setdiff(annotate_cluster_list,c('PB','Bcells'))
for (cluster_name in annotate_cluster_list)
{
  graphics.off()
  print(cluster_name)
  temp_Bcell_cluster <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$clusters == cluster_name),]
  # plot <- ggplot(temp_Bcell_cluster,aes(x=days, y=percentage,color = age_group)) +  ggtitle(cluster_name) +
  #   facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
  #   geom_point(size = 1)+ geom_line(aes(group = donor_ID)) + 
  #   scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
  #   theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  #   ylab('Bcell %') +
  #   # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  #   theme_bw()
  # print(plot)
  # dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_lognorm_LAIV_age_group_BCL6_pos_',cluster_name,'_condition_ratio2.pdf',sep = ''),width = 5, height = 4)
  plot <- ggplot(temp_Bcell_cluster,aes(x=days, y=percentage,color = age_group,fill = age_group,shape = age_group, linetype = age_group)) +  ggtitle(cluster_name) +
    geom_point(size = 1.5) + #,aes(shape = donor_ID)
    geom_smooth(method = "loess",se = F, linewidth=0.5) +
    ylim(0,max(temp_Bcell_cluster$percentage)) +
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('% in B cells') +
    scale_shape_manual(values = c('02-03yrs' = 25,'07-09yrs'= 23,'27-39yrs' = 22)) +
    scale_linetype_manual(values = c('02-03yrs' = 'dotted','07-09yrs'= 'longdash','27-39yrs' = 'solid')) +
    scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
    # geom_rect(aes(xmin=4, xmax=10, ymin=0, ymax=Inf),alpha = 0.1,fill = 'gray',color = 'gray') +
    theme_bw()
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_BCL6_pos_',gsub('\\/','or',cluster_name),'_condition_ratio.pdf',sep = ''),width = 2.6, height = 1.5)
  temp_Bcell_cluster <- temp_Bcell_cluster[temp_Bcell_cluster$donor_ID != '07yrs M IMD170',]
  # plot <- ggplot(temp_Bcell_cluster,aes(x=days, y=cluster_count,color = age_group)) +  ggtitle(cluster_name) +
  #   facet_wrap( ~ age_group, scales = "free_x",ncol = 3) +
  #   geom_point(size = 1)+ geom_line(aes(group = donor_ID)) + 
  #   scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
  #   theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  #   ylab('#Bcells') +
  #   # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  #   theme_bw()
  # print(plot)
  # dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_BCL6_pos_',cluster_name,'_condition_count2.pdf',sep = ''),width = 5, height = 4)
  plot <- ggplot(temp_Bcell_cluster,aes(x=days, y=cluster_count,color = age_group,fill = age_group,shape = age_group)) +  ggtitle(cluster_name) +
    geom_point(size = 2) + #,aes(shape = donor_ID)
    geom_smooth(method = "loess",se = F) +
    theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
    ylab('#Bcells') +
    ylim(0,max(temp_Bcell_cluster$cluster_count)) +
    scale_shape_manual(values = c('02-03yrs' = 25,'07-09yrs'= 23,'27-39yrs' = 22)) +
    scale_colour_manual(values = c('02-03yrs' = '#F16C23','07-09yrs'= '#1b7c3d','27-39yrs' = '#2b6a99')) +
    # geom_rect(aes(xmin=4, xmax=10, ymin=0, ymax=Inf),alpha = 0.1,fill = 'gray',color = 'gray') +
    theme_bw()
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_BCL6_pos_',gsub('\\/','or',cluster_name),'_condition_count.pdf',sep = ''),width = 3, height = 2)
}
#### split by age ################ 
# mean of each cell count
library(purrr) 
library(ggalluvial)
temp <- Bcell_cluster_condition_donor_count
temp <- temp[!(temp$donor_ID %in% c('02yrs M IMD085','07yrs M IMD170')),]
temp <- temp[temp$day_group %in% c('day0','day04','day06-07','day10','day12-14'),]
temp <- temp %>% group_by(age_group,clusters,day_group) %>% mutate(avg_cluster_count = mean(cluster_count),avg_fraction = mean(percentage))
temp <- unique(temp[,c('age_group','clusters','day_group','avg_cluster_count','avg_fraction')])
temp$days <- 0
temp$days[temp$day_group == 'day04'] <- 4
temp$days[temp$day_group == 'day06-07'] <- 7
temp$days[temp$day_group == 'day10'] <- 10
temp$days[temp$day_group == 'day12-14'] <- 14

p <- DimPlot(Bcell_LAIV.subset, group.by = 'BCL6_pos', label = TRUE, reduction = "umap")
pbuild <- ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # Pull the data used for the plot
pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
cols_group <- unique(pdata$colour) # Get a vector of unique colors
ordered_BCL6_pos_list <- c('naive','cytokine-rich','CXCR4+act.','CD23+act.','CD83+folli.','prolif. CD83+folli.','preGC/folli.',
                                    'switched GC','GC LZ','GC-like','memory','switched memory','IgA+/IgM+PB','IgG+PB')
# ordered_BCL6_pos_list <- c('naive','act.-like folli.','CXCR4+','act.-like','switching folli./memory','folli. memory',
#                                     'GC-like','switched memory','switched GC','GC','GC LZ','CD83+CXCR4+GC','switched act. memory','unswitched act.','prolif. PB-like GC','IgM+PB','IgA+PB','IgG+PB')

names(cols_group) <- ordered_BCL6_pos_list
# # Idents(TB_shrink.seurat_update) <- TB_shrink.seurat_update@meta.data$sample_tag
temp$clusters <- factor(temp$clusters,levels = ordered_BCL6_pos_list)
temp1 <- temp[temp$clusters %in% c('switched memory','memory','GC-like'),]
temp1 <- temp[temp$clusters %in% c('switched GC','GC LZ','preGC/folli.','prolif. CD83+folli.','CD83+folli.'),]
temp1 <- temp[temp$clusters %in% c('CD23+act.','CXCR4+act.','cytokine-rich','naive'),]

ggplot(temp1, aes(x = days, stratum = age_group, alluvium = clusters, y = avg_cluster_count, label = clusters,fill = clusters)) +#, 
  geom_alluvium(aes(fill = clusters), width = 1,alpha = 1,color = 'black',linewidth = 0.1) +
  facet_grid(age_group ~ .,scales = "fixed") +
  # geom_stratum(width = 1/12, fill = "black",size = 0.5) +
  scale_fill_manual(values = cols_group) +
  ggtitle('') + ylab('#cells') + #ylim(0,300) + xlim(0,14)
theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_BCL6_pos_count.pdf',sep = ''),width = 4, height = 5)

ggplot(temp1, aes(x = days, stratum = age_group, alluvium = clusters, y = avg_fraction, label = clusters,fill = clusters)) +#, 
  geom_alluvium(aes(fill = clusters), width = 0.5,alpha = 1,color = 'black',size = 0.1) +
  facet_grid(age_group ~ .,scales = "fixed") +
  # geom_stratum(width = 1/12, fill = "black") +
  scale_fill_manual(values = cols_group) +
  ggtitle('') + ylab('% in B cells') + #ylim(0,300) + xlim(0,14)
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_BCL6_pos_fraction.pdf',sep = ''),width = 4, height = 5)

######## fraction: compare each age_group pairs ###########################
Bcell_LAIV_cluster_ttest_age_group_paired <- expand.grid(clusters = annotate_cluster_list, age_group1 = age_group_list[1], age_group2 = age_group_list[3],
                                                  day_group = LAIV_day_group_list)
Bcell_LAIV_cluster_ttest_age_group_paired <- transform(Bcell_LAIV_cluster_ttest_age_group_paired, age_group1 = as.character(age_group1),
                                                age_group2 = as.character(age_group2))
Bcell_LAIV_cluster_ttest_age_group_paired <- Bcell_LAIV_cluster_ttest_age_group_paired[Bcell_LAIV_cluster_ttest_age_group_paired$age_group1 != Bcell_LAIV_cluster_ttest_age_group_paired$age_group2,]
Bcell_LAIV_cluster_ttest_age_group_paired$diff <- NaN
Bcell_LAIV_cluster_ttest_age_group_paired$pval <- NaN
Bcell_LAIV_cluster_ttest_age_group_paired$n1 <- NaN
Bcell_LAIV_cluster_ttest_age_group_paired$n2 <- NaN
Bcell_LAIV_cluster_ttest_age_group_paired$mean1 <- NaN
Bcell_LAIV_cluster_ttest_age_group_paired$mean2 <- NaN

age_pair_list <- data.frame(age_group1 = c('02-03yrs','02-03yrs','07-09yrs'),
                            age_group2 = c('07-09yrs','27-39yrs','27-39yrs'))

library(dplyr)
for (cluster_name in annotate_cluster_list) {
  print(cluster_name)
  graphics.off()
  temp <- Bcell_cluster_condition_donor_count[Bcell_cluster_condition_donor_count$clusters == cluster_name,]
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
        row_select <- (Bcell_LAIV_cluster_ttest_age_group_paired$clusters == cluster_name) &
          (Bcell_LAIV_cluster_ttest_age_group_paired$age_group1 == temp_age_group1) &
          (Bcell_LAIV_cluster_ttest_age_group_paired$age_group2 == temp_age_group2) &
          (Bcell_LAIV_cluster_ttest_age_group_paired$day_group == day_group_name)
        Bcell_LAIV_cluster_ttest_age_group_paired$n1[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == temp_age_group1]))
        Bcell_LAIV_cluster_ttest_age_group_paired$n2[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == temp_age_group2]))
        Bcell_LAIV_cluster_ttest_age_group_paired$mean1[row_select] <- ttest_result[["estimate"]][["mean of x"]]
        Bcell_LAIV_cluster_ttest_age_group_paired$mean2[row_select] <- ttest_result[["estimate"]][["mean of y"]]
        Bcell_LAIV_cluster_ttest_age_group_paired$pval[row_select] <- ttest_result$p.value
        Bcell_LAIV_cluster_ttest_age_group_paired$sd[row_select] <- sd(temp1$percentage[(temp1$age_group == temp_age_group1) |
                                                                                          (temp1$age_group == temp_age_group2)])
      }
    }

    # }

  }
}

Bcell_LAIV_cluster_ttest_age_group_paired <- Bcell_LAIV_cluster_ttest_age_group_paired %>% arrange(pval)
Bcell_LAIV_cluster_ttest_age_group_paired$diff <- Bcell_LAIV_cluster_ttest_age_group_paired$mean2 - 
  Bcell_LAIV_cluster_ttest_age_group_paired$mean1

Bcell_LAIV_cluster_ttest_age_group_paired$fdr <- p.adjust(Bcell_LAIV_cluster_ttest_age_group_paired$pval,method = 'fdr')
Bcell_LAIV_cluster_ttest_age_group_paired <- Bcell_LAIV_cluster_ttest_age_group_paired %>% arrange(fdr)
Bcell_LAIV_cluster_ttest_age_group_paired$FC <- Bcell_LAIV_cluster_ttest_age_group_paired$diff/Bcell_LAIV_cluster_ttest_age_group_paired$sd
Bcell_LAIV_cluster_ttest_age_group_paired <- Bcell_LAIV_cluster_ttest_age_group_paired[Bcell_LAIV_cluster_ttest_age_group_paired$cluster != 'PB',]
write.xlsx(Bcell_LAIV_cluster_ttest_age_group_paired,'tonsil_LAIV_120a_s_Bcell_LAIV_age_BCL6_pos_age_group_ttest_paired.xlsx')


###### prepare for the correlation network ############################
library(reshape2)
temp <- Bcell_cluster_condition_donor_count[!((Bcell_cluster_condition_donor_count$donor_ID == '02yrs F IMD030') & (Bcell_cluster_condition_donor_count$days == 12)),]
temp <- temp[!((temp$donor_ID == '33yrs M VIP024') & (temp$days == 12)),]
temp <- temp[temp$donor_ID != '02yrs M IMD085',]
temp <- dcast(temp[,c('clusters','day_group','donor_ID','percentage')], clusters + day_group ~ donor_ID, value.var="percentage")
temp <- temp[temp$day_group %in% LAIV_day_group_list,]

temp_final <- merge(temp,Bcell_LAIV_cluster_ttest_age_group_paired,by = c('clusters','day_group'))
temp_final <- temp_final[(!is.na(temp_final$n1)) & (!is.na(temp_final$n2)) & (!is.na(temp_final$pval)),]
temp_final <- temp_final[rowSums(temp_final[,3:11] == 0) < 9,]
temp_final$cluster <- as.character(temp_final$cluster)
temp_final$fdr <- p.adjust(temp_final$pval,method = 'fdr')
temp_final <- temp_final %>% arrange(fdr)
temp_final$sig <- 'no'
temp_final$sig[temp_final$pval <= 0.05] <- 'yes'
rownames(temp_final) <- paste(temp_final$day_group,temp_final$clusters,sep = '_')

write.xlsx(temp_final,'tonsil_LAIV_120a_s_Bcell_log1p_cca_CycleRegressOut_cluster_annotated_LAIV_fraction_cleaned_CytoScape.xlsx')


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
Bcell_LAIV_cluster_ttest_age_group <- expand.grid(cluster = annotate_cluster_list, age_group = age_group_list, day_group = LAIV_day_group_list)
Bcell_LAIV_cluster_ttest_age_group$diff <- NaN
Bcell_LAIV_cluster_ttest_age_group$pval <- NaN
Bcell_LAIV_cluster_ttest_age_group$n1 <- NaN
Bcell_LAIV_cluster_ttest_age_group$n2 <- NaN
Bcell_LAIV_cluster_ttest_age_group$mean1 <- NaN
Bcell_LAIV_cluster_ttest_age_group$mean2 <- NaN
Bcell_LAIV_cluster_correlate_age <- expand.grid(cluster = annotate_cluster_list, day_group = LAIV_day_group_list)
Bcell_LAIV_cluster_correlate_age$r <- NaN
Bcell_LAIV_cluster_correlate_age$pval <- NaN
Bcell_LAIV_cluster_correlate_age$n <- NaN

library(dplyr)
for (cluster_name in annotate_cluster_list) {
  print(cluster_name)
  graphics.off()
  temp <- Bcell_cluster_condition_donor_count[Bcell_cluster_condition_donor_count$clusters == cluster_name,]
  rownames(temp) <- temp$condition
  # if ((!all(temp$percentage == 0)) & (!all(temp$percentage == 100))) {
  for (day_group_name in LAIV_day_group_list) {
    temp1 <- temp[temp$days %in% LAIV_day_group_dict[[day_group_name]],]
    if (!all(temp1$percentage == 0)) {
      for (age_group_index in c(1:3)) {
        age_group_name <- age_group_list[age_group_index]
        ttest_result <- t.test(temp1$percentage[temp1$age_group == age_group_name],
                               temp1$percentage[temp1$age_group != age_group_name])
        row_select <- (Bcell_LAIV_cluster_ttest_age_group$cluster == cluster_name) & 
          (Bcell_LAIV_cluster_ttest_age_group$age_group == age_group_list[age_group_index]) &
          (Bcell_LAIV_cluster_ttest_age_group$day_group == day_group_name)
        Bcell_LAIV_cluster_ttest_age_group$n1[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == age_group_name]))
        Bcell_LAIV_cluster_ttest_age_group$n2[row_select] <- length(unique(temp1$donor_ID[temp1$age_group != age_group_name]))
        Bcell_LAIV_cluster_ttest_age_group$mean1[row_select] <- ttest_result[["estimate"]][["mean of x"]]
        Bcell_LAIV_cluster_ttest_age_group$mean2[row_select] <- ttest_result[["estimate"]][["mean of y"]]
        Bcell_LAIV_cluster_ttest_age_group$pval[row_select] <- ttest_result$p.value
      }
      row_select <- (Bcell_LAIV_cluster_correlate_age$cluster == cluster_name) &
        (Bcell_LAIV_cluster_correlate_age$day_group == day_group_name)
      cor_test <- cor.test(temp1$percentage, temp1$age)
      Bcell_LAIV_cluster_correlate_age$n[row_select] <- length(unique(temp$donor_ID))
      Bcell_LAIV_cluster_correlate_age$r[row_select] <- cor_test[["estimate"]][["cor"]]
      Bcell_LAIV_cluster_correlate_age$pval[row_select] <- cor_test$p.value  
    }  
    
    # }
    
  }
}

Bcell_LAIV_cluster_ttest_age_group <- Bcell_LAIV_cluster_ttest_age_group %>% arrange(pval)
Bcell_LAIV_cluster_correlate_age <- Bcell_LAIV_cluster_correlate_age %>% arrange(pval)
Bcell_LAIV_cluster_ttest_age_group$diff <- Bcell_LAIV_cluster_ttest_age_group$mean1 - Bcell_LAIV_cluster_ttest_age_group$mean2

# Bcell_LAIV_cluster_ttest_age_group$pval_adj <- p.adjust(Bcell_LAIV_cluster_ttest_age_group$pval,method = 'fdr')
# Bcell_LAIV_cluster_correlate_age$pval_adj <- p.adjust(Bcell_LAIV_cluster_correlate_age$pval,method = 'fdr')
Bcell_LAIV_cluster_ttest_age_group <- Bcell_LAIV_cluster_ttest_age_group %>% arrange(cluster,pval)
Bcell_LAIV_cluster_correlate_age <- Bcell_LAIV_cluster_correlate_age %>% arrange(cluster,desc(r))

Bcell_LAIV_cluster_ttest_age_group <- Bcell_LAIV_cluster_ttest_age_group[Bcell_LAIV_cluster_ttest_age_group$age_group != '07-09yrs',]
Bcell_LAIV_cluster_correlate_age <- Bcell_LAIV_cluster_correlate_age[Bcell_LAIV_cluster_correlate_age$cluster != 'PB',]

write.xlsx(Bcell_LAIV_cluster_ttest_age_group,'tonsil_LAIV_120a_s_Bcell_lognorm_LAIV_age_cluster_annotated_age_group_ttest.xlsx')
write.xlsx(Bcell_LAIV_cluster_correlate_age,'tonsil_LAIV_120a_s_Bcell_lognorm_LAIV_age_cluster_annotated_age_corr.xlsx')


####### fraction: visualization ####################################
cluster_name <- 'AbCD11c+TBX21+B'
day_group_list <- c(0,4,10)
day_group_name <- 'day0'
age_group_name <- '27-39yrs'
marker_name <- 'CXCR5'
x_lab <- 2.2
plot_corr <- 0
graphics.off()
temp <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$clusters == cluster_name),]
temp <- temp[temp$days != 8,]
temp <- temp[temp$donor_ID != "02yrs M IMD085",]
temp$day_group <- temp$days
temp$day_group[temp$days %in% c(0)] <- 'day0'
temp$day_group[temp$days %in% c(4)] <- 'day04'
temp$day_group[temp$days %in% c(6,7)] <- 'day06-07'
temp$day_group[temp$days %in% c(10)] <- 'day10'
temp$day_group[temp$days %in% c(12,14)] <- 'day12-14'

temp_donorID_batch_list <- unlist(batch_donorID_list[which(sapply(batch_marker_list, function(x) (marker_name %in% x)))])
for (day_group_name in setdiff(LAIV_day_group_list,c('day07-08','day07'))) {
  Bcell_cluster_condition_donor_count$day_group[Bcell_cluster_condition_donor_count$days %in% LAIV_day_group_dict[[day_group_name]]] <-
    day_group_name
}

# for (cluster_name in cluster_name){#annotate_cluster_list){
#   print(cluster_name)
#   temp <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$clusters == cluster_name) &
#                                                 # (Bcell_cluster_condition_donor_count$days %in% LAIV_day_group_dict[[day_group_name]]) &
#                                                 (Bcell_cluster_condition_donor_count$day_group %in% setdiff(LAIV_day_group_list,c('day07-08','day07'))) &
#                                                 Bcell_cluster_condition_donor_count$donor_ID %in% temp_donorID_batch_list,]
# 
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
    # ylim(0,max(temp$percentage) + 0.5*abs(max(temp$percentage))) +
    RotatedAxis()
  print(plot)
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_AbCD11cTBET_pos_age_group_ttest_paired.pdf',sep = ''),width = 6, height = 3)

# }

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
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_cluster_annotated_age_corr_',cluster_name,'.pdf',sep = ''),width = 6, height = 3)

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

dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_cluster_annotated_age_corr_',cluster_name,'_',day_group_name,'.pdf',sep = ''),width = 3, height = 3)


######## count: compare each age_group pairs ###########################
p <- match(Bcell_cluster_condition_donor_count$condition, flow_count_table$condition)
temp <- flow_count_table[p,]
Bcell_cluster_condition_donor_count$Bcell_count <- temp$B_cells
Bcell_cluster_condition_donor_count$cluster_count <- Bcell_cluster_condition_donor_count$percentage*Bcell_cluster_condition_donor_count$Bcell_count/100
Bcell_cluster_condition_donor_count_cleaned <- Bcell_cluster_condition_donor_count[Bcell_cluster_condition_donor_count$donor_ID != '07yrs M IMD170',]
Bcell_LAIV_cluster_ttest_age_group_paired <- expand.grid(cluster = annotate_cluster_list, age_group1 = age_group_list[1:2], age_group2 = age_group_list[2:3],
                                                         day_group = LAIV_day_group_list)
Bcell_LAIV_cluster_ttest_age_group_paired <- transform(Bcell_LAIV_cluster_ttest_age_group_paired, age_group1 = as.character(age_group1),
                                                       age_group2 = as.character(age_group2))
Bcell_LAIV_cluster_ttest_age_group_paired <- Bcell_LAIV_cluster_ttest_age_group_paired[Bcell_LAIV_cluster_ttest_age_group_paired$age_group1 != Bcell_LAIV_cluster_ttest_age_group_paired$age_group2,]
# Bcell_LAIV_cluster_ttest_age_group_paired$diff <- NaN
# Bcell_LAIV_cluster_ttest_age_group_paired$pval <- NaN
# Bcell_LAIV_cluster_ttest_age_group_paired$n1 <- NaN
# Bcell_LAIV_cluster_ttest_age_group_paired$n2 <- NaN
# Bcell_LAIV_cluster_ttest_age_group_paired$mean1 <- NaN
# Bcell_LAIV_cluster_ttest_age_group_paired$mean2 <- NaN

age_pair_list <- data.frame(age_group1 = c('02-03yrs','02-03yrs','07-09yrs'),
                            age_group2 = c('07-09yrs','27-39yrs','27-39yrs'))
# 
# library(dplyr)
# for (cluster_name in annotate_cluster_list) {
#   print(cluster_name)
#   graphics.off()
#   temp <- Bcell_cluster_condition_donor_count_cleaned[Bcell_cluster_condition_donor_count_cleaned$clusters == cluster_name,]
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
#         row_select <- (Bcell_LAIV_cluster_ttest_age_group_paired$cluster == cluster_name) & 
#           (Bcell_LAIV_cluster_ttest_age_group_paired$age_group1 == temp_age_group1) &
#           (Bcell_LAIV_cluster_ttest_age_group_paired$age_group2 == temp_age_group2) &
#           (Bcell_LAIV_cluster_ttest_age_group_paired$day_group == day_group_name)
#         Bcell_LAIV_cluster_ttest_age_group_paired$n1[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == temp_age_group1]))
#         Bcell_LAIV_cluster_ttest_age_group_paired$n2[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == temp_age_group2]))
#         Bcell_LAIV_cluster_ttest_age_group_paired$mean1[row_select] <- ttest_result[["estimate"]][["mean of x"]]
#         Bcell_LAIV_cluster_ttest_age_group_paired$mean2[row_select] <- ttest_result[["estimate"]][["mean of y"]]
#         Bcell_LAIV_cluster_ttest_age_group_paired$pval[row_select] <- ttest_result$p.value
#       }
#     }  
#     
#     # }
#     
#   }
# }
# 
# Bcell_LAIV_cluster_ttest_age_group_paired <- Bcell_LAIV_cluster_ttest_age_group_paired %>% arrange(pval)
# Bcell_LAIV_cluster_ttest_age_group_paired$diff <- Bcell_LAIV_cluster_ttest_age_group_paired$mean1 - Bcell_LAIV_cluster_ttest_age_group_paired$mean2
# 
# # Bcell_LAIV_cluster_ttest_age_group_paired$pval_adj <- p.adjust(Bcell_LAIV_cluster_ttest_age_group_paired$pval,method = 'fdr')
# Bcell_LAIV_cluster_ttest_age_group_paired <- Bcell_LAIV_cluster_ttest_age_group_paired %>% arrange(cluster,pval)
# 
# Bcell_LAIV_cluster_ttest_age_group_paired <- Bcell_LAIV_cluster_ttest_age_group_paired[Bcell_LAIV_cluster_ttest_age_group_paired$cluster != 'PB',]
# write.xlsx(Bcell_LAIV_cluster_ttest_age_group_paired,'tonsil_LAIV_120a_s_Bcell_LAIV_age_cluster_annotated_age_group_ttest_count_paired.xlsx')
# 

######## count: all cell FC diff visualization ###########################
Bcell_LAIV_cluster_ttest_age_group <- expand.grid(cluster = annotate_cluster_list, age_group = age_group_list, day_group = LAIV_day_group_list)
Bcell_LAIV_cluster_ttest_age_group$diff <- NaN
Bcell_LAIV_cluster_ttest_age_group$pval <- NaN
Bcell_LAIV_cluster_ttest_age_group$n1 <- NaN
Bcell_LAIV_cluster_ttest_age_group$n2 <- NaN
Bcell_LAIV_cluster_ttest_age_group$mean1 <- NaN
Bcell_LAIV_cluster_ttest_age_group$mean2 <- NaN
Bcell_LAIV_cluster_correlate_age <- expand.grid(cluster = annotate_cluster_list, day_group = LAIV_day_group_list)
Bcell_LAIV_cluster_correlate_age$r <- NaN
Bcell_LAIV_cluster_correlate_age$pval <- NaN
Bcell_LAIV_cluster_correlate_age$n <- NaN

library(dplyr)
for (cluster_name in annotate_cluster_list) {
  print(cluster_name)
  graphics.off()
  temp <- Bcell_cluster_condition_donor_count_cleaned[Bcell_cluster_condition_donor_count_cleaned$clusters == cluster_name,]
  rownames(temp) <- temp$condition
  # if ((!all(temp$cluster_count == 0)) & (!all(temp$cluster_count == 100))) {
  for (day_group_name in LAIV_day_group_list) {
    temp1 <- temp[temp$days %in% LAIV_day_group_dict[[day_group_name]],]
    if (!all(temp1$cluster_count == 0)) {
      for (age_group_index in c(1:3)) {
        age_group_name <- age_group_list[age_group_index]
        ttest_result <- t.test(temp1$cluster_count[temp1$age_group == age_group_name],
                               temp1$cluster_count[temp1$age_group != age_group_name])
        row_select <- (Bcell_LAIV_cluster_ttest_age_group$cluster == cluster_name) & 
          (Bcell_LAIV_cluster_ttest_age_group$age_group == age_group_list[age_group_index]) &
          (Bcell_LAIV_cluster_ttest_age_group$day_group == day_group_name)
        Bcell_LAIV_cluster_ttest_age_group$n1[row_select] <- length(unique(temp1$donor_ID[temp1$age_group == age_group_name]))
        Bcell_LAIV_cluster_ttest_age_group$n2[row_select] <- length(unique(temp1$donor_ID[temp1$age_group != age_group_name]))
        Bcell_LAIV_cluster_ttest_age_group$mean1[row_select] <- ttest_result[["estimate"]][["mean of x"]]
        Bcell_LAIV_cluster_ttest_age_group$mean2[row_select] <- ttest_result[["estimate"]][["mean of y"]]
        Bcell_LAIV_cluster_ttest_age_group$pval[row_select] <- ttest_result$p.value
      }
      row_select <- (Bcell_LAIV_cluster_correlate_age$cluster == cluster_name) &
        (Bcell_LAIV_cluster_correlate_age$day_group == day_group_name)
      cor_test <- cor.test(temp1$cluster_count, temp1$age)
      Bcell_LAIV_cluster_correlate_age$n[row_select] <- length(unique(temp$donor_ID))
      Bcell_LAIV_cluster_correlate_age$r[row_select] <- cor_test[["estimate"]][["cor"]]
      Bcell_LAIV_cluster_correlate_age$pval[row_select] <- cor_test$p.value  
    }  
    
    # }
    
  }
}

Bcell_LAIV_cluster_ttest_age_group <- Bcell_LAIV_cluster_ttest_age_group %>% arrange(pval)
Bcell_LAIV_cluster_correlate_age <- Bcell_LAIV_cluster_correlate_age %>% arrange(pval)
Bcell_LAIV_cluster_ttest_age_group$diff <- Bcell_LAIV_cluster_ttest_age_group$mean1 - Bcell_LAIV_cluster_ttest_age_group$mean2

# Bcell_LAIV_cluster_ttest_age_group$pval_adj <- p.adjust(Bcell_LAIV_cluster_ttest_age_group$pval,method = 'fdr')
# Bcell_LAIV_cluster_correlate_age$pval_adj <- p.adjust(Bcell_LAIV_cluster_correlate_age$pval,method = 'fdr')
Bcell_LAIV_cluster_ttest_age_group <- Bcell_LAIV_cluster_ttest_age_group %>% arrange(cluster,pval)
Bcell_LAIV_cluster_correlate_age <- Bcell_LAIV_cluster_correlate_age %>% arrange(cluster,desc(r))

Bcell_LAIV_cluster_ttest_age_group <- Bcell_LAIV_cluster_ttest_age_group[Bcell_LAIV_cluster_ttest_age_group$cluster != 'PB',]
Bcell_LAIV_cluster_correlate_age <- Bcell_LAIV_cluster_correlate_age[Bcell_LAIV_cluster_correlate_age$cluster != 'PB',]

Bcell_LAIV_cluster_ttest_age_group <- Bcell_LAIV_cluster_ttest_age_group[Bcell_LAIV_cluster_ttest_age_group$age_group != '07-09yrs',]
write.xlsx(Bcell_LAIV_cluster_ttest_age_group,'tonsil_LAIV_120a_s_Bcell_LAIV_age_cluster_annotated_age_group_ttest_count.xlsx')
write.xlsx(Bcell_LAIV_cluster_correlate_age,'tonsil_LAIV_120a_s_Bcell_LAIV_age_cluster_annotated_age_corr_count.xlsx')
####### count: visualization ####################################
cluster_name <- 'IGHA1-secreted'
day_group_list <- 10
day_group_name <- 'day10'
age_group_name <- '27-39yrs'
marker_name <- 'CXCR5'
x_lab <- 2.2
plot_corr <- 0
graphics.off()
temp_donorID_batch_list <- unlist(batch_donorID_list[which(sapply(batch_marker_list, function(x) (marker_name %in% x)))])
for (day_group_name in setdiff(LAIV_day_group_list,c('day07-08','day07'))) {
  Bcell_cluster_condition_donor_count_cleaned$day_group[Bcell_cluster_condition_donor_count_cleaned$days %in% LAIV_day_group_dict[[day_group_name]]] <-
    day_group_name
}
for (cluster_name in cluster_name){#annotate_cluster_list){
  print(cluster_name)
  temp <- Bcell_cluster_condition_donor_count_cleaned[(Bcell_cluster_condition_donor_count_cleaned$clusters == cluster_name) &
                                                # (Bcell_cluster_condition_donor_count_cleaned$days %in% LAIV_day_group_dict[[day_group_name]]) &
                                                (Bcell_cluster_condition_donor_count_cleaned$day_group %in% setdiff(LAIV_day_group_list,c('day07-08','day07'))) &
                                                Bcell_cluster_condition_donor_count_cleaned$donor_ID %in% temp_donorID_batch_list,]
  
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
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_PB_subisotype_age_group_count_',cluster_name,'.pdf',sep = ''),width = 6, height = 3)
}

cluster_name <- 'IGHG2-secreted'
day_group_list <- 10
day_group_name <- 'day10'
x_lab <- 2.2
graphics.off()

for (cluster_name in cluster_name){#annotate_cluster_list){
  print(cluster_name)
  temp <- Bcell_cluster_condition_donor_count_cleaned[(Bcell_cluster_condition_donor_count_cleaned$clusters == cluster_name) &
                                                        (Bcell_cluster_condition_donor_count_cleaned$days %in% LAIV_day_group_dict[[day_group_name]])]
  ttest_result <- t.test(temp$cluster_count[temp$age_group %in% '02-03yrs'],
                         temp$cluster_count[temp$age_group %in% '27-39yrs'])
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
  dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_PB_subisotype_age_group_count_',cluster_name,'.pdf',sep = ''),width = 6, height = 3)
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
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_cluster_annotated_age_corr_',cluster_name,'.pdf',sep = ''),width = 6, height = 3)

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

dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_cluster_annotated_age_corr_',cluster_name,'_',day_group_name,'.pdf',sep = ''),width = 3, height = 3)


##### within PB cells #####################
PB_secreted <- data.frame(t(as.matrix(Bcell.subset@assays$raw@data)),check.names = F)[,rownames(Bcell.subset)[grepl('-secreted',rownames(Bcell.subset))]]
Bcell.subset$max <- colnames(PB_secreted)[apply(PB_secreted,1,which.max)]

Bcell.database <- data.frame(cbind(Bcell.subset$donor_ID,as.character(Bcell.subset$max),Bcell.subset$condition,Bcell.subset$days,Bcell.subset$stimulation,Bcell.subset$age_group,Bcell.subset$day_group,Bcell.subset$age))
colnames(Bcell.database) <- c('donor_ID','clusters','condition','days','stimulation','age_group','day_group','age')
Bcell.database <- Bcell.database[(Bcell.subset$major_clusters == 'PB') & (Bcell.subset$stimulation != 'LAIV-'),]
# sapply(Bcell.database,class)
Bcell_cluster_condition_donor_count <- dplyr::count(Bcell.database, donor_ID, clusters,condition, days,stimulation,age_group,day_group,age)

Bcell_cluster_condition_donor_count <- Bcell_cluster_condition_donor_count %>% group_by(donor_ID, condition, days,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
Bcell_cluster_condition_donor_count$days <- as.numeric(substring(Bcell_cluster_condition_donor_count$days,first = 4,last = 5))
Bcell_cluster_condition_donor_count$percentage <- Bcell_cluster_condition_donor_count$n/Bcell_cluster_condition_donor_count$total*100
Bcell_cluster_condition_donor_count$age <- as.numeric(Bcell_cluster_condition_donor_count$age)

for (donor_name in unique(Bcell_LAIV.subset$donor_ID)){
  temp_condition_list <- sort(unique(Bcell_LAIV.subset$condition[Bcell_LAIV.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$donor_ID == donor_name) & (Bcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Bcell_cluster_condition_donor_count$clusters)){
      if_row <- ((Bcell_cluster_condition_donor_count$donor_ID == donor_name) & (Bcell_cluster_condition_donor_count$clusters == cluster_name) & (Bcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount$clusters <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        Bcell_cluster_condition_donor_count[nrow(Bcell_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(Bcell.database)

annotate_cluster_list <- unique(Bcell.subset$max)

temp <- Bcell_cluster_condition_donor_count[Bcell_cluster_condition_donor_count$donor_ID != '02yrs M IMD085',]
# temp <- temp[temp$clusters != 'IGHG3-secreted',]
ggplot(temp,aes(x=days, y=percentage,color = clusters)) +
  facet_wrap( ~ donor_ID, scales = "free_x",ncol = 3) +
  geom_point(size = 1)+ geom_line(aes(group = interaction(donor_ID,clusters))) + 
  theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  ylab('PB %') +
  # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  theme_bw()

dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_PB_LAIV_age_group_subisotype_donorID_ratio.pdf',sep = ''),width = 6, height = 5)

# mean of each cell count
library(purrr) 
library(ggalluvial)
temp <- Bcell_cluster_condition_donor_count
temp <- temp[!(temp$donor_ID %in% c('02yrs M IMD085')),]
temp <- temp[temp$days %in% c(0,4,6,14,10,7,12),]
temp <- temp %>% group_by(age_group,clusters,days) %>% mutate(avg_cluster_count = mean(n),avg_fraction = mean(percentage))
temp <- unique(temp[,c('age_group','clusters','days','avg_cluster_count','avg_fraction')])
temp$day_group <- temp$days
temp$day_group[temp$days %in% c(6,7)] <- 7
temp <- temp[(temp$days != 12) | !(temp$donor_ID %in% c('02yrs F IMD030','33yrs M VIP024')),]
temp$day_group[temp$days %in% c(12,14)] <- 14

p <- DimPlot(Bcell_LAIV.subset, group.by = 'max', label = TRUE, reduction = "umap")
pbuild <- ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # Pull the data used for the plot
pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
cols_group <- unique(pdata$colour) # Get a vector of unique colors
ordered_max_list <- c('naive','cytokine-rich','CXCR4+act.','CD23+act.','CD83+folli.','prolif. CD83+folli.','preGC/folli.',
                                    'switched GC','GC LZ','GC-like','memory','switched memory','IgA+/IgM+PB','IgG+PB')
# ordered_max_list <- c('naive','act.-like folli.','CXCR4+','act.-like','switching folli./memory','folli. memory',
#                                     'GC-like','switched memory','switched GC','GC','GC LZ','CD83+CXCR4+GC','switched act. memory','unswitched act.','prolif. PB-like GC','IgM+PB','IgA+PB','IgG+PB')

names(cols_group) <- random_max_list
# Idents(TB_shrink.seurat_update) <- TB_shrink.seurat_update@meta.data$sample_tag
temp$clusters <- factor(temp$clusters,levels = ordered_max_list)
ggplot(temp, aes(x = days, stratum = age_group, alluvium = clusters, y = avg_cluster_count, label = clusters,fill = clusters)) +#, 
  geom_alluvium(aes(fill = clusters), width = 1,alpha = 1,color = 'black',size = 0.1) +
  facet_grid(age_group ~ .,scales = "fixed") +
  # geom_stratum(width = 1/12, fill = "black",size = 0.5) +
  # scale_fill_manual(values = cols_group) +
  ggtitle('') + ylab('#cells') + #ylim(0,300) + xlim(0,14)
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_cluster_annotated_count.pdf',sep = ''),width = 5, height = 5)

ggplot(temp, aes(x = days, stratum = age_group, alluvium = clusters, y = avg_fraction, label = clusters,fill = clusters)) +#, 
  geom_alluvium(aes(fill = clusters), width = 0.5,alpha = 1,color = 'black',size = 0.1) +
  facet_grid(age_group ~ .,scales = "fixed") +
  # geom_stratum(width = 1/12, fill = "black") +
  # scale_fill_manual(values = cols_group) +
  ggtitle('') + ylab('% in B cells') + #ylim(0,300) + xlim(0,14)
  theme(text = element_text(size = 10),plot.title = element_text(size = 10, face = "bold")) 
dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_max_fraction.pdf',sep = ''),width = 5, height = 5)

######## PB cell fraction comparison visualization ###########################
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
  temp <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$clusters == cluster_name) & 
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
cluster_correlate_age <- cluster_correlate_age %>% arrange(pval)
write.xlsx(cluster_ttest_age_group,'tonsil_LAIV_120a_s_PB_subisotype_age_group_ttest.xlsx')
write.xlsx(cluster_correlate_age,'tonsil_LAIV_120a_s_PB_subisotype_age_corr.xlsx')

########## cell number #####################
flow_count_table <- read.xlsx('tonsil_LAIV_120a_s_flow_cell_count.xlsx', sheetName = "Sheet1")
flow_count_table <- flow_count_table[,!colnames(flow_count_table) == '...1']
flow_count_table$condition <- gsub('day0 LAIV-','day0',flow_count_table$condition)
flow_count_table$B_cells[grepl('IMD170',flow_count_table$donor_ID)] <- 0

Bcell.database <- data.frame(cbind(Bcell_LAIV.subset$donor_ID,as.character(Bcell_LAIV.subset$cluster_annotated),Bcell_LAIV.subset$condition,Bcell_LAIV.subset$days,Bcell_LAIV.subset$stimulation,Bcell_LAIV.subset$age_group,Bcell_LAIV.subset$day_group,Bcell_LAIV.subset$age))
colnames(Bcell.database) <- c('donor_ID','clusters','condition','days','stimulation','age_group','day_group','age')
# sapply(Bcell.database,class)
Bcell.database <- Bcell.database[Bcell.database$clusters == 'memory B',]
Bcell_cluster_condition_donor_count <- dplyr::count(Bcell.database, donor_ID, clusters,condition, days,stimulation,age_group,day_group,age)
Bcell_cluster_condition_donor_count <- Bcell_cluster_condition_donor_count %>% group_by(donor_ID, condition, days,stimulation,age_group,day_group,age) %>% mutate(total = sum(n))
Bcell_cluster_condition_donor_count$days <- as.numeric(substring(Bcell_cluster_condition_donor_count$days,first = 4,last = 5))
Bcell_cluster_condition_donor_count$percentage <- Bcell_cluster_condition_donor_count$n/Bcell_cluster_condition_donor_count$total*100
Bcell_cluster_condition_donor_count$age <- as.numeric(Bcell_cluster_condition_donor_count$age)
# Filling zero values for the condition that has no counts!!
for (donor_name in unique(Bcell_LAIV.subset$donor_ID)){
  temp_condition_list <- sort(unique(Bcell_LAIV.subset$condition[Bcell_LAIV.subset$donor_ID == donor_name]))
  for (condition_name in temp_condition_list){
    temp_Tcount <- Bcell_cluster_condition_donor_count[(Bcell_cluster_condition_donor_count$donor_ID == donor_name) & (Bcell_cluster_condition_donor_count$condition == condition_name),][1,]
    for (cluster_name in unique(Bcell_LAIV.subset$cluster_annotated)){
      if_row <- ((Bcell_cluster_condition_donor_count$donor_ID == donor_name) & (Bcell_cluster_condition_donor_count$clusters == cluster_name) & (Bcell_cluster_condition_donor_count$condition == condition_name))
      if (sum(if_row) == 0){
        temp_Tcount$clusters <- cluster_name
        temp_Tcount$n <- 0
        temp_Tcount$percentage <- 0
        Bcell_cluster_condition_donor_count[nrow(Bcell_cluster_condition_donor_count) + 1,] <- temp_Tcount
      }
    }
  }
}
rm(Bcell.database)

Bcell_cluster_condition_donor_count_cleaned <- Bcell_cluster_condition_donor_count 
p <- match(Bcell_cluster_condition_donor_count$condition, flow_count_table$condition)
temp <- flow_count_table[p,]
Bcell_cluster_condition_donor_count$Bcell_count <- temp$B_cells
Bcell_cluster_condition_donor_count$cluster_count <- Bcell_cluster_condition_donor_count$percentage*Bcell_cluster_condition_donor_count$Bcell_count/100
temp <- Bcell_cluster_condition_donor_count[Bcell_cluster_condition_donor_count$donor_ID != '02yrs M IMD085',]
temp <- temp[temp$clusters != 'IGHG3-secreted',]
temp <- temp[temp$clusters != 'B cells',]
ggplot(temp,aes(x=days, y=cluster_count,color = clusters)) +
  facet_wrap( ~ donor_ID, scales = "free_x",ncol = 3) +
  geom_point(size = 1)+ geom_line(aes(group = interaction(donor_ID,clusters))) + 
  theme(text = element_text(size = 15),plot.title = element_text(size = 18, face = "bold")) + RotatedAxis() +
  ylab('#Bcells') +
  # scale_colour_manual(values = c("LAIV+" = "red", "LAIV-" = "darkgrey")) +
  theme_bw()

dev.print(pdf, paste('tonsil_LAIV_120a_s_Bcell_LAIV_age_group_cluster_annotated_donorID_count.pdf',sep = ''),width = 6, height = 5)
