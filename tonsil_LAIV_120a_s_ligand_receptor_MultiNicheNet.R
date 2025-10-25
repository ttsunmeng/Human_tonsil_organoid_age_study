library(Seurat)
library(dplyr)
library(ggplot2)
library(xlsx)
library(ggcorrplot)
`%notin%` <- Negate(`%in%`)
library(ggrepel)
library(plotly)
library(htmlwidgets)
library(forcats)
library(SingleCellExperiment)
library(multinichenetr)

# devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
library(tidyverse)
# devtools::install_github("saeyslab/multinichenetr")

Bcell_LAIV.subset <- readRDS('tonsil_LAIV_120a_s_lognorm_Bcell_LAIV.rds')
Bcell_LAIV.subset$subsets <- Bcell_LAIV.subset$major_clusters2
Bcell_LAIV.subset$celltypes <- 'Bcells'
CD4_LAIV.subset <- readRDS('tonsil_LAIV_120a_s_lognorm_CD4_LAIV.rds')
CD4_LAIV.subset <- subset(CD4_LAIV.subset,days != 'day02')
CD4_LAIV.subset <- subset(CD4_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
CD4_LAIV.subset$subsets <- CD4_LAIV.subset$major_cluster
CD4_LAIV.subset$celltypes <- 'CD4 Tcells'
CD8_LAIV.subset <- readRDS('tonsil_LAIV_120a_s_CD8.rds')
CD8_LAIV.subset <- subset(CD8_LAIV.subset, stimulation != 'LAIV-')
CD8_LAIV.subset <- subset(CD8_LAIV.subset,days != 'day02')
CD8_LAIV.subset <- subset(CD8_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
CD8_LAIV.subset[['RNA']] <- NULL
CD8_LAIV.subset[['ADT']] <- NULL
CD8_LAIV.subset[['raw']] <- NULL
CD8_LAIV.subset[['integrated']] <- NULL
CD8_LAIV.subset$subsets <- 'CD8'
gdT_LAIV.subset <- readRDS('tonsil_LAIV_120a_s_gdT.rds')
gdT_LAIV.subset <- subset(gdT_LAIV.subset, stimulation != 'LAIV-')
gdT_LAIV.subset <- subset(gdT_LAIV.subset,days != 'day02')
gdT_LAIV.subset <- subset(gdT_LAIV.subset,(condition != '33yrs M VIP024 day12 LAIV+') & (condition != '02yrs F IMD030 day12 LAIV+'))
gdT_LAIV.subset[['RNA']] <- NULL
gdT_LAIV.subset[['ADT']] <- NULL
gdT_LAIV.subset[['raw']] <- NULL
gdT_LAIV.subset[['integrated']] <- NULL
gdT_LAIV.subset$subsets <- 'gdT'

Bcell_Tcell_LAIV.subset <- merge(merge(merge(Bcell_LAIV.subset, y = CD4_LAIV.subset),CD8_LAIV.subset),gdT_LAIV.subset)

# tutorial from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
# These models are the human data
# Need to download from zenodo manually!
# ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
ligand_target_matrix = readRDS("/Users/mengsun/Library/CloudStorage/GoogleDrive-ttmsun@stanford.edu/My\ Drive/Stanford/RNA-seq/analysis_papers/ligand_receptor/ligand_target_matrix_nsga2r_final.rds")
colnames(ligand_target_matrix) = colnames(ligand_target_matrix) %>% make.names()
rownames(ligand_target_matrix) = rownames(ligand_target_matrix) %>% make.names()
lr_network = readRDS("/Users/mengsun/Library/CloudStorage/GoogleDrive-ttmsun@stanford.edu/My\ Drive/Stanford/RNA-seq/analysis_papers/ligand_receptor/lr_network_human_21122021.rds")
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor) %>% mutate(ligand = make.names(ligand), receptor = make.names(receptor))

# Prepare SingleCellExperiment Objects
Bcell_Tcell_LAIV.subset <- subset(Bcell_Tcell_LAIV.subset, (day_group != 'day0') & (age_group != '07-09yrs'))
Bcell_Tcell_LAIV.subset$groups <- 'O'
Bcell_Tcell_LAIV.subset$groups[Bcell_Tcell_LAIV.subset$age_group == '02-03yrs'] <- 'Y'
Bcell_Tcell_LAIV.subset$groups[Bcell_Tcell_LAIV.subset$age_group == '27-39yrs'] <- 'A'
Bcell_Tcell_LAIV.subset$subsets <- make.names(Bcell_Tcell_LAIV.subset$subsets)
Bcell_Tcell_LAIV.subset$donor_ID <- make.names(Bcell_Tcell_LAIV.subset$donor_ID)
Bcell_Tcell_LAIV.subset <- subset(Bcell_Tcell_LAIV.subset,subsets != "early.activated.CD4")
Bcell_Tcell_LAIV.subset <- subset(Bcell_Tcell_LAIV.subset,subsets != "switched.GC")

scaled_integrated_assay <- CreateAssayObject(counts = expm1(as.matrix(Bcell_Tcell_LAIV.subset@assays[["log1p"]]@data)))
Bcell_Tcell_LAIV.subset[["RNA"]] <- scaled_integrated_assay
Bcell_Tcell_LAIV.subset[['integrated_scale']] <- NULL

sce = Seurat::as.SingleCellExperiment(Bcell_Tcell_LAIV.subset, assay = "RNA")
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()
# temp = readRDS(url("https://zenodo.org/record/8010790/files/sce_subset_misc.rds"))
# SummarizedExperiment::colData(temp)[,"MIS.C.AgeTier"] %>% unique()
# Define metadata columns
sample_id = "donor_ID"
group_id = "groups"
celltype_id = "subsets"
covariates = NA
batches = NA

senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()

sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% c(senders_oi, receivers_oi)]

# Define the contrasts and covariates of interest for the DE analysis
min_cells = 10
# contrasts_oi = c("'Y-(O+A)/2','O-(Y+A)/2','A-(Y+O)/2'")

contrasts_oi = c("'Y-A','A-Y'")
contrast_tbl = tibble(
  contrast = c("Y-A","A-Y"), 
  group = c("A","Y"))
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% contrast_tbl$group]

# Define the parameters for the NicheNet ligand activity analysis
logFC_threshold = 0.50
p_val_threshold = 0.05
fraction_cutoff = 0.05
p_val_adj = FALSE 
empirical_pval = FALSE
top_n_target = 250
cores_system = 8
n.cores = min(cores_system, union(senders_oi, receivers_oi) %>% length()) # use one core per receiver cell type

# Define the weights of the prioritization of both expression, differential expression and NicheNet activity information
prioritizing_weights_DE = c("de_ligand" = 1,
                            "de_receptor" = 1)
prioritizing_weights_activity = c("activity_scaled" = 2)

prioritizing_weights_expression_specificity = c("exprs_ligand" = 2,
                                                "exprs_receptor" = 2)

prioritizing_weights_expression_sufficiency = c("frac_exprs_ligand_receptor" = 1)

prioritizing_weights_relative_abundance = c( "abund_sender" = 0,
                                             "abund_receiver" = 0)
prioritizing_weights = c(prioritizing_weights_DE, 
                         prioritizing_weights_activity, 
                         prioritizing_weights_expression_specificity,
                         prioritizing_weights_expression_sufficiency, 
                         prioritizing_weights_relative_abundance)

multinichenet_output = multi_nichenet_analysis(sce = sce, celltype_id = celltype_id, sample_id = sample_id, group_id = group_id, 
                                               lr_network = lr_network, ligand_target_matrix = ligand_target_matrix, contrasts_oi = contrasts_oi, contrast_tbl = contrast_tbl, batches = batches, covariates = covariates,
                                               prioritizing_weights = prioritizing_weights, min_cells = min_cells, logFC_threshold = logFC_threshold, p_val_threshold = p_val_threshold,  
                                               fraction_cutoff = fraction_cutoff, p_val_adj = p_val_adj, empirical_pval = empirical_pval, top_n_target = top_n_target, n.cores = n.cores, sender_receiver_separate = FALSE, verbose = TRUE)

multinichenet_output$prioritization_tables$group_prioritization_tbl %>% head()

prioritized_tbl_oi_all = get_top_n_lr_pairs(multinichenet_output$prioritization_tables, 50, rank_per_group = FALSE)

prioritized_tbl_oi = multinichenet_output$prioritization_tables$group_prioritization_tbl %>%
  filter(id %in% prioritized_tbl_oi_all$id) %>%
  distinct(id, sender, receiver, ligand, receptor, group) %>% left_join(prioritized_tbl_oi_all)
prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] = 0

senders_receivers = union(prioritized_tbl_oi$sender %>% unique(), prioritized_tbl_oi$receiver %>% unique()) %>% sort()

colors_sender = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)
colors_receiver = RColorBrewer::brewer.pal(n = length(senders_receivers), name = 'Spectral') %>% magrittr::set_names(senders_receivers)

circos_list = make_circos_group_comparison(prioritized_tbl_oi, colors_sender, colors_receiver)
dev.print(pdf, paste('tonsil_Rhapsody_LAIV_120a_s_raw_noday0_MultiNicheNet_Tcell_Bcells_Adults.pdf',sep = ''),width = 5, height = 5)

