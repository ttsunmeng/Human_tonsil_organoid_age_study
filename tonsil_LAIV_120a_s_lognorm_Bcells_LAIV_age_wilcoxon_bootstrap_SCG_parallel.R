library(dplyr)
library(Seurat)
library(gridExtra)
library(grid)
library(gtable)
# ########## all calculation ###################################################
Bcell.subset <- readRDS('tonsil_LAIV_120a_s_Bcell_LAIV.rds')
# Bcell.subset <- subset(Bcell.subset,days != 'day02')
batch_marker_list <- list()
batch_marker_list[[1]] <- readRDS('tonsil_LAIV_120a.rds')
batch_marker_list[[2]] <- readRDS('tonsil_LAIV_120b.rds')
batch_marker_list[[3]] <- readRDS('tonsil_LAIV_120c.rds')
batch_marker_list[[4]] <- readRDS('tonsil_LAIV_120j.rds')
all.markers <- rownames(Bcell.subset)
age_group_list <- sort(unique(Bcell.subset$age_group))
LAIV_day_group_dict <- list()
LAIV_day_group_dict[['day0']] <- c(0)
LAIV_day_group_dict[['day04']] <- c(4)
# LAIV_day_group_dict[['day07']] <- c(7)
LAIV_day_group_dict[['day10']] <- c(10)
# LAIV_day_group_dict[['day14']] <- c(14)
LAIV_day_group_dict[['day06-07']] <- c(6,7)
LAIV_day_group_dict[['day07-08']] <- c(7,8)
LAIV_day_group_dict[['day12-14']] <- c(12,14)

LAIV_day_group_list <- names(LAIV_day_group_dict)

Bcells_dataframe_lognorm <- data.frame(t(as.matrix(Bcell.subset@assays[["normalized"]]@data)),check.names = F)
Bcells_dataframe_lognorm$condition <- Bcell.subset$condition

# initializing output tables
n_cells <- 300#floor(exp(mean(log(sort(table(Bcells_dataframe_lognorm$condition))))))
print(n_cells)
# takes about 4 hours
n_bootstrap <- 50

age_pair_list <- data.frame(age_group1 = c('02-03yrs','02-03yrs','07-09yrs'),
                            age_group2 = c('07-09yrs','27-39yrs','27-39yrs'))

batch_donorID_list <- list()
batch_donorID_list[[1]] <- unique(Bcell.subset$donor_ID[Bcell.subset$batch == 'batch1'])
batch_donorID_list[[2]] <- unique(Bcell.subset$donor_ID[Bcell.subset$batch == 'batch2'])
batch_donorID_list[[3]] <- unique(Bcell.subset$donor_ID[Bcell.subset$batch == 'batch3'])
batch_donorID_list[[4]] <- unique(Bcell.subset$donor_ID[Bcell.subset$batch == 'batch4'])

library(foreach)
# library(doSNOW)
# cl <- makeCluster(6)
# registerDoSNOW(cl)
# Bcells <- txtProgressBar(max = n_bootstrap, style = 3)
# progress <- function(n) setTxtProgressBar(Bcells, n)
# opts <- list(progress = progress)

my.cluster <- parallel::makeCluster(16)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
temp_condition_list <- unique(Bcell.subset$condition)
Bcells_diff_wilcoxon <- foreach (index = 1:n_bootstrap,.combine = "rbind") %dopar% {#,.options.snow = opts
  temp_lognorm_gene_data_all <- data.frame()
  temp_output <- expand.grid(marker = all.markers, day_group = LAIV_day_group_list,
                             age_group1 = age_group_list[1:2], age_group2 = age_group_list[2:3])
  temp_output <- transform(temp_output, age_group1 = as.character(age_group1),
                           age_group2 = as.character(age_group2))
  temp_output <- temp_output[temp_output$age_group1 != temp_output$age_group2,]
  temp_output$pval <- NaN
  temp_output$n1 <- NaN
  temp_output$n2 <- NaN
  temp_output$mean1 <- NaN
  temp_output$mean2 <- NaN
  temp_output$pct1 <- NaN
  temp_output$pct2 <- NaN
  temp_output$std1 <- NaN
  temp_output$std2 <- NaN
  temp_output$log2FC <- NaN
  temp_output$bootstrap_index <- index
  for (condition_name in temp_condition_list){
    temp <- which(Bcells_dataframe_lognorm$condition == condition_name)
    donorname <- substring(condition_name, first = 1, last = 14)
    temp <- sample(temp, size = n_cells,replace = T)
    temp_lognorm_sample_diff <- data.frame(condition = condition_name,index = temp,
                                         age_group = Bcell.subset$age_group[Bcell.subset$donor_ID == donorname][1],
                                         donor_ID = donorname,
                                         days = as.numeric(substring(condition_name, first = 19, last = 20)))
    temp_lognorm_gene_data_all <- rbind(temp_lognorm_gene_data_all,temp_lognorm_sample_diff)
  }
  temp_lognorm_gene_data_all <- temp_lognorm_gene_data_all[temp_lognorm_gene_data_all$donor_ID != '02yrs M IMD085',]
  temp_lognorm_gene_data_all <- temp_lognorm_gene_data_all[(temp_lognorm_gene_data_all$days != 12) | (temp_lognorm_gene_data_all$donor_ID != '33yrs M VIP024'),]
  temp_lognorm_gene_data_all <- temp_lognorm_gene_data_all[(temp_lognorm_gene_data_all$days != 12) | (temp_lognorm_gene_data_all$donor_ID != '02yrs F IMD030'),]
  for (day_group_name in LAIV_day_group_list) {
    print(day_group_name)
    temp_lognorm_gene_data <- temp_lognorm_gene_data_all[temp_lognorm_gene_data_all$days %in% LAIV_day_group_dict[[day_group_name]],]
    for (marker_name in all.markers) {
      # print(marker_name)
      # check what are the donorID list of a particular marker
      temp_donorID_batch_list <- unlist(batch_donorID_list[which(sapply(batch_marker_list, function(x) (marker_name %in% x)))])
      temp_LAIV <- temp_lognorm_gene_data[(temp_lognorm_gene_data$donor_ID %in% temp_donorID_batch_list),]
      if (length(unique(temp_LAIV$donor_ID)) > 1) {
        temp_LAIV$expr <- Bcells_dataframe_lognorm[temp_LAIV$index,][,marker_name]
        for (age_group_index in c(1:3)) {
          temp_age_group1 <- age_pair_list$age_group1[age_group_index]
          temp_age_group2 <- age_pair_list$age_group2[age_group_index]
          wilcox_result <- wilcox.test(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group1],
                                       temp_LAIV$expr[temp_LAIV$age_group == temp_age_group2])
          row_select <- (temp_output$marker == marker_name) & 
            (temp_output$age_group1 == temp_age_group1) &
            (temp_output$age_group2 == temp_age_group2) &
            (temp_output$day_group == day_group_name)
          temp_output$n1[row_select] <- 
            length(unique(temp_LAIV$donor_ID[temp_LAIV$age_group == temp_age_group1]))
          temp_output$n2[row_select] <- 
            length(unique(temp_LAIV$donor_ID[temp_LAIV$age_group == temp_age_group2]))
          temp_output$mean1[row_select] <- mean(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group1])
          temp_output$mean2[row_select] <- mean(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group2])
          temp_output$std1[row_select] <- sd(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group1])
          temp_output$std2[row_select] <- sd(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group2])
          temp_output$pval[row_select] <- wilcox_result$p.value
          temp_output$pct1[row_select] <- sum(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group1] > 0)/sum(temp_LAIV$age_group == temp_age_group1)
          temp_output$pct2[row_select] <- sum(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group2] > 0)/sum(temp_LAIV$age_group == temp_age_group2)
        }
      } 
    } 
  }
  return(temp_output)
}
saveRDS(Bcells_diff_wilcoxon,'tonsil_LAIV_120a_s_Bcells_LAIV_lognorm_age_wilcoxon_output_parallel_n300.rds')
parallel::stopCluster(cl = my.cluster)
# stopCluster(cl)

