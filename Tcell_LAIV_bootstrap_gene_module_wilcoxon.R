library(dplyr)
library(Seurat)
library(gridExtra)
library(grid)
library(gtable)

cor_cluster_table <- readRDS('Tcell_LAIV_module_cor_cluster_table.rds')
Tcell_module_dataframe <- readRDS('Tcell_LAIV_gene_module_table.rds')
batch_marker_list <- list()
batch_marker_list[[1]] <- readRDS('tonsil_LAIV_120a.rds')
batch_marker_list[[2]] <- readRDS('tonsil_LAIV_120b.rds')
batch_marker_list[[3]] <- readRDS('tonsil_LAIV_120c.rds')
batch_marker_list[[4]] <- readRDS('tonsil_LAIV_120j.rds')
minimal_marker_list <- Reduce(intersect,batch_marker_list[2:4])
module_list <- unique(cor_cluster_table$name[cor_cluster_table$name != ''])
  
age_group_list <- sort(unique(Tcell_module_dataframe$age_group))
LAIV_day_group_list <- c(0,4,7,10,14)
age_pair_list <- data.frame(age_group1 = c('02-03yrs','02-03yrs','07-09yrs'),
                            age_group2 = c('07-09yrs','27-39yrs','27-39yrs'))
n_cells <- 300 #from each donor
n_bootstrap <- 50

condition_list <- unique(Tcell_module_dataframe$condition)

library(foreach)
my.cluster <- parallel::makeCluster(16)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
Tcells_diff_wilcoxon <- foreach (index = 1:n_bootstrap,.combine = "rbind") %dopar% {#,.options.snow = opts
  temp_lognorm_gene_data_all <- data.frame()
  temp_output <- expand.grid(module = module_list, day_group = LAIV_day_group_list,
                             age_group1 = age_group_list[1:2], age_group2 = age_group_list[2:3])
  temp_output <- transform(temp_output, age_group1 = as.character(age_group1),
                           age_group2 = as.character(age_group2))
  temp_output <- temp_output[temp_output$age_group1 != temp_output$age_group2,]
  temp_output$n_feature <- NaN
  temp_output$n_min_feature <- NaN
  temp_output$pval <- NaN
  temp_output$mean1 <- NaN
  temp_output$mean2 <- NaN
  temp_output$std <- NaN
  temp_output$min_pval <- NaN
  temp_output$min_mean1 <- NaN
  temp_output$min_mean2 <- NaN
  temp_output$min_std <- NaN
  temp_output$bootstrap_index <- index
  for (condition_name in condition_list){
    temp <- which(Tcell_module_dataframe$condition == condition_name)
    donorname <- substring(condition_name, first = 1, last = 14)
    temp <- sample(temp, size = n_cells,replace = T)
    temp_lognorm_sample_diff <- data.frame(condition = condition_name,index = temp)
    temp_lognorm_gene_data_all <- rbind(temp_lognorm_gene_data_all,temp_lognorm_sample_diff)
  }
  for (module_name in module_list) {
    feature_list <- cor_cluster_table$marker[cor_cluster_table$name == module_name]
    min_feature_list <- feature_list[feature_list %in% minimal_marker_list]
    row_select <- (temp_output$module == module_name)
    if (length(min_feature_list) > 0){
      temp_lognorm_gene_data_all$min_expr <- Tcell_module_dataframe[,paste('min_',module_name,sep = '')][temp_lognorm_gene_data_all$index]
      temp_output$n_min_feature[row_select] <- length(min_feature_list)
    } else {
      temp_output$n_min_feature[row_select] <- 0
    }
    temp_lognorm_gene_data_all$expr <- Tcell_module_dataframe[,module_name][temp_lognorm_gene_data_all$index]
    temp_lognorm_gene_data_all$age_group <- Tcell_module_dataframe$age_group[temp_lognorm_gene_data_all$index]
    temp_lognorm_gene_data_all$day_group <- Tcell_module_dataframe$day_group[temp_lognorm_gene_data_all$index]
    temp_output$n_feature[row_select] <- length(feature_list)
    for (day_name in LAIV_day_group_list) {
      # print(day_name)
      temp_LAIV <- temp_lognorm_gene_data_all[temp_lognorm_gene_data_all$day_group == day_name,]
      for (age_group_index in c(1:3)) {
        temp_age_group1 <- age_pair_list$age_group1[age_group_index]
        temp_age_group2 <- age_pair_list$age_group2[age_group_index]
        row_select <- (temp_output$module == module_name) & 
          (temp_output$age_group1 == temp_age_group1) &
          (temp_output$age_group2 == temp_age_group2) &
          (temp_output$day_group == day_name)
        wilcox_result <- wilcox.test(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group1],
                                     temp_LAIV$expr[temp_LAIV$age_group == temp_age_group2])
        temp_output$pval[row_select] <- wilcox_result$p.value
        temp_output$mean1[row_select] <- mean(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group1])
        temp_output$mean2[row_select] <- mean(temp_LAIV$expr[temp_LAIV$age_group == temp_age_group2])
        temp_output$std[row_select] <- sd(temp_LAIV$expr[(temp_LAIV$age_group == temp_age_group1) |
                                                           (temp_LAIV$age_group == temp_age_group2)])
        if (length(min_feature_list) > 0){
          wilcox_result <- wilcox.test(temp_LAIV$min_expr[temp_LAIV$age_group == temp_age_group1],
                                       temp_LAIV$min_expr[temp_LAIV$age_group == temp_age_group2])
          temp_output$min_pval[row_select] <- wilcox_result$p.value
          temp_output$min_mean1[row_select] <- mean(temp_LAIV$min_expr[temp_LAIV$age_group == temp_age_group1])
          temp_output$min_mean2[row_select] <- mean(temp_LAIV$min_expr[temp_LAIV$age_group == temp_age_group2])
          temp_output$min_std[row_select] <- sd(temp_LAIV$min_expr[(temp_LAIV$age_group == temp_age_group1) |
                                                                     (temp_LAIV$age_group == temp_age_group2)])
        }
      }
    } 
  }
  return(temp_output)
}
saveRDS(Tcells_diff_wilcoxon,'tonsil_LAIV_120a_s_Tcell_LAIV_lognorm_scale_module_age_wilcoxon_output_parallel_n300.rds')
parallel::stopCluster(cl = my.cluster)
# stopCluster(cl)