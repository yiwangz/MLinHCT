rm(list=ls())

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if(length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("data.table","fastDummies","caret","tidyverse")
ipak(packages)

############################################################################################
############################################################################################
load("~/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/MLinHCT/Publications/JCO_2022/Data/dt_sj_imp10.Rdata")
load("~/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/MLinHCT/Publications/JCO_2022/Data/dt_msk_imp10.Rdata")

m <- length(unique(dt_sj_imp10$.imp))
n_sj <- length(unique(dt_sj_imp10$entity_id))
n_msk <- length(unique(dt_msk_imp10$entity_id))

dt_sj_imp10 <- dt_sj_imp10[,.SD,.SDcol=which(colnames(dt_sj_imp10) %in% colnames(dt_msk_imp10))]
colnames_index <- sapply(1:ncol(dt_sj_imp10), function(i) which(colnames(dt_sj_imp10) == colnames(dt_msk_imp10)[i]))
dt_sj_imp10 <- dt_sj_imp10[,.SD,.SDcol=colnames_index]

#############################################
## death info
dt_sj_imp10[,obs_time_100day:=ifelse(num_t_surv<100, num_t_surv, 100)]
dt_sj_imp10[,delta_100day:=ifelse(cat_death=="Deceased" & num_t_surv<100, "Deceased", "Alive")]
dt_sj_imp10[,delta_100day:=factor(delta_100day, levels = c("Alive","Deceased"))]

dt_sj_imp10[,obs_time_1year:=ifelse(num_t_surv<365, num_t_surv, 365)]
dt_sj_imp10[,delta_1year:=ifelse(cat_death=="Deceased" & num_t_surv<365, "Deceased", "Alive")]
dt_sj_imp10[,delta_1year:=factor(delta_1year, levels = c("Alive","Deceased"))]

dt_sj_imp10[,obs_time_2year:=ifelse(num_t_surv<365*2, num_t_surv, 365*2)]
dt_sj_imp10[,delta_2year:=ifelse(cat_death=="Deceased" & num_t_surv<365*2, "Deceased", "Alive")]
dt_sj_imp10[,delta_2year:=factor(delta_2year, levels = c("Alive","Deceased"))]

dt_sj_imp10[,cat_dt_trans:=ifelse(as.numeric(substr(num_dt_trans,1,4))>2013,1,0)]
dt_sj_imp10[,num_dt_trans:=NULL]

table(dt_sj_imp10$delta_100day[seq(1,nrow(dt_sj_imp10),m)])  ## 60 death
table(dt_sj_imp10$delta_1year[seq(1,nrow(dt_sj_imp10),m)])  ## 187 death
table(dt_sj_imp10$delta_2year[seq(1,nrow(dt_sj_imp10),m)])  ## 233 death

############################################################################################
## recode aGVHD
agvhd_recode <- function(vec){
  vec_recode <- rep("0", length(vec))
  for(i in 1:length(vec)){
    if(vec[i] %in% c("1","2")){
      vec_recode[i] <- "1"
    }else if (vec[i] %in% c("3","4")){
      vec_recode[i] <- "2"
    }  
  }
  return(vec_recode)
}
  
dt_sj_imp10[,("max_cat_agvhd_2"):=lapply(.SD, function(i) agvhd_recode(i)), .SDcols="max_cat_agvhd_2"]

############################################################################################
## feature encoding
dt_sj_imp10_agg <- dt_sj_imp10 ##if no averaging is performed for multiple imputed datasets

feature_names <- setdiff(colnames(dt_sj_imp10_agg), c("entity_id",".id",".imp","cat_death","num_t_surv",
                                                      "obs_time_100day","delta_100day",
                                                      "obs_time_1year","delta_1year",
                                                      "obs_time_2year","delta_2year"))
cat_names <- stringr::str_subset(feature_names, "cat_")
lgl_names <- stringr::str_subset(feature_names, "lgl")
num_names <- setdiff(feature_names, c(cat_names, lgl_names))

dt_sj_imp10_dum_cat <- do.call(cbind, lapply(1:length(cat_names), function(i) fastDummies::dummy_cols(dt_sj_imp10_agg[,.SD, .SDcol=cat_names[i]])[,-1]))

# dt_all_agg <- rbind(dt_sj_imp10, dt_msk_imp10)
# num_means <- dt_all_agg[, sapply(.SD, function(i) mean(i)), .SDcol=num_names]
# num_sds <- dt_all_agg[, sapply(.SD, function(i) sd(i)), .SDcol=num_names]
# dt_sj_imp10_scale_num <- do.call(cbind, lapply(1:length(num_names), function(i) (dt_sj_imp10_agg[,.SD,.SDcol=num_names[i]]-num_means[i])/num_sds[i]))
dt_sj_imp10_scale_num <- dt_sj_imp10_agg[, lapply(.SD, function(i) scale(i, scale = TRUE, center = TRUE)), .SDcol=num_names]
colnames(dt_sj_imp10_scale_num) <- num_names

dt_sj_imp10_agg[, (lgl_names):=lapply(.SD, as.logical), .SDcols=lgl_names]

dt_sj_imp10_scale <- cbind(dt_sj_imp10_dum_cat, dt_sj_imp10_scale_num, dt_sj_imp10_agg[,.SD,.SDcol=lgl_names])

############################################################################################
## check Zero- and Near Zero-Variance Predictors
## since all the categories of a categorical variable should be forced into the model, this step is omitted for categorical variables
## but we need to perform for numeric variables
zeroVar_check_num <- list()
zeroVar_names_num <- list()

for(i in 1:m){
  zero_var_check_num <- caret::nearZeroVar(dt_sj_imp10_scale[seq(i,nrow(dt_sj_imp10_scale),m),.SD,.SDcol=num_names], saveMetrics= TRUE, uniqueCut=1)
  zeroVar_check_num[[i]] <- zero_var_check_num
  zeroVar_names_num[[i]] <- rownames(zero_var_check_num[which(zero_var_check_num$zeroVar==TRUE | zero_var_check_num$nzv==TRUE),])
}

unique(unlist(zeroVar_names_num))
# a total of 0 variables should be deleted

# zeroVar_mat_num <- do.call(rbind, lapply(1:m, function(i) zeroVar_check_num[[i]][rownames(zeroVar_check_num[[i]]) %in% unique(unlist(zeroVar_names_num)),1:2]))
# zeroVar_mat_num$var <- rownames(zeroVar_mat_num)[1:length(unique(unlist(zeroVar_names_num)))]
# zeroVar_mat_summary_num <- zeroVar_mat_num %>% group_by(var) %>% dplyr::summarise(freqRatio_mean=mean(freqRatio),
#                                                                                   percentUnique_mean=mean(percentUnique))

## identify correlated predictors
cat_names <- stringr::str_subset(colnames(dt_sj_imp10_scale), "cat_")
lgl_names <- stringr::str_subset(colnames(dt_sj_imp10_scale), "lgl")
num_names <- setdiff(colnames(dt_sj_imp10_scale), c(cat_names,lgl_names))  

high_num_names <- list()

for(i in 1:m){
  descrCor <- cor(dt_sj_imp10_scale[seq(i,nrow(dt_sj_imp10_scale),m),.SD,.SDcols=num_names])
  # summary(descrCor[upper.tri(descrCor)])
  highlyCorDescr <- caret::findCorrelation(descrCor, cutoff = .75)  ##columns to remove to reduce correlation
  high_num_names[[i]] <- num_names[highlyCorDescr]    
}

## although categorical variables will not be deleted because of zero/near-zero variance
## we still want to check which reference categorical should be deleted
zeroVar_check_cat <- list()

for(i in 1:m){
  zero_var_check_cat <- caret::nearZeroVar(dt_sj_imp10_scale[seq(i,nrow(dt_sj_imp10_scale),m),.SD,.SDcol=cat_names], saveMetrics= TRUE, uniqueCut=1)
  zeroVar_check_cat[[i]] <- zero_var_check_cat
}

zeroVar_mat_cat <- do.call(rbind, lapply(1:m, function(i) zeroVar_check_cat[[i]][,1:2]))
zeroVar_mat_cat$var <- rownames(zeroVar_mat_cat)[1:length(cat_names)]
zeroVar_mat_summary_cat <- zeroVar_mat_cat %>% group_by(var) %>% dplyr::summarise(freqRatio_mean=mean(freqRatio),
                                                                                  percentUnique_mean=mean(percentUnique))

## identify linear dependencies
## similar reason as we skip the zeroVar_names step
# linear_den_cat_names <- list()
# 
# for(i in 1:m){
#   comboInfo <- caret::findLinearCombos(data_scale[((i-1)*n+1):(i*n),.SD,.SDcols=c(cat_names, lgl_names)])
#   linear_den_cat_names[[i]] <- c(cat_names, lgl_names)[comboInfo$remove]  
# }

remove_names <- unique(c(unique(unlist(high_num_names)),
                         "cat_degree_match10_8",
                         "cat_disease_status_at_trans_Active (Malignant)",
                         "cat_donor_sex_Male",
                         "cat_dx_grp_Other Leukemias",
                         "cat_ethnicity_Not Hispanic or Latino",
                         "cat_matched_related_mismatched_unrelated",
                         "cat_prep_type_Non-Myeloablative",
                         "cat_product_type_PBSC",
                         "cat_race_Other",
                         "cat_sex_Male",
                         "max_cat_agvhd_2_2",
                         "cat_sex_dir_Female->Female",
                         "cat_dt_trans_0"))
dt_sj_imp10_scale_final <- dt_sj_imp10_scale[,.SD,.SDcols=colnames(dt_sj_imp10_scale)[!(colnames(dt_sj_imp10_scale)%in%remove_names)]]

############################################################################################
baseline_names <- setdiff(stringr::str_subset(colnames(dt_sj_imp10_scale_final), "cat_|lgl_|age|trans|sex"),
                          c("max_cat_agvhd_2_0","max_cat_agvhd_2_1"))
longitudinal_names <- setdiff(colnames(dt_sj_imp10_scale_final), baseline_names)

dt_sj_imp10_scale_final <- cbind(dt_sj_imp10_agg[,c(".id","cat_death","num_t_surv",
                                                    "obs_time_100day","delta_100day",
                                                    "obs_time_1year","delta_1year",
                                                    "obs_time_2year","delta_2year")], 
                                 dt_sj_imp10_scale_final[,.SD,.SDcol=c(baseline_names,longitudinal_names)])
colnames(dt_sj_imp10_scale_final) <- gsub("/|->| |,|-|)|\\[|\\]", "_", colnames(dt_sj_imp10_scale_final))
colnames(dt_sj_imp10_scale_final)[15] <- "cat_disease_status_at_trans_Active_Non_Malignant"

############################################################################################
setwd("/Users/yzhou42/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/MLinHCT_YZ/ensemble_learner/results_JCO/SJ_results")
save(dt_sj_imp10_scale_final, file = "dt_sj_imp10_scale_final.Rdata")





############################################################################################
############################################################################################
############################################################################################
## process MSKCC data
## death info
dt_msk_imp10[,obs_time_100day:=ifelse(num_t_surv<100, num_t_surv, 100)]
dt_msk_imp10[,delta_100day:=ifelse(cat_death=="Deceased" & num_t_surv<100, "Deceased", "Alive")]
dt_msk_imp10[,delta_100day:=factor(delta_100day, levels = c("Alive","Deceased"))]

dt_msk_imp10[,obs_time_1year:=ifelse(num_t_surv<365, num_t_surv, 365)]
dt_msk_imp10[,delta_1year:=ifelse(cat_death=="Deceased" & num_t_surv<365, "Deceased", "Alive")]
dt_msk_imp10[,delta_1year:=factor(delta_1year, levels = c("Alive","Deceased"))]

dt_msk_imp10[,obs_time_2year:=ifelse(num_t_surv<365*2, num_t_surv, 365*2)]
dt_msk_imp10[,delta_2year:=ifelse(cat_death=="Deceased" & num_t_surv<365*2, "Deceased", "Alive")]
dt_msk_imp10[,delta_2year:=factor(delta_2year, levels = c("Alive","Deceased"))]

dt_msk_imp10[,cat_dt_trans:=ifelse(as.numeric(substr(num_dt_trans,1,4))>2013,1,0)]
dt_msk_imp10[,num_dt_trans:=NULL]

table(dt_msk_imp10$delta_100day[seq(1,nrow(dt_msk_imp10),m)])  ## 8 death
table(dt_msk_imp10$delta_1year[seq(1,nrow(dt_msk_imp10),m)])  ## 25 death
table(dt_msk_imp10$delta_2year[seq(1,nrow(dt_msk_imp10),m)])  ## 43 death

#############################################
## recode aGVHD
dt_msk_imp10[,("max_cat_agvhd_2"):=lapply(.SD, function(i) agvhd_recode(i)), .SDcols="max_cat_agvhd_2"]

############################################################################################
## feature encoding
dt_msk_imp10_agg <- dt_msk_imp10 ##if no averaging is performed for multiple imputed datasets

feature_names <- setdiff(colnames(dt_msk_imp10_agg), c("entity_id",".id",".imp","cat_death","num_t_surv",
                                                       "obs_time_100day","delta_100day",
                                                       "obs_time_1year","delta_1year",
                                                       "obs_time_2year","delta_2year"))
cat_names <- stringr::str_subset(feature_names, "cat_")
lgl_names <- stringr::str_subset(feature_names, "lgl")
num_names <- setdiff(feature_names, c(cat_names, lgl_names))

dt_msk_imp10_dum_cat <- do.call(cbind, lapply(1:length(cat_names), function(i) fastDummies::dummy_cols(dt_msk_imp10_agg[,.SD, .SDcol=cat_names[i]])[,-1]))
dt_msk_imp10_dum_cat <- dt_msk_imp10_dum_cat[,.SD,.SDcols = which(sapply(1:ncol(dt_msk_imp10_dum_cat),function(i) length(which(dt_msk_imp10_dum_cat[,.SD,.SDcol=i]==1)))>0)]


# num_means_sj <- dt_sj_imp10_agg[, sapply(.SD, function(i) mean(i)), .SDcol=num_names]
# num_sds_sj <- dt_sj_imp10_agg[, sapply(.SD, function(i) sd(i)), .SDcol=num_names]
# num_means_msk <- dt_msk_imp10_agg[, sapply(.SD, function(i) mean(i)), .SDcol=num_names]
# num_sds_msk <- dt_msk_imp10_agg[, sapply(.SD, function(i) sd(i)), .SDcol=num_names]
# dt_msk_imp10_scale_num <- do.call(cbind, lapply(1:length(num_names), function(i) (dt_msk_imp10_agg[,.SD,.SDcol=num_names[i]]-num_means[i])/num_sds[i]))
dt_msk_imp10_scale_num <- dt_msk_imp10_agg[, lapply(.SD, function(i) scale(i, scale = TRUE, center = TRUE)), .SDcol=num_names]
colnames(dt_msk_imp10_scale_num) <- num_names

dt_msk_imp10_agg[, (lgl_names):=lapply(.SD, as.logical), .SDcols=lgl_names]

dt_msk_imp10_scale <- cbind(dt_msk_imp10_dum_cat, dt_msk_imp10_scale_num, dt_msk_imp10_agg[,.SD,.SDcol=lgl_names])

############################################################################################
dt_msk_imp10_scale_final <- dt_msk_imp10_scale[,.SD,.SDcols=colnames(dt_msk_imp10_scale)[!(colnames(dt_msk_imp10_scale)%in%remove_names)]]

############################################################################################
baseline_names <- setdiff(stringr::str_subset(colnames(dt_msk_imp10_scale_final), "cat_|lgl_|age|trans|sex"),
                          c("max_cat_agvhd_2_0","max_cat_agvhd_2_1"))
longitudinal_names <- setdiff(colnames(dt_msk_imp10_scale_final), baseline_names)


dt_msk_imp10_scale_final <- cbind(dt_msk_imp10_agg[,c(".id","cat_death","num_t_surv",
                                                      "obs_time_100day","delta_100day",
                                                      "obs_time_1year","delta_1year",
                                                      "obs_time_2year","delta_2year")], 
                                  dt_msk_imp10_scale_final[,.SD,.SDcol=c(baseline_names,longitudinal_names)])
colnames(dt_msk_imp10_scale_final) <- gsub("/|->| |,|-|)|\\[|\\]", "_", colnames(dt_msk_imp10_scale_final))
colnames(dt_msk_imp10_scale_final)[15] <- "cat_disease_status_at_trans_Active_Non_Malignant"

############################################################################################
setwd("/Users/yzhou42/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/MLinHCT_YZ/ensemble_learner/results_JCO/MSKCC_results")
save(dt_msk_imp10_scale_final, file = "dt_msk_imp10_scale_final.Rdata")






