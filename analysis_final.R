rm(list=ls())

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if(length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

packages <- c("data.table","e1071","randomForest","naivebayes","corrplot","adabag","caret","caretEnsemble","glmnet",
              "RSNNS","mlbench","keras","tensorflow","DMwR","ROSE","pROC","GA","qs","dplyr","dm","mlr3","mice",
              "fastAdaboost","stringr","fastshap","fastDummies","tidyverse","pROC","corrplot","caTools","ipred",
              "plyr","rpart","earth","mda","nnet","gbm","kernlab","ggplot2","survival","stepPlr","randomGLM",
              "mboost","ada","deepboost","monmlp","cowplot","survminer","shapr","fastshap") #"RWeka","rJava", 
ipak(packages)

#fastAdaboost
#DMwR
#deepboost
#adabag

############################################################################################
############################################################################################
## function to perform univariate logistic regression (dimension reduction)
uni_logit <- function(data, outcome, i){
  
  formula_uni <- as.formula(paste0(outcome,"~", colnames(data)[c(baseline_index,longitudinal_index)[i]]))
  model_logit <- glm(formula_uni, data=data, family = "binomial")
  pvalue <- coef(summary(model_logit))[2,4]  
  
  return(pvalue)
}

## function to force all categories of a selected categorical variable into the final model
force_cat <- function(var_screen, cat_names, data_pvalue){
  
  var_cat <- var_screen[which(var_screen$name %in% cat_names),] 
  var_noncat <- var_screen[which(!(var_screen$name %in% cat_names)),] 
  
  var_cat_all <- NULL
  
  if(nrow(var_cat)>0){
    if(("cat_degree_match10_5" %in% var_cat$name) | ("cat_degree_match10_6" %in% var_cat$name) |
       ("cat_degree_match10_7" %in% var_cat$name) | ("cat_degree_match10_9" %in% var_cat$name) |
       ("cat_degree_match10_10" %in% var_cat$name)){
      var_cat_all <- rbind(var_cat_all, data_pvalue[which(data_pvalue$name %in% c("cat_degree_match10_5","cat_degree_match10_6",
                                                                                  "cat_degree_match10_7","cat_degree_match10_9",
                                                                                  "cat_degree_match10_10")),])
    }
    
    if(("cat_disease_status_at_trans_Active_Non_Malignant" %in% var_cat$name) | ("cat_disease_status_at_trans_Remission" %in% var_cat$name)){
      var_cat_all <- rbind(var_cat_all, data_pvalue[which(data_pvalue$name %in% c("cat_disease_status_at_trans_Active_Non_Malignant","cat_disease_status_at_trans_Remission")),])
    }
    
    if("cat_donor_sex_Female" %in% var_cat$name){
      var_cat_all <- rbind(var_cat_all, data_pvalue[data_pvalue$name=="cat_donor_sex_Female",])
    }
    
    if(("cat_dx_grp_ALL" %in% var_cat$name) | ("cat_dx_grp_AML" %in% var_cat$name) | ("cat_dx_grp_Aplastic_Anemia" %in% var_cat$name) |
       ("cat_dx_grp_Sickle_Cell_Anemia" %in% var_cat$name) | ("cat_dx_grp_Lymphoma" %in% var_cat$name) | ("cat_dx_grp_Myelodysplastic_Syndrome" %in% var_cat$name) |
       ("cat_dx_grp_CML" %in% var_cat$name) | ("cat_dx_grp_Other_Non_Malignancy" %in% var_cat$name)){
      var_cat_all <- rbind(var_cat_all, data_pvalue[which(data_pvalue$name %in% c("cat_dx_grp_ALL", "cat_dx_grp_AML","cat_dx_grp_Aplastic_Anemia",
                                                                                  "cat_dx_grp_Sickle_Cell_Anemia","cat_dx_grp_Lymphoma","cat_dx_grp_Myelodysplastic_Syndrome",
                                                                                  "cat_dx_grp_CML","cat_dx_grp_Other_Non_Malignancy")),])
    }
    
    if("cat_ethnicity_Hispanic_or_Latino" %in% var_cat$name){
      var_cat_all <- rbind(var_cat_all, data_pvalue[data_pvalue$name=="cat_ethnicity_Hispanic_or_Latino",])
    }
    
    if(("cat_matched_related_matched_related" %in% var_cat$name) | ("cat_matched_related_mismatched_related" %in% var_cat$name) | 
       ("cat_matched_related_matched_unrelated" %in% var_cat$name)){
      var_cat_all <- rbind(var_cat_all, data_pvalue[which(data_pvalue$name %in% c("cat_matched_related_matched_related", 
                                                                                  "cat_matched_related_mismatched_related",
                                                                                  "cat_matched_related_matched_unrelated")),])
    }
    
    if(("cat_prep_type_Myeloablative" %in% var_cat$name) | ("cat_prep_type_Reduced_Intensity" %in% var_cat$name)){
      var_cat_all <- rbind(var_cat_all, data_pvalue[which(data_pvalue$name %in% c("cat_prep_type_Myeloablative", "cat_prep_type_Reduced_Intensity")),])
    }
    
    if("cat_product_type_Marrow" %in% var_cat$name){
      var_cat_all <- rbind(var_cat_all, data_pvalue[data_pvalue$name=="cat_product_type_Marrow",])
    }
    
    if(("cat_race_Black" %in% var_cat$name) | ("cat_race_White" %in% var_cat$name)){
      var_cat_all <- rbind(var_cat_all, data_pvalue[which(data_pvalue$name %in% c("cat_race_Black", "cat_race_White")),])
    }
    
    if("cat_sex_Female" %in% var_cat$name){
      var_cat_all <- rbind(var_cat_all, data_pvalue[data_pvalue$name=="cat_sex_Female",])
    }
    
    if(("cat_sex_dir_Female_Male" %in% var_cat$name) | ("cat_sex_dir_Male_Female" %in% var_cat$name) |
       ("cat_sex_dir_Male_Male" %in% var_cat$name)){
      var_cat_all <- rbind(var_cat_all, data_pvalue[which(data_pvalue$name %in% c("cat_sex_dir_Female_Male", "cat_sex_dir_Male_Female",
                                                                                  "cat_sex_dir_Male_Male")),])
    }
    
    if("cat_dt_trans_1" %in% var_cat$name){
      var_cat_all <- rbind(var_cat_all, data_pvalue[data_pvalue$name=="cat_dt_trans_1",])
    }
    
    if("lgl_donor_related" %in% var_cat$name){
      var_cat_all <- rbind(var_cat_all, data_pvalue[data_pvalue$name=="lgl_donor_related",])
    }
    
    if("lgl_malignant" %in% var_cat$name){
      var_cat_all <- rbind(var_cat_all, data_pvalue[data_pvalue$name=="lgl_malignant",])
    }
    
    if(("max_cat_agvhd_2_0" %in% var_cat$name) | ("max_cat_agvhd_2_1" %in% var_cat$name)){
      var_cat_all <- rbind(var_cat_all, data_pvalue[which(data_pvalue$name %in% c("max_cat_agvhd_2_0", "max_cat_agvhd_2_1")),])
    }
    
  }
  
  var_screen_final <- rbind(var_cat_all, var_noncat)
  
  return(var_screen_final)
}

## function to perform bootstrap on testing data
boot_pred <- function(data, outcome, i, models, stack.rf=NULL){
  set.seed(1234+i)
  
  data_boot <- data[sample(1:nrow(data), size=nrow(data), replace = TRUE),]
  
  if(is.null(stack.rf)){
    model_preds_boot <- predict(models, newdata=data_boot, type="prob")[,2]
    model_preds_boot <- data.frame(model_preds_boot)
  }else{
    model_preds_boot <- lapply(models, predict, newdata=data_boot, type="prob")
    model_preds_boot <- lapply(model_preds_boot, function(x) x[,2])
    model_preds_boot <- data.frame(model_preds_boot)
    ens_preds_boot <- predict(stack.rf, newdata=data_boot, type="prob")
    model_preds_boot$ensemble <- ens_preds_boot  
  }
  
  return(caTools::colAUC(model_preds_boot, unlist(data_boot[,.SD,.SDcols=outcome])))
}

# ## function to build ensemble learner
# ensemble_builder <- function(data_train, data_test, outcome, var_screen, 
#                              algorithmList=c("naive_bayes","rf"), #c("knn","glm","naive_bayes","fda","mlp","gbm","xgbTree")
#                              model, seed=20221010, nboot=500, kfold=5, repeats=10){
#   
#   ##############################################
#   ## determine predictors
#   if(model=="full"){
#     var_names <- var_screen[,name]
#   }else if(model=="base_only"){
#     var_names <- var_screen[index%in%baseline_index,name]
#   }else if(model=="base_origin"){
#     var_names <- colnames(data_train)[baseline_index]
#   }else if(model=="long_only"){
#     var_names <- var_screen[index%in%longitudinal_index,name]
#   }
#   
#   ##############################################
#   control <- caret::trainControl(method="repeatedcv", 
#                                  number=kfold, 
#                                  repeats=repeats,
#                                  savePredictions=TRUE, 
#                                  classProbs=TRUE,
#                                  summaryFunction=twoClassSummary)
#   
#   
#   set.seed(seed)
#   formula_model <- as.formula(paste0(outcome, "~."))
#   models_fit <- caretEnsemble::caretList(formula_model, 
#                                      data=data_train[,.SD,.SDcols=c(outcome,var_names)],
#                                      trControl=control,
#                                      methodList=algorithmList,
#                                      metric ="ROC")
#   
#   # tuneList=list(
#   #   knn=caretModelSpec(method="knn", tuneGrid=data.frame(k=seq(2,10,1))),
#   #   glm=caretModelSpec(method="glm"),
#   #   naive_bayes=caretModelSpec(method="naive_bayes", tuneGrid=expand.grid(usekernel=c(TRUE,FALSE),
#   #                                                                         laplace=c(0,0.5,1), 
#   #                                                                         adjust=c(0.75,1,1.25,1.5))),
#   #   fda=caretModelSpec(method="fda", tuneGrid=expand.grid(degree=1, nprune=(1:10)*3)),
#   #   mlp=caretModelSpec(method="mlp", tuneGrid=expand.grid(size=seq(1:5))),
#   #   xgbTree=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(nrounds=(1:4)*50, 
#   #                                                                 max_depth=2:5, 
#   #                                                                 eta=c(0.3,0.4), 
#   #                                                                 gamma=0,
#   #                                                                 colsample_bytree=c(0.6,0.8,1.0), 
#   #                                                                 min_child_weight=2,
#   #                                                                 subsample=c(0.5,0.75,1.0))))
#   
#   
#   
#   results <- caret::resamples(models_fit)
#   summary(results)
#   
#   #write.table(results$values, file=paste0("ROC_sen_spec_values_",model,".txt"), col.names = TRUE, row.names = FALSE, sep = ",")
#   # png(paste0("ROC_sec_spec_",model,".png"), height = 10, width = 10, res = 300, units = "in")
#   # dotplot(results)
#   # dev.off()
# 
#   #write.table(modelCor(results), file=paste0("ROC_cor_",model,".txt"), col.names = TRUE, row.names = FALSE, sep = ",")
#   # png(paste("ROC_cor_",model,".png"), height = 10, width = 10, res = 300, units = "in")
#   # splom(results)
#   # dev.off()
# 
#   stackControl <- caret::trainControl(method="repeatedcv", 
#                                       number=kfold, 
#                                       repeats=repeats, 
#                                       savePredictions=TRUE, 
#                                       classProbs=TRUE,
#                                       summaryFunction=twoClassSummary)
#   
#   stack.rf <- caretEnsemble::caretStack(models_fit, 
#                                         method="rf", 
#                                         trControl=stackControl,
#                                         metric="ROC")
#   print(stack.rf)
#   
#   ##############################################
#   model_preds_prob <- lapply(models_fit, predict, newdata=data_test, type="prob")
#   model_preds_prob <- lapply(model_preds_prob, function(x) x[,2])
#   model_preds_prob <- data.frame(model_preds_prob)
#   ens_preds_prob <- predict(stack.rf, newdata=data_test, type="prob")
#   model_preds_prob$ensemble <- ens_preds_prob
#   
#   ##############################################
#   ## bootstrap on testing data
#   pred_boot_auc <- t(sapply(1:nboot, function(i) boot_pred(data_test, outcome, i, models_fit, stack.rf)))
#   colnames(pred_boot_auc) <- c(algorithmList, "ensemble")
#   
#   # pred_boot_auc_summary <- as.data.frame(t(do.call(cbind,lapply(1:ncol(pred_boot_auc), function(i) quantile(pred_boot_auc[,i], probs=c(0.025,0.5,0.975))))))
#   # pred_boot_auc_summary$method <- c(algorithmList, "ensemble")
#   
#   ##############################################
#   return(list("pred_train"=results$values,
#               "pred_test_prob"=model_preds_prob,
#               "pred_test_auc"=caTools::colAUC(model_preds_prob, unlist(data_test[,.SD,.SDcols=outcome])),
#               "pred_boot_auc"=pred_boot_auc#,
#               #"pred_boot_auc_summary"=pred_boot_auc_summary
#               ))
# }

## function to build base-learner (only one machine learning algorithm)
base_builder <- function(data_train, data_test, data_test_ex, outcome, var_screen_all, 
                         method, model, seed=20221010, nboot=500, kfold=5, repeats=10){
  
  ##############################################
  ## determine predictors
  if(model=="full"){
    var_names <- var_screen_all[,name]
  }else if(model=="base_only"){
    var_names <- var_screen_all[index%in%baseline_index,name]
  }else if(model=="base_origin"){
    var_names <- colnames(data_train)[baseline_index]
  }else if(model=="long_only"){
    var_names <- var_screen_all[index%in%longitudinal_index,name]
  }
  
  ##############################################
  control <- caret::trainControl(method="repeatedcv", 
                                 number=kfold, 
                                 repeats=repeats,
                                 savePredictions=TRUE, 
                                 classProbs=TRUE,
                                 summaryFunction=twoClassSummary)
  
  
  set.seed(seed)
  formula_model <- as.formula(paste0(outcome, "~."))
  model_fit <- caret::train(formula_model, 
                            data=data_train[,.SD,.SDcols=c(outcome,var_names)],
                            trControl=control,
                            method = method,
                            metric ="ROC")
  
  model_fit
  
  ##############################################
  model_preds_prob <- predict(model_fit, newdata=data_test, type="prob")[,2]
  model_preds_prob <- data.frame(model_preds_prob)
  
  model_preds_prob_ex <- predict(model_fit, newdata=data_test_ex, type="prob")[,2]
  model_preds_prob_ex <- data.frame(model_preds_prob_ex)
  
  ##############################################
  ## bootstrap on testing data
  pred_boot_auc <- sapply(1:nboot, function(i) boot_pred(data_test, outcome, i, model_fit))
  pred_boot_auc_ex <- sapply(1:nboot, function(i) boot_pred(data_test_ex, outcome, i, model_fit))
  
  ##############################################
  ## Shapley Additive Explanations
  pred_wrapper <- function(object, newdata) predict(object, newdata=newdata, type="prob")[,2]
  
  # cl = parallel::makeCluster(4L)
  # doParallel::registerDoParallel(cl, cores = length(cl))
  set.seed(seed)
  shap <- fastshap::explain(
    model_fit,
    X = data_train[,.SD,.SDcols=var_names],
    pred_wrapper = pred_wrapper,
    nsim = 10L,
    adjust = TRUE,
    newdata = data_test[,.SD,.SDcols=var_names]
  )
  # doParallel::stopImplicitCluster()
  shap_values <- shapviz::shapviz(shap, X = data_test[,.SD,.SDcols=var_names])
  
  ##############################################
  return(list("pred_train"=model_fit$results,
              "pred_test_prob"=model_preds_prob,
              "pred_test_prob_ex"=model_preds_prob_ex,
              "pred_test_auc"=caTools::colAUC(model_preds_prob, unlist(data_test[,.SD,.SDcols=outcome])),
              "pred_test_auc_ex"=caTools::colAUC(model_preds_prob_ex, unlist(data_test_ex[,.SD,.SDcols=outcome])),
              "pred_boot_auc"=pred_boot_auc,
              "pred_boot_auc_ex"=pred_boot_auc_ex,
              "shap_values"=shap_values))
}

#function for final prediction
MLinHCT_pred <- function(data_scale_final, analysis_index, outcome, ratio=0.7, seed,
                         baseline_index, longitudinal_index, n_var=50, method, cat_names, train_index, data_test_ex){
  
  data_analysis <- data_scale_final[analysis_index,]
  
  data_train <- data_analysis[train_index,]
  data_test <- data_analysis[-train_index,]
  
  # if(outcome=="delta_100day"){
  #   data_train <- SMOTE("delta_100day" ~ ., data  = data_train)
  # }
  
  ############################################################################################
  ## dimension reduction (logistic regression)
  pvalue_vec <- sapply(1:length(c(baseline_index, longitudinal_index)), function(i) uni_logit(data_train, outcome, i))
  
  data_pvalue <- data.table(index=c(baseline_index, longitudinal_index),
                            name=colnames(data_train)[c(baseline_index, longitudinal_index)],
                            pvalue=pvalue_vec)
  
  data_pvalue <- data_pvalue[order(pvalue)]
  
  ############################################################################################
  var_screen <- data_pvalue[1:n_var,]
  
  var_screen_all <- force_cat(var_screen, cat_names = cat_names, data_pvalue=data_pvalue)
  
  ############################################################################################
  full_results <- base_builder(data_train, data_test, data_test_ex, outcome, var_screen_all, method="naive_bayes", model="full", seed=seed)
  #base_only_results <- base_builder(data_train, data_test, data_test_ex, outcome, var_screen_all, method="naive_bayes", model="base_only", seed=seed)
  base_origin_results <- base_builder(data_train, data_test, data_test_ex, outcome, var_screen_all, method="naive_bayes", model="base_origin", seed=seed)
  long_only_results <- base_builder(data_train, data_test, data_test_ex, outcome, var_screen_all, method="naive_bayes", model="long_only", seed=seed)
  
  return(list("data_train"=data_train,
              "data_test"=data_test,
              "data_pvalue"=data_pvalue,
              "var_screen_all"=var_screen_all,
              "full_results"=full_results,
              #"base_only_results"=base_only_results,
              "base_origin_results"=base_origin_results,
              "long_only_results"=long_only_results))
}


############################################################################################
############################################################################################
############################### different model construction ###############################
############################################################################################
load("/Users/yzhou42/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/MLinHCT_YZ/ensemble_learner/results_JCO/SJ_results/dt_sj_imp10_scale_final.Rdata")
load("/Users/yzhou42/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/MLinHCT_YZ/ensemble_learner/results_JCO/MSKCC_results/dt_msk_imp10_scale_final.Rdata")

n_sj <- length(unique(dt_sj_imp10_scale_final[,.id]))
n_msk <- length(unique(dt_msk_imp10_scale_final[,.id]))
m <- max(table(dt_sj_imp10_scale_final[,.id]))

baseline_names <- setdiff(stringr::str_subset(colnames(dt_sj_imp10_scale_final), "cat_|lgl_|age|trans|sex"),
                          c("max_cat_agvhd_2_0","max_cat_agvhd_2_1","cat_death"))
longitudinal_names <- setdiff(colnames(dt_sj_imp10_scale_final), c(baseline_names, 
                                                                   ".id","cat_death","num_t_surv",
                                                                   "obs_time_100day","delta_100day",
                                                                   "obs_time_1year","delta_1year",
                                                                   "obs_time_2year","delta_2year"))
baseline_index <- which(colnames(dt_sj_imp10_scale_final) %in% baseline_names)  
longitudinal_index <- which(colnames(dt_sj_imp10_scale_final) %in% longitudinal_names) 


feature_names <- setdiff(colnames(dt_sj_imp10_scale_final), c(".id","cat_death","num_t_surv",
                                                              "obs_time_100day","delta_100day",
                                                              "obs_time_1year","delta_1year",
                                                              "obs_time_2year","delta_2year"))
cat_names <- stringr::str_subset(feature_names, "cat_|lgl_")

############################################################################################
############################################################################################
## overall survival
set.seed(20221010)

data_analysis <- dt_sj_imp10_scale_final[seq(1,nrow(dt_sj_imp10_scale_final),by=m),]
train_index <- caret::createDataPartition(unlist(data_analysis[,.SD,.SDcols="delta_100day"]), p = 0.7, list = FALSE, times = 1)

results_100day <- lapply(1:m, function(rep) MLinHCT_pred(data_scale_final=dt_sj_imp10_scale_final,
                                                         analysis_index=seq(rep,nrow(dt_sj_imp10_scale_final),by=m),
                                                         outcome="delta_100day",
                                                         ratio=0.7,
                                                         seed=20221010,
                                                         baseline_index,
                                                         longitudinal_index,
                                                         n_var=50,
                                                         method="naive_bayes",
                                                         cat_names=cat_names,
                                                         train_index=train_index,
                                                         data_test_ex = dt_msk_imp10_scale_final))

results_1year <- lapply(1:m, function(rep) MLinHCT_pred(data_scale_final=dt_sj_imp10_scale_final,
                                                        analysis_index=seq(rep,nrow(dt_sj_imp10_scale_final),by=m),
                                                        outcome="delta_1year",
                                                        ratio=0.7,
                                                        seed=20221010,
                                                        baseline_index,
                                                        longitudinal_index,
                                                        n_var=50,
                                                        method="naive_bayes",
                                                        cat_names=cat_names,
                                                        train_index=train_index,
                                                        data_test_ex = dt_msk_imp10_scale_final))

results_2year <- lapply(1:m, function(rep) MLinHCT_pred(data_scale_final=dt_sj_imp10_scale_final,
                                                        analysis_index=seq(rep,nrow(dt_sj_imp10_scale_final),by=m),
                                                        outcome="delta_2year",
                                                        ratio=0.7,
                                                        seed=20221010,
                                                        baseline_index,
                                                        longitudinal_index,
                                                        n_var=50,
                                                        method="naive_bayes",
                                                        cat_names=cat_names,
                                                        train_index=train_index,
                                                        data_test_ex = dt_msk_imp10_scale_final))

############################################################################################
############################################################################################
setwd("/Users/yzhou42/Library/CloudStorage/OneDrive-St.JudeChildren'sResearchHospital/MLinHCT_YZ/ensemble_learner/results_JCO/SJ_results")
save(results_100day, file = "results_100day.Rdata")
save(results_1year, file = "results_1year.Rdata")
save(results_2year, file = "results_2year.Rdata")

