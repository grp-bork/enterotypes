## R script to construct prediction models for enterotyping

## load library
library(tidyverse)
library(caret)
library(doMC)
library(pROC)
library(scales)
registerDoMC(cores = 4)

## read training data
d <- read_rds("motus3Matrix_GTDBGenuslevel_adults_filtered.rds") %>% t() %>% data.frame(check.names = F)
md <- read.delim("metadata_adults_filtered_withEnterotypes.tsv", header = T, row.names = 1)
md <- md[order(rownames(md)), ]
d <- d[order(rownames(d)), ]
identical(rownames(d), rownames(md))

keep <- !is.na(md$studies)
d <- d[keep, ]
md <- md[keep, ]

## exclude low abundant and  "unassigned" species
keep <- apply(d, 2, mean) > 1E-4
d <- d[, keep]
d <- d %>% select(-contains("-1"))
dim(d)

ord <- order(rownames(d))
d <- d[ord, ]
md <- md[ord, ]

## read validation data
d2 <- read_rds("motus3Matrix_GTDBGenuslevel_adults_test.rds") %>% t()
md2 <- read.delim("metadata_tt_withEnterotypes.tsv", header = T, row.names = 1)
d2 <- d2[order(rownames(d2)), ]
md2 <- md2[order(rownames(md2)), ]

keep1 <- rownames(d2) %in% rownames(md2)
keep2 <- rownames(md2) %in% rownames(d2)
d2 <- d2[keep1, ]
md2 <- md2[keep2, ]

md2$pam_k3 <- md2$pam_k3 %>% str_replace("Firmicutes_1", "Firmicutes") ## modify label

## set parameters
study.int <- md$studies %>% as.numeric()
study.int <- groupKFold(md$studies, 10)
train_control <- trainControl(method = "cv", savePredictions = T, index = study.int)
tune_grid <- expand.grid( nrounds=c(100, 200, 300), max_depth = c(3:5), eta = c(0.01), gamma = c(0.01), colsample_bytree = c(0.75), subsample = 1, min_child_weight = 1)

## function for fuzzy clustering
make_model_for_fuzzy <- function(enterotype, prob, target, type){
  print(target)
  
  ## construct model
  model <- caret::train(d, enterotype, trControl = train_control, method = "xgbTree", metric = "Rsquared", tuneGrid = tune_grid)
  out <- paste0("results/model/model.", type, "_", target, ".rds")
  saveRDS(model, file = out)
  
  ## evaluate model  
  pred <- predict(model, d2)
  
  ## plot AUC
  df <- data.frame(id = rownames(d2), prob = prob, pred = pred)
  p <- ggplot(df, aes(x = prob, y = pred)) +
    theme_bw() +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm") +
    xlab("Probability") +
    ylab("Predicted value") +
    ggtitle(target) +
    ggpubr::stat_cor()
  p
  
  out <- paste0("results/model_motus3/scatter.validation.", type, "_", target, ".pdf")
  ggsave(p, filename = out, width = 3, height = 3)
  
  out <- paste0("results/model_motus3/scatter.validation.", type, "_", target, ".rds")
  saveRDS(df, file = out)
}



## construct models (fuzzy 3) #####################################
## for Bact_Phoc
enterotype <- md$prob_Bact_Phoc
prob <- md2$prob_Bact_Phoc
make_model_for_fuzzy(enterotype, prob, "Bact_Phoc", "f3")

## for Firmicutes
enterotype <- md$prob_Firmicutes
prob <- md2$prob_Firmicutes
make_model_for_fuzzy(enterotype, prob, "Firmicutes", "f3")

## for Prevotella
enterotype <- md$prob_Prevotella
prob <- md2$prob_Prevotella
make_model_for_fuzzy(enterotype, prob, "Prevotella", "f3")


## construct models (fuzzy 2)
## for Bact_Phoc
enterotype <- md$fkm_k2_prob_Bact_Phoc
prob <- md2$fkm_k2_prob_Bact_Phoc
make_model_for_fuzzy(enterotype, prob, "Bact_Phoc", "f2")

## for Prevotella
enterotype <- md$fkm_k2_prob_Prevotella
prob <- md2$fkm_k2_prob_Prevotella
make_model_for_fuzzy(enterotype, prob, "Prevotella", "f2")
###################################################################



## function for pam clustering
make_model_for_pam <- function(enterotype, target){
  print(target)
  
  ## construct model
  model <- caret::train(d, enterotype, trControl = train_control, method = "glmnet", metric = "roc", tuneGrid = expand.grid(alpha = 1, lambda = 10 ^ (1:10 * -1))) ## LASSO
  out <- paste0("results/model/model.", target, ".rds")
  saveRDS(model, file = out)
  
  ## evaluate model  
  pred <- predict(model, d2, type = "prob")
  rownames(pred) <- rownames(d2)
  
  ## save result
  prob.max <- apply(pred, 1, which.max)
  res <- data.frame(enterotype = md2[, target], enterotype_pred = colnames(pred)[prob.max], pred) %>% rownames_to_column()
  colnames(res)[1] <- "sample_id"
  colnames(res)[4:ncol(res)] <- colnames(res)[4:ncol(res)] %>% paste0(., "_pred")
  res <- res %>% arrange(sample_id)
  
  out <- paste0("results/model_evaluation_motus3/external_validation.", str_replace(target, "k", ""), ".tsv")
  write.table(res, file = out, sep = "\t", col.names = NA)
  
  ## calc AUC for each enterotype
  df <- c()
  auc.vec <- c()
  title <- c()
  enterotype <- factor(enterotype)

  for(i in levels(enterotype)){
    print(i)  
    lab <- ifelse(md2[, target] == i, 1, 0)
    
    roc <- roc(lab, pred[, i])
    df.temp <- data.frame(fpr = 1 - roc$specificities, tpr = roc$sensitivities, enterotype = i)
    auc <- round(roc$auc, digit = 3) %>% format(digits = 3)
    
    df <- rbind(df, df.temp)
    temp <- paste0(i, ": ", auc, "\n")
    title <- paste0(title, temp)
  }
  
  ## plot AUC
  roc <- ggplot(df, aes(x = fpr, y = tpr, col = enterotype)) + 
    theme_light() +
    geom_path() +
    xlab("False positive rate") +
    ylab("True positive rate") +
    ggtitle(target) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "gray") +
    annotate(geom = "text", x = 0.8, y = 0.05, label = title, size = 2) +
    theme(plot.title = element_text(size = 9))
  roc
  
  out <- paste0("results/model_motus3/roc.validation.", target, ".pdf")
  ggsave(roc, filename = out, width = 4, height = 3)
  
  out <- paste0("results/model_motus3/roc.validation.", target, ".rds")
  saveRDS(df, file = out)
}


## construct models (pam 3) ######################################
enterotype <- md$pam_k3
target <- "pam_k3"
make_model_for_pam(enterotype, "pam_k3")

## construct models (pam 2)
enterotype <- md$pam_k2
target <- "pam_k2"
make_model_for_pam(enterotype, "pam_k2")

## construct models (pam 4)
enterotype <- md$pam_k4
target <- "pam_k4"
make_model_for_pam(enterotype, "pam_k4")
###################################################################


