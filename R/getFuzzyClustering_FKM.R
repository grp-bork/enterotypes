# FKM FUZZY Clustering with gtdb  
# author: Marisa Metzger
# ######################

# we want to use the fuzzy clustering approach for soft clustering 

library(here)
library(fclust) # for FKM clustering
library(tidyverse)
library(futile.logger)


flog.info("Adults Filtered COHORT")
flog.info("Loading data")

outDir <- here("results/")
figDir <- here("figures/")
dataDir <- here("data/")

# load data
data_used <- readRDS(here("data/motusMatrix_Genuslevel_gtdb_adults_filtered.rds"))
data_used <- na.omit(data_used)

# remove rows with only zeros.
data_used <- t(data_used[rowSums(data_used[, -1])>0,])
dataset<- "motus2.5_adults_filtered"

# FKM from the clust package
library(fclust)
#
flog.info("fuzzy clustering")
# fuzzy k-means algorithm
data_fkm <- fclust::FKM(data_used, k=1:15, index = "SIL")
# get the best clutsernumber
clusterNumber <- as.data.frame(data_fkm$clus) %>%
  select(Cluster) %>%
  unique(.) %>%
  nrow(.)

saveRDS(data_fkm, paste0(outDir, "FKM_cluster_",Sys.Date(),"_", dataset, ".Rds"))

print(summary(data_fkm))




