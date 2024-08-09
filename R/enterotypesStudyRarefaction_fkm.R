set.seed(35215)
library(tidyverse)
library(here)
library(futile.logger)
library(fclust)


flog.info("load data")
# load data ------------------------------
dataMatrix <- readRDS(here("data/motusMatrix_Genuslevel_gtdb_adults_filtered.rds"))
dataMatrix <- dataMatrix[rowSums(is.na(dataMatrix)) == 0, ]

metadata <- readRDS(here("data/metadata.rds"))  %>% filter(!is.na(studies))
metadata <- metadata %>% filter(internal_id %in% colnames(dataMatrix))


samplesStudiesDf <- metadata %>% dplyr::select(studies, internal_id) %>% magrittr::set_colnames(c("studyID", "sampleID"))
studyIDs <- samplesStudiesDf %>% pull(studyID) %>% unique()
sampleIDs <- samplesStudiesDf %>% pull(sampleID) %>% unique()


flog.info("make rarefraction iterations")
# make rarefaction iterations --------------------------
rarefactionLevels <- c(seq(from=0.3,to=0.9,by=0.1),1) # or could go up to 1 if clustering is not deterministic
nReps <-  20 # number of repetitions per rarefaction level


nStudiesTotal <- unique(studyIDs) %>% length()
getIDs <- function(propStudies){
  nStudiesForIter <- round(nStudiesTotal*propStudies)
  selectedStudies <- sample(x = studyIDs,size = nStudiesForIter)
  selectedSamples <- samplesStudiesDf %>% filter(studyID %in% selectedStudies) %>% pull(sampleID)
  nSamplesForIter <- length(selectedSamples)
  return(list("selectedStudies"=selectedStudies,"selectedSamples"=selectedSamples,"propStudies"=propStudies))
}

x <- rarefactionLevels %>% rep(times=nReps) %>% sort()

iterations <- map(x,getIDs)

res <- list()
# do clustering on rarefied data ----------------------

doSelectedClustering <- function(selections){
  selectedSamples <- selections$selectedSamples
  selectedStudies <- selections$selectedStudies
  propStudies <- selections$propStudies
  print(paste0("Clustering ",propStudies," prop of studies: ",length(selectedStudies)," studies, ",length(selectedSamples)," samples..."))
 
  # filter the dataMatrix for samples
  dataMatrixSelected <- dataMatrix[,selectedSamples]
  
  # do clustering
  SW <- getBestClusterNumber(dataMatrixSelected)
  
  res$propStudies <- propStudies
  res$nStudies <- length(selectedStudies)
  res$nSamples <- length(selectedSamples)
  res$silhouette.vector <- SW 
  #best_k <- res$k
  return(res)
}

getBestClusterNumber <- function(dataMatrixSelected ){
  #Sys.sleep(1)# for testing parallelisation
  # do clustering for each value of k to be tested
  # gives list of cluster results, ready to calculate CH index from
 
  ############################################################
  # do clustering
  #############################################################
  flog.info("run fuzzy clustering")
  data_used <- t(dataMatrixSelected[rowSums(dataMatrixSelected)>0,])
  data_fkm <- fclust::FKM(data_used, k=2:7, index = "SIL")
  #data_fkm.f <- fclust::FKM(data_used, k=2:7, index = "SIL.F")
  silhouette.width <- data_fkm$criterion %>% as.vector()
 
  return(silhouette.width)
 }

# run the clustering analysis --------------------------------------
flog.info("run clustering")
# serial execution
#res <- map_dfr(iterations, doSelectedClustering)

# parallel execution
#install.packages("furrr")
#install.packages("future")
options(future.globals.maxSize = 3 * 1024^3)
future::plan(future::multisession, workers = 6) # change 6 to your number of cores/processes
resL <- furrr::future_map(.x = iterations, 
                          .f = doSelectedClustering,
                          .progress = T,
                          .options = furrr::furrr_options(seed = TRUE)) #map_dfr not in furrr package yet
#res <- map_dfr(resL,bind_rows)

flog.info("save the results")
# save the results
saveRDS(resL, here("results", Sys.Date(),"_motus_fuzzy_rarefraction_results.rds"))

