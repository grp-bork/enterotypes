set.seed(35215)
library(tidyverse)
library(here)

library(futile.logger)
require(cluster)
require(clusterSim)

flog.info("load data")
# load data ------------------------------
dataMatrix <- readRDS(here("data/motusMatrix_Genuslevel_gtdb_adults_filtered.rds"))
# remove the -1 'species'
dataMatrix <- dataMatrix[row.names(dataMatrix)!="-1",]
dataMatrix <- dataMatrix[rowSums(is.na(dataMatrix)) == 0, ]
distObj <- readRDS(here("data/motusMatrix_Genuslevel_dist_JSD.rds"))
distMatrix <- as.matrix(distObj)
distMatrix <- distMatrix[rownames(distMatrix) %in% colnames(dataMatrix), colnames(distMatrix) %in% colnames(dataMatrix)]
metadata <- readRDS(here("data/metadata.rds"))  %>% filter(!is.na(studies))
metadata <- metadata %>% filter(internal_id %in% colnames(distMatrix))


samplesStudiesDf <- metadata %>% dplyr::select(studies, internal_id) %>% magrittr::set_colnames(c("studyID", "sampleID"))
studyIDs <- samplesStudiesDf %>% pull(studyID) %>% unique()
sampleIDs <- samplesStudiesDf %>% pull(sampleID) %>% unique()


flog.info("make rarefraction iterations")
# make rarefaction iterations --------------------------
rarefactionLevels <- seq(from=0.3,to=1,by=0.1) # or could go up to 1 if clustering is not deterministic
nReps <- 20 # number of repetitions per rarefaction level


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


# do clustering on rarefied data ----------------------

doSelectedClustering <- function(selections){
  selectedSamples <- selections$selectedSamples
  selectedStudies <- selections$selectedStudies
  propStudies <- selections$propStudies
  print(paste0("Clustering ",propStudies," prop of studies: ",length(selectedStudies)," studies, ",length(selectedSamples)," samples..."))
  # create rarefied distance matrix
  distMatrixSelected <- distMatrix[selectedSamples,selectedSamples]
  
  # filter the dataMatrix for samples
  dataMatrixSelected <- dataMatrix[,rownames(as.matrix(distMatrixSelected))]
  
  # do clustering
  res <- getBestClusterNumber(distMatrixSelected, dataMatrixSelected)
  res$propStudies <- propStudies
  res$nStudies <- length(selectedStudies)
  res$nSamples <- length(selectedSamples)
  #best_k <- res$k
  return(res)
}

# pam clustering function 
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

getBestClusterNumber <- function(distMatrixSelected,dataMatrixSelected ){
  #Sys.sleep(1)# for testing parallelisation
  # do clustering for each value of k to be tested
  # gives list of cluster results, ready to calculate CH index from
  
    # calculate cluster performance metrics for each clustering attempt
  ch.index=NULL
  silhouette.width=NULL
  data.dist <- as.matrix(distMatrixSelected)
  for (k in 2:10){
    if (k==1) {
      ch.index[k]=NA
      silhouette.width[k]=NA
    } else {
      flog.info("Calculating PAM clusters for k = %s", k)
      data.cluster_temp=pam.clustering(data.dist, k)
      flog.info("Calculating CH index for k = %s", k)
      ch.index[k]=clusterSim::index.G1(t(dataMatrixSelected),
                                       data.cluster_temp,
                                       d = data.dist,
                                       centrotypes = "medoids")
      flog.info("Calculating Silhouette width for k = %s", k)
      silhouette.width[k]=mean(silhouette(data.cluster_temp,
                                          data.dist)[,3])
    }
  }
  
  ############################################################
  # pick best cluster number based on silhouette
  ############################################################
  bestK_sil_index <-  which.max(silhouette.width) # since vector starts at k = 2
  
  # % distance between best and second best
  secondBestK_sil <-  sort(silhouette.width)[length(silhouette.width) -2] # -2 because ot the NA
  diffProp_sil <- (silhouette.width[bestK_sil_index] - secondBestK_sil) / silhouette.width[bestK_sil_index]
  
  bestK_sil <- bestK_sil_index  
  bestK_sil_value <- silhouette.width[bestK_sil_index] # in case you need it
  bestK_sil_value2nd <- secondBestK_sil # in case you need it
  
  ############################################################
  # pick best cluster number based on ch index 
  ############################################################
  bestK_ch_index <- which.max(ch.index)
  
  # % distance between best and second best
  secondBestK_ch <- sort(ch.index)[length(ch.index) -2]
  diffProb_ch <- (ch.index[bestK_ch_index] - secondBestK_ch) / ch.index[bestK_ch_index]
  
  bestK_ch <- bestK_ch_index 
  bestK_ch_value <- ch.index[bestK_ch_index]
  bestK_ch_value2nd <- secondBestK_ch 
  
  ############################################################
  # get all values within 10% of best cluster number 
  ############################################################
  #threshold = 0.1
  #near_max_ch <- ch.index[ch.index >= bestK_ch_value * (1 - threshold) & ch.index <= bestK_ch_value * (1 + threshold)]
  #near_max_sil <- silhouette.width[silhouette.width >= bestK_sil_value * (1 - threshold) & silhouette.width <= bestK_sil_value * (1 + threshold)]
  
  return(list(k_best_sil=bestK_sil,
              k_sil_value=bestK_sil_value,
              k_sil_value_2ndBest=bestK_sil_value2nd,
              propDiffBetweenBestAndSecondBest_sil=diffProp_sil,
              k_best_ch=bestK_ch, 
              k_ch_value=bestK_ch_value, 
              k_ch_value_2ndBest = bestK_ch_value2nd, 
              probDiffBetweenBestAndSecondBest_ch=diffProb_ch))
}

# run the clustering analysis --------------------------------------
flog.info("run clustering")
# serial execution
#res <- map_dfr(iterations, doSelectedClustering)

# parallel execution
options(future.globals.maxSize = 3 * 1024^3)
future::plan(future::multisession, workers = 6) # change 6 to your number of cores/processes
resL <- furrr::future_map(.x = iterations, 
                          .f = doSelectedClustering,
                          .progress = T,
                          .options = furrr::furrr_options(seed = TRUE)) #map_dfr not in furrr package yet
res <- map_dfr(resL,bind_rows)

flog.info("save the results")
# save the results
saveRDS(res, here("results", Sys.Date(),"_pam_rarefraction_results.rds"))


  
