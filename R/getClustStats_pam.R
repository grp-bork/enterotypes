# following protocol from 2011 paper, as described in 2017 paper 
# here: https://enterotype.embl.de/enterotypes.html

library(here)
library(futile.logger)
require(cluster)
require(clusterSim)
library(tidyverse)

flog.threshold(DEBUG)


flog.info("Loading data")

outDir <- here("results/")
figDir <- here("figures/")
#Load data on motus level 
dataMatrix <- readRDS(here("data/motusMatrix_Genuslevel_gtdb_adults_filtered.rds"))
dataMatrix <- dataMatrix[rowSums(dataMatrix[, -1])>0,]
distObj <- readRDS(here("data/motusMatrix_Genuslevel_dist_JSD.rds"))
distMethod <- "distJSD"
cohort <- "adults_filtered"

# filter the dataMatrix for samples
dataMatrix <- dataMatrix[,colnames(dataMatrix) %in% rownames(as.matrix(distObj))]

# remove NA rows
dataMatrix <- dataMatrix[rowSums(is.na(dataMatrix)) == 0, ]


pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

ch.index=NULL
silhouette.width=NULL
data.dist <- as.matrix(distObj)
# filter the distance matrix for samples
data.dist <- data.dist[rownames(data.dist) %in% colnames(dataMatrix), colnames(data.dist) %in% colnames(dataMatrix)]

for (k in 1:20){
  if (k==1) {
    ch.index[k]=NA
    silhouette.width[k]=NA
  } else {
    flog.info("Calculating PAM clusters for k = %s", k)
    data.cluster_temp=pam.clustering(data.dist, k)
    flog.info("Calculating CH index for k = %s", k)
    ch.index[k]=clusterSim::index.G1(t(dataMatrix),
                                     data.cluster_temp,
                                     d = data.dist,
                                     centrotypes = "medoids")
    flog.info("Calculating Silhouette width for k = %s", k)
    silhouette.width[k]=mean(silhouette(data.cluster_temp,
                                        data.dist)[,3])
  }
}

flog.info("Saving results")
saveRDS(file = paste0(outDir, "/chIndexVector_", Sys.Date(),distMethod, "_",cohort,  ".Rds"),object = ch.index)
saveRDS(file = paste0(outDir, "/silhouetteWidthVector_", Sys.Date(),distMethod, "_",cohort,  ".Rds"),object = silhouette.width)

ch.index <- readRDS(file = paste0(outDir, "/chIndexVector_", Sys.Date(),distMethod,"_", cohort,  ".Rds"))
silhouette.width <- readRDS(file = paste0(outDir, "/silhouetteWidthVector_",Sys.Date(), distMethod,"_", cohort, ".Rds"))

flog.info("Plotting results")

png(paste0(figDir,"/chIndexPlot_", Sys.Date() ,distMethod, "_",cohort,  ".png"))
plot(ch.index, type="h", xlab="k clusters", ylab="CH index")
dev.off()

png(paste0(figDir, "/silhouetteWidthPlot_", Sys.Date(),distMethod,  "_",cohort, ".png"))
plot(silhouette.width, type="h", xlab="k clusters", ylab="Silhouette width")
dev.off()


 
