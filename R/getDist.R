library(futile.logger)
flog.threshold(DEBUG)
library(here)
setwd(here("data/"))

# Load motus Matrix
flog.debug("Loading motusMatrix")
prefix <- "motusMatrix_Genuslevel"
dM <- readRDS(here("data/motusMatrix_Genuslevel_gtdb_adults_filtered.rds"))
# remove the -1 column 
dM <- dM[!rownames(dM) %in% c("-1"),]
source("../../../../scripts/distJSD_fast.R")

dM <- as.matrix(dM)


distMethods <- c("JSD")

# function to get distance 
distJSD_fmarotta <- function(inMatrix, pseudocount=0.000001, ...) {
  inMatrix[inMatrix == 0] <- pseudocount
  KLD_matrix <- sapply(as.data.frame(inMatrix), function(x) {
    crossprod(x, log(2 * x / (x + inMatrix)))
  })
  resultsMatrix <- sqrt(.5 * (KLD_matrix + t(KLD_matrix)))
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix)
}

# get JSD distance matrix

for(distMeth in distMethods){
  
  flog.debug("Calculating distance matrix with %s distances",distMeth)
  
  if(distMeth == "JSD"){
    distObj <- distJSD_fmarotta(dM)
   
  }else{
    # samples need to be rows; species are columns
    distObj <- dist(t(dM),method = tolower(distMeth))
  }
  
  flog.debug("Saving distance matrix")
  saveRDS(distObj,file = paste0(prefix,"_dist_",distMeth,".rds"))
}
