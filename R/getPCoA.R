# following protocol from 2011 paper, as described in 2017 paper 
# here: https://enterotype.embl.de/enterotypes.html

library(here)
library(futile.logger)
require(ade4)
require(cluster)
library(scales) 



## here costumize the cluster number --------------------------------------
numClusters <- 3
rep <- 1
# -------------------------------------------------------------------------
al<-0.7
 defaultCols <- c(alpha("deeppink3", al),
                  alpha("dodgerblue", al),
                  alpha("springgreen4", al),
                  alpha("chocolate", al),
                  alpha("orchid3", al),
                  alpha("gray50", al),
		alpha("olivedrab2", al))


flog.threshold(DEBUG)

outdir <- here("results/")
figdir <- here("figures/")
flog.info("Loading data")
# data
dataMatrix <- readRDS(here("data/motusMatrix_Genuslevel_gtdb_adults_filtered.rds"))
dataMatrix <- dataMatrix[!rownames(dataMatrix) %in% c("-1"),]

distObj <- readRDS(here("data/motusMatrix_Genuslevel_dist_JSD.rds"))
#---------------------------------------------------------------------------------------------------
# filter the dataMatrix for samples
dataMatrix <- dataMatrix[,colnames(dataMatrix) %in% colnames(as.matrix(distObj))]

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

# remove taxa if their mean abundance across all samples was below 0.01%
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

data.denoized=noise.removal(dataMatrix, percent=0.01)
# remove samples with less than 80% of original community represented
#data.denoized<-data.denoized[,colSums(data.denoized)>0.8]

samplesToUse <- colnames(data.denoized)
data.dist <- as.matrix(distObj)[samplesToUse,samplesToUse]

# # PCoA ----------
flog.info("PCoA")
obs.pcoa=dudi.pco(as.dist(data.dist), scannf=F, nf=10,full = F)
saveRDS(obs.pcoa,paste0(outdir,"pcoa.rds"))

png(paste0(figdir,"pcoa_scree_.png"))
screeplot(obs.pcoa,npcs = 10)
dev.off()

#s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F)
png(paste0(figdir,"pcoa_k.png"))
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F,
        cell=0, cstar=0,
        col = defaultCols)

dev.off()

