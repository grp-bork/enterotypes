
# DMM CLUSTERING

#based on this tutorial: https://microbiome.github.io/tutorials/DMM.html 

# The data for this analysis contains a sequence count motus Matrix (not relative abundance) on genus level (GDTB Taxonomy) for all adults (age > 5 years)

# SET UP ----- 

library(DirichletMultinomial)
library(magrittr)
library(dplyr)
library(here)

outDir <- paste0(here(), "results/")
dataDir <- paste0(here(), "data/")

# LOAD THE DATA ---- 

data_adults <- readRDS(paste0(dataDir, "motusMatrix_Genuslevel_gtdb_adults_filtered.rds"))
data_adults <- data_adults[row.names(data_adults)!="-1",]
data_adults <- na.omit(data_adults)
# load metadata 
md <- readr::read_tsv(paste0(dataDir, "metadata.rds"))


data_use <- data_adults
# remove rows with only zero 
data_use <- data_adults[rowSums(data_use[, -1])>0,]

# FOR RELATIVE DATA: TRANSFORMATION 
da2 <- round(data_use * 10000, 0)

# square root of the abundance 
da2_s <- sqrt(da2)

# convert data to samples x taxa format 
count <- as.matrix(t(da2))
count <- count[rowSums(count) > 0,]

# DMM ----
print(Sys.time())
fit <- lapply(4, dmn, count = count, verbose = TRUE)
#fit <- mclapp1y(1:50, dmn, count = count, verbose = TRUE, mc.cores = 10)
saveRDS(fit, file = paste0(outDir, Sys.Date(), "_DMM_adults_filtered_gtdb_k4.Rds"))
#save.image(file = paste0(outDir, "DMM_rarefy_sqrt_gtdb_healhty.Rdata"))
print(Sys.time())
