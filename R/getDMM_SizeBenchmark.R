library(DirichletMultinomial)
library(magrittr)
library(dplyr)
library(here)

outDir <- paste0(here(), "results/")
dataDir <- paste0(here(), "data/")

set.seed(132940)

# LOAD THE DATA ---- 
# load metadata 
md <- readr::read_tsv(paste0(dataDir, "metadata_.rds"))

data<- readRDS(paste0(dataDir, "motusMatrix_Genuslevel_adults_filtered.rds"))
data <- data[row.names(data)!="-1",]
data <- na.omit(data)


## PREPARE THE MOTUS-DATA -------------------------------
# remove rows with only zero 
data <- data[rowSums(data[, -1])>0,]

# normalization - transform to relative abundance, multiply with the lowest read number (sum by sample) or someth.
# and round to a real number.  this was also suggested by Sara Vieira -Silva. as done in the metacardis paper (Jeroen Raes lab)

da2 <- round(data * 10000, 0)

# square root of the abundance 
da2_s <- sqrt(da2)

# convert data to samples x taxa format 
count <- as.matrix(t(da2_s))
count <- count[rowSums(count) > 0,]


# PREPARE DIFFERENT SIZES BUT FULL HETEROGENITY OF THE DATASET 

size.total <- md %>% pull(internal_id)

size.half.1 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/2))
size.half.2<- md %>% pull(internal_id) %>% sample(., round(length(size.total)/2))
size.half.3 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/2))

size.quarter.1 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/4))
size.quarter.2 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/4))
size.quarter.3 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/4))

size.eighth.1 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/8))
size.eighth.2 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/8))
size.eighth.3 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/8))

size.sixteenth.1 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/16))
size.sixteenth.2 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/16))
size.sixteenth.3 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/16))

size.demisemiq.1 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/32))
size.demisemiq.2 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/32))
size.demisemiq.3 <- md %>% pull(internal_id) %>% sample(., round(length(size.total)/32))

dmm_sizecohorts <- list("size.total" = size.total, 
                       # "size.threequarter.1"=size.threequarter.1, "size.threequarter.2"=size.threequarter.2, "size.threequarter.3"=size.threequarter.3,
                        "size.half.1"=size.half.1, "size.half.2"=size.half.2, "size.half.3"=size.half.3, 
                        "size.quarter.1"=size.quarter.1, "size.quarter.2"=size.quarter.2, "size.quarter.3"=size.quarter.3,
                        "size.eighth.1"=size.eighth.1, "size.eighth.2"=size.eighth.2, "size.eighth.3"=size.eighth.3,
                        "size.sixteenth.1"=size.sixteenth.1, "size.sixteenth.2"=size.sixteenth.2, "size.sixteenth.3"=size.sixteenth.3,
                       "size.demisemiq.1" = size.demisemiq.1, "size.demisemiq.2" = size.demisemiq.2, "size.demisemiq.3" = size.demisemiq.3)

saveRDS(dmm_sizecohorts, here("fulldata_v2/data/processed/adults_filtered_gtdb/dmm_benchmark_SizeSubcohorts_internal_ids.rds"))
## DMM FOR THE SIZE SUBCOHORTS ---- 

for(i in 1:length(dmm_sizecohorts)){
    i=1
  tmp <- count[rownames(count) %in% dmm_sizecohorts[[i]], ]
  print(Sys.time())
  print(paste0("Calculating DMM for ", names(dmm_sizecohorts)[[i]]))
  fit <- lapply(1:15, dmn, count = tmp, verbose = TRUE)
  saveRDS(fit, file = paste0(outDir, "DMM_SizeBenchmark_", names(dmm_sizecohorts)[[i]], ".Rds"))
}
