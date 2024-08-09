# PERMANOVA ANALYSIS

library(tidyverse)
library(here)
library(vegan)
library(futile.logger)
library(parallel)

# data 
flog.info("loading the data")
outDir <- here("results/")
dataDir <- here("data/")

dataMatrix <- readRDS(here("data/motusMatrix_Genuslevel_gtdb_adults_filtered.rds"))
dist_full <- readRDS(here("data/motusMatrix_Genuslevel_dist_JSD.rds"))
metadata <-  readRDS(here("data/metadata.rds"))
#----------------------------------------------------------------------------------------
flog.info("Country Effect")
md <- metadata %>%
  mutate(subject_id = case_when(is.na(subject_id) ~ internal_id,
                                TRUE ~ subject_id)) %>%
  filter(!is.na(geographic_location)) %>%
  ungroup()

dist <- as.matrix(dist_full)
dist <- dist[rownames(dist) %in% md$internal_id, colnames(dist) %in% md$internal_id]
print(dim(dist))
md <- md %>% filter(internal_id %in% colnames(dist))
dist <- as.dist(dist)

flog.info("calculate PERMANOVA")

ad_country <-   adonis2(dist ~ geographic_location, data = md, by = "margin", permutations = 999, parallel = 20)

flog.info("save PERMANOVA")
saveRDS(ad_country, paste0(outDir, "permanova_dist_marginal_999perm_country.rds"))
#----------------------------------------------------------------------------------------
flog.info("Study Effect")
md <- metadata %>%
  mutate(subject_id = case_when(is.na(subject_id) ~ internal_id,
                                TRUE ~ subject_id)) %>%
  filter(!is.na(studies)) %>%
  ungroup()

dist <- as.matrix(dist_full)
dist <- dist[rownames(dist) %in% md$internal_id, colnames(dist) %in% md$internal_id]
print(dim(dist))
md <- md %>% filter(internal_id %in% colnames(dist))
dist <- as.dist(dist)

flog.info("calculate PERMANOVA")


ad_studies <-    adonis2(dist ~ studies, data = md, by = "margin", permutations = 999, parallel = 20)

flog.info("save PERMANOVA")
saveRDS(ad_studies, paste0(outDir, "permanova_dist_marginal_999perm_studies.rds"))
#----------------------------------------------------------------------------------------
flog.info("Disease Effect")
md <- metadata %>%
  mutate(subject_id = case_when(is.na(subject_id) ~ internal_id,
                                TRUE ~ subject_id)) %>%
  filter(!is.na(disease)) %>%
  filter(!is.na(pam_k3)) %>% ungroup()

dist <- as.matrix(dist_full)
dist <- dist[rownames(dist) %in% md$internal_id, colnames(dist) %in% md$internal_id]
md <- md %>% filter(internal_id %in% colnames(dist))
print(dim(dist))
dist <- as.dist(dist)

flog.info("calculate PERMANOVA")

ad_disease <- adonis2(dist ~ disease, data = md, by = "margin", permutations = 999, parallel = 20)

flog.info("save PERMANOVA")
saveRDS(ad_disease, paste0(outDir, "permanova_dist_marginal_500perm_disease.rds"))
#----------------------------------------------------------------------------------------
flog.info("subject_id Effect")
md <- metadata %>%
  filter(!is.na(subject_id)) %>%
  filter(!is.na(pam_k3)) %>% ungroup()

dist <- as.matrix(dist_full)
dist <- dist[rownames(dist) %in% md$internal_id, colnames(dist) %in% md$internal_id]
md <- md %>% filter(internal_id %in% colnames(dist))
print(dim(dist))
dist <- as.dist(dist)

flog.info("calculate PERMANOVA")

ad_subject <-adonis2(dist ~ subject_id, data = md, by = "margin", permutations = 999, parallel = 20)

flog.info("save PERMANOVA")
saveRDS(ad_subject, paste0(outDir, "permanova_dist_marginal_500perm_subject.rds"))
#----------------------------------------------------------------------------------------
flog.info("age Effect")
md <- metadata %>%
  mutate(subject_id = case_when(is.na(subject_id) ~ internal_id,
                                TRUE ~ subject_id)) %>%
  filter(!is.na(age_years)) %>%
  filter(!is.na(pam_k3)) %>% ungroup()

dist <- as.matrix(dist_full)
dist <- dist[rownames(dist) %in% md$internal_id, colnames(dist) %in% md$internal_id]
md <- md %>% filter(internal_id %in% colnames(dist))
print(dim(dist))
dist <- as.dist(dist)

flog.info("calculate PERMANOVA")

ad_age <- adonis2(dist ~  age_years, data = md, by = "margin", permutations = 999, parallel = 20)

flog.info("save PERMANOVA")
saveRDS(ad_age, paste0(outDir, "permanova_dist_marginal_500perm_age.rds"))
#----------------------------------------------------------------------------------------
flog.info("bmi Effect")
md <- metadata %>%
  mutate(subject_id = case_when(is.na(subject_id) ~ internal_id,
                                TRUE ~ subject_id)) %>%
  filter(!is.na(bmi)) %>%
  filter(!is.na(pam_k3)) %>% ungroup()

dist <- as.matrix(dist_full)
dist <- dist[rownames(dist) %in% md$internal_id, colnames(dist) %in% md$internal_id]
md <- md %>% filter(internal_id %in% colnames(dist))
print(dim(dist))
dist <- as.dist(dist)

flog.info("calculate PERMANOVA")

ad_bmi <- adonis2(dist ~ bmi, data = md, by = "margin", permutations = 999, parallel = 20)

flog.info("save PERMANOVA")
saveRDS(ad_bmi, paste0(outDir, "permanova_dist_marginal_500perm_bmi.rds"))
#----------------------------------------------------------------------------------------
flog.info("bristol stool Effect")
md <- metadata %>%
  mutate(subject_id = case_when(is.na(subject_id) ~ internal_id,
                                TRUE ~ subject_id)) %>%
  filter(!is.na(bristol_stool_scale))

dist <- as.matrix(dist_full)
dist <- dist[rownames(dist) %in% md$ena_ers_sample_id, colnames(dist) %in% md$ena_ers_sample_id]
md <- md %>% filter(ena_ers_sample_id %in% colnames(dist))
print(dim(dist))
dist <- as.dist(dist)

flog.info("calculate PERMANOVA")

ad_bristol <- adonis2(dist ~ bristol_stool_scale, data = md, by = "margin", permutations = 999, parallel = 20)

flog.info("save PERMANOVA")
saveRDS(ad_bristol, paste0(outDir, "permanova_dist_marginal_500perm_bristol.rds"))
#----------------------------------------------------------------------------------------
flog.info("Smoker")
md <- metadata %>%
  mutate(subject_id = case_when(is.na(subject_id) ~ internal_id,
                                TRUE ~ subject_id)) %>%
  filter(!is.na(smoker))

dist <- as.matrix(dist_full)
dist <- dist[rownames(dist) %in% md$ena_ers_sample_id, colnames(dist) %in% md$ena_ers_sample_id]
md <- md %>% filter(ena_ers_sample_id %in% colnames(dist))
print(dim(dist))
dist <- as.dist(dist)

flog.info("calculate PERMANOVA")

ad_smoker <- adonis2(dist ~ smoker, data = md, by = "margin", permutations = 999, parallel = 20)

flog.info("save PERMANOVA")
saveRDS(ad_smoker,paste0(outDir, "permanova_dist_marginal_500perm_smoker.rds"))
