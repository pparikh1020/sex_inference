# XIST ENSG00000229807 -- y chromosome
# EIF1AY ENSG00000198692 -- y chromosome
# KDM5D ENSG00000012817 -- y chromosome
# UTY ENSG00000183878 -- y chromosome
# DDX3Y ENSG00000067048 -- y chromosome
# RPS4Y1 ENSG00000129824 -- y chromosome

library(tximport)
library(GenomicFeatures)
library(tidyr)
library(tidyverse)
library(dplyr)
library(caret)
library(glmnet)
library(janitor)

old_ei <- read.csv("https://hpc.nih.gov/~parikhpp/EiaD/2022_metadata.csv")

#Imported from the create_sex_inference_metadata.R script output
sex_inf_metadata <- 
  read.csv("/Users/parikhpp/NEI Projects/sex_inference/sex_inference_metadata.csv")

# Download gene counts files from the gene_counts folder in the EiaD_build repository
gene_counts_recount3 <- 
  vroom::vroom("/Users/parikhpp/git/EiaD_build/gene_counts/recount3_transformed_counts.csv")
gene_counts_local_additions <- 
  vroom::vroom("/Users/parikhpp/git/EiaD_build/gene_counts/local_data_additions_transformed_counts.csv")

# Combine recount3 and local samples for sex inference
combined_gene_counts <- gene_counts_local_additions %>% left_join(gene_counts_recount3, by = "gene_id")

genes <- c("ENSG00000229807.11", "ENSG00000198692.9", "ENSG00000012817.15", "ENSG00000183878.15",
           "ENSG00000067048.16", "ENSG00000129824.15")

combined_gene_counts_subset <- combined_gene_counts %>% filter(gene_id %in% genes)

# Log normalization
log_normalized_counts <- combined_gene_counts_subset %>% select(-gene_id) %>% +1 %>%  log()

# Only use recount3 samples included in eyeIntegration
log_normalized_counts_keep <- intersect(names(log_normalized_counts), sex_inf_metadata$run_accession)
combined_gene_counts_subset <- combined_gene_counts_subset %>% select(one_of(c("gene_id", log_normalized_counts_keep)))
#Only use samples which pass quality check
qc_keep <- old_ei %>% filter(Kept == "Kept") %>% pull(run_accession)
combined_gene_counts_subset_quality <- combined_gene_counts_subset %>% select(one_of(c("gene_id", qc_keep)))

#Transform data long-ways

gene_counts_long <- t(combined_gene_counts_subset_quality) %>% as.data.frame()
gene_counts_long <- gene_counts_long %>% row_to_names(row_number = 1)
gene_counts_long <- tibble::rownames_to_column(gene_counts_long, "run_accession")

# Add sex data

regression_data <- gene_counts_long %>% left_join(sex_inf_metadata %>% select(run_accession, sex), by = "run_accession")
colnames(regression_data) = c("run_accession", "XIST","EIF1AY","KDM5D","UTY","DDX3Y","RPS4Y1","sex")

# Add mapping data and percentages

recount3_mapping_information <- vroom::vroom("/Users/parikhpp/git/EiaD_build/mapping_data/recount3_mapping_information.csv")
local_data_additions_mapping_information <- vroom::vroom("/Users/parikhpp/git/EiaD_build/mapping_data/local_data_additions_mapping_information.csv")

mapping_data <- bind_rows(recount3_mapping_information, local_data_additions_mapping_information)

regression_data_final <- regression_data %>% 
  left_join(mapping_data %>% 
              select(recount_qc.aligned_reads..chrx, recount_qc.aligned_reads..chry, external_id), 
            by = c("run_accession" = "external_id"))

write.csv(regression_data_final[,2:8], "inference_regression_data.csv", row.names = FALSE)

# Plotting data

# X alignment
ggplot(regression_data_final, aes(x = sex, y = recount_qc.aligned_reads..chrx, color = sex)) +
  geom_boxplot() + theme_bw()
# Y alignment
ggplot(regression_data_final, aes(x = sex, y = recount_qc.aligned_reads..chry, color = sex)) +
  geom_boxplot() + theme_bw()
# X/Y ratio
ggplot(regression_data_final, 
       aes(x = sex, y = (recount_qc.aligned_reads..chrx/recount_qc.aligned_reads..chry), color = sex)) +
  geom_boxplot() + theme_bw()

# Which samples are overlapping?

regression_data_final <- mutate(regression_data_final, 
                                x_y_ratio = (recount_qc.aligned_reads..chrx/recount_qc.aligned_reads..chry))
regression_data_final %>% filter(sex == "female") %>% filter(x_y_ratio < 100)

# SRR2895372 - SRP064956
# SRR2895378 - SRP064956
# SRR2895380 - SRP064956
# SRR2895381 - SRP064956
# SRR4252548 - SRP090027
# SRR5535471 - SRP090027
# SRR5535778 - SRP090027
# SRR5535934 - SRP090027

# The above samples are in the same tail range as the male plot.

# Investigating whether these are the same outliers as for the female only plot

ggplot(regression_data_final, aes(x = sex, y = recount_qc.aligned_reads..chry, color = sex)) +
  geom_boxplot() + theme_bw()
regression_data_final %>% filter(sex == "female") %>% filter(recount_qc.aligned_reads..chry > 0.01)

# It seems as if though the same 8 samples, including an additional ninth overlap between male and female samples

# SRR7461319 - SRP151763
# SRR2895372
# SRR2895378
# SRR2895380
# SRR2895381
# SRR4252548
# SRR5535471
# SRR5535778
# SRR5535934

# All samples were confirmed female.
