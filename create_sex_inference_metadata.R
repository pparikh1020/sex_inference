setwd("/Users/parikhpp/git/sex_inference")

# Import Data
eyeIntegration22 <- 
  read.csv("/Users/parikhpp/NEI Projects/metadata_creation/eyeIntegration22_meta_2022_10_25.csv")
# Downloaded from the SRA
gender_data <- read.csv("gender_metadata.csv")

# Only keep samples in EyeIntegration
eyeIntegration22_gtex <- eyeIntegration22 %>% filter(study_accession == "SRP012682")
gender_data_gtex <- gender_data %>% filter(SRA.Study == "SRP012682")

eyeIntegration22_not_gtex <- eyeIntegration22 %>% filter(study_accession != "SRP012682")
gender_data_not_gtex <- gender_data %>% filter(SRA.Study != "SRP012682")

gtex_keep <- intersect(eyeIntegration22_gtex$gtex_accession, gender_data_gtex$Sample.Name)
gtex_data <- gender_data_gtex %>% filter(Sample.Name %in% gtex_keep)

not_gtex_keep <- intersect(eyeIntegration22_not_gtex$run_accession, gender_data_not_gtex$Run)
not_gtex_data <- gender_data_not_gtex %>% filter(Run %in% not_gtex_keep)

gender_metadata <- bind_rows(not_gtex_data, gtex_data)

##### Subset for ocular samples where gender is given

sex_inf <- gender_metadata %>% 
  filter(SRA.Study!="SRP012682") %>% 
  filter(sex != "") %>% 
  filter(sex != "not applicable") %>% 
  filter(sex != "not determined")

sex_inf <- rename(sex_inf, run_accession = Run)

# Pull in old metadata
old_ei <- read.csv("https://hpc.nih.gov/~parikhpp/EiaD/2022_metadata.csv")

# Combine data

sex_inf_meta <- old_ei %>% left_join(sex_inf, by = "run_accession")
sex_inf_meta <- sex_inf_meta %>% filter(!(is.na(sex)))

write.csv(sex_inf_meta, "sex_inference_metadata.csv")
