#Data clean and imputation
##################################################################################
#Load packages
library(tidyverse)
library(softImpute)

#Load datasets
data <- readRDS("master_processed_data_v20240825.RDS")
data_raw <- readRDS("master_processed_data_v20240825.RDS")
data_har <- readRDS("master_harmonized_data_v20240825.RDS")

#Separate training and challenge data
challenge_id <- data_har$challenge$subject_specimen$specimen_id
train_id <- data_har$training$subject_specimen$specimen_id

train_tbl <- data$subject_specimen %>% filter(specimen_id %in% train_id)
challenge_tbl <- data$subject_specimen %>% filter(specimen_id %in% challenge_id)

#Features
data$plasma_ab_titer_raw <- data$plasma_ab_titer$raw_data
data$plasma_ab_titer <- data$plasma_ab_titer$batchCorrected_data
data$plasma_cytokine_concentrations_by_olink <- data$plasma_cytokine_concentrations_by_olink$batchCorrected_data
data$plasma_cytokine_concentrations_by_legendplex <- data$plasma_cytokine_concentrations_by_legendplex$normalized_data
data$pbmc_cell_frequency <- data$pbmc_cell_frequency$batchCorrected_data
data$t_cell_polarization <- data$t_cell_polarization$raw_data
data$t_cell_activation <- data$t_cell_activation$raw_data
data$pbmc_gene_expression_raw_count <- data$pbmc_gene_expression$raw_count$batchCorrected_data
data$pbmc_gene_expression_tpm <- data$pbmc_gene_expression$tpm$batchCorrected_data
data$pbmc_gene_expression <- NULL

#IgG titres data (raw data)
ab_titer_raw <- data$plasma_ab_titer_raw %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$plasma_ab_titer_raw)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  mutate(type = paste0(type, "_raw")) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(ab_titer_raw, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(ab_titer_raw, by = "specimen_id")

#IgG titres data (batch correlated data)
ab_titer <- data$plasma_ab_titer %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$plasma_ab_titer)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(ab_titer, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(ab_titer, by = "specimen_id")

#cytokine data (olink)
cytokine_olink <- data$plasma_cytokine_concentrations_by_olink %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$plasma_cytokine_concentrations_by_olink)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  mutate(type = paste0("cytokine_olink_",type)) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(cytokine_olink, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(cytokine_olink, by = "specimen_id")

#cytokine data (lp)
cytokine_lp <- data$plasma_cytokine_concentrations_by_legendplex %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$plasma_cytokine_concentrations_by_legendplex)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  mutate(type = paste0("cytokine_lp_",type)) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(cytokine_lp, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(cytokine_lp, by = "specimen_id")

#cell frequency data
cell_freq <- data$pbmc_cell_frequency %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$pbmc_cell_frequency)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(cell_freq, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(cell_freq, by = "specimen_id")

#T-cell data (Polar)
tcell_polar <- data$t_cell_polarization %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$t_cell_polarization)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(tcell_polar, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(tcell_polar, by = "specimen_id")

#T-cell data (act)
tcell_act <- data$t_cell_activation %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$t_cell_activation)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(tcell_act, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(tcell_act, by = "specimen_id")

#Gene RNAseq data (raw count)
gene_count <- data$pbmc_gene_expression_raw_count %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$pbmc_gene_expression_raw_count)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  mutate(type = paste0("raw_count_",type)) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(gene_count, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(gene_count, by = "specimen_id")

#Gene RNAseq data (tpm)
gene_tpm <- data$pbmc_gene_expression_tpm %>% as_tibble(row_name = FALSE) %>%
  mutate(type = row.names(data$pbmc_gene_expression_tpm)) %>%
  pivot_longer(cols = !type, names_to = "specimen_id", values_to = "value") %>%
  mutate(type = paste0("tpm_",type)) %>%
  pivot_wider(names_from = type, values_from = value) %>%
  mutate(specimen_id = as.integer(specimen_id))

train_tbl <- train_tbl %>% left_join(gene_tpm, by = "specimen_id")
challenge_tbl <- challenge_tbl %>% left_join(gene_tpm, by = "specimen_id")

#save the cleaned dataset
save(train_tbl, challenge_tbl, file = "data_cleaned.RData")

#Clean demographic data
demo_training <- data_har$training$subject_specimen %>% filter(timepoint == 0) %>%
  select(subject_id, specimen_type, infancy_vac, biological_sex, ethnicity, race,
         year_of_birth, date_of_boost, dataset)

demo_challenge <- data_har$challenge$subject_specimen %>% filter(timepoint == 0) %>%
  select(subject_id, specimen_type, infancy_vac, biological_sex, ethnicity, race,
         year_of_birth, date_of_boost, dataset)

#save the demographic data
save(demo_training, demo_challenge, file = "demo_data.RData")

##################################################################################
##Feature selection, and data imputation

train_tbl %>% filter(timepoint == 14) %>% select(subject_id, IgG_PT_raw) %>%
  filter(!is.na(IgG_PT_raw)) -> task1

train_tbl %>% filter(timepoint == 1) %>% select(subject_id, Monocytes) %>%  
  filter(!is.na(Monocytes)) -> task2

train_tbl %>% filter(timepoint == 1) %>% select(subject_id, tpm_ENSG00000277632.1) %>% 
  filter(!is.na(tpm_ENSG00000277632.1)) -> task3 

train_day0 <- train_tbl %>% filter(timepoint == 0)

train_y <- train_day0 %>% select(subject_id, dataset, IgG_PT_raw) %>%
  rename(IgG_t0 = IgG_PT_raw) %>%
  left_join(task1, by = "subject_id") %>%
  left_join(task2, by = "subject_id") %>%
  left_join(task3, by = "subject_id")

#Summarize the missing data in features
train_day0
var_na <- as.numeric(colSums(is.na(train_day0)))
var_na

missing_var <- tibble(variable = names(train_day0), missing = var_na,
                      prop = missing/nrow(train_day0))

missing_var %>% View()

#Summarize the missing data in subjects
sub_na <- as.numeric(rowSums(is.na(train_day0)))
missing_sub <- tibble(variable = train_day0$subject_id, missing = sub_na,
                      prop = missing/ncol(train_day0))
missing_sub %>% View()

train_day0 %>% select(subject_id, starts_with("P"), TT) %>% View()

challenge_day0 <- challenge_tbl %>% filter(timepoint == 0)

tibble(variable = names(challenge_day0), missing = as.numeric(colSums(is.na(challenge_day0))),
       prop = missing/nrow(challenge_day0)) %>% View()

tibble(variable = challenge_day0$subject_id, missing = as.numeric(rowSums(is.na(challenge_day0))),
       prop = missing/ncol(challenge_day0)) %>% View()


#EXCLUDE features and subjects with high proportion of missing data
train_day0_sub <- train_day0 %>% select(-starts_with("cytokine"), -starts_with("PHA"),
                                        -starts_with("PT"), -TT, -starts_with("raw"),
                                        -ends_with("_raw"))
tibble(variable = names(train_day0_sub), missing = as.numeric(colSums(is.na(train_day0_sub))),
       prop = missing/nrow(train_day0_sub)) %>% View()
tibble(variable = train_day0_sub$subject_id, missing = as.numeric(rowSums(is.na(train_day0_sub))),
       prop = missing/ncol(train_day0_sub)) %>% View()

train_day0_sub <- train_day0_sub[as.numeric(rowSums(is.na(train_day0_sub)))/ncol(train_day0_sub)<=0.5,]

challenge_day0_sub <- challenge_day0 %>% select(-starts_with("cytokine"), -starts_with("PHA"),
                                                -starts_with("PT"), -TT, -starts_with("raw"),
                                                -ends_with("_raw"))

load("demo_data.RData")
demo_training <- demo_training %>% mutate(age = year(date_of_boost) - year(year_of_birth)) %>%
  mutate(infancy_vac = as.numeric(infancy_vac == "wP"),
         biological_sex = as.numeric(biological_sex == "Male")) %>% 
  select(-specimen_type, -ethnicity, -year_of_birth,
         -date_of_boost, -dataset, -race)
demo_challenge <- demo_challenge %>% mutate(age = year(date_of_boost) - year(year_of_birth)) %>%
  mutate(infancy_vac = as.numeric(infancy_vac == "wP"),
         biological_sex = as.numeric(biological_sex == "Male")) %>% 
  select(-specimen_type, -ethnicity, -year_of_birth,
         -date_of_boost, -dataset, -race)

demo_training_sub <- demo_training %>% filter(subject_id %in% train_day0_sub$subject_id)
train_cleaned <- cbind(train_day0_sub[,2:3], demo_training_sub[,2:4],
                       train_day0_sub[,8:6680]) %>% as_tibble()
challenge_cleaned <- cbind(challenge_day0_sub[,2:3], demo_challenge[,2:4],
                           challenge_day0_sub[,8:6680]) %>% as_tibble()

names(train_cleaned)[c(12,19,26,33)] <- c("IgG1_FIM2_3", "IgG2_FIM2_3", "IgG3_FIM2_3",
                                          "IgG4_FIM2_3")
names(challenge_cleaned)[c(12,19,26,33)] <- c("IgG1_FIM2_3", "IgG2_FIM2_3", "IgG3_FIM2_3",
                                              "IgG4_FIM2_3")

cell_freq_index <- paste0("cell_freq_", c(1:39))
names(train_cleaned)[37:75] <- cell_freq_index
names(challenge_cleaned)[37:75] <- cell_freq_index

save(train_cleaned, challenge_cleaned, train_y, 
     file = "train_challenge_cleaned.RData")

load("train_challenge_cleaned.RData")

all_cleaned <- rbind(train_cleaned, challenge_cleaned)

#softImpute data imputation
data_ori <- as.matrix(all_cleaned[,6:6678])
imp_soft <- softImpute(data_ori, rank=10, lambda=30)
data_imp <- complete(data_ori, imp_soft)

all_cleaned[,6:6678] <- data_imp

train_softimp <- all_cleaned %>% filter(subject_id %in% train_cleaned$subject_id)
challenge_softimp <- all_cleaned %>% filter(subject_id %in% challenge_cleaned$subject_id)

save(train_softimp, challenge_softimp, train_y, file = "imputed_softImpute.RData")


