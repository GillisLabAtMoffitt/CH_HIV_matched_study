# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "CH and HIV")

clinical_data <- 
  # readxl::read_xlsx(paste0(path, "/raw data/Sample_Sheet_withEpiAge_20230316.xlsx"),
  # readxl::read_xlsx(paste0("/Users/colinccm/Documents/GitHub/Gillis", "/Sample_Sheet_withEpiAge_Methylation_Analysis_Cohort_20230321.xlsx")#,
  readxl::read_xlsx(paste0(here::here(), "/hiv_data.xlsx")#,
                                      # sheet = "Sample_Sheet_20220811"
                    ) %>% 
  janitor::clean_names() %>% 
  filter(study_id != "P30_S3_009" | is.na(study_id))

dictionary <- 
  readxl::read_xlsx(paste0(here::here(), "/Table 1_new variables-modified.xlsx"),
                    sheet = "new variables-modified"
  )

treatment <- 
  readxl::read_xlsx(paste0(here::here(), "/CH_HIV_data_treatment status_NG_03.31.23.xlsx")
  ) %>% 
  select(chip_barcode, treatment_status)

CH_status <- 
  # readxl::read_xlsx(paste0(path, "/raw data/CICPT2987-combined_calls_03.17.23_FINALwU2AF1.xlsx")) %>%
  readxl::read_xlsx(paste0(here::here(), "/CICPT2987-combined_calls_03.17.23_FINALwU2AF1.xlsx")) %>%
  rename(CH_status = CH) %>% 
  mutate(MRN = as.character(MRN))

Linked_IDs <- 
  # readxl::read_xlsx(paste0(path, "/raw data/HIV_Participant_ID_update_batch C update.xlsx")) %>% 
  readxl::read_xlsx(paste0(here::here(), "/dat C update_copy.xlsx")) %>% 
  janitor::clean_names() %>% 
  mutate(across(contains("mrn"), ~ as.character(.))) %>% 
  select(chip_barcode, batch, newmrn, pipeline_file_name,
         sample_id_from_anders, sample_id_fo_chip_from_du, 
         study_id_correct
         ) %>% 
  filter(!is.na(newmrn))


################################################################################# II ### Data cleaning
CH_status <- CH_status %>% 
  inner_join(., Linked_IDs,
            by = c("MRN" = "newmrn", 
                   "patient_id" = "pipeline_file_name")) %>% 
  select(MRN, patient_id, batch, everything())

clinical_data_B <- clinical_data %>% 
  filter(batch == "B") %>% 
  left_join(., CH_status,
             by = c("batch", 
                    "study_id" = "study_id_correct")) %>% 
  rename(chip_barcode = chip_barcode.x) %>% 
  select(-c(chip_barcode.y, sample_id_from_anders, sample_id_fo_chip_from_du))
  
clinical_data_A <- clinical_data %>% 
  filter(batch == "A") %>% 
  left_join(., CH_status,
            by = c("chip_barcode", "batch", 
                   "sample_id_from_ab" = "sample_id_from_anders",
                   "sample_id_for_chip" = "sample_id_fo_chip_from_du")) %>% 
  select(-study_id_correct)

clinical_data_C <- clinical_data %>% 
  filter(batch == "C") %>% 
  left_join(., CH_status,
            by = c("chip_barcode", "batch", 
                   "sample_id_from_ab" = "sample_id_from_anders",
                   "sample_id_for_chip" = "sample_id_fo_chip_from_du",
                   "study_id" = "study_id_correct"))

CH_HIV_data <- 
  bind_rows(clinical_data_A, clinical_data_B, clinical_data_C) %>% 
  select(chip_barcode, sample_id_for_chip, sample_id_from_ab,CH_status, MRN, patient_id, batch, everything()) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "Race") %>% select(Category, new_race = `New Category`),
            by= c("race" = "Category")) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "Ethnicity") %>% select(Category, new_ethnicity = `New Category`),
            by= c("ethnicity" = "Category")) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "Primary Site") %>% select(Category, new_primary_site = `New Category`),
            by= c("primary_site" = "Category")) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "TNM Stage") %>% select(Category, new_tnm_stage = `New Category`),
            by= c("tnm_stage" = "Category")) %>% 
  left_join(., treatment, 
            by= c("chip_barcode")) %>% 
  mutate(treatment_status = factor(treatment_status, levels = c("before", "after start"))) %>% 
  mutate(os_event = case_when(
    vital_status == 1               ~ 0,
    vital_status == 2               ~ 1
  )) %>% 
  mutate(month_at_os = interval(start = treatment_date, end = vital_status_date)/
           duration(n=1, units = "months")) %>% 
  mutate(new_ethnicity = factor(new_ethnicity, levels = c("Non-Hispanic", "Hispanic"))) %>% 
  mutate(race = factor(race, levels = c("White", "Black", "Other"))) %>% 
  mutate(new_tnm_stage = case_when(
    new_tnm_stage == "Unknown"      ~ NA_character_,
    TRUE                            ~ new_tnm_stage
  )) %>% 
  mutate(clinical_stage = case_when(
    tnm_stage == "99 UNK"           ~ NA_character_,
    str_detect(tnm_stage, "1|2")    ~ "I-II",
    str_detect(tnm_stage, "3|4")    ~ "III-IV"
  )) %>% 
  mutate(CH_status = factor(CH_status, levels = c("NO", "CH")))
  
# ################################################################################# III ### Merge data
# CH_HIV_data <- clinical_data %>% 
#   full_join(.,
#             CH_status1 %>% 
#               select(MRN, patient_id, chip_barcode, sample_id_from_ab, new_sample_id_for_chip,
#                      everything()),
#             by = c("chip_barcode", "sample_id_from_ab"))








write_rds(CH_HIV_data, "CH_HIV_data.rds")
write_csv(CH_HIV_data, "CH_HIV_data.csv")

# End cleaning