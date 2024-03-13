# Import Library
library(tidyverse)
library(lubridate)


################################################################################# I ### Load data
path <- fs::path("", "Volumes", "Gillis_Research","Christelle Colin-Leitzinger", "CH and HIV")

id_to_remove <- read_csv(paste0(here::here(), "/id_to_remove_from_analysis.csv"))
                         
clinical_data <- 
  # readxl::read_xlsx(paste0(here::here(), "/Sample_Sheet_withEpiAge_20230316.xlsx")#, # For full data
  # readxl::read_xlsx(paste0(here::here(), "/hiv_data.xlsx")#, # For short data
  readxl::read_xlsx(paste0(here::here(), "/allSample_wClocks_M4Madded_20230509.xlsx")#, # For new clocks
                    ) %>% 
  janitor::clean_names() %>% 
  filter(study_id != id_to_remove$id_to_remove | is.na(study_id))

dictionary <- 
  readxl::read_xlsx(paste0(here::here(), "/Table 1_new variables-modified04.04.23.xlsx"),
                    sheet = "new variables-modified"
  )

smoking_CR <- readxl::read_xlsx(paste0(here::here(), "/HIVreport_for_smoking.xlsx"),
                             sheet = "Smoking_from_CR"
)
smoking_CERNER <- readxl::read_xlsx(paste0(here::here(), "/HIVreport_for_smoking.xlsx"),
                                sheet = "Smoking_from_Cerner"
)
smoking_patientquestionnaire <- readxl::read_xlsx(paste0(here::here(), "/HIVreport_for_smoking.xlsx"),
                                sheet = "Summary_for_smoking"
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

smoking <- smoking_CR %>% 
  select(PATIENT_ID, DERIVED_TOBACCO_SMOKING_STATUS_DESC) %>% 
  distinct() %>% 
  arrange(PATIENT_ID, DERIVED_TOBACCO_SMOKING_STATUS_DESC) %>% 
  distinct(PATIENT_ID, .keep_all = TRUE) %>% 
  full_join(., smoking_CERNER %>% 
              select(MRN, PATIENT_ID, RECENTLY_REPORTED_CIGARETTE_SMOKING_SV2_DESC),
            by = "PATIENT_ID") %>% 
  full_join(., smoking_patientquestionnaire %>% 
              select(PATIENT_ID, RECENTLY_REPORTED_CIGARETTE_SMOKING_SV2_DESC),
            by = "PATIENT_ID") %>% 
  mutate(DERIVED_TOBACCO_SMOKING_STATUS_DESC = case_when(
    DERIVED_TOBACCO_SMOKING_STATUS_DESC == "Unknown"             ~ NA_character_,
    DERIVED_TOBACCO_SMOKING_STATUS_DESC == "Current smoker"      ~ "Ever",
    DERIVED_TOBACCO_SMOKING_STATUS_DESC == "Former smoker"       ~ "Ever",
    DERIVED_TOBACCO_SMOKING_STATUS_DESC == "Never smoker"        ~ "Never"
  )) %>% 
  mutate(across(c("RECENTLY_REPORTED_CIGARETTE_SMOKING_SV2_DESC.x",
                  "RECENTLY_REPORTED_CIGARETTE_SMOKING_SV2_DESC.y"), 
                ~ case_when(
                  . == "Missing"  ~ NA_character_,
                  TRUE                                                         ~ .
                )
                  )) %>% 
  mutate(across(-c("MRN", "PATIENT_ID"), ~factor(., levels = c("Never", "Ever"))))


data_change <- read_csv(paste0(here::here(), "/id_to_change_data_tx.csv"))
data_change2 <- read_csv(paste0(here::here(), "/id_to_change_data_tx_after.csv"))
data_change3 <- read_csv(paste0(here::here(), "/id_to_change_data_date.csv"))

CH_HIV_data <- 
  bind_rows(clinical_data_A, clinical_data_B, clinical_data_C) %>% 
  select(chip_barcode, sample_id_for_chip, sample_id_from_ab,CH_status, MRN, patient_id, batch, everything()) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "Race") %>% select(Category, new_race = `New Category`),
            by= c("race" = "Category")) %>% 
  mutate(new_race = factor(new_race, levels= c("White", "Black", "Other"))) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "Ethnicity") %>% select(Category, new_ethnicity = `New Category`),
            by= c("ethnicity" = "Category")) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "Primary Site") %>% select(Category, new_primary_site = `New Category`),
            by= c("primary_site" = "Category")) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "TNM Stage") %>% select(Category, new_tnm_stage = `New Category`),
            by= c("tnm_stage" = "Category")) %>% 
  left_join(.,
            dictionary %>% filter(Variable == "Tx Summary") %>% select(Category, new_tx_summary = `New Category`),
            by= c("tx_summary" = "Category")) %>% 
  left_join(., treatment, 
            by= c("chip_barcode")) %>% 
  left_join(., smoking,
            by= "MRN") %>% 
  mutate(sex = factor(sex, levels= c("MALE", "FEMALE"))) %>% 
  mutate(treatment_status = case_when(
    MRN == data_change$id_change_data      ~ "before",
    MRN == data_change2$id_change_data     ~ "after start",
    TRUE                                   ~ treatment_status
  )) %>% 
  mutate(treatment_date = case_when(
    MRN == data_change3$id_change_data     ~ "2019-05-21",
    TRUE                                   ~ treatment_date
  )) %>% 

  mutate(treatment_status = factor(treatment_status, levels = c("before", "after start"))) %>% 
  mutate(os_event = case_when(
    vital_status == 1               ~ 0,
    vital_status == 2               ~ 1
  )) %>% 
  mutate(across(c("treatment_date", "specimen_to_surgery_lag", "specimen_to_treatment_lag"), ~ na_if(., "NA"))) %>% 
  mutate(across(c("specimen_to_surgery_lag", "specimen_to_treatment_lag"), ~ as.numeric(.))) %>% 
  # mutate(treatment_date = na_if(treatment_date, "NA")) %>% 
  mutate(across(c("dob", "specimen_date", "treatment_date", "vital_status_date"), ~as.Date(.))) %>% 
  mutate(age_at_specimen = case_when(
    !is.na(age_at_specimen)           ~ age_at_specimen,
    is.na(age_at_specimen)            ~ round(interval(start = dob, end = specimen_date)/
                                     duration(n=1, units = "years"), 0)
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
  mutate(tnm_stage = case_when(
    tnm_stage == "99 UNK"           ~ NA_character_,
    str_detect(tnm_stage, "1")      ~ "I",
    str_detect(tnm_stage, "2")      ~ "II",
    str_detect(tnm_stage, "3")      ~ "III",
    str_detect(tnm_stage, "4")      ~ "IV"
  )) %>% 
  mutate(primary_site = case_when(
    str_detect(primary_site, "TONGUE")                   ~ "Tongue",
    str_detect(primary_site, "CHEEK MUCOSA")             ~ "Mouth",
    str_detect(primary_site, "INTRAHEPATIC BILE DUC")    ~ "Liver",
    str_detect(primary_site, "ESOPHAGUS")                ~ "Esophagus",
    
    str_detect(primary_site, "GASTRIC")                  ~ "Stomach",
    str_detect(primary_site, "STOMACH")                  ~ "Stomach",
    str_detect(primary_site, "DUODENUM")                 ~ "Stomach",
    
    str_detect(primary_site, "COLON")                    ~ "Colorectal",
    str_detect(primary_site, "APPENDIX")                 ~ "Colorectal",
    str_detect(primary_site, "RECTUM")                   ~ "Colorectal",
    str_detect(primary_site, "ANUS")                     ~ "Colorectal",
    str_detect(primary_site, "ANAL")                     ~ "Colorectal",
    
    str_detect(primary_site, "PANCREAS")                 ~ "Pancreas",
    str_detect(primary_site, "LUNG")                     ~ "Lung",
    str_detect(primary_site, "SKIN")                     ~ "Skin",
    str_detect(primary_site, "BREAST")                   ~ "Breast",
    TRUE                                                 ~ primary_site
  )) %>% 
  mutate(CH_status = factor(CH_status, levels = c("NO", "CH")))
  
# ################################################################################# III ### Merge data
# CH_HIV_data <- clinical_data %>% 
#   full_join(.,
#             CH_status1 %>% 
#               select(MRN, patient_id, chip_barcode, sample_id_from_ab, new_sample_id_for_chip,
#                      everything()),
#             by = c("chip_barcode", "sample_id_from_ab"))








write_rds(CH_HIV_data, "CH_HIV_data_full_03112024.rds")
write_csv(CH_HIV_data, "CH_HIV_data_03112024.csv")

# End cleaning