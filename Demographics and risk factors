#PAD demographics
PAD_events_exposed_link_covars_geno <- readRDS("PAD_events_exposed_link_covars_geno")

#age
shapiro.test(PAD_events_exposed_link_covars_geno$age_first_med) #age is normally distributed

PAD_events_exposed_link_covars_geno2 <- PAD_events_exposed_link_covars_geno[PI == "PI"]
PAD_events_exposed_link_covars_geno2[, .(
  mean_age = mean(age_first_med, na.rm = TRUE),
  sd_age = sd(age_first_med, na.rm = TRUE)
             )]

PAD_events_exposed_link_covars_geno3 <- PAD_events_exposed_link_covars_geno[PI == "other"]
PAD_events_exposed_link_covars_geno3[, .(
  mean_age = mean(age_first_med, na.rm = TRUE),
  sd_age = sd(age_first_med, na.rm = TRUE)
)]

t.test(age_first_med ~ PI, data = PAD_events_exposed_link_covars_geno)

#adding co-morbidity
PAD_events_exposed_link_covars_geno <- readRDS("PAD_events_exposed_link_covars_geno")

#link with CHD
CHD = read_csv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/Coronary_heart_disease/Coronary_heart_disease_summary_report.csv")
CHD <- as.data.table(CHD)
CHD[, CHD := 1]
PAD_BMI_DM_CKD_CHD = PAD_BMI_DM_CKD %>% left_join(CHD, by="nhs_number")

#link with CVD
CVD = read_csv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/Cerebrovascular_disease/Cerebrovascular_disease_summary_report.csv")
CVD <- as.data.table(CVD)
CVD[, CVD := 1]
PAD_BMI_DM_CKD_CVD = PAD_events_exposed_link_covars_geno %>% left_join(CVD, by="nhs_number")

#link with HTN
HTN = read_csv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/Hypertension/Hypertension_summary_report.csv")
HTN <- as.data.table(HTN)
HTN[, HTN := 1]
PAD_BMI_DM_CKD_HTN = PAD_events_exposed_link_covars_geno %>% left_join(HTN, by="nhs_number")
