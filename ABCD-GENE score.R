#ABCD Gene Score
PAD_events_exposed_link_covars_geno <- readRDS("PAD_events_exposed_link_covars_geno")

#link with BMI
BMI = read_csv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/Obesity/Obesity_summary_report.csv")
BMI <- as.data.table(BMI)
BMI[, BMI := 1]
PAD_BMI = PAD_events_exposed_link_covars_geno %>% left_join(BMI, by="nhs_number")

#link with DM
T1DM = read_csv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/Type_1_Diabetes/Type_1_Diabetes_summary_report.csv")
T1DM <- as.data.table(T1DM)
T1DM[, T1DM := 1]
PAD_BMI_DM = PAD_BMI%>% left_join(T1DM, by="nhs_number")

T2DM = read_csv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/Type_2_Diabetes/Type_2_Diabetes_summary_report.csv")
T2DM <- as.data.table(T2DM)
T2DM[, T2DM := 1]
T2DM <- T2DM %>% rename(nhs_number="nhs_number")
PAD_BMI_DM = PAD_BMI_DM%>% left_join(T2DM, by="nhs_number")

unspecDM = read_csv("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/Unspecified_or_Rare_Diabetes/Unspecified_or_Rare_Diabetes_summary_report.csv")
unspecDM <- as.data.table(unspecDM)
unspecDM[, unspecDM := 1]
PAD_BMI_DM = PAD_BMI_DM%>% left_join(unspecDM, by="nhs_number")

#link with CKD
CKD_data = fread("/genesandhealth/library-red/genesandhealth/phenotypes_rawdata/DSA__BartsHealth_NHS_Trust/2024_09_ResearchDataset/RDE_APC_DIAGNOSIS.ascii.redacted.tab")

CKD_data$CKD= ifelse(CKD_data$ICD_Diagnosis_Cd=='N183'|CKD_data$ICD_Diagnosis_Cd=='N184'|CKD_data$ICD_Diagnosis_Cd=='N185',
                       '1','0')

#rename nhs number column
CKD_data <- CKD_data %>% rename(nhs_number="PseudoNHS_2024-07-10")

#only keep CKD patients
CKD_data <- CKD_data %>%
  filter(CKD ==1)

#return to one row per patient
CKD_data <- CKD_data[, .SD[1], by = nhs_number]

#link CKD to dataset
PAD_BMI_DM_CKD = PAD_BMI_DM %>% left_join(CKD_data, by = "nhs_number")

#calculate total ABCD Gene Score per person
PAD_BMI_DM_CKD[, score :=
                 (age_first_med >75) * 4 +
                 (BMI == 1 & !is.na(BMI)) *4 +
                 (CKD ==1 & !is.na(CKD)) *3 +
                 (((T1DM ==1 & !is.na(T1DM)) | (T2DM ==1 & !is.na(T2DM)) | (unspecDM ==1 & !is.na(unspecDM)))) *3 +
                 (!is.na(IM) & IM == "IM") *6 +
                 (!is.na(PM) & PM == "PM") *24
               ]

#split into individual sections
PAD_BMI_DM_CKD[, DM_score := as.integer(T1DM ==1 | T2DM ==1 | unspecDM ==1)]
PAD_BMI_DM_CKD[, age_score := as.integer(age_first_med >75)]

PAD_BMI_DM_CKD[is.na(BMI), BMI :=0]
PAD_BMI_DM_CKD[is.na(CKD), CKD :=0]
PAD_BMI_DM_CKD[is.na(DM_score), DM_score :=0]
PAD_BMI_DM_CKD[is.na(age_score), age_score :=0]

#make binary group where score >= 10
PAD_BMI_DM_CKD[, risk_score := as.integer(score>=10)]

saveRDS(PAD_BMI_DM_CKD, "PAD_BMI_DM_CKD")
