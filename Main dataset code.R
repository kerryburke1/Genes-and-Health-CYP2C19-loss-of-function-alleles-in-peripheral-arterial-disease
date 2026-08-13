library(tidyverse)
library(data.table)

setwd("/home/ivm/")

#all patients with PAD (up to 05/02/2024)
PAD= fread("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/Peripheral_arterial_disease//Peripheral_arterial_disease_summary_report.csv")

#add in prescription data
Meds<-fread("/genesandhealth/library-red/genesandhealth/phenotypes_rawdata/DSA__Discovery_7CCGs/2024_07_Discovery/gh3_medication_ord.csv")

Meds [, clopi := ifelse (grepl("clopidogrel|Clopidogrel|plavix|Plavix", original_term,ignore.case=TRUE)
                         ,"1", "0")]

#rename nhs number column
Meds <- Meds %>% rename(nhs_number="pseudo_nhs_number")

#remove other drug prescriptions
Meds2 <- Meds[clopi == 1] 

#merge clopi prescriptions into episode per patient
setorder(Meds2, nhs_number, clinical_effective_date)
Meds2[, prescription_date := as.Date(clinical_effective_date)] #make date
Meds2[, quantity_value := as.numeric(as.character(quantity_value))] #make numeric

cutoff <- as.Date("2024-02-05") 
Meds2<- Meds2[prescription_date <= cutoff] #remove prescriptions after 5th Feb 2024

episodes <- Meds2[, .(
  first_med = first(prescription_date),
  last_med = last(prescription_date) + last(quantity_value)
), by = nhs_number]

saveRDS(episodes, "rx_episodes.rds")
episodes <- readRDS("rx_episodes.rds")

#link prescription data to PAD
PAD_clopi = PAD%>% left_join(episodes, by = "nhs_number")
nrow(PAD_clopi)

#remove patients with no clopi prescription
PAD_clopi <- PAD_clopi[!is.na(first_med)]
nrow(PAD_clopi)

saveRDS(PAD_clopi, "PAD_clopi.rds")
PAD_clopi <- readRDS("PAD_clopi.rds")

#add in limb events (MALE = major adverse limb event)
OPCSBarts= fread("/genesandhealth/library-red/genesandhealth/phenotypes_rawdata/DSA__BartsHealth_NHS_Trust/2024_09_ResearchDataset/RDE_APC_OPCS.ascii.redacted.tab")

OPCSBarts$MALE= ifelse(OPCSBarts$OPCS_Proc_Cd=='X091'|OPCSBarts$OPCS_Proc_Cd=='X092'|OPCSBarts$OPCS_Proc_Cd=='X093'
                       |OPCSBarts$OPCS_Proc_Cd=='X094'|OPCSBarts$OPCS_Proc_Cd=='X095'|OPCSBarts$OPCS_Proc_Cd=='X098'
                       |OPCSBarts$OPCS_Proc_Cd=='X099'|OPCSBarts$OPCS_Proc_Cd=='X101'|OPCSBarts$OPCS_Proc_Cd=='X102'
                       |OPCSBarts$OPCS_Proc_Cd=='X103'|OPCSBarts$OPCS_Proc_Cd=='X104'|OPCSBarts$OPCS_Proc_Cd=='X108'
                       |OPCSBarts$OPCS_Proc_Cd=='X109'|OPCSBarts$OPCS_Proc_Cd=='X111'|OPCSBarts$OPCS_Proc_Cd=='X112'
                       |OPCSBarts$OPCS_Proc_Cd=='X118'|OPCSBarts$OPCS_Proc_Cd=='X119'
                       
                       |OPCSBarts$OPCS_Proc_Cd=='L161'|OPCSBarts$OPCS_Proc_Cd=='L162'|OPCSBarts$OPCS_Proc_Cd=='L163'
                       |OPCSBarts$OPCS_Proc_Cd=='L168'
                       |OPCSBarts$OPCS_Proc_Cd=='L169'|OPCSBarts$OPCS_Proc_Cd=='L204'|OPCSBarts$OPCS_Proc_Cd=='L205'
                       |OPCSBarts$OPCS_Proc_Cd=='L206'|OPCSBarts$OPCS_Proc_Cd=='L214'|OPCSBarts$OPCS_Proc_Cd=='L215'
                       |OPCSBarts$OPCS_Proc_Cd=='L216'
                       |OPCSBarts$OPCS_Proc_Cd=='L222'|OPCSBarts$OPCS_Proc_Cd=='L223'
                       |OPCSBarts$OPCS_Proc_Cd=='L251'|OPCSBarts$OPCS_Proc_Cd=='L252'|OPCSBarts$OPCS_Proc_Cd=='L253'
                       |OPCSBarts$OPCS_Proc_Cd=='L261'|OPCSBarts$OPCS_Proc_Cd=='L262'|OPCSBarts$OPCS_Proc_Cd=='L263'
                       |OPCSBarts$OPCS_Proc_Cd=='L265'|OPCSBarts$OPCS_Proc_Cd=='L266'|OPCSBarts$OPCS_Proc_Cd=='L267'
                       |OPCSBarts$OPCS_Proc_Cd=='L268'|OPCSBarts$OPCS_Proc_Cd=='L269'
                       
                       
                       |OPCSBarts$OPCS_Proc_Cd=='L501'
                       |OPCSBarts$OPCS_Proc_Cd=='L502'|OPCSBarts$OPCS_Proc_Cd=='L503'|OPCSBarts$OPCS_Proc_Cd=='L504'
                       |OPCSBarts$OPCS_Proc_Cd=='L505'|OPCSBarts$OPCS_Proc_Cd=='L506'|OPCSBarts$OPCS_Proc_Cd=='L508'
                       |OPCSBarts$OPCS_Proc_Cd=='L509'|OPCSBarts$OPCS_Proc_Cd=='L511'|OPCSBarts$OPCS_Proc_Cd=='L512'
                       |OPCSBarts$OPCS_Proc_Cd=='L513'|OPCSBarts$OPCS_Proc_Cd=='L514'|OPCSBarts$OPCS_Proc_Cd=='L515'
                       |OPCSBarts$OPCS_Proc_Cd=='L516'|OPCSBarts$OPCS_Proc_Cd=='L518'|OPCSBarts$OPCS_Proc_Cd=='L519'
                       |OPCSBarts$OPCS_Proc_Cd=='L521'|OPCSBarts$OPCS_Proc_Cd=='L522'|OPCSBarts$OPCS_Proc_Cd=='L528'
                       |OPCSBarts$OPCS_Proc_Cd=='L529'|OPCSBarts$OPCS_Proc_Cd=='L531'|OPCSBarts$OPCS_Proc_Cd=='L532'
                       |OPCSBarts$OPCS_Proc_Cd=='L541'|OPCSBarts$OPCS_Proc_Cd=='L542'|OPCSBarts$OPCS_Proc_Cd=='L544'
                       |OPCSBarts$OPCS_Proc_Cd=='L548'|OPCSBarts$OPCS_Proc_Cd=='L549'
                       |OPCSBarts$OPCS_Proc_Cd=='L581'|OPCSBarts$OPCS_Proc_Cd=='L582'|OPCSBarts$OPCS_Proc_Cd=='L583'
                       |OPCSBarts$OPCS_Proc_Cd=='L584'|OPCSBarts$OPCS_Proc_Cd=='L585'|OPCSBarts$OPCS_Proc_Cd=='L586'
                       |OPCSBarts$OPCS_Proc_Cd=='L587'|OPCSBarts$OPCS_Proc_Cd=='L588'|OPCSBarts$OPCS_Proc_Cd=='L589'
                       |OPCSBarts$OPCS_Proc_Cd=='L591'|OPCSBarts$OPCS_Proc_Cd=='L592'|OPCSBarts$OPCS_Proc_Cd=='L593'
                       |OPCSBarts$OPCS_Proc_Cd=='L594'|OPCSBarts$OPCS_Proc_Cd=='L595'|OPCSBarts$OPCS_Proc_Cd=='L596'
                       |OPCSBarts$OPCS_Proc_Cd=='L597'|OPCSBarts$OPCS_Proc_Cd=='L598'|OPCSBarts$OPCS_Proc_Cd=='L599'
                       |OPCSBarts$OPCS_Proc_Cd=='L601'|OPCSBarts$OPCS_Proc_Cd=='L602'|OPCSBarts$OPCS_Proc_Cd=='L603'
                       |OPCSBarts$OPCS_Proc_Cd=='L604'|OPCSBarts$OPCS_Proc_Cd=='L608'|OPCSBarts$OPCS_Proc_Cd=='L609'
                    
                       |OPCSBarts$OPCS_Proc_Cd=='L621'|OPCSBarts$OPCS_Proc_Cd=='L622'|OPCSBarts$OPCS_Proc_Cd=='L631'
                       |OPCSBarts$OPCS_Proc_Cd=='L632'|OPCSBarts$OPCS_Proc_Cd=='L633'
                       |OPCSBarts$OPCS_Proc_Cd=='L635'|OPCSBarts$OPCS_Proc_Cd=='L638'|OPCSBarts$OPCS_Proc_Cd=='L639'
                       |OPCSBarts$OPCS_Proc_Cd=='L651'|OPCSBarts$OPCS_Proc_Cd=='L652'|OPCSBarts$OPCS_Proc_Cd=='L653'
                       |OPCSBarts$OPCS_Proc_Cd=='L658'|OPCSBarts$OPCS_Proc_Cd=='L659'
                       
                       |OPCSBarts$OPCS_Proc_Cd=='L661'|OPCSBarts$OPCS_Proc_Cd=='L662'|OPCSBarts$OPCS_Proc_Cd=='L665'
                       |OPCSBarts$OPCS_Proc_Cd=='L667'|OPCSBarts$OPCS_Proc_Cd=='L668'|OPCSBarts$OPCS_Proc_Cd=='L669',
                       
                       '1','0')

#only keep MALE events
OPCSBarts <- OPCSBarts %>%
  filter(MALE ==1)
nrow(OPCSBarts)

#rename nhs number column
OPCSBarts <- OPCSBarts %>% rename (nhs_number = `PseudoNHS_2024-07-10`)

#remove duplicate events recorded on same day
OPCSBarts[, OPCS_Proc_Dt := as.Date(OPCS_Proc_Dt, format = "%d/%m/%Y %H:%M")] #make date
dup_dates <- OPCSBarts[, .N, by = .(nhs_number, OPCS_Proc_Dt)][N>1]
view(OPCSBarts[dup_dates, on = .(nhs_number, OPCS_Proc_Dt)]) #view surgeries for same patient on same day

setorder(OPCSBarts, nhs_number, OPCS_Proc_Dt)
OPCSBarts <- OPCSBarts[, .SD[1], by=.(nhs_number, OPCS_Proc_Dt)] #remove duplicates, keeping first event per day

#link MALE to PAD_clopi
PAD_clopi_MALE = PAD_clopi %>% left_join(OPCSBarts, by = "nhs_number")

#copy of dataset for MALE analysis
PAD_clopi_MALE2 <- PAD_clopi_MALE

#find MALE events during clopi episode period
PAD_clopi_MALE2 [, OPCS_Proc_Dt := as.Date(OPCS_Proc_Dt, format = "%d/%m/%Y %H:%M")] #make date

PAD_clopi_MALE2[, exposed := fifelse(OPCS_Proc_Dt > first_med &
                                       OPCS_Proc_Dt <= last_med, 1L, 0L)]

#merge into main MALE dataset
PAD_clopi_MALE[PAD_clopi_MALE2, on = "nhs_number", exposed := i.exposed]

#make MALE event dates not linked to clopi episode "NA"
PAD_clopi_MALE[, OPCS_Proc_Dt := fifelse(!is.na(OPCS_Proc_Dt) & OPCS_Proc_Dt > first_med & OPCS_Proc_Dt <= last_med, OPCS_Proc_Dt, as.Date(NA))]

#return to one row per patient, with multiple MALE events moved into extra columns
setorder(PAD_clopi_MALE, nhs_number, OPCS_Proc_Dt)
PAD_clopi_MALE[, event_rank := seq_len(.N), by = nhs_number]

PAD_clopi_MALE3 <- dcast(
  PAD_clopi_MALE,
  nhs_number ~ event_rank,
  value.var = c("OPCS_Proc_Dt", "Proc_Desc")
)

setnames(
  PAD_clopi_MALE3,
  old = names(PAD_clopi_MALE3)[-1],
  new = paste0("OPCS_Proc_Dt", names(PAD_clopi_MALE3)[-1])
)

PAD_clopi_MALE4 <- PAD_clopi_MALE[, .SD[1], by = nhs_number]

PAD_clopi_MALE5 <- PAD_clopi_MALE4[PAD_clopi_MALE3, on = "nhs_number"]

#add stroke events
all_stroke <- readRDS("all_stroke.rds")

PAD_clopi_MALE_stroke = PAD_clopi_MALE5 %>% left_join(all_stroke, by = "nhs_number")

#find stroke events during clopi episode period
PAD_clopi_MALE_stroke[, date.y := as.Date(date.y)] #make dates
PAD_clopi_MALE_stroke[, Activity_date1 := as.Date(Activity_date1)] 
PAD_clopi_MALE_stroke[, Activity_date2 := as.Date(Activity_date2)]
PAD_clopi_MALE_stroke[, Activity_date3 := as.Date(Activity_date3)]
PAD_clopi_MALE_stroke[, Activity_date4 := as.Date(Activity_date4)]
PAD_clopi_MALE_stroke[, Activity_date5 := as.Date(Activity_date5)]
PAD_clopi_MALE_stroke[, Activity_date6 := as.Date(Activity_date6)]

PAD_clopi_MALE_stroke[, exposed2 := as.integer(
  (!is.na(date.y) & date.y > first_med & date.y <= last_med)|
    (!is.na(Activity_date1) & Activity_date1 > first_med & Activity_date1 <= last_med)|
    (!is.na(Activity_date2) & Activity_date2 > first_med & Activity_date2 <= last_med)|
    (!is.na(Activity_date3) & Activity_date3 > first_med & Activity_date3 <= last_med)|
    (!is.na(Activity_date4) & Activity_date4 > first_med & Activity_date4 <= last_med)|
    (!is.na(Activity_date5) & Activity_date5 > first_med & Activity_date5 <= last_med)|
    (!is.na(Activity_date6) & Activity_date6 > first_med & Activity_date6 <= last_med)
)]

#make stroke dates outside the clopi window NA
PAD_clopi_MALE_stroke[, date.y := fifelse(!is.na(date.y) & date.y > first_med & date.y <= last_med, date.y, as.Date(NA))]
PAD_clopi_MALE_stroke[, Activity_date1  := fifelse(!is.na(Activity_date1) & Activity_date1 > first_med & Activity_date1 <= last_med, Activity_date1, as.Date(NA))]
PAD_clopi_MALE_stroke[, Activity_date2  := fifelse(!is.na(Activity_date2) & Activity_date2 > first_med & Activity_date2 <= last_med, Activity_date2, as.Date(NA))]
PAD_clopi_MALE_stroke[, Activity_date3  := fifelse(!is.na(Activity_date3) & Activity_date3 > first_med & Activity_date3 <= last_med, Activity_date3, as.Date(NA))]
PAD_clopi_MALE_stroke[, Activity_date4  := fifelse(!is.na(Activity_date4) & Activity_date4 > first_med & Activity_date4 <= last_med, Activity_date4, as.Date(NA))]
PAD_clopi_MALE_stroke[, Activity_date5  := fifelse(!is.na(Activity_date5) & Activity_date5 > first_med & Activity_date5 <= last_med, Activity_date5, as.Date(NA))]
PAD_clopi_MALE_stroke[, Activity_date6  := fifelse(!is.na(Activity_date6) & Activity_date6 > first_med & Activity_date6 <= last_med, Activity_date6, as.Date(NA))]

#find MI events
MI= fread("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/CustomPhenotypes/custom_phenotypes/GNH0005_MyocardialInfarction_extended/GNH0005_MyocardialInfarction_extended_summary_report.csv")  

MI<-MI %>%  mutate(MI=1)
MI[is.na(date)] #no NA dates

#find recurrent MI
rMI= fread("/genesandhealth/library-red/genesandhealth/phenotypes_curated/version008_2024_02/3digitICD10/3-digit-ICD/I22/I22_summary_report.csv")  

rMI<-rMI %>%  mutate(rMI=1)
rMI[is.na(date)] #no NA dates

#link MI to rMI
allMI = MI %>% left_join(rMI, by = "nhs_number")

#link all MI to PAD_clopi_MALE_stroke
PAD_clopi_MALE_stroke_MI = PAD_clopi_MALE_stroke %>% left_join(allMI, by = "nhs_number")

#find MI events during clopi episode period
PAD_clopi_MALE_stroke_MI [, date.x.y := as.Date(date.x.y)] #make date
PAD_clopi_MALE_stroke_MI [, date.y.y := as.Date(date.y.y)] #make date

PAD_clopi_MALE_stroke_MI[date.x.y == date.y.y, date.y.y := as.Date(NA)] #make dates NA for any recurrent MIs already counted as first MI (events on same day)

PAD_clopi_MALE_stroke_MI[, exposed3 := as.integer(
  (!is.na(date.x.y) & date.x.y > (first_med) & date.x.y <= last_med)|
    (!is.na(date.y.y) & date.y.y > (first_med) & date.y.y <= last_med)
)] #find MI events in clopi window

#make MI dates outside the clopi window NA
PAD_clopi_MALE_stroke_MI[, date.x.y := fifelse(!is.na(date.x.y) & date.x.y > (first_med) & date.x.y <= last_med, date.x.y, as.Date(NA))]
PAD_clopi_MALE_stroke_MI[, date.y.y := fifelse(!is.na(date.y.y) & date.y.y > (first_med) & date.y.y <= last_med, date.y.y, as.Date(NA))]

#find deaths
death= fread("/genesandhealth/library-red/genesandhealth/phenotypes_rawdata/DSA__Discovery_7CCGs/2024_07_Discovery/gh3_demographicsDataset_forTRE.csv")  
death[,deceased := fifelse(YearOfDeath == "NULL", 0, 1)]
death <- death %>%
  filter(deceased ==1) #only keep deaths
death[is.na(YearOfDeath)] #no NA dates

#rename nhs number column
death <- death %>% rename (nhs_number = "PseudoNHSNumber")

#link death to PAD_clopi_MALE_stroke_MI
PAD_clopi_MALE_stroke_MI_death = PAD_clopi_MALE_stroke_MI %>% left_join(death, by = "nhs_number")

#find deaths during clopi episode period
PAD_clopi_MALE_stroke_MI_death [, YearOfDeath := as.Date(YearOfDeath, format = "%d/%m/%Y")] #make date

PAD_clopi_MALE_stroke_MI_death[, exposed4 := fifelse(YearOfDeath > first_med &
                                                       YearOfDeath <= last_med, 1L, 0L)]

#make death dates not linked to clopi episode "NA"
PAD_clopi_MALE_stroke_MI_death[, YearOfDeath := fifelse(!is.na(YearOfDeath) & YearOfDeath > first_med & YearOfDeath <= last_med, YearOfDeath, as.Date(NA))]

#add column to specify if patient had any MALE or MACE event during clopi prescriptions
PAD_clopi_MALE_stroke_MI_death[, any_event := as.integer(rowSums(.SD ==1, na.rm = TRUE) >0),
                               .SDcols = c("exposed", "exposed2", "exposed3", "exposed4")]

#add column with total number of clopi-related MALE or MACE events per person
date_cols <- c("OPCS_Proc_DtOPCS_Proc_Dt_1","OPCS_Proc_DtOPCS_Proc_Dt_2","OPCS_Proc_DtOPCS_Proc_Dt_3",
               "OPCS_Proc_DtOPCS_Proc_Dt_4","OPCS_Proc_DtOPCS_Proc_Dt_5","OPCS_Proc_DtOPCS_Proc_Dt_6",
               "OPCS_Proc_DtOPCS_Proc_Dt_7","OPCS_Proc_DtOPCS_Proc_Dt_8","OPCS_Proc_DtOPCS_Proc_Dt_9",
               "OPCS_Proc_DtOPCS_Proc_Dt_10","OPCS_Proc_DtOPCS_Proc_Dt_11", "YearOfDeath", "date.y.x",
               "date.x.y", "date.y.y", "Activity_date1","Activity_date2","Activity_date3","Activity_date4",
               "Activity_date5","Activity_date6")
PAD_clopi_MALE_stroke_MI_death[, n_events :=  as.integer(rowSums(!is.na(.SD))), .SDcols = date_cols]

#full events dataset
PAD_events_exposed <- PAD_clopi_MALE_stroke_MI_death

saveRDS(PAD_events_exposed, "PAD_events_exposed.rds")
PAD_events_exposed <- readRDS("PAD_events_exposed.rds")

#join with link file
link_file = read_csv("/genesandhealth/library-red/genesandhealth/2025_02_10__MegaLinkage_forTRE.csv")
link_file <-link_file%>% rename('nhs_number' = 'pseudonhs_2024-07-10')
link_file <- link_file%>%
  distinct(nhs_number, .keep_all = TRUE)
PAD_events_exposed_link = PAD_events_exposed %>% left_join(link_file,by="nhs_number")

#join with co-variates 
covars = read_tsv("/genesandhealth/library-red/genesandhealth/GSAv3EAMD/Jan2024_51k_TOPMED-r3_Imputation_b38/GNH.51170.noEthnicOutliers.covariates.20PCs.tab")
PAD_events_exposed_link_covars = PAD_events_exposed_link %>% left_join(covars,by="OrageneID")

#read in genotype
CYP2C19_genos = read_tsv("/genesandhealth/red/Emma_Maia_Trios/CYP2C19.raw")

CYP2C19_genos = CYP2C19_genos %>% mutate(PM = 
                                           ifelse(`10_94781859_G_A_G`<0.5|`10_94780653_G_A_G`<0.5 |((`10_94781859_G_A_G`<1.5 & `10_94781859_G_A_G`>0.5 ) & (`10_94780653_G_A_G`<1.5 & `10_94780653_G_A_G`>0.5)),"PM","other"))


CYP2C19_genos= CYP2C19_genos %>% mutate(PI = ifelse(`10_94781859_G_A_G`<1.5|`10_94780653_G_A_G`<1.5,"PI","other"))

CYP2C19_genos = CYP2C19_genos %>% mutate(IM = ifelse((`10_94781859_G_A_G`<1.5 & `10_94781859_G_A_G`>0.5 & `10_94780653_G_A_G`>1.5)|(`10_94780653_G_A_G`<1.5 &`10_94780653_G_A_G`>0.5 & `10_94781859_G_A_G`>1.5),"IM","other"))

CYP2C19_genos = CYP2C19_genos %>%mutate(UR = ifelse(`10_94761900_C_T_C`<0.5 & `10_94780653_G_A_G`>1.5 & `10_94781859_G_A_G`>1.5,"UR","other"))

CYP2C19_genos = CYP2C19_genos %>%mutate(RM = ifelse(`10_94761900_C_T_C`<1.5 & `10_94780653_G_A_G`>1.5 & `10_94781859_G_A_G`>1.5,"RM","other"))

#join genotype with covars 
PAD_events_exposed_link_covars_geno = PAD_events_exposed_link_covars %>% left_join(CYP2C19_genos,by="IID")

#remove patients with no genotype data
PAD_events_exposed_link_covars_geno <- PAD_events_exposed_link_covars_geno[PM=="PM"|PM=="other"] 

#join with questionnaire data
q_file = read_csv("/genesandhealth/library-red/genesandhealth/phenotypes_rawdata/QMUL__Stage1Questionnaire/2025_04_25__S1QSTredacted.csv")
q_file <- q_file %>% rename (OrageneID = S1QST_Oragene_ID)

q_file <- q_file%>%
  distinct(OrageneID, .keep_all = TRUE)

PAD_events_exposed_link_covars_geno = PAD_events_exposed_link_covars_geno %>% left_join(q_file,by="OrageneID")

#add column for age at first prescription
library(lubridate)
PAD_events_exposed_link_covars_geno[, dob := as.Date(paste0("01-", `S1QST_MM-YYYY_ofBirth`), format = "%d-%m-%Y")] #convert to dob
PAD_events_exposed_link_covars_geno[, age_first_med := floor(time_length(interval(dob, first_med), "years"))]

saveRDS(PAD_events_exposed_link_covars_geno, "PAD_events_exposed_link_covars_geno")
PAD_events_exposed_link_covars_geno <- readRDS("PAD_events_exposed_link_covars_geno")
