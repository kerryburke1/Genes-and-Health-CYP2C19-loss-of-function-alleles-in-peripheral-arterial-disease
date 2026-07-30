#initial and recurrent event analysis
PAD_BMI_DM_CKD <- readRDS("PAD_BMI_DM_CKD")

#binary lof score
PAD_BMI_DM_CKD[, lof :=
                 fifelse(PI == "other", 0,
                         fifelse(IM == "IM", 1,
                                 fifelse(PM == "PM", 2, NA_integer_
                                 )))]

#logistic regression - for any event - lof only (age included as con-founder)
logmet <- glm(any_event ~ PI
              + as.numeric(age_first_med)
              + as.factor(S1QST_Gender.y)
              + as.numeric(PC1)
              + as.numeric(PC2),
              family = binomial,
              data = PAD_BMI_DM_CKD
)
summary(logmet)

cbind(
  OR = exp(coef(logmet)),
  exp(confint(logmet))
)

fit <- glm(any_event ~ as.factor(lof) + age_first_med + S1QST_Gender.y + PC1 + PC2,    #with parametised lof
              family = binomial,
              data = PAD_BMI_DM_CKD
)
summary(fit) 

#logistic regression for any event - ABCD-GENE score (no age as cofounder, as included in score)
logmet <- glm(any_event ~ score
              + as.factor(S1QST_Gender.y)
              + as.numeric(PC1)
              + as.numeric(PC2),
              family = binomial,
              data = PAD_BMI_DM_CKD
)
summary(logmet)

cbind(
  OR = exp(coef(logmet)),
  exp(confint(logmet))
)

logmet <- glm(any_event ~ risk_score   #using >=10 threshold
              + as.factor(S1QST_Gender.y)
              + as.numeric(PC1)
              + as.numeric(PC2),
              family = binomial,
              data = PAD_BMI_DM_CKD
)
summary(logmet)

cbind(
  OR = exp(coef(logmet)),
  exp(confint(logmet))
)

#logistic regression - for second (recurrent) event
PAD_BMI_DM_CKD$event_date <- pmin(PAD_BMI_DM_CKD$date.y.y,
                                  PAD_BMI_DM_CKD$date.x.y,
                                  PAD_BMI_DM_CKD$date.y.x,
                                  PAD_BMI_DM_CKD$Activity_date1,
                                  PAD_BMI_DM_CKD$Activity_date2,
                                  PAD_BMI_DM_CKD$Activity_date3,
                                  PAD_BMI_DM_CKD$Activity_date4,
                                  PAD_BMI_DM_CKD$Activity_date5,
                                  PAD_BMI_DM_CKD$Activity_date6,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_1,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_2,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_3,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_4,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_5,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_6,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_7,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_8,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_9,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_10,
                                  PAD_BMI_DM_CKD$OPCS_Proc_DtOPCS_Proc_Dt_11,
                                  na.rm = TRUE)     #death excluded from first event

PAD_BMI_DM_CKD2 <- subset(PAD_BMI_DM_CKD, !is.na(event_date)) #limit data set to those with first event

PAD_BMI_DM_CKD2 [, event_date := as.Date(event_date)] #make date

#find second events (death now included)
PAD_BMI_DM_CKD2[, event_date2 :=
                  apply(.SD,1,function(x) {
                    x<- sort (x[!is.na(x)])
                    if(length(x) >=2) x[2] else as.Date(NA)
                  }),
                .SDcols = c("OPCS_Proc_DtOPCS_Proc_Dt_1","OPCS_Proc_DtOPCS_Proc_Dt_2","OPCS_Proc_DtOPCS_Proc_Dt_3",
                            "OPCS_Proc_DtOPCS_Proc_Dt_4","OPCS_Proc_DtOPCS_Proc_Dt_5","OPCS_Proc_DtOPCS_Proc_Dt_6",
                            "OPCS_Proc_DtOPCS_Proc_Dt_7","OPCS_Proc_DtOPCS_Proc_Dt_8","OPCS_Proc_DtOPCS_Proc_Dt_9",
                            "OPCS_Proc_DtOPCS_Proc_Dt_10","OPCS_Proc_DtOPCS_Proc_Dt_11", "YearOfDeath", "date.y.x",
                            "date.x.y", "date.y.y", "Activity_date1","Activity_date2","Activity_date3","Activity_date4",
                            "Activity_date5","Activity_date6")]

PAD_BMI_DM_CKD2 [, event_date2 := as.Date(event_date2)] #make date

#PAD_BMI_DM_CKD2[event_date2 > (event_date + 365), event_date2 := as.Date(NA)] #make dates not within a year of first event n/a (code added as needed for sensitivity analysis)

PAD_BMI_DM_CKD2$both_dates <- as.integer(!is.na(PAD_BMI_DM_CKD2$event_date) & !is.na(PAD_BMI_DM_CKD2$event_date2))

#run for lof
logmet_lof <- glm(both_dates ~ as.factor(lof)
                  + as.numeric(age_first_med)
                  + as.factor(S1QST_Gender.y)
                  + as.numeric(PC1)
                  + as.numeric(PC2),
                  family = binomial,
                  data = PAD_BMI_DM_CKD2
)
summary(logmet_lof)

cbind(
  OR = exp(coef(logmet_lof)),
  exp(confint(logmet_lof))
)

fit.lof <- glm(both_dates ~ as.factor(lof) + age_first_med + S1QST_Gender.y + PC1 + PC2,   #with parametised lof
           family = binomial,
           data = PAD_BMI_DM_CKD2
)
summary(fit.lof)

cbind(
  OR = exp(coef(fit.lof)),
  exp(confint(fit.lof))
)

#run for ABCD-GENE score
logmet_ABCD <- glm(both_dates ~ score
                   + as.factor(S1QST_Gender.y)
                   + as.numeric(PC1)
                   + as.numeric(PC2),
                   family = binomial,
                   data = PAD_BMI_DM_CKD2
)
summary(logmet_ABCD)

cbind(
  OR = exp(coef(logmet_ABCD)),
  exp(confint(logmet_ABCD))
)

logmet_ABCD_binary <- glm(both_dates ~ risk_score             #using >=10 threshold
                   + as.factor(S1QST_Gender.y)
                   + as.numeric(PC1)
                   + as.numeric(PC2),
                   family = binomial,
                   data = PAD_BMI_DM_CKD2
)
summary(logmet_ABCD_binary)

cbind(
  OR = exp(coef(logmet_ABCD_binary)),
  exp(confint(logmet_ABCD_binary))
)

#compare PM vs ABCD-GENE score for recurrent events
logmet_PM <- glm(both_dates ~ PM        #PM vs not-PM, no age co-variate for comparison purposes
                  + as.factor(S1QST_Gender.y)
                  + as.numeric(PC1)
                  + as.numeric(PC2),
                  family = binomial,
                  data = PAD_BMI_DM_CKD2
)
summary(logmet_PM) 

cbind(
  OR = exp(coef(logmet_PM)),
  exp(confint(logmet_PM))
)

logmet_both <- glm(both_dates ~ PM + score
                   + as.factor(S1QST_Gender.y)
                   + as.numeric(PC1)
                   + as.numeric(PC2),
                   family = binomial,
                   data = PAD_BMI_DM_CKD2
)
summary(logmet_both)

cbind(
  OR = exp(coef(logmet_both)),
  exp(confint(logmet_both))
)

AIC(logmet_PM, logmet_ABCD, logmet_both)

#logistic regression for individual sections of ABCD-GENE score (component analysis)
logmet_sections <- glm(both_dates ~ as.factor(lof)
                       + as.factor(age_score)
                       + as.factor(BMI)
                       + as.factor(CKD)
                       + as.factor(DM_score),
                       family = binomial,
                       data = PAD_BMI_DM_CKD2
)
summary(logmet_sections)

cbind(
  OR = exp(coef(logmet_sections)),
  exp(confint(logmet_sections))
)

#calculate optimum ABCD-GENE score threshold for this cohort
library(pROC)

roc_obj <- roc(PAD_BMI_DM_CKD2$both_dates, PAD_BMI_DM_CKD2$score)

coords(roc_obj,
       x = "best",
       best.method = "youden",
       ret = c("threshold", "sensitivity", "specificity"))

#AUC
auc(roc_obj)

ci.auc(roc_obj)
