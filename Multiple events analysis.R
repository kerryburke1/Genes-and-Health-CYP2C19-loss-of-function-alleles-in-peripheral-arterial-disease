#MCF for multiple events
library(tidyverse)
library(data.table)

library(survival)
library(survminer)

PAD_BMI_DM_CKD <- readRDS("PAD_BMI_DM_CKD")

PAD_BMI_DM_CKD[, first_med := as.Date(first_med)] #make date
PAD_BMI_DM_CKD[, last_med := as.Date(last_med)] #make date

PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_1 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_1)] #make date
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_2 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_2)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_3 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_3)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_4 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_4)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_5 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_5)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_6 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_6)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_7 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_7)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_8 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_8)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_9 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_9)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_10 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_10)] 
PAD_BMI_DM_CKD[, OPCS_Proc_DtOPCS_Proc_Dt_11 := as.Date(OPCS_Proc_DtOPCS_Proc_Dt_11)] 
PAD_BMI_DM_CKD [, date.y.y := as.Date(date.y.y)] #make date
PAD_BMI_DM_CKD [, Activity_date1 := as.Date(Activity_date1)]
PAD_BMI_DM_CKD [, Activity_date2 := as.Date(Activity_date2)]
PAD_BMI_DM_CKD [, Activity_date3 := as.Date(Activity_date3)]
PAD_BMI_DM_CKD [, Activity_date4 := as.Date(Activity_date4)]
PAD_BMI_DM_CKD [, Activity_date5 := as.Date(Activity_date5)]
PAD_BMI_DM_CKD [, Activity_date6 := as.Date(Activity_date6)]
PAD_BMI_DM_CKD [, date.x.y := as.Date(date.x.y)] #make date
PAD_BMI_DM_CKD [, date.y.x := as.Date(date.y.x)] #make date
PAD_BMI_DM_CKD [, YearOfDeath := as.Date(YearOfDeath, format = "%d/%m/%Y")]

#multiple events negative binomial regression
PAD_BMI_DM_CKD$risk_window <- as.numeric(PAD_BMI_DM_CKD$last_med
                                - PAD_BMI_DM_CKD$first_med)
library(MASS)

#Lof binary score
PAD_BMI_DM_CKD[, lof :=
                 fifelse(PI == "other", 0,
                         fifelse(IM == "IM", 1,
                                 fifelse(PM == "PM", 2, NA_integer_
                                 )))]

#run for multiple events and lof (negative binomial regression)
fit <- glm.nb(n_events ~ PI
              + as.numeric(PAD_BMI_DM_CKD$age_first_med)
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2)
              + offset(log(risk_window)),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

fit <- glm.nb(n_events ~ as.factor(lof)     #with parametised lof
              + as.numeric(PAD_BMI_DM_CKD$age_first_med)
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2)
              + offset(log(risk_window)),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

#run for multiple events and ABCD-GENE score (negative binomial regression)
fit <- glm.nb(n_events ~ score 
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2)
              + offset(log(risk_window)),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

fit <- glm.nb(n_events ~ risk_score 
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2)
              + offset(log(risk_window)),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

#make column for n_events within 2 years only
date_cols <- c("OPCS_Proc_DtOPCS_Proc_Dt_1","OPCS_Proc_DtOPCS_Proc_Dt_2","OPCS_Proc_DtOPCS_Proc_Dt_3",
               "OPCS_Proc_DtOPCS_Proc_Dt_4","OPCS_Proc_DtOPCS_Proc_Dt_5","OPCS_Proc_DtOPCS_Proc_Dt_6",
               "OPCS_Proc_DtOPCS_Proc_Dt_7","OPCS_Proc_DtOPCS_Proc_Dt_8","OPCS_Proc_DtOPCS_Proc_Dt_9",
               "OPCS_Proc_DtOPCS_Proc_Dt_10","OPCS_Proc_DtOPCS_Proc_Dt_11", "YearOfDeath", "date.y.x",
               "date.x.y", "date.y.y", "Activity_date1","Activity_date2","Activity_date3","Activity_date4",
               "Activity_date5","Activity_date6")
PAD_BMI_DM_CKD[, n_events2yr := rowSums(sapply(.SD, function(x) {
  !is.na(x) & x > first_med & x <= first_med +730
  }), na.rm = TRUE),
  .SDcols = date_cols]

#retest association of lof and binary ABCD-GENE score at 2 years
fit <- glm.nb(n_events2yr ~ risk_score
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

fit <- glm.nb(n_events2yr ~ PI
              + as.numeric(PAD_BMI_DM_CKD$age_first_med)
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

fit <- glm.nb(n_events2yr ~ as.factor(lof)
              + as.numeric(PAD_BMI_DM_CKD$age_first_med)
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

#make column for n_events within 5 years only
date_cols <- c("OPCS_Proc_DtOPCS_Proc_Dt_1","OPCS_Proc_DtOPCS_Proc_Dt_2","OPCS_Proc_DtOPCS_Proc_Dt_3",
               "OPCS_Proc_DtOPCS_Proc_Dt_4","OPCS_Proc_DtOPCS_Proc_Dt_5","OPCS_Proc_DtOPCS_Proc_Dt_6",
               "OPCS_Proc_DtOPCS_Proc_Dt_7","OPCS_Proc_DtOPCS_Proc_Dt_8","OPCS_Proc_DtOPCS_Proc_Dt_9",
               "OPCS_Proc_DtOPCS_Proc_Dt_10","OPCS_Proc_DtOPCS_Proc_Dt_11", "YearOfDeath", "date.y.x",
               "date.x.y", "date.y.y", "Activity_date1","Activity_date2","Activity_date3","Activity_date4",
               "Activity_date5","Activity_date6")

PAD_BMI_DM_CKD[, n_events5yr := rowSums(sapply(.SD, function(x) {
  !is.na(x) & x > first_med & x <= first_med +1825
}), na.rm = TRUE),
.SDcols = date_cols]

#retest association of lof and binary ABCD-GENE score at 5 years
fit <- glm.nb(n_events5yr ~ risk_score
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

fit <- glm.nb(n_events5yr ~ PI
              + as.numeric(PAD_BMI_DM_CKD$age_first_med)
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

fit <- glm.nb(n_events5yr ~ as.factor(lof)
              + as.numeric(PAD_BMI_DM_CKD$age_first_med)
              + as.factor(PAD_BMI_DM_CKD$S1QST_Gender.y)
              + as.numeric(PAD_BMI_DM_CKD$PC1)
              + as.numeric(PAD_BMI_DM_CKD$PC2),
              data = PAD_BMI_DM_CKD)
summary(fit)

cbind(
  IRR = coef(fit),
  confint(fit)
)

#multiple events mean cumulative function graph for lof

#wide to long
event_cols <- c("OPCS_Proc_DtOPCS_Proc_Dt_1","OPCS_Proc_DtOPCS_Proc_Dt_2","OPCS_Proc_DtOPCS_Proc_Dt_3",
                "OPCS_Proc_DtOPCS_Proc_Dt_4","OPCS_Proc_DtOPCS_Proc_Dt_5","OPCS_Proc_DtOPCS_Proc_Dt_6",
                "OPCS_Proc_DtOPCS_Proc_Dt_7","OPCS_Proc_DtOPCS_Proc_Dt_8","OPCS_Proc_DtOPCS_Proc_Dt_9",
                "OPCS_Proc_DtOPCS_Proc_Dt_10","OPCS_Proc_DtOPCS_Proc_Dt_11", "YearOfDeath", "date.y.x",
                "date.x.y", "date.y.y", "Activity_date1","Activity_date2","Activity_date3","Activity_date4",
                "Activity_date5","Activity_date6")

PAD_BMI_DM_CKD [, first_med := as.Date(first_med)] #make date

#final date is last med or death
PAD_BMI_DM_CKD$final_date <- pmin(PAD_BMI_DM_CKD$last_med,
                                  PAD_BMI_DM_CKD$YearOfDeath,
                                  na.rm = TRUE)

PAD_BMI_DM_CKD [, final_date := as.Date(final_date)] #make date

PAD_BMI_DM_CKD3 <- melt(
  PAD_BMI_DM_CKD,
  id.vars = c("nhs_number", "PI", "first_med", "final_date" ),
  measure.vars = event_cols,
  value.name = "event_date"
)

PAD_BMI_DM_CKD3 <- PAD_BMI_DM_CKD3[!is.na(event_date)]

PAD_BMI_DM_CKD3[, time := as.numeric(event_date - first_med)]

PAD_BMI_DM_CKD[, censor_time := as.numeric(final_date - first_med)]

PAD_BMI_DM_CKD3 <- merge(
  PAD_BMI_DM_CKD3,
  PAD_BMI_DM_CKD[, .(nhs_number, censor_time)],
  by = "nhs_number",
  all.x = TRUE
)

PAD_BMI_DM_CKD3 <- PAD_BMI_DM_CKD3[time <= censor_time]

PAD_BMI_DM_CKD_MCF <- PAD_BMI_DM_CKD3[, .N, by= .(time, PI)]
setorder(PAD_BMI_DM_CKD_MCF, PI, time)

n_group <- PAD_BMI_DM_CKD[, .N, by = PI]
PAD_BMI_DM_CKD_MCF <- merge(PAD_BMI_DM_CKD_MCF, n_group, by = "PI")
PAD_BMI_DM_CKD_MCF[, mcf := cumsum(N.x) / N.y, by = PI]

PAD_BMI_DM_CKD_MCF$time <- PAD_BMI_DM_CKD_MCF$time/365.25 #change time to years not days

plot(NULL,
     xlab = "Time since first clopidogrel prescription (years)",
     ylab = "Mean cummulative number of events",
     xlim = range(PAD_BMI_DM_CKD_MCF$time),
     ylim = range(PAD_BMI_DM_CKD_MCF$mcf))

cols <- c("blue", "red")
i <- 1

for (g in unique(PAD_BMI_DM_CKD_MCF$PI)) {
  lines(PAD_BMI_DM_CKD_MCF[PI == g]$time,
        PAD_BMI_DM_CKD_MCF[PI == g]$mcf,
        col = cols[i],
        lwd = 2)
  i <- i +1
}

legend("bottomright",
       legend = expression(italic("CYP2C19") ~ "LoF alleles", "No LoF alleles"),
       col = cols,
       lwd = 2)

#run for ABCD-GENE score

#wide to long
event_cols <- c("OPCS_Proc_DtOPCS_Proc_Dt_1","OPCS_Proc_DtOPCS_Proc_Dt_2","OPCS_Proc_DtOPCS_Proc_Dt_3",
                "OPCS_Proc_DtOPCS_Proc_Dt_4","OPCS_Proc_DtOPCS_Proc_Dt_5","OPCS_Proc_DtOPCS_Proc_Dt_6",
                "OPCS_Proc_DtOPCS_Proc_Dt_7","OPCS_Proc_DtOPCS_Proc_Dt_8","OPCS_Proc_DtOPCS_Proc_Dt_9",
                "OPCS_Proc_DtOPCS_Proc_Dt_10","OPCS_Proc_DtOPCS_Proc_Dt_11", "YearOfDeath", "date.y.x",
                "date.x.y", "date.y.y", "Activity_date1","Activity_date2","Activity_date3","Activity_date4",
                "Activity_date5","Activity_date6")

PAD_BMI_DM_CKD [, first_med := as.Date(first_med)] #make date

#end date is last med or death
PAD_BMI_DM_CKD$final_date <- pmin(PAD_BMI_DM_CKD$last_med,
                                  PAD_BMI_DM_CKD$YearOfDeath,
                                  na.rm = TRUE)

PAD_BMI_DM_CKD [, final_date := as.Date(final_date)] #make date

PAD_BMI_DM_CKD4 <- melt(
  PAD_BMI_DM_CKD,
  id.vars = c("nhs_number", "risk_score", "first_med", "final_date" ),
  measure.vars = event_cols,
  value.name = "event_date"
)

PAD_BMI_DM_CKD4 <- PAD_BMI_DM_CKD4[!is.na(event_date)]

PAD_BMI_DM_CKD4[, time := as.numeric(event_date - first_med)]

PAD_BMI_DM_CKD[, censor_time := as.numeric(final_date - first_med)]

PAD_BMI_DM_CKD4 <- merge(
  PAD_BMI_DM_CKD4,
  PAD_BMI_DM_CKD[, .(nhs_number, censor_time)],
  by = "nhs_number",
  all.x = TRUE
)

PAD_BMI_DM_CKD4 <- PAD_BMI_DM_CKD4[time <= censor_time]

PAD_BMI_DM_CKD_MCF <- PAD_BMI_DM_CKD4[, .N, by= .(time, risk_score)]
setorder(PAD_BMI_DM_CKD_MCF, risk_score, time)

n_group <- PAD_BMI_DM_CKD[, .N, by = risk_score]
PAD_BMI_DM_CKD_MCF <- merge(PAD_BMI_DM_CKD_MCF, n_group, by = "risk_score")
PAD_BMI_DM_CKD_MCF[, mcf := cumsum(N.x) / N.y, by = risk_score]

PAD_BMI_DM_CKD_MCF$time <- PAD_BMI_DM_CKD_MCF$time/365.25 #change time to years not days

plot(NULL,
     xlab = "Time since first clopidogrel prescription (years)",
     ylab = "Mean cummulative number of events",
     xlim = range(PAD_BMI_DM_CKD_MCF$time),
     ylim = range(PAD_BMI_DM_CKD_MCF$mcf))

cols <- c("red","blue")
i <- 1

for (g in unique(PAD_BMI_DM_CKD_MCF$risk_score)) {
  lines(PAD_BMI_DM_CKD_MCF[risk_score == g]$time,
        PAD_BMI_DM_CKD_MCF[risk_score == g]$mcf,
        col = cols[i],
        lwd = 2)
  i <- i +1
}

legend("bottomright",
       legend = c( "ABCD-GENE score <10", "ABCD-GENE score \u2265 10"),
       col = cols,
       lwd = 2)
