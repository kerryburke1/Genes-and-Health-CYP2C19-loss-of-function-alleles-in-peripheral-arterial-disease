#time on clopidogrel
PAD_events_exposed_link_covars_geno <- readRDS("PAD_events_exposed_link_covars_geno")

PAD_events_exposed_link_covars_geno$risk_window <- as.numeric(PAD_events_exposed_link_covars_geno$last_med
                                                              - PAD_events_exposed_link_covars_geno$first_med)

#binary lof score
PAD_events_exposed_link_covars_geno[, lof :=
                 fifelse(PI == "other", 0,
                         fifelse(IM == "IM", 1,
                                 fifelse(PM == "PM", 2, NA_integer_
                                 )))]

#time on clopi stats                                                           
PAD_events_exposed_link_covars_geno[, .(
  n = .N,
  median = median(risk_window, na.rm = TRUE),
  IQR = IQR(risk_window, na.rm = TRUE)
), by = PI]                

PAD_events_exposed_link_covars_geno[, {
  s <-shapiro.test(risk_window)
  .(
    W = unname(s$statistic),
    p.value = s$p.value
  )
}, by = lof] #not normal distributions

wilcox.test(risk_window ~ PI,
            data = PAD_events_exposed_link_covars_geno)

PAD_events_exposed_link_covars_geno2 <- PAD_events_exposed_link_covars_geno[PM=="other" | IM =="IM"]

wilcox.test(risk_window ~ IM,
       data = PAD_events_exposed_link_covars_geno2)    #test for IM vs NM only

#switching of clopi within 30 days of first event (death excluded)
PAD_events_exposed_link_covars_geno <- readRDS("PAD_events_exposed_link_covars_geno")

PAD_events_exposed_link_covars_geno$first_event <- pmin(PAD_events_exposed_link_covars_geno$date.y.y,
                                                        PAD_events_exposed_link_covars_geno$date.x.y,
                                                        PAD_events_exposed_link_covars_geno$date.y.x,
                                                        PAD_events_exposed_link_covars_geno$Activity_date1,
                                                        PAD_events_exposed_link_covars_geno$Activity_date2,
                                                        PAD_events_exposed_link_covars_geno$Activity_date3,
                                                        PAD_events_exposed_link_covars_geno$Activity_date4,
                                                        PAD_events_exposed_link_covars_geno$Activity_date5,
                                                        PAD_events_exposed_link_covars_geno$Activity_date6,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_1,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_2,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_3,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_4,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_5,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_6,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_7,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_8,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_9,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_10,
                                                        PAD_events_exposed_link_covars_geno$OPCS_Proc_DtOPCS_Proc_Dt_11,
                                                        na.rm = TRUE)
PAD_events_exposed_link_covars_geno [, first_event := as.Date(first_event)] #make date

#find second event (death excluded)
PAD_events_exposed_link_covars_geno[, second_event :=
                                      apply(.SD,1,function(x) {
                                        x<- sort (x[!is.na(x)])
                                        if(length(x) >=2) x[2] else as.Date(NA)
                                      }),
                                    .SDcols = c("OPCS_Proc_DtOPCS_Proc_Dt_1","OPCS_Proc_DtOPCS_Proc_Dt_2","OPCS_Proc_DtOPCS_Proc_Dt_3",
                                                "OPCS_Proc_DtOPCS_Proc_Dt_4","OPCS_Proc_DtOPCS_Proc_Dt_5","OPCS_Proc_DtOPCS_Proc_Dt_6",
                                                "OPCS_Proc_DtOPCS_Proc_Dt_7","OPCS_Proc_DtOPCS_Proc_Dt_8","OPCS_Proc_DtOPCS_Proc_Dt_9",
                                                "OPCS_Proc_DtOPCS_Proc_Dt_10","OPCS_Proc_DtOPCS_Proc_Dt_11", "date.y.x",
                                                "date.x.y", "date.y.y", "Activity_date1","Activity_date2","Activity_date3","Activity_date4",
                                                "Activity_date5","Activity_date6")
]

PAD_events_exposed_link_covars_geno [, second_event := as.Date(second_event)] #make date

PAD_events_exposed_link_covars_geno[, last_med := as.Date(last_med)] #make date

PAD_events_exposed_link_covars_geno [ last_med > first_event &
                                        last_med  <= first_event + 30, .N]

PAD_events_exposed_link_covars_geno [ last_med > second_event &
                                        last_med <= second_event + 30, .N]

PAD_events_exposed_link_covars_geno [ last_med > first_event &
                                        last_med <= first_event + 30 &
                                        PI=="other", .N]

PAD_events_exposed_link_covars_geno [ last_med > second_event &
                                        last_med <= second_event + 30 &
                                        PI=="other", .N]
