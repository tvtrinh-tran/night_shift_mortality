##########Main points of this script:
# R libraries
# Data loading
# Create variables for exposures, outcomes, covariates
# Check for Cox models' assumptions


#####Necessary R library and data loading for the Night Shift project
library(here)
library(tidyverse)
library(gtsummary)
library(survival)
library(tidymodels)
library(purrr)
library(rlang)
library(survminer)
library(gt)
#library(broom)
#library(broom.helpers)
#library(splines)
# library(Gmisc, quietly = TRUE)
# library(glue)
# library(htmlTable)
# library(grid)
# library(magrittr)

# library(EValue)
# library(meta)
# library(metafor)
# library(ggplot2)
# library(gridExtra)

directory="C:/Users/trant4/OneDrive - National Institutes of Health/Trinh/Postdoc NCI - Risk factors and thyroid cancer/Night_Shift_Mortality"

# label: data
night_shift_all =
  read.csv(paste(directory,"/Night_Shift_Mortality_Data/nightshift_cohort_20241008.csv",sep = ""),
           header = TRUE,
           stringsAsFactors = FALSE,
           fileEncoding = "latin1"
  )

night_shift_popu = night_shift_all|>select(!starts_with(c("COMB_","TYW_","PRIORITY", "LQ2_","QX3_SMK_CUR","QX3_SMK_QTY","QX3_MARSTAT")))|> replace_values_with_na()
night_shift_rad_dose = night_shift_all|>select(starts_with(c("COHORTID","COMB_","TYW_","PRIORITY", "LQ2_","QX3_SMK_CUR","QX3_SMK_QTY","QX3_MARSTAT")))  
str(night_shift_popu)
summary(night_shift_popu)

#Create labels for existing categorical variables
night_shift_popu$COHORTID = factor(night_shift_popu$COHORTID)
night_shift_popu$QX3_EDUC = factor(night_shift_popu$QX3_EDUC, labels=c("1-8 years (grade school)","9-12 years (high school)","Other (e.g. vocational)","2-year hospital rad tech program","1-4 years college","Graduate school"))|>set_unknown_level()
night_shift_popu$QX3_INCOME = factor(night_shift_popu$QX3_INCOME, labels=c("Less than $25,000","$25,000-$49,999","$50,000-$74,999","$75,000-$99,999","$100,000 or more"))|>set_unknown_level()
night_shift_popu$QX4_SMK_STATUS = factor(night_shift_popu$QX4_SMK_STATUS, labels=c("Never smoked","Former smoker","Current smoker"))|>set_unknown_level()
night_shift_popu$QX4_FH_LUNG_EV = factor(night_shift_popu$QX4_FH_LUNG_EV, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_FH_BREAST_EV = factor(night_shift_popu$QX4_FH_BREAST_EV, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_CHILD_EV = factor(night_shift_popu$QX4_CHILD_EV, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_MENOPAUSE_STATUS = factor(night_shift_popu$QX4_MENOPAUSE_STATUS, labels=c("No, still having periods","Yes","No, menstrual periods are irregular or using hormones","Never menstruated"))|>set_unknown_level()
night_shift_popu$QX4_HRT_EV = factor(night_shift_popu$QX4_HRT_EV, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_HRT_CUR = factor(night_shift_popu$QX4_HRT_CUR, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_EV_U30 = factor(night_shift_popu$QX4_NIGHT_PERM_EV_U30, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_EV_30 = factor(night_shift_popu$QX4_NIGHT_PERM_EV_30, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_EV_40 = factor(night_shift_popu$QX4_NIGHT_PERM_EV_40, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_EV_50 = factor(night_shift_popu$QX4_NIGHT_PERM_EV_50, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_DUR_U30 = factor(night_shift_popu$QX4_NIGHT_PERM_DUR_U30, labels=c("1","2-3","4-5","6-7","8+"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_DUR_30 = factor(night_shift_popu$QX4_NIGHT_PERM_DUR_30, labels=c("1","2-3","4-5","6-7","8+"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_DUR_40 = factor(night_shift_popu$QX4_NIGHT_PERM_DUR_40, labels=c("1","2-3","4-5","6-7","8+"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_DUR_50 = factor(night_shift_popu$QX4_NIGHT_PERM_DUR_50, labels=c("1","2-3","4-5","6-7","8+"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_TIMES_U30 = factor(night_shift_popu$QX4_NIGHT_PERM_TIMES_U30, labels=c("1-3","4-5","6-9","10-14","15-19","20+"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_TIMES_30 = factor(night_shift_popu$QX4_NIGHT_PERM_TIMES_30, labels=c("1-3","4-5","6-9","10-14","15-19","20+"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_TIMES_40 = factor(night_shift_popu$QX4_NIGHT_PERM_TIMES_40, labels=c("1-3","4-5","6-9","10-14","15-19","20+"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_PERM_TIMES_50 = factor(night_shift_popu$QX4_NIGHT_PERM_TIMES_50, labels=c("1-3","4-5","6-9","10-14","15-19","20+"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_ROTATING_EV_U30 = factor(night_shift_popu$QX4_NIGHT_ROTATING_EV_U30, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_ROTATING_EV_30 = factor(night_shift_popu$QX4_NIGHT_ROTATING_EV_30, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_ROTATING_EV_40 = factor(night_shift_popu$QX4_NIGHT_ROTATING_EV_40, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_NIGHT_ROTATING_EV_50 = factor(night_shift_popu$QX4_NIGHT_ROTATING_EV_50, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_EXERCISE_WALK = factor(night_shift_popu$QX4_EXERCISE_WALK, labels=c("None","0.5","1","1.5","2-3","4-6","7-10","11+"))|>set_unknown_level()
night_shift_popu$QX4_EXERCISE_MOD = factor(night_shift_popu$QX4_EXERCISE_MOD, labels=c("None","0.5","1","1.5","2-3","4-6","7-10","11+"))|>set_unknown_level()
night_shift_popu$QX4_EXERCISE_STREN = factor(night_shift_popu$QX4_EXERCISE_STREN, labels=c("None","0.5","1","1.5","2-3","4-6","7-10","11+"))|>set_unknown_level()
night_shift_popu$QX4_EXERCISE_WEIGHT = factor(night_shift_popu$QX4_EXERCISE_WEIGHT, labels=c("None","0.5","1","1.5","2-3","4-6","7-10","11+"))|>set_unknown_level()
night_shift_popu$QX4_SLEEP_HOURS_WEEKDAY = factor(night_shift_popu$QX4_SLEEP_HOURS_WEEKDAY, labels=c("1-4","5","6","7","8","9","10+"))|>set_unknown_level()
night_shift_popu$QX4_SLEEP_HOURS_WEEKEND = factor(night_shift_popu$QX4_SLEEP_HOURS_WEEKEND, labels=c("1-4","5","6","7","8","9","10+"))|>set_unknown_level()
night_shift_popu$QX4_SLEEP_ADVERSE = factor(night_shift_popu$QX4_SLEEP_ADVERSE, labels=c("None","1","2-3","4-5","6-7","8+"))|>set_unknown_level()
night_shift_popu$QX4_SLEEP_LIGHT = factor(night_shift_popu$QX4_SLEEP_LIGHT, labels=c("Bright light","Some light","Completely dark"))|>set_unknown_level()
night_shift_popu$QX4_SLEEP_MIDNIGHT_EV = factor(night_shift_popu$QX4_SLEEP_MIDNIGHT_EV, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$QX4_SLEEP_MIDNIGHT_BEDTIME = factor(night_shift_popu$QX4_SLEEP_MIDNIGHT_BEDTIME, labels=c("12:00-1:00 am","1:00-2:00 am","2:00-3:00 am","After 3:00 am"))|>set_unknown_level()
night_shift_popu$QX4_SLEEP_MIDNIGHT_TIMES = factor(night_shift_popu$QX4_SLEEP_MIDNIGHT_TIMES, labels=c("1-4","5-8","9-15","16+"))|>set_unknown_level()
night_shift_popu$QX4_PERSON_TYPE = factor(night_shift_popu$QX4_PERSON_TYPE, labels=c("Morning person","Evening person","Neither","Both"))|>set_unknown_level()
night_shift_popu$PREV_CANCER = factor(night_shift_popu$PREV_CANCER, labels=c("No","Yes","Unknown"))|>set_unknown_level()
night_shift_popu$PREV_LUNG = factor(night_shift_popu$PREV_LUNG, labels=c("No","Yes","Unknown"))|>set_unknown_level()
night_shift_popu$PREV_BREAST = factor(night_shift_popu$PREV_BREAST, labels=c("No","Yes","Unknown","n/a"))|>set_unknown_level()
night_shift_popu$PREV_CVD = factor(night_shift_popu$PREV_CVD, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$FH_LUNG = factor(night_shift_popu$FH_LUNG, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$FH_BREAST = factor(night_shift_popu$FH_BREAST, labels=c("No","Yes"))|>set_unknown_level()
night_shift_popu$EXIT_STATUS = factor(night_shift_popu$EXIT_STATUS, labels=c("Alive", "Dead"))|>set_unknown_level()
night_shift_popu$RACE = factor(night_shift_popu$RACE, labels=c("white","black","asian/pi","ai/an","other","unknown"))

#convert continuous variable from character to numeric:
night_shift_popu[sapply(night_shift_popu, is.character)] = lapply(night_shift_popu[sapply(night_shift_popu, is.character)], as.numeric)


#####Create variables for exposures, outcomes, and covariates
night_shift_popu = night_shift_popu |>
  mutate(
    
    #####EXPOSURES
    
    #Ever night shift: permanent
    night_perm_ev_all = factored_case_when(QX4_NIGHT_PERM_EV_U30=="No" & QX4_NIGHT_PERM_EV_30=="No" & QX4_NIGHT_PERM_EV_40=="No" & QX4_NIGHT_PERM_EV_50=="No"~"No",
                                           QX4_NIGHT_PERM_EV_U30=="Yes"|QX4_NIGHT_PERM_EV_30=="Yes"|QX4_NIGHT_PERM_EV_40=="Yes"|QX4_NIGHT_PERM_EV_50=="Yes"~"Yes",
                                           TRUE ~ "Unknown status"),
    
    #Ever night shift: rotating
    night_rotating_ev_all = factored_case_when(QX4_NIGHT_ROTATING_EV_U30=="No" & QX4_NIGHT_ROTATING_EV_30=="No" & QX4_NIGHT_ROTATING_EV_40=="No" & QX4_NIGHT_ROTATING_EV_50=="No"~"No",
                                               QX4_NIGHT_ROTATING_EV_U30=="Yes"|QX4_NIGHT_ROTATING_EV_30=="Yes"|QX4_NIGHT_ROTATING_EV_40=="Yes"|QX4_NIGHT_ROTATING_EV_50=="Yes"~"Yes",
                                               TRUE ~ "Unknown status"),
    
    #Ever night shift: all
    night_ev_all = factored_case_when(night_perm_ev_all=="No" & night_rotating_ev_all=="No"~"Never worked permanent or rotating night shifts",
                                      night_perm_ev_all=="Yes" ~ "Permanent night shift work",
                                      night_perm_ev_all=="No" & night_rotating_ev_all=="Yes"~"Rotating night shift work",
                                      TRUE ~ "Unknown status"),
    night_ev_all_2cat = factor(night_ev_all,
                               labels = c("Known night shift work status",
                                          "Known night shift work status",
                                          "Known night shift work status",
                                          "Unknown night shift work status")),
    
    
    #Age first start night shift work:

    night_perm_age_start_5cat = factored_case_when(night_ev_all=="Never worked permanent or rotating night shifts"~"Never worked permanent or rotating night shifts",
                                              QX4_NIGHT_PERM_EV_U30=="Yes"~"Under 30",
                                              QX4_NIGHT_PERM_EV_U30=="No"&QX4_NIGHT_PERM_EV_30=="Yes"~"Age 30-39",
                                              QX4_NIGHT_PERM_EV_U30=="No"&QX4_NIGHT_PERM_EV_30=="No"&
                                                (QX4_NIGHT_PERM_EV_40=="Yes"|(QX4_NIGHT_PERM_EV_40=="No"&QX4_NIGHT_PERM_EV_50=="Yes"))~"Age 40-59",
                                              night_ev_all=="Rotating night shift work" ~ "Rotating night shift work",
                                              TRUE ~ "Unknown age start"),
    night_perm_age_start_5cat = factor(night_perm_age_start_5cat,levels=c("Never worked permanent or rotating night shifts",
                                                                          "Under 30",
                                                                          "Age 30-39",
                                                                          "Age 40-59",
                                                                          "Unknown age start",
                                                                          "Rotating night shift work")),
    night_perm_age_start_4cat = factor(night_perm_age_start_5cat,
                                  labels = c("Never worked permanent or rotating night shifts",
                                             "Under 30",
                                             "After 30",
                                             "After 30",
                                             "Unknown age start",
                                             "Rotating night shift work")),
    night_perm_age_start_5cat_nightworkers = factor(night_perm_age_start_5cat, exclude = "Never worked permanent or rotating night shifts"),
    night_perm_age_start_4cat_nightworkers = factor(night_perm_age_start_4cat, exclude = "Never worked permanent or rotating night shifts"),
    
    #Duration of night shift (years): total number of years employed in permanent night-shift work
    night_perm_dur_u30 = case_when(QX4_NIGHT_PERM_EV_U30=="No" ~ 0,
                                   QX4_NIGHT_PERM_EV_U30=="Unknown status" | QX4_NIGHT_PERM_DUR_U30=="Unknown status" ~ NA,
                                   QX4_NIGHT_PERM_DUR_U30=="1" ~ 1,
                                   QX4_NIGHT_PERM_DUR_U30=="2-3" ~ 2.5,
                                   QX4_NIGHT_PERM_DUR_U30=="4-5" ~ 4.5,
                                   QX4_NIGHT_PERM_DUR_U30=="6-7"~6.5,
                                   QX4_NIGHT_PERM_DUR_U30=="8+"~9),
    
    night_perm_dur_30 = case_when(QX4_NIGHT_PERM_EV_30=="No" ~ 0,
                                  QX4_NIGHT_PERM_EV_30=="Unknown status" | QX4_NIGHT_PERM_DUR_30=="Unknown status" ~ NA,
                                  QX4_NIGHT_PERM_DUR_30=="1" ~ 1,
                                  QX4_NIGHT_PERM_DUR_30=="2-3" ~ 2.5,
                                  QX4_NIGHT_PERM_DUR_30=="4-5" ~ 4.5,
                                  QX4_NIGHT_PERM_DUR_30=="6-7"~6.5,
                                  QX4_NIGHT_PERM_DUR_30=="8+"~9),
    
    night_perm_dur_40 = case_when(QX4_NIGHT_PERM_EV_40=="No" ~ 0,
                                  QX4_NIGHT_PERM_EV_40=="Unknown status" | QX4_NIGHT_PERM_DUR_40=="Unknown status" ~ NA,
                                  QX4_NIGHT_PERM_DUR_40=="1" ~ 1,
                                  QX4_NIGHT_PERM_DUR_40=="2-3" ~ 2.5,
                                  QX4_NIGHT_PERM_DUR_40=="4-5" ~ 4.5,
                                  QX4_NIGHT_PERM_DUR_40=="6-7"~6.5,
                                  QX4_NIGHT_PERM_DUR_40=="8+"~9),
    
    night_perm_dur_50 = case_when(QX4_NIGHT_PERM_EV_50=="No" ~ 0,
                                  QX4_NIGHT_PERM_EV_50=="Unknown status" | QX4_NIGHT_PERM_DUR_50=="Unknown status" ~ NA,
                                  QX4_NIGHT_PERM_DUR_50=="1" ~ 1,
                                  QX4_NIGHT_PERM_DUR_50=="2-3" ~ 2.5,
                                  QX4_NIGHT_PERM_DUR_50=="4-5" ~ 4.5,
                                  QX4_NIGHT_PERM_DUR_50=="6-7"~6.5,
                                  QX4_NIGHT_PERM_DUR_50=="8+"~9),
    
    night_perm_dur_all_cont = night_perm_dur_u30 + night_perm_dur_30 + night_perm_dur_40 + night_perm_dur_50,#will generate NA value if missing any value for dur_u30, dur_30, dur_40, dur_50
    
    night_perm_dur_all_6cat = factored_case_when(night_ev_all=="Never worked permanent or rotating night shifts"~"Never worked permanent or rotating night shifts",
                                                night_perm_ev_all=="Yes"&5>night_perm_dur_all_cont & night_perm_dur_all_cont>0~"Less than 5 years",
                                                night_perm_ev_all=="Yes"&10>night_perm_dur_all_cont & night_perm_dur_all_cont>=5~"5 to 9 years",
                                                night_perm_ev_all=="Yes"&night_perm_dur_all_cont & night_perm_dur_all_cont>=10~"10 years or more",
                                                night_perm_ev_all=="Yes"&is.na(night_perm_dur_all_cont)~"Yes with unknown duration",
                                                TRUE~"Unknown status"),
    
    night_perm_dur_all_5cat = factor(night_perm_dur_all_6cat,
                                     labels = c("Never worked permanent or rotating night shifts",
                                                "Less than 5 years",
                                                "More than 5 years",
                                                "More than 5 years",
                                                "Yes with unknown duration",
                                                "Unknown status")),
    night_perm_dur_all_6cat_nightworkers = factor(night_perm_dur_all_6cat, exclude = c("Never worked permanent or rotating night shifts","Yes with unknown duration","Unknown status")),
    night_perm_dur_all_5cat_nightworkers = factor(night_perm_dur_all_5cat, exclude = c("Never worked permanent or rotating night shifts","Yes with unknown duration","Unknown status")),
    
    #Intensity of night shift: mean number of night shifts per MONTH averaged over the duration of all jobs involving permanent night shifts
    night_perm_intensity_u30 = case_when(QX4_NIGHT_PERM_EV_U30=="No" ~ 0,
                                         QX4_NIGHT_PERM_EV_U30=="Unknown status" | QX4_NIGHT_PERM_TIMES_U30=="Unknown status" ~ NA,
                                         QX4_NIGHT_PERM_TIMES_U30=="1-3" ~ 2,
                                         QX4_NIGHT_PERM_TIMES_U30=="4-5" ~ 4.5,
                                         QX4_NIGHT_PERM_TIMES_U30=="6-9" ~ 7.5,
                                         QX4_NIGHT_PERM_TIMES_U30=="10-14"~12,
                                         QX4_NIGHT_PERM_TIMES_U30=="15-19"~17,
                                         QX4_NIGHT_PERM_TIMES_U30=="20+"~25),
    
    night_perm_intensity_30 = case_when(QX4_NIGHT_PERM_EV_30=="No" ~ 0,
                                        QX4_NIGHT_PERM_EV_30=="Unknown status" | QX4_NIGHT_PERM_TIMES_30=="Unknown status" ~ NA,
                                        QX4_NIGHT_PERM_TIMES_30=="1-3" ~ 2,
                                        QX4_NIGHT_PERM_TIMES_30=="4-5" ~ 4.5,
                                        QX4_NIGHT_PERM_TIMES_30=="6-9" ~ 7.5,
                                        QX4_NIGHT_PERM_TIMES_30=="10-14"~12,
                                        QX4_NIGHT_PERM_TIMES_30=="15-19"~17,
                                        QX4_NIGHT_PERM_TIMES_30=="20+"~25),
    
    night_perm_intensity_40 = case_when(QX4_NIGHT_PERM_EV_40=="No" ~ 0,
                                        QX4_NIGHT_PERM_EV_40=="Unknown status" | QX4_NIGHT_PERM_TIMES_40=="Unknown status" ~ NA,
                                        QX4_NIGHT_PERM_TIMES_40=="1-3" ~ 2,
                                        QX4_NIGHT_PERM_TIMES_40=="4-5" ~ 4.5,
                                        QX4_NIGHT_PERM_TIMES_40=="6-9" ~ 7.5,
                                        QX4_NIGHT_PERM_TIMES_40=="10-14"~12,
                                        QX4_NIGHT_PERM_TIMES_40=="15-19"~17,
                                        QX4_NIGHT_PERM_TIMES_40=="20+"~25),
    
    night_perm_intensity_50 = case_when(QX4_NIGHT_PERM_EV_50=="No" ~ 0,
                                        QX4_NIGHT_PERM_EV_50=="Unknown status" | QX4_NIGHT_PERM_TIMES_50=="Unknown status" ~ NA,
                                        QX4_NIGHT_PERM_TIMES_50=="1-3" ~ 2,
                                        QX4_NIGHT_PERM_TIMES_50=="4-5" ~ 4.5,
                                        QX4_NIGHT_PERM_TIMES_50=="6-9" ~ 7.5,
                                        QX4_NIGHT_PERM_TIMES_50=="10-14"~12,
                                        QX4_NIGHT_PERM_TIMES_50=="15-19"~17,
                                        QX4_NIGHT_PERM_TIMES_50=="20+"~25),
    
    night_perm_cumulative_all_cont = case_when(night_ev_all=="Never worked permanent or rotating night shifts"~0,
                                               .default = (night_perm_intensity_u30*night_perm_dur_u30 + night_perm_intensity_30*night_perm_dur_30 + 
                                                             night_perm_intensity_40*night_perm_dur_40 + night_perm_intensity_50*night_perm_dur_50)*12),
    #multiple with 12 months to have the total number of shifts
    
    night_perm_cumulative_all_6cat = factored_case_when(night_ev_all=="Never worked permanent or rotating night shifts"~"Never worked permanent or rotating night shifts",
                                                       night_perm_ev_all=="Yes"&500>night_perm_cumulative_all_cont & night_perm_cumulative_all_cont>0~"Less than 500 shifts",
                                                       night_perm_ev_all=="Yes"&1500>night_perm_cumulative_all_cont & night_perm_cumulative_all_cont>=500~"500 to 1499 shifts",
                                                       night_perm_ev_all=="Yes"&night_perm_cumulative_all_cont>=1500~"1500 shifts or more",
                                                       night_perm_ev_all=="Yes"&is.na(night_perm_cumulative_all_cont)~"Yes with unknown cumulative exposure",
                                                       TRUE~"Unknown status"),
    night_perm_cumulative_all_5cat = factor(night_perm_cumulative_all_6cat,
                                            labels = c("Never worked permanent or rotating night shifts",
                                                       "Less than 500 shifts",
                                                       "More than 500 shifts",
                                                       "More than 500 shifts",
                                                       "Yes with unknown cumulative exposure",
                                                       "Unknown status")),
    
    night_perm_cumulative_all_6cat_nightworkers = factor(night_perm_cumulative_all_6cat, exclude = c("Never worked permanent or rotating night shifts","Yes with unknown cumulative exposure","Unknown status")),
    night_perm_cumulative_all_5cat_nightworkers = factor(night_perm_cumulative_all_5cat, exclude = c("Never worked permanent or rotating night shifts","Yes with unknown cumulative exposure","Unknown status")),
    
    night_perm_intensity_all_cont = night_perm_cumulative_all_cont/(12*night_perm_dur_all_cont),
    
    night_perm_intensity_all_5cat = factored_case_when(night_ev_all=="Never worked permanent or rotating night shifts"~"Never worked permanent or rotating night shifts",
                                                      night_perm_ev_all=="Yes"&15>night_perm_intensity_all_cont & night_perm_intensity_all_cont>0~"Less than 15 shifts per month",
                                                      night_perm_ev_all=="Yes"&night_perm_intensity_all_cont>=15~"15 shifts per month or more",
                                                      night_perm_ev_all=="Yes"&is.na(night_perm_intensity_all_cont)~"Yes with unknown intensity",
                                                      TRUE~"Unknown status"),
    night_perm_intensity_all_3cat = factor(night_perm_intensity_all_5cat, exclude = c("Yes with unknown intensity","Unknown status")),
    night_perm_intensity_all_5cat_nightworkers = factor(night_perm_intensity_all_5cat, exclude = c("Never worked permanent or rotating night shifts","Yes with unknown intensity","Unknown status")),
    
    
    ######Sleep patterns during one year before the completion of the 2012-2014 questionnaire
    sleep_hours_weekday_cont = case_when(QX4_SLEEP_HOURS_WEEKDAY=="1-4" ~ 2.5  ,
                                         QX4_SLEEP_HOURS_WEEKDAY=="5" ~ 5,
                                         QX4_SLEEP_HOURS_WEEKDAY=="6" ~ 6,
                                         QX4_SLEEP_HOURS_WEEKDAY=="7" ~ 7,
                                         QX4_SLEEP_HOURS_WEEKDAY=="8" ~ 8,
                                         QX4_SLEEP_HOURS_WEEKDAY=="9" ~ 9,
                                         QX4_SLEEP_HOURS_WEEKDAY=="10+" ~ 10,
                                         QX4_SLEEP_HOURS_WEEKDAY=="Unknown status"~NA),
    
    sleep_hours_weekend_cont = case_when(QX4_SLEEP_HOURS_WEEKEND=="1-4" ~ 2.5  ,
                                         QX4_SLEEP_HOURS_WEEKEND=="5" ~ 5,
                                         QX4_SLEEP_HOURS_WEEKEND=="6" ~ 6,
                                         QX4_SLEEP_HOURS_WEEKEND=="7" ~ 7,
                                         QX4_SLEEP_HOURS_WEEKEND=="8" ~ 8,
                                         QX4_SLEEP_HOURS_WEEKEND=="9" ~ 9,
                                         QX4_SLEEP_HOURS_WEEKEND=="10+" ~ 10,
                                         QX4_SLEEP_HOURS_WEEKEND=="Unknown status"~NA),
    
    sleep_hours_all_cont = (sleep_hours_weekday_cont*5 + sleep_hours_weekend_cont*2)/7, #missing data in either sleep_hours_weekday_cont or sleep_hours_weekend_cont will result in missing data here
    sleep_hours_all_cat  = factored_case_when(6>=sleep_hours_all_cont & sleep_hours_all_cont>0~"6 hours or less",
                                              (sleep_hours_all_cont==7 | sleep_hours_all_cont==8)~"7-8 hours",
                                              sleep_hours_all_cont>=9~"9 hours or more",
                                              TRUE ~ "Unknown status"),
    
    time_per_month_midnight = factored_case_when(QX4_SLEEP_MIDNIGHT_EV=="No"~"None",
                                                 QX4_SLEEP_MIDNIGHT_EV=="Yes"&QX4_SLEEP_MIDNIGHT_TIMES %in% c("1-4","5-8","9-15")~"1-15 times",
                                                 QX4_SLEEP_MIDNIGHT_EV=="Yes"&QX4_SLEEP_MIDNIGHT_TIMES == "16+"~"16 times or more",
                                                 TRUE ~ "Unknown status"),
    
    
    morning_chronotype = factor(QX4_PERSON_TYPE, levels=c("Morning person","Evening person","Neither","Both","Unknown" ),
                                labels=c("Morning person","Evening person","Neither","Both","Unknown status")),
    
#####Outcomes
    MOR_CANCER_MODIF = case_when(MOR_CANCER==2~0,.default = MOR_CANCER),
    MOR_CVD_MODIF = case_when(MOR_CVD==2~0,.default = MOR_CVD),
    MOR_LUNG_MODIF = case_when(MOR_LUNG==2~0,.default = MOR_LUNG),
    MOR_BREAST_MODIF = case_when(MOR_BREAST==2~0,.default = MOR_BREAST),
    
###### Covariates
    sex = factor(SEX, labels=c("Male", "Female")),
    race_modif = factor(RACE, levels=c("white","black","asian/pi","ai/an","other","unknown"),
                        labels = c("White","Black","Other/Unknown status","Other/Unknown status","Other/Unknown status","Other/Unknown status")),
    firstwork_year_cat = factored_case_when(1963>FIRSTWORK_YEAR &  FIRSTWORK_YEAR>0~"Before 1963",
                                            1969>FIRSTWORK_YEAR &  FIRSTWORK_YEAR>=1963~"1963–1969",
                                            1974>FIRSTWORK_YEAR &  FIRSTWORK_YEAR>1970~"1970–1974",
                                            FIRSTWORK_YEAR>=1975~"After 1975",
                                            TRUE ~ "Unknown status"),
    bmi_cont = (QX4_WEIGHT*0.453592)/(QX4_HEIGHT*0.0254)**2,#convert weight from pounds to kg and height from inches to meter
    bmi_cat = factored_case_when(25>bmi_cont & bmi_cont>0~"<25",
                                     30>bmi_cont & bmi_cont>=25~"25-29.9",
                                 bmi_cont>=30~"30 or higher",
                                 TRUE ~ "Unknown status"),
    smk_status = factor(QX4_SMK_STATUS, labels=c("Never smoked","Former smoker","Current smoker","Unknown status")),
    
    #calculate MET-hours per week, follow the steps in this paper (https://ascopubs.org/doi/10.1200/JCO.19.02407#suppl-html-sec-1)
    exercise_walk_hours_perweek = case_when(QX4_EXERCISE_WALK=="None" | QX4_EXERCISE_WALK=="Unknown" ~ 0  ,
                              QX4_EXERCISE_WALK=="0.5" ~ 0.5,
                              QX4_EXERCISE_WALK=="1" ~ 1,
                              QX4_EXERCISE_WALK=="1.5" ~ 1.5,
                              QX4_EXERCISE_WALK=="2-3" ~ 2.5,
                              QX4_EXERCISE_WALK=="4-6"~5,
                              QX4_EXERCISE_WALK=="7-10"~8.5,
                              QX4_EXERCISE_WALK=="11+"~11),
    exercise_mod_hours_perweek = case_when(QX4_EXERCISE_MOD=="None" | QX4_EXERCISE_WALK=="Unknown" ~ 0  ,
                              QX4_EXERCISE_MOD=="0.5" ~ 0.5,
                              QX4_EXERCISE_MOD=="1" ~ 1,
                              QX4_EXERCISE_MOD=="1.5" ~ 1.5,
                              QX4_EXERCISE_MOD=="2-3" ~ 2.5,
                              QX4_EXERCISE_MOD=="4-6"~5,
                              QX4_EXERCISE_MOD=="7-10"~8.5,
                              QX4_EXERCISE_MOD=="11+"~11),
    exercise_stren_hours_perweek = case_when(QX4_EXERCISE_STREN=="None" | QX4_EXERCISE_WALK=="Unknown" ~ 0  ,
                              QX4_EXERCISE_STREN=="0.5" ~ 0.5,
                              QX4_EXERCISE_STREN=="1" ~ 1,
                              QX4_EXERCISE_STREN=="1.5" ~ 1.5,
                              QX4_EXERCISE_STREN=="2-3" ~ 2.5,
                              QX4_EXERCISE_STREN=="4-6"~5,
                              QX4_EXERCISE_STREN=="7-10"~8.5,
                              QX4_EXERCISE_STREN=="11+"~11),
    exercise_weight_hours_perweek = case_when(QX4_EXERCISE_WEIGHT=="None" | QX4_EXERCISE_WALK=="Unknown" ~ 0  ,
                              QX4_EXERCISE_WEIGHT=="0.5" ~ 0.5,
                              QX4_EXERCISE_WEIGHT=="1" ~ 1,
                              QX4_EXERCISE_WEIGHT=="1.5" ~ 1.5,
                              QX4_EXERCISE_WEIGHT=="2-3" ~ 2.5,
                              QX4_EXERCISE_WEIGHT=="4-6"~5,
                              QX4_EXERCISE_WEIGHT=="7-10"~8.5,
                              QX4_EXERCISE_WEIGHT=="11+"~11),
    #https://www.physio-pedia.com/images/c/c7/Quidelines_for_interpreting_the_IPAQ.pdf convert MET-minutes/week to MET-hours/week, consider weight training and resistance exercises as vigorous exercise
    METhours_perweek = (3.3*exercise_walk_hours_perweek + 4.0*exercise_mod_hours_perweek + 8.0*exercise_stren_hours_perweek + 8.0*exercise_weight_hours_perweek),
    
    phy_act = factored_case_when(METhours_perweek==0~"0 MET-hours/week",
                                 15>METhours_perweek & METhours_perweek>0~"0.1-14.9 MET-hours/week",
                                 METhours_perweek>=15~"≥15 MET-hours/week",
                                 TRUE ~ "Unknown status"),
    edu_level = factor(night_shift_popu$QX3_EDUC, levels = c("1-8 years (grade school)",
                                                             "9-12 years (high school)",
                                                             "2-year hospital rad tech program",
                                                             "1-4 years college",
                                                             "Graduate school",
                                                             "Other (e.g. vocational)",
                                                             "Unknown"),
                       labels=c("Radiotechnology technologist program, high school or lower",
                                "Radiotechnology technologist program, high school or lower",
                                "Radiotechnology technologist program, high school or lower",
                                "1-4 years college",
                                "Graduate school",
                                "Other/Unknown status",
                                "Other/Unknown status")),
    income = factor(QX3_INCOME, levels = c("Less than $25,000",
                                           "$25,000-$49,999",
                                           "$50,000-$74,999",
                                           "$75,000-$99,999",
                                           "$100,000 or more",
                                           "Unknown"),
                    labels=c("Less than $50,000",
                             "Less than $50,000",
                             "$50,000-$74,999",
                             "$75,000-$99,999",
                             "$100,000 or more",
                             "Unknown status")))|>
  mutate(night_perm_dur_u30=NULL,
         night_perm_dur_30=NULL,
         night_perm_dur_40=NULL,
         night_perm_dur_50=NULL,
         night_perm_intensity_u30=NULL,
         night_perm_intensity_30=NULL,
         night_perm_intensity_40=NULL,
         night_perm_intensity_50=NULL,
         sleep_hours_weekday_cont=NULL,
         sleep_hours_weekend_cont=NULL,
         exercise_walk_hours_perweek=NULL,
         exercise_mod_hours_perweek=NULL,
         exercise_stren_hours_perweek=NULL,
         exercise_weight_hours_perweek=NULL,
         METhours_perweek=NULL)
night_shift_popu_before_exclude_unknown = night_shift_popu

night_shift_popu=night_shift_popu[night_shift_popu$night_ev_all!="Unknown status",]|>
  mutate(night_ev_all = droplevels(night_ev_all, exclude = "Unknown status"))
night_shift_popu_cvd=night_shift_popu[night_shift_popu$PREV_CVD!="Yes",]
night_shift_popu_cancer=night_shift_popu[night_shift_popu$PREV_CANCER!="Yes",]
night_shift_popu_lung=night_shift_popu[night_shift_popu$PREV_LUNG!="Yes",]
night_shift_popu_breast=night_shift_popu[night_shift_popu$sex=="Female"&night_shift_popu$PREV_BREAST!="Yes",]

#Test Cox proportional hazard assumption for all outcomes: https://bookdown.org/rwnahhas/RMPH/survival-phassumption.html
#No violation found
# night_shift_popu_test=night_shift_popu
# categories = c("night_ev_all","night_perm_ev_all","night_rotating_ev_all","night_perm_dur_all_cat",
#                "night_perm_intensity_all_cat","night_perm_cumulative_all_cat")
# 
# cox_models_assumption(data=night_shift_popu_test,categories=categories,response_formula = "Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)")
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)~night_perm_dur_all_cont+ tt(night_perm_dur_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)~night_perm_intensity_all_cont+ tt(night_perm_intensity_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)~night_perm_cumulative_all_cont+ tt(night_perm_cumulative_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# 
# cox_models_assumption(data=night_shift_popu_test,categories=categories,response_formula = "Surv(ENTRY_AGE, EXIT_AGE, MOR_CANCER_MODIF)")
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_CANCER)~night_perm_dur_all_cont+ tt(night_perm_dur_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_CANCER)~night_perm_intensity_all_cont+ tt(night_perm_intensity_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_CANCER)~night_perm_cumulative_all_cont+ tt(night_perm_cumulative_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# 
# cox_models_assumption(data=night_shift_popu_test,categories=categories,response_formula = "Surv(ENTRY_AGE, EXIT_AGE, MOR_CVD_MODIF)")
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_CVD)~night_perm_dur_all_cont+ tt(night_perm_dur_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_CVD)~night_perm_intensity_all_cont+ tt(night_perm_intensity_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_CVD)~night_perm_cumulative_all_cont+ tt(night_perm_cumulative_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)

# cox_models_assumption(data=night_shift_popu_test[night_shift_popu_test$sex=="female",],categories=categories,response_formula = "Surv(ENTRY_AGE, EXIT_AGE, MOR_BREAST_MODIF)")
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_BREAST)~night_perm_dur_all_cont+ tt(night_perm_dur_all_cont), data = night_shift_popu_test[night_shift_popu_test$sex=="female",], tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_BREAST)~night_perm_intensity_all_cont+ tt(night_perm_intensity_all_cont), data = night_shift_popu_test[night_shift_popu_test$sex=="female",], tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_BREAST)~night_perm_cumulative_all_cont+ tt(night_perm_cumulative_all_cont), data = night_shift_popu_test[night_shift_popu_test$sex=="female",], tt = function(x, t, ...) x * t)

# cox_models_assumption(data=night_shift_popu_test,categories=categories,response_formula = "Surv(ENTRY_AGE, EXIT_AGE, MOR_LUNG_MODIF)")
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_LUNG)~night_perm_dur_all_cont+ tt(night_perm_dur_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_LUNG)~night_perm_intensity_all_cont+ tt(night_perm_intensity_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_LUNG)~night_perm_cumulative_all_cont+ tt(night_perm_cumulative_all_cont), data = night_shift_popu_test, tt = function(x, t, ...) x * t)

#Test for log-linear assumption
#No evidence of non-linearity found
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE) ~ night_perm_dur_all_cont + log(night_perm_dur_all_cont) + sqrt(night_perm_dur_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_dur_all_cont)>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE) ~ night_perm_intensity_all_cont + log(night_perm_intensity_all_cont) + sqrt(night_perm_intensity_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_intensity_all_cont)>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE) ~ night_perm_cumulative_all_cont + log(night_perm_cumulative_all_cont) + sqrt(night_perm_cumulative_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_cumulative_all_cont)>0,])
# 
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_CANCER_MODIF) ~ night_perm_dur_all_cont + log(night_perm_dur_all_cont) + sqrt(night_perm_dur_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_dur_all_cont)>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_CANCER_MODIF) ~ night_perm_intensity_all_cont + log(night_perm_intensity_all_cont) + sqrt(night_perm_intensity_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_intensity_all_cont)>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_CANCER_MODIF) ~ night_perm_cumulative_all_cont + log(night_perm_cumulative_all_cont) + sqrt(night_perm_cumulative_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_cumulative_all_cont)>0,])
# 
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_CVD_MODIF) ~ night_perm_dur_all_cont + log(night_perm_dur_all_cont) + sqrt(night_perm_dur_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_dur_all_cont)>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_CVD_MODIF) ~ night_perm_intensity_all_cont + log(night_perm_intensity_all_cont) + sqrt(night_perm_intensity_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_intensity_all_cont)>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_CVD_MODIF) ~ night_perm_cumulative_all_cont + log(night_perm_cumulative_all_cont) + sqrt(night_perm_cumulative_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_cumulative_all_cont)>0,])
# 
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_LUNG_MODIF) ~ night_perm_dur_all_cont + log(night_perm_dur_all_cont) + sqrt(night_perm_dur_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_dur_all_cont)>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_LUNG_MODIF) ~ night_perm_intensity_all_cont + log(night_perm_intensity_all_cont) + sqrt(night_perm_intensity_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_intensity_all_cont)>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_LUNG_MODIF) ~ night_perm_cumulative_all_cont + log(night_perm_cumulative_all_cont) + sqrt(night_perm_cumulative_all_cont), 
#                            data = night_shift_popu_test[(night_shift_popu_test$night_perm_cumulative_all_cont)>0,])
# 
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_BREAST_MODIF) ~ night_perm_dur_all_cont + log(night_perm_dur_all_cont) + sqrt(night_perm_dur_all_cont), 
#                 data = night_shift_popu_test[night_shift_popu_test$sex=="female"&night_shift_popu_test$night_perm_cumulative_all_cont>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_BREAST_MODIF) ~ night_perm_intensity_all_cont + log(night_perm_intensity_all_cont) + sqrt(night_perm_intensity_all_cont), 
#                 data = night_shift_popu_test[night_shift_popu_test$sex=="female"&night_shift_popu_test$night_perm_cumulative_all_cont>0,])
# survminer::ggcoxfunctional(Surv(ENTRY_AGE, EXIT_AGE, MOR_BREAST_MODIF) ~ night_perm_cumulative_all_cont + log(night_perm_cumulative_all_cont) + sqrt(night_perm_cumulative_all_cont), 
#                 data = night_shift_popu_test[night_shift_popu_test$sex=="female"&night_shift_popu_test$night_perm_cumulative_all_cont>0,])
# 
