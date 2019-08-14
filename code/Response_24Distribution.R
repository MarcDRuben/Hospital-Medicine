
#Characterize time-of-day pattern of response to hydralazine treatment for all inpatients who received the drug at CCHMC from 2011-17

#Required R packages
library(tidyverse)
library(lubridate)
library(cosinor2)
library(broom)

#Load patient data. Create unique identifiers for each patient-order and patient-dose (note: an order may have multiple doses associated)
dat <- read_csv("./data/CCHMC_EMR_Hydralazine_Deidentified.csv") %>%
  mutate(PtOrdID = paste(PatID, ORDER_START_TIME, sep = "_")) %>%
  mutate(PtDosID = paste(PatID, MEDS_TAKEN, sep = "_"))

###FIX THIS###
filter(!is.na(DOSE_MG_KG)) %>% #optional remove doses without information on concentration

#Vector specifying time bins for analysis; 2 h intervals used
times = c("12AM-2AM", "2AM-4AM", "4AM-6AM", "6AM-8AM", "8AM-10AM", "10AM-12PM", "12PM-2PM", "2PM-4PM", "4PM-6PM", "6PM-8PM", "8PM-10PM", "10PM-12AM") 

#Clean & transform data for analysis
dat_ct <- dat %>%
  #optional filter on ROUTE
  filter(ROUTE == "Intravenous") %>%
  #convert all doses to mg/kg
  mutate(DOSE_MG_KG = ifelse(UNIT == "mg", round(HV_DISCRETE_DOSE/WEIGHT_KG, digits = 2), HV_DISCRETE_DOSE)) %>%
  #optional remove doses without concentration 
  filter(!is.na(DOSE_MG_KG)) %>%
  #add variable that specifies order and taken TIMES but NOT date.
  mutate_at(vars(ORDER_START_NODATE = ORDER_START_TIME, MEDS_TAKEN_NODATE = MEDS_TAKEN),
            .funs = funs(format(ymd_hms(.), "%H:%M:%S"))) %>% 
  mutate(ORDER_START_NODATE = as.POSIXct(ORDER_START_NODATE, format = "%H:%M:%S")) %>%
  mutate(MEDS_TAKEN_NODATE = as.POSIXct(MEDS_TAKEN_NODATE, format = "%H:%M:%S")) %>%
  mutate_at(vars(ORDER_TIME_BIN = ORDER_START_NODATE, MEDS_TAKEN_BIN = MEDS_TAKEN_NODATE),
            .funs = funs(cut(., breaks = "2 hours", labels = times))) %>% #order and dose times into 2 h bins
  #patients into age bins
  mutate(AGE_BIN = case_when(
    CURRENT_AGE >= 1 & CURRENT_AGE < 5 ~ "1-5y",
    CURRENT_AGE >= 5 & CURRENT_AGE < 7 ~ "5-7y",
    CURRENT_AGE >= 7 & CURRENT_AGE < 10 ~ "7-10y",
    CURRENT_AGE >= 10 & CURRENT_AGE < 15 ~ "10-15y",
    CURRENT_AGE >= 15 & CURRENT_AGE < 20 ~ "15-20y",
    TRUE ~ ">20y")) %>% 
  mutate(AGE_BIN = factor(AGE_BIN, levels = c("1-5y", "5-7y", "7-10y", "10-15y", "15-20y", ">20y"))) %>%
  #add dose bin info
  mutate(dose_bin = case_when(
    DOSE_MG_KG < 0.1 ~ 1,
    DOSE_MG_KG >= 0.1 & DOSE_MG_KG <= 0.2 ~ 2,
    DOSE_MG_KG > 0.2 & DOSE_MG_KG <= 0.3 ~ 3,
    DOSE_MG_KG > 0.3 & DOSE_MG_KG <= 0.4 ~ 4,
    DOSE_MG_KG > 0.4 & DOSE_MG_KG <= 0.5 ~ 5,
    DOSE_MG_KG > 0.5 & DOSE_MG_KG <= 0.6 ~ 6,
    DOSE_MG_KG > 0.6 & DOSE_MG_KG <= 0.7 ~ 7,
    DOSE_MG_KG > 0.7 & DOSE_MG_KG <= 0.8 ~ 8,
    DOSE_MG_KG > 0.8 & DOSE_MG_KG <= 0.9 ~ 9,
    DOSE_MG_KG > 0.9 & DOSE_MG_KG <= 1 ~ 10,
    DOSE_MG_KG > 1 ~ 11)) 

#Filter blood pressure (BP) data, require that doses be separated by at least 4 h in cases where patient rec'd > 1 dose as part of an order 
dat_ct2 <- dat_ct %>%
  group_by(PatID) %>%
  #add variable defining hrs from previous dose for each patient-order
  mutate(hfpd = difftime(MEDS_TAKEN, lag(MEDS_TAKEN, default=first(MEDS_TAKEN)), units = "hours")) %>%
  mutate(hfpd = round(hfpd, digits = 3)) %>%
  #remove any dose that was given less than 4 h from a prior dose
  group_by(PtDosID) %>%
  filter(!any(hfpd > 0 & hfpd < 4)) %>%
  ungroup(PatID, PtDosID) %>%
  #hours from dose for each bp measure
  mutate(BP_HRS_FROM_DOSE = round(difftime(BP_TIME, MEDS_TAKEN, units = c("hours")), digits = 2)) %>%
  separate(BP_VALUE, c("sbp", "dbp"), sep = "/") %>%
  mutate_at(vars(sbp, dbp), .funs = funs(as.numeric(.))) %>%
  select(PatID, PtDosID, DOSE_MG_KG, dose_bin, MEDS_TAKEN_BIN, BP_TIME, BP_HRS_FROM_DOSE, sbp, dbp)
  
#Additional filtering on BP data, require that doses be flanked closely by BP records, then compute treatment response measure
bpR <- dat_ct2 %>%
  group_by(PtDosID) %>%
  #filter doses with >= 1 bp record in the half-hour BEFORE dose
  filter(any(BP_HRS_FROM_DOSE < 0 & BP_HRS_FROM_DOSE >= -0.5)) %>%
  #filter doses with >= 1 bp record in the window from 0.5h to 5h AFTER dose
  filter(any(BP_HRS_FROM_DOSE >= 0 & BP_HRS_FROM_DOSE <= 2)) %>%
  #add variable that computes the MEAN SBP & DBP over the three hours AFTER DOSE
  mutate(sbp_mean = round(mean(sbp[which(BP_HRS_FROM_DOSE >= 0 & BP_HRS_FROM_DOSE <= 3)]), digits = 1)) %>%
  mutate(dbp_mean = round(mean(dbp[which(BP_HRS_FROM_DOSE >= 0 & BP_HRS_FROM_DOSE <= 3)]), digits = 1)) %>%
  #add variable to specify the dosing BP (BP measure recorded just prior to dose)
  mutate(dosing_sbp = sbp[max(which(BP_HRS_FROM_DOSE < 0))]) %>%
  mutate(dosing_dbp = dbp[max(which(BP_HRS_FROM_DOSE < 0))]) %>%
  #add var to specify max pct. change in systolic & diastolic bp -- from dosing bp to mean response over period set above  
  mutate(sbp_Pchg = sbp_mean/dosing_sbp - 1, dbp_Pchg = dbp_mean/dosing_dbp - 1) %>%
  mutate(sbp_Pchg = round(sbp_Pchg*100, digits = 2), dbp_Pchg = round(dbp_Pchg*100, digits = 2)) %>%
  distinct(PtDosID, .keep_all = TRUE) %>% ungroup(PtDosID) %>%
  #remove BP response outliers; lower & upper 3rd percentiles
  filter(sbp_Pchg < quantile(sbp_Pchg, probs = 0.97) & sbp_Pchg > quantile(sbp_Pchg, probs = 0.03)) %>%
  filter(dbp_Pchg < quantile(dbp_Pchg, probs = 0.97) & dbp_Pchg > quantile(dbp_Pchg, probs = 0.03)) %>%
  #remove lowest and highest doses
  filter(dose_bin == 2) %>% mutate(dose_bin = factor(dose_bin)) %>%
  select(PtDosID, MEDS_TAKEN_BIN, dose_bin, dosing_dbp, dbp_Pchg, dosing_sbp, sbp_Pchg)

#Model (anova) function, add covariate of dosing BP (systolic or diastolic) to the model 
model_aov<-function(df) {
  aov(sbp_Pchg ~ dosing_sbp + MEDS_TAKEN_BIN, data=df)
}

#Fit ANCOVA model & tidy
fit_aov <- bpR %>%
  select(dose_bin, MEDS_TAKEN_BIN, dosing_dbp, dbp_Pchg) %>% nest(-dose_bin) %>%
  mutate(model = map(data, model_aov)) %>%
  mutate(tidied = map(model, tidy)) %>% unnest(tidied)
write.csv(tukey_out, "./your directory.filename.csv")

#Post-hoc pairwise comparison tests; #Model (anova) function, add covariate of dosing BP (systolic or diastolic) to the model 
model_aov <- aov(dbp_Pchg ~ dosing_dbp + MEDS_TAKEN_BIN, data = bpR)
summary(model_aov)
#tukey test
model_tukey <- TukeyHSD(model_aov, "dose_time_bin")
tukey_out <- as.data.frame(model_tukey$MEDS_TAKEN_BIN)
write.csv(tukey_out, "./your directory.filename.csv")

#model stats for plot
plt_mod <- mstr_mod %>%
  filter(term == "dose_time_bin")

#avgs for line plot 
avgs <- mstr %>%
  group_by(dose_bin, dose_time_bin) %>%
  summarise(mean_sysC = round(mean(sbp_Pchg), digits = 3), med_sysC = median(sbp_Pchg),
            mean_diaC = round(mean(dbp_Pchg), digits = 3), med_diaC = median(dbp_Pchg)) %>%
  ungroup()

#dose counts for plot
count <- mstr %>%
  group_by(dose_bin) %>%
  summarise(cnt = n())


#boxplot single plots (diastolic or systolic) response vs time
p <- ggplot(mstr, aes(x = dose_time_bin, y = sbp_Pchg)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.5, outlier.alpha = 0.75) +
  geom_line(aes(x = dose_time_bin, y = med_sysC), data = avgs, group = 1, color = "red") +
  labs(x = "Dosing time", y = "Systolic response (% change)") +
  geom_text(aes(x = Inf, y = Inf, label = paste0("ANCOVA, P = ", signif(p.value, digits = 1)), hjust = 1.1, vjust = 1.2, fontface = "bold"),  data = plt_mod, size = 1.85) +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9, hjust = 1.3))
ggsave(p, filename = "./figs/HydralazineResponse.DoseBin2.Systolic.AncovaStats.pdf", width = 1.65, height = 2.25)


#boxplot response vs time for all dosage bins
#rename facet labels
fct_lbs <- c(
  `1` = "< 0.1 mg/kg",
  `3` = "0.2 >= mg/kg <= 0.3", 
  `4` = "0.3 >= mg/kg <= 0.4")

p <- ggplot(mstr, aes(x = dose_time_bin, y = sbp_Pchg)) +
  geom_boxplot(outlier.shape = NA, notch = TRUE) +
  geom_line(aes(x = dose_time_bin, y = med_sysC), data = avgs, group = 1, color = "red") +
  facet_wrap(~dose_bin, scales = "free_y", ncol = 3, labeller = as_labeller(fct_lbs)) +
  labs(x = "Hydralazine dose time", y = "BP response (% change in systolic)" ) +
  geom_text(aes(x = Inf, y = Inf, label = paste0("anova = ", signif(p.value, digits = 1)), hjust = 1.1, vjust = 1),  data = plt_mod, size = 2.25) +
  #geom_text(aes(x = Inf, y = Inf, label = paste0("n = ", cnt), hjust = 3, vjust = 1),  data = count, size = 2.25) +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 7), axis.title.y = element_text(size = 8, hjust = 1),
        axis.title.x = element_text(size = 8), strip.text.x = element_text(size = 8, margin = margin(0, 0, 0, 0)))
ggsave(p, filename = "./figs/HydralazineResponse.Systolic.DoseBins1.3.4.AncovaStats.pdf", width = 4.25, height = 2.15)





