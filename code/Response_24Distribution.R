#clean & prep BP data flanking hydralazine doses @ CCHMC from 2011-17
library(tidyverse)
library(readxl)
library(lubridate)
library(broom)

#all sheet names from the excel worbook 
sheets <- excel_sheets("./data/Req9050_Smith_Hydralazine_YG_01252019.xlsx")
#BP data into single df
bp_sheets <- sheets[grepl("BP", sheets)]
bp_df_in <- map_dfr(bp_sheets, ~read_excel(path = "./data/Req9050_Smith_Hydralazine_YG_01252019.xlsx", sheet = .x))



#add dose bin info
mutate(dose_bin = case_when(
  dose_mg_kg < 0.1 ~ 1,
  dose_mg_kg >= 0.1 & dose_mg_kg <= 0.2 ~ 2,
  dose_mg_kg > 0.2 & dose_mg_kg <= 0.3 ~ 3,
  dose_mg_kg > 0.3 & dose_mg_kg <= 0.4 ~ 4,
  dose_mg_kg > 0.4 & dose_mg_kg <= 0.5 ~ 5,
  dose_mg_kg > 0.5 & dose_mg_kg <= 0.6 ~ 6,
  dose_mg_kg > 0.6 & dose_mg_kg <= 0.7 ~ 7,
  dose_mg_kg > 0.7 & dose_mg_kg <= 0.8 ~ 8,
  dose_mg_kg > 0.8 & dose_mg_kg <= 0.9 ~ 9,
  dose_mg_kg > 0.9 & dose_mg_kg <= 1 ~ 10,
  dose_mg_kg > 1 ~ 11)) %>%

#add patient-dose IDs & clean
bp_df <- bp_df_in %>%
  mutate(PtDosID = paste(PAT_MRN_ID, Meds_Taken_Time, sep = "_")) %>%
  #filter doses that are separated by at least 6 h in cases where patient rec'd > 1 dose 
  group_by(PAT_MRN_ID) %>%
  #add variable for hrs from previous dose for each patient
  mutate(hfpd = difftime(Meds_Taken_Time, lag(Meds_Taken_Time, default=first(Meds_Taken_Time)), units = "hours")) %>%
  mutate(hfpd = round(hfpd, digits = 3)) %>%
  #remove any dose to close to prior dose (hfpd < 6 h)
  group_by(PtDosID) %>%
  filter(!any(hfpd > 0 & hfpd < 4)) %>%
  ungroup(PAT_MRN_ID, PtDosID) %>%
  #hours from dose for each bp measure
  mutate(bp_hrs_from_dose = round(difftime(BP_Recorded_Time, Meds_Taken_Time, units = c("hours")), digits = 2)) %>%
  separate(meas_value, c("sbp", "dbp"), sep = "/") %>%
  mutate_at(vars(sbp, dbp), .funs = funs(as.numeric(.))) %>%
  select(PAT_MRN_ID, PtDosID, dose_time = Meds_Taken_Time, bp_time = BP_Recorded_Time, bp_hrs_from_dose, sbp, dbp)
  
#filter doses with flanking bp records, compute response measure
bpr_df <- bp_df %>%
  group_by(PtDosID) %>%
  #filter doses with >= 1 bp record in the half-hour BEFORE dose
  filter(any(bp_hrs_from_dose < 0 & bp_hrs_from_dose >= -0.5)) %>%
  #filter doses with >= 1 bp record in the window from 0.5h to 5h AFTER dose
  filter(any(bp_hrs_from_dose >= 0 & bp_hrs_from_dose <= 2)) %>%
  #mean sbp & dbp over the five hours AFTER dose
  mutate(sbp_mean = round(mean(sbp[which(bp_hrs_from_dose >= 0 & bp_hrs_from_dose <= 3)]), digits = 1)) %>%
  mutate(dbp_mean = round(mean(dbp[which(bp_hrs_from_dose >= 0 & bp_hrs_from_dose <= 3)]), digits = 1)) %>%
  #dosing BP
  mutate(dosing_sbp = sbp[max(which(bp_hrs_from_dose < 0))]) %>%
  mutate(dosing_dbp = dbp[max(which(bp_hrs_from_dose < 0))]) %>%
  #var for max pct. change in systolic & diastolic bp -- from dosing bp to mean response over period set above  
  mutate(sbp_Pchg = sbp_mean/dosing_sbp - 1, dbp_Pchg = dbp_mean/dosing_dbp - 1) %>%
  #mutate(sbp_Pchg = -(dosing_sbp - sbp_mean) / dosing_sbp, dbp_Pchg = -(dosing_dbp - dbp_mean) / dosing_dbp) %>%
  mutate(sbp_Pchg = round(sbp_Pchg*100, digits = 2), dbp_Pchg = round(dbp_Pchg*100, digits = 2)) %>%
  #set BP changes above zero to zero
  #mutate(sbp_Pchg = ifelse(sbp_Pchg > 0, 0, sbp_Pchg), dbp_Pchg = ifelse(dbp_Pchg > 0, 0, dbp_Pchg)) %>%
  #add variable for no clinical response (NCR) ideal clinical resonse (ICR), excess response (ER)
  mutate(sbp_ClinResp = case_when(
    sbp_Pchg > -10 ~ "ncr",
    sbp_Pchg <= -10 & sbp_Pchg > -25 ~ "icr", 
    TRUE ~ "er")) %>%
  mutate(dbp_ClinResp = case_when(
    dbp_Pchg > -10 ~ "ncr",
    dbp_Pchg <= -10 & dbp_Pchg > -25 ~ "icr", 
    TRUE ~ "er")) %>%
  ungroup(PtDosID) %>%
  mutate(dosing_sbp_bin = as.numeric(cut_number(dosing_sbp, 8))) %>%
  mutate(dosing_dbp_bin = as.numeric(cut_number(dosing_dbp, 8))) %>%
  distinct(PtDosID, .keep_all = TRUE) %>%
  #remove BP response outliers; lower & upper 3rd percentiles
  #filter(dosing_sbp < quantile(dosing_sbp, probs = 0.99) & dosing_sbp > quantile(dosing_sbp, probs = 0.01)) %>%
  #filter(dosing_dbp < quantile(dosing_dbp, probs = 0.99) & dosing_dbp > quantile(dosing_dbp, probs = 0.01)) %>%
  filter(sbp_Pchg < quantile(sbp_Pchg, probs = 0.97) & sbp_Pchg > quantile(sbp_Pchg, probs = 0.03)) %>%
  filter(dbp_Pchg < quantile(dbp_Pchg, probs = 0.97) & dbp_Pchg > quantile(dbp_Pchg, probs = 0.03))
 
#patient demos master df, grab columns for joining with BP data
labs = c("12AM-2AM", "2AM-4AM", "4AM-6AM", "6AM-8AM", "8AM-10AM", "10AM-12PM","12PM-2PM", "2PM-4PM", "4PM-6PM", "6PM-8PM", "8PM-10PM", "10PM-12AM") 
labs_hires = c("12AM-1AM", "1AM-2AM", "2AM-3AM", "3AM-4AM", "4AM-5AM", "5AM-6AM","6AM-7AM", "7AM-8AM", "8AM-9AM", "9AM-10AM", "10AM-11AM", "11AM-12PM",
               "12PM-1PM", "1PM-2PM", "2PM-3PM", "3PM-4PM", "4PM-5PM", "5PM-6PM","6PM-7PM", "7PM-8PM", "8PM-9PM", "9PM-10PM", "10PM-11PM", "11PM-12AM") 
ages <- c("age 1-5", "age 5-7", "age 7-10", "age 10-15", "age 15-20", "age > 20")
pt_info <- read_csv("./data/Hydralazine_PatientDemos.csv") %>%
  mutate(age_bin = factor(age_bin, levels = ages)) %>%
  mutate(order_time_bin_hires = factor(order_time_bin_hires, levels = labs_hires)) %>%
  mutate(dose_time_bin = factor(dose_time_bin, levels = labs)) %>%
  mutate(dose_time_bin_hires = factor(dose_time_bin_hires, levels = labs_hires)) %>%
  #remove higher doses
  #filter(dose_bin == 2) %>%
  #filter(dose_bin == 1 | dose_bin == 3 | dose_bin == 4) %>%
  mutate(dose_bin = factor(dose_bin)) %>%
  select(PtDosID, PtOrdID, sex, age, age_bin, race, order_time, dose_time2 = dose_time, dose_time_nodate, dose_time_bin, dose_time_bin_hires,
         dose_mg_kg, dose_bin, freq, dept)

#join bp data with patient demo for master df
mstr <- inner_join(bpr_df, pt_info, by = "PtDosID") %>%
  select(PtDosID, sex, age_bin, race, dose_mg_kg, dose_bin, freq, dept, dose_time, dose_time_nodate,
         dose_time_bin, dose_time_bin_hires, dosing_sbp, dosing_sbp_bin, sbp_Pchg, sbp_ClinResp, dosing_dbp, dosing_dbp_bin, dbp_Pchg, dbp_ClinResp)
#write.table(mstr, "./data/Hydralazine_Filtered_BPresponses.csv", sep = "," , row.names = FALSE)  

#vars for linear modeling (adjust for diastolic vs. systolic)
mstr_c <- mstr %>%
  select(dose_bin, dose_time_bin, dosing_sbp, sbp_Pchg) %>%
  nest(-dose_bin)

#model (anova) function, add covariate of dosing BP (systolic or diastolic) to the model 
model_aov<-function(df) {
  aov(sbp_Pchg ~ dosing_sbp + dose_time_bin, data=df)
}

#fit model & tidy
mstr_mod <- mstr_c %>%
  mutate(model = map(data, model_aov)) %>%
  mutate(tidied = map(model, tidy)) %>%
  unnest(tidied)

#post-hoc pairwise comparison test, dose_bin #2 only --where multivariate analysis showed differences
posthoc <- mstr %>% filter(dose_bin == 2)
model_aov <- aov(dbp_Pchg ~ dosing_dbp + dose_time_bin, data = posthoc)
summary(model_aov)
#tukey test
model_tukey <- TukeyHSD(model_aov, "dose_time_bin")
tukey_out <- as.data.frame(model_tukey$dose_time_bin)
write.csv(tukey_out, "./data/HydralazineResponse.DoseBin2.Diastolic.TukeyStats.csv")

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


#violin plot dosa44ge distribution ()
vio <- ggplot(mstr, aes(x = "", y = dose_mg_kg)) +
  geom_violin() +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.1)) +
  scale_x_discrete(expand = c(50, 50)) +
  geom_rect(fill = "#009E73", alpha = 0.005, xmin = 0, xmax = Inf,  ymin = 0.1, ymax = 0.2) +
  #geom_rect(fill = "#E69F00", alpha = 0.005, xmin = 0, xmax = Inf,  ymin = 0, ymax = 0.1) +
  #geom_rect(fill = "#56B4E9", alpha = 0.005, xmin = 0, xmax = Inf,  ymin = 0.2, ymax = 0.3) +
  #geom_rect(fill = "#D55E00", alpha = 0.005, xmin = 0, xmax = Inf,  ymin = 0.3, ymax = 0.4) +
  labs(x = "Hydralazine doses", y = "Dosage (mg/kg)") +
  theme_bw() +
  theme(panel.border = element_blank(), axis.title = element_text(size = 9), axis.text.y = element_text(size = 9),
        axis.text.x = element_blank(), panel.grid = element_blank(), axis.ticks = element_blank())
ggsave(vio, filename = "./figs/ViolinPlot.Hydralazine.DosageBin2.Distribution.pdf", width = 1.75, height = 2.25)



#histogram clinical response by time bin, facet by dosing bin
clinR <- mstr %>%
  group_by(dosing_sbp_bin, dose_time_bin, sbp_ClinResp) %>%
  summarise(cnt = n()) %>%
  mutate(percent = round(cnt / sum(cnt), digits = 2)) %>%
  ungroup()
p <- ggplot(clinR, aes(x = dose_time_bin, y = percent)) +
  geom_bar(aes(fill = sbp_ClinResp), stat = "identity") +
  facet_wrap(~dosing_sbp_bin)

#scratch - continuous time data
p2 <- ggplot(mstr, aes(x = dose_time_nodate, y = sbp_Pchg)) +
  geom_point(alpha = 0.4) +
  #scale_y_continuous(limits = quantile(mstr_dos$sbp_chg, c(0.1, 0.9))) +
  facet_wrap(~ dosing_sbp_bin, scales = "free_y", ncol = 5) +
  theme(axis.text.x = element_text(angle = 90))





