
#Characterize time-of-day pattern of response to IV-hydralazine treatment

#Required R packages
library(tidyverse)
library(broom)

input <- read_csv("./data/CCHMC_EMR_Hydralazine_Deidentified.csv")

tBin = c("12AM-2AM", "2AM-4AM", "4AM-6AM", "6AM-8AM", "8AM-10AM", "10AM-12PM", "12PM-2PM", "2PM-4PM", "4PM-6PM", "6PM-8PM", "8PM-10PM", "10PM-12AM") 

dat <- input %>%
  #add a patient-order-dose identifier variable
  mutate(PT_ORD_DOS = paste(PATIENT, ORDER_NUM, DOSE_NUM, sep = "_")) %>%
  #one row per dose
  #group_by(PT_ORD_DOS) %>% distinct(PT_ORD_DOS, .keep_all = TRUE) %>%
  #IV hydralazine only, remove doses w/o concentration reported
  filter(ROUTE == "Intravenous" & !is.na(DOSE_MG_KG)) %>%
  #add dose bins
  mutate(DOSAGE_BIN = case_when(
    DOSE_MG_KG < 0.1 ~ 1,
    DOSE_MG_KG >= 0.1 & DOSE_MG_KG <= 0.2 ~ 2,
    DOSE_MG_KG > 0.2 & DOSE_MG_KG <= 0.3 ~ 3,
    TRUE ~ 4
  )) %>%
  mutate(DOSE_TAKEN = factor(DOSE_TAKEN, levels = tBin)) %>%
  ungroup()

#note: the data above was prefiltered to remove any doses that were administered < 4 h from a prior dose in that patient (using the "HOURS_FROM_PRIOR_DOSE" variable)
#additional filtering: remove any doses that are not flanked closely by BP records
dat_f <- dat %>%
  group_by(PT_ORD_DOS) %>%
  #filter doses with >= 1 bp record in the half-hour BEFORE dose
  filter(any(BP_TIME_HRS_FROM_DOSE < 0 & BP_TIME_HRS_FROM_DOSE >= -0.5)) %>%
  #filter doses with >= 1 bp record in the window from 0.5h to 5h AFTER dose
  filter(any(BP_TIME_HRS_FROM_DOSE >= 0 & BP_TIME_HRS_FROM_DOSE <= 2)) %>%
  ungroup()
  
#compute measure of treatment response
dat_rsp <- dat_f %>%
  separate(BP_VALUE, c("sbp", "dbp"), sep = "/") %>%
  mutate_at(vars(sbp, dbp), .funs = funs(as.numeric(.))) %>%
  #compute mean SBP & DBP across the 3 h AFTER a dose
  group_by(PT_ORD_DOS) %>%
  mutate(sbp_mean = round(mean(sbp[which(BP_TIME_HRS_FROM_DOSE >= 0 & BP_TIME_HRS_FROM_DOSE <= 3)]), digits = 1)) %>%
  mutate(dbp_mean = round(mean(dbp[which(BP_TIME_HRS_FROM_DOSE >= 0 & BP_TIME_HRS_FROM_DOSE <= 3)]), digits = 1)) %>%
  #add variable to specify the dosing BP (BP measure recorded just prior to dose)
  mutate(DOSING_SBP = sbp[max(which(BP_TIME_HRS_FROM_DOSE < 0))]) %>%
  mutate(DOSING_DBP = dbp[max(which(BP_TIME_HRS_FROM_DOSE < 0))]) %>%
  #compute percent change in SBP -- from dosing BP to mean response over post-dose window set above (3 h) 
  mutate(SBP_CHG = sbp_mean/DOSING_SBP - 1, DBP_CHG = dbp_mean/DOSING_DBP - 1) %>%
  mutate(SBP_CHG = round(SBP_CHG*100, digits = 2), DBP_CHG = round(DBP_CHG*100, digits = 2)) %>%
  distinct(PT_ORD_DOS, .keep_all = TRUE) %>% ungroup()

#post-analysis filtering: remove response outliers; lower & upper 3rd percentiles
dat_pf <- dat_rsp %>%
  filter(SBP_CHG < quantile(SBP_CHG, probs = 0.97) & SBP_CHG > quantile(SBP_CHG, probs = 0.03) &
  DBP_CHG < quantile(DBP_CHG, probs = 0.97) & DBP_CHG > quantile(DBP_CHG, probs = 0.03)) %>%
  #remove lowest and highest doses
  mutate(DOSAGE_BIN = factor(DOSAGE_BIN)) %>%
  select(PT_ORD_DOS, DOSE_TAKEN, DOSAGE_BIN, DOSING_DBP, DBP_CHG, DOSING_SBP, SBP_CHG)

#tidy prep to test linear model ANOVA
nst <- dat_pf %>%
  select(DOSAGE_BIN, DOSE_TAKEN, DOSING_SBP, DOSING_DBP, SBP_CHG, DBP_CHG) %>%
  nest(-c(DOSAGE_BIN))

#model function, add covariate of dosing BP to the model 
aov_Sbp <-function(df) {
  aov(SBP_CHG ~ DOSING_SBP + DOSE_TAKEN, data=df)
}
aov_Dbp <-function(df) {
  aov(DBP_CHG ~ DOSING_DBP + DOSE_TAKEN, data=df)
}

#fit  model & tidy
out_sbp <- nst %>% #sbp
  mutate(model = map(data, aov_Sbp)) %>%
  mutate(tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(RESPONSE = "SBP")

out_dbp <- nst %>% #dbp
  mutate(model = map(data, aov_Dbp)) %>%
  mutate(tidied = map(model, tidy)) %>%
  unnest(tidied) %>%
  mutate(RESPONSE = "DBP")

ancova_comb <- bind_rows(out_sbp, out_dbp) %>%
  filter(term == "DOSE_TAKEN")

#transform data to narrow format for boxplot
resp_nf <- dat_pf %>% select(DOSAGE_BIN, DOSE_TAKEN, SBP=SBP_CHG, DBP=DBP_CHG) %>%
  gather(RESPONSE, PCT_CHG, SBP:DBP)

#median responses for smooth curve in boxplot
resp_med <- resp_nf %>% group_by(DOSAGE_BIN, RESPONSE, DOSE_TAKEN) %>%
  summarise(PCT_CHG = round(median(PCT_CHG), digits = 3))

#facet labels
  fct_lbs1 <- c(
    `DBP` = "DBP",
    `SBP` = "SBP")

#boxplot response; dose bin 2 SBP and DBP
p <- ggplot(resp_nf, aes(x = DOSE_TAKEN, y = PCT_CHG)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.5) +
  geom_line(aes(x = DOSE_TAKEN, y = PCT_CHG), data = resp_med, group = 1, color = "red") +
  scale_x_discrete(labels = c("12", "", "4", "", "8", "", "12", "", "4", "", "8", "")) +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(cols = vars(DOSAGE_BIN), rows = vars(RESPONSE), scales = "free_y") +
  labs(x = "Dosing time", y = "Response (% change)") +
  geom_label(aes(x = 5.75, y = 3, label = paste0("ANCOVA, p = ", signif(p.value, digits = 1))),
             data = ancova_comb, size = 2) +
  theme(axis.text.x = element_text(size = 8, hjust = 0.5, vjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 9),
        panel.background = element_blank(),
        strip.background.y = element_blank(),
        strip.text.x = element_text(size = 9, margin = margin(0, 0, 0, 0)))
#ggsave(p2, filename = "./figs/HydralazineResponse.Boxplot.AncovaStats.pdf", width = 5, height = 4)
