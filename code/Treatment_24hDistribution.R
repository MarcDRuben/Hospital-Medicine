
#Characterize the 24 distribution of hydralazine treatment (drug order AND administration) for all inpatients who received the drug at CCHMC from 2011-17

#Required R packages
library(tidyverse)
library(readxl)
library(lubridate)
library(cosinor2)
library(patchwork)

#Load patient data. Create unique identifiers for each patient-order and patient-dose (note: an order may have multiple doses associated)
dat <- read_csv("./data/CCHMC_EMR_Hydralazine_Deidentified.csv") %>%
  mutate(PtOrdID = paste(PatID, ORDER_START_TIME, sep = "_")) %>%
  mutate(PtDosID = paste(PatID, MEDS_TAKEN, sep = "_"))
  
#Vector specifying time bins for analysis; 2 h intervals used
times = c("12AM-2AM", "2AM-4AM", "4AM-6AM", "6AM-8AM", "8AM-10AM", "10AM-12PM", "12PM-2PM", "2PM-4PM", "4PM-6PM", "6PM-8PM", "8PM-10PM", "10PM-12AM") 

#Clean & transform data for analysis / visualization 
dat_ct <- dat %>%
  filter(ROUTE == "Intravenous") %>%  #optional filter on ROUTE
  distinct(PtDosID, .keep_all = TRUE) %>% 
  mutate(DOSE_MG_KG = ifelse(UNIT == "mg", round(HV_DISCRETE_DOSE/WEIGHT_KG, digits = 2), HV_DISCRETE_DOSE)) %>%  #convert all doses to mg/kg
  filter(!is.na(DOSE_MG_KG)) %>% #optional remove doses without information on concentration 
  mutate_at(vars(ORDER_START_TIME_NODATE = ORDER_START_TIME, MEDS_TAKEN_NODATE = MEDS_TAKEN),
            .funs = funs(format(ymd_hms(.), "%H:%M:%S"))) %>% #add variable for order and taken times without date
  mutate(ORDER_START_TIME_NODATE = as.POSIXct(ORDER_START_TIME_NODATE, format = "%H:%M:%S")) %>%
  mutate(MEDS_TAKEN_NODATE = as.POSIXct(MEDS_TAKEN_NODATE, format = "%H:%M:%S")) %>%
  mutate_at(vars(ORDER_TIME_BIN = ORDER_START_TIME_NODATE, DOSE_TIME_BIN = MEDS_TAKEN_NODATE),
              .funs = funs(cut(., breaks = "2 hours", labels = times))) %>% #order and dose times into 2 h bins
  mutate(AGE_BIN = case_when(
  CURRENT_AGE >= 1 & CURRENT_AGE < 5 ~ "1-5y",
  CURRENT_AGE >= 5 & CURRENT_AGE < 7 ~ "5-7y",
  CURRENT_AGE >= 7 & CURRENT_AGE < 10 ~ "7-10y",
  CURRENT_AGE >= 10 & CURRENT_AGE < 15 ~ "10-15y",
  CURRENT_AGE >= 15 & CURRENT_AGE < 20 ~ "15-20y",
  TRUE ~ ">20y")) %>%   #patients into age bins
  mutate(AGE_BIN = factor(AGE_BIN, levels = c("1-5y","5-7y","7-10y","10-15y", "15-20y", ">20y")))

#Subset the first-doses administered from each order (note: in some cases the first dose is the only dose in an order)
dat_first <- dat_ct %>%
  group_by(PtOrdID) %>%
  filter(MEDS_TAKEN == min(MEDS_TAKEN)) %>%
  gather(key = "STAGE", value = "TIME", ORDER_TIME_BIN:DOSE_TIME_BIN, factor_key = TRUE) %>% #convert data to narrow format for plotting
  ungroup(PtOrdID) %>%
  mutate(TIME = factor(TIME, levels = times), STAGE = factor(STAGE, levels = c("ORDER_TIME_BIN", "DOSE_TIME_BIN"))) %>%
  select(PatID, AGE_BIN, PtOrdID, PtDosID, STAGE, TIME) #optional select whichever additional variables of interest for analysis

#Plot histogram of orders and first-doses by time-of-day
p <- ggplot(dat_first, aes(x = TIME, group = STAGE, fill = STAGE)) +
  geom_bar(stat = "count", position = position_dodge(width = -0.35), colour = "black") +
  scale_fill_manual(values = c("#D55E00", "#999999"), labels = c("Orders", "First doses")) +
  scale_y_continuous(limits = c(0, 625), expand = c(0, 0)) +
  #facet_wrap(~age, scales = "free") +
  theme(axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), legend.position = "top",
        legend.direction = "horizontal", legend.key.size = unit(0.75, "line"), legend.title = element_blank(),
        legend.text = element_text(size = 11, margin = margin(r = 10, unit = "pt")),
        legend.margin = margin(0, 0, -7, 0, unit = "pt"))

#Cosinor regression to test for 24 h rhythm in orders and/or first-doses
cr_dat_first <- dat_first %>%
  group_by(STAGE, TIME) %>%
  summarise(COUNT = n()) %>% #tally the number of orders and first-doses in each time bin
  mutate(TIME_INDX = seq(1, 12, by = 1)) #add index for time bins
cr_dat_first <- data.frame(cr_dat_first) #convert to data.frame for cosinor regression function

var <- c("ORDER_TIME_BIN", "DOSE_TIME_BIN")
fitD <- statsD <- NULL
for (i in 1:length(var)) {
  tmp <- cr_dat_first %>% filter(STAGE == var[i])
  mod <- cosinor.lm(COUNT ~ time(TIME_INDX), data = tmp, period = 12) 
  fit <- data.frame(fit = mod$fit$fitted.values, TIME_INDX = tmp$TIME_INDX) %>% #df with fitted values for cosinor curve
    mutate(STAGE = var[i])
  fitD <- bind_rows(fitD, fit)
  stats <- data.frame(pval = signif(summary(mod)$transformed.table$p.value[2], digits = 1),
                      pval2 = signif(cosinor.detect(mod)[4], digits = 2), 
                      acr = abs(summary(mod)$transformed.table$estimate[3]),
                      acr2 = abs(correct.acrophase(mod)),
                      rsq = signif(cosinor.PR(mod)[2], digits = 2), STAGE = var[i])
  statsD <- bind_rows(statsD, stats)
  }
statsD$STAGE <- factor(statsD$STAGE, levels = var)
fitD$STAGE <- factor(fitD$STAGE, levels = var)
  
#facet labels for plot
facet_labs <- c(
  `DOSE_TIME_BIN` = "First doses",
  `ORDER_TIME_BIN` = "Orders")

p2 <- ggplot(fitD) +
  geom_point(aes(x = TIME_INDX, y = COUNT, color = STAGE), data = cr_dat_first, shape = 1, stroke = 0.75, size = 2.3) +
  scale_color_manual(values = c("#D55E00", "#999999")) +
  geom_smooth(aes(x = TIME_INDX, y = fit), data = fitD, colour = "black", size = 0.5) +
  scale_x_continuous(breaks = seq(1, 12, by = 1), labels = times) +
  geom_text(aes(x = 8, y = 350), data = statsD, label = paste0("P = ", statsD$pval), size = 2.75, fontface = "bold", color = "black") +
  facet_wrap(~STAGE, ncol = 2, labeller = as_labeller(facet_labs)) +
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 9, margin = margin(0, 0, 0, 0)), 
        strip.background = element_rect(fill = "lightgrey", color = "lightgrey"))

combo_plot <- p + p2 + plot_layout(ncol = 1, heights = c(1, 0.75))
#ggsave(patch, filename = "./figs/Hydralazine.OrdersVsDoses.pdf", width = 3, height = 5.5)
