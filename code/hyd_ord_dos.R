#clean & prep data on all hydralazine doses & flanking BP measures CCHMC from 2011-17
library(tidyverse)
library(readxl)
library(lubridate)
library(patchwork)
library(cosinor2)

#all sheet names from the excel worbook 
sheets <- excel_sheets("./data/Req9050_Smith_Hydralazine_YG_01252019.xlsx")

#patient info w/body weight into df (weight taken at order_start)
pt_wt_in <- read_excel("./data/Req9050_Smith_Hydralazine_YG_01252019.xlsx", sheet = "Weight")
pt_wt <- pt_wt_in %>%
  select(PAT_MRN_ID, ORDER_START_TIME, weight_kg = WEIGHT_KG)

#patient info w/all other demo info into df
pt_demo_in <- read_excel("./data/Req9050_Smith_Hydralazine_YG_01252019.xlsx", sheet = "Hydralazine")
pt_demo <- pt_demo_in %>%
  mutate(PtDosID = paste(PAT_MRN_ID, taken_time, sep = "_")) %>%
  mutate(PtOrdID = paste(PAT_MRN_ID, ORDER_START_TIME, sep = "_"))
  
#join all patient info, remove dups, convert doses to mg/kg, and clean up for plotting
labs = c("12AM-2AM", "2AM-4AM", "4AM-6AM", "6AM-8AM", "8AM-10AM", "10AM-12PM","12PM-2PM", "2PM-4PM", "4PM-6PM", "6PM-8PM", "8PM-10PM", "10PM-12AM") 
labs_hires = c("12AM-1AM", "1AM-2AM", "2AM-3AM", "3AM-4AM", "4AM-5AM", "5AM-6AM","6AM-7AM", "7AM-8AM", "8AM-9AM", "9AM-10AM", "10AM-11AM", "11AM-12PM",
               "12PM-1PM", "1PM-2PM", "2PM-3PM", "3PM-4PM", "4PM-5PM", "5PM-6PM","6PM-7PM", "7PM-8PM", "8PM-9PM", "9PM-10PM", "10PM-11PM", "11PM-12AM") 

pt_info <- left_join(pt_demo, pt_wt, by = c("PAT_MRN_ID", "ORDER_START_TIME")) %>%
  filter(ROUTE == "Intravenous") %>%
  distinct(PtDosID, .keep_all = TRUE) %>% 
  mutate(dose_mg_kg = ifelse(UNIT == "mg", round(HV_DISCRETE_DOSE/weight_kg, digits = 2), HV_DISCRETE_DOSE)) %>%
  filter(!is.na(dose_mg_kg)) %>%
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
  select(PtDosID, PtOrdID, PtID = PAT_MRN_ID, sex = SEX, age = CurrentAge, race, order_time = ORDER_START_TIME,
         dose_time = taken_time, dose_mg_kg, dose_bin, freq = FREQUENCY, dept = DEPARTMENT_NAME) %>%
  mutate_at(vars(order_time_nodate = order_time, dose_time_nodate = dose_time),
            .funs = funs(format(ymd_hms(.), "%H:%M:%S"))) %>%
  mutate(order_time_nodate = as.POSIXct(order_time_nodate, format = "%H:%M:%S")) %>%
  mutate(dose_time_nodate = as.POSIXct(dose_time_nodate, format = "%H:%M:%S")) %>%
  mutate_at(vars(order_time_bin = order_time_nodate, dose_time_bin = dose_time_nodate),
              .funs = funs(cut(., breaks = "2 hours", labels = labs))) %>%
  mutate_at(vars(order_time_bin_hires = order_time_nodate, dose_time_bin_hires = dose_time_nodate),
            .funs = funs(cut(., breaks = "1 hours", labels = labs_hires))) %>%
  #add age bins
  mutate(age_bin = case_when(
  age >= 1 & age < 5 ~ "age 1-5",
  age >= 5 & age < 7 ~ "age 5-7",
  age >= 7 & age < 10 ~ "age 7-10",
  age >= 10 & age < 15 ~ "age 10-15",
  age >= 15 & age < 20 ~ "age 15-20",
  TRUE ~ "age > 20")) %>%
  mutate(age_bin = factor(age_bin, levels = c("age 1-5","age 5-7","age 7-10","age 10-15", "age 15-20", "age > 20")))
write.table(pt_info, "./data/Hydralazine_PatientDemos.csv", sep = "," , row.names = FALSE)

#first (or only) doses of an order
first_dose <- pt_info %>%
  group_by(PtOrdID) %>%
  mutate(doses_per_order = ifelse(n() == 1, "one", "multiple")) %>%
  filter(dose_time == min(dose_time)) %>%
  gather(key = "order_dose", value = "time", order_time_bin:dose_time_bin, factor_key = TRUE) %>%
  ungroup(PtOrdID) %>%
  mutate(time = factor(time, levels = labs), order_dose = factor(order_dose,
                                                                 levels = c("order_time_bin", "dose_time_bin")))

#plot first (or only) doses of an order
p <- ggplot(first_dose, aes(x = time, group = order_dose, fill = order_dose)) +
  geom_bar(stat = "count", position = position_dodge(width = -0.35), colour = "black") +
  scale_fill_manual(values = c("#D55E00", "#999999"), labels = c("Orders", "First doses")) +
  scale_y_continuous(limits = c(0, 625), expand = c(0, 0)) +
  #facet_wrap(~age, scales = "free") +
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 9),
        #axis.title.y = element_text(size = 13),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), legend.position = "top",
        legend.direction = "horizontal", legend.key.size = unit(0.75, "line"), legend.title = element_blank(),
        legend.text = element_text(size = 10, margin = margin(r = 10, unit = "pt")),
        legend.margin = margin(0, 0, -7, 0, unit = "pt"))

#cosinor regression; tally number of first (or only) doses of an order at each time
cr_first_dose <- first_dose %>%
  group_by(order_dose, time) %>%
  summarise(number = n()) %>%
  mutate(time_adj = seq(1, 12, by = 1))
cr_first_dose <- data.frame(cr_first_dose)

var <- c("order_time_bin", "dose_time_bin")
fitD <- statsD <- NULL
for (i in 1:length(var)) {
  tmp <- cr_first_dose %>% filter(order_dose == var[i])
  mod <- cosinor.lm(number ~ time(time_adj), data = tmp, period = 12) 
  fit <- data.frame(fit = mod$fit$fitted.values, time_adj = tmp$time_adj) %>% #df with fitted values for cosinor curve
    mutate(order_dose = var[i])
  fitD <- bind_rows(fitD, fit)
  stats <- data.frame(pval = signif(summary(mod)$transformed.table$p.value[2], digits = 1),
                      pval2 = signif(cosinor.detect(mod)[4], digits = 2), 
                      acr = abs(summary(mod)$transformed.table$estimate[3]),
                      acr2 = abs(correct.acrophase(mod)),
                      rsq = signif(cosinor.PR(mod)[2], digits = 2), order_dose = var[i])
  statsD <- bind_rows(statsD, stats)
  }
statsD$order_dose <- factor(statsD$order_dose, levels = var)
fitD$order_dose <- factor(fitD$order_dose, levels = var)
  
#facet labels
facet_labs <- c(
  `dose_time_bin` = "First doses",
  `order_time_bin` = "Orders")

p2 <- ggplot(fitD) +
  geom_point(aes(x = time_adj, y = number, color = order_dose), data = cr_first_dose, shape = 1, stroke = 0.75, size = 2.3) +
  scale_color_manual(values = c("#D55E00", "#999999")) +
  geom_smooth(aes(x = time_adj, y = fit), data = fitD, colour = "black", size = 0.5) +
  scale_x_continuous(breaks = seq(1, 12, by = 1), labels = labs) +
  geom_text(aes(x = 8, y = 350), data = statsD, label = paste0("P = ", statsD$pval), size = 2.75, fontface = "bold", color = "black") +
  facet_wrap(~order_dose, ncol = 2, labeller = as_labeller(facet_labs)) +
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 9, margin = margin(0, 0, 0, 0)), 
        strip.background = element_rect(fill = "lightgrey", color = "lightgrey"))

patch <- p + p2 + plot_layout(ncol = 1, heights = c(1, 0.75))
ggsave(patch, filename = "./figs/Hydralazine.OrdersVsDoses.pdf", width = 3, height = 5.5)

##

#plot single doses only
p3 <- ggplot(multi_dose, aes(x = time, group = order_dose, fill = order_dose)) +
  geom_bar(stat = "count", position = position_dodge(width = -0.35), colour = "black") +
  scale_fill_manual(values = c("#D55E00", "#999999"), labels = c("taken (first dose)", "ordered")) +
  scale_y_continuous(limits = c(0, 350), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 11),
        #axis.title.y = element_text(size = 13),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), legend.position = "top",
        legend.direction = "horizontal", legend.key.size = unit(1, "line"), legend.title = element_blank(),
        legend.text = element_text(size = 10, margin = margin(r = 10, unit = "pt")))



