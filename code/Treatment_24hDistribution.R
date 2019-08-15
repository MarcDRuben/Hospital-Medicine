
#Characterize the 24 distribution of hydralazine (orders AND first-doses)

#Required R packages
library(tidyverse)
library(cosinor2)
library(patchwork)

#Load patient data
input <- read_csv("./data/CCHMC_EMR_Hydralazine_Deidentified.csv")
  
dat <- input %>%
  #add a patient-order-dose id var
  mutate(PT_ORD_DOS = paste(PATIENT, ORDER_NUM, DOSE_NUM, sep = "_")) %>%
  #first-doses only
  filter(DOSE_NUM == 1) %>%
  distinct(PT_ORD_DOS, .keep_all = TRUE)
  
#transform to narrow format for plotting
tBin = c("12AM-2AM", "2AM-4AM", "4AM-6AM", "6AM-8AM", "8AM-10AM", "10AM-12PM", "12PM-2PM", "2PM-4PM", "4PM-6PM", "6PM-8PM", "8PM-10PM", "10PM-12AM") 
dat_nf <- dat %>%  
  gather(key = "ORDER_DOSE", value = "TIME", ORDER_START:DOSE_TAKEN, factor_key = TRUE) %>% 
  mutate(TIME = factor(TIME, levels = tBin), ORDER_DOSE = factor(ORDER_DOSE, levels = c("ORDER_START", "DOSE_TAKEN"))) %>%
  select(PT_ORD_DOS, ORDER_DOSE, TIME)

#Plot histogram of orders and first-doses by time-of-day
p <- ggplot(dat_nf, aes(x = TIME, group = ORDER_DOSE, fill = ORDER_DOSE)) +
  geom_bar(stat = "count", position = position_dodge(width = -0.35), colour = "black") +
  scale_fill_manual(values = c("#D55E00", "#999999"), labels = c("Orders", "First doses")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = 11, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 11),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), legend.position = "top",
        legend.direction = "horizontal", legend.key.size = unit(0.75, "line"), legend.title = element_blank(),
        legend.text = element_text(size = 11, margin = margin(r = 10, unit = "pt")),
        legend.margin = margin(0, 0, -7, 0, unit = "pt"))

#Cosinor regression to test for 24 h rhythm in orders and/or first-doses
cr_frst <- dat_nf %>%
  group_by(ORDER_DOSE, TIME) %>%
  summarise(COUNT = n()) %>% #tally orders and first-doses in each time bin
  mutate(TIME_INDX = seq(1, 12, by = 1)) #index time bins
cr_frst <- data.frame(cr_frst) #convert to data.frame

var <- c("ORDER_START", "DOSE_TAKEN")
fitD <- statsD <- NULL
for (i in 1:length(var)) {
  tmp <- cr_frst %>% filter(ORDER_DOSE == var[i])
  mod <- cosinor.lm(COUNT ~ time(TIME_INDX), data = tmp, period = 12) 
  fit <- data.frame(fit = mod$fit$fitted.values, TIME_INDX = tmp$TIME_INDX) %>% #df with fitted values for cosinor curve
    mutate(ORDER_DOSE = var[i])
  fitD <- bind_rows(fitD, fit)
  stats <- data.frame(pval = signif(summary(mod)$transformed.table$p.value[2], digits = 1),
                      #acr = abs(summary(mod)$transformed.table$estimate[3]),
                      phase = abs(correct.acrophase(mod)),
                      rsq = signif(cosinor.PR(mod)[2], digits = 2), ORDER_DOSE = var[i])
  statsD <- bind_rows(statsD, stats)
  }
statsD$ORDER_DOSE <- factor(statsD$ORDER_DOSE, levels = var)
fitD$ORDER_DOSE <- factor(fitD$ORDER_DOSE, levels = var)
  
#facet labels for plot
facet_labs <- c(
  `DOSE_TAKEN` = "First doses",
  `ORDER_START` = "Orders")

p2 <- ggplot(fitD) +
  geom_point(aes(x = TIME_INDX, y = COUNT, color = ORDER_DOSE), data = cr_frst, shape = 1, stroke = 0.75, size = 2.3) +
  scale_color_manual(values = c("#D55E00", "#999999")) +
  geom_smooth(aes(x = TIME_INDX, y = fit), data = fitD, colour = "black", size = 0.5) +
  scale_x_continuous(breaks = seq(1, 12, by = 1), labels = tBin) +
  geom_text(aes(x = 8, y = 350), data = statsD, label = paste0("P = ", statsD$pval), size = 2.75, fontface = "bold", color = "black") +
  facet_wrap(~ORDER_DOSE, ncol = 2, labeller = as_labeller(facet_labs)) +
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_blank(),
        legend.position = "none",
        strip.text.x = element_text(size = 9, margin = margin(0, 0, 0, 0)), 
        strip.background = element_rect(fill = "lightgrey", color = "lightgrey"))

combo_plot <- p + p2 + plot_layout(ncol = 1, heights = c(1, 0.75))
#ggsave(combo_plot, filename = "./figs/Hydralazine.OrdersVsDoses.pdf", width = 3, height = 5.5)
