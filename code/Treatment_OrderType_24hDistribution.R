
#Characterize the 24 distribution of hydralazine doses by (1) order type ("once", "as needed", or "scheduled") and
#and (2) dose type ("first-dose" vs "all others")

#Required R packages
library(tidyverse)
library(cosinor2)

#Load patient data
input <- read_csv("./data/CCHMC_EMR_Hydralazine_Deidentified.csv")

tBin = c("12AM-2AM", "2AM-4AM", "4AM-6AM", "6AM-8AM", "8AM-10AM", "10AM-12PM", "12PM-2PM", "2PM-4PM", "4PM-6PM", "6PM-8PM", "8PM-10PM", "10PM-12AM") 
frqs <- c("Scheduled", "As Needed", "Once")
doses <- c("First-dose", "All other doses")

dat <- input %>%
  #add a patient-order-dose identifier variable
  mutate(PT_ORD_DOS = paste(PATIENT, ORDER_NUM, DOSE_NUM, sep = "_")) %>%
  #first-doses only
  group_by(PT_ORD_DOS) %>% distinct(PT_ORD_DOS, .keep_all = TRUE) %>%
  #variable for first- vs. all other-doses
  mutate(DOSE_NUM_BIN = ifelse(DOSE_NUM == 1, "First-dose", "All other doses")) %>%
  mutate(DOSE_NUM_BIN = factor(DOSE_NUM_BIN, levels = doses)) %>%
  ungroup(PT_ORD_DOS) %>%
  #bin order types into "once", "as needed", or "scheduled"
  mutate(FREQ = case_when(
    grepl("as needed", FREQUENCY, ignore.case = TRUE) ~ "As Needed",
    grepl("once", FREQUENCY, ignore.case = TRUE) ~ "Once",
    TRUE ~ "Scheduled")) %>%
  mutate(FREQ = factor(FREQ, levels = frqs)) %>%
  mutate(DOSE_TAKEN = factor(DOSE_TAKEN, levels = tBin)) %>%
  select(PT_ORD_DOS, DOSE_NUM_BIN, FREQ, DOSE_TAKEN)

#compute percent of doses in each time bin, index time bins for cosinor regression
dos <- dat %>%
  group_by(FREQ, DOSE_NUM_BIN, DOSE_TAKEN) %>%
  summarise(count = n()) %>%
  mutate(PCT_DOS = round(count/sum(count), digits = 3)) %>%
  #remove two doses from "once" orders that are not first-doses
  filter(!(FREQ == "Once" & DOSE_NUM_BIN == "All other doses")) %>%
  #index times 
  mutate(index = seq(1, 12, by = 1)) %>%
  ungroup(FREQ, DOSE_NUM_BIN, DOSE_TAKEN)
#filter(FREQ != "Once")
dos <- data.frame(dos)

#dummy dataframe for cosinor regression
dum <- dos %>% filter(FREQ == "Once" & DOSE_NUM_BIN == "First-dose")
dum$DOSE_NUM_BIN <- "All other doses"
dum_bnd <- bind_rows(dos, dum)
dos <- data.frame(dum_bnd)
dos$FREQ <- factor(dos$FREQ, levels = frqs)
dos$DOSE_NUM_BIN <- factor(dos$DOSE_NUM_BIN, levels = doses)

#cosinor regression
fitD <- statsD <- NULL
for (i in 1:length(frqs)) {
  tmp <- dos %>% filter(FREQ == frqs[i])
  for (j in 1:length(doses))  {
    tmp2 <- tmp %>% filter(DOSE_NUM_BIN == doses[j])
    mod <- cosinor.lm(PCT_DOS ~ time(index), data = tmp2, period = 12) 
    fit <- data.frame(fit = mod$fit$fitted.values, index = tmp2$index) %>% mutate(DOSE_NUM_BIN = doses[j], FREQ = frqs[i])
    fitD <- bind_rows(fitD, fit)
    stats <- data.frame(pval = signif(summary(mod)$transformed.table$p.value[2], digits = 1), 
                        phase = abs(correct.acrophase(mod)),
                        rsq = signif(cosinor.PR(mod)[2], digits = 2), 
                        amp = signif(summary(mod)$transformed.table$estimate[2], digits = 2),
                        DOSE_NUM_BIN = doses[j], FREQ = frqs[i])
    statsD <- bind_rows(statsD, stats)
    }
  }

#add fit means to statsD, and add p val thresholding for labels
statsD <- statsD %>%
  mutate(pval = signif(pval, digits = 2)) %>%
  mutate(pval_cut = case_when(
    pval >= 0.05 ~ paste("p = ", as.character(pval), sep = ""),
    pval < 0.001 ~ "p < 0.001",
    pval < 0.05 ~ paste("p = ", as.character(pval), sep = ""))) %>%
  filter(!(FREQ == "Once" & DOSE_NUM_BIN == "All other doses"))

statsD$FREQ <- factor(statsD$FREQ, levels = frqs)
statsD$DOSE_NUM_BIN <- factor(statsD$DOSE_NUM_BIN, levels = doses)

fitD$FREQ <- factor(fitD$FREQ, levels = frqs)
fitD$DOSE_NUM_BIN <- factor(fitD$DOSE_NUM_BIN, levels = doses)  

#helper dataframes for plot labels
help <- dat %>%
  group_by(FREQ, DOSE_NUM_BIN) %>%
  count() %>%
  filter(!(FREQ == "Once" & DOSE_NUM_BIN == "All other doses"))

#max percent val for each bin
help2 <- dos %>%
  group_by(FREQ, DOSE_NUM_BIN) %>%
  filter(PCT_DOS == max(PCT_DOS))

##histogram first- and all other-doses by time of day for "once", "as needed", and "scheduled" orders
p <- ggplot(fitD %>% filter(!(FREQ == "Once" & DOSE_NUM_BIN == "All other doses"))) +
  geom_bar(aes(x = DOSE_TAKEN, y = PCT_DOS), data = dos %>% filter(!(FREQ == "Once" & DOSE_NUM_BIN == "All other doses")),
           stat = "identity", position = position_dodge(width = -0.35),
           color = "black", fill = "#999999") +
  geom_smooth(aes(x = index, y = fit), colour = "black", size = 0.5) +
  geom_text(aes(x = 3.25, y = 0.105), data = help, label = paste("n=", help$n), size = 2.5, color = "black") +
  geom_label(aes(x = 8, y = 0.04), data = statsD, label = statsD$pval_cut, size = 2.5, color = "black") +
  geom_label(aes(x = 8, y = 0.015), data = statsD, label = paste0("R^2 == ", statsD$Percent.rhythm, sep = ""), parse = TRUE,
             size = 2.5, color = "black") +
  scale_x_discrete(labels = c("12", "", "4", "", "8", "", "12", "", "4", "", "8", "")) +
  
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(cols = vars(FREQ), rows = vars(DOSE_NUM_BIN), scales = "free_y") +
  labs(y = "Percent of doses") +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 9, vjust = 1, hjust = 0.5),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 10, margin = margin(0, 0, 0, 0)))
ggsave(p, filename = "./figs/Hydralazine.DoseTiming.Different.OrderFreqs.pdf", width = 5, height = 3)
