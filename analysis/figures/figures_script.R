## This script produces the figures used in
## "Interleukin-6 is a mediator of therapeutic efficacy in acute respiratory 
## distress syndrome: an individual patient data meta-analysis of RCTs"

# Author: L Kramer
# Contact: l.kramer1@amsterdamumc.nl

library(tidyr)
library(dplyr)
library(ggplot2)

# read in direct, indirect, total effect estimates
res <- bind_rows(cbind(readRDS("ALVEOLI_res.rds"), 
                       study = c("ALVEOLI")),
                 cbind(readRDS("ARMA_res.rds"), 
                       study = c("ARMA")),
                 cbind(readRDS("SAVEMORE_res.rds"), 
                       study = c("SAVE-MORE")),
                 cbind(readRDS("COUNTERCOVID_res.rds"), 
                       study = c("COUNTER-COVID")),
                 cbind(readRDS("HARP2_res.rds"), 
                       study = c("HARP-2"))) %>% 
  mutate(class = if_else(is.na(class), "All", class))

# read in interaction effect estimates
beta_est <- bind_rows(cbind(readRDS("ALVEOLI_beta_est.rds"), 
                            study = c("ALVEOLI")),
                      cbind(readRDS("ARMA_beta_est.rds"), 
                            study = c("ARMA")),
                      cbind(readRDS("SAVEMORE_beta_est.rds"), 
                            study = c("SAVE-MORE")),
                      cbind(readRDS("COUNTERCOVID_beta_est.rds"), 
                            study = c("COUNTER-COVID")),
                      cbind(readRDS("HARP2_beta_est.rds"), 
                            study = c("HARP-2"))) %>% 
  mutate(class = if_else(is.na(class), "All", class))

# read in association estimates
alpha_est <- bind_rows(
  cbind(readRDS("ALVEOLI_alpha_est.rds"), 
        study = c("ALVEOLI")),
  cbind(readRDS("ARMA_alpha_est.rds"), 
        study = c("ARMA")),
  cbind(readRDS("SAVEMORE_alpha_est.rds"), 
        study = c("SAVE-MORE")),
  cbind(readRDS("COUNTERCOVID_alpha_est.rds"), 
        study = c("COUNTER-COVID")),
  cbind(readRDS("HARP2_alpha_est.rds"),
        study = c("HARP-2")))%>% 
  mutate(class = if_else(is.na(class), "All", class))

res <- res %>% group_by(study, endpoint, class)%>% 
  select(est, CI_lower, CI_upper, effect)%>% 
  mutate(id = case_when(study == "COUNTER-COVID" ~ "A",
                        study == "SAVE-MORE" ~ "B",
                        study == "ARMA" ~"C",
                        study == "ALVEOLI" & class == "All" ~"D",
                        study == "ALVEOLI" & class == "Hypo-inflammatory" ~"E",
                        study == "ALVEOLI" & class == "Hyper-inflammatory" ~"F",
                        study == "HARP-2" & class == "All" ~"G",
                        study == "HARP-2" & class == "Hypo-inflammatory" ~"H",
                        study == "HARP-2" & class == "Hyper-inflammatory" ~"I",
                        study == "Pooled" ~"J"))%>%
  mutate(id = factor(id, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")),
         effect = factor(effect),
         class = factor(class),
         endpoint = factor(endpoint))

## Figure 1 -------------------------------------------------------------------

# facet labels
facet.labs <- c(paste("Imatinib", sep = ""), 
                paste("Anakinra", sep = ""), 
                paste("Low TV", sep = ""), 
                paste("High PEEP", sep = ""),
                paste("Hypo-inflammatory", sep = ""),
                paste("Hyper-inflammatory", sep = ""),
                paste("Simvastatin", sep = ""),
                paste("Hypo-inflammatory", sep = ""),
                paste("Hyper-inflammatory", sep = ""), 
                "Pooled")

# for grey shades of hypo and hyper
a <- res %>% filter(effect == "total (Cox-PH)" & endpoint == "28-day endpoint") %>% 
  mutate(a =
           ifelse(class== "All", "black", "grey60")) %>% arrange(desc(id))

a <- a$a
names(facet.labs) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

p1 <- res%>% 
  filter(effect == "total (Cox-PH)" & endpoint == "28-day endpoint") %>%
  ggplot(aes(y = id,fill = class, color = class))+
  scale_color_manual(values=c("black","grey60", "grey60"))+
  theme_classic() +
  geom_point(aes(x=exp(est)), size=3, position=position_dodge(width = 1)) +
  geom_linerange(aes(xmin=exp(CI_lower), xmax=exp(CI_upper)), position=position_dodge(width = 1)) +
  geom_vline(xintercept = 1, linetype="dotted") +
  labs(x="HR, 95% CI", y= "Intervention", fill = "", col ="" )+
  guides(fill = "none", col = "none")+ 
  scale_x_continuous(trans = "log10", limits = c(0.2, 6), breaks = c(0.25, 0.5, 1, 2))+
  scale_y_discrete(labels = facet.labs, limits = rev)+
  #coord_fixed(ratio=.9)+
  ggtitle(expression("Total effect " *italic((Intervention %->% Survival)), ", HR (95%-CI)"))+
  theme(text = element_text(size = 20),
        axis.text.y = element_text(color = a),
        plot.title = element_text(size = 20,face = "bold"),
        plot.margin=unit(c(0.5,0.3, 0, 0.1), "cm") # TRouBLe
  )+
  geom_text(data = res %>% 
              filter(effect == "total (Cox-PH)" &
                       endpoint == "28-day endpoint"),
            aes(label = paste0(format(round(exp(est),2), 
                                      nsmall = 2),
                               " (", format(round(exp(CI_lower), 2),
                                            nsmall =2), ", ",
                               format(round(exp(CI_upper), 2), 
                                      nsmall = 2), ")"), x = 3.8),size = 3.5)+ 
  annotate(geom = 'segment', y = Inf, yend = Inf, x = 0, 
           color = "black", xend = Inf, size = 1) 


p2 <- res %>% 
  filter(effect == "direct" & endpoint == "28-day endpoint") %>%
  ggplot(aes(y = id,fill = class, color = class))+ 
  scale_color_manual(values=c("black","grey60", "grey60"))+
  theme_classic() +
  geom_point(aes(x=exp(est)), size=3, position=position_dodge(width = 1)) +
  geom_linerange(aes(xmin=exp(CI_lower), xmax=exp(CI_upper)), position=position_dodge(width = 1)) +
  geom_vline(xintercept = 1, linetype="dotted") +
  labs(x="HR, 95% CI", y= "", fill = "", col ="" )+
  guides(fill = "none", col = "none")+ 
  scale_x_continuous(trans = "log10", limits = c(0.23, 10), breaks = c(0.25, 0.5, 1, 2, 4))+
  scale_y_discrete(labels = facet.labs, limits = rev)+
  scale_shape_manual(values=c(5,17 ,18 ,19.5)) +
  ggtitle(expression("Direct effect " *italic((Intervention %->% Survival)), ", HR (95%-CI)"))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        plot.margin=unit(c(0.5,0.1, 0, -0.5), "cm"),
        text = element_text(size = 20)) +
  geom_text(data = res%>% 
              filter(effect == "direct" & endpoint == "28-day endpoint"),
            aes(label = paste0(format(round(exp(est),2), 
                                      nsmall = 2), " (",
                               format(round(exp(CI_lower), 2), 
                                      nsmall = 2), ", ",
                               format(round(exp(CI_upper), 2), 
                                      nsmall = 2), ")"), x = 8),
            size = 3.5)+ 
  annotate(geom = 'segment', y = Inf, yend = Inf, x = 0, 
           color = "black", xend = Inf, size = 1) 




p3 <- res %>% 
  filter(effect == "indirect" & endpoint == "28-day endpoint") %>%
  ggplot(aes(y = id,fill = class, color = class))+ 
  scale_color_manual(values=c("black","grey60",
                              "grey60"))+
  theme_classic() +
  geom_point(aes(x=round(exp(est),2)), size=3,
             position=position_dodge(width = 1)) +
  geom_linerange(aes(xmin=round(exp(CI_lower),2), 
                     xmax= round(exp(CI_upper),2)),
                 position=position_dodge(width = 1)) +
  geom_vline(xintercept = 1, linetype="dotted") +
  labs(x="HR, 95% CI", y= "", fill = "", col ="" )+
  guides(fill = "none", col = "none")+ 
  scale_x_continuous(trans = "log10",limits = c(0.849, 1.15), 
                     breaks = c(0.9, 1, 1.05))+
  scale_y_discrete(labels = facet.labs, limits = rev)+
  scale_shape_manual(values=c(5,17 ,18 ,19.5)) +
  ggtitle(expression("Indirect effect " *italic((Intervention %->% IL6 %->% Survival)), ", HR (95%-CI)"))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        plot.margin=unit(c(0.5,0.3, 0, 0.1), "cm"),
        text = element_text(size = 20)) +
  geom_text(data = res %>% 
              filter(effect == "indirect" & endpoint == "28-day endpoint"),
            aes(
              label = paste0(format(round(exp(est),2), nsmall = 2), " (",
                             format(round(exp(CI_lower), 2), nsmall = 2), ", ",
                             format(round(exp(CI_upper), 2), nsmall =2), ")"),
              x = 1.1), size = 3.5)+ 
  annotate(geom = 'segment', y = Inf, yend = Inf, x = 0, 
           color = "black", xend = Inf, size = 1) 


library(patchwork)

figure1 <- p1 + p2 + p3
figure1

#ggsave("./figure1.pdf", plot = figure1B, width = 500, height = 200, units = "mm")

## Figure 2 --------------------------------------------------------------------

### Panel A --------------------------------------------------------------------

# load long data frames with biomarker concentrations
countercovid_long <-readRDS("countercovid_long.rds")
savemore_long <-readRDS("savemore_long.rds") %>% filter(biomarker == "il6")
arma_long<-readRDS("arma_long.rds")
arma_surv <- readRDS("arma_surv.rds") 
alveoli_long<- readRDS("alveoli_long.rds")
alveoli_surv <- readRDS("alveoli_surv.rds") 
harp2_long<- readRDS("harp2_long.rds") %>% filter(biomarker == "IL_6")

library(labelled)
savemore_long <- savemore_long %>% remove_val_labels()

arma_long <- merge(arma_long, arma_surv[c("record.id", "death_d90", "time_mort90", "death_d28", "time_mort28")], by = "record.id") 

alveoli_long <- merge(alveoli_long, alveoli_surv[c("record.id", "death_d90", "time_mort90", "death_d28", "time_mort28")], by = "record.id") 

countercovid <- countercovid_long %>%select( c(record.id, randomized_group, day, conc_log10,death_d28, time_mort28)) %>%
  mutate(study = "COUNTER-COVID",
         day = as.numeric(day),
         class = "All",
         record.id = as.factor(record.id),
         death_d28 = as.factor(death_d28))
savemore <- savemore_long %>%select(c(record.id, randomized_group,day, conc_log10, death_d28, time_mort28))%>%
  mutate(study = "SAVE-MORE",
         class = "All",
         day = as.numeric(day),
         death_d28 = as.factor(death_d28))
arma <- arma_long %>%select(c(record.id,randomized_group, day, conc_log10,death_d28, time_mort28))%>%
  mutate(study ="ARMA",
         class = "All",
         day = as.numeric(day),
         record.id = as.factor(record.id),
         death_d28 = as.factor(death_d28)) 
alveoli <- alveoli_long %>%select( c(record.id,randomized_group, day, conc_log10,death_d28,time_mort28, class))%>%
  mutate(study = "ALVEOLI",
         class = case_when(class == "1" ~ "Hypo-inflammatory",
                           class == "2" ~ "Hyper-inflammatory"),
         day = as.numeric(day),
         record.id = as.factor(record.id),
         death_d28 = as.factor(death_d28))

harp2 <- harp2_long%>%
  select(c(record.id, randomized_group,day, conc_log10,death_d28, time_mort28, LCA_class))%>%
  mutate(study = "HARP-2",
         class = LCA_class,
         day = as.numeric(day),
         record.id = as.factor(record.id),
         death_d28 = as.factor(death_d28)) %>% select(!c(LCA_class))


savemore <- savemore %>% mutate(day = day - 1,
                                death_d28 = as.factor(death_d28))

data_all <- rbind(countercovid, savemore, arma, alveoli, harp2) %>% mutate(study = as.factor(study)) %>%
  filter(!is.na(conc_log10)  & !is.na(death_d28) & day<7) %>% 
  filter(day == 0 | day == 2 | day == 3 | day == 4)

# saveRDS(data_all, "data_all.rds")


# read in IL-6 concentrations in long format for all studies
# data_all <- readRDS("data_all.rds")

# plot deceased vs. alive IL-6 values
panelA <- data_all %>% filter(day == 0 | day == 3) %>%
  ggplot(aes(x = as.factor(day), y = 10^conc_log10, fill = death_d28))+
  geom_boxplot(outliers = FALSE)+
  theme_bw() +
  geom_jitter(width = .35, alpha = 0.05)+
  ggsci::scale_fill_lancet(alpha = 0.4, name = "Patient status after 28 days:", labels = c("Alive", "Deceased"))+
  scale_y_log10(breaks= scales::trans_breaks("log10", function(x) 10^x), labels =scales::trans_format("log10", scales::math_format(10^.x)))+
  theme_bw()+
  annotation_logticks(sides = "l")+
  #  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE)+
  labs(x = "Days from randomization", 
       y = "IL-6 concentration")+
  theme(text = element_text(size = 20),
        legend.position = c(0.75,0.9))+theme(text = element_text(size = 20))

panelA

### Panel B --------------------------------------------------------------------

# META ANALYZE THE ASSOCIATION ESTIMATE

set.seed(15)

# split by endpoint
alpha_meta_28 <- alpha_est %>% 
  filter(class == "All" & endpoint == "28-day endpoint")

alpha_meta_90 <- alpha_est %>% 
  filter(class == "All" & endpoint == "90-day endpoint")

#install.packages("metafor")
library(metafor)

## 28-day endpoint

# Meta-analysis with a random-effects model
meta_result28 <- rma(yi = alpha_meta_28$est, 
                     sei = alpha_meta_28$se, 
                     method = "REML")  


## 90-day endpoint 

# Meta-analysis with a random-effects model
meta_result90 <- rma(yi = alpha_meta_90$est, 
                     sei = alpha_meta_90$se, 
                     method = "REML")  

### add results to data
alpha_est <- alpha_est %>% rbind(data.frame(model = c("Association"),
                                            endpoint = c("28-day endpoint"),
                                            class = c("All"),
                                            est = c(meta_result28$beta),
                                            se = c(meta_result28$se),
                                            CI_lower = c(
                                              meta_result28$beta - 
                                                1.96*meta_result28$se),
                                            CI_upper = c(
                                              meta_result28$beta + 
                                                1.96*meta_result28$se),
                                            study = c("Pooled"))) 

alpha_est <- alpha_est %>% rbind(data.frame(model = c("Association"),
                                            endpoint = c("90-day endpoint"),
                                            class = c("All"),
                                            est = c(meta_result90$beta),
                                            se = c(meta_result90$se),
                                            CI_lower = c(
                                              meta_result90$beta - 
                                                1.96*meta_result90$se),
                                            CI_upper = c(
                                              meta_result90$beta + 
                                                1.96*meta_result90$se),
                                            study = c("Pooled")))

alpha_est$study <- factor(alpha_est$study, levels = 
                            c("Pooled", "HARP-2", "ALVEOLI", 
                              "ARMA", "SAVE-MORE", "COUNTER-COVID"))


alpha_est <- alpha_est %>% mutate(class = factor(class),
                                  endpoint =factor(endpoint))%>%
  group_by(study, endpoint, class) %>% 
  select(est, CI_lower, CI_upper) %>%
  mutate(id = case_when(study == "COUNTER-COVID" ~ "A",
                        study == "SAVE-MORE" ~ "B",
                        study == "ARMA" ~"C",
                        study == "ALVEOLI" & 
                          class == "All" ~"D",
                        study == "ALVEOLI" & 
                          class == "Hypo-inflammatory" ~"E",
                        study == "ALVEOLI" & 
                          class == "Hyper-inflammatory" ~"F",
                        study == "HARP-2" & class == "All" ~"G",
                        study == "HARP-2" & 
                          class == "Hypo-inflammatory" ~"H",
                        study == "HARP-2" & 
                          class == "Hyper-inflammatory" ~"I",
                        study == "Pooled" ~"J"))%>%
  mutate(id = factor(id, 
                     levels = 
                       c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")))

facet.labs <- c("COUNTER-COVID", "SAVE-MORE", "ARMA",
                "ALVEOLI", "Hypo-inflammatory", "Hyper-inflammatory",
                "HARP-2", "Hypo-inflammatory", "Hyper-inflammatory", 
                "Pooled")

names(facet.labs) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")

a <- alpha_est %>% filter(endpoint == "28-day endpoint") %>% 
  mutate(a =
           ifelse(class== "All", "black", "grey60")) %>% arrange(desc(id))

a <- a$a 

# log x axis:
pB <- alpha_est %>% filter(endpoint == "28-day endpoint") %>%
  ggplot(aes(y = id, fill = class, col = class))+
  scale_color_manual(values=c("black","grey60", "grey60"))+
  theme_bw() +
  geom_point(aes(x=exp(est)), size=3, position=position_dodge(width = 1)) +
  geom_linerange(aes(xmin=exp(CI_lower), xmax=exp(CI_upper)), 
                 position=position_dodge(width = 1)) +
  geom_vline(xintercept = 1, linetype="dotted") +
  labs(x="", y= "", fill = "", col ="" )+
  guides(fill = "none", col = "none")+ 
  scale_x_continuous(trans = "log10", breaks =c(1, 2, 8, 32))+
  scale_y_discrete(limits = rev,
                   labels= facet.labs)+
  theme(text = element_text(size = 10),
        axis.text.y = element_text(color = a))+#, angle = 45))+ 
  annotate(geom = 'segment', y = Inf, yend = Inf, x = 0, color = "black", 
           xend = Inf, size = 1) 

panelB <- pB+
  geom_text(data = alpha_est %>%
              filter(endpoint == "28-day endpoint"),
            aes(label = paste0(format(round(exp(est),2), nsmall = 2), " (", 
                               format(round(exp(CI_lower), 2), nsmall =2), ", ", format(round(exp(CI_upper), 2), nsmall = 2), ")"), x = exp(est), y =id), vjust = 2, size =6)+ 
  annotate(geom = 'segment', y = Inf, yend = Inf, x = 0,  color = "black", 
           xend = Inf, size = 1) +theme(text = element_text(size = 20))


# panelB
#ggsave("./association.pdf")


### Panel C --------------------------------------------------------------------

# Aim: plot for change in IL-6 between day 0 and 3 and mortality risk

# Strategy: predict mortality risk after day 3 using the joint models for each
# trial

## load fitted joint model objects
jointfit.savemore <- readRDS("./jointfit.savemore28.rds")
jointfit.countercovid <- readRDS("./jointfit_countercovid_28.rds")
jointfit.arma <- readRDS("./jointfit.arma28.rds")
jointfit.alveoli <- readRDS("./jointfit_alveoli_28.rds")
jointfit.harp2 <- readRDS("./jointfit_harp2_28_d3.rds")
jointfit.alveoli_hyper <- readRDS("./jointfit_alveoli_hyper_28.rds")
jointfit.alveoli_hypo <- readRDS("./jointfit_alveoli_hypo_28.rds")
jointfit.harp2_hyper <- readRDS("./jointfit_harp2_hyper_28_d3.rds")
jointfit.harp2_hypo <- readRDS("./jointfit_harp2_hypo_28_d3.rds")

# reshape all data as needed
data_all <- data_all %>% 
  mutate(
    time_mort28 = 3, # set last time of observation to 3
    death_d28 = 0) # pretend all are alive at day 3

data_all <- data_all %>% # unique subgroup identfiers
  mutate(id = case_when(study == "COUNTER-COVID" ~ "A",
                        study == "SAVE-MORE" ~ "B",
                        study == "ARMA" ~"C",
                        study == "ALVEOLI" & class == "All" ~"D", 
                        study == "ALVEOLI" & class == "Hyper-inflammatory" ~"E",
                        study == "ALVEOLI" & class == "Hypo-inflammatory" ~"F",
                        study == "HARP-2" & class == "All" ~ "G",
                        study == "HARP-2" & class == "Hyper-inflammatory" ~"H",
                        study == "HARP-2" & class == "Hypo-inflammatory" ~"I"))

data_all <- data_all[order(data_all$record.id, data_all$day), ] # sort by patient and day

## get joint model predictions

# countercovid 
pred_surv_cc <- predict(jointfit.countercovid,
                        newdata =data_all[data_all$study == "COUNTER-COVID",], 
                        process = "event", # gets the cumulative risk probabilities (CIFs)
                        times = 28,
                        return_newdata = T) %>% 
  filter(day == 28)%>% # save results for risk at day 28
  select(c(record.id, pred_CIF)) 

pred_surv_cc$id <- "A" # unique identifyer

# savemore
pred_surv_sm <- predict(jointfit.savemore,
                        newdata =data_all[data_all$study == "SAVE-MORE",], 
                        process = "event", 
                        times = 28,
                        return_newdata = T) %>%
  filter(day == 28)%>% 
  select(c(record.id, pred_CIF)) 

pred_surv_sm$id <- "B"

# arma 
pred_surv_ar <- predict(jointfit.arma, 
                        newdata =data_all[data_all$study == "ARMA",], 
                        process = "event",
                        times = 28,
                        return_newdata = T) %>%
  filter(day == 28)%>% 
  select(c(record.id, pred_CIF))

pred_surv_ar$id <- "C"

# alveoli all
pred_surv_al <- predict(jointfit.alveoli, 
                        newdata =data_all[data_all$study == "ALVEOLI",], 
                        process = "event", 
                        times = 28,
                        return_newdata = T) %>% 
  filter(day == 28)%>% 
  select(c(record.id, pred_CIF))

pred_surv_al$id <- "D"

# alveoli  hyper
pred_surv_al_y <- predict(jointfit.alveoli_hyper, 
                          newdata =data_all[c(
                            data_all$study == "ALVEOLI" & 
                              data_all$class == "Hyper-inflammatory"
                          ),], 
                          process = "event",
                          times = 28,
                          return_newdata = T) %>%
  filter(day == 28)%>% 
  select(c(record.id, pred_CIF))

pred_surv_al_y$id <- "E"


# alveoli  hypo
pred_surv_al_o <- predict(jointfit.alveoli_hypo, 
                          newdata =data_all[c(
                            data_all$study == "ALVEOLI" &
                              data_all$class == "Hypo-inflammatory"
                          ),], 
                          process = "event",  times = 28,
                          return_newdata = T) %>% 
  filter(day == 28)%>% 
  select(c(record.id, pred_CIF)) 

pred_surv_al_o$id <- "F"

# harp2 
pred_surv_ha <- predict(jointfit.harp2, 
                        newdata =data_all[data_all$study == "HARP-2",], 
                        process = "event", 
                        times = 28,
                        return_newdata = T) %>% 
  filter(day == 28)%>% 
  select(c(record.id, pred_CIF)) 

pred_surv_ha$id <- "G"

# harp2 hyper
pred_surv_ha_y <- predict(jointfit.harp2_hyper, 
                          newdata =data_all[c(
                            data_all$study == "HARP-2" & 
                              data_all$class == "Hyper-inflammatory"),] %>%
                            mutate(randomized_group = factor(
                              randomized_group,
                              levels = c("Placebo", "Simvastatin"))), 
                          process = "event", 
                          times = 28,
                          return_newdata = T) %>% 
  filter(day == 28)%>% 
  select(c(record.id, pred_CIF))

pred_surv_ha_y$id <- "H"

# harp2 hypo
pred_surv_ha_o <- predict(jointfit.harp2_hypo, 
                          newdata =data_all[c(
                            data_all$study == "HARP-2" & 
                              data_all$class == "Hypo-inflammatory"),], 
                          process = "event",
                          times = 28,
                          return_newdata = T) %>%
  filter(day == 28)%>% 
  select(c(record.id, pred_CIF)) 

pred_surv_ha_o$id <- "I"

# combine CIFs
CIF <- rbind(pred_surv_cc,pred_surv_sm, pred_surv_ar, pred_surv_al, pred_surv_al_y, pred_surv_al_o, pred_surv_ha, pred_surv_ha_y, pred_surv_ha_o)

data_all_CIF <-  CIF %>% left_join(data_all %>% select(-id), by = c("record.id"))

# saveRDS(data_all_CIF, "data_all_CIF.rds")

# data_all_CIF <- readRDS("data_all_CIF.rds")

## Next: get change in IL-6 between day 0 and day 3 for each patient where both measures are available

# pivot wider so each patient's `conc_log10` at each day is in a separate column
data_all_wide <- data_all_CIF %>% 
  group_by(record.id, day, death_d28, time_mort28, 
           study, pred_CIF, class, id, randomized_group) %>%
  summarise(conc_log10 = mean(conc_log10, na.rm = TRUE), .groups = "drop")%>%
  pivot_wider(
    names_from = day,
    values_from = conc_log10,
    names_prefix = "conc_log10_day"
  ) %>%
  mutate(id = as.factor(id))


# subset: all studies except COUNTER-COVID uses day 0 and 3
data_all_wide_b <- data_all_wide %>%
  filter(study != "COUNTER-COVID",
         conc_log10_day0 != "NULL",
         conc_log10_day3 != "NULL") %>%
  mutate(
    conc_log10_day0 = as.numeric(conc_log10_day0),
    conc_log10_day3 = as.numeric(conc_log10_day3),
    change_il6 = conc_log10_day3 - conc_log10_day0
  )

# subset: COUNTER-COVID (uses day 2)
data_all_wide_a <- data_all_wide %>%
  filter(study == "COUNTER-COVID",
         conc_log10_day0 != "NULL",
         conc_log10_day2 != "NULL") %>%
  mutate(
    conc_log10_day0 = as.numeric(conc_log10_day0),
    conc_log10_day2 = as.numeric(conc_log10_day2),
    change_il6 = conc_log10_day2 - conc_log10_day0
  )

# Combine back into one dataset
data_all_wide <- bind_rows(data_all_wide_b, data_all_wide_a)

library(ggsci)

# facet labels
facet.labs <- c(paste("COUNTER-COVID (Imatinib)", sep = ""), 
                paste("SAVE-MORE (Anakinra)", sep = ""), 
                paste("ARMA (Low VT)", sep = ""), 
                paste("ALVEOLI (High PEEP)", sep = ""),
                paste("ALVEOLI: hyper", sep = ""),
                paste("ALVEOLI: hypo", sep = ""),
                paste("HARP-2 (Simvastatin)", sep = ""),
                paste("HARP-2: hyper", sep = ""),
                paste("HARP-2: hypo", sep = ""))


names(facet.labs) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

# set intervention and control for all
data_all_wide <- data_all_wide %>%
  mutate(randomized_group = case_when(
    randomized_group == "Imatinib" ~ "Intervention",
    randomized_group == "Placebo" ~ "Control",
    randomized_group == "Anakinra" ~ "Intervention",
    randomized_group == "Randomized: 12 ml/kg" ~ "Control",
    randomized_group == "lower PEEP" ~ "Control",
    randomized_group == "Randomized: 6 ml/kg" ~ "Intervention",
    randomized_group == "higher PEEP" ~ "Intervention",
    randomized_group == "Simvastatin" ~ "Intervention"))

# make ggplot
panelC <- ggplot(data = data_all_wide,
                 aes(x = change_il6, y = pred_CIF
                 ))+
  geom_smooth(se = T)+
  geom_jitter(alpha = .05)+
  labs(x = "log10 change in IL-6 from baseline", 
       y = "28-day mortality risk")+
  geom_vline(xintercept =0, linetype = "dotted")+
  theme_bw() +
  facet_wrap(~id, labeller = labeller(id = facet.labs))+
  theme(text = element_text(size = 20))

#panelC


### Panel D --------------------------------------------------------------------

# data frame with total effects on y, time*treat on x

plotB_dat <- beta_est %>% mutate(time_x_treat = est) %>%
  left_join(res %>% 
              filter(effect == "total (Cox-PH)" & 
                       endpoint == "28-day endpoint") %>%
              mutate(total = est,
                     t_lower = CI_lower,
                     t_upper = CI_upper) %>%
              select(c(id, study, class, total, t_lower, t_upper)),
            by = c("study", "class"))

plotB_dat$label<- c("High PEEP","High PEEP (hypo-inflammatory)",
                    "High PEEP (hyper-inflammatory)", "Low TV",
                    "Anakinra","Imatinib","Simvastatin",
                    "Simvastatin (hypo-inflammatory)", 
                    "Simvastatin (hyper-inflammatory)") %>% as.factor()


# facet labels
facet.labs <- c(paste("Imatinib", sep = ""), 
                paste("Anakinra", sep = ""), 
                paste("Low VT", sep = ""), 
                paste("High PEEP", sep = ""),
                paste("High PEEP (hypo-inflammatory)", sep = ""),
                paste("High PEEP (hyper-inflammatory)", sep = ""),
                paste("Simvastatin", sep = ""),
                paste("Simvastatin (hypo-inflammatory)", sep = ""),
                paste("Simvastatin (hyper-inflammatory)", sep = ""))

names(facet.labs) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

panelD <- plotB_dat%>%
  ggplot(aes(x = as.numeric(time_x_treat), y = exp(total),  color = id))+
  geom_hline(yintercept = 1, linetype = "dotted", size = 1)+
  geom_vline(xintercept = 0, linetype = "dotted", size = 1)+
  geom_point(size = 3)+
  scale_y_continuous(limits = c(0, 2), oob = scales::rescale_none)+
  theme_bw()+ 
  #  annotate("text", x = -0.05, y = 1.8, label = paste0("r = ", round(cor_test$estimate, 2)), size = 6)+
  xlab("Intervention effect on IL-6 over time (95% CI)")+
  ylab("Total effect for 28-day endpoint (HR, 95% CI)")+
  geom_errorbar(aes(ymin = exp(t_lower), ymax = exp(t_upper)))+
  geom_errorbar(aes(xmin =CI_lower, xmax = CI_upper))+
  ggsci::scale_color_lancet(name = "Intervention:", labels =facet.labs)+
  annotate("segment", x =0, y =0.05,xend = -0.1, yend = 0.05, size = 7, 
           linejoin= "mitre",color = "darkgray",
           arrow = arrow(type = "closed", length = unit(0.02,"npc")))+
  annotate("text", x = -.05,y =.05,label = "Reduced IL-6", size = 4, 
           color = "white")+
  annotate("segment", x =-0.12, y =1,xend = -0.12, yend = 0.20, size = 7, 
           linejoin= "mitre", color = "darkgray",
           arrow = arrow(type = "closed", length = unit(0.01,"npc")))+
  annotate("text", x = -.12,y =.56,label = "Reduced mortality", size = 4, color = "white", angle = 90)+
  theme(legend.position = c(0.25, 0.75))+
  theme(text = element_text(size = 20))

#panelD


### Figure 2: combined ---------------------------------------------------------


library(patchwork)
library(ggtext)

# custom labels for axes

panelB <- panelB +
  labs(x = "Association between IL-6 and the hazard of death over 28
       days<br>(Joint model estimate per log<sub>10</sub> unit increase: 
       HR, 95%
       CI)</br>", size = 20)+ 
  theme(axis.title.x = element_markdown())

panelA <- panelA + 
  labs(y = "IL-6<br>log<sub>10</sub>(pg/ml)", size = 20)+ 
  theme(axis.title.y = element_markdown(), 
        text = element_text(size = 20))

panelC <- panelC +
  labs(x = "<br>log<sub>10</sub> IL-6 change within 3 days after randomization",
       y = "28-day mortality risk", size = 20)+ 
  theme(axis.title.x = element_markdown(), 
        text = element_text(size = 20))

panelD <- panelD + 
  labs(x = "Intervention effect on IL-6 per day\n(Mixed effects estimate:
       time*treat, 95% CI)",
       y = "Total effect over 28 days\n(Cox-PH estimate: HR, 95% CI)", 
       size =20)

# together
figure2 <- (
  (panelA + plot_spacer()+ panelB + plot_layout(widths = c(1, 0.05, 1))) /
    plot_spacer() /
    (panelC + plot_spacer()+ panelD + plot_layout(widths = c(1, 0.05, 1)))
) +
  plot_layout(heights = c(1, 0.1, 1)) +
  plot_annotation(tag_levels = "A") +
  theme(
    plot.tag = element_text(size = 20, face = "bold"), 
    text = element_text(size = 20) 
  )

figure2


## Supplemental figures --------------------------------------------------------

### Figure S3 ------------------------------------------------------------------

res_c <- readRDS("C:/Users/linak/Documents/lungs/lungs/res_c.rds")


# facet labels
facet.labs <- c(paste("Imatinib", sep = ""), 
                paste("Anakinra", sep = ""), 
                paste("Low TV", sep = ""), 
                paste("High PEEP", sep = ""),
                paste("Hypo-inflammatory", sep = ""),
                paste("Hyper-inflammatory", sep = ""),
                paste("Simvastatin", sep = ""),
                paste("Hypo-inflammatory", sep = ""),
                paste("Hyper-inflammatory", sep = ""))

res_c$LCA_class <- res_c$LCA

names(facet.labs) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

p1 <- res_c %>% filter(effect == "total (Cox-PH)" ) %>%
  ggplot(aes(y = id,fill = LCA_class, color = LCA_class))+
  scale_color_manual(values=c("black","grey60", "grey60"))+
  theme_classic() +
  geom_point(aes(x=exp(est)),size=3, position=position_dodge(width = 1)) +
  geom_linerange(aes(xmin=exp(CI_lower), xmax=exp(CI_upper)), position=position_dodge(width = 1)) +
  geom_vline(xintercept = 1, linetype="dotted") +
  labs(x="HR, 95% CI", y= "", fill = "Intervention", col ="" )+
  guides(fill = "none", col = "none")+ 
  scale_x_continuous(trans = "log2", limits = c(0.1, 6), breaks = c(0.25, 0.5, 1, 2))+
  scale_y_discrete(labels = facet.labs, limits = rev)+
  #coord_fixed(ratio=.9)+
  ggtitle("Total effect")+
  theme(text = element_text(size = 20),
        axis.text.y = element_text(color = a),
        plot.title = element_text(size = 20, face = "bold"),
        plot.margin=unit(c(0.5,0.3, 0, 0.1), "cm") # TRouBLe
  )+
  geom_text(data = res_c %>% filter(effect == "total (Cox-PH)"),
            aes(label = paste0(format(round(exp(est),2), nsmall = 2), " (", format(round(exp(CI_lower), 2), nsmall =2), ", ", format(round(exp(CI_upper), 2), nsmall = 2), ")"), x = 3.8), size = 3.5)+ 
  annotate(geom = 'segment', y = Inf, yend = Inf, x = 0, color = "black", xend = Inf, size = 1) 



p2 <- res_c %>% filter(effect == "direct") %>%
  ggplot(aes(y = id,fill = LCA_class, color = LCA_class))+ 
  scale_color_manual(values=c("black","grey60", "grey60"))+
  theme_classic() +
  geom_point(aes(x=exp(est)), size=3, position=position_dodge(width = 1)) +
  geom_linerange(aes(xmin=exp(CI_lower), xmax=exp(CI_upper)), position=position_dodge(width = 1)) +
  geom_vline(xintercept = 1, linetype="dotted") +
  labs(x="HR, 95% CI", y= "", fill = "", col ="" )+
  guides(fill = "none", col = "none")+ 
  scale_x_continuous(trans = "log2", limits = c(0.2, 8), breaks = c(0.25, 0.5, 1, 2, 4))+
  scale_y_discrete(labels = facet.labs, limits = rev)+
  #    coord_fixed(ratio=.6)+
  scale_shape_manual(values=c(5,17 ,18 ,19.5)) +
  ggtitle("Direct effect")+
  # theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        plot.margin=unit(c(0.5,0.1, 0, -0.5), "cm"),
        text = element_text(size = 20)) + 
  geom_text(data = res_c %>% filter(effect == "direct"),
            aes(
              label = paste0(format(round(exp(est),2), nsmall = 2), " (", format(round(exp(CI_lower), 2), nsmall = 2), ", ", format(round(exp(CI_upper), 2), nsmall =2), ")"), x = 6), size = 3.5)+ 
  annotate(geom = 'segment', y = Inf, yend = Inf, x = 0, color = "black", xend = Inf, size = 1) 



p3 <- res_c %>% filter(effect == "indirect") %>%
  ggplot(aes(y = id,fill = LCA_class, color = LCA_class))+ 
  scale_color_manual(values=c("black","grey60",
                              "grey60"))+
  theme_classic() +
  geom_point(aes(x=round(exp(est),2)), size=3, position=position_dodge(width = 1)) +
  geom_linerange(aes(xmin=round(exp(CI_lower),2), xmax= round(exp(CI_upper),2)), position=position_dodge(width = 1)) +
  geom_vline(xintercept = 1, linetype="dotted") +
  labs(x="HR, 95% CI", y= "", fill = "", col ="" )+
  guides(fill = "none", col = "none")+ 
  scale_x_continuous(trans = "log2",limits = c(0.86, 1.15), breaks = c(0.9, 1, 1.05))+
  scale_y_discrete(labels = facet.labs, limits = rev)+
  #    coord_fixed(ratio=.5)+
  scale_shape_manual(values=c(5,17 ,18 ,19.5)) +
  ggtitle("Indirect effect")+
  theme(axis.text.y = element_blank(),
        axis.ticks.y= element_blank(),
        axis.line.y = element_blank(),
        plot.title = element_text(size = 20, face = "bold"),
        #   axis.line.y = element_blank(),
        plot.margin=unit(c(0.5,0.3, 0, 0.1), "cm"),
        text = element_text(size = 20)) +
  geom_text(data = res_c %>% filter(effect == "indirect"),
            aes(
              label = paste0(format(round(exp(est),2), nsmall = 2), " (", format(round(exp(CI_lower), 2), nsmall = 2), ", ", format(round(exp(CI_upper), 2), nsmall =2), ")"), x = 1.1), size = 3.5)+ 
  annotate(geom = 'segment', y = Inf, yend = Inf, x = 0, color = "black", xend = Inf, size = 1) 



figure3S <- p1 + p2 + p3

figure3S
