## Script to analyse all datasets together combined
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

## Load packages ----
library(tidyverse)
library(plyr)
library(emmeans) 
library(papaja)

## Source functions ----
source("funcs/theme.R")
source("funcs/rainclouds.R")
source("dataset_a.R")
source("dataset_b.R")
source("dataset_c.R")

## Add NAs for dataset C; does not have finger tapping data
dataset_c$CONS_change <- rep(NA, nrow(dataset_c) )

## Get colnumnnames to match on
colnames_a <- colnames(dataset_a)
colnames_b <- colnames(dataset_b)
colnames_c <- colnames(dataset_c)

common_colnames <- intersect(colnames_a,  colnames_b)
common_colnames <- intersect(common_colnames,  colnames_c)

## Filter only the columns that are common
dataset_a <- dataset_a %>% dplyr::select(common_colnames)
dataset_b <- dataset_b %>% dplyr::select(common_colnames)
dataset_c <- dataset_c %>% dplyr::select(common_colnames)

dataset_a$dataset <- "A"
dataset_b$dataset <- "B"
dataset_c$dataset <- "C"

dataset_a_comb <- dataset_a
dataset_b_comb <- dataset_b
dataset_c_comb <- dataset_c

## Add power Delta data
select_vars <- c("1","1.2","1.4","1.6","1.8",
                 "2","2.2","2.4","2.6","2.8",
                 "3","3.2","3.4","3.6","3.8",
                 "4","4.2","4.4","4.6","4.8")
mergedData_wide_dataset_a <- mergedData_wide_dataset_a %>% mutate(av_delta = rowMeans(dplyr::select(., select_vars)))
mergedData_wide_dataset_b <- mergedData_wide_dataset_b %>% mutate(av_delta = rowMeans(dplyr::select(., select_vars)))
mergedData_wide_dataset_c <- mergedData_wide_dataset_c %>% mutate(av_delta = rowMeans(dplyr::select(., select_vars)))

pow_a  <- mergedData_wide_dataset_a %>% dplyr::select(pnum, datasetnum, av_delta)
pow_b  <- mergedData_wide_dataset_b %>% dplyr::select(pnum, datasetnum, av_delta)
pow_c  <- mergedData_wide_dataset_c %>% dplyr::select(pnum, datasetnum, av_delta)

dataset_a_comb <- merge(dataset_a_comb,pow_a, by = c("pnum", "datasetnum"))
dataset_b_comb <- merge(dataset_b_comb,pow_b, by = c("pnum", "datasetnum"))
dataset_c_comb <- merge(dataset_c_comb,pow_c, by = c("pnum", "datasetnum"))

## Rename to combine 
dataset_a_comb$group<- revalue(dataset_a_comb$group, c("Patients" = "Medicated"))
dataset_b_comb$group<- revalue(dataset_b_comb$group, c("control" = "Controls",
                                             "patient_unmed" = "Unmedicated",
                                             "patient_7d" = "Medicated"))
dataset_c_comb$group<- revalue(dataset_c_comb$group, c("control" = "Controls",
                                             "patient_28" = "Medicated_long",
                                             "patient_7" = "Medicated"))

dataset_comb <- rbind(dataset_a_comb, dataset_b_comb, dataset_c_comb)

## This leaves out Dataset C 28d repeated patients; change to create other combinations
dataset_comb <- dataset_comb %>% filter(group == "Controls" |  group == "Medicated") 

## Summary statistics and ttests ----
## Sleep architecture ----
dataset_comb_hypno_long <- dataset_comb %>%
  gather(sleepstage, score, S1_percent,S2_percent,S3andS4_perc,REM_percent,NREM_percent,
         Wake_after_sleep_onset_percent,Total_sleep_time_min, Sleep_Onset_min, SWS_onset_min, REM_onset_min)

dataset_comb_hypno_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

hypnogram_ttests_perc_dataset_comb <- lapply(dataset_comb[,c("age", "S1_percent", "S2_percent", "S3andS4_perc",  
                                                             "REM_percent", "NREM_percent", "Wake_after_sleep_onset_percent", 
                                                             "Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min" 
                                                             )], function(x) t.test(x ~ dataset_comb$group, var.equal = FALSE))
hypnogram_S4_ttests_perc_dataset_comb <- t.test(dataset_comb$S4_percent ~ dataset_comb$group) 

lm_waso_dataset <-(lm(Wake_after_sleep_onset_percent~group*dataset, dataset_comb))
emmeans(lm_waso_dataset, specs = pairwise ~ dataset)


## Spindle data ----
dataset_comb_spindle_long <- dataset_comb %>%
  gather(spindle, score,  spin_dens, spin_count,spin_amp, spin_dur, spin_freq)

dataset_comb_spindle_long %>%
  dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%
  get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_comb <- lapply(dataset_comb[,c("spin_dens","spin_dur",
                                                      "spin_amp", "spin_count",  "spin_freq" 
                                                      )], function(x) t.test(x ~ dataset_comb$group, var.equal = FALSE))

## SO data ----
dataset_comb_so_long <- dataset_comb %>%
  gather(so, score,  so_dens,so_amp, so_dur, so_freq, so_count)

dataset_comb_so_long %>%
  dplyr::select(so, score, group) %>%
  group_by(group, so) %>%
  get_summary_stats(score, type = "mean_se")

so_ttests_dataset_comb <- lapply(dataset_comb[,c("so_dens","so_dur","so_amp" ,
                                                 "so_freq", "so_count"
                                                 )], function(x) t.test(x ~ dataset_comb$group, var.equal = FALSE))

## Coupling data ----
dataset_comb_coupling_long <- dataset_comb %>%
  gather(coupling, score,  coupling_occur,mean_delay,sd_delay)

dataset_comb_coupling_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_ttests_dataset_comb <- lapply(dataset_comb[,c("coupling_occur","mean_delay","sd_delay" 
                                                       )], function(x) t.test(x ~ dataset_comb$group, var.equal = FALSE))

bf_meandelay_dataset_comb  <- ttestBF(formula = mean_delay ~ group, data = dataset_comb)
bf_sddelay_dataset_comb  <- ttestBF(formula = sd_delay ~ group, data = dataset_comb)


## Plots ----
## Spindles 
plot_rain_spinamp_dataset_comb <- ggplot(dataset_comb,aes(x=group,y=spin_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle amplitude')+ xlab('')+ ylim(0,100) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_comb) + 
  scale_fill_manual(values=color_dataset_comb) + 
  scale_x_discrete(labels = c("Controls (comb)", "Medicated (comb)")) +
  ggtitle("Group difference in spindle amplitude")

plot_rain_spidens_dataset_comb <- ggplot(dataset_comb,aes(x=group,y=spin_dens,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_dens),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle density')+ xlab('')+ ylim(0,4) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_comb) + 
  scale_fill_manual(values=color_dataset_comb) + 
  scale_x_discrete(labels = c("Controls (comb)", "Medicated (comb)")) +
  ggtitle("Group difference in spindle density") 


## SOs
plot_rain_soamp_dataset_comb <- ggplot(dataset_comb,aes(x=group,y=so_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO amplitude')+ xlab('')+ ylim(0,300) +
  raincloud_theme + my_theme + 
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_comb) + 
  scale_fill_manual(values=color_dataset_comb) + 
  scale_x_discrete(labels = c("Controls (comb)", "Medicated (comb)")) +
  ggtitle("Group difference in SO amplitude") +
  geom_signif(annotation=formatC("***", digits=1), y_position=280,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

plot_rain_sodur_dataset_comb <- ggplot(dataset_comb,aes(x=group,y=so_dur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_dur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO duration')+ xlab('')+ ylim(1,1.6) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_comb) + 
  scale_fill_manual(values=color_dataset_comb) +
  scale_x_discrete(labels = c("Controls (comb)", "Medicated (comb)")) +
  ggtitle("Group difference in SO duration")+
  geom_signif(annotation=formatC("***", digits=1), y_position=1.55,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

## Coupling
plot_rain_amount_coup_dataset_comb <- ggplot(dataset_comb,aes(x=group,y=coupling_occur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=coupling_occur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO-fast spindles')+ xlab('')+ ylim(0,320) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_comb) + 
  scale_fill_manual(values=color_dataset_comb) + 
  scale_x_discrete(labels = c("Controls (comb)", "Medicated (comb)")) +
  ggtitle("Group difference in amount of SO-fast spindle couplings")+
  geom_signif(annotation=formatC("*", digits=1), y_position=300,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

plot_rain_mean_delay_dataset_comb <- ggplot(dataset_comb,aes(x=group,y=mean_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=mean_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Mean delay')+ xlab('')+ ylim(0.3,0.85) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_comb) + 
  scale_fill_manual(values=color_dataset_comb) + 
  scale_x_discrete(labels = c("Controls (comb)", "Medicated (comb)")) +
  ggtitle("Group difference in mean delay") 

plot_rain_sd_delay_dataset_comb <- ggplot(dataset_comb,aes(x=group,y=sd_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=sd_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Delay dispersion (SD)')+ xlab('')+ ylim(0.15,0.42) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_comb) + 
  scale_fill_manual(values=color_dataset_comb) + 
  scale_x_discrete(labels = c("Controls (comb)", "Medicated (comb)")) +
  ggtitle("Group difference in delay dispersion (SD)")+
  geom_signif(annotation=formatC("***", digits=1), y_position=0.4,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

## Add the plots on top of eachother
spin_plot_dataset_comb     <- ggarrange(plot_rain_spidens_dataset_comb, plot_rain_spinamp_dataset_comb,nrow = 2, labels = c("A", "B") )
so_plot_dataset_comb       <- ggarrange(plot_rain_soamp_dataset_comb, plot_rain_sodur_dataset_comb,nrow = 2, labels = c("A", "B") )
coupling_plot_dataset_comb <- ggarrange(plot_rain_amount_coup_dataset_comb, plot_rain_mean_delay_dataset_comb, plot_rain_sd_delay_dataset_comb, nrow = 3, labels = c("A", "B", "C") )

## Save the images
# ggsave("figs/spin_plot_dataset_comb_v2.png", spin_plot_dataset_comb, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/so_plot_dataset_comb_v2.png", so_plot_dataset_comb, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/coupling_plot_dataset_comb_v2.png", coupling_plot_dataset_comb, device = "png", dpi = 300, width = 7, height = 10)

## Depression severity ----
## Take Hamd_week 0 for dataset A
## Take Hamd week 1 for Dataset B
## Take Hamd week 1 for Dataset C
dataset_a_hamd_0     <- dataset_comb %>% filter(dataset=="A") %>%dplyr::select(hamd_0)
dataset_b_hamd_week1 <- dataset_comb %>% filter(dataset=="B") %>%dplyr::select(hamd_week1)
dataset_c_hamd_week1 <- dataset_comb %>% filter(dataset=="C") %>%dplyr::select(hamd_week1)

colnames(dataset_a_hamd_0)     <- "hamd_comb"
colnames(dataset_b_hamd_week1) <- "hamd_comb"
colnames(dataset_c_hamd_week1) <- "hamd_comb"
hamd_comb    <- rbind(dataset_a_hamd_0, dataset_b_hamd_week1, dataset_c_hamd_week1)
dataset_comb <- cbind(dataset_comb, hamd_comb)

dataset_comb_pat<- dataset_comb %>% filter(group == "Medicated")
  
## Depression severity and sleep parameters ----
lm_hamds_spindens_all  <- summary(lm(spin_dens ~ hamd_comb, dataset_comb))
lm_hamds_spincount_all <- summary(lm(spin_count ~ hamd_comb, dataset_comb))
lm_hamds_so_amp_all    <- summary(lm(so_amp ~ hamd_comb, dataset_comb))
lm_hamds_so_dur_all    <- summary(lm(so_dur ~ hamd_comb, dataset_comb))
lm_hamds_sd_delay_all  <- summary(lm(sd_delay ~ hamd_comb, dataset_comb))
lm_hamds_cons_all      <- summary(lm(CONS_change ~ hamd_comb, dataset_comb))

lm_hamd_dataset_all <- lapply(dataset_comb[,c("spin_dens","spin_count","spin_amp", 
                                              "spin_dur", "spin_freq","so_dens","so_dur",
                                              "so_amp", "so_count", "so_freq" , "coupling_occur",
                                              "mean_delay","sd_delay")], function(x) summary(lm(x ~ dataset_comb$hamd_comb)))

## Correlation table of depression severity, number of episodes and response outcome
dataset_comb_controls <- dataset_comb  %>% dplyr::filter(group == "Controls")
dataset_comb_patients <- dataset_comb  %>% dplyr::filter(group == "Medicated")

lm_hamd_dataset_all_corr_pat <- lapply(dataset_comb_patients[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq","so_dens","so_dur","so_amp", "so_count", "so_freq" , "coupling_occur","mean_delay","sd_delay", "av_delta")], function(x) cor.test(x, dataset_comb_patients$hamd_comb))
lm_eps_dataset_all_corr_pat  <- lapply(dataset_comb_patients[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq","so_dens","so_dur","so_amp", "so_count", "so_freq" , "coupling_occur","mean_delay","sd_delay", "av_delta")], function(x) cor.test(x, dataset_comb_patients$no_episodes))
lm_resp_dataset_all_corr_pat <- lapply(dataset_comb_patients[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq","so_dens","so_dur","so_amp", "so_count", "so_freq" , "coupling_occur","mean_delay","sd_delay", "av_delta")], function(x) cor.test(x, dataset_comb_patients$responder))

## Moderation effect of age
dataset_comb$hamd_comb_center <- scale(dataset_comb$hamd_comb, scale = F)
dataset_comb$age_center <- scale(dataset_comb$age, scale = F)

dataset_comb$age_center_high <- dataset_comb$age_center + sd(dataset_comb$age_center)
dataset_comb$age_center_low  <- dataset_comb$age_center - sd(dataset_comb$age_center)

lm_hamds_spindens_all_high    <- summary(lm(spin_dens ~ hamd_comb*age_center_high, dataset_comb))
lm_hamds_spindens_all_low     <- summary(lm(spin_dens ~ hamd_comb*age_center_low, dataset_comb))
lm_hamds_spindens_all_center  <- summary(lm(spin_dens ~ hamd_comb*age_center, dataset_comb))

lm_hamds_spincount_all <- summary(lm(spin_count ~ hamd_comb*age, dataset_comb))
lm_hamds_age_all       <- summary(lm(hamd_comb ~ age, dataset_comb))

## depression severity and age
lm_hamds_spin_dens_age_all  <- summary(lm(spin_dens ~ age*hamd_comb, dataset_comb))
lm_hamds_spincount_all <- summary(lm(CONS_change ~ no_episodes, dataset_comb))
lm_hamds_cons_all      <- summary(lm(CONS_change ~ hamd_comb, dataset_comb))

## Number of episodes and sleep parameters
lm_no_episodes_spifreq_all   <- (lm(spin_freq ~ no_episodes, dataset_comb))
lm_no_episodes_spindens_all  <- (lm(spin_dens ~ no_episodes, dataset_comb))
lm_no_episodes_spincount_all <- (lm(spin_count ~ no_episodes, dataset_comb))
lm_no_episodes_spinamp_all   <- (lm(spin_amp ~ no_episodes, dataset_comb))
lm_no_episodes_so_amp_all    <- (lm(so_amp ~ no_episodes, dataset_comb))
lm_no_episodes_so_dur_all    <- (lm(so_dur ~ no_episodes, dataset_comb))
lm_no_episodes_sd_delay_all  <- (lm(sd_delay ~ no_episodes, dataset_comb))
lm_no_episodes_cons_all      <- (lm(CONS_change ~ no_episodes, dataset_comb))
lm_no_episodes_age_cons_all  <- (lm(no_episodes ~ age, dataset_comb))

cor.test(dataset_comb$spin_dens, dataset_comb$no_episodes)
cor.test(dataset_comb$spin_count, dataset_comb$no_episodes)
cor.test(dataset_comb$spin_amp, dataset_comb$no_episodes)
cor.test(dataset_comb$spin_freq, dataset_comb$no_episodes)
cor.test(dataset_comb$spin_dur, dataset_comb$no_episodes)
cor.test(dataset_comb$so_dens, dataset_comb$no_episodes)
cor.test(dataset_comb$so_count, dataset_comb$no_episodes)
cor.test(dataset_comb$so_amp, dataset_comb$no_episodes)
cor.test(dataset_comb$so_freq, dataset_comb$no_episodes)
cor.test(dataset_comb$so_dur, dataset_comb$no_episodes)
cor.test(dataset_comb$coupling_occur, dataset_comb$no_episodes)
cor.test(dataset_comb$mean_delay, dataset_comb$no_episodes)
cor.test(dataset_comb$sd_delay, dataset_comb$no_episodes)

dataset_comb_no_episodes_bf <- dataset_comb_pat %>% drop_na(no_episodes)
bf_lm_noepi_spindens_comb = regressionBF(spin_dens ~ no_episodes, data = dataset_comb_no_episodes_bf)
bf_lm_noepi_so_amp_comb = regressionBF(so_amp ~ no_episodes, data = dataset_comb_no_episodes_bf)

lm_no_episodes_age_all  <- summary(lm(no_episodes ~ age, dataset_comb))
lm_hamds_spindens_all  <- summary(lm(spin_dens ~ no_episodes*age, dataset_comb_pat))

## Responder / non-responder (Not available in Dataset A)
dataset_comb_pat_resp <- dataset_comb_pat %>% filter(dataset != "A")

## Outcome response and sleep parameters
lm_resp_rem_all       <- summary(lm( responder ~ REM_percent, dataset_comb_pat_resp))
lm_resp_spindens_all  <- summary(lm( responder ~ spin_dens, dataset_comb_pat_resp))
lm_resp_spincount_all <- summary(lm( responder ~ spin_count, dataset_comb_pat_resp))
lm_resp_so_amp_all    <- summary(lm( responder ~ so_amp, dataset_comb_pat_resp))
lm_resp_so_dur_all    <- summary(lm( responder ~ so_dur, dataset_comb_pat_resp))
lm_resp_sd_delay_all  <- summary(lm( responder ~ sd_delay, dataset_comb_pat_resp))
lm_resp_sd_delay_all  <- summary(lm( responder ~ coupling_occur, dataset_comb_pat_resp))
lm_resp_cons_all      <- summary(lm( responder ~ CONS_change, dataset_comb_pat_resp))
lm_resp_rem_all       <- summary(lm( responder ~ hamd_comb, dataset_comb_pat_resp))

## Effects of medication ----
dataset_b_med <- dataset_comb_pat %>% pivot_wider(names_from = med_group, values_from = med_group)
dataset_a_med <- read.xlsx("dataset_a_medication.xlsx", sheetIndex = 1)
dataset_c_med <- read.xlsx("dataset_c_medication.xlsx", sheetIndex = 1)

dataset_b_med <- dataset_b_med %>%  
 dplyr::filter(dataset == "B" & group == "Medicated" ) %>% 
 dplyr::select(data_name,datasetnum,pnum,group,med_type, SSRI, SNRI, TCA, NaSSA, NDRI, NaSSRI, NARI)

## Get counts of type
dataset_b_med %<>% mutate_if(is.character,as.factor)
dataset_b_med$SSRI   <- recode( dataset_b_med$SSRI,  SSRI = "1")
dataset_b_med$SNRI   <- recode( dataset_b_med$SNRI,  SNRI = "1")
dataset_b_med$TCA    <- recode( dataset_b_med$TCA,   TCA = "1")
dataset_b_med$NaSSA  <- recode( dataset_b_med$NaSSA, NaSSA = "1")
dataset_b_med$NDRI   <- recode( dataset_b_med$NDRI,  NDRI = "1")
dataset_b_med$NaSSRI <- recode( dataset_b_med$NaSSRI,NaSSRI = "1")
dataset_b_med$NARI   <- recode( dataset_b_med$NARI,  NARI = "1")
dataset_b_med$Rest   <- dataset_b_med$NDRI
dataset_b_med$sleep  <- NA
dataset_b_med %<>% mutate_at(vars(SSRI:Rest), funs(as.numeric))
dataset_b_med[is.na(dataset_b_med)] <- 0
dataset_b_med$antidepressant <- dataset_b_med$med_type

dataset_a_med  <- dataset_a_med %>% filter(Group == "Patients")
dataset_c_med  <- dataset_c_med %>% filter(group == "patient_7")

dataset_a_med <-   dataset_a_med%>% 
  dplyr::select(data_name, antidepressant,SSRI, SNRI, TCA, NaSSRI, Rest, sleep)
dataset_b_med <-   dataset_b_med%>% 
 dplyr::select(data_name,antidepressant, SSRI, SNRI, TCA, NaSSRI, Rest, sleep)
dataset_c_med <-   dataset_c_med%>% 
  dplyr::select(data_name, antidepressant, SSRI, SNRI, TCA, NaSSRI, Rest, sleep)

dataset_a_med$dataset <-"A"
dataset_b_med$dataset <-"B"
dataset_c_med$dataset <-"C"

dataset_comb_med <-  rbind(dataset_a_med, dataset_b_med, dataset_c_med)
dataset_comb_med$antidepressant <- str_to_lower(dataset_comb_med$antidepressant)

## What are the most common?
names(sort(summary(as.factor(dataset_comb_med$antidepressant)), decreasing=T)[1:6])
length(grep("venlafaxin", dataset_comb_med$antidepressant)) # 24
length(grep("mirtazapin", dataset_comb_med$antidepressant)) # 19
length(grep("citalopram", dataset_comb_med$antidepressant))
length(grep("trimipramin", dataset_comb_med$antidepressant)) #18

dataset_comb_med$venlafaxin <- 0
dataset_comb_med[grep("venlafaxin", dataset_comb_med$antidepressant),]$venlafaxin <- 1
dataset_comb_med$mirtazapin <- 0
dataset_comb_med[grep("mirtazapin", dataset_comb_med$antidepressant),]$mirtazapin <- 1
dataset_comb_med$trimipramin <- 0
dataset_comb_med[grep("trimipramin", dataset_comb_med$antidepressant),]$trimipramin <- 1
dataset_comb_med$citalopram <- 0
dataset_comb_med[grep("citalopram", dataset_comb_med$antidepressant),]$citalopram <- 1

## merge and only take the medicated patients
dataset_comb_med_all <- merge(dataset_comb, dataset_comb_med, by = "data_name")
dataset_comb_med_all<- dataset_comb_med_all%>% filter(group ==  "Medicated")

dataset_comb_med_all_na <- dataset_comb_med_all
dataset_comb_med_all_na <- dataset_comb_med_all_na %>% drop_na(SSRI)

## statistics on medication type
summary(lm(spin_dens ~ SSRI, dataset_comb_med_all))
summary(lm(spin_count ~ SSRI, dataset_comb_med_all))
summary(lm(so_amp ~ SSRI, dataset_comb_med_all))
summary(lm(so_dur ~ SSRI, dataset_comb_med_all))
summary(lm(NREM_percent ~ SSRI, dataset_comb_med_all))
summary(lm(coupling_occur ~ SSRI, dataset_comb_med_all))
summary(lm(mean_delay ~ SSRI, dataset_comb_med_all))
summary(lm(sd_delay ~ SSRI, dataset_comb_med_all))
summary(lm(CONS_change ~ SSRI, dataset_comb_med_all))
summary(lm(spin_dens ~ SNRI, dataset_comb_med_all))
summary(lm(spin_count ~ SNRI, dataset_comb_med_all))
summary(lm(so_amp ~ SNRI, dataset_comb_med_all))
summary(lm(so_dur ~ SNRI, dataset_comb_med_all))
summary(lm(NREM_percent ~ SNRI, dataset_comb_med_all))
summary(lm(REM_percent ~ SNRI, dataset_comb_med_all))
summary(lm(coupling_occur ~ SNRI, dataset_comb_med_all))
summary(lm(mean_delay ~ SNRI, dataset_comb_med_all))
summary(lm(sd_delay ~ SNRI, dataset_comb_med_all))
summary(lm(CONS_change ~ SNRI, dataset_comb_med_all))

summary(lm(spin_dens ~ sleep, dataset_comb_med_all))
summary(lm(spin_count ~ sleep, dataset_comb_med_all))
summary(lm(spin_amp ~ sleep, dataset_comb_med_all))
summary(lm(so_amp ~ sleep, dataset_comb_med_all)) #
summary(lm(so_dur ~ sleep, dataset_comb_med_all))
summary(lm(coupling_occur ~ sleep, dataset_comb_med_all))
summary(lm(mean_delay ~ sleep, dataset_comb_med_all))#
summary(lm(sd_delay ~ sleep, dataset_comb_med_all))
summary(lm(NREM_percent ~ sleep, dataset_comb_med_all))
summary(lm(CONS_change ~ sleep, dataset_comb_med_all))

summary(lm(spin_dens ~ TCA, dataset_comb_med_all))
summary(lm(spin_count ~ TCA, dataset_comb_med_all))
summary(lm(so_amp ~ TCA, dataset_comb_med_all))#
summary(lm(so_dur ~ TCA, dataset_comb_med_all))
summary(lm(coupling_occur ~ TCA, dataset_comb_med_all))
summary(lm(mean_delay ~ TCA, dataset_comb_med_all))
summary(lm(sd_delay ~ TCA, dataset_comb_med_all))
summary(lm(NREM_percent ~ TCA, dataset_comb_med_all))
summary(lm(CONS_change ~ TCA, dataset_comb_med_all))

## count as rest
summary(lm(spin_dens ~ NaSSRI, dataset_comb_med_all))  #
summary(lm(spin_count ~ NaSSRI, dataset_comb_med_all))#
summary(lm(so_amp ~ NaSSRI, dataset_comb_med_all))
summary(lm(so_dur ~ NaSSRI, dataset_comb_med_all))#
summary(lm(coupling_occur ~ NaSSRI, dataset_comb_med_all))
summary(lm(mean_delay ~ NaSSRI, dataset_comb_med_all))
summary(lm(sd_delay ~ NaSSRI, dataset_comb_med_all))
summary(lm(NREM_percent ~ NaSSRI, dataset_comb_med_all))
summary(lm(CONS_change ~ NaSSRI, dataset_comb_med_all))

summary(lm(spin_dens ~ Rest, dataset_comb_med_all))
summary(lm(spin_count ~ Rest, dataset_comb_med_all))
summary(lm(so_amp ~ Rest, dataset_comb_med_all))
summary(lm(so_dur ~ Rest, dataset_comb_med_all))
summary(lm(coupling_occur ~ Rest, dataset_comb_med_all))
summary(lm(mean_delay ~ Rest, dataset_comb_med_all))
summary(lm(sd_delay ~ Rest, dataset_comb_med_all))
summary(lm(NREM_percent ~ Rest, dataset_comb_med_all))
summary(lm(CONS_change ~ Rest, dataset_comb_med_all))

## seperate drugs
summary(lm(spin_dens ~ venlafaxin, dataset_comb_med_all))#
summary(lm(spin_count ~ venlafaxin, dataset_comb_med_all))
summary(lm(so_amp ~ venlafaxin, dataset_comb_med_all))
summary(lm(so_dur ~ venlafaxin, dataset_comb_med_all))
summary(lm(coupling_occur ~ venlafaxin, dataset_comb_med_all))
summary(lm(mean_delay ~ venlafaxin, dataset_comb_med_all))
summary(lm(sd_delay ~ venlafaxin, dataset_comb_med_all))
summary(lm(NREM_percent ~ venlafaxin, dataset_comb_med_all))
summary(lm(REM_percent ~ venlafaxin, dataset_comb_med_all))#
summary(lm(CONS_change ~ venlafaxin, dataset_comb_med_all))

summary(lm(spin_dens ~ mirtazapin, dataset_comb_med_all))
summary(lm(spin_count ~ mirtazapin, dataset_comb_med_all))
summary(lm(so_amp ~ mirtazapin, dataset_comb_med_all))
summary(lm(so_dur ~ mirtazapin, dataset_comb_med_all))
summary(lm(coupling_occur ~ mirtazapin, dataset_comb_med_all))
summary(lm(mean_delay ~ mirtazapin, dataset_comb_med_all))
summary(lm(sd_delay ~ mirtazapin, dataset_comb_med_all))
summary(lm(NREM_percent ~ mirtazapin, dataset_comb_med_all))
summary(lm(REM_percent ~ mirtazapin, dataset_comb_med_all))
summary(lm(CONS_change ~ mirtazapin, dataset_comb_med_all))

summary(lm(spin_dens ~ trimipramin, dataset_comb_med_all))
summary(lm(spin_count ~ trimipramin, dataset_comb_med_all))
summary(lm(so_amp ~ trimipramin, dataset_comb_med_all))
summary(lm(so_dur ~ trimipramin, dataset_comb_med_all))
summary(lm(coupling_occur ~ trimipramin, dataset_comb_med_all))
summary(lm(mean_delay ~ trimipramin, dataset_comb_med_all))
summary(lm(sd_delay ~ trimipramin, dataset_comb_med_all))
summary(lm(NREM_percent ~ trimipramin, dataset_comb_med_all))
summary(lm(REM_percent ~ trimipramin, dataset_comb_med_all))
summary(lm(CONS_change ~ trimipramin, dataset_comb_med_all))

cor.test(dataset_comb_med_all$sleep, dataset_comb_med_all$so_amp)
cor.test(dataset_comb_med_all$TCA, dataset_comb_med_all$so_amp)
cor.test(dataset_comb_med_all$SNRI, dataset_comb_med_all$REM_percent)
cor.test(dataset_comb_med_all$SNRI, dataset_comb_med_all$spin_dens)
cor.test(dataset_comb_med_all$SSRI, dataset_comb_med_all$REM_percent)
cor.test(dataset_comb_med_all$SSRI, dataset_comb_med_all$spin_dens)
cor.test(dataset_comb_med_all$venlafaxin, dataset_comb_med_all$spin_dens)
cor.test(dataset_comb_med_all$venlafaxin, dataset_comb_med_all$REM_percent)
  

## Specific for Dataset B
dataset_b_unmed  <- dataset_b_conpat %>%   filter(group  == "patient_unmed")
summary(lm(hamd_0 ~  spin_amp, dataset_b_unmed))
summary(lm(spin_amp  ~ no_episodes , dataset_b_unmed))

cor.test(dataset_b_unmed$spin_amp, dataset_b_unmed$no_episodes)
dataset_b_unmed_na <- dataset_b_unmed
dataset_b_unmed_na$no_episodes[34] <- NA #outlier
cor.test(dataset_b_unmed_na$spin_amp, dataset_b_unmed_na$no_episodes)


## Finger tapping -----
## Source again to get clean datasets
source("dataset_a.R")
source("dataset_b.R")

dataset_a$dataset <- "A"
dataset_b$dataset <- "B"

dataset_a <- dataset_a %>% dplyr::rename(BASE = START)

colnames_a <- colnames(dataset_a)
colnames_b <- colnames(dataset_b)

common_colnames <- intersect(colnames_a,  colnames_b)

dataset_a <- dataset_a %>% dplyr::select(common_colnames)
dataset_b <- dataset_b %>% dplyr::select(common_colnames)

dataset_a_comb <- dataset_a
dataset_b_comb <- dataset_b

dataset_a_comb$group<- revalue(dataset_a_comb$group, c("Patients" = "Medicated"))
dataset_b_comb$group<- revalue(dataset_b_comb$group, c("control" = "Controls",
                                                       "patient_unmed" = "Unmedicated",
                                                       "patient_7d" = "Medicated"))
dataset_comb_beh<- rbind(dataset_a_comb, dataset_b_comb)

dataset_comb_beh_med <- dataset_comb_beh %>% filter(group!= "Unmedicated")

## Get summary statistics 
dataset_comb_beh_med_long <- dataset_comb_beh_med %>% gather(memory, score, BASE,TRAIN_change, CONS_change)
dataset_comb_beh_med_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")

## T- test Group differences in behaviour
beh_BASE_ttest_dataset_comb_beh <- t.test(BASE ~ group, data = dataset_comb_beh_med)
beh_train_ttest_dataset_comb_beh<- t.test(TRAIN_change ~ group, data = dataset_comb_beh_med)
beh_cons_ttest_dataset_comb_beh <- t.test(CONS_change  ~ group, data = dataset_comb_beh_med)

# Plots 
ggplot(dataset_comb_beh_med,aes(x=group,y=BASE,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=BASE),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Correctly tapped')+ xlab('')+ my_theme +   raincloud_theme + 
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Baseline")+ scale_fill_manual(values = color_dataset_comb) +  scale_colour_manual(values = color_dataset_comb)  +
  geom_signif(annotation=formatC("*", digits=1), y_position=19.5,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK") 

ggplot(dataset_comb_beh_med,aes(x=group,y=TRAIN_change,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=TRAIN_change),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Change performance')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Training")+ scale_fill_manual(values = color_dataset_comb) +  scale_colour_manual(values = color_dataset_comb) 

ggplot(dataset_comb_beh_med,aes(x=group,y=CONS_change,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=CONS_change),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Change performance')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Consolidation")+ scale_fill_manual(values = color_dataset_comb) +  scale_colour_manual(values = color_dataset_comb) +
  geom_signif(annotation=formatC("***", digits=1), y_position=0.9,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK") 

# remove outliers 3sd from the mean ---- 
dataset_comb_beh_med_rm <- dataset_comb_beh_med
which_outl(dataset_comb_beh_med, "BASE")        #1
which_outl(dataset_comb_beh_med, "TRAIN_change")  #3
which_outl(dataset_comb_beh_med, "CONS_change") #1
which_outl(dataset_comb_beh_med, "CONS_change") #1 <- removes 2 outliers
dataset_comb_beh_med_rm <- rm_outl(dataset_comb_beh_med_rm, "BASE")        #1
dataset_comb_beh_med_rm <- rm_outl(dataset_comb_beh_med_rm, "TRAIN_change")  #3
# dataset_comb_beh_med_rm$CONS_change[139] <-NA
# dataset_comb_beh_med_rm$CONS_change[146] <-NA

beh_BASE_ttest_dataset_comb_beh_rm  <- t.test(BASE ~ group, data = dataset_comb_beh_med_rm)
beh_train_ttest_dataset_comb_beh_rm <- t.test(TRAIN_change ~ group, data = dataset_comb_beh_med_rm)
beh_cons_ttest_dataset_comb_beh_rm  <- t.test(CONS_change  ~ group, data = dataset_comb_beh_med_rm)

dataset_comb_beh_med_long_rm <- dataset_comb_beh_med_rm %>% gather(memory, score, BASE,TRAIN_change, CONS_change)
dataset_comb_beh_med_long_rm %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")

dataset_comb_beh_med_long <- dataset_comb_beh_med %>% gather(memory, score, BASE,TRAIN_change, CONS_change)
dataset_comb_beh_med_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")


# Plots 
rainplot_baseline_dataset_comb_beh <- ggplot(dataset_comb_beh_med_rm,aes(x=group,y=BASE,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=BASE),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Baseline correctly tapped')+ xlab('')+ my_theme +   raincloud_theme + 
  guides(fill=FALSE,colour=FALSE)+ theme(axis.title.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.line.x = element_blank(),
                                         axis.text.x = element_blank()) + ylim(0,22) +
  scale_fill_manual(values = color_dataset_comb) +  scale_colour_manual(values = color_dataset_comb)  +
  geom_signif(annotation=formatC("*", digits=1), y_position=20,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK") 

rainplot_training_dataset_comb_beh <- ggplot(dataset_comb_beh_med_rm,aes(x=group,y=TRAIN_change*100,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=TRAIN_change*100),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Training \npercent change [%]')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+ theme(axis.title.x = element_blank(),
                                         axis.ticks.x = element_blank(),
                                         axis.line.x = element_blank(),
                                         axis.text.x = element_blank()) +
  scale_fill_manual(values = color_dataset_comb) +  scale_colour_manual(values = color_dataset_comb) 

rainplot_consolidation_dataset_comb_beh  <- ggplot(dataset_comb_beh_med_rm,aes(x=group,y=CONS_change*100,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=CONS_change*100),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Consolidation \npercent change [%]')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+   theme(axis.title.x = element_blank(),
                                           axis.ticks.x = element_blank(),
                                           axis.line.x = element_blank(),
                                           axis.text.x = element_blank()) + ylim(-100, 110)+
  scale_fill_manual(values = color_dataset_comb) +  scale_colour_manual(values = color_dataset_comb) +
  geom_signif(annotation=formatC("***", digits=1), y_position=105,xmin=1.2, xmax=2.2, 
              textsize = 8, colour="BLACK") 

beh_plot_dataset_comb_beh     <- ggarrange(rainplot_baseline_dataset_comb_beh, rainplot_training_dataset_comb_beh, 
                                           rainplot_consolidation_dataset_comb_beh, nrow = 3, ncol = 1)
# ggsave("figs/beh_plot_dataset_comb_beh.png", beh_plot_dataset_comb_beh, device = "png", dpi = 300, width = 8, height = 12)
# ggsave("figs/beh_plot_dataset_comb_beh.pdf", beh_plot_dataset_comb_beh, device = "pdf", dpi = 300, width = 8, height = 10)


###### Behaviour and sleep parameters
summary(lm(CONS_change ~ spin_dens*group, dataset_comb_beh_med_rm))

dataset_comb_beh_med_rm_singles<- dataset_comb_beh_med
dataset_comb_beh_med_rm_singles$CONS_change[139]<- NA
dataset_comb_beh_med_rm_singles$CONS_change[146]<- NA

lapply(dataset_comb_beh_med_rm[,c("spin_dens","spin_count","spin_amp","spin_dur", 
                                  "spin_freq", "so_dens","so_dur","so_amp" ,"so_count","so_freq", 
                                  "coupling_occur","mean_delay","sd_delay"
                                  )], function(x) summary(lm(dataset_comb_beh_med_rm$CONS_change ~ x * dataset_comb_beh_med_rm$group)))

dataset_comb_beh_med_rm_singles_pat <- dataset_comb_beh_med_rm_singles %>% filter(group== "Medicated")
dataset_comb_beh_med_rm_singles_con <- dataset_comb_beh_med_rm_singles %>% filter(group== "Controls")

cor.test(dataset_comb_beh_med_rm_singles_pat$CONS_change, dataset_comb_beh_med_rm_singles_pat$spin_dens)
cor.test(dataset_comb_beh_med_rm_singles_con$CONS_change, dataset_comb_beh_med_rm_singles_con$spin_dens)

