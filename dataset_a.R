## Script to analyse Dataset A
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

## Load packages ----
options(warn = -1)
library(xlsx)
library(readxl)
library(tidyverse)
library(ggpubr)
library(plotrix)
library(ggsignif)
library(rstatix)
library(foreign)
library(BayesFactor)
library(rmatio)
options(warn = 1)

## Source functions ----
source("funcs/theme.R")
source("funcs/rainclouds.R")
source("funcs/outliers.R")

## Main script  ----
## Part I: Read in and prepare data ----
## Read in behavioural data
dataset_a <- read_excel("dataset_a_subject_info.xlsx") 
dataset_a$group <- as.factor(dataset_a$group)

## Order by datasetname
dataset_a <- dataset_a[order(dataset_a$datasetnum),] 

## Read in sleep table
hypval_dataset_a  <- read.delim("dataset_a_hypvals.csv", sep=",",dec= ".") 

## Select variables
coupling_dens_info <- hypval_dataset_a %>%
  dplyr::select(datasetnum, S2_without_MA_min,
                S3_without_MA_min,S4_without_MA_min)

hypval_dataset_a <- hypval_dataset_a %>% 
  dplyr::select(datasetnum, Total_sleep_time_min:SWS_onset_min, 
                S2_onset_min:Wake_after_sleep_onset_min,SWS_min, 
                S1_percent:Wake_after_sleep_onset_percent) 

## Combine with subject data
dataset_a <- merge(dataset_a, hypval_dataset_a, by = "datasetnum")

# Combine S3 and S4 to fit AASM scoring
dataset_a$S3andS4_min   <- dataset_a$S3_min     + dataset_a$S4_min
dataset_a$S3andS4_perc  <- dataset_a$S3_percent + dataset_a$S4_percent
dataset_a$NREM_percent  <- dataset_a$S2_percent + dataset_a$S3andS4_perc

## Read in spindle data
spindle_dataset_a <- read.delim("dataset_a_spindles.csv", sep=",",dec= ".")   

## Recode to match all different datasets
spindle_dataset_a <- mutate(spindle_dataset_a, channel = fct_recode(channel,"EEGC3" = "C3A2 EEG","EEGC4" = "C4A1 EEG"))

## Split by channel to merge average below
spindle_C3 <- spindle_dataset_a[spindle_dataset_a$channel == "EEGC3", ]
spindle_C4 <- spindle_dataset_a[spindle_dataset_a$channel == "EEGC4", ]

## Create spindle parameters by averaging over C3 and C4 
spindle_C3$spin_dens       <- (spindle_C3$density_per_epoch + spindle_C4$density_per_epoch)/2 
spindle_C3$spin_centerfreq <- (spindle_C3$used_center_freq+spindle_C4$used_center_freq)/2
spindle_C3$spin_count      <- (spindle_C3$count+spindle_C4$count)/2  
spindle_C3$spin_amp        <- (spindle_C3$mean_amplitude_trough2peak_potential+spindle_C4$mean_amplitude_trough2peak_potential)/2  
spindle_C3$spin_dur        <- (spindle_C3$mean_duration_seconds+spindle_C4$mean_duration_seconds)/2  
spindle_C3$spin_freq       <- (spindle_C3$mean_frequency_by_mean_pk_trgh_cnt_per_dur+spindle_C4$mean_frequency_by_mean_pk_trgh_cnt_per_dur)/2  
spindle_C3$spin_act        <-  spindle_C3$spin_amp * spindle_C3$spin_dur

spindle_params <- spindle_C3 %>% dplyr::select(datasetnum, spin_dens,spin_centerfreq, spin_count, spin_amp, spin_dur, spin_freq, spin_act)
dataset_a <- merge(dataset_a, spindle_params, by = "datasetnum")

## Read in SO data 
so_dataset_a  <- read.delim("dataset_a_so.csv", sep=",",dec= ".")   

## Recode to match all different datasets
so_dataset_a <- mutate(so_dataset_a, channel = fct_recode(channel,"EEGC3" = "C3A2 EEG","EEGC4" = "C4A1 EEG"))

## Split by channel to merge average below
so_C3 <- so_dataset_a[so_dataset_a$channel == "EEGC3", ]
so_C4 <- so_dataset_a[so_dataset_a$channel == "EEGC4", ]

## Create spindle parameters by averaging over C3 and C4
so_C3$so_count       <- (so_C3$count + so_C4$count)/2
so_C3$so_dens        <- (so_C3$density_per_epoch + so_C4$density_per_epoch)/2
so_C3$so_dur         <- (so_C3$mean_duration_seconds + so_C4$mean_duration_seconds)/2
so_C3$so_amp         <- (so_C3$mean_amplitude_peak2trough_potential+so_C4$mean_amplitude_peak2trough_potential)/2
so_C3$so_freq        <- (so_C3$mean_frequency_by_duration+so_C4$mean_frequency_by_duration)/2
so_C3$so_act         <-  so_C3$so_amp * so_C3$so_dur
so_C3$so_down_slope  <- (so_C3$mean_slope_to_trough_min_potential_per_second+so_C4$mean_slope_to_trough_min_potential_per_second)/2
so_C3$so_up_slope    <- (so_C3$mean_slope_trough_to_up_max_potential_per_second+so_C4$mean_slope_trough_to_up_max_potential_per_second)/2

so_params <- so_C3 %>% dplyr::select(datasetnum, so_count,so_dens, so_dur, so_amp, so_freq, so_act, so_down_slope, so_up_slope)
dataset_a <- merge(dataset_a, so_params, by = "datasetnum")
 
## Read in spindle-SO coupling data
coupling_dataset_a  <- read.delim("dataset_a_coupling.csv", sep=",",dec= ".") 

## Recode to match all different datasets
coupling_dataset_a <- mutate(coupling_dataset_a, test_channel = fct_recode(test_channel,"EEGC3" = "C3A2 EEG","EEGC4" = "C4A1 EEG"))
coupling_dataset_a <- mutate(coupling_dataset_a, target_channel = fct_recode(target_channel,"EEGC3" = "C3A2 EEG","EEGC4" = "C4A1 EEG"))

## Split by channel to merge average below
coupling_C3 <- coupling_dataset_a[coupling_dataset_a$test_channel == "EEGC3", ]
coupling_C4 <- coupling_dataset_a[coupling_dataset_a$test_channel == "EEGC4", ]

## Count the amount of co-occurences
coupling_occur_c3 <- table(coupling_C3$event_files_number)
coupling_occur_c4 <- table(coupling_C4$event_files_number)

## Get sum of the amount of co-occurences
coupling_occur     <- (coupling_occur_c3+coupling_occur_c4)/2  ##Calculate the mean of the two channels
coupling_occur_sum <- coupling_occur_c3+coupling_occur_c4  

coupling_occur     <- as.numeric(coupling_occur)
coupling_occur_sum <- as.numeric(coupling_occur_sum)

## Add to dataframe
dataset_a$coupling_occur     <- coupling_occur
dataset_a$coupling_occur_sum <- coupling_occur_sum

## Create the delay variable by subtracting the trough of the SO from the spindle
coupling_C3$delay_c3 <- coupling_C3$test_seconds_trough_max - coupling_C3$target_seconds_trough_max
coupling_C4$delay_c4 <- coupling_C4$test_seconds_trough_max - coupling_C4$target_seconds_trough_max

## Calculate mean and sd of delay per channel
mean_delay_c3 <- aggregate(delay_c3 ~ event_files_number, data=coupling_C3, FUN=mean)
sd_delay_c3   <- aggregate(delay_c3 ~ event_files_number, data=coupling_C3, FUN=sd)
mean_delay_c4 <- aggregate(delay_c4 ~ event_files_number, data=coupling_C4, FUN=mean)
sd_delay_c4   <- aggregate(delay_c4 ~ event_files_number, data=coupling_C4, FUN=sd)

## Average over channels
mean_delay     <- (mean_delay_c3+mean_delay_c4)/2
sd_delay       <- (sd_delay_c3+sd_delay_c4)/2  

# Add to main dataframe
dataset_a$mean_delay     <- mean_delay$delay_c3
dataset_a$sd_delay       <- sd_delay$delay_c3

## Supplemental analysis: get parameter information of matching versus mismatching events
## Based on MATLAB script (spisop_coupling_match_vs_mismatch.m) output; a .mat file containing matching and mismatching events
## Loop will fill create dataframe containing this info per participant
spin_dur_match    <- NULL
spin_dur_mismatch <- NULL
spin_amp_match    <- NULL
spin_amp_mismatch <- NULL
spin_freq_match   <- NULL
spin_freq_mismatch<- NULL
spin_perc_match   <- NULL
spin_perc_mismatch<- NULL
so_dur_match      <- NULL
so_dur_mismatch   <- NULL
so_amp_match      <- NULL
so_amp_mismatch   <- NULL
so_perc_match     <- NULL
so_perc_mismatch  <- NULL

i = 1
for (subnr in 1:80){ 
coupling_subj <- read.mat(paste(".../matching_mat/coupling_v1_sub",subnr, ".mat",sep=""))

spin_dur_match [i]   <- unlist(coupling_subj$coupling$match$spindle$duration)
spin_dur_mismatch[i] <- unlist(coupling_subj$coupling$mismatch$spindle$duration)
spin_amp_match[i]    <- unlist(coupling_subj$coupling$match$spindle$amplitude)
spin_amp_mismatch[i] <- unlist(coupling_subj$coupling$mismatch$spindle$amplitude)
spin_freq_match[i]   <- unlist(coupling_subj$coupling$match$spindle$frequency)
spin_freq_mismatch[i]<- unlist(coupling_subj$coupling$mismatch$spindle$frequency)
spin_perc_match[i]   <- unlist(coupling_subj$coupling$match$spindle$percentage)

so_dur_match[i]   <- unlist(coupling_subj$coupling$match$so$duration)
so_dur_mismatch[i]<- unlist(coupling_subj$coupling$mismatch$so$duration)
so_amp_match[i]   <- unlist(coupling_subj$coupling$match$so$amplitude)
so_amp_mismatch[i]<- unlist(coupling_subj$coupling$mismatch$so$amplitude)
so_perc_match[i]  <- unlist(coupling_subj$coupling$match$so$percentage)
  i = i +1     
} 

# Add to main dataframe
dataset_a$spin_dur_match    <- spin_dur_match
dataset_a$spin_dur_mismatch <- spin_dur_mismatch
dataset_a$spin_amp_match    <- spin_amp_match
dataset_a$spin_amp_mismatch <- spin_amp_mismatch
dataset_a$spin_freq_match   <- spin_freq_match
dataset_a$spin_freq_mismatch<- spin_freq_mismatch
dataset_a$spin_perc_match   <- spin_perc_match
dataset_a$so_dur_match      <- so_dur_match
dataset_a$so_dur_mismatch   <- so_dur_mismatch
dataset_a$so_amp_match      <- so_amp_match
dataset_a$so_amp_mismatch   <- so_amp_mismatch
dataset_a$so_perc_match     <- so_perc_match

dataset_a$spin_dur_diff  <- spin_dur_match - spin_dur_mismatch
dataset_a$spin_amp_diff  <- spin_amp_match - spin_amp_mismatch
dataset_a$spin_freq_diff <- spin_freq_match - spin_freq_mismatch
dataset_a$so_dur_diff    <- so_dur_match - so_dur_mismatch
dataset_a$so_amp_diff    <- so_amp_match - so_amp_mismatch

## Read in power data 
power_dataset_a  <- read.delim("dataset_a_pow.csv", sep=",",dec= ".") 

## change to numeric
power_dataset_a$freq <- as.numeric(as.character(power_dataset_a$freq))
power_dataset_a$mean_powerDensity_over_segments <- as.numeric(as.character(power_dataset_a$mean_powerDensity_over_segments))

power_dataset_a <- mutate(power_dataset_a, channel = fct_recode(channel,"EEGC3" = "C3A2 EEG","EEGC4" = "C4A1 EEG"))

# Merge with data
mergedData_long_dataset_a <- merge(dataset_a, power_dataset_a, by=c("datasetnum"))

## Reshape the data, per frequency
power_dataset_a_freq <- power_dataset_a %>%
  dplyr::select(datasetnum, channel, freq, mean_powerDensity_over_segments) %>%
  spread(key = freq, value = mean_powerDensity_over_segments)

## Average over the two channels
power_dataset_a_freq_reshape <- aggregate(. ~ datasetnum, data = power_dataset_a_freq, mean)

## Merge with dataframe in wide shape
mergedData_wide_dataset_a <- merge(dataset_a, power_dataset_a_freq_reshape, by = ("datasetnum"))

# Transform to decibel
mergedData_long_dataset_a$mean_densDB <- mergedData_long_dataset_a$mean_powerDensity_over_segments + 1
mergedData_long_dataset_a$mean_densDB <- 10*log10(mergedData_long_dataset_a$mean_densDB)

## Calculate t.test for each of the frequencies between groups CONTROLS AND PATIENTS
column_dataset_a_06  <- which(colnames(mergedData_wide_dataset_a)=="0.6")
listtests   <- lapply(mergedData_wide_dataset_a[,c(column_dataset_a_06:(column_dataset_a_06+147))], function(x) t.test(x ~ mergedData_wide_dataset_a$group, var.equal = FALSE))

## Extract the p-value and make it a dataframe format
listpvalues <- sapply(listtests, function(x) x[3])
listpvalues <- as.data.frame(listpvalues)

## Reshape into long format
listpvalues <- listpvalues %>% gather(key = freq , value = pvalue, starts_with("X")  )  

new_names        <- seq(0.6, 30, by=0.2)
listpvalues$Freq <- paste0(new_names)
listpvalues$Freq <- as.numeric(listpvalues$Freq)
listpvalues      <- within(listpvalues, rm(freq))

data_long_pvalues <- merge(mergedData_long_dataset_a, listpvalues, by.x=c("freq"), by.y=c("Freq"))

###### Permutation test -- patients controls
freq_columns_start   <- which(colnames(mergedData_wide_dataset_a)=="0.6")
freq_columns         <-  freq_columns_start:(freq_columns_start+147)

pb = txtProgressBar(min = 0, max = length(freq_columns), initial = 0) 
perm_pvalue = rep(NA,length(freq_columns))
index = 1

for(i in freq_columns_start:(freq_columns_start+147)){
  ## Make new dataset
  perm_data <- mergedData_wide_dataset_a %>% dplyr::select(i, group) 
  
  # Split the data
  controls_freq <- perm_data[perm_data$group == "Controls",]
  patients_freq <- perm_data[perm_data$group == "Patients",]
  
  # Make seperate lists
  combined_groups <- c(controls_freq$group, patients_freq$group)
  combined_freqs  <- c(controls_freq[,1], patients_freq[,1])
  
  # Number of simulations 
  nsims <- 10000
  
  diff_obs <- mean(controls_freq[,1])- mean(patients_freq[,1])
  diffs <- rep(NA, nsims)
  
  ## Permutation for loop
  for(k in 1:nsims){
    shuffled_labels <- sample(combined_groups, replace = FALSE)
    diffs[k] <- mean(combined_freqs[shuffled_labels == 1]) - mean(combined_freqs[shuffled_labels == 2])
    
  }
  
  ## p-value = (number of more extreme differences than diff_obs/nsims)
  length(diffs[abs(diffs)>= abs(diff_obs)])/nsims
  perm_pvalue[index] <- mean(abs(diffs) > abs(diff_obs))
  index = index + 1
  
  setTxtProgressBar(pb,i)
  print
}  

## Add them to the other dataframe
perm_pvalue <- as.data.frame(perm_pvalue)

## Catch loop in case a value is zero (will result in errors)
for (i in 1:nrow(perm_pvalue)) {
  if ( perm_pvalue[i,1] == 0.000){
    perm_pvalue[i,1] <- 0.0001
  }
}

new_names <- seq(0.6, 30, by=0.2)
perm_pvalue$Freq <- paste0(new_names)
perm_pvalue$Freq<- as.numeric(perm_pvalue$Freq)

dataset_a_long_pvalues_perm <- merge(data_long_pvalues, perm_pvalue, by.x=c("freq"), by.y=c("Freq"))

############################# Split datasets #############################

dataset_a <-  mutate(dataset_a, x = fct_recode(group,  "control" = "Controls", "Patients" = "patient"))

## Split into seperate dataframes:
## controls 
dataset_a_con <- dataset_a %>% 
  dplyr::filter(group == "Controls")
dataset_a_con$group <- droplevels(dataset_a_con$group)
dataset_a_con$group <- factor(dataset_a_con$group)

## patients 
dataset_a_pat <- dataset_a %>% 
  dplyr::filter(group == "Patients")
dataset_a_pat$group <- droplevels(dataset_a_pat$group)
dataset_a_pat$group <- factor(dataset_a_pat$group)

############################# Part II: Summary statistics and ttests #############################

## Demographics ----
dataset_a_demo_long <- dataset_a %>%
  gather(demo, score, age, gender, hamd_0, no_episodes)

dataset_a_demo_long %>%
  dplyr::select(demo, score, group) %>%
  group_by(group, demo) %>%
  get_summary_stats(score, type = "mean_se")

## Sleep architecture ---- 
dataset_a_hypno_long <- dataset_a %>%
  gather(sleepstage, score, S1_percent,S2_percent,S3andS4_perc,REM_percent,NREM_percent,Wake_after_sleep_onset_percent,Total_sleep_time_min, Sleep_Onset_min, SWS_onset_min, REM_onset_min)

dataset_a_hypno_min_long <- dataset_a %>%
  gather(sleepstage, score, S1_min,S2_min,S3_min,REM_min)

dataset_a_hypno_min_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

dataset_a_hypno_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

hypnogram_ttests_perc_dataset_a    <- lapply(dataset_a[,c("age", "S1_percent", "S2_percent", "S3andS4_perc",  "REM_percent", "NREM_percent", "Wake_after_sleep_onset_percent", "Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min" )], function(x) t.test(x ~ dataset_a$group, var.equal = FALSE))
hypnogram_S4_ttests_perc_dataset_a <- t.test(dataset_a$S4_percent ~ dataset_a$group) 

## Calculate bayesfactors
ttestBF(formula = S1_percent ~ group, data = dataset_a)
ttestBF(formula = S2_percent ~ group, data = dataset_a)
ttestBF(formula = S3andS4_perc ~ group, data = dataset_a)
ttestBF(formula = NREM_percent ~ group, data = dataset_a)
ttestBF(formula = REM_percent ~ group, data = dataset_a)
ttestBF(formula = Wake_after_sleep_onset_percent ~ group, data = dataset_a)
ttestBF(formula = Total_sleep_time_min ~ group, data = dataset_a)
ttestBF(formula = Sleep_Onset_min ~ group, data = dataset_a)
ttestBF(formula = SWS_onset_min ~ group, data = dataset_a)
ttestBF(formula = REM_onset_min ~ group, data = dataset_a)

## Correlation between REM and S1 in patients
cor.test(dataset_a_pat$S1_percent, dataset_a_pat$REM_percent)

## Spindle data ---- 
dataset_a_spindle_long <- dataset_a %>%
  gather(spindle, score,  spin_dens, spin_count,spin_amp, spin_dur, spin_freq)

dataset_a_spindle_long %>%
  dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%
  get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_a <- lapply(dataset_a[,c("spin_dens","spin_dur","spin_amp", "spin_count", "spin_centerfreq", "spin_freq" )], function(x) t.test(x ~ dataset_a$group, var.equal = FALSE))
 
## without outlier - one control value below the mean
dataset_a_spin_dens_rm <- rm_outl(dataset_a, "spin_dens")

dataset_a_spindle_rm_long <- dataset_a_spin_dens_rm %>%
  gather(spindle, score,  spin_dens, spin_count,spin_amp, spin_dur)

dataset_a_spindle_rm_long %>% dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_a_rm <- lapply(dataset_a_spin_dens_rm[,c("spin_dens","spin_dur","spin_amp", "spin_count" )], function(x) t.test(x ~ dataset_a_spin_dens_rm$group, var.equal = FALSE))

## Calculate bayes factor
ttestBF(formula = spin_dens ~ group, data = dataset_a)
ttestBF(formula = spin_count ~ group, data = dataset_a)
ttestBF(formula = spin_amp ~ group, data = dataset_a)
ttestBF(formula = spin_freq ~ group, data = dataset_a)
ttestBF(formula = spin_dur ~ group, data = dataset_a)

## SO data ----  
dataset_a_so_long <- dataset_a %>%
  gather(so, score,  so_dens,so_amp, so_dur, so_count, so_freq)

dataset_a_so_long %>%
  dplyr::select(so, score, group) %>%
  group_by(group, so) %>%
  get_summary_stats(score, type = "mean_se")

so_ttests_dataset_a <- lapply(dataset_a[,c("so_dens","so_dur","so_amp" ,"so_count","so_freq")], function(x) t.test(x ~ dataset_a$group, var.equal = FALSE))

ttestBF(formula = so_dens ~ group, data = dataset_a)
ttestBF(formula = so_count ~ group, data = dataset_a)
ttestBF(formula = so_amp ~ group, data = dataset_a)
ttestBF(formula = so_freq ~ group, data = dataset_a)
ttestBF(formula = so_dur ~ group, data = dataset_a)

## Coupling data ---- 
dataset_a_coupling_dens_info <- merge(dataset_a, coupling_dens_info, by = "datasetnum")
dataset_a_coupling_dens_info$NREM_min <- (dataset_a_coupling_dens_info$S2_without_MA_min +
                                              dataset_a_coupling_dens_info$S3_without_MA_min) /2 
dataset_a_coupling_dens_info$NREM_min <- (dataset_a_coupling_dens_info$NREM_min +
                                            dataset_a_coupling_dens_info$S4_without_MA_min) /2 

dataset_a$coupling_dens<- dataset_a_coupling_dens_info$coupling_occur  /  dataset_a_coupling_dens_info$NREM_min

dataset_a_coupling_long <- dataset_a %>%
  gather(coupling, score,  coupling_occur,mean_delay,sd_delay, coupling_dens)

dataset_a_coupling_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_ttests_dataset_a <- lapply(dataset_a[,c("coupling_occur","mean_delay","sd_delay", "coupling_dens" )], function(x) t.test(x ~ dataset_a$group, var.equal = FALSE))

ttestBF(formula = coupling_occur ~ group, data = dataset_a)
ttestBF(formula = mean_delay ~ group, data = dataset_a)
ttestBF(formula = sd_delay ~ group, data = dataset_a)

## Coupling matching vs nonmatching spindles
dataset_a_coupling_matchspin_long <- dataset_a %>%
  gather(coupling, score, spin_amp_match, spin_dur_match,spin_freq_match,spin_perc_match,
         so_amp_match, so_dur_match,so_perc_match)

dataset_a_coupling_matchspin_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

dataset_a_coupling_diffspin_long <- dataset_a %>%
  gather(coupling, score, spin_amp_diff, spin_dur_diff,spin_freq_diff, so_amp_diff, so_dur_diff)

dataset_a_coupling_diffspin_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_diff_ttests_dataset_a <- lapply(dataset_a[,c("spin_amp_diff","spin_dur_diff","spin_freq_diff","spin_perc_match","so_amp_diff", "so_dur_diff", "so_perc_match" )], function(x) t.test(x ~ dataset_a$group, var.equal = FALSE))
coupling_match_ttests_dataset_a <- lapply(dataset_a[,c("spin_amp_match","spin_dur_match","spin_freq_match","so_amp_match", "so_dur_match")], function(x) t.test(x ~ dataset_a$group, var.equal = FALSE))

ttestBF(formula = so_dur_match ~ group, data = dataset_a)
ttestBF(formula = so_amp_match ~ group, data = dataset_a)
ttestBF(formula = spin_amp_match ~ group, data  = dataset_a)
ttestBF(formula = spin_dur_match ~ group, data  = dataset_a)
ttestBF(formula = spin_freq_match ~ group, data = dataset_a)
ttestBF(formula = so_amp_match ~ group, data    = dataset_a)
ttestBF(formula = so_dur_match ~ group, data    = dataset_a)
ttestBF(formula = spin_amp_diff ~ group, data  = dataset_a)
ttestBF(formula = spin_dur_diff ~ group, data  = dataset_a)
ttestBF(formula = spin_freq_diff ~ group, data = dataset_a)
ttestBF(formula = so_amp_diff ~ group, data    = dataset_a)
ttestBF(formula = so_dur_diff ~ group, data    = dataset_a)


## Part III: Plots ---- 
## Power plot
freq_spec_dataset_a <- dataset_a_long_pvalues_perm %>%
  ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
  theme_classic()+ my_theme +
  stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
  scale_color_manual(values=color_dataset_a, labels = c("Controls (A)", "Medicated (A)")) +  
  ggtitle("Frequency spectrum")+
  ylab("Power (dB)") +
  theme(plot.title  = element_text(hjust = 0.5)) +                                
  theme(      
    axis.text.y  = element_text( color="#000000"),
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.9, 0.4)) + 
  NULL

perm_plot_dataset_a <- ggplot(dataset_a_long_pvalues_perm, aes(x=freq, y=abs(log10(perm_pvalue)), colour="#000000")) +
  theme_classic()+ my_theme +
  stat_summary(fun.y="mean", geom="bar", colour="black")+ 
  xlab("\nFrequency") + ylab("P-values") +
  theme(       
    axis.text.y = element_text( color="#000000"),
    axis.text.x = element_text(color="#000000"), 
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 4),labels=c("0",  "10e-1","10e-2", "10e-3","10e-4"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) 
NULL

## Spindles
plot_rain_spinamp_dataset_a <- ggplot(dataset_a,aes(x=group,y=spin_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle amplitude')+ xlab('')+ ylim(0,100) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) + 
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in spindle amplitude")

dataset_a_spinamp_rm <- rm_outl(dataset_a, "spin_amp")

plot_rain_spinamp_dataset_a <- ggplot(dataset_a_spinamp_rm,aes(x=group,y=spin_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle amplitude')+ xlab('')+ ylim(0,55) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) + 
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in spindle amplitude")

plot_rain_spidens_dataset_a <- ggplot(dataset_a,aes(x=group,y=spin_dens,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_dens),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle density')+ xlab('')+ ylim(0,4) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) + 
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in spindle density") +
  geom_signif(annotation=formatC("#", digits=1), y_position=3.6,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

##  SOs
plot_rain_soamp_dataset_a <- ggplot(dataset_a,aes(x=group,y=so_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO amplitude')+ xlab('')+ ylim(0,300) +
  raincloud_theme + my_theme + 
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) + 
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in SO amplitude") +
  geom_signif(annotation=formatC("***", digits=1), y_position=280,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

plot_rain_sodur_dataset_a <- ggplot(dataset_a,aes(x=group,y=so_dur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_dur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO duration')+ xlab('')+ ylim(1,1.6) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) +
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in SO duration")+
  geom_signif(annotation=formatC("***", digits=1), y_position=1.5,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

plot_rain_sodens_dataset_a <- ggplot(dataset_a,aes(x=group,y=so_dens,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_dens),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO density')+ xlab('')+ ylim(0,3) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) + 
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in SO density")

## Coupling
plot_rain_amount_coup_dataset_a <- ggplot(dataset_a,aes(x=group,y=coupling_occur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=coupling_occur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO-fast spindles')+ xlab('')+ ylim(0,300) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) + 
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in amount of SO-fast spindle couplings")

plot_rain_mean_delay_dataset_a <- ggplot(dataset_a,aes(x=group,y=mean_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=mean_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Mean delay')+ xlab('')+ ylim(0.2,0.85) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) + 
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in mean delay") 

plot_rain_sd_delay_dataset_a <- ggplot(dataset_a,aes(x=group,y=sd_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=sd_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Delay dispersion (SD)')+ xlab('')+ ylim(0.15,0.45) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_a) + 
  scale_fill_manual(values=color_dataset_a) + 
  scale_x_discrete(labels = c("Controls (A)", "Medicated (A)")) +
  ggtitle("Group difference in delay dispersion (SD)")+
  geom_signif(annotation=formatC("**", digits=1), y_position=0.40,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

## Add the plots on top of eachother
freq_plot_dataset_a     <- ggarrange(freq_spec_dataset_a, perm_plot_dataset_a, heights = c(2, 1), ncol = 1, nrow = 2)
spin_plot_dataset_a     <- ggarrange(plot_rain_spidens_dataset_a, plot_rain_spinamp_dataset_a,nrow = 2, labels = c("A", "B") )
so_plot_dataset_a       <- ggarrange(plot_rain_soamp_dataset_a, plot_rain_sodur_dataset_a,nrow = 2, labels = c("A", "B") )
coupling_plot_dataset_a <- ggarrange(plot_rain_amount_coup_dataset_a, plot_rain_mean_delay_dataset_a, plot_rain_sd_delay_dataset_a, nrow = 3, labels = c("A", "B", "C") )

## Save the images
# ggsave("figs/freq_plot_dataset_a.png", freq_plot_dataset_a, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/spin_plot_dataset_a_inverse_edf_142.png", spin_plot_dataset_a, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/so_plot_dataset_a_inverse_edf_142.png", so_plot_dataset_a, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/coupling_plot_dataset_a_inverse_edf_142.png", coupling_plot_dataset_a, device = "png", dpi = 300, width = 7, height = 10)


## Part IV: Additional analyses ----
##  Finger tapping ----
## Read in data
tapping  <- read.delim("dataset_a_tapping.csv", sep=",",dec= ".") 
tapping  <- subset(tapping, select=c("START", "END", "TEST", "TRAIN", "CONS",  "name")) 

dataset_a <- merge(dataset_a, tapping, by.x = "data_name", by.y = "name" )

dataset_a$TRAIN_change <- ((dataset_a$END -  dataset_a$START)/ dataset_a$START)
dataset_a$CONS_change  <- ((dataset_a$TEST - dataset_a$END)  / dataset_a$END)

## Get summary statistics 
dataset_a_long <- dataset_a %>% gather(memory, score, START, END, TEST, TRAIN_change, CONS_change)
dataset_a_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")

## T- test group differences in behaviour
beh_start_ttest_dataset_a <- t.test(START ~ group, data = dataset_a)
beh_train_ttest_dataset_a <- t.test(TRAIN ~ group, data = dataset_a)
beh_cons_ttest_dataset_a  <- t.test(CONS  ~ group, data = dataset_a) # This is TEST/END --> percentage change 
beh_cons_ttest_dataset_a  <- t.test(CONS_change  ~ group, data = dataset_a)  

## remove outliers 3sd from the mean
dataset_a_rm <- dataset_a

dataset_a_rm <- rm_outl(dataset_a_rm, "START")        #1
dataset_a_rm <- rm_outl(dataset_a_rm, "TRAIN_change")  #1
dataset_a_rm <- rm_outl(dataset_a_rm, "CONS_change") #1

## Get summary statistics - after outlier removal
dataset_a_rm_long <- dataset_a_rm %>% gather(memory, score, START, END, TEST, TRAIN_change, CONS_change)
dataset_a_rm_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")

## T- test Group differences in behaviour - after outlier removal
beh_start_ttest_dataset_a_rm <- t.test(START ~ group, data = dataset_a_rm)
beh_train_ttest_dataset_a_rm <- t.test(TRAIN ~ group, data = dataset_a_rm)
beh_cons_ttest_dataset_a_rm <- t.test(CONS_change  ~ group, data = dataset_a_rm)  

## bayes statistics
bayes_start_ttest <- ttestBF(formula = START ~ group, data = dataset_a)
bayes_train_ttest <- ttestBF(formula = TRAIN ~ group, data = dataset_a)
bayes_cons_ttest  <- ttestBF(formula = CONS_change  ~ group, data = dataset_a)

## plots 
rainplot_baseline_dataset_a <- ggplot(dataset_a_rm,aes(x=group,y=START,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=START),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Correctly tapped')+ xlab('')+ my_theme +   raincloud_theme + 
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Baseline")+ scale_fill_manual(values = color_dataset_a) +  scale_colour_manual(values = color_dataset_a) +
  geom_signif(annotation=formatC("***", digits=1), y_position=15,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

rainplot_training_dataset_a  <- ggplot(dataset_a_rm,aes(x=group,y=TRAIN_change,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=TRAIN_change),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Change performance')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Training")+ scale_fill_manual(values = color_dataset_a) +  scale_colour_manual(values = color_dataset_a) 

rainplot_consolidation_dataset_a  <- ggplot(dataset_a_rm,aes(x=group,y=CONS_change,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=CONS_change),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Change performance')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Consolidation")+ scale_fill_manual(values = color_dataset_a) +  scale_colour_manual(values = color_dataset_a) +
  geom_signif(annotation=formatC("***", digits=1), y_position=0.5,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")

beh_plot_dataset_a     <- ggarrange(rainplot_baseline_dataset_a, rainplot_training_dataset_a, 
                                     rainplot_consolidation_dataset_a, ncol = 3, nrow = 1)

ggsave("figs/beh_plot_dataset_a.png", beh_plot_dataset_a, device = "png", dpi = 300, width = 14, height = 5)

######################################## Relation between sleep and memory consolidation ----

## sleep architecture and cons
lm_CONS_sleeptable_dataset_a <- lapply(dataset_a[,c("S1_percent","S2_percent","S3_percent","REM_percent","Wake_after_sleep_onset_percent","Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min")], function(x) summary(lm(dataset_a$CONS_change ~ x * dataset_a$group)))
lm_CONS_sleeptable_dataset_a_rm <- lapply(dataset_a_rm[,c("S1_percent","S2_percent","S3_percent","REM_percent","Wake_after_sleep_onset_percent","Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min")], function(x) lm(dataset_a_rm$CONS_change ~ x * dataset_a_rm$group))

## power and cons
select_vars <- c("1","1.2","1.4","1.6","1.8",
                 "2","2.2","2.4","2.6","2.8",
                 "3","3.2","3.4","3.6","3.8",
                 "4","4.2","4.4","4.6","4.8")

mergedData_wide_dataset_a <- mergedData_wide_dataset_a %>% mutate(av_delta = rowMeans(dplyr::select(., select_vars)))
mergedData_wide_dataset_a <- merge(mergedData_wide_dataset_a, tapping, by.x = "data_name", by.y = "name")
mergedData_wide_dataset_a$CONS_change  <- ((mergedData_wide_dataset_a$TEST - mergedData_wide_dataset_a$END)  / mergedData_wide_dataset_a$END)

mergedData_wide_dataset_a$av_delta <- mergedData_wide_dataset_a$av_delta + 1
mergedData_wide_dataset_a$av_delta <- 10*log10(mergedData_wide_dataset_a$av_delta)
 
## spindle params and cons
lm_CONS_spindle_dataset_a <- lapply(dataset_a[,c("spin_dens","spin_count","spin_amp","spin_dur", "spin_freq")], function(x) lm(dataset_a$CONS_change ~ x * dataset_a$group))
lm_CONS_spindle_dataset_a_rm <- lapply(dataset_a_rm[,c("spin_dens","spin_count","spin_amp","spin_dur", "spin_freq")], function(x) lm(dataset_a_rm$CONS_change ~ x * dataset_a_rm$group))

lm_cons_dens_dataset_a <- lm(CONS_change ~ spin_dens * group, dataset_a)
lm_cons_dens_dataset_a_rm <- lm(CONS_change ~ spin_dens * group, dataset_a_rm)

## coupling and cons
lm_CONS_coupling_dataset_a <- lapply(dataset_a[,c("coupling_occur","mean_delay","sd_delay")], function(x) lm(dataset_a$CONS_change ~ x * dataset_a$group))
lm_CONS_coupling_dataset_a_rm <- lapply(dataset_a_rm[,c("coupling_occur","mean_delay","sd_delay")], function(x) lm(dataset_a_rm$CONS_change ~ x * dataset_a_rm$group))

dataset_a_rm %>% 
  group_by(group) %>%
  filter(!is.na(CONS_change)) %>%
  dplyr::summarize(cor(CONS_change, coupling_occur))

dataset_a_rm %>% 
  group_by(group) %>%
  filter(!is.na(CONS_change)) %>%
  dplyr::summarize(cor(CONS_change, sd_delay))
 
## Depression severity ----
dataset_a_hamd_long <- dataset_a %>% filter(group=="Patients") %>% gather(hamd, score, hamd_0)
dataset_a_hamd_long %>% dplyr::select(hamd, score, group) %>% group_by(hamd) %>%  get_summary_stats(score, type = "mean_se")
