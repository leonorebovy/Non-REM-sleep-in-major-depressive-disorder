## Script to analyse Dataset C
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

## NOTE: one subject is missing (19_2 or E0000787)
## In addition, 5 controls and 2 7day patients were double with dataset A:
## - B0001465, C0001233, C0001238, C0001251, D0001406
## And C0000913, D0001049

## This means, to be deleted are:
# - 43, 44, 45, 59, 61
# - 19_2 and to match also 19_1
# - 17_1 and to match also 17_2
# - 27_1 and to match also 27_2

## Load packages ----
options(warn = -1)
library(xlsx)
library(readxl)
library(tidyverse)
library(ggpubr)
library(plotrix)
library(ggsignif)
library(rstatix)
library(coin)
library(BayesFactor)
library(rmatio)
options(warn = 1)

## Source functions ----
source("funcs/theme.R")
source("funcs/rainclouds.R")
source("funcs/outliers.R")
source("funcs/combine_pow_output_dataset_c.R")

## Main script  ----
## Part I: Read in and prepare data ----
## Read in behavioural data
dataset_c <- read_excel("dataset_c_subject_info.xlsx") 
dataset_c$group <- as.factor(dataset_c$group)
  
## Order by datasetname:
dataset_c <- dataset_c[order(dataset_c$datasetnum),] 

## Read in sleep table
hypval_dataset_c  <- read.delim("dataset_c_hypvals.csv", sep=",",dec= ".") 

coupling_dens_info <- hypval_dataset_c %>%
  dplyr::select(datasetnum, S2_without_MA_min,
                S3_without_MA_min,S4_without_MA_min)

hypval_dataset_c <- hypval_dataset_c %>% 
  dplyr::select(datasetnum, Total_sleep_time_min:SWS_onset_min, 
                S2_onset_min:Wake_after_sleep_onset_min,SWS_min, 
                S1_percent:Wake_after_sleep_onset_percent) 

dataset_c <- merge(dataset_c, hypval_dataset_c, by = "datasetnum", all.x = T)

# Combine S3 and S4 
dataset_c$S3andS4_min   <- dataset_c$S3_min     + dataset_c$S4_min
dataset_c$S3andS4_perc  <- dataset_c$S3_percent + dataset_c$S4_percent
dataset_c$NREM_percent  <- dataset_c$S2_percent + dataset_c$S3andS4_perc

## Read in spindle data 
spindle_dataset_c <- read.delim("dataset_c_spindles.csv", sep=",",dec= ".")   

## Recode to match all different datasets
spindle_dataset_c <- mutate(spindle_dataset_c, channel = fct_recode(channel, "EEGC3" = "EEG_-C3/A2", 
                                                                    "EEGC4" = "EEG_-C4/A1",
                                                                    "EEGC3" = "EEG C3",
                                                                    "EEGC4" = "EEG C4"))
# Split by channel to merge average below
spindle_C3 <- spindle_dataset_c[spindle_dataset_c$channel == "EEGC3", ]
spindle_C4 <- spindle_dataset_c[spindle_dataset_c$channel == "EEGC4", ]

# Create spindle parameters by averaging over C3 and C4 
spindle_C3$spin_dens       <- (spindle_C3$density_per_epoch + spindle_C4$density_per_epoch)/2 
spindle_C3$spin_centerfreq <- (spindle_C3$used_center_freq+spindle_C4$used_center_freq)/2
spindle_C3$spin_count      <- (spindle_C3$count+spindle_C4$count)/2  
spindle_C3$spin_amp        <- (spindle_C3$mean_amplitude_trough2peak_potential+spindle_C4$mean_amplitude_trough2peak_potential)/2  
spindle_C3$spin_dur        <- (spindle_C3$mean_duration_seconds+spindle_C4$mean_duration_seconds)/2  
spindle_C3$spin_freq       <- (spindle_C3$mean_frequency_by_mean_pk_trgh_cnt_per_dur+spindle_C4$mean_frequency_by_mean_pk_trgh_cnt_per_dur)/2  
spindle_C3$spin_act        <- spindle_C3$spin_amp * spindle_C3$spin_dur

spindle_params <- spindle_C3 %>% dplyr::select(datasetnum, spin_dens,spin_centerfreq, spin_count, spin_amp, spin_dur, spin_freq, spin_act)
dataset_c <- merge(dataset_c, spindle_params, by = "datasetnum")

## Read in SO data 
so_dataset_c  <- read.delim("dataset_c_so.csv", sep=",",dec= ".")
so_dataset_c <- mutate(so_dataset_c, channel = fct_recode(channel, "EEGC3" = "EEG_-C3/A2", 
                                                                    "EEGC4" = "EEG_-C4/A1",
                                                                    "EEGC3" = "EEG C3",
                                                                    "EEGC4" = "EEG C4"))
# Split by channel to merge average below
so_C3 <- so_dataset_c[so_dataset_c$channel == "EEGC3", ]
so_C4 <- so_dataset_c[so_dataset_c$channel == "EEGC4", ]

# Create spindle paramters by averageing over C3 and C4 
so_C3$so_count       <- (so_C3$count + so_C4$count)/2
so_C3$so_dens        <- (so_C3$density_per_epoch + so_C4$density_per_epoch)/2
so_C3$so_dur         <- (so_C3$mean_duration_seconds + so_C4$mean_duration_seconds)/2  
so_C3$so_amp         <- (so_C3$mean_amplitude_peak2trough_potential+so_C4$mean_amplitude_peak2trough_potential)/2
so_C3$so_freq        <- (so_C3$mean_frequency_by_duration+so_C4$mean_frequency_by_duration)/2  
so_C3$so_act         <-  so_C3$so_amp * so_C3$so_dur
so_C3$so_down_slope  <- (so_C3$mean_slope_to_trough_min_potential_per_second+so_C4$mean_slope_to_trough_min_potential_per_second)/2 
so_C3$so_up_slope    <- (so_C3$mean_slope_trough_to_up_max_potential_per_second+so_C4$mean_slope_trough_to_up_max_potential_per_second)/2 

so_params <- so_C3 %>% dplyr::select(datasetnum, so_count,so_dens, so_dur, so_amp, so_freq, so_act, so_down_slope, so_up_slope)
dataset_c <- merge(dataset_c, so_params, by = "datasetnum")

## Read in spindle-SO coupling data 
coupling_dataset_c  <- read.delim("dataset_c_coupling.csv", sep=",",dec= ".")

coupling_dataset_c <- mutate(coupling_dataset_c, test_channel = fct_recode(test_channel, "EEGC3" = "EEG_-C3/A2", 
                                                          "EEGC4" = "EEG_-C4/A1",
                                                          "EEGC3" = "EEG C3",
                                                          "EEGC4" = "EEG C4"))
coupling_dataset_c <- mutate(coupling_dataset_c, target_channel = fct_recode(target_channel, "EEGC3" = "EEG_-C3/A2", 
                                                                      "EEGC4" = "EEG_-C4/A1",
                                                                      "EEGC3" = "EEG C3",
                                                                      "EEGC4" = "EEG C4"))

# Split by channel to merge average below
coupling_C3 <- coupling_dataset_c[coupling_dataset_c$test_channel == "EEGC3", ]
coupling_C4 <- coupling_dataset_c[coupling_dataset_c$test_channel == "EEGC4", ]

## Count the amount of co-occurences
coupling_occur_c3 <- table(coupling_C3$event_files_number)
coupling_occur_c4 <- table(coupling_C4$event_files_number)

## Get sum of the amount of co-occurences
coupling_occur     <- (coupling_occur_c3+coupling_occur_c4)/2
coupling_occur_sum <- coupling_occur_c3+coupling_occur_c4  
coupling_occur     <- as.numeric(coupling_occur)
coupling_occur_sum <- as.numeric(coupling_occur_sum)

## Add to dataframe
dataset_c$coupling_occur     <- coupling_occur
dataset_c$coupling_occur_sum <- coupling_occur_sum

# Create the delay variable by subtracting the trough of the SO from the spindle
coupling_C3$delay_c3 <- coupling_C3$test_seconds_trough_max - coupling_C3$target_seconds_trough_max
coupling_C4$delay_c4 <- coupling_C4$test_seconds_trough_max - coupling_C4$target_seconds_trough_max

# Calculate mean and sd of delay per channel
mean_delay_c3 <- aggregate(delay_c3 ~ event_files_number, data=coupling_C3, FUN=mean)
sd_delay_c3   <- aggregate(delay_c3 ~ event_files_number, data=coupling_C3, FUN=sd)
mean_delay_c4 <- aggregate(delay_c4 ~ event_files_number, data=coupling_C4, FUN=mean)
sd_delay_c4   <- aggregate(delay_c4 ~ event_files_number, data=coupling_C4, FUN=sd)

# Average over channels
mean_delay     <- (mean_delay_c3+mean_delay_c4)/2 
sd_delay       <- (sd_delay_c3+sd_delay_c4)/2

# Add to main dataframe
dataset_c$mean_delay     <- mean_delay$delay_c3
dataset_c$sd_delay       <- sd_delay$delay_c3

## Supplemental analysis: get parameter information of matching versus mismatching events
## Based on MATLAB script (spisop_coupling_match_vs_mismatch.m) output; a .mat file containing matching and mismatching events
## Mismatching events
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

for (subnr in 1:98){ 
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
  
  i = i + 1     
} 


# Add to main dataframe
dataset_c$spin_dur_match    <- spin_dur_match
dataset_c$spin_dur_mismatch <- spin_dur_mismatch
dataset_c$spin_amp_match    <- spin_amp_match
dataset_c$spin_amp_mismatch <- spin_amp_mismatch
dataset_c$spin_freq_match   <- spin_freq_match
dataset_c$spin_freq_mismatch<- spin_freq_mismatch
dataset_c$spin_perc_match   <- spin_perc_match
dataset_c$so_dur_match      <- so_dur_match
dataset_c$so_dur_mismatch   <- so_dur_mismatch
dataset_c$so_amp_match      <- so_amp_match
dataset_c$so_amp_mismatch   <- so_amp_mismatch
dataset_c$so_perc_match     <- so_perc_match

dataset_c$spin_dur_diff  <- spin_dur_match - spin_dur_mismatch
dataset_c$spin_amp_diff  <- spin_amp_match - spin_amp_mismatch
dataset_c$spin_freq_diff <- spin_freq_match - spin_freq_mismatch
dataset_c$so_dur_diff    <- so_dur_match - so_dur_mismatch
dataset_c$so_amp_diff    <- so_amp_match - so_amp_mismatch

## IMPORTANT: Filter data ----
dataset_c <- dataset_c %>% dplyr::filter(dataset_c$pp_id != "19_1" &
                                         dataset_c$pp_id != "17_1" &  dataset_c$pp_id != "17_2" &
                                         dataset_c$pp_id != "27_1" &  dataset_c$pp_id != "27_2" & 
                                         dataset_c$pp_id != "43"   &  dataset_c$pp_id != "44"   & 
                                         dataset_c$pp_id != "45"   &  dataset_c$pp_id != "59"   &    
                                         dataset_c$pp_id != "61")

## Read in power data
# Merge with data
mergedData_long_dataset_c <- merge(dataset_c, power_dataset_c, by=c("datasetnum"))
 
 ## Reshape the data, per frequency
power_dataset_c_freq <- power_dataset_c %>%
   dplyr::select(datasetnum, channel, freq, mean_powerDensity_over_segments) %>%
   spread(key = freq, value = mean_powerDensity_over_segments)
 
 ## Average over the two channels
 power_dataset_c_freq_reshape <- aggregate(. ~ datasetnum, data = power_dataset_c_freq, mean)

 ## Merge with data
 mergedData_wide_dataset_c <- merge(dataset_c, power_dataset_c_freq_reshape, by = ("datasetnum"))

 # Transform to decibel
 mergedData_long_dataset_c$mean_densDB <- mergedData_long_dataset_c$mean_powerDensity_over_segments + 1
 mergedData_long_dataset_c$mean_densDB <- 10*log10(mergedData_long_dataset_c$mean_densDB)

 ## Create subgroups:
 ## control patient wide
 mergedData_wide_dataset_c_control_patient <- mergedData_wide_dataset_c %>%
   dplyr::filter(group == "control" | group == "patient_7")

 mergedData_wide_dataset_c_control_patient$group <- droplevels(mergedData_wide_dataset_c_control_patient$group)
 mergedData_wide_dataset_c_control_patient$group <- factor(mergedData_wide_dataset_c_control_patient$group)

 ## patient wide
 mergedData_wide_dataset_c_patients <- mergedData_wide_dataset_c %>%
   dplyr::filter(group == "patient_7" | group == "patient_28")

 mergedData_wide_dataset_c_patients$group <- droplevels(mergedData_wide_dataset_c_patients$group)
 mergedData_wide_dataset_c_patients$group <- factor(mergedData_wide_dataset_c_patients$group)

 ## control patient 28 wide
 mergedData_wide_dataset_c_conpat_med <- mergedData_wide_dataset_c %>%
   dplyr::filter(group == "control" | group == "patient_28")
 
 mergedData_wide_dataset_c_conpat_med$group <- droplevels(mergedData_wide_dataset_c_conpat_med$group)
 mergedData_wide_dataset_c_conpat_med$group <- factor(mergedData_wide_dataset_c_conpat_med$group)
 
 ## control patient long
 mergedData_long_dataset_c_control_patient <- mergedData_long_dataset_c %>%
   dplyr::filter(group == "control" | group == "patient_7")

 mergedData_long_dataset_c_control_patient$group <- droplevels(mergedData_long_dataset_c_control_patient$group)
 mergedData_long_dataset_c_control_patient$group <- factor(mergedData_long_dataset_c_control_patient$group)

 ## patient longs
 mergedData_long_dataset_c_patients <- mergedData_long_dataset_c %>%
   dplyr::filter(group == "patient_7" | group == "patient_28")

 mergedData_long_dataset_c_patients$group <- droplevels(mergedData_long_dataset_c_patients$group)
 mergedData_long_dataset_c_patients$group <- factor(mergedData_long_dataset_c_patients$group)
 
 ##  control patient 28 wide
 mergedData_long_dataset_c_conpat_med <- mergedData_long_dataset_c %>%
   dplyr::filter(group == "control" | group == "patient_28")
 
 mergedData_long_dataset_c_conpat_med$group <- droplevels(mergedData_long_dataset_c_conpat_med$group)
 mergedData_long_dataset_c_conpat_med$group <- factor(mergedData_long_dataset_c_conpat_med$group)

 ## get patient data in paired format: remove patient 19
 mergedData_wide_dataset_c_patients_pair <- mergedData_wide_dataset_c_patients %>% dplyr::filter(mergedData_wide_dataset_c_patients$pp_id != "19_1")
 mergedData_long_dataset_c_patients_pair <- mergedData_long_dataset_c_patients %>% dplyr::filter(mergedData_long_dataset_c_patients$pp_id != "19_1")


 ## Calculate t.test for each of the frequencies between groups CONTROLS AND PATIENTS
 column_dataset_c_06       <- which(colnames(mergedData_wide_dataset_c_control_patient)=="0.6")
 column_dataset_c_pats_06  <- which(colnames(mergedData_wide_dataset_c_patients_pair)=="0.6")
 column_dataset_c_conpat_med_06  <- which(colnames(mergedData_wide_dataset_c_conpat_med)=="0.6")
 
 listtests   <- lapply(mergedData_wide_dataset_c_control_patient[,c(column_dataset_c_06:(column_dataset_c_06+147))], function(x) t.test(x ~ mergedData_wide_dataset_c_control_patient$group, var.equal = FALSE))
 listtests_pats   <- lapply(mergedData_wide_dataset_c_patients_pair[,c(column_dataset_c_pats_06:(column_dataset_c_pats_06+147))], function(x) t.test(x ~ mergedData_wide_dataset_c_patients_pair$group, var.equal = TRUE, paired = TRUE))
 listtests_conpat_med   <- lapply(mergedData_wide_dataset_c_conpat_med[,c(column_dataset_c_conpat_med_06:(column_dataset_c_conpat_med_06+147))], function(x) t.test(x ~ mergedData_wide_dataset_c_conpat_med$group, var.equal = F, paired = F))
 
 ## Extract the p-value and make it a dataframe format
 listpvalues <- sapply(listtests, function(x) x[3])
 listpvalues <- as.data.frame(listpvalues)

 listpvalues_pats <- sapply(listtests_pats, function(x) x[3])
 listpvalues_pats <- as.data.frame(listpvalues_pats)
 
 listpvalues_conpat_med <- sapply(listtests_conpat_med, function(x) x[3])
 listpvalues_conpat_med <- as.data.frame(listpvalues_conpat_med)

 ## Reshape into long format
 listpvalues <- listpvalues %>%
   gather(key = freq , value = pvalue, starts_with("X")  )

 listpvalues_pats <- listpvalues_pats %>%
   gather(key = freq , value = pvalue, starts_with("X")  )
 
 listpvalues_conpat_med <- listpvalues_conpat_med %>%
   gather(key = freq , value = pvalue, starts_with("X")  )
 
 new_names <- seq(0.6, 30, by=0.2)
 listpvalues$Freq <- paste0(new_names)
 listpvalues$Freq<- as.numeric(listpvalues$Freq)
 listpvalues <- within(listpvalues, rm(freq))

 listpvalues_pats$Freq <- paste0(new_names)
 listpvalues_pats$Freq<- as.numeric(listpvalues_pats$Freq)
 listpvalues_pats <- within(listpvalues_pats, rm(freq))
 
 listpvalues_conpat_med$Freq <- paste0(new_names)
 listpvalues_conpat_med$Freq<- as.numeric(listpvalues_conpat_med$Freq)
 listpvalues_conpat_med <- within(listpvalues_conpat_med, rm(freq))

 data_long_pvalues <- merge(mergedData_long_dataset_c_control_patient, listpvalues, by.x=c("freq"), by.y=c("Freq"))
 data_long_pvalues_pats <- merge(mergedData_long_dataset_c_patients_pair, listpvalues_pats, by.x=c("freq"), by.y=c("Freq"))
 data_long_pvalues_conpat_med <- merge(mergedData_long_dataset_c_conpat_med, listpvalues_conpat_med, by.x=c("freq"), by.y=c("Freq"))
 
 ###### Permutation test -- patients controls
 freq_columns_start   <- which(colnames(mergedData_wide_dataset_c_control_patient)=="0.6")
 freq_columns <-  freq_columns_start:(freq_columns_start+147)

 pb = txtProgressBar(min = 0, max = length(freq_columns), initial = 0)
 perm_pvalue = rep(NA,length(freq_columns))
 index = 1


 for(i in freq_columns_start:(freq_columns_start+147)){
   ## Make new dataset
   perm_data <- mergedData_wide_dataset_c_control_patient %>% dplyr::select(i, group)

   # Split the data
   controls_freq <- perm_data[perm_data$group == "control",]
   patients_freq <- perm_data[perm_data$group == "patient_7",]

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

 ## Catch loop in cass a value is zero (will result in errors)
 for (i in 1:nrow(perm_pvalue)) {
   if ( perm_pvalue[i,1] == 0.000){
     perm_pvalue[i,1] <- 0.0001
   }
 }

 new_names <- seq(0.6, 30, by=0.2)
 perm_pvalue$Freq <- paste0(new_names)
 perm_pvalue$Freq<- as.numeric(perm_pvalue$Freq)

 dataset_c_long_pvalues_perm <- merge(data_long_pvalues, perm_pvalue, by.x=c("freq"), by.y=c("Freq"))

 ###### Permutation test -- patients paired
 freq_columns_start   <- which(colnames(mergedData_wide_dataset_c_patients_pair)=="0.6")
 freq_columns <-  freq_columns_start:(freq_columns_start+147)

 pb = txtProgressBar(min = 0, max = length(freq_columns), initial = 0)
 perm_pvalue_pats = rep(NA,length(freq_columns))
 index = 1


 for(i in freq_columns_start:(freq_columns_start+147)){
   ## Make new dataset
   perm_data <- mergedData_wide_dataset_c_patients_pair %>% dplyr::select(i, group)

   # Split the data
   controls_freq <- perm_data[perm_data$group == "patient_7",]
   patients_freq <- perm_data[perm_data$group == "patient_28",]

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
   perm_pvalue_pats[index] <- mean(abs(diffs) > abs(diff_obs))
   index = index + 1

   setTxtProgressBar(pb,i)
   print
 }

 ## Add them to the other dataframe
 perm_pvalue_pats <- as.data.frame(perm_pvalue_pats)

 ## Catch loop in cass a value is zero (will result in errors)
 for (i in 1:nrow(perm_pvalue_pats)) {
   if ( perm_pvalue_pats[i,1] == 0.000){
     perm_pvalue_pats[i,1] <- 0.0001
   }
 }

 new_names <- seq(0.6, 30, by=0.2)
 perm_pvalue_pats$Freq <- paste0(new_names)
 perm_pvalue_pats$Freq<- as.numeric(perm_pvalue_pats$Freq)

 dataset_c_long_pvalues_perm_pats <- merge(data_long_pvalues_pats, perm_pvalue_pats, by.x=c("freq"), by.y=c("Freq"))

 dataset_c_long_pvalues_perm_pats_10HZ <- dataset_c_long_pvalues_perm_pats %>% dplyr::filter(freq == 10)
 dataset_c_long_pvalues_perm_pats_10HZ %>% group_by(group) %>% dplyr::summarise(mean = mean(mean_densDB, na.rm = T))

 ## Alternative: use coin package to check calculation by hand:: is the same! :-)
 
 freq_columns_start   <- which(colnames(mergedData_wide_dataset_c_patients_pair)=="0.6")
 freq_columns <-  freq_columns_start:(freq_columns_start+147)
 
 pb = txtProgressBar(min = 0, max = length(freq_columns), initial = 0)
 perm_pvalue_pats = rep(NA,length(freq_columns))
 perm_pvalue_pats_ttest = rep(NA,length(freq_columns))
 perm_pvalue_pats_coin = rep(NA,length(freq_columns))
 perm_pvalue_pats_coin_paired = rep(NA,length(freq_columns))
 
 mergedData_wide_dataset_c_patients_pair$pp_id_match <- str_remove(mergedData_wide_dataset_c_patients_pair$pp_id, "_1")
 mergedData_wide_dataset_c_patients_pair$pp_id_match <- str_remove(mergedData_wide_dataset_c_patients_pair$pp_id_match, "_2")
 mergedData_wide_dataset_c_patients_pair$pp_id_match <- as.factor(mergedData_wide_dataset_c_patients_pair$pp_id_match)
 index = 1
 
 for(i in freq_columns_start:(freq_columns_start+147)){
   ## Make new dataset
   perm_data <- mergedData_wide_dataset_c_patients_pair %>% dplyr::select(i, group, pp_id_match)
   perm_data$pp_id_match <- as.factor(perm_data$pp_id_match)
   
   ## Get pvalue from coin package function
   a <- independence_test(perm_data[[1]] ~ group, data = perm_data) 
   b <- symmetry_test(perm_data[[1]] ~ group | pp_id_match, data = perm_data) 
   
   perm_pvalue_pats_coin[index] <- a@distribution@pvalue(a@statistic@teststatistic)
   perm_pvalue_pats_coin_paired[index] <- b@distribution@pvalue(b@statistic@teststatistic)
   index = index + 1
 }
 
 new_names <- seq(0.6, 30, by=0.2)
 perm_pvalue_pats_coin_paired <- as.data.frame(perm_pvalue_pats_coin_paired)
 perm_pvalue_pats_coin_paired$freq <- paste0(new_names)
 perm_pvalue_pats_coin_paired$freq<- as.numeric(perm_pvalue_pats_coin_paired$freq)
 
 dataset_c_long_pvalues_perm_pats <- merge(data_long_pvalues_pats, perm_pvalue_pats_coin_paired, by=c("freq"))
 
 
 ###### Permutation test --  controls patients 28d
 freq_columns_start   <- which(colnames(mergedData_wide_dataset_c_conpat_med)=="0.6")
 freq_columns <-  freq_columns_start:(freq_columns_start+147)
 
 pb = txtProgressBar(min = 0, max = length(freq_columns), initial = 0)
 perm_pvalue_conpat_med = rep(NA,length(freq_columns))
 index = 1
 
 
 for(i in freq_columns_start:(freq_columns_start+147)){
   ## Make new dataset
   perm_data <- mergedData_wide_dataset_c_conpat_med %>% dplyr::select(i, group)
   
   # Split the data
   controls_freq <- perm_data[perm_data$group == "control",]
   patients_freq <- perm_data[perm_data$group == "patient_28",]
   
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
   perm_pvalue_conpat_med[index] <- mean(abs(diffs) > abs(diff_obs))
   index = index + 1
   
   setTxtProgressBar(pb,i)
   print
 }
 
 ## Add them to the other dataframe
 perm_pvalue_conpat_med <- as.data.frame(perm_pvalue_conpat_med)
 
 ## Catch loop in cass a value is zero (will result in errors)
 for (i in 1:nrow(perm_pvalue_conpat_med)) {
   if ( perm_pvalue_conpat_med[i,1] == 0.000){
     perm_pvalue_conpat_med[i,1] <- 0.0001
   }
 }
 
 new_names <- seq(0.6, 30, by=0.2)
 perm_pvalue_conpat_med$Freq <- paste0(new_names)
 perm_pvalue_conpat_med$Freq<- as.numeric(perm_pvalue_conpat_med$Freq)
 
 dataset_c_long_pvalues_perm_conpat_med <- merge(data_long_pvalues_conpat_med, perm_pvalue_conpat_med, by.x=c("freq"), by.y=c("Freq"))
 
 
############################# Split datasets #############################
# reorder factor levels
dataset_c$group <- factor(dataset_c$group, levels = c("control", "patient_7", "patient_28"))

dataset_c_coupling_dens_info <- merge(dataset_c, coupling_dens_info, by = "datasetnum")
dataset_c_coupling_dens_info$NREM_min <- (dataset_c_coupling_dens_info$S2_without_MA_min +
                                            dataset_c_coupling_dens_info$S3_without_MA_min) /2 
dataset_c_coupling_dens_info$NREM_min <- (dataset_c_coupling_dens_info$NREM_min +
                                            dataset_c_coupling_dens_info$S4_without_MA_min) /2 

dataset_c$coupling_dens<- dataset_c_coupling_dens_info$coupling_occur  /  dataset_c_coupling_dens_info$NREM_min

## Split into three groups:
## controls and patients
dataset_c_conpat <- dataset_c %>% 
  dplyr::filter(group == "control" | group == "patient_7")
dataset_c_conpat$group <- droplevels(dataset_c_conpat$group)
dataset_c_conpat$group <- factor(dataset_c_conpat$group)

## patients 7 and 28
dataset_c_patients <- dataset_c %>% 
  dplyr::filter(group == "patient_7" | group == "patient_28")
dataset_c_patients$group <- droplevels(dataset_c_patients$group)
dataset_c_patients$group <- factor(dataset_c_patients$group)

## patients 7 and 28
dataset_c_conpat_med <- dataset_c %>% 
  dplyr::filter(group == "control" | group == "patient_28")
dataset_c_patients$group <- droplevels(dataset_c_patients$group)
dataset_c_patients$group <- factor(dataset_c_patients$group)

############################# Part II: Summary statistics and ttests #############################

## Demographics ----
dataset_c_demo_long <- dataset_c %>%
  gather(demo, score, age, gender, hamd_0, no_episodes)

dataset_c_demo_long %>%
  dplyr::select(demo, score, group) %>%
  group_by(group, demo) %>%
  get_summary_stats(score, type = "mean_se")

## Sleep architecture control - patient7d ---- 
dataset_c_conpat_hypno_long <- dataset_c_conpat %>%
  gather(sleepstage, score, S1_percent,S2_percent,S3andS4_perc,REM_percent,NREM_percent,Wake_after_sleep_onset_percent,Total_sleep_time_min, Sleep_Onset_min, SWS_onset_min, REM_onset_min)

dataset_c_conpat_hypno_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

hypnogram_ttests_dataset_c_conpat <- lapply(dataset_c_conpat[,c("age", "S1_percent", "S2_percent", "S3andS4_perc",  "REM_percent", "NREM_percent", "Wake_after_sleep_onset_percent", "Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min" )], function(x) t.test(x ~ dataset_c_conpat$group, var.equal = FALSE))
hypnogram_S4_ttests_dataset_c_conpat <- t.test(dataset_c_conpat$S4_percent ~ dataset_c_conpat$group) 

ttestBF(formula = S1_percent ~ group, data = dataset_c_conpat)
ttestBF(formula = S2_percent ~ group, data = dataset_c_conpat)
ttestBF(formula = S3andS4_perc ~ group, data = dataset_c_conpat)
ttestBF(formula = NREM_percent ~ group, data = dataset_c_conpat)
ttestBF(formula = REM_percent ~ group, data = dataset_c_conpat)
ttestBF(formula = Wake_after_sleep_onset_percent ~ group, data = dataset_c_conpat)
ttestBF(formula = Total_sleep_time_min ~ group, data = dataset_c_conpat)
ttestBF(formula = Sleep_Onset_min ~ group, data = dataset_c_conpat)
ttestBF(formula = SWS_onset_min ~ group, data = dataset_c_conpat)
ttestBF(formula = REM_onset_min ~ group, data = dataset_c_conpat)

## Sleep architecture patient7d - patient28d ---- 
dataset_c_patients_hypno_long <- dataset_c_patients %>%
  gather(sleepstage, score, S1_percent,S2_percent,S3andS4_perc,REM_percent,NREM_percent,Wake_after_sleep_onset_percent,Total_sleep_time_min, Sleep_Onset_min, SWS_onset_min, REM_onset_min)

dataset_c_patients_hypno_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

hypnogram_ttests_dataset_c_patients <- lapply(dataset_c_patients[,c("age", "S1_percent", "S2_percent", "S3andS4_perc",  "REM_percent", "NREM_percent", "Wake_after_sleep_onset_percent", "Total_sleep_time_min", "Sleep_Onset_min")], function(x) t.test(x ~ dataset_c_patients$group, paired = T))
hypnogram_S4_ttests_dataset_c_patients <- t.test(dataset_c_patients$S4_percent ~ dataset_c_patients$group, paired = T) 

## since there are NAs in SWS and REM onset: remove extra to get full pairs
dataset_c_patients_SWS_onset <- dataset_c_patients %>% drop_na(SWS_onset_min) ## Subject 4
dataset_c_patients_SWS_onset <- dataset_c_patients_SWS_onset %>% dplyr::filter(dataset_c_patients_SWS_onset$pp_id != "4_1" )
dataset_c_patients_REM_onset <- dataset_c_patients %>% drop_na(REM_onset_min) ## Subject 6
dataset_c_patients_REM_onset <- dataset_c_patients_REM_onset %>% dplyr::filter(dataset_c_patients_REM_onset$pp_id != "6_1" )

ttests_dataset_c_patients_SWS_onset <- t.test(dataset_c_patients_SWS_onset$SWS_onset_min ~dataset_c_patients_SWS_onset$group, paired = T)
ttests_dataset_c_patients_REM_onset <- t.test(dataset_c_patients_REM_onset$REM_onset_min ~dataset_c_patients_REM_onset$group, paired = T)

ttestBF(formula = S1_percent ~ group, data = dataset_c_patients)
ttestBF(formula = S2_percent ~ group, data = dataset_c_patients)
ttestBF(formula = S3andS4_perc ~ group, data = dataset_c_patients)
ttestBF(formula = NREM_percent ~ group, data = dataset_c_patients)
ttestBF(formula = REM_percent ~ group, data = dataset_c_patients)
ttestBF(formula = Wake_after_sleep_onset_percent ~ group, data = dataset_c_patients)
ttestBF(formula = Total_sleep_time_min ~ group, data = dataset_c_patients)
ttestBF(formula = Sleep_Onset_min ~ group, data = dataset_c_patients)
ttestBF(formula = SWS_onset_min ~ group, data = dataset_c_patients_SWS_onset)
ttestBF(formula = REM_onset_min ~ group, data = dataset_c_patients_REM_onset)


## Sleep architecture control - patient28d ---- 
dataset_c_conpat_med_hypno_long <- dataset_c_conpat_med %>%
  gather(sleepstage, score, S1_percent,S2_percent,S3andS4_perc,REM_percent,NREM_percent,Wake_after_sleep_onset_percent,Total_sleep_time_min, Sleep_Onset_min, SWS_onset_min, REM_onset_min)

dataset_c_conpat_med_hypno_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

hypnogram_ttests_dataset_c_conpat_med <- lapply(dataset_c_conpat_med[,c("age", "S1_percent", "S2_percent", "S3andS4_perc",  "REM_percent", "NREM_percent", "Wake_after_sleep_onset_percent", "Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min")], function(x) t.test(x ~ dataset_c_conpat_med$group, paired = F))
hypnogram_S4_ttests_dataset_c_conpat_med <- t.test(dataset_c_conpat_med$S4_percent ~ dataset_c_conpat_med$group) 

ttestBF(formula = S1_percent ~ group, data = dataset_c_conpat_med)
ttestBF(formula = S2_percent ~ group, data = dataset_c_conpat_med)
ttestBF(formula = S3andS4_perc ~ group, data = dataset_c_conpat_med)
ttestBF(formula = NREM_percent ~ group, data = dataset_c_conpat_med)
ttestBF(formula = REM_percent ~ group, data = dataset_c_conpat_med)
ttestBF(formula = Wake_after_sleep_onset_percent ~ group, data = dataset_c_conpat_med)
ttestBF(formula = Total_sleep_time_min ~ group, data = dataset_c_conpat_med)
ttestBF(formula = Sleep_Onset_min ~ group, data = dataset_c_conpat_med)

dataset_c_conpat_med_SWS <- dataset_c_conpat_med%>% drop_na(SWS_onset_min)
ttestBF(formula = SWS_onset_min ~ group, data = dataset_c_conpat_med_SWS)
dataset_c_conpat_med_REM <- dataset_c_conpat_med%>% drop_na(REM_onset_min)
ttestBF(formula = REM_onset_min ~ group, data = dataset_c_conpat_med_REM)

dataset_c_pat_med7 <- dataset_c_patients %>% filter(group ==  "patient_7")
cor.test(dataset_c_pat_med7$S1_percent, dataset_c_pat_med7$REM_percent)

dataset_c_pat_med28 <- dataset_c_patients %>% filter(group ==  "patient_28")
cor.test(dataset_c_pat_med28$S1_percent, dataset_c_pat_med28$REM_percent)

## Spindle data control - patient7d ---- 
dataset_c_conpat_spindle_long <- dataset_c_conpat %>%
  gather(spindle, score, spin_dens, spin_count,spin_amp, spin_dur, spin_freq)

dataset_c_conpat_spindle_long %>%
  dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%
  get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_c_conpat <- lapply(dataset_c_conpat[,c("spin_dens","spin_dur","spin_amp", "spin_count", "spin_centerfreq" , "spin_freq")], function(x) t.test(x ~ dataset_c_conpat$group, var.equal = FALSE))

## Spindle data patient7d - patient28d ---- 
dataset_c_patients_spindle_long <- dataset_c_patients %>%
  gather(spindle, score, spin_dens, spin_count,spin_amp, spin_dur, spin_freq)

dataset_c_patients_spindle_long %>%
  dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%
  get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_c_patients <- lapply(dataset_c_patients[,c("spin_dens","spin_dur","spin_amp", "spin_count", "spin_centerfreq" , "spin_freq")], function(x) t.test(x ~ dataset_c_patients$group, paired = T))

## Spindle data control - patient28d ---- 
dataset_c_conpat_med_spindle_long <- dataset_c_conpat_med %>%
  gather(spindle, score, spin_dens, spin_count,spin_amp, spin_dur,  spin_freq)

dataset_c_conpat_med_spindle_long %>%
  dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%
  get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_c_conpat_med <- lapply(dataset_c_conpat_med[,c("spin_dens","spin_dur","spin_amp", "spin_count", "spin_centerfreq", "spin_freq" )], function(x) t.test(x ~ dataset_c_conpat_med$group, paired = F))

## Calculate bayes factor
ttestBF(formula = spin_dens ~ group, data = dataset_c_conpat)
ttestBF(formula = spin_count ~ group, data = dataset_c_conpat)
ttestBF(formula = spin_amp ~ group, data = dataset_c_conpat)
ttestBF(formula = spin_freq ~ group, data = dataset_c_conpat)
ttestBF(formula = spin_dur ~ group, data = dataset_c_conpat)
ttestBF(formula = spin_dens ~ group, data = dataset_c_patients)
ttestBF(formula = spin_count ~ group, data = dataset_c_patients)
ttestBF(formula = spin_amp ~ group, data = dataset_c_patients)
ttestBF(formula = spin_freq ~ group, data = dataset_c_patients)
ttestBF(formula = spin_dur ~ group, data = dataset_c_patients)
ttestBF(formula = spin_dens ~ group, data = dataset_c_conpat_med)
ttestBF(formula = spin_count ~ group, data = dataset_c_conpat_med)
ttestBF(formula = spin_amp ~ group, data = dataset_c_conpat_med)
ttestBF(formula = spin_freq ~ group, data = dataset_c_conpat_med)
ttestBF(formula = spin_dur ~ group, data = dataset_c_conpat_med)

## SO data control - patient7d ---- 
dataset_c_conpat_so_long <- dataset_c_conpat %>%
  gather(so, score, so_dens, so_count, so_amp, so_dur, so_freq)

dataset_c_conpat_so_long %>%
  dplyr::select(so, score, group) %>%
  group_by(group, so) %>%
  get_summary_stats(score, type = "mean_se")

so_ttests_dataset_c_conpat <- lapply(dataset_c_conpat[,c("so_dens","so_dur","so_amp","so_count" ,"so_freq")], function(x) t.test(x ~ dataset_c_conpat$group, var.equal = FALSE))

## SO data patient7d - patient28d ---- 
dataset_c_patients_so_long <- dataset_c_patients %>%
  gather(so, score, so_dens,so_count,so_amp, so_dur, so_freq)

dataset_c_patients_so_long %>%
  dplyr::select(so, score, group) %>%
  group_by(group, so) %>%
  get_summary_stats(score, type = "mean_se")

so_ttests_dataset_c_patients <- lapply(dataset_c_patients[,c("so_dens","so_dur","so_amp" ,"so_count","so_freq")], function(x) t.test(x ~ dataset_c_patients$group, paired = T))

## SO data control - patient28d ---- 
dataset_c_conpat_med_so_long <- dataset_c_conpat_med %>%
  gather(so, score, so_dens,so_count,so_amp, so_dur, so_freq)

dataset_c_conpat_med_so_long %>%
  dplyr::select(so, score, group) %>%
  group_by(group, so) %>%
  get_summary_stats(score, type = "mean_se")

so_ttests_dataset_c_conpat_med <- lapply(dataset_c_conpat_med[,c("so_dens","so_dur","so_amp","so_count","so_freq" )], function(x) t.test(x ~ dataset_c_conpat_med$group, paired = F))

ttestBF(formula = so_dens ~ group, data = dataset_c_conpat)
ttestBF(formula = so_count ~ group, data = dataset_c_conpat)
ttestBF(formula = so_amp ~ group, data = dataset_c_conpat)
ttestBF(formula = so_freq ~ group, data = dataset_c_conpat)
ttestBF(formula = so_dur ~ group, data = dataset_c_conpat)
ttestBF(formula = so_dens ~ group, data = dataset_c_patients)
ttestBF(formula = so_count ~ group, data = dataset_c_patients)
ttestBF(formula = so_amp ~ group, data = dataset_c_patients)
ttestBF(formula = so_freq ~ group, data = dataset_c_patients)
ttestBF(formula = so_dur ~ group, data = dataset_c_patients)
ttestBF(formula = so_dens ~ group, data = dataset_c_conpat_med)
ttestBF(formula = so_count ~ group, data = dataset_c_conpat_med)
ttestBF(formula = so_amp ~ group, data = dataset_c_conpat_med)
ttestBF(formula = so_freq ~ group, data = dataset_c_conpat_med)
ttestBF(formula = so_dur ~ group, data = dataset_c_conpat_med)

## Coupling data  control - patient7d ---- 
dataset_c_conpat_coupling_long <- dataset_c_conpat %>%
  gather(coupling, score, coupling_occur,mean_delay,sd_delay)

dataset_c_conpat_coupling_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_ttests_dataset_c_conpat <- lapply(dataset_c_conpat[,c("coupling_occur","mean_delay","sd_delay","coupling_dens" )], function(x) t.test(x ~ dataset_c_conpat$group, var.equal = FALSE))

## Coupling data patient7d - patient28d  ---- 
dataset_c_patients_coupling_long <- dataset_c_patients %>%
  gather(coupling, score, coupling_occur,mean_delay,sd_delay)

dataset_c_patients_coupling_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_ttests_dataset_c_patients <- lapply(dataset_c_patients[,c("coupling_occur","mean_delay","sd_delay","coupling_dens"   )], function(x) t.test(x ~ dataset_c_patients$group, paired = T))

## Coupling data  control - patient28d ---- 
dataset_c_conpat_med_coupling_long <- dataset_c_conpat_med %>%
  gather(coupling, score, coupling_occur,mean_delay,sd_delay)

dataset_c_conpat_med_coupling_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_ttests_dataset_c_conpat_med <- lapply(dataset_c_conpat_med[,c("coupling_occur","mean_delay","sd_delay","coupling_dens"   )], function(x) t.test(x ~ dataset_c_conpat_med$group, paired = F))

## Calculate bayes factor
ttestBF(formula = coupling_occur ~ group, data = dataset_c_conpat)
ttestBF(formula = mean_delay ~ group, data = dataset_c_conpat)
ttestBF(formula = sd_delay ~ group, data = dataset_c_conpat)
ttestBF(formula = coupling_occur ~ group, data = dataset_c_patients)
ttestBF(formula = mean_delay ~ group, data = dataset_c_patients)
ttestBF(formula = sd_delay ~ group, data = dataset_c_patients)
ttestBF(formula = coupling_occur ~ group, data = dataset_c_conpat_med)
ttestBF(formula = mean_delay ~ group, data = dataset_c_conpat_med)
ttestBF(formula = sd_delay ~ group, data = dataset_c_conpat_med)

## Coupling matching vs nonmatching controls - patients
dataset_c_conpat_coupling_matchspin_long <- dataset_c_conpat %>%
  gather(coupling, score, spin_amp_match, spin_dur_match,spin_freq_match,spin_perc_match,
         spin_amp_mismatch, spin_dur_mismatch,spin_freq_mismatch)

dataset_c_conpat_coupling_matchspin_long <- dataset_c_conpat %>%
  gather(coupling, score, spin_dur_match, spin_dur_mismatch, spin_dur_diff, so_dur_match, so_dur_mismatch, so_dur_diff)

dataset_c_conpat_coupling_matchspin_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

dataset_c_conpat_coupling_matchspin_long <- dataset_c_conpat  %>%
  gather(coupling, score, spin_amp_match, spin_dur_match,spin_freq_match,spin_perc_match,
         so_amp_match, so_dur_match,so_perc_match)

dataset_c_conpat_coupling_matchspin_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

dataset_c_conpat_med_coupling_matchspin_long <- dataset_c_conpat_med  %>%
  gather(coupling, score, spin_amp_match, spin_dur_match,spin_freq_match,spin_perc_match,
         so_amp_match, so_dur_match,so_perc_match)

dataset_c_conpat_med_coupling_matchspin_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

lapply(dataset_c_conpat[,c("spin_amp_diff","spin_dur_diff","spin_freq_diff","so_amp_diff", "so_dur_diff" )], function(x) t.test(x, mu = 0, var.equal = FALSE))

coupling_diff_ttests_dataset_c_conpat <- lapply(dataset_c_conpat[,c("spin_amp_diff","spin_dur_diff","spin_freq_diff","spin_perc_match","so_amp_diff", "so_dur_diff", "so_perc_match" )], function(x) t.test(x ~ dataset_c_conpat$group, var.equal = FALSE))
coupling_diff_ttests_dataset_c_patients <- lapply(dataset_c_patients[,c("spin_amp_diff","spin_dur_diff","spin_freq_diff","spin_perc_match","so_amp_diff", "so_dur_diff", "so_perc_match" )], function(x) t.test(x ~ dataset_c_patients$group, var.equal = FALSE, paired = T))
coupling_diff_ttests_dataset_c_conpat_med <- lapply(dataset_c_conpat_med[,c("spin_amp_diff","spin_dur_diff","spin_freq_diff","spin_perc_match","so_amp_diff", "so_dur_diff", "so_perc_match" )], function(x) t.test(x ~ dataset_c_conpat_med$group, var.equal = FALSE))

coupling_match_ttests_dataset_c_conpat <- lapply(dataset_c_conpat[,c("spin_amp_match","spin_dur_match","spin_freq_match","so_amp_match", "so_dur_match")], function(x) t.test(x ~ dataset_c_conpat$group, var.equal = FALSE))
coupling_match_ttests_dataset_c_patients <- lapply(dataset_c_patients[,c("spin_amp_match","spin_dur_match","spin_freq_match","so_amp_match", "so_dur_match")], function(x) t.test(x ~ dataset_c_patients$group, var.equal = FALSE, paired =T))
coupling_match_ttests_dataset_c_conpat_med <- lapply(dataset_c_conpat_med[,c("spin_amp_match","spin_dur_match","spin_freq_match","so_amp_match", "so_dur_match" )], function(x) t.test(x ~ dataset_c_conpat_med$group, var.equal = FALSE))

ttestBF(formula = spin_amp_match ~ group, data  = dataset_c_conpat)
ttestBF(formula = spin_dur_match ~ group, data  = dataset_c_conpat)
ttestBF(formula = spin_freq_match ~ group, data = dataset_c_conpat)
ttestBF(formula = so_amp_match ~ group, data    = dataset_c_conpat)
ttestBF(formula = so_dur_match ~ group, data    = dataset_c_conpat)
ttestBF(formula = spin_amp_match ~ group, data  = dataset_c_patients)
ttestBF(formula = spin_dur_match ~ group, data  = dataset_c_patients)
ttestBF(formula = spin_freq_match ~ group, data = dataset_c_patients)
ttestBF(formula = so_amp_match ~ group, data    = dataset_c_patients)
ttestBF(formula = so_dur_match ~ group, data    = dataset_c_patients)
ttestBF(formula = spin_amp_match ~ group, data  = dataset_c_conpat_med)
ttestBF(formula = spin_dur_match ~ group, data  = dataset_c_conpat_med)
ttestBF(formula = spin_freq_match ~ group, data = dataset_c_conpat_med)
ttestBF(formula = so_amp_match ~ group, data    = dataset_c_conpat_med)
ttestBF(formula = so_dur_match ~ group, data    = dataset_c_conpat_med)

ttestBF(formula = spin_amp_diff ~ group, data  = dataset_c_conpat)
ttestBF(formula = spin_dur_diff ~ group, data  = dataset_c_conpat)
ttestBF(formula = spin_freq_diff ~ group, data = dataset_c_conpat)
ttestBF(formula = so_amp_diff ~ group, data    = dataset_c_conpat)
ttestBF(formula = so_dur_diff ~ group, data    = dataset_c_conpat)
ttestBF(formula = spin_amp_diff ~ group, data  = dataset_c_patients)
ttestBF(formula = spin_dur_diff ~ group, data  = dataset_c_patients)
ttestBF(formula = spin_freq_diff ~ group, data = dataset_c_patients)
ttestBF(formula = so_amp_diff ~ group, data    = dataset_c_patients)
ttestBF(formula = so_dur_diff ~ group, data    = dataset_c_patients)
ttestBF(formula = spin_amp_diff ~ group, data  = dataset_c_conpat_med)
ttestBF(formula = spin_dur_diff ~ group, data  = dataset_c_conpat_med)
ttestBF(formula = spin_freq_diff ~ group, data = dataset_c_conpat_med)
ttestBF(formula = so_amp_diff ~ group, data    = dataset_c_conpat_med)
ttestBF(formula = so_dur_diff ~ group, data    = dataset_c_conpat_med)

dataset_c_conpat_med_long_spin_dur <- dataset_c_conpat_med %>%
  gather(coupling, score, spin_dur_match, spin_dur_mismatch, spin_dur_diff)

dataset_c_conpat_med_long_spin_dur %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

dataset_c_conpat%>%
  gather(match, score, spin_amp_diff, spin_dur_diff,spin_freq_diff, so_amp_diff, so_dur_diff) %>%
  dplyr::select(match, score, group) %>%
  group_by(group, match) %>% get_summary_stats(score, type = "mean_se")

dataset_c_patients %>%
  gather(match, score, spin_amp_diff, spin_dur_diff,spin_freq_diff, so_amp_diff, so_dur_diff) %>%
  dplyr::select(match, score, group) %>%
  group_by(group, match) %>% get_summary_stats(score, type = "mean_se")

dataset_c_conpat_med %>%
  gather(match, score, spin_amp_diff, spin_dur_diff,spin_freq_diff, so_amp_diff, so_dur_diff) %>%
  dplyr::select(match, score, group) %>%
  group_by(group, match) %>% get_summary_stats(score, type = "mean_se")


## Part III: Plots ---- 
## Combine data for plots
dataset_c_long_pvalues_perm_select <- dataset_c_long_pvalues_perm  %>% dplyr::select(pnum, datasetnum, group, freq, mean_densDB, pvalue)
dataset_c_long_pvalues_perm_pats_select   <- dataset_c_long_pvalues_perm_pats   %>% dplyr::select(pnum, datasetnum,  group, freq, mean_densDB, pvalue)
dataset_c_long_pvalues_perm_conpat_med_select   <- dataset_c_long_pvalues_perm_conpat_med  %>% dplyr::select(pnum, datasetnum, group, freq, mean_densDB, pvalue)

dataset_c_long_pvalues_perm_control <- dataset_c_long_pvalues_perm_select %>% filter(group == "control")
dataset_c_long_pvalues_perm_pat7    <- dataset_c_long_pvalues_perm_select %>% filter(group == "patient_7")
dataset_c_long_pvalues_perm_pat28   <- dataset_c_long_pvalues_perm_pats_select %>% filter(group == "patient_28")

dataset_c_long_pvalues_perm_three <- rbind(dataset_c_long_pvalues_perm_control, dataset_c_long_pvalues_perm_pat7, dataset_c_long_pvalues_perm_pat28)

## Power plot control - patient7d
 freq_spec_dataset_c_conpat <- dataset_c_long_pvalues_perm %>%
   ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
   theme_classic()+ my_theme +
   stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
   scale_color_manual(values=color_dataset_c, labels = c("Controls (C)", "Medicated (C)")) +  
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
 
 perm_plot_dataset_c_conpat <- ggplot(dataset_c_long_pvalues_perm, aes(x=freq, y=abs(log10(perm_pvalue)), colour="#000000")) +
   theme_classic()+ my_theme +
   stat_summary(fun.y="mean", geom="bar", colour="black")+ 
   xlab("\nFrequency") + ylab("P-values (con - 7d)") +
   theme(       
     axis.text.y = element_text( color="#000000"),
     axis.title.x = element_blank(),
     axis.line.x  = element_blank(),
     axis.text.x  = element_blank(),
     axis.ticks.x = element_blank(),
     legend.position='none')+ 
   geom_hline(yintercept = abs(log10(0.05)))+
   scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
   scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 
 NULL
 
 ## Power plot patient7d - patient28d
  freq_spec_dataset_c_pats <- dataset_c_long_pvalues_perm_pats %>%
   ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
   theme_classic()+my_theme +
   stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
   scale_color_manual(values=color_dataset_c_pat, labels = c("Medicated - 7days (C)", "Medicated - 28days (C)")) +  
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
 
 perm_plot_dataset_c_pats <- ggplot(dataset_c_long_pvalues_perm_pats, aes(x=freq, y=abs(log10(perm_pvalue_pats_coin_paired)), colour="#000000")) +
   theme_classic()+ my_theme +
   stat_summary(fun.y="mean", geom="bar", colour="black")+ 
   xlab("\nFrequency") + ylab("P-values (7d - 28d)") + 
   theme(       
     axis.text.y = element_text( color="#000000"),
     #axis.text.x = element_text(color="#000000"),
     axis.title.x = element_blank(),
     axis.line.x  = element_blank(),
     axis.text.x  = element_blank(),
     axis.ticks.x = element_blank(), 
     legend.position='none')+ 
   geom_hline(yintercept = abs(log10(0.05)))+
   scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
   scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 
 NULL
 
 ## Power plot control - patient7d - patient28d
 freq_spec_dataset_c_three <- dataset_c_long_pvalues_perm_three %>%
   ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
   theme_classic()+my_theme +
   stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
   scale_color_manual(values=color_dataset_c_all, labels = c("Controls - 7days (C)", "Medicated - 7days (C)", "Medicated - 28days (C)")) +  
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
 
 perm_plot_dataset_c_conpat_med <- ggplot(dataset_c_long_pvalues_perm_conpat_med, aes(x=freq, y=abs(log10(perm_pvalue_conpat_med)), colour="#000000")) +
   theme_classic()+ my_theme +
   stat_summary(fun.y="mean", geom="bar", colour="black")+ 
   xlab("\nFrequency") + ylab("P-values (con - 28d)") +
   theme(       
     axis.text.y = element_text( color="#000000"),
     axis.text.x = element_text(color="#000000"), 
     legend.position='none')+ 
   geom_hline(yintercept = abs(log10(0.05)))+
   scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
   scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 

## Spindles plot control - patient7d - patient28d
plot_rain_spidens_dataset_c <- ggplot(dataset_c,aes(x=group,y=spin_dens,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_dens),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle density')+ xlab('')+ ylim(0,4) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_c_all) + 
  scale_fill_manual(values=color_dataset_c_all) + 
  scale_x_discrete(labels = c("Controls (C)", "Medicated 7d (C)", "Medicated 28d (C)")) +
  ggtitle("Group difference in spindle density") +
  geom_signif(annotation=formatC("*", digits=1), y_position=3.6,xmin=1.2, xmax=3.2, 
              textsize = 6, colour="BLACK") 

plot_rain_spinamp_dataset_c <- ggplot(dataset_c,aes(x=group,y=spin_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle amplitude')+ xlab('')+ ylim(0,55) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_c_all) + 
  scale_fill_manual(values=color_dataset_c_all) + 
  scale_x_discrete(labels = c("Controls (C)", "Medicated 7d (C)", "Medicated 28d (C)")) +
  ggtitle("Group difference in spindle amplitude")

## SO plot control - patient7d - patient28d
plot_rain_soamp_dataset_c <- ggplot(dataset_c,aes(x=group,y=so_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO amplitude')+ xlab('')+ ylim(0,300) +
  raincloud_theme + my_theme + 
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_c_all) + 
  scale_fill_manual(values=color_dataset_c_all) + 
  scale_x_discrete(labels = c("Controls (C)", "Medicated 7d (C)", "Medicated 28d (C)")) +
  ggtitle("Group difference in SO amplitude")  +
  geom_signif(annotation=formatC("*", digits=1), y_position=270,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK")  +
  geom_signif(annotation=formatC("*", digits=1), y_position=290,xmin=1.2, xmax=3.2, 
              textsize = 6, colour="BLACK") 

plot_rain_sodur_dataset_c <- ggplot(dataset_c,aes(x=group,y=so_dur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_dur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO duration')+ xlab('')+ ylim(1,1.6) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_c_all) + 
  scale_fill_manual(values=color_dataset_c_all) + 
  scale_x_discrete(labels = c("Controls (C)", "Medicated 7d (C)", "Medicated 28d (C)")) +
  ggtitle("Group difference in SO duration") +
  geom_signif(annotation=formatC("*", digits=1), y_position=1.58,xmin=1.2, xmax=3.2, 
              textsize = 6, colour="BLACK") 

plot_rain_sodens_dataset_c <- ggplot(dataset_c,aes(x=group,y=so_dens,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_dens),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO density')+ xlab('')+ ylim(0,3) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_c_all) + 
  scale_fill_manual(values=color_dataset_c_all) + 
  scale_x_discrete(labels = c("Controls (C)", "Medicated 7d (C)", "Medicated 28d (C)")) +
  ggtitle("Group difference in SO density")

## Coupling plot control - patient7d - patient28d
plot_rain_amount_coup_dataset_c <- ggplot(dataset_c,aes(x=group,y=coupling_occur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=coupling_occur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO-fast spindles')+ xlab('')+ ylim(0,300) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_c_all) + 
  scale_fill_manual(values=color_dataset_c_all) + 
  scale_x_discrete(labels = c("Controls (C)", "Medicated 7d (C)", "Medicated 28d (C)")) +
  ggtitle("Group difference in amount of SO-fast spindle couplings")

plot_rain_mean_delay_dataset_c <- ggplot(dataset_c,aes(x=group,y=mean_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=mean_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Mean delay')+ xlab('')+ ylim(0.2,0.85) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_c_all) + 
  scale_fill_manual(values=color_dataset_c_all) + 
  scale_x_discrete(labels = c("Controls (C)", "Medicated 7d (C)", "Medicated 28d (C)")) +
  ggtitle("Group difference in mean delay") 

plot_rain_sd_delay_dataset_c <- ggplot(dataset_c,aes(x=group,y=sd_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=sd_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Delay dispersion (SD)')+ xlab('')+ ylim(0.15,0.45) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_c_all) + 
  scale_fill_manual(values=color_dataset_c_all) + 
  ggtitle("Group difference in delay dispersion (SD)") +
  scale_x_discrete(labels = c("Controls (C)", "Medicated 7d (C)", "Medicated 28d (C)")) +
  geom_signif(annotation=formatC("*", digits=1), y_position=0.4,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK") +
  geom_signif(annotation=formatC("*", digits=1), y_position=0.43,xmin=1.2, xmax=3.2, 
              textsize = 6, colour="BLACK") 

## Add the plots on top of eachother
freq_plot_dataset_c       <- ggarrange(freq_spec_dataset_c_conpat, perm_plot_dataset_c_conpat, heights = c(2, 1), ncol = 1, nrow = 2)
freq_plot_pats_dataset_c  <- ggarrange(freq_spec_dataset_c_pats, perm_plot_dataset_c_pats, heights = c(2, 1), ncol = 1, nrow = 2)
freq_plot_three_dataset_c  <- ggarrange(freq_spec_dataset_c_three, perm_plot_dataset_c_conpat, perm_plot_dataset_c_pats, perm_plot_dataset_c_conpat_med, heights = c(3, 1, 1, 1), ncol = 1, nrow = 4)

spin_plot_dataset_c     <- ggarrange(plot_rain_spidens_dataset_c, plot_rain_spinamp_dataset_c,nrow = 2, labels = c("A", "B") )
so_plot_dataset_c       <- ggarrange(plot_rain_soamp_dataset_c, plot_rain_sodur_dataset_c,nrow = 2, labels = c("A", "B") )
coupling_plot_dataset_c <- ggarrange(plot_rain_amount_coup_dataset_c, plot_rain_mean_delay_dataset_c, plot_rain_sd_delay_dataset_c, nrow = 3, labels = c("A", "B", "C") )
 
## Save the images
#ggsave("figs/freq_plot_dataset_c.png", freq_plot_dataset_c, device = "png", dpi = 300, width = 7, height = 10)
#ggsave("figs/freq_plot_pats_dataset_c.png", freq_plot_pats_dataset_c, device = "png", dpi = 300, width = 7, height = 10)
#ggsave("figs/freq_plot_three_dataset_c.png", freq_plot_three_dataset_c, device = "png", dpi = 300, width = 7, height = 15)
#ggsave("figs/spin_plot_dataset_c.png", spin_plot_dataset_c, device = "png", dpi = 300, width = 7, height = 10)
#ggsave("figs/so_plot_dataset_c_v2.png", so_plot_dataset_c, device = "png", dpi = 300, width = 7, height = 10)
#ggsave("figs/coupling_plot_dataset_c_v2.png", coupling_plot_dataset_c, device = "png", dpi = 300, width = 7, height = 10)

## Part IV: Additional analyses ----
## Depression severity ----
dataset_c_hamd_long <- dataset_c %>% filter(group=="patient_7") %>% gather(hamd, score, hamd_0, hamd_week1, hamd_week4)
dataset_c_hamd_long %>% dplyr::select(hamd, score, group) %>% group_by(hamd) %>%  get_summary_stats(score, type = "mean_se")

## split into parts
dataset_c_hamd_long_base_7  <- dataset_c_hamd_long %>% filter(hamd == "hamd_0" | hamd == "hamd_week1")
dataset_c_hamd_long_7_28    <- dataset_c_hamd_long %>% filter(hamd == "hamd_week1" | hamd == "hamd_week4")
dataset_c_hamd_long_base_7_ttest  <- t.test(score ~ hamd, data = dataset_c_hamd_long_base_7)
dataset_c_hamd_long_7_28_ttest <- t.test(score ~ hamd, data = dataset_c_hamd_long_7_28)
dataset_c_hamd_long$hamd <- as.factor(dataset_c_hamd_long$hamd)



