## Script to analyse Dataset B
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

## IMPORTANT: 
## We have data for 40 controls, 40 unmedicated and 38 medicated
## The analyses were initially run on all 118 datasets, but subsequent statistics will be performed with excluding the two matching unmedicated patients as well
## Missing are LP_104_2 and LP_107_2  //  remove -> LP_104_1 and LP_107_1

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
source("funcs/combine_pow_output_dataset_b.R")
source("funcs/get_hipp_size_values.R")

## Main script  ----
## Part I: Read in and prepare data ----
## Read in behavioural data
dataset_b <- read_excel("dataset_b_subject_info.xlsx") 
dataset_b$group <- as.factor(dataset_b$group)
dataset_b$PSQI_sum <- as.numeric(as.character(dataset_b$PSQI_sum))

## Order by datasetname
dataset_b <- dataset_b[order(dataset_b$datasetnum),] 

## Read in sleep table
hypval_dataset_b  <- read.delim("dataset_b_hypvals.csv", sep=",",dec= ".") 

coupling_dens_info <- hypval_dataset_b %>%
  dplyr::select(datasetnum, S2_without_MA_min,
                S3_without_MA_min,S4_without_MA_min)

hypval_dataset_b <- hypval_dataset_b %>% 
  dplyr::select(datasetnum, Total_sleep_time_min:SWS_onset_min, 
                S2_onset_min:Wake_after_sleep_onset_min,SWS_min, 
                S1_percent:Wake_after_sleep_onset_percent) 

dataset_b <- merge(dataset_b, hypval_dataset_b, by = "datasetnum", all.x = T)

# Combine S3 and S4 to fit AASM scoring
dataset_b$S3andS4_min   <- dataset_b$S3_min     + dataset_b$S4_min
dataset_b$S3andS4_perc  <- dataset_b$S3_percent + dataset_b$S4_percent
dataset_b$NREM_percent  <- dataset_b$S2_percent + dataset_b$S3andS4_perc

## Read in spindle data 
spindle_dataset_b <- read.delim("dataset_b_spindles.csv", sep=",",dec= ".")   

## Recode to match all different datasets
spindle_dataset_b <- mutate(spindle_dataset_b, channel = fct_recode(channel,"EEGC3" = "C3-TP10", "EEGC4" = "C4-TP9"))

# Split by channel to merge average below
spindle_C3 <- spindle_dataset_b[spindle_dataset_b$channel == "EEGC3", ]
spindle_C4 <- spindle_dataset_b[spindle_dataset_b$channel == "EEGC4", ]

# Create spindle parameters by averaging over C3 and C4 
spindle_C3$spin_dens       <- (spindle_C3$density_per_epoch + spindle_C4$density_per_epoch)/2 
spindle_C3$spin_centerfreq <- (spindle_C3$used_center_freq+spindle_C4$used_center_freq)/2
spindle_C3$spin_count      <- (spindle_C3$count+spindle_C4$count)/2  
spindle_C3$spin_amp        <- (spindle_C3$mean_amplitude_trough2peak_potential+spindle_C4$mean_amplitude_trough2peak_potential)/2  
spindle_C3$spin_dur        <- (spindle_C3$mean_duration_seconds+spindle_C4$mean_duration_seconds)/2  
spindle_C3$spin_freq       <- (spindle_C3$mean_frequency_by_mean_pk_trgh_cnt_per_dur+spindle_C4$mean_frequency_by_mean_pk_trgh_cnt_per_dur)/2  
spindle_C3$spin_act        <- spindle_C3$spin_amp * spindle_C3$spin_dur

spindle_params <- spindle_C3 %>% dplyr::select(datasetnum, spin_dens,spin_centerfreq, spin_count, spin_amp, spin_dur, spin_freq, spin_act)
dataset_b <- merge(dataset_b, spindle_params, by = "datasetnum")

## Read in SO data 
so_dataset_b  <- read.delim("dataset_b_so.csv", sep=",",dec= ".")   
so_dataset_b <- mutate(so_dataset_b, channel = fct_recode(channel,"EEGC3" = "C3-TP10","EEGC4" = "C4-TP9"))

# Split by channel to merge average below
so_C3 <- so_dataset_b[so_dataset_b$channel == "EEGC3", ]
so_C4 <- so_dataset_b[so_dataset_b$channel == "EEGC4", ]

# Create spindle parameters by averaging over C3 and C4
so_C3$so_count       <- (so_C3$count + so_C4$count)/2
so_C3$so_dens        <- (so_C3$density_per_epoch + so_C4$density_per_epoch)/2
so_C3$so_dur         <- (so_C3$mean_duration_seconds + so_C4$mean_duration_seconds)/2
so_C3$so_amp         <- (so_C3$mean_amplitude_peak2trough_potential+so_C4$mean_amplitude_peak2trough_potential)/2
so_C3$so_freq        <- (so_C3$mean_frequency_by_duration+so_C4$mean_frequency_by_duration)/2
so_C3$so_act         <-  so_C3$so_amp * so_C3$so_dur
so_C3$so_down_slope  <- (so_C3$mean_slope_to_trough_min_potential_per_second+so_C4$mean_slope_to_trough_min_potential_per_second)/2
so_C3$so_up_slope    <- (so_C3$mean_slope_trough_to_up_max_potential_per_second+so_C4$mean_slope_trough_to_up_max_potential_per_second)/2

so_params <- so_C3 %>% dplyr::select(datasetnum, so_count,so_dens, so_dur, so_amp, so_freq, so_act, so_down_slope, so_up_slope)
dataset_b <- merge(dataset_b, so_params, by = "datasetnum")

## Read in spindle-SO coupling data 
coupling_dataset_b  <- read.delim("dataset_b_coupling.csv", sep=",",dec= ".")
coupling_dataset_b  <- mutate(coupling_dataset_b, test_channel = fct_recode(test_channel,"EEGC3" = "C3-TP10", "EEGC4" = "C4-TP9"))
coupling_dataset_b  <- mutate(coupling_dataset_b, target_channel = fct_recode(target_channel,"EEGC3" = "C3-TP10","EEGC4" = "C4-TP9"))

# Split by channel to merge average below
coupling_C3        <- coupling_dataset_b[coupling_dataset_b$test_channel == "EEGC3", ]
coupling_C4        <- coupling_dataset_b[coupling_dataset_b$test_channel == "EEGC4", ]

## Count the amount of co-occurences
coupling_occur_c3  <- table(coupling_C3$event_files_number)
coupling_occur_c4  <- table(coupling_C4$event_files_number)

## Get sum of the amount of co-occurences
coupling_occur     <- (coupling_occur_c3+coupling_occur_c4)/2  
coupling_occur_sum <- coupling_occur_c3+coupling_occur_c4
coupling_occur     <- as.numeric(coupling_occur)
coupling_occur_sum <- as.numeric(coupling_occur_sum)

## Add to dataframe
dataset_b$coupling_occur     <- coupling_occur
dataset_b$coupling_occur_sum <- coupling_occur_sum

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
dataset_b$mean_delay     <- mean_delay$delay_c3
dataset_b$sd_delay       <- sd_delay$delay_c3

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
for (subnr in 1:118){ 
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
dataset_b$spin_dur_match    <- spin_dur_match
dataset_b$spin_dur_mismatch <- spin_dur_mismatch
dataset_b$spin_amp_match    <- spin_amp_match
dataset_b$spin_amp_mismatch <- spin_amp_mismatch
dataset_b$spin_freq_match   <- spin_freq_match
dataset_b$spin_freq_mismatch<- spin_freq_mismatch
dataset_b$spin_perc_match   <- spin_perc_match
dataset_b$so_dur_match      <- so_dur_match
dataset_b$so_dur_mismatch   <- so_dur_mismatch
dataset_b$so_amp_match      <- so_amp_match
dataset_b$so_amp_mismatch   <- so_amp_mismatch
dataset_b$so_perc_match     <- so_perc_match

dataset_b$spin_dur_diff  <- spin_dur_match - spin_dur_mismatch
dataset_b$spin_amp_diff  <- spin_amp_match - spin_amp_mismatch
dataset_b$spin_freq_diff <- spin_freq_match - spin_freq_mismatch
dataset_b$so_dur_diff    <- so_dur_match - so_dur_mismatch
dataset_b$so_amp_diff    <- so_amp_match - so_amp_mismatch

## IMPORTANT: Filter data ----
## We have data for 40 controls, 40 unmedicated and 38 medicated
## The analyses were initially run on all 118 datasets, but subsequent statistics will be performed with excluding the two matching unmedicated patients as well
## Missing are LP_104_2 and LP_107_2  //  remove -> LP_104_1 and LP_107_1

dataset_b <- dataset_b %>% dplyr::filter(dataset_b$id_loreta != "104" & dataset_b$id_loreta != "104_2" &
                                         dataset_b$id_loreta != "107" & dataset_b$id_loreta != "107_2")

dataset_b_coupling_dens_info <- merge(dataset_b, coupling_dens_info, by = "datasetnum")
dataset_b_coupling_dens_info$NREM_min <- (dataset_b_coupling_dens_info$S2_without_MA_min +
                                            dataset_b_coupling_dens_info$S3_without_MA_min) /2 
dataset_b_coupling_dens_info$NREM_min <- (dataset_b_coupling_dens_info$NREM_min +
                                            dataset_b_coupling_dens_info$S4_without_MA_min) /2 

dataset_b$coupling_dens<- dataset_b_coupling_dens_info$coupling_occur  /  dataset_b_coupling_dens_info$NREM_min

############################# Split datasets #############################

## Reorder factor levels
dataset_b$group <- factor(dataset_b$group, levels = c("control", "patient_unmed", "patient_7d"))

## Split into three groups:
## controls and patients unmed
dataset_b_conpat <- dataset_b %>% 
  dplyr::filter(group == "control" | group == "patient_unmed")
dataset_b_conpat$group <- droplevels(dataset_b_conpat$group)
dataset_b_conpat$group <- factor(dataset_b_conpat$group)

## controls and patients med
dataset_b_conpat_med <- dataset_b %>% 
  dplyr::filter(group == "control" | group == "patient_7d")
dataset_b_conpat_med$group <- droplevels(dataset_b_conpat_med$group)
dataset_b_conpat_med$group <- factor(dataset_b_conpat_med$group)

## patients 7 and 28
dataset_b_patients <- dataset_b %>% 
  dplyr::filter(group == "patient_unmed" | group == "patient_7d")
dataset_b_patients$group <- droplevels(dataset_b_patients$group)
dataset_b_patients$group <- factor(dataset_b_patients$group)

## Read in power data
# Merge with data
mergedData_long_dataset_b <- merge(dataset_b, power_dataset_b, by=c("datasetnum"))

## Reshape the data, per frequency
power_dataset_b_freq <- power_dataset_b %>%
  dplyr::select(datasetnum, channel, freq, mean_powerDensity_over_segments) %>%
  spread(key = freq, value = mean_powerDensity_over_segments)

## Average over the two channels
power_dataset_b_freq_reshape <- aggregate(. ~ datasetnum, data = power_dataset_b_freq, mean)

## Merge with data
mergedData_wide_dataset_b <- merge(dataset_b, power_dataset_b_freq_reshape, by = ("datasetnum"))

# Transform to decibel
mergedData_long_dataset_b$mean_densDB <- mergedData_long_dataset_b$mean_powerDensity_over_segments + 1
mergedData_long_dataset_b$mean_densDB <- 10*log10(mergedData_long_dataset_b$mean_densDB)

## Create subgroups:
## control patient wide
mergedData_wide_dataset_b_control_patient <- mergedData_wide_dataset_b %>%
  dplyr::filter(group == "control" | group == "patient_unmed")

mergedData_wide_dataset_b_control_patient$group <- droplevels(mergedData_wide_dataset_b_control_patient$group)
mergedData_wide_dataset_b_control_patient$group <- factor(mergedData_wide_dataset_b_control_patient$group)

## patient wide
mergedData_wide_dataset_b_patients <- mergedData_wide_dataset_b %>%
  dplyr::filter(group == "patient_unmed" | group == "patient_7d")

mergedData_wide_dataset_b_patients$group <- droplevels(mergedData_wide_dataset_b_patients$group)
mergedData_wide_dataset_b_patients$group <- factor(mergedData_wide_dataset_b_patients$group)

## control patient 28 wide
mergedData_wide_dataset_b_conpat_med <- mergedData_wide_dataset_b %>%
  dplyr::filter(group == "control" | group == "patient_7d")

mergedData_wide_dataset_b_conpat_med$group <- droplevels(mergedData_wide_dataset_b_conpat_med$group)
mergedData_wide_dataset_b_conpat_med$group <- factor(mergedData_wide_dataset_b_conpat_med$group)

## control patient long
mergedData_long_dataset_b_control_patient <- mergedData_long_dataset_b %>%
  dplyr::filter(group == "control" | group == "patient_unmed")

mergedData_long_dataset_b_control_patient$group <- droplevels(mergedData_long_dataset_b_control_patient$group)
mergedData_long_dataset_b_control_patient$group <- factor(mergedData_long_dataset_b_control_patient$group)

## patient longs
mergedData_long_dataset_b_patients <- mergedData_long_dataset_b %>%
  dplyr::filter(group == "patient_unmed" | group == "patient_7d")

mergedData_long_dataset_b_patients$group <- droplevels(mergedData_long_dataset_b_patients$group)
mergedData_long_dataset_b_patients$group <- factor(mergedData_long_dataset_b_patients$group)

##  control patient 28 wide
mergedData_long_dataset_b_conpat_med <- mergedData_long_dataset_b %>%
  dplyr::filter(group == "control" | group == "patient_7d")

mergedData_long_dataset_b_conpat_med$group <- droplevels(mergedData_long_dataset_b_conpat_med$group)
mergedData_long_dataset_b_conpat_med$group <- factor(mergedData_long_dataset_b_conpat_med$group)

mergedData_wide_dataset_b_patients_pair <- mergedData_wide_dataset_b_patients
mergedData_long_dataset_b_patients_pair <- mergedData_long_dataset_b_patients

## Calculate t.test for each of the frequencies between groups CONTROLS AND PATIENTS
column_dataset_b_06       <- which(colnames(mergedData_wide_dataset_b_control_patient)=="0.6")
column_dataset_b_pats_06  <- which(colnames(mergedData_wide_dataset_b_patients_pair)=="0.6")
column_dataset_b_conpat_med_06  <- which(colnames(mergedData_wide_dataset_b_conpat_med)=="0.6")

listtests   <- lapply(mergedData_wide_dataset_b_control_patient[,c(column_dataset_b_06:(column_dataset_b_06+147))], function(x) t.test(x ~ mergedData_wide_dataset_b_control_patient$group, var.equal = FALSE))
listtests_pats   <- lapply(mergedData_wide_dataset_b_patients_pair[,c(column_dataset_b_pats_06:(column_dataset_b_pats_06+147))], function(x) t.test(x ~ mergedData_wide_dataset_b_patients_pair$group, var.equal = TRUE, paired = TRUE))
listtests_conpat_med   <- lapply(mergedData_wide_dataset_b_conpat_med[,c(column_dataset_b_conpat_med_06:(column_dataset_b_conpat_med_06+147))], function(x) t.test(x ~ mergedData_wide_dataset_b_conpat_med$group, var.equal = F, paired = F))

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

data_long_pvalues <- merge(mergedData_long_dataset_b_control_patient, listpvalues, by.x=c("freq"), by.y=c("Freq"))
data_long_pvalues_pats <- merge(mergedData_long_dataset_b_patients_pair, listpvalues_pats, by.x=c("freq"), by.y=c("Freq"))
data_long_pvalues_conpat_med <- merge(mergedData_long_dataset_b_conpat_med, listpvalues_conpat_med, by.x=c("freq"), by.y=c("Freq"))

###### Permutation test -- patients controls
freq_columns_BASE   <- which(colnames(mergedData_wide_dataset_b_control_patient)=="0.6")
freq_columns <-  freq_columns_BASE:(freq_columns_BASE+147)

pb = txtProgressBar(min = 0, max = length(freq_columns), initial = 0)
perm_pvalue = rep(NA,length(freq_columns))
index = 1


for(i in freq_columns_BASE:(freq_columns_BASE+147)){
  ## Make new dataset
  perm_data <- mergedData_wide_dataset_b_control_patient %>% dplyr::select(i, group)
  
  # Split the data
  controls_freq <- perm_data[perm_data$group == "control",]
  patients_freq <- perm_data[perm_data$group == "patient_unmed",]
  
  # Make seperate lists
  combined_groups <- c(controls_freq$group, patients_freq$group)
  combined_freqs  <- c(controls_freq[,1], patients_freq[,1])
  
  # Number of simulations
  nsims <- 10
  
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

dataset_b_long_pvalues_perm <- merge(data_long_pvalues, perm_pvalue, by.x=c("freq"), by.y=c("Freq"))




###### Permutation test -- patients paired 
freq_columns_BASE   <- which(colnames(mergedData_wide_dataset_b_patients_pair)=="0.6")
freq_columns <-  freq_columns_BASE:(freq_columns_BASE+147)

pb = txtProgressBar(min = 0, max = length(freq_columns), initial = 0)
perm_pvalue_pats = rep(NA,length(freq_columns))
perm_pvalue_pats_ttest = rep(NA,length(freq_columns))
perm_pvalue_pats_coin = rep(NA,length(freq_columns))
perm_pvalue_pats_coin_paired = rep(NA,length(freq_columns))
mergedData_wide_dataset_b_patients_pair$id_loreta_match <- str_remove(mergedData_wide_dataset_b_patients_pair$id_loreta, "_2")
mergedData_wide_dataset_b_patients_pair$id_loreta_match <- as.factor(mergedData_wide_dataset_b_patients_pair$id_loreta_match)
index = 1

## Alternative: use coin package to check calculation by hand:: is the same! :-)
for(i in freq_columns_BASE:(freq_columns_BASE+147)){
  ## Make new dataset
  perm_data <- mergedData_wide_dataset_b_patients_pair %>% dplyr::select(i, group, id_loreta_match)
  perm_data$id_loreta_match <- as.factor(perm_data$id_loreta_match)
  
  ## Get pvalue from coin package function
  a <- independence_test(perm_data[[1]] ~ group, data = perm_data) 
  b <- symmetry_test(perm_data[[1]] ~ group | id_loreta_match, data = perm_data) 
  
  perm_pvalue_pats_coin[index] <- a@distribution@pvalue(a@statistic@teststatistic)
  perm_pvalue_pats_coin_paired[index] <- b@distribution@pvalue(b@statistic@teststatistic)
  index = index + 1
}
  
index = 1

for(i in freq_columns_BASE:(freq_columns_BASE+147)){
  ## Make new dataset
  perm_data <- mergedData_wide_dataset_b_patients_pair %>% dplyr::select(i, group)
  
  # Split the data
  controls_freq <- perm_data[perm_data$group == "patient_unmed",]
  patients_freq <- perm_data[perm_data$group == "patient_7d",]
  
  # Make seperate lists
  combined_groups <- c(controls_freq$group, patients_freq$group)
  combined_freqs  <- c(controls_freq[,1], patients_freq[,1])
  
  # Number of simulations
  nsims <- 10
  
  diff_obs <- mean(controls_freq[,1])- mean(patients_freq[,1])
  diffs_ttest_obs <- t.test(perm_data[[1]] ~ group, data = perm_data, paired = TRUE)$p.value
  
  diffs <- rep(NA, nsims)
  diffs_ttest_pvalues <- rep(NA, nsims)
  
  
  ## Permutation for loop
  for(k in 1:nsims){
    shuffled_labels <- sample(combined_groups, replace = FALSE)
    perm_data$group_shuffled <- shuffled_labels
    
    
    diffs[k] <- mean(combined_freqs[shuffled_labels == 1]) - mean(combined_freqs[shuffled_labels == 2])
    diffs_ttest_pvalues[k] <- t.test(perm_data[[1]] ~ group_shuffled, data = perm_data, paired = FALSE)$p.value
    
  }
  
  ## p-value = (number of more extreme differences than diff_obs/nsims)
  length(diffs[abs(diffs)>= abs(diff_obs)])/nsims
  perm_pvalue_pats[index] <- mean(abs(diffs) > abs(diff_obs))
  
  index = index + 1
  
  length(diffs_ttest_pvalues[abs(perm_pvalue_pats)<= abs(diff_obs)])/nsims
  perm_pvalue_pats_ttest[index] <- mean(abs(diffs_ttest_pvalues) < abs(diffs_ttest_obs))  
  
  setTxtProgressBar(pb,i)
  print
}

## Add them to the other dataframe
perm_pvalue_pats <- as.data.frame(perm_pvalue_pats)
perm_pvalue_pats_ttest      <- perm_pvalue_pats_ttest[-1]
perm_pvalue_pats_ttest <- as.data.frame(perm_pvalue_pats_ttest)


## Catch loop in cass a value is zero (will result in errors)
for (i in 1:nrow(perm_pvalue_pats)) {
  if ( perm_pvalue_pats[i,1] == 0.000){
    perm_pvalue_pats[i,1] <- 0.0001
  }
}

for (i in 1:nrow(perm_pvalue_pats_ttest)) {
  if ( perm_pvalue_pats_ttest[i,1] == 0.000){
    perm_pvalue_pats_ttest[i,1] <- 0.0001
  }
}

new_names <- seq(0.6, 30, by=0.2)
perm_pvalue_pats$Freq <- paste0(new_names)
perm_pvalue_pats$Freq<- as.numeric(perm_pvalue_pats$Freq)
perm_pvalue_pats_ttest$Freq <- paste0(new_names)
perm_pvalue_pats_ttest$Freq<- as.numeric(perm_pvalue_pats_ttest$Freq)

dataset_b_long_pvalues_perm_pats <- merge(data_long_pvalues_pats, perm_pvalue_pats, by.x=c("freq"), by.y=c("Freq"))
dataset_b_long_pvalues_perm_pats_ttest <- merge(data_long_pvalues_pats, perm_pvalue_pats_ttest, by.x=c("freq"), by.y=c("Freq"))


dataset_b_long_pvalues_perm_pats_10HZ <- dataset_b_long_pvalues_perm_pats %>% dplyr::filter(freq == 10)
dataset_b_long_pvalues_perm_pats_10HZ %>% group_by(group) %>% dplyr::summarise(mean = mean(mean_densDB, na.rm = T))


c <- cbind(perm_pvalue_pats$perm_pvalue_pats, perm_pvalue_pats_ttest$perm_pvalue_pats_ttest ,perm_pvalue_pats_coin, perm_pvalue_pats_coin_paired)


###### Permutation test --  controls patients 28d
freq_columns_BASE   <- which(colnames(mergedData_wide_dataset_b_conpat_med)=="0.6")
freq_columns <-  freq_columns_BASE:(freq_columns_BASE+147)

pb = txtProgressBar(min = 0, max = length(freq_columns), initial = 0)
perm_pvalue_conpat_med = rep(NA,length(freq_columns))
index = 1


for(i in freq_columns_BASE:(freq_columns_BASE+147)){
  ## Make new dataset
  perm_data <- mergedData_wide_dataset_b_conpat_med %>% dplyr::select(i, group)
  
  # Split the data
  controls_freq <- perm_data[perm_data$group == "control",]
  patients_freq <- perm_data[perm_data$group == "patient_7d",]
  
  # Make seperate lists
  combined_groups <- c(controls_freq$group, patients_freq$group)
  combined_freqs  <- c(controls_freq[,1], patients_freq[,1])
  
  # Number of simulations
  nsims <- 10
  
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

dataset_b_long_pvalues_perm_conpat_med <- merge(data_long_pvalues_conpat_med, perm_pvalue_conpat_med, by.x=c("freq"), by.y=c("Freq"))


############################# Part II: Summary statistics and ttests #############################

## Demographics ----
dataset_b_demo_long <- dataset_b %>%
  gather(demo, score, age, gender, hamd_0, no_episodes)

dataset_b %>%
  group_by(group) %>%
  summarise_at(vars(age, gender, hamd_0, no_episodes), funs(mean(., na.rm=TRUE)))
dataset_b %>%
  group_by(group) %>%
  summarise_at(vars(age, gender, hamd_0, no_episodes), funs(std.error(., na.rm=TRUE)))

table(dataset_b$group, dataset_b$gender)

## Sleep architecture control - unmed ---- 
dataset_b_conpat_hypno_long <- dataset_b_conpat %>%
  gather(sleepstage, score, S1_percent,S2_percent,S3andS4_perc,REM_percent,NREM_percent,Wake_after_sleep_onset_percent,Total_sleep_time_min, Sleep_Onset_min, SWS_onset_min, REM_onset_min)

dataset_b_conpat_hypno_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

hypnogram_ttests_dataset_b_conpat <- lapply(dataset_b_conpat[,c("age", "S1_percent", "S2_percent", "S3andS4_perc",  "REM_percent", "NREM_percent", "Wake_after_sleep_onset_percent", "Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min" )], function(x) t.test(x ~ dataset_b_conpat$group, var.equal = FALSE))
  
ttestBF(formula = S1_percent ~ group, data = dataset_b_conpat)
ttestBF(formula = S2_percent ~ group, data = dataset_b_conpat)
ttestBF(formula = S3andS4_perc ~ group, data = dataset_b_conpat)
ttestBF(formula = NREM_percent ~ group, data = dataset_b_conpat)
ttestBF(formula = REM_percent ~ group, data = dataset_b_conpat)
ttestBF(formula = Wake_after_sleep_onset_percent ~ group, data = dataset_b_conpat)
ttestBF(formula = Total_sleep_time_min ~ group, data = dataset_b_conpat)
ttestBF(formula = Sleep_Onset_min ~ group, data = dataset_b_conpat)
ttestBF(formula = SWS_onset_min ~ group, data = dataset_b_conpat)
ttestBF(formula = REM_onset_min ~ group, data = dataset_b_conpat)

## Sleep architecture control - med ---- 
dataset_b_conpat_med_hypno_long <- dataset_b_conpat_med %>%
  gather(sleepstage, score, S1_percent,S2_percent,S3andS4_perc,REM_percent,NREM_percent,Wake_after_sleep_onset_percent,Total_sleep_time_min, Sleep_Onset_min, SWS_onset_min, REM_onset_min)

dataset_b_conpat_med_hypno_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

hypnogram_ttests_dataset_b_conpat_med <- lapply(dataset_b_conpat_med[,c("age", "S1_percent", "S2_percent", "S3andS4_perc",  "REM_percent", "NREM_percent", "Wake_after_sleep_onset_percent", "Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min" )], function(x) t.test(x ~ dataset_b_conpat_med$group, var.equal = FALSE))

ttestBF(formula = S1_percent ~ group, data = dataset_b_conpat_med)
ttestBF(formula = S2_percent ~ group, data = dataset_b_conpat_med)
ttestBF(formula = S3andS4_perc ~ group, data = dataset_b_conpat_med)
ttestBF(formula = NREM_percent ~ group, data = dataset_b_conpat_med)
ttestBF(formula = REM_percent ~ group, data = dataset_b_conpat_med)
ttestBF(formula = Wake_after_sleep_onset_percent ~ group, data = dataset_b_conpat_med)
ttestBF(formula = Total_sleep_time_min ~ group, data = dataset_b_conpat_med)
ttestBF(formula = Sleep_Onset_min ~ group, data = dataset_b_conpat_med)
ttestBF(formula = REM_onset_min ~ group, data = dataset_b_conpat_med)

dataset_b_conpat_med_SWS <- dataset_b_conpat_med %>% drop_na(SWS_onset_min) ## Subject 113_2
ttestBF(formula = SWS_onset_min ~ group, data = dataset_b_conpat_med_SWS)

## Sleep architecture unmed - med ---- 
dataset_b_patients_hypno_long <- dataset_b_patients %>%
  gather(sleepstage, score, S1_percent,S2_percent,S3andS4_perc,REM_percent,NREM_percent,Wake_after_sleep_onset_percent,Total_sleep_time_min, Sleep_Onset_min, SWS_onset_min, REM_onset_min)

dataset_b_patients_hypno_long %>%
  dplyr::select(sleepstage, score, group) %>%
  group_by(group, sleepstage) %>%
  get_summary_stats(score, type = "mean_se")

hypnogram_ttests_dataset_b_patients <- lapply(dataset_b_patients[,c("age", "S1_percent", "S2_percent", "S3andS4_perc",  "REM_percent", "NREM_percent", "Wake_after_sleep_onset_percent", "Total_sleep_time_min", "Sleep_Onset_min", "REM_onset_min")], function(x) t.test(x ~ dataset_b_patients$group, paired = T))

ttestBF(formula = S1_percent ~ group, data = dataset_b_patients)
ttestBF(formula = S2_percent ~ group, data = dataset_b_patients)
ttestBF(formula = S3andS4_perc ~ group, data = dataset_b_patients)
ttestBF(formula = NREM_percent ~ group, data = dataset_b_patients)
ttestBF(formula = REM_percent ~ group, data = dataset_b_patients)
ttestBF(formula = Wake_after_sleep_onset_percent ~ group, data = dataset_b_patients)
ttestBF(formula = Total_sleep_time_min ~ group, data = dataset_b_patients)
ttestBF(formula = Sleep_Onset_min ~ group, data = dataset_b_patients)
ttestBF(formula = REM_onset_min ~ group, data = dataset_b_patients)

## since there are NAs in SWS: remove extra to get full pairs
dataset_b_patients_SWS_onset <- dataset_b_patients %>% drop_na(SWS_onset_min) ## Subject 113_2
dataset_b_patients_SWS_onset <- dataset_b_patients_SWS_onset %>% dplyr::filter(dataset_b_patients_SWS_onset$id_loreta != "113" )

ttests_dataset_b_patients_SWS_onset <- t.test(dataset_b_patients_SWS_onset$SWS_onset_min ~dataset_b_patients_SWS_onset$group, paired = T)
ttestBF(formula = SWS_onset_min ~ group, data = dataset_b_patients_SWS_onset)

dataset_b_pat_med <- dataset_b_patients %>% filter(group ==  "patient_7d")
cor.test(dataset_b_pat_med$S1_percent, dataset_b_pat_med$REM_percent)


## Spindle data control - unmed ----  
dataset_b_conpat_spindle_long <- dataset_b_conpat %>%
  gather(spindle, score, spin_dens, spin_count,spin_amp, spin_dur, spin_freq)

dataset_b_conpat_spindle_long %>%
  dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%
  get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_b_conpat <- lapply(dataset_b_conpat[,c("spin_dens","spin_dur","spin_amp", "spin_count", "spin_centerfreq", "spin_freq"  )], function(x) t.test(x ~ dataset_b_conpat$group, var.equal = FALSE))

## Spindle data unmed -  med ----  
dataset_b_patients_spindle_long <- dataset_b_patients %>%
  gather(spindle, score, spin_dens, spin_count,spin_amp, spin_dur,spin_freq)

dataset_b_patients_spindle_long %>%
  dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%
  get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_b_patients <- lapply(dataset_b_patients[,c("spin_dens","spin_dur","spin_amp", "spin_count" , "spin_centerfreq", "spin_freq" )], function(x) t.test(x ~ dataset_b_patients$group, paired = T))

## Spindle data control -  med ----  
dataset_b_conpat_med_spindle_long <- dataset_b_conpat_med %>%
  gather(spindle, score, spin_dens, spin_count,spin_amp, spin_dur, spin_freq)

dataset_b_conpat_med_spindle_long %>%
  dplyr::select(spindle, score, group) %>%
  group_by(group, spindle) %>%
  get_summary_stats(score, type = "mean_se")

spindle_ttests_dataset_b_conpat_med <- lapply(dataset_b_conpat_med[,c("spin_dens","spin_dur","spin_amp", "spin_count", "spin_centerfreq" , "spin_freq" )], function(x) t.test(x ~ dataset_b_conpat_med$group, paired = F))

## Calculate bayes factor
ttestBF(formula = spin_dens ~ group, data = dataset_b_conpat)
ttestBF(formula = spin_count ~ group, data = dataset_b_conpat)
ttestBF(formula = spin_amp ~ group, data = dataset_b_conpat)
ttestBF(formula = spin_freq ~ group, data = dataset_b_conpat)
ttestBF(formula = spin_dur ~ group, data = dataset_b_conpat)
ttestBF(formula = spin_dens ~ group, data = dataset_b_patients)
ttestBF(formula = spin_count ~ group, data = dataset_b_patients)
ttestBF(formula = spin_amp ~ group, data = dataset_b_patients)
ttestBF(formula = spin_freq ~ group, data = dataset_b_patients)
ttestBF(formula = spin_dur ~ group, data = dataset_b_patients)
ttestBF(formula = spin_dens ~ group, data = dataset_b_conpat_med)
ttestBF(formula = spin_count ~ group, data = dataset_b_conpat_med)
ttestBF(formula = spin_amp ~ group, data = dataset_b_conpat_med)
ttestBF(formula = spin_freq ~ group, data = dataset_b_conpat_med)
ttestBF(formula = spin_dur ~ group, data = dataset_b_conpat_med)

## SO data control - unmed ---- 
dataset_b_conpat_so_long <- dataset_b_conpat %>%
  gather(so, score, so_dens,so_amp, so_dur, so_count,so_freq)

dataset_b_conpat_so_long %>%
  dplyr::select(so, score, group) %>%
  group_by(group, so) %>%
  get_summary_stats(score, type = "mean_se")

so_ttests_dataset_b_conpat <- lapply(dataset_b_conpat[,c("so_dens","so_dur","so_amp", "so_count", "so_freq" )], function(x) t.test(x ~ dataset_b_conpat$group, var.equal = FALSE))

## SO data unmed - med ---- 
dataset_b_patients_so_long <- dataset_b_patients %>%
  gather(so, score, so_dens,so_amp, so_dur, so_count,so_freq)

dataset_b_patients_so_long %>%
   dplyr::select(so, score, group) %>%
   group_by(group, so) %>%
   get_summary_stats(score, type = "mean_se")

so_ttests_dataset_b_patients <- lapply(dataset_b_patients[,c("so_dens","so_dur","so_amp", "so_count" , "so_freq" )], function(x) t.test(x ~ dataset_b_patients$group, paired = T))

## SO data control - med ---- 
dataset_b_conpat_med_so_long <- dataset_b_conpat_med %>%
  gather(so, score, so_dens,so_amp, so_dur,so_count,so_freq)

dataset_b_conpat_med_so_long %>%
  dplyr::select(so, score, group) %>%
  group_by(group, so) %>%
  get_summary_stats(score, type = "mean_se")

so_ttests_dataset_b_conpat_med <- lapply(dataset_b_conpat_med[,c("so_dens","so_dur","so_amp", "so_count" , "so_freq")], function(x) t.test(x ~ dataset_b_conpat_med$group, paired = F))

## Calculate bayes factor
ttestBF(formula = so_dens ~ group, data = dataset_b_conpat)
ttestBF(formula = so_count ~ group, data = dataset_b_conpat)
ttestBF(formula = so_amp ~ group, data = dataset_b_conpat)
ttestBF(formula = so_freq ~ group, data = dataset_b_conpat)
ttestBF(formula = so_dur ~ group, data = dataset_b_conpat)
ttestBF(formula = so_dens ~ group, data = dataset_b_patients)
ttestBF(formula = so_count ~ group, data = dataset_b_patients)
ttestBF(formula = so_amp ~ group, data = dataset_b_patients)
ttestBF(formula = so_freq ~ group, data = dataset_b_patients)
ttestBF(formula = so_dur ~ group, data = dataset_b_patients)
ttestBF(formula = so_dens ~ group, data = dataset_b_conpat_med)
ttestBF(formula = so_count ~ group, data = dataset_b_conpat_med)
ttestBF(formula = so_amp ~ group, data = dataset_b_conpat_med)
ttestBF(formula = so_freq ~ group, data = dataset_b_conpat_med)
ttestBF(formula = so_dur ~ group, data = dataset_b_conpat_med)


## Coupling data  control - unmed ---- 
dataset_b_conpat_coupling_long <- dataset_b_conpat %>%
  gather(coupling, score, coupling_occur,mean_delay,sd_delay)

dataset_b_conpat_coupling_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_ttests_dataset_b_conpat <- lapply(dataset_b_conpat[,c("coupling_occur","mean_delay","sd_delay", "coupling_dens" )], function(x) t.test(x ~ dataset_b_conpat$group, var.equal = FALSE))

## Coupling data unmed - med ---- 
dataset_b_patients_coupling_long <- dataset_b_patients %>%
  gather(coupling, score, coupling_occur,mean_delay,sd_delay)

dataset_b_patients_coupling_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_ttests_dataset_b_patients <- lapply(dataset_b_patients[,c("coupling_occur","mean_delay","sd_delay","coupling_dens"  )], function(x) t.test(x ~ dataset_b_patients$group, paired = T))

## Coupling data  control - med ---- 
dataset_b_conpat_med_coupling_long <- dataset_b_conpat_med %>%
  gather(coupling, score, coupling_occur,mean_delay,sd_delay)

dataset_b_conpat_med_coupling_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

coupling_ttests_dataset_b_conpat_med <- lapply(dataset_b_conpat_med[,c("coupling_occur","mean_delay","sd_delay", "coupling_dens")], function(x) t.test(x ~ dataset_b_conpat_med$group, paired = F))

ttestBF(formula = coupling_occur ~ group, data = dataset_b_conpat)
ttestBF(formula = mean_delay ~ group, data = dataset_b_conpat)
ttestBF(formula = sd_delay ~ group, data = dataset_b_conpat)
ttestBF(formula = coupling_occur ~ group, data = dataset_b_patients)
ttestBF(formula = mean_delay ~ group, data = dataset_b_patients)
ttestBF(formula = sd_delay ~ group, data = dataset_b_patients)
ttestBF(formula = coupling_occur ~ group, data = dataset_b_conpat_med)
ttestBF(formula = mean_delay ~ group, data = dataset_b_conpat_med)
ttestBF(formula = sd_delay ~ group, data = dataset_b_conpat_med)


## Coupling matching vs nonmatching
dataset_b_conpat_med_coupling_matchspin_long <- dataset_b_conpat_med  %>%
  gather(coupling, score, spin_amp_match, spin_dur_match,spin_freq_match,spin_perc_match,
         so_amp_match, so_dur_match,so_perc_match)

dataset_b_conpat_med_coupling_matchspin_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

dataset_b_conpat_coupling_matchspin_long <- dataset_b_conpat  %>%
  gather(coupling, score, spin_amp_match, spin_dur_match,spin_freq_match,spin_perc_match,
         so_amp_match, so_dur_match,so_perc_match)

dataset_b_conpat_coupling_matchspin_long %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

dataset_b_conpat_med %>%
    gather(match, score, spin_freq_match,spin_freq_mismatch) %>%
    dplyr::select(match, score, group) %>%
    summary(lm(score~match * group))

dataset_b_conpat%>%
  gather(match, score, spin_amp_diff, spin_dur_diff,spin_freq_diff, so_amp_diff, so_dur_diff) %>%
  dplyr::select(match, score, group) %>%
  group_by(group, match) %>% get_summary_stats(score, type = "mean_se")

dataset_b_patients %>%
  gather(match, score, spin_amp_diff, spin_dur_diff,spin_freq_diff, so_amp_diff, so_dur_diff) %>%
  dplyr::select(match, score, group) %>%
  group_by(group, match) %>% get_summary_stats(score, type = "mean_se")

dataset_b_conpat_med %>%
  gather(match, score, spin_amp_diff, spin_dur_diff,spin_freq_diff, so_amp_diff, so_dur_diff) %>%
  dplyr::select(match, score, group) %>%
  group_by(group, match) %>% get_summary_stats(score, type = "mean_se")

coupling_diff_ttests_dataset_b_conpat <- lapply(dataset_b_conpat[,c("spin_amp_diff","spin_dur_diff","spin_freq_diff","spin_perc_match","so_amp_diff", "so_dur_diff", "so_perc_match" )], function(x) t.test(x ~ dataset_b_conpat$group, var.equal = FALSE))
coupling_diff_ttests_dataset_b_patients <- lapply(dataset_b_patients[,c("spin_amp_diff","spin_dur_diff","spin_freq_diff","spin_perc_match","so_amp_diff", "so_dur_diff", "so_perc_match" )], function(x) t.test(x ~ dataset_b_patients$group, var.equal = FALSE, paired = T))
coupling_diff_ttests_dataset_b_conpat_med <- lapply(dataset_b_conpat_med[,c("spin_amp_diff","spin_dur_diff","spin_freq_diff","spin_perc_match","so_amp_diff", "so_dur_diff", "so_perc_match" )], function(x) t.test(x ~ dataset_b_conpat_med$group, var.equal = FALSE))

coupling_match_ttests_dataset_b_conpat <- lapply(dataset_b_conpat[,c("spin_amp_match","spin_dur_match","spin_freq_match","so_amp_match", "so_dur_match")], function(x) t.test(x ~ dataset_b_conpat$group, var.equal = FALSE))
coupling_match_ttests_dataset_b_patients <- lapply(dataset_b_patients[,c("spin_amp_match","spin_dur_match","spin_freq_match","so_amp_match", "so_dur_match")], function(x) t.test(x ~ dataset_b_patients$group, var.equal = FALSE, paired = T))
coupling_match_ttests_dataset_b_conpat_med <- lapply(dataset_b_conpat_med[,c("spin_amp_match","spin_dur_match","spin_freq_match","so_amp_match", "so_dur_match" )], function(x) t.test(x ~ dataset_b_conpat_med$group, var.equal = FALSE))

dataset_b_conpat_med_long_spin_diff <- dataset_b_conpat_med %>%
  gather(coupling, score, spin_freq_diff, spin_dur_diff)

dataset_b_conpat_med_long_spin_diff %>%
  dplyr::select(coupling, score, group) %>%
  group_by(group, coupling) %>%
  get_summary_stats(score, type = "mean_se")

ttestBF(formula = spin_amp_match ~ group, data  = dataset_b_conpat)
ttestBF(formula = spin_dur_match ~ group, data  = dataset_b_conpat)
ttestBF(formula = spin_freq_match ~ group, data = dataset_b_conpat)
ttestBF(formula = so_amp_match ~ group, data    = dataset_b_conpat)
ttestBF(formula = so_dur_match ~ group, data    = dataset_b_conpat)
ttestBF(formula = spin_amp_match ~ group, data  = dataset_b_patients)
ttestBF(formula = spin_dur_match ~ group, data  = dataset_b_patients)
ttestBF(formula = spin_freq_match ~ group, data = dataset_b_patients)
ttestBF(formula = so_amp_match ~ group, data    = dataset_b_patients)
ttestBF(formula = so_dur_match ~ group, data    = dataset_b_patients)
ttestBF(formula = spin_amp_match ~ group, data  = dataset_b_conpat_med)
ttestBF(formula = spin_dur_match ~ group, data  = dataset_b_conpat_med)
ttestBF(formula = spin_freq_match ~ group, data = dataset_b_conpat_med)
ttestBF(formula = so_amp_match ~ group, data    = dataset_b_conpat_med)
ttestBF(formula = so_dur_match ~ group, data    = dataset_b_conpat_med)

ttestBF(formula = spin_amp_diff ~ group, data  = dataset_b_conpat)
ttestBF(formula = spin_dur_diff ~ group, data  = dataset_b_conpat)
ttestBF(formula = spin_freq_diff ~ group, data = dataset_b_conpat)
ttestBF(formula = so_amp_diff ~ group, data    = dataset_b_conpat)
ttestBF(formula = so_dur_diff ~ group, data    = dataset_b_conpat)
ttestBF(formula = spin_amp_diff ~ group, data  = dataset_b_patients)
ttestBF(formula = spin_dur_diff ~ group, data  = dataset_b_patients)
ttestBF(formula = spin_freq_diff ~ group, data = dataset_b_patients)
ttestBF(formula = so_amp_diff ~ group, data    = dataset_b_patients)
ttestBF(formula = so_dur_diff ~ group, data    = dataset_b_patients)
ttestBF(formula = spin_amp_diff ~ group, data  = dataset_b_conpat_med)
ttestBF(formula = spin_dur_diff ~ group, data  = dataset_b_conpat_med)
ttestBF(formula = spin_freq_diff ~ group, data = dataset_b_conpat_med)
ttestBF(formula = so_amp_diff ~ group, data    = dataset_b_conpat_med)
ttestBF(formula = so_dur_diff ~ group, data    = dataset_b_conpat_med)


## Part III: Plots ---- 
## Combine data for plots
dataset_b_long_pvalues_perm_select <- dataset_b_long_pvalues_perm  %>% dplyr::select(pnum, id_loreta, datasetnum, group, freq, mean_densDB, pvalue)
dataset_b_long_pvalues_perm_pats_select   <- dataset_b_long_pvalues_perm_pats   %>% dplyr::select(pnum, id_loreta, datasetnum,  group, freq, mean_densDB, pvalue)
dataset_b_long_pvalues_perm_conpat_med_select   <- dataset_b_long_pvalues_perm_conpat_med  %>% dplyr::select(pnum, id_loreta, datasetnum, group, freq, mean_densDB, pvalue)

dataset_b_long_pvalues_perm_control <- dataset_b_long_pvalues_perm_select %>% filter(group == "control")
dataset_b_long_pvalues_perm_pat_unmed    <- dataset_b_long_pvalues_perm_select %>% filter(group == "patient_unmed")
dataset_b_long_pvalues_perm_pat_med   <- dataset_b_long_pvalues_perm_pats_select %>% filter(group == "patient_7d")

dataset_b_long_pvalues_perm_three <- rbind(dataset_b_long_pvalues_perm_control, dataset_b_long_pvalues_perm_pat_unmed, dataset_b_long_pvalues_perm_pat_med)

## Power plot control - unmed
freq_spec_dataset_b_conpat <- dataset_b_long_pvalues_perm %>%
  ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
  theme_classic()+ my_theme +
  stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
  scale_color_manual(values=color_dataset_b, labels = c("Controls (B)", "Medicated (B)")) +  
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

perm_plot_dataset_b_conpat <- ggplot(dataset_b_long_pvalues_perm, aes(x=freq, y=abs(log10(perm_pvalue)), colour="#000000")) +
  theme_classic()+ my_theme +
  stat_summary(fun.y="mean", geom="bar", colour="black")+ 
  xlab("\nFrequency") + ylab("P-values (con - unmed)") +
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

## Power plot unmed - med
freq_spec_dataset_b_pats <- dataset_b_long_pvalues_perm_pats_ttest %>%
  ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
  theme_classic()+my_theme +
  stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
  scale_color_manual(values=color_dataset_b_pat, labels = c("Medicated - 7days (B)", "Medicated - 28days (B)")) +  
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

perm_plot_dataset_b_pats <- ggplot(dataset_b_long_pvalues_perm_pats_ttest, aes(x=freq, y=abs(log10(perm_pvalue_pats_ttest)), colour="#000000")) +
  theme_classic()+ my_theme +
  stat_summary(fun.y="mean", geom="bar", colour="black")+ 
  xlab("\nFrequency") + ylab("P-values (unmed - med)") +
  theme(       
    axis.text.y = element_text( color="#000000"),
    axis.title.x = element_blank(),
    axis.line.x  = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(), 
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 4),labels=c("0",  "10e-1","10e-2", "10e-3", "10e-4"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 
  NULL

## Power plot control - unmed - med
freq_spec_dataset_b_three <- dataset_b_long_pvalues_perm_three %>%
  ggplot(aes(x = freq, y = mean_densDB, colour = group)) + 
  theme_classic()+my_theme +
  stat_summary(fun.y="mean", geom = 'smooth', se = FALSE)+              
  scale_color_manual(values=color_dataset_b_all, labels = c("Controls - 7days (B)", "Medicated - 7days (B)", "Medicated - 28days (B)")) +  
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

perm_plot_dataset_b_conpat_med <- ggplot(dataset_b_long_pvalues_perm_conpat_med, aes(x=freq, y=abs(log10(perm_pvalue_conpat_med)), colour="#000000")) +
  theme_classic()+ my_theme +
  stat_summary(fun.y="mean", geom="bar", colour="black")+ 
  xlab("\nFrequency") + ylab("P-values (con - med)") +
  theme(       
    axis.text.y = element_text( color="#000000"),
    axis.text.x = element_text(color="#000000"), 
    legend.position='none')+ 
  geom_hline(yintercept = abs(log10(0.05)))+
  scale_y_continuous(limits = c(0, 3),labels=c("0",  "10e-1","10e-2", "10e-3"))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) 

## Spindles control - unmed - med
plot_rain_spidens_dataset_b <- ggplot(dataset_b,aes(x=group,y=spin_dens,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_dens),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle density')+ xlab('')+ ylim(0,4) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values= color_dataset_b_all) + 
  scale_x_discrete(labels = c("Controls (B)", "Unmedicated (B)", "Medicated (B)")) +
  scale_fill_manual(values=color_dataset_b_all) + 
  ggtitle("Group difference in spindle density") 

plot_rain_spinamp_dataset_b <- ggplot(dataset_b,aes(x=group,y=spin_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=spin_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Spindle amplitude')+ xlab('')+ ylim(0,55) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+ 
  scale_color_manual(values=color_dataset_b_all) + 
  scale_fill_manual(values=color_dataset_b_all) + 
  scale_x_discrete(labels = c("Controls (B)", "Unmedicated (B)", "Medicated (B)")) +
  ggtitle("Group difference in spindle amplitude") +
  geom_signif(annotation=formatC("*", digits=1), y_position=42,xmin=1.2, xmax=2.2, 
              textsize = 6, colour="BLACK") +
  geom_signif(annotation=formatC("*", digits=1), y_position=47,xmin=1.2, xmax=3.2, 
              textsize = 6, colour="BLACK") 

## SO control - unmed - med
plot_rain_soamp_dataset_b <- ggplot(dataset_b,aes(x=group,y=so_amp,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_amp),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO amplitude')+ xlab('')+ ylim(0,300) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+
  scale_color_manual(values=color_dataset_b_all) + 
  scale_fill_manual(values=color_dataset_b_all) + 
  scale_x_discrete(labels = c("Controls (B)", "Unmedicated (B)", "Medicated (B)")) +
  ggtitle("Group difference in SO amplitude")  

plot_rain_sodur_dataset_b <- ggplot(dataset_b,aes(x=group,y=so_dur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=so_dur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO duration')+ xlab('')+ ylim(1,1.6) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+
  scale_color_manual(values=color_dataset_b_all) + 
  scale_fill_manual(values=color_dataset_b_all) + 
  scale_x_discrete(labels = c("Controls (B)", "Unmedicated (B)", "Medicated (B)")) +
  ggtitle("Group difference in SO duration") +
  geom_signif(annotation=formatC("*", digits=1), y_position=1.55,xmin=2.2, xmax=3.2,
              textsize = 6, colour="BLACK")

## Coupling control - unmed - med
plot_rain_amount_coup_dataset_b <- ggplot(dataset_b,aes(x=group,y=coupling_occur,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=coupling_occur),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('SO-fast spindles')+ xlab('')+ ylim(0,300) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+
  scale_color_manual(values=color_dataset_b_all) + 
  scale_fill_manual(values=color_dataset_b_all) + 
  scale_x_discrete(labels = c("Controls (B)", "Unmedicated (B)", "Medicated (B)")) +  ggtitle("Group difference in amount of SO-fast spindle couplings")

plot_rain_mean_delay_dataset_b <- ggplot(dataset_b,aes(x=group,y=mean_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=mean_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Mean delay')+ xlab('')+ ylim(0.2,0.85) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+
  scale_color_manual(values=color_dataset_b_all) + 
  scale_fill_manual(values=color_dataset_b_all) + 
  scale_x_discrete(labels = c("Controls (B)", "Unmedicated (B)", "Medicated (B)")) +  ggtitle("Group difference in mean delay")

plot_rain_sd_delay_dataset_b <- ggplot(dataset_b,aes(x=group,y=sd_delay,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=2.5, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=sd_delay),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Delay dispersion (SD)')+ xlab('')+ ylim(0.15,0.45) +
  raincloud_theme + my_theme +
  guides(fill=FALSE,colour=FALSE)+
  scale_color_manual(values=color_dataset_b_all) + 
  scale_fill_manual(values=color_dataset_b_all) + 
  scale_x_discrete(labels = c("Controls (B)", "Unmedicated (B)", "Medicated (B)")) +  ggtitle("Group difference in delay dispersion (SD)") +
  geom_signif(annotation=formatC("*", digits=1), y_position=0.4,xmin=1.2, xmax=3.2,
              textsize = 6, colour="BLACK")


## Add the plots on top of eachother
freq_plot_dataset_b       <- ggarrange(freq_spec_dataset_b_conpat, perm_plot_dataset_b_conpat, heights = c(2, 1), ncol = 1, nrow = 2)
freq_plot_pats_dataset_b  <- ggarrange(freq_spec_dataset_b_pats, perm_plot_dataset_b_pats, heights = c(2, 1), ncol = 1, nrow = 2)
freq_plot_three_dataset_b  <- ggarrange(freq_spec_dataset_b_three, perm_plot_dataset_b_conpat, perm_plot_dataset_b_pats, perm_plot_dataset_b_conpat_med, heights = c(3, 1, 1, 1.5), ncol = 1, nrow = 4)

spin_plot_dataset_b     <- ggarrange(plot_rain_spidens_dataset_b, plot_rain_spinamp_dataset_b,nrow = 2, labels = c("A", "B") )
so_plot_dataset_b       <- ggarrange(plot_rain_soamp_dataset_b, plot_rain_sodur_dataset_b,nrow = 2, labels = c("A", "B") )
coupling_plot_dataset_b <- ggarrange(plot_rain_amount_coup_dataset_b, plot_rain_mean_delay_dataset_b, plot_rain_sd_delay_dataset_b, nrow = 3, labels = c("A", "B", "C") )

## Save the images
# ggsave("figs/freq_plot_dataset_b.png", freq_plot_dataset_b, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/freq_plot_pats_dataset_b.png", freq_plot_pats_dataset_b, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/freq_plot_three_dataset_b.png", freq_plot_three_dataset_b, device = "png", dpi = 300, width = 7, height = 15)
# ggsave("figs/spin_plot_dataset_b.png", spin_plot_dataset_b, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/so_plot_dataset_b.png", so_plot_dataset_b, device = "png", dpi = 300, width = 7, height = 10)
# ggsave("figs/coupling_plot_dataset_b.png", coupling_plot_dataset_b, device = "png", dpi = 300, width = 7, height = 10)

## Part IV: Additional analyses ----
## Hippocampal size ----
## Merge with data from sourced function
dataset_b_conpat_hipp <- merge(dataset_b_conpat, hipp_all, by.x=c("id_loreta"), by.y=c("loreta_id"), all.x = T)

## Group differences in size left/right
t.test(dataset_b_conpat_hipp$hipp_size_left  ~ dataset_b_conpat_hipp$group)
t.test(dataset_b_conpat_hipp$hipp_size_right ~ dataset_b_conpat_hipp$group)

## Combine
hipp_all_long <- dataset_b_conpat_hipp %>% gather(hemisphere, size, hipp_size_left, hipp_size_right)

Anova(lm(size ~ group * hemisphere, hipp_all_long))

group_by(hipp_all_long, group, hemisphere) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(size, na.rm = TRUE),
    sd = sd(size, na.rm = TRUE)
  )



##  Finger tapping ----
## Read in data
tapping  <- read_excel("tapping_info.xlsx")
tapping  <- subset(tapping, select=c("BASE", "END", "TEST", "TRAIN", "CONS",  "id_loreta", "datasetnum")) 

dataset_b <- merge(dataset_b, tapping, by = c("datasetnum", "id_loreta"), all.x = T )
cols<- c("BASE", "END", "TEST", "TRAIN", "CONS")
dataset_b[cols] <- sapply(dataset_b[cols],as.numeric)

dataset_b$TRAIN_change <- ((dataset_b$END -  dataset_b$BASE)/ dataset_b$BASE)
dataset_b$CONS_change  <- ((dataset_b$TEST - dataset_b$END)  / dataset_b$END)
is.na(dataset_b) <- sapply(dataset_b, is.infinite)

## Split into three groups:
## controls and patients unmed
dataset_b_conpat <- dataset_b %>% 
  dplyr::filter(group == "control" | group == "patient_unmed")
dataset_b_conpat$group <- droplevels(dataset_b_conpat$group)
dataset_b_conpat$group <- factor(dataset_b_conpat$group)

## controls and patients med
dataset_b_conpat_med <- dataset_b %>% 
  dplyr::filter(group == "control" | group == "patient_7d")
dataset_b_conpat_med$group <- droplevels(dataset_b_conpat_med$group)
dataset_b_conpat_med$group <- factor(dataset_b_conpat_med$group)

## patients 7 and 28
dataset_b_patients <- dataset_b %>% 
  dplyr::filter(group == "patient_unmed" | group == "patient_7d")
dataset_b_patients$group <- droplevels(dataset_b_patients$group)
dataset_b_patients$group <- factor(dataset_b_patients$group)

## Get summary statistics 
dataset_b_long <- dataset_b %>% gather(memory, score, BASE,TRAIN_change, CONS_change)
dataset_b_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")

## T- test Group differences in behaviour
beh_BASE_ttest_dataset_b  <- t.test(BASE ~ group, data = dataset_b_conpat)
beh_train_ttest_dataset_b <- t.test(TRAIN_change ~ group, data = dataset_b_conpat)
beh_cons_ttest_dataset_b  <- t.test(CONS_change  ~ group, data = dataset_b_conpat)  
beh_BASE_ttest_dataset_b_conpat_med  <- t.test(BASE ~ group, data = dataset_b_conpat_med)
beh_train_ttest_dataset_b_conpat_med <- t.test(TRAIN_change ~ group, data = dataset_b_conpat_med)
beh_cons_ttest_dataset_b_conpat_med  <- t.test(CONS_change  ~ group, data = dataset_b_conpat_med)  
beh_BASE_ttest_dataset_b_patients  <- t.test(BASE ~ group, data = dataset_b_patients)
beh_train_ttest_dataset_b_patients <- t.test(TRAIN_change ~ group, data = dataset_b_patients)
beh_cons_ttest_dataset_b_patients  <- t.test(CONS_change  ~ group, data = dataset_b_patients) 

dataset_b_conpat_base_na   <- dataset_b_conpat %>% drop_na(BASE)
dataset_b_conpat_train_na  <- dataset_b_conpat %>% drop_na(TRAIN_change)
dataset_b_conpat_cons_na   <- dataset_b_conpat %>% drop_na(CONS_change)

## bayes statistics
bayes_BASE_ttest  <- ttestBF(formula = BASE ~ group, data = dataset_b_conpat_base_na)
bayes_train_ttest <- ttestBF(formula = TRAIN_change ~ group, data = dataset_b_conpat_train_na)
bayes_cons_ttest  <- ttestBF(formula = CONS_change  ~ group, data = dataset_b_conpat_cons_na)

## remove outliers 3sd from the mean
dataset_b_conpat_rm <- dataset_b_conpat
dataset_b_conpat_med_rm <- dataset_b_conpat_med
dataset_b_patients_rm <- dataset_b_patients
dataset_b_rm <- dataset_b

dataset_b_conpat_rm <- rm_outl(dataset_b_conpat_rm, "BASE")        #1
dataset_b_conpat_rm <- rm_outl(dataset_b_conpat_rm, "TRAIN_change")  #1
dataset_b_conpat_rm <- rm_outl(dataset_b_conpat_rm, "CONS_change") #1
dataset_b_conpat_med_rm <- rm_outl(dataset_b_conpat_med_rm, "BASE")        #1
dataset_b_conpat_med_rm <- rm_outl(dataset_b_conpat_med_rm, "TRAIN_change")  #1
dataset_b_conpat_med_rm <- rm_outl(dataset_b_conpat_med_rm, "CONS_change") #1
#dataset_b_patients_rm <- rm_outl(dataset_b_patients_rm, "BASE")        #1
dataset_b_patients_rm <- rm_outl(dataset_b_patients_rm, "TRAIN_change")  #1
dataset_b_patients_rm <- rm_outl(dataset_b_patients_rm, "CONS_change") #1
dataset_b_rm <- rm_outl(dataset_b_rm, "BASE")        #1
dataset_b_rm <- rm_outl(dataset_b_rm, "TRAIN_change")  #1
dataset_b_rm <- rm_outl(dataset_b_rm, "CONS_change") #1

## Get summary statistics 
dataset_b_conpat_rm_long <- dataset_b_conpat_rm %>% gather(memory, score, BASE, TRAIN_change, CONS_change)
dataset_b_conpat_rm_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")
dataset_b_conpat_med_rm_long <- dataset_b_conpat_med_rm %>% gather(memory, score, BASE, TRAIN_change, CONS_change)
dataset_b_conpat_med_rm_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")
dataset_b_patients_rm_long <- dataset_b_patients_rm %>% gather(memory, score, BASE, TRAIN_change, CONS_change)
dataset_b_patients_rm_long %>% dplyr::select(memory, score, group) %>% group_by(group, memory) %>%  get_summary_stats(score, type = "mean_se")

## plots 
rainplot_baseline_dataset_b <- ggplot(dataset_b_rm,aes(x=group,y=BASE,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=BASE),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Correctly tapped')+ xlab('')+ my_theme +   raincloud_theme + 
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Baseline")+ scale_fill_manual(values = color_dataset_b_all) +  scale_colour_manual(values = color_dataset_b_all) 

rainplot_training_dataset_b  <- ggplot(dataset_b_rm,aes(x=group,y=TRAIN_change,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=TRAIN_change),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Change performance')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Training")+ scale_fill_manual(values = color_dataset_b_all) +  scale_colour_manual(values = color_dataset_b_all) 

rainplot_consolidation_dataset_b  <- ggplot(dataset_b_rm,aes(x=group,y=CONS_change,fill=group,colour=group))+
  geom_flat_violin(position=position_nudge(x=.25,y=0),adjust=2)+
  geom_point(position=position_jitter(width=.05, height = 0, seed = 18),size=4, alpha = 0.5)+
  geom_boxplot(aes(x=as.numeric(group)+0.25,y=CONS_change),
               outlier.shape=NA,alpha=0.3,width=.1,colour="BLACK", size = 1.2) +
  stat_summary(fun.y=mean, geom="point", colour="black", fill="yellow", shape=23, stroke = 1,size=3) +
  ylab('Change performance')+ xlab('')+ my_theme +   raincloud_theme + geom_hline(yintercept=0, linetype="dashed")+
  guides(fill=FALSE,colour=FALSE)+ theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  ggtitle("Consolidation")+ scale_fill_manual(values = color_dataset_b_all) +  scale_colour_manual(values = color_dataset_b_all)

beh_plot_dataset_b     <- ggarrange(rainplot_baseline_dataset_b, rainplot_training_dataset_b, 
                                    rainplot_consolidation_dataset_b, ncol = 3, nrow = 1)

ggsave("figs/beh_plot_dataset_b.png", beh_plot_dataset_b, device = "png", dpi = 300, width = 14, height = 5)


######################################## Relation between sleep and memory consolidation ----

lm_CONS_sleeptable_dataset_b <- lapply(dataset_b_conpat[,c("S1_percent","S2_percent","S3_percent","REM_percent","Wake_after_sleep_onset_percent","Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min")], function(x) lm(dataset_b_conpat$CONS_change ~ x * dataset_b_conpat$group))
lm_CONS_sleeptable_dataset_b_rm <- lapply(dataset_b_rm[,c("S1_percent","S2_percent","S3_percent","REM_percent","Wake_after_sleep_onset_percent","Total_sleep_time_min", "Sleep_Onset_min", "SWS_onset_min", "REM_onset_min")], function(x) lm(dataset_b_rm$CONS_change ~ x * dataset_b_rm$group))

## power and cons
select_vars <- c("1","1.2","1.4","1.6","1.8",
                 "2","2.2","2.4","2.6","2.8",
                 "3","3.2","3.4","3.6","3.8",
                 "4","4.2","4.4","4.6")

mergedData_wide_dataset_b_patients_pair <- mergedData_wide_dataset_b_patients_pair %>% mutate(av_delta = rowMeans(dplyr::select(., select_vars)))
mergedData_wide_dataset_b_patients_pair <- merge(mergedData_wide_dataset_b_patients_pair, tapping, by = "datasetnum")
mergedData_wide_dataset_b_patients_pair %<>% mutate_if(is.character,as.numeric)
mergedData_wide_dataset_b_patients_pair$CONS_change  <- ((mergedData_wide_dataset_b_patients_pair$TEST - mergedData_wide_dataset_b_patients_pair$END)  / mergedData_wide_dataset_b_patients_pair$END)

mergedData_wide_dataset_b_patients_pair$av_delta <- mergedData_wide_dataset_b_patients_pair$av_delta + 1
mergedData_wide_dataset_b_patients_pair$av_delta <- 10*log10(mergedData_wide_dataset_b_patients_pair$av_delta)
is.na(mergedData_wide_dataset_b_patients_pair) <- sapply(mergedData_wide_dataset_b_patients_pair, is.infinite)

select_vars <- c("13.6","13.8","14","14.2","14.4","14.6")

mergedData_wide_dataset_b_conpat_med <- mergedData_wide_dataset_b_conpat_med %>% mutate(av_sigma = rowMeans(dplyr::select(., select_vars)))
mergedData_wide_dataset_b_conpat_med <- merge(mergedData_wide_dataset_b_conpat_med, tapping, by = "datasetnum")
mergedData_wide_dataset_b_conpat_med %<>% mutate_if(is.character,as.numeric)
mergedData_wide_dataset_b_conpat_med$CONS_change  <- ((mergedData_wide_dataset_b_conpat_med$TEST - mergedData_wide_dataset_b_conpat_med$END)  / mergedData_wide_dataset_b_conpat_med$END)

mergedData_wide_dataset_b_conpat_med$av_sigma <- mergedData_wide_dataset_b_conpat_med$av_sigma + 1
mergedData_wide_dataset_b_conpat_med$av_sigma <- 10*log10(mergedData_wide_dataset_b_conpat_med$av_sigma)
is.na(mergedData_wide_dataset_b_conpat_med) <- sapply(mergedData_wide_dataset_b_conpat_med, is.infinite)

## spindle params and cons
lm_CONS_spindle_dataset_b_conpat <- lapply(dataset_b_conpat[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq")], function(x) lm(dataset_b_conpat$CONS_change ~ x * dataset_b_conpat$group))
lm_CONS_spindle_dataset_b_conpat_rm <- lapply(dataset_b_conpat_rm[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq")], function(x) lm(dataset_b_conpat_rm$CONS_change ~ x * dataset_b_conpat_rm$group))

lm_CONS_spindle_dataset_b_conpat_med <- lapply(dataset_b_conpat_med[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq")], function(x) summary(lm(dataset_b_conpat_med$CONS_change ~ x * dataset_b_conpat_med$group)))
lm_CONS_spindle_dataset_b_patients <- lapply(dataset_b_patients[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq")], function(x) summary(lm(dataset_b_patients$CONS_change ~ x * dataset_b_patients$group)))

## coupling and cons
lm_CONS_coupling_dataset_b <- lapply(dataset_b_conpat[,c("coupling_occur","mean_delay","sd_delay")], function(x) lm(dataset_b_conpat$CONS_change ~ x * dataset_b_conpat$group))
lm_CONS_coupling_dataset_b_rm <- lapply(dataset_b_conpat_rm[,c("coupling_occur","mean_delay","sd_delay")], function(x) lm(dataset_b_conpat_rm$CONS_change ~ x * dataset_b_conpat_rm$group))

dataset_b_rm %>% 
  group_by(group) %>%
  filter(!is.na(CONS_change)) %>%
  dplyr::summarize(cor(CONS_change, coupling_occur))

dataset_b_rm %>% 
  group_by(group) %>%
  filter(!is.na(CONS_change)) %>%
  dplyr::summarize(cor(CONS_change, sd_delay))


## Correlate hippocampal size with CONS
dataset_b_conpat_hipp_CONS <- merge(dataset_b_rm, hipp_all, by.x=c("id_loreta"), by.y=c("loreta_id"), all.x = T)
dataset_b_pat_hipp_CONS <- dataset_b_conpat_hipp_CONS %>% dplyr::filter(group == "patient_unmed")
dataset_b_con_hipp_CONS <- dataset_b_conpat_hipp_CONS %>% dplyr::filter(group == "control")

dataset_b_pat_hipp_CONS$hipp_size_av <- (dataset_b_pat_hipp_CONS$hipp_size_left + dataset_b_pat_hipp_CONS$hipp_size_right)/2
dataset_b_con_hipp_CONS$hipp_size_av <- (dataset_b_con_hipp_CONS$hipp_size_left + dataset_b_con_hipp_CONS$hipp_size_right)/2
dataset_b_conpat_hipp_CONS$hipp_size_av <- (dataset_b_conpat_hipp_CONS$hipp_size_left + dataset_b_conpat_hipp_CONS$hipp_size_right)/2

## Any interaction effects for all sleep parameters?
lapply(dataset_b_pat_hipp_CONS[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq","so_dens","so_dur","so_amp", "so_count", "so_freq" , "coupling_occur","mean_delay","sd_delay")], function(x) summary(lm(x ~ dataset_b_pat_hipp_CONS$hipp_size_av)))
lapply(dataset_b_con_hipp_CONS[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq","so_dens","so_dur","so_amp", "so_count", "so_freq" , "coupling_occur","mean_delay","sd_delay")], function(x) summary(lm(x ~ dataset_b_con_hipp_CONS$hipp_size_av)))
lapply(dataset_b_conpat_hipp_CONS[,c("spin_dens","spin_count","spin_amp", "spin_dur", "spin_freq","so_dens","so_dur","so_amp", "so_count", "so_freq" , "coupling_occur","mean_delay","sd_delay")], function(x) summary(lm(x ~ dataset_b_conpat_hipp_CONS$hipp_size_av*dataset_b_conpat_hipp_CONS$group)))

## specific correlation per group
cor.test(dataset_b_con_hipp_CONS$spin_freq, dataset_b_con_hipp_CONS$hipp_size_av)
cor.test(dataset_b_pat_hipp_CONS$spin_freq, dataset_b_pat_hipp_CONS$hipp_size_av)

cor.test(dataset_b_con_hipp_CONS$mean_delay, dataset_b_con_hipp_CONS$hipp_size_av)
cor.test(dataset_b_pat_hipp_CONS$mean_delay, dataset_b_pat_hipp_CONS$hipp_size_av)

## Adjust for plotting
dataset_b_conpat_hipp_CONS  <- dataset_b_conpat_hipp_CONS %>% dplyr::filter(group == "control" | group == "patient_unmed")
dataset_b_conpat_hipp_CONS$group<- revalue(dataset_b_conpat_hipp_CONS$group, c("control" = "Controls",
                                                       "patient_unmed" = "Unmedicated"))
## plot
plot_hipp_spinfreq <- ggplot(dataset_b_conpat_hipp_CONS, aes(x = hipp_size_av , y = spin_freq, colour = group)) + 
  geom_point( size = 4) + geom_smooth(method = "lm", se = FALSE, fullrange=TRUE, size = 1.2) + my_theme +
  ylab("Spindle frequency [Hz]") + 
  scale_color_manual(values=color_dataset_b) + xlim(3000, 5000)+
  scale_fill_manual(values=color_dataset_b) + guides(fill=FALSE,colour=FALSE) +
  labs(x = expression ("Average hippocampal volume in"~mm^3 ))
 
plot_hipp_meandelay <-  ggplot(dataset_b_conpat_hipp_CONS, aes(x = hipp_size_av , y = mean_delay, colour = group)) + 
  geom_point( size = 4) + geom_smooth(method = "lm", se = FALSE, fullrange=TRUE, size = 1.2) + my_theme +
  ylab("Mean delay [s]") + 
  scale_color_manual(values=color_dataset_b, name = "Groups") + 
  scale_fill_manual(values=color_dataset_b) + xlim(3000, 5000)+
  theme( legend.position = c(0.82, 0.1),
         legend.title=element_text(size=13),
         legend.text =element_text(size=13)) +
  labs(x = expression ("Average hippocampal volume in"~mm^3 ))

plot_hipp <- ggarrange(plot_hipp_spinfreq, plot_hipp_meandelay, ncol = 2,labels = c("A", "B"))

# ggsave("figs/plot_hipp.png", plot_hipp, device = "png", dpi = 300, width = 10, height = 7)
# ggsave("figs/plot_hipp.pdf", plot_hipp, device = "pdf", dpi = 300, width = 10, height = 7)


## Depression severity ----
dataset_b_hamd_long <- dataset_b %>% filter(group=="patient_unmed") %>% gather(hamd, score, hamd_0, hamd_week1, hamd_week4)
dataset_b_hamd_long %>% dplyr::select(hamd, score, group) %>% group_by(hamd) %>%  get_summary_stats(score, type = "mean_se")

dataset_b_unmed<- dataset_b  %>% filter(group=="patient_unmed")

## PSQI ----
dataset_b_conpat %>%
  dplyr::select(PSQI_sum, group) %>%
  group_by(group) %>% get_summary_stats(PSQI_sum, type = "mean_se")

psqi_ttest_dataset_b_conpat <- t.test(dataset_b_conpat$PSQI_sum ~ dataset_b_conpat$group, var.equal = FALSE)
dataset_b_conpat_psqi <- dataset_b_conpat %>% drop_na(PSQI_sum) 

dataset_b_conpat_med %>%
  dplyr::select(PSQI_sum, group) %>%
  group_by(group) %>% get_summary_stats(PSQI_sum, type = "mean_se")

psqi_dataset_b_conpat_med <- t.test(dataset_b_conpat_med$PSQI_sum ~ dataset_b_conpat_med$group, var.equal = FALSE)
dataset_b_conpat_med_psqi <- dataset_b_conpat_med %>% drop_na(PSQI_sum) 


dataset_b_patients %>%
  dplyr::select(PSQI_sum, group) %>%
  group_by(group) %>% get_summary_stats(PSQI_sum, type = "mean_se")

dataset_b_patients_NA <- dataset_b_patients_NA[complete.cases(dataset_b_patients_NA),]
dataset_b_patients_NA <- dataset_b_patients_NA[complete.cases(dataset_b_patients_NA[c(63, 58),]),]

t.test(dataset_b_patients_NA$PSQI_sum ~ dataset_b_patients_NA$group, var.equal = FALSE, paired = T)

psqi_dataset_b_patients <- t.test(dataset_b_patients$PSQI_sum ~ dataset_b_patients$group, var.equal = FALSE)
dataset_b_conpat_patients_psqi <- dataset_b_patients %>% drop_na(PSQI_sum) 
bf_dataset_b_conpat_patients_psqi  <- ttestBF(formula = PSQI_sum ~ group, data = dataset_b_conpat_patients_psqi)

cor.test(dataset_b$PSQI_sum, dataset_b$Wake_after_sleep_onset_percent)
cor.test(dataset_b$PSQI_sum, dataset_b$Sleep_Onset_min)
cor.test(dataset_b$PSQI_sum, dataset_b$Total_sleep_time_min)
cor.test(dataset_b$PSQI_sum, dataset_b$SWS_min)

