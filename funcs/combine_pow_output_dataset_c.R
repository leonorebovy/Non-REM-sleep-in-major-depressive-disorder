## extra script to combine the spisop power analysis into one csv file 
## it was too large  (= too much RAM needed) to do on it's own
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

library(plyr)
library(tidyverse)

power  <- read.delim(".../dataset_c_pow_full_channels_datanum_1.csv", sep=",",dec= ".") 
power_colnames <- colnames(power)

no_subjects <- 2:118

for(subnr in no_subjects){
  power_subj <- read.delim(paste(".../dataset_c_run_1_pow_full_channels_datanum_",subnr,".csv",sep=""), sep=",", dec= ".")
  power <- rbind(power, power_subj)
}

## drop unused levels
power[] <- lapply(power, function(x) if(is.factor(x)) factor(x) else x)

## rename 
power_dataset_b <- power

## change to numeric
power_dataset_b$freq <- as.numeric(as.character(power_dataset_b$freq))
power_dataset_b$mean_powerDensity_over_segments <- as.numeric(as.character(power_dataset_b$mean_powerDensity_over_segments))

## Rename the eeg channels
power_dataset_b <- mutate(power_dataset_b, channel = fct_recode(channel, "EEGC3" = "C3-TP10", "EEGC4" = "C4-TP9"))

## clean up
rm(power, power_subj, no_subjects, power_colnames, subnr)
