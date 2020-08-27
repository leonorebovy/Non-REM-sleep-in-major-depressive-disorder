## extra script to calculate hippocampal sizes from anatomical MRI data
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

library(tidyverse)

### Read in hippocampal size values

## Controls ----

options(digits=10)

hipp_env <- new.env()
with(hipp_env, { 
    ## Specify the amount of folders containing Loreta subjects 
    controls <- list.files(path = ".../raw_mri/healthy/", pattern = "^LK_")
    
    amount_controls_subj <- length(controls)
    
    ## Initialize empty vector
    hipp_control_left_vector = 0 
    hipp_control_right_vector = 0 
    
    ## For loop over the folder content to make a list of file names
    for (subjn in 1:amount_controls_subj){ #last Loreta ID number
      hipp_control_left <- read.delim(paste(".../raw_mri/healthy/",
                                           controls[subjn],
                                           "/FIRST/volume_Hipp_left.txt",sep = ""),header = FALSE, sep = " ")
      
      hipp_control_right <- read.delim(paste(".../raw_mri/healthy/",
                                            controls[subjn],
                                            "/FIRST/volume_Hipp_right.txt",sep = ""),header = FALSE, sep = " ")
      
      hipp_control_left_vector <- append(hipp_control_left_vector, hipp_control_left[2])
      hipp_control_right_vector <- append(hipp_control_right_vector, hipp_control_right[2])
    }
    
    ## Remove the first inital 0 again
    hipp_control_left_vector <- hipp_control_left_vector[-1]
    hipp_control_right_vector <- hipp_control_right_vector[-1]
    
    hipp_controls_left        <- as.data.frame(hipp_control_left_vector)
    hipp_controls_left        <- as.data.frame(t(hipp_controls_left))
    hipp_controls_right       <- as.data.frame(hipp_control_right_vector)
    hipp_controls_right       <- as.data.frame(t(hipp_controls_right))
    
    loreta_controls   <- as.data.frame(controls)
    
    hipp_controls        <- cbind(loreta_controls,hipp_controls_left,hipp_controls_right)
    
    rownames(hipp_controls) <- NULL
    names(hipp_controls)[1] <- "loreta_id"
    names(hipp_controls)[2] <- "hipp_size_left"
    names(hipp_controls)[3] <- "hipp_size_right"
    
    ## save the output
    #write.table(hipp_controls, file = "hipp_size_controls.txt", sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)
    
    
## Patients ----
    
    ## Specify the amount of folders containing Loreta subjects 
    patients <- list.files(path = ".../raw_mri/patients/", pattern = "^LP_")
    
    amount_patients_subj <- length(patients)
    
    ## Initialize empty vector
    hipp_patients_left_vector = 0 
    hipp_patients_right_vector = 0 
    
    ## For loop over the folder content to make a list of file names
    for (subjn in 1:amount_patients_subj){ #last Loreta ID number
      hipp_patients_left <- read.delim(paste(".../raw_mri/patients/",
                                            patients[subjn],
                                            "/FIRST/volume_Hipp_left.txt",sep = ""),header = FALSE, sep = " ")
      
      hipp_patients_right <- read.delim(paste(".../raw_mri/patients/",
                                             patients[subjn],
                                             "/FIRST/volume_Hipp_right.txt",sep = ""),header = FALSE, sep = " ")
      
      hipp_patients_left_vector <- append(hipp_patients_left_vector, hipp_patients_left[1,2])
      hipp_patients_right_vector <- append(hipp_patients_right_vector, hipp_patients_right[1,2])
    }
    
    ## Remove the first inital 0 again
    hipp_patients_left_vector <- hipp_patients_left_vector[-1]
    hipp_patients_right_vector <- hipp_patients_right_vector[-1]
    
    hipp_patients_left        <- as.data.frame(hipp_patients_left_vector)
    #hipp_patients_left        <- as.data.frame(t(hipp_patients_left))
    hipp_patients_right       <- as.data.frame(hipp_patients_right_vector)
    #hipp_patients_right       <- as.data.frame(t(hipp_patients_right))
    
    loreta_patients   <- as.data.frame(patients)
    
    hipp_patients        <- cbind(loreta_patients,hipp_patients_left,hipp_patients_right)
    
    rownames(hipp_patients) <- NULL
    names(hipp_patients)[1] <- "loreta_id"
    names(hipp_patients)[2] <- "hipp_size_left"
    names(hipp_patients)[3] <- "hipp_size_right"
    
    ## save the output
    #write.table(hipp_patients, file = "hipp_size_patients.txt", sep = ",", dec = ".", row.names = TRUE, col.names = TRUE)
    
    
    # hipp_patients$group <- "patient"
    # hipp_controls$group <- "control"
     
    hipp_all <- rbind(hipp_patients, hipp_controls)
  
      ## Create matched pair values
      hipp_pairs <- hipp_all[FALSE,]
      names(hipp_pairs)[1] <- "loreta_id"
      names(hipp_pairs)[2] <- "hipp_size_left"
      names(hipp_pairs)[3] <- "hipp_size_right"
     
     hipp_all_pairs <- hipp_all %>% dplyr::filter(loreta_id  == "LK_14" | loreta_id  == "LP_06" |
                                                  loreta_id  == "LK_15" | loreta_id  == "LP_38" |
                                                  loreta_id  == "LK_18" | loreta_id  == "LP_41" |
                                                  loreta_id  == "LK_20" | loreta_id  == "LP_05" |
                                                  loreta_id  == "LK_31" | loreta_id  == "LP_02" |
                                                  loreta_id  == "LK_33" | loreta_id  == "LP_58" |
                                                  loreta_id  == "LK_35" | loreta_id  == "LP_81" |
                                                  loreta_id  == "LK_48" | loreta_id  == "LP_29" |
                                                  loreta_id  == "LK_55" | loreta_id  == "LP_23" |
                                                  loreta_id  == "LK_63" | loreta_id  == "LP_52" |
                                                  loreta_id  == "LK_64" | loreta_id  == "LP_22" |
                                                  loreta_id  == "LK_67" | loreta_id  == "LP_74" |
                                                  loreta_id  == "LK_68" | loreta_id  == "LP_42" |
                                                  loreta_id  == "LK_85" | loreta_id  == "LP_03"  )
     
     hipp_all_pairs$pairing <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14, 5, 14, 4, 1, 11, 9, 8, 2, 3, 13, 10, 6, 12, 7)

     })

hipp_all <- with(hipp_env, hipp_all)
hipp_all_pairs <- with(hipp_env, hipp_all_pairs)
 
hipp_all$loreta_id <- str_remove_all(hipp_all$loreta_id, "[LKP_]")
hipp_all$loreta_id <- substr(hipp_all$loreta_id,regexpr("[^0]",hipp_all$loreta_id),nchar(hipp_all$loreta_id))

hipp_all_pairs$loreta_id <- str_remove_all(hipp_all_pairs$loreta_id, "[LKP_]")
hipp_all_pairs$loreta_id <- substr(hipp_all_pairs$loreta_id,regexpr("[^0]",hipp_all_pairs$loreta_id),nchar(hipp_all_pairs$loreta_id))

