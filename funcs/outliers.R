## Functions to check for outliers and to remove them (replace with NA)
## Writen for >3sd
## Project: Non-REM sleep in major depressive disorder
## Author: Leonore Bovy

##input is dataframe and variable of interest

check_outl  <- function(dataframe, x) {
  nrow(dataframe[which(dataframe[[x]] < mean(dataframe[[x]] , na.rm = T) - 3*sd(dataframe[[x]] , na.rm = T) | dataframe[[x]]  > mean(dataframe[[x]] , na.rm = T) + 3*sd(dataframe[[x]]  ,na.rm = T)),])}

which_outl  <- function(dataframe, x) {
  dataframe[which(dataframe[[x]] < mean(dataframe[[x]] , na.rm = T) - 3*sd(dataframe[[x]] , na.rm = T) | dataframe[[x]]  > mean(dataframe[[x]] , na.rm = T) + 3*sd(dataframe[[x]]  ,na.rm = T)),]}

rm_outl  <- function(dataframe, x) {
  dataframe[which(dataframe[[x]] < mean(dataframe[[x]] , na.rm = T) - 3*sd(dataframe[[x]] , na.rm = T) | dataframe[[x]]  > mean(dataframe[[x]] , na.rm = T) + 3*sd(dataframe[[x]]  ,na.rm = T)),][[x]] <- NA
  return(dataframe)
}