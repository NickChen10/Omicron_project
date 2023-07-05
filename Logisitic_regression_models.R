#Description: Logistic regression modeling for periods of variant co-circulation 
#Author: Nicholas Chen

#Load required packages
library(tidyverse)
library(lme4)
library(lubridate)
library(ggeffects)
library(snakecase)
library(cowplot)



##### DATA FORMATTING #####

# Read in new metadata sheet and subset to fields of interest 
metadata <- read.csv("C:/Users/nicho/Box Sync/Yale Virology sample intake/vaccination_metadata_database/vaccine_metadata_clean.csv")
metadata <- metadata %>% select(-c(X,patient_class,country,city,zip,tube_label,original_sample_type,report_ct,extraction_pcr_date,yale_n1_fam:yale_n1_ge_ml, collection_date_gsheet)) %>% rename(collection_date = collection_date_virology)
metadata[,c("dob","collection_date","vax_1_date","vax_2_date","booster_1_date","booster_2_date","series_complete_date")] <- 
  lapply(metadata[,c("dob","collection_date","vax_1_date","vax_2_date","booster_1_date","booster_2_date","series_complete_date")], as.Date,format="%m/%d/%Y")

# Apply filters (Clinical virology, Connecticut residents, non-missing MPI/MRN/Lineage info/dob, first instance of duplicate MRNs, consistent vaccination records 
metadata <- metadata %>% filter(source == "Yale Clinical Virology Lab" & state == "Connecticut" & filter == "connecticut" 
                                & metadata == T & inconsistent_vax_info == F & dup_mpi_gsheet==F & miscoded_mpi == F & !is.na(lineage) 
                                & !is.na(dob)) %>% filter(!duplicated(mrn))

# Create interval column 
metadata$Interval <- ""

# Categorize vaccination status
metadata <- metadata %>% 
  mutate(vax_status = case_when(
    difftime(collection_date, vax_1_date, units = "days") < 14 ~ "Non_Breakthrough",
    difftime(collection_date, vax_1_date, units = "days") >= 14 ~ "One_Dose_Breakthrough",
    difftime(collection_date, vax_2_date, units = "days") >= 14 ~ "Two_Dose_Breakthrough",
    difftime(collection_date, booster_1_date, units = "days") >= 14 ~ "Three_Dose_Breakthrough",
    difftime(collection_date, booster_2_date, units = "days") >= 14 ~ "Four_Dose_Breakthrough",
    TRUE ~ "Non_Breakthrough"
  ))
  

#Time since vax (also shitty old code)
metadata$time_since_vax <- ""
for(i in 1:nrow(metadata)){
  if(metadata[i,"vax_status"]=='One_Dose_Breakthrough'){
    metadata[i,"time_since_vax"] <- difftime(metadata[i,"collection_date"],metadata[i,"vax_1_date"],units='days')
  } else if(metadata[i,"vax_status"]=='Two_Dose_Breakthrough'){
    metadata[i,"time_since_vax"] <- difftime(metadata[i,"collection_date"],metadata[i,"vax_2_date"],units='days')
  } else if(metadata[i,"vax_status"]=='Three_Dose_Breakthrough'){
    metadata[i,"time_since_vax"] <- difftime(metadata[i,"collection_date"],metadata[i,"booster_1_date"],units='days')
  } else if(metadata[i,"vax_status"]=='Four_Dose_Breakthrough'){
    metadata[i,"time_since_vax"] <- difftime(metadata[i,"collection_date"],metadata[i,"booster_2_date"],units='days')
  } else { 
    NULL
  }
}

#LEFT OFF HERE

#Time since vax as numeric
metadata$time_since_vax <- as.numeric(metadata$time_since_vax)

#Updated vaccination categories (finally some good code)
for(i in 1:nrow(metadata)){
  if(metadata[i,"vax_status"] == "Non_Breakthrough"){
    metadata[i,"vax_status_revised"]  <- "non_breakthrough"
  } else if(metadata[i,"vax_status"] == "One_Dose_Breakthrough"){
    metadata[i,"vax_status_revised"]  <- "one_dose_breakthrough"
  } else if(metadata[i,"vax_status"] %in% c("Two_Dose_Breakthrough","Three_Dose_Breakthrough","Four_Dose_Breakthrough") & metadata[i,"time_since_vax"] <150){
    metadata[i,"vax_status_revised"]  <- "complete_breakthrough_<5_months"
  } else if(metadata[i,"vax_status"] %in% c("Two_Dose_Breakthrough","Three_Dose_Breakthrough","Four_Dose_Breakthrough") & metadata[i,"time_since_vax"] >=150){
    metadata[i,"vax_status_revised"]  <- "complete_breakthrough_>=5_months"
  }
}