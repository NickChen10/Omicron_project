#Description: Logistic regression modeling for periods of variant co-circulation 
#Author: Nicholas Chen

#Load required packages
library(tidyverse)
library(lme4)
library(lubridate)
library(ggeffects)
library(snakecase)
library(cowplot)
library(table1)




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
    difftime(collection_date, booster_2_date, units = "days") >= 14 ~ "Four_Dose_Breakthrough",
    difftime(collection_date, booster_1_date, units = "days") >= 14 ~ "Three_Dose_Breakthrough",
    difftime(collection_date, vax_2_date, units = "days") >= 14 ~ "Two_Dose_Breakthrough",
    difftime(collection_date, vax_1_date, units = "days") >= 14 ~ "One_Dose_Breakthrough",
    TRUE ~ "Non_Breakthrough"
  ))
  

# Calculate time from most recent vaccination to collection date
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

metadata$time_since_vax <- as.numeric(metadata$time_since_vax)

# Update vaccination statuses
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

# Need to combine this data with an older dataset, so restrict this dataset to the last date of the previous dataset 
metadata <- metadata %>% filter(collection_date >= "2022-09-01")

# Recode lineages 
metadata <- metadata %>%
  mutate(ShortLin = case_when(
    endsWith(lineage, "2.12.1") ~ "BA.2",
    startsWith(lineage, "BG") ~ "BA.2",
    startsWith(lineage, "BA.2.75") ~ "BA.2",
    startsWith(lineage, "BA.1.1") ~ "BA.1",
    startsWith(lineage, "BA.1") ~ "BA.1",
    startsWith(lineage, "BA.2") ~ "BA.2",
    startsWith(lineage, "BK") ~ "BA.2",
    startsWith(lineage, "BP") ~ "BA.2",
    startsWith(lineage, "BL") ~ "BA.2",
    startsWith(lineage, "BH") ~ "BA.2",
    startsWith(lineage, "BM") ~ "BA.2",
    startsWith(lineage, "CH") ~ "BA.2",
    startsWith(lineage, "BJ") ~ "BA.2",
    startsWith(lineage, "BR") ~ "BA.2",
    startsWith(lineage, "BA.4.6") ~ "BA.4",
    startsWith(lineage, "BA.4") ~ "BA.4",
    startsWith(lineage, "BA.5") ~ "BA.5",
    startsWith(lineage, "BU") ~ "BA.5",
    startsWith(lineage, "BV") ~ "BA.5",
    startsWith(lineage, "BE") ~ "BA.5",
    startsWith(lineage, "BT") ~ "BA.5",
    startsWith(lineage, "BF.7") ~ "BA.5", #BF.7 as BA.5 now
    startsWith(lineage, "BF") ~ "BA.5",
    startsWith(lineage, "BQ.1.1") ~ "BA.5", #BQ.1.1 as BA.5 now 
    startsWith(lineage, "BQ.2") ~ "BA.5",
    startsWith(lineage, "BQ.3") ~ "BA.5",
    startsWith(lineage, "CN") ~ "BA.5",
    startsWith(lineage, "CD") ~ "BA.5",
    startsWith(lineage, "CE") ~ "BA.5",
    startsWith(lineage, "CF") ~ "BA.5",
    startsWith(lineage, "CL") ~ "BA.5",
    startsWith(lineage, "CG") ~ "BA.5",
    startsWith(lineage, "CK") ~ "BA.5",
    startsWith(lineage, "BZ") ~ "BA.5",
    startsWith(lineage, "XBB.1.5") ~ "XBB.1",
    startsWith(lineage, "XBB.1") ~ "XBB.1", 
    TRUE ~ "Other"
  ))

# Format approximate age
birth_year<-as.numeric(format(as.Date(metadata$dob, format="%Y-%m-%d"),"%Y"))
collection_year<-as.numeric(format(as.Date(metadata$collection_date, format="%Y-%m-%d"),"%Y"))
approx_age<-collection_year -
  birth_year
metadata$age <- approx_age
rm(approx_age,birth_year, collection_year,i)


##### COMBINE DATASETS #####
# Append the current dataset to a previous dataset that extended further back in time


# Read in old dataset and format
old_dataset <- read.csv("C:/Users/nicho/OneDrive/Documents/Yale Docs/Thesis_revised/Data/VBT_complete_data.csv")
names(old_dataset) <- to_any_case(names(old_dataset), case = "snake")

# Retain only the fields of interest for both datasets 
metadata_clean<- metadata %>% select(sample_id,collection_date,mpi, age,gender,city_gsheet,ShortLin,vax_status_revised) %>% 
  rename(town = city_gsheet, variant = ShortLin)
old_dataset_clean<- old_dataset %>% select(sample_id,collection_date,mpi, age,gender,city,lineage_name,vax_status_revised) %>% 
  rename(town = city, variant = lineage_name)
old_dataset_clean$collection_date <- as.Date(old_dataset_clean$collection_date, format="%Y-%m-%d")

# Combine datsets  
data <- rbind(old_dataset_clean, metadata_clean)
rm(old_dataset,old_dataset_clean,metadata,metadata_clean)


##### MODELING ######

for(i in 1:5){ #Run for loop once to generate 5 models at the 5 week duration
  data$Interval <- NA
  
  #Pre-Delta to Delta date intervals
  interval_1_dates <- data.frame(start=ymd(c("2021-06-09","2021-06-07","2021-06-03","2021-05-29","2021-05-26","2021-05-21")),
                                 end=ymd(c("2021-06-30","2021-07-05","2021-07-08","2021-07-10","2021-07-14","2021-07-16")))
  
  interval_1 <- interval(interval_1_dates[3,1], interval_1_dates[3,2]) #Change row number to date interval of interest (5 weeks = 3)
  data$Interval[which(data$collection_date %within% interval_1)] <- 1
  
  #Delta to BA.1 date intervals
  interval_2_dates <- data.frame(start=ymd(c("2021-12-09","2021-12-05","2021-11-30","2021-11-25","2021-11-20","2021-11-18")),
                                 end=ymd(c("2021-12-30","2022-01-02","2022-01-04","2022-01-06","2022-01-08","2022-01-13")))
  
  interval_2 <- interval(interval_2_dates[3,1], interval_2_dates[3,2]) #Change row number to date interval of interest (5 weeks = 3)
  data$Interval[which(data$collection_date %within% interval_2)] <- 2
  
  #BA.1 to BA.2 date intervals
  interval_3_dates <- data.frame(start=ymd(c("2022-03-03","2022-02-28","2022-02-23","2022-02-20","2022-02-19","2022-02-15")),
                                 end=ymd(c("2022-03-24","2022-03-28","2022-03-30","2022-04-03","2022-04-09","2022-04-12")))
  
  interval_3 <- interval(interval_3_dates[3,1], interval_3_dates[3,2]) #Change row number to date interval of interest (5 weeks = 3)
  data$Interval[which(data$collection_date %within% interval_3)] <- 3
  
  #BA.2 to BA.4&5 date intervals
  interval_4_dates <- data.frame(start=ymd(c("2022-06-14","2022-06-11","2022-06-09","2022-06-05","2022-06-02","2022-05-30")),
                                 end=ymd(c("2022-07-05","2022-07-09","2022-07-14","2022-07-17","2022-07-21","2022-07-26")))
  
  interval_4 <- interval(interval_4_dates[3,1], interval_4_dates[3,2]) #Change row number to date interval of interest (5 weeks = 3)
  data$Interval[which(data$collection_date %within% interval_4)] <- 4
  
  
  #Others to XBB.1.5 date intervals
  interval_5_dates <- data.frame(start=ymd(c("2022-12-12","2022-12-08","2022-12-07","2022-12-03","2022-11-30","2022-11-27")),
                                 end=ymd(c("2023-01-02","2023-01-05","2023-01-11","2023-01-14","2023-01-18","2023-01-22")))
  
  interval_5 <- interval(interval_5_dates[3,1], interval_5_dates[3,2]) #Change row number to date interval of interest (5 weeks = 3)
  data$Interval[which(data$collection_date %within% interval_5)] <- 5
  
  #Selects intervals sequentially
  interval <- data %>% filter(Interval==i)
  interval <- interval %>% filter(age >=5) #Restrict to those over 5 years old 
  
  #Covariate: age 
  age_5_17 <- as.numeric((interval$age <= 17) & (interval$age >= 5))
  age_18_39 <- as.numeric((interval$age <= 39) & (interval$age >= 18))
  age_40_64 <- as.numeric((interval$age <= 64) & (interval$age >= 40))
  age_65_plus <- as.numeric(interval$age >= 65)
  age_cat<-cbind(age_5_17,age_40_64,age_65_plus) 
  
  
  #Dichotomize outcome variable (variant)
  if(unique(interval$Interval)==1){
    variant<-0+as.numeric(interval$variant == "Delta")
  } else if (unique(interval$Interval)==2){
    variant<-1-as.numeric(interval$variant == "Delta")
  } else if (unique(interval$Interval)==3){
    variant<-1-as.numeric(interval$variant == "Omicron_BA.1")
  } else if (unique(interval$Interval)==4){
    variant<-1-as.numeric(interval$variant == "Omicron_BA.2")
  } else if (unique(interval$Interval)==5){
    variant<-0+as.numeric(interval$variant == "XBB.1")
  }
  
  #Primary predictor variable: vaccination status
  vax_status <-model.matrix(~0 + interval$vax_status_revised)
  vax_status <- vax_status[,c("interval$vax_status_revisednon_breakthrough","interval$vax_status_revisedone_dose_breakthrough","interval$vax_status_revisedcomplete_breakthrough_>=5_months","interval$vax_status_revisedcomplete_breakthrough_<5_months")]
  colnames(vax_status) <- c("non_breakthrough","one_dose_breakthrough","complete_breakthrough_>=5_months","complete_breakthrough_<5_months")
  vax_status<-vax_status[,-1]
  
  
  #Covariate: time
  time<-scale(as.numeric(as.Date(interval$collection_date)))
  
  #Covariate: sex
  gender <- interval$gender
  gender[which(!(gender %in% c("M","F")))] <- NA
  gender[which(gender=="M")] <- 0
  gender[which(gender=="F")] <- 1
  
  #Random effect covariate: town of residence
  town <-as.factor(interval$town)
  
  #Mixed effect multivariate logistic regression model
  if(unique(interval$Interval)==4){
    model<-glmer(variant ~ vax_status + 
                   time + 
                   age_cat +  
                   #gender + #Removed for interval 4 due to missing data
                   (1|town), 
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa"))
  } else{
    model<-glmer(variant ~ vax_status + 
                   time + 
                   age_cat +  
                   gender+ #removed for interval 4
                   (1|town), 
                 family = binomial(link = "logit"),
                 control = glmerControl(optimizer = "bobyqa"))
  }
  
  
  #Pull exponentiated model coefficients for plotting 
  for(x in 2:4){
    est<-exp(summary(model)$coef[x])
    ll<-exp(summary(model)$coef[x]-(1.96*summary(model)$coef[x,2]))
    ul<-exp(summary(model)$coef[x]+(1.96*summary(model)$coef[x,2]))
    assign(paste0("est_",x),est)
    assign(paste0("ll_",x),ll)
    assign(paste0("ul_",x),ul)
  }
  
  #Format plot data
  if(unique(interval$Interval)==1){
    plot_data <- data.frame(
      Index = c(1:3),
      label = c("1 dose",'2+ doses \n\u2265 5 months','2+ doses \n <5 months'),
      OR = c(NA,est_3,est_4), #Interval 1 could not provide an estimate for the partial vaccination variable level, define as NA
      LL = c(NA,ll_3,ll_4),
      UL = c(NA,ul_3,ul_4))
  } else{
    plot_data <- data.frame(
      Index = c(1:3),
      label = c("1 dose",'2+ doses \n\u2265 5 months','2+ doses \n <5 months'),
      OR = c(est_2,est_3,est_4),
      LL = c(ll_2,ll_3,ll_4),
      UL = c(ul_2,ul_3,ul_4))
  }
  
  
  
  #Plotting names
  p_names <- list("Delta Emergence","BA.1 Emergence","BA.2 Emergence","BA.4/5 Emergence","XBB.1 Emergence")
  

  #Plot colors 
  
  cols<- c("#ff6e40","#184e27","#CE4E50","#7294D4","#892F2E")
  
  #Forest plots
  p <- ggplot(plot_data, aes(y = Index, x = OR)) +
    geom_point(shape = 18, size = 5,color=cols[[i]]) +  
    geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25,size=1, color=cols[[i]]) +
    geom_vline(xintercept = 1, color = "black", linetype = "dashed", cex = 1, alpha = 0.5) +
    scale_y_continuous(name = "", breaks=1:3, labels = plot_data$label, trans = "reverse") +
    scale_x_log10()+
    xlab("Odds Ratio (95% CI)") + 
    ylab("") +
    ggtitle(paste0(p_names[[i]]))+
    theme_bw() +
    theme(plot.title = element_text(size=14, hjust = 0.5, face="bold"),
          panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black",face="bold"),
          axis.text.x.bottom = element_text(size = 12, colour = "black",face="bold"),
          axis.title.x = element_text(size = 12, colour = "black",face="bold")) 
  
  assign(paste0(p_names[[i]],"_plot"),p)  
  assign(paste0(p_names[[i]],"_model"),model)  
  
  
} #Ignore warnings


# Combine plots 
plot_grid(`Delta Emergence_plot`+xlab(""),
          `BA.1 Emergence_plot`+xlab(""),
          `BA.2 Emergence_plot`,
          `BA.4/5 Emergence_plot`,
          `XBB.1 Emergence_plot`)
#ggsave("forest_plots.png",width=30, height=20, units = "cm")



##### TABLE #####
table <- data
table <- table %>% filter(age>=5)

for(i in 1:nrow(table)){
  if(table$age[i] %in% c(5:17)){
    table$age_cat[i] <- "5-17"
  } else if(table$age[i] %in% c(18:39)){
    table$age_cat[i] <- "18-39"
  } else if(table$age[i] %in% c(40:64)){
    table$age_cat[i] <- "40-64"
  } else if(table$age[i]>=65){
    table$age_cat[i] <- "65+"
  }
}

table$gender[which(!(table$gender %in% c("F","M", NA)))] <- "Other"
table$gender <- factor(table$gender, levels=c("F","M","Other"))
table$age_cat <- factor(table$age_cat, levels=c("5-17","18-39","40-64","65+"))
table$vax_status_revised <- factor(table$vax_status_revised, levels=c("non_breakthrough","one_dose_breakthrough","complete_breakthrough_>=5_months","complete_breakthrough_<5_months"))

table1(~age_cat + gender | vax_status_revised,data=table)


models<- list(`Delta Emergence_model`,`BA.1 Emergence_model`,`BA.2 Emergence_model`,`BA.4/5 Emergence_model`,`XBB.1 Emergence_model`)




for(x in 1:9){
    model <- `Delta Emergence_model`
    est<-exp(summary(model)$coef[x])
    ll<-exp(summary(model)$coef[x]-(1.96*summary(model)$coef[x,2]))
    ul<-exp(summary(model)$coef[x]+(1.96*summary(model)$coef[x,2]))
    assign(paste0(x,"_est"),est)
    assign(paste0(x,"_ll"),ll)
    assign(paste0(x,"_ul"),ul)  
}
 
model_results <- data.frame(name=c("int","one dose", ">5 months", "<5 months",
                                   "time","5-17","40-64","65+","gender"),
                            est=c(`1_est`,`2_est`,`3_est`,`4_est`,`5_est`,
                                  `6_est`,`7_est`,`8_est`,`9_est`), 
                            ll=c(`1_ll`,`2_ll`,`3_ll`,`4_ll`,`5_ll`,
                                 `6_ll`,`7_ll`,`8_ll`,`9_ll`),
                            ul=c(`1_ul`,`2_ul`,`3_ul`,`4_ul`,`5_ul`,
                                 `6_ul`,`7_ul`,`8_ul`,`9_ul`))
#It's getting caught up on the one that has a missing exp cat?
