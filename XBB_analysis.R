#Description: Assessing factors that influence XBB Ct values
#Author: Nicholas Chen

#Load packages
library(tidyverse)
library(readxl)
library(xlsx)
library(snakecase)
library(plotly)
library(ggsignif)
library(cowplot)
library(patchwork)
library(lmtest)
library(ggplot2)
library(ggforce)
library(ggpubr)




# Load data and subset to only XBB lineages
data <- read_excel("C:/Users/nicho/GLab Dropbox/Nicholas Chen/Quantiative_Epi/XBB Project/GLab_SC2_sequencing_data.xlsx", 
                   col_types = c("text", "text", "numeric", 
                                 "text", "date", "text", "text", "date", 
                                 "text", "text", "numeric", "text", 
                                 "text", "text", "text", "text", "text", 
                                 "text", "text", "text", "numeric"))
data <- data %>% filter(grepl("XBB", Lineage))
names(data) <- to_any_case(names(data),case="snake")


##### CT VALUES OVER TIME #####

# Fitting a simple linear model to determine if the Ct values are changing over time 
mod1 <- lm(yale_n_1_fam~collection_date, data)
summary(mod1) # Insignificant

# Visualize 
ggplot(data, aes(x=collection_date,y=yale_n_1_fam))+
  geom_point()+
  theme_classic()+
  labs(
    x="collection date",
    y="n1"
  )

# Test for homoscedasticity and model the residuals
bptest(mod1, studentize = T) # Not significant = homosecasticity 
data$predicted <- predict(mod1)
data$residuals <- resid(mod1)
ggplot(data = data, aes(x = predicted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, col = "red") + # add horizontal line at 0
  labs(title = "Residuals vs. Fitted Values",
       x = "Fitted Values", y = "Residuals")+
  theme_classic()


ggplot(data = data, aes(x = residuals)) +
  geom_histogram(bins = 15,
                 col = "black",
                 fill = "blue",
                 alpha = 0.2,
                 closed = "left",
                 na.rm = TRUE) +
  labs(title = "Frequency Histogram of Model Residuals",
       x = "Residuals", y = "Count")+
  theme_classic()

ggplot(data, aes(sample = residuals)) +
  geom_qq() +
  geom_qq_line(col = "blue") +
  labs(title = "Normal Q-Q Plot of Model Residuals",
       x = "Theoretical Quantiles", y = "Sample Quantiles")+
  theme_classic()

# No evidence for the ct values changing over time 


###### COVARIATES #####

# Read in patient information
setwd('C:/Users/nicho/Box Sync/Yale Clinical Virology Lab Samples')
data$original_id <- as.numeric(data$original_id, scientific = F)
vir_1 <- read_excel('Grubaugh COVID Sequencing Manifest Metadata FINAL 2022.10.27.xlsx')
vir_2 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2022.11.03.xlsx")
vir_3 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2022.11.10.xlsx")
vir_4 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata DRAFT 2022.11.17.xlsx")
vir_4 <- vir_4[,c(1:20)]
vir_5 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2022.11.22.xlsx")
vir_6 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2022.12.01.xlsx")
vir_7 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2022.12.08.xlsx")
vir_8 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2022.12.15.xlsx")
vir_9 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2022.12.22.xlsx")
vir_10 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2023.01.04.xlsx")
vir_10b <- read_excel("Grubaugh 2 COVID Sequencing Manifest Metadata FINAL 2023.01.04.xlsx")
vir_11 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2023.01.11.xlsx")
vir_12 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2023.01.18.xlsx", 
                     col_types = c("text", "text", "numeric", 
                                   "text", "date", "text", "text", "text", 
                                   "date", "date", "text", "text", "text", 
                                   "text", "numeric", "numeric", "date", 
                                   "date", "date", "date"))
vir_13 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2023.01.25.xlsx",
                     col_types = c("text", "text", "numeric", 
                                   "text", "date", "text", "text", "text", 
                                   "date", "date", "text", "text", "text", 
                                   "text", "numeric", "numeric", "date", 
                                   "date", "date", "date"))
vir_14 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2023.02.01.xlsx",
                     col_types = c("text", "text", "numeric", 
                                   "text", "date", "text", "text", "text", 
                                   "date", "date", "text", "text", "text", 
                                   "text", "numeric", "numeric", "date", 
                                   "date", "date", "date"))
vir_15 <- read_excel("Grubaugh COVID Sequencing Manifest Metadata FINAL 2023.02.08.xlsx",
                     col_types = c("text", "text", "numeric", 
                                   "text", "date", "text", "text", "text", 
                                   "date", "date", "text", "text", "text", 
                                   "text", "numeric", "numeric", "date", 
                                   "date", "date", "date"))
virology <- rbind(vir_1,vir_2,vir_3,vir_4,vir_5,vir_6,vir_7,vir_8,vir_9,vir_10, vir_10b, vir_11,vir_12,vir_13,vir_14,vir_15) %>% 
  select(MPI,Gender,`Collection Date`,`Patient Class`,days_since_latest_vax,DOB)


data$MPI <- data$original_id
merge <- merge(data, virology, by="MPI") %>% select(-original_id)


# Covariate: SEX
gender <- merge %>% select(yale_n_1_fam, Gender) %>% filter(Gender %in% c("M","F"))
male <- merge %>% select(yale_n_1_fam, Gender) %>% filter(Gender == "M")
female <- merge %>% select(yale_n_1_fam, Gender) %>% filter(Gender == "F")

# Assess normalcy
shapiro.test(male$yale_n_1_fam) #fail to reject = normal 
shapiro.test(female$yale_n_1_fam) #reject = not normal 

# Assess variance
var.test(male$yale_n_1_fam,female$yale_n_1_fam) #equal variances 

# Assumption 2 of t-test is met, assumption 1 is circumvented on account of n>30
t.test(male$yale_n_1_fam,female$yale_n_1_fam, var.equal = T) #no difference

# Plot
ggboxplot(gender, x = "Gender", y = 'yale_n_1_fam', 
          color = "Gender", palette = c("#00AFBB", "#E7B800"),
          ylab = "Yale-N1(FAM)", xlab = "Gender")

merge <- merge %>% 
  mutate(sex = case_when(
    startsWith(Gender,"F") ~ "Female",
    startsWith(Gender,"M") ~ "Male",
    TRUE ~ "Other"
  )) 
merge <- merge %>% filter(sex !="Other")


xbb_sex <- ggplot(merge,aes(x=sex,y=yale_n_1_fam, color=sex))+
  geom_hline(yintercept = 25.2,colour = "lightgrey")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Ct Value by Sex",x=NULL, y = "Ct (N)")+
  scale_color_manual(values=c("#00A08A","#FF0000"))+
  geom_signif(comparisons = list(c("Female","Male")),map_signif_level=T,color="black")+ 
  coord_fixed(ratio = 1/14.1, clip = 'off', expand = TRUE, ylim = c(40,10))+
  scale_y_reverse(breaks = seq(10,40,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=10, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  #stat_summary(fun.data = median.n, geom = "label", colour = "black",
  #             size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(sex) %>% count(), 
            aes(label = paste0("n=",n), x = sex), y = -45, size = 3, fontface = c("italic"),
            check_overlap = TRUE,color="black")

setwd("C:/Users/nicho/Documents/GitHub/Omicron_Project/")
#ggsave(filename="xbb_sex.png", units = "cm")


# Covariate: PATIENT CLASS 
patient_class <- merge %>% select(yale_n_1_fam, `Patient Class`) %>% rename(class = `Patient Class`)

# Assess normalcy
with(patient_class, shapiro.test(yale_n_1_fam[class == "Outpatient"])) #not normal
with(patient_class, shapiro.test(yale_n_1_fam[class == "Inpatient"])) #normal 
with(patient_class, shapiro.test(yale_n_1_fam[class == "Emergency"])) #normal

# Assess variance
bartlett.test(yale_n_1_fam ~ class, data = patient_class) #unequal variance

#Approximately normal but with unequal variances = Welch ANOVA
oneway.test(yale_n_1_fam ~ class, data = patient_class, var.equal = FALSE) # Significant 
anova <- aov(yale_n_1_fam ~ class, data = patient_class) # Confirm result with standard ANOVA
summary(anova) # Concordance

# Post hoc pairwise comparisons 
TukeyHSD(anova, conf.level=.95) # Significant difference between out-patient and emergency

patient_class$class <- factor(patient_class$class, levels= c("Inpatient","Outpatient","Emergency"))

xbb_class <- ggplot(patient_class,aes(x=class,y=yale_n_1_fam, color=class))+
  geom_hline(yintercept = 25.2,colour = "lightgrey")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Ct Value by Patient Class",x=NULL, y = "Ct (N)")+
  scale_color_manual(values=c("#46ACC8", "#E58601", "#B40F20"))+
  geom_signif(comparisons = list(c("Outpatient","Emergency")),map_signif_level=T,color="black")+ 
  coord_fixed(ratio = 1/9.7, clip = 'off', expand = TRUE, ylim = c(40,10))+
  scale_y_reverse(breaks = seq(10,40,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=10, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  #stat_summary(fun.data = median.n, geom = "label", colour = "black",
  #             size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(class) %>% count(), 
            aes(label = paste0("n=",n), x = class), y = -45, size = 3, fontface = c("italic"),
            check_overlap = TRUE,color="black")



# Covariate: AGE
# Approximate age
merge$age <- round(as.numeric(merge$collection_date-merge$DOB)/365.25,digits=0)

# Age categories
merge <- merge %>% mutate(agec = 
                            if_else(age>=0 & age<18,1,
                                    if_else(age>=18 & age<50,2,
                                            if_else(age>=50 & age<70,3,
                                                    if_else(age>=70,4,9999)))))
merge$age_cat <- factor(merge$agec, levels= c(1,2,3,4),labels=c("<18","18-49","50-69","70+"))
merge<- merge%>% filter(!(is.na(age_cat)))

# Assess normalcy
with(merge, shapiro.test(yale_n_1_fam[age_cat == "<18"])) #normal
with(merge, shapiro.test(yale_n_1_fam[age_cat == "18-49"])) #normal 
with(merge, shapiro.test(yale_n_1_fam[age_cat == "50-69"])) #normal 
with(merge, shapiro.test(yale_n_1_fam[age_cat == "70+"])) #normal

# Assess variance
bartlett.test(yale_n_1_fam ~ age_cat, data = merge) #equal variances 

# ANOVA
anova_2 <- aov(yale_n_1_fam ~ age_cat, data = merge)
summary(anova_2) #Significant

# Pairwise post-hoc 
TukeyHSD(anova_2, conf.level=.95)
# Significant difference between 70+ and 18-49


xbb_age <- ggplot(merge,aes(x=age_cat,y=yale_n_1_fam, color=age_cat))+
  geom_hline(yintercept = 25.2,colour = "lightgrey")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Ct Value by Age Category",x=NULL, y = "Ct (N)")+
  scale_color_manual(values=c("#0B775E", "#354823" ,"#F2300F","#273046"))+
  geom_signif(comparisons = list(c("70+","18-49")),map_signif_level=T,color="black")+ 
  coord_fixed(ratio = 1/7.4, clip = 'off', expand = TRUE, ylim = c(40,10))+
  scale_y_reverse(breaks = seq(10,40,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=10, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  #stat_summary(fun.data = median.n, geom = "label", colour = "black",
  #             size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(age_cat) %>% count(), 
            aes(label = paste0("n=",n), x = age_cat), y = -45, size = 3, fontface = c("italic"),
            check_overlap = TRUE,color="black")




# Covariate = VACCINATION STATUS
# Read in vaccination data and merge with Ct value data 
vax_data <- read.csv("C:/Users/nicho/Box Sync/Yale Virology sample intake/vaccination_metadata_database/vaccine_metadata_clean.csv")
vax_data <- vax_data %>% select(-c(X,patient_class,country,city,zip,tube_label,original_sample_type,report_ct,extraction_pcr_date,yale_n1_fam:yale_n1_ge_ml, collection_date_gsheet)) %>% 
  rename(collection_date = collection_date_virology, MPI = mpi)
vax_data[,c("dob","collection_date","vax_1_date","vax_2_date","booster_1_date","booster_2_date","series_complete_date")] <- 
  lapply(vax_data[,c("dob","collection_date","vax_1_date","vax_2_date","booster_1_date","booster_2_date","series_complete_date")], as.Date,format="%m/%d/%Y")

# Apply filters (Clinical virology, Connecticut residents, non-missing MPI/MRN/Lineage info/dob, first instance of duplicate MRNs, consistent vaccination records 
vax_data <- vax_data %>% filter(source == "Yale Clinical Virology Lab" & state == "Connecticut" & filter == "connecticut" 
                                & metadata == T & inconsistent_vax_info == F & dup_mpi_gsheet==F & miscoded_mpi == F & !is.na(lineage) 
                                & !is.na(dob) & duplicate_mrn == F) %>% filter(!duplicated(mrn)) 


xbb_vax_data <- merge(merge[,c("MPI","yale_n_1_fam")], vax_data, by="MPI")


# Create interval column 
xbb_vax_data$Interval <- ""

# Categorize vaccination status 
xbb_vax_data <- xbb_vax_data %>% 
  mutate(vax_status = case_when(
    difftime(collection_date, vax_1_date, units = "days") < 14 ~ "Non_Breakthrough",
    difftime(collection_date, booster_2_date, units = "days") >= 14 ~ "Four_Dose_Breakthrough",
    difftime(collection_date, booster_1_date, units = "days") >= 14 ~ "Three_Dose_Breakthrough",
    difftime(collection_date, vax_2_date, units = "days") >= 14 ~ "Two_Dose_Breakthrough",
    difftime(collection_date, vax_1_date, units = "days") >= 14 ~ "One_Dose_Breakthrough",
    TRUE ~ "Non_Breakthrough"
  ))


# Calculate time from most recent vaccination to collection date
xbb_vax_data$time_since_vax <- ""
for(i in 1:nrow(xbb_vax_data)){
  if(xbb_vax_data[i,"vax_status"]=='One_Dose_Breakthrough'){
    xbb_vax_data[i,"time_since_vax"] <- difftime(xbb_vax_data[i,"collection_date"],xbb_vax_data[i,"vax_1_date"],units='days')
  } else if(xbb_vax_data[i,"vax_status"]=='Two_Dose_Breakthrough'){
    xbb_vax_data[i,"time_since_vax"] <- difftime(xbb_vax_data[i,"collection_date"],xbb_vax_data[i,"vax_2_date"],units='days')
  } else if(xbb_vax_data[i,"vax_status"]=='Three_Dose_Breakthrough'){
    xbb_vax_data[i,"time_since_vax"] <- difftime(xbb_vax_data[i,"collection_date"],xbb_vax_data[i,"booster_1_date"],units='days')
  } else if(xbb_vax_data[i,"vax_status"]=='Four_Dose_Breakthrough'){
    xbb_vax_data[i,"time_since_vax"] <- difftime(xbb_vax_data[i,"collection_date"],xbb_vax_data[i,"booster_2_date"],units='days')
  } else { 
    NULL
  }
}

xbb_vax_data$time_since_vax <- as.numeric(xbb_vax_data$time_since_vax)

# Update vaccination statuses
for(i in 1:nrow(xbb_vax_data)){
  if(xbb_vax_data[i,"vax_status"] == "Non_Breakthrough"){
    xbb_vax_data[i,"vax_status_revised"]  <- "0 Doses"
  } else if(xbb_vax_data[i,"vax_status"] == "One_Dose_Breakthrough"){
    xbb_vax_data[i,"vax_status_revised"]  <- "1 Dose"
  } else if(xbb_vax_data[i,"vax_status"] %in% c("Two_Dose_Breakthrough","Three_Dose_Breakthrough","Four_Dose_Breakthrough") & xbb_vax_data[i,"time_since_vax"] <150){
    xbb_vax_data[i,"vax_status_revised"]  <- "2 Doses <5 Months"
  } else if(xbb_vax_data[i,"vax_status"] %in% c("Two_Dose_Breakthrough","Three_Dose_Breakthrough","Four_Dose_Breakthrough") & xbb_vax_data[i,"time_since_vax"] >=150){
    xbb_vax_data[i,"vax_status_revised"]  <- "2 Doses >=5 Months"
  }
}


# Assess normalcy
with(xbb_vax_data, shapiro.test(yale_n_1_fam[vax_status_revised == "0 Doses"])) #normal
with(xbb_vax_data, shapiro.test(yale_n_1_fam[vax_status_revised == "1 Dose"])) #normal
with(xbb_vax_data, shapiro.test(yale_n_1_fam[vax_status_revised == "2 Doses <5 Months"])) #normal
with(xbb_vax_data, shapiro.test(yale_n_1_fam[vax_status_revised == "2 Doses >=5 Months"])) #normal


# Assess variance
bartlett.test(yale_n_1_fam ~ vax_status_revised, data = xbb_vax_data) #equal variances 

# ANOVA
anova_3 <- aov(yale_n_1_fam ~ vax_status_revised, data = xbb_vax_data)
summary(anova_3) #Not significant

# Plotting
xbb_vax <- ggplot(xbb_vax_data,aes(x=vax_status_revised,y=yale_n_1_fam, color=vax_status_revised))+
  geom_hline(yintercept = 25.2,colour = "lightgrey")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Ct Value by Vaccination Status",x=NULL, y = "Ct (N)")+
  scale_color_manual(values=c("#046C9A", "#899DA4" , "black","#02401B"))+
  geom_signif(comparisons = list(c("0 Doses","2 Doses >=5 Months")),map_signif_level=T,color="black")+ 
  coord_fixed(ratio = 1/7.4, clip = 'off', expand = TRUE, ylim = c(40,10))+
  scale_y_reverse(breaks = seq(10,40,by = 5))+
  scale_x_discrete(labels = c("0 Doses","1 Dose", "2 Doses \n<5 Months", "2 Doses \n>=5 Months"))+
  theme_classic()+theme(legend.position="none",axis.text.y = element_text(size = 10, colour = "black",face="bold"),
                        axis.text.x.bottom = element_text(size = 8, colour = "black",face="bold"),title=element_text(size=10, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  #stat_summary(fun.data = median.n, geom = "label", colour = "black",
  #             size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(vax_status_revised) %>% count(), 
            aes(label = paste0("n=",n), x = vax_status_revised), y = -46, size = 3, fontface = c("italic"),
            check_overlap = TRUE,color="black")

plot_grid(xbb_sex,xbb_class+labs(y=""),xbb_age+labs(y=""),xbb_vax+labs(y=""),nrow=1)

ggsave(filename="xbb_wrap.png", height=10,width=32,units = "cm")
