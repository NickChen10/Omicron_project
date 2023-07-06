#Description: Satewide and Yale specific variant trends
#Author: Nicholas Chen


#Load required packages
library(tidyverse)
library(ggplot2)
library(readxl)
library(snakecase)



##### DATA FORMATTING #####
path=getwd()

# Read weekly frequency data from CovidTracker.CT
data <- read_csv(paste0(path,"/GitHub/Omicron_Project/Data/SARS-CoV-2 variant frequency in Connecticut.csv"))

# Format dates and recode lineages
data <- data %>% filter(between(
  week, as.Date('2022-01-01'), 
  as.Date('2023-01-30'))) 
data <- data %>% mutate(ShortLin = `Variant names`) %>%
  mutate(ShortLin = case_when(
    endsWith(ShortLin, "Delta") ~ "Others",
    startsWith(ShortLin, "Omicron (BA.1)") ~ "BA.1",
    startsWith(ShortLin, "Omicron (BA.2)") ~ "BA.2",
    startsWith(ShortLin, "Omicron (BA.4)") ~ "BA.4",
    endsWith(ShortLin, "Other") ~ "Others",
    startsWith(ShortLin, "Omicron (BA.5)") ~ "BA.5",
    startsWith(ShortLin, "Omicron (XBB)") ~ "XBB",
    endsWith(ShortLin, "") ~ "Others"
  ))
data$ShortLin[is.na(data$ShortLin)]<- "Others"


# Read data from Yale SARS-CoV-2 Genomic Surveillance Initiative
yale_data <-  read_excel("C:/Users/nicho/GLab Dropbox/Nicholas Chen/Quantiative_Epi/XBB Project/GLab_SC2_sequencing_data.xlsx", 
                    col_types = c("text", "text", "numeric", 
                                  "text", "date", "text", "text", "date", 
                                  "text", "text", "numeric", "text", 
                                  "text", "text", "text", "text", "text", 
                                  "text", "text", "text", "numeric"))
# Format names and dates 
names(yale_data) <- to_any_case(names(yale_data),case="snake")
yale_data[,c("collection_date","extraction_pcr_date")] <- lapply(yale_data[,c("collection_date","extraction_pcr_date")],as.Date,format="%m/%d/%Y")

# Filter to only sequenced data and dates of interest
data_clean <- data_clean %>% filter(grepl('connecticut', Filter)) %>% 
  filter(between(collection, 
                 as.Date('2022-01-01'), as.Date('2023-01-30'))) 

# Recode lineages
yale_data <- yale_data %>%
  mutate(ShortLin = case_when(
    endsWith(lineage, "2.12.1") ~ "BA.2.12.1",
    startsWith(lineage, "Q") ~ "Other",
    startsWith(lineage, "B.1.1.7") ~ "Other",
    startsWith(lineage, "AY") ~ "Other",
    startsWith(lineage, "B.1.617") ~ "Other",
    startsWith(lineage, "BG") ~ "BA.2",
    startsWith(lineage, "BA.2.75") ~ "BA.2",
    startsWith(lineage, "BA.1.1") ~ "BA.1.1",
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
    startsWith(lineage, "BA.4.6") ~ "BA.4.6",
    startsWith(lineage, "BA.4") ~ "BA.4",
    startsWith(lineage, "BA.5") ~ "BA.5",
    startsWith(lineage, "BU") ~ "BA.5",
    startsWith(lineage, "BV") ~ "BA.5",
    startsWith(lineage, "BE") ~ "BA.5",
    startsWith(lineage, "BT") ~ "BA.5",
    startsWith(lineage, "BF.7") ~ "BA.5",
    startsWith(lineage, "BF") ~ "BA.5", #BF.7 as BA.5
    startsWith(lineage, "BQ.1.1") ~ "BQ.1.1",  
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
    startsWith(lineage, "XBB.1.5") ~ "XBB.1.5",
    startsWith(lineage, "XBB.1") ~ "XBB.1", 
    TRUE ~ "Other"
  ))
yale_data <- yale_data %>% filter(ShortLin != "Other")

##### PLOTTING #####
#Weekly variants across CT
OmiCols <- c("BA.1" = "#184E27", 
             "BA.1.1"= "#74875B",
             "BA.2" = "#CE4E50",
             "BA.2.12.1" ="#FB7069" ,
             "BA.4" = "#CAC6EF",
             "BA.4.6" = "#CAC6EF",
             "BA.5" = "#7294D4",
             "BQ.1.1" = "#184294",
             "XBB" = "#892F2E",
             "XBB.1" = "#892F2E",
             "XBB.1.5" = "#892F2E",
             "Others" = "grey"
)

ggplot(data, aes(x = as.Date(`week`), y = `n`, fill = as.factor(`ShortLin`), alluvium = as.factor(`ShortLin`))) + 
  geom_col(width = 5, show.legend = FALSE) +
  theme_classic() + 
  scale_x_date(expand=c(0,0)) + 
  scale_y_continuous(expan=c(0,0))+
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black",face ="bold"),
        legend.title = element_blank(),
        plot.title = element_text(size = 14, color = "black",face ="bold"))+
  ggtitle("Samples per week Connecticut") + 
  ylab("number of samples") +
  scale_fill_manual(values=c(OmiCols)) +
  scale_color_manual(values=c(OmiCols))+
  coord_fixed(ratio = 0.075, xlim = c(as.Date("2022-01-01"),as.Date("2023-01-20")))

ggsave(filename="cases_CT.png", units = "cm")



# Weekly variants sequenced at Yale
# Assign week
lineages <- yale_data %>% mutate(CopyDate = collection_date)
lineages <- lineages %>% group_by(`week` = cut(`CopyDate`, "week"))
lineages$week <- as.Date(lineages$week)
lineages <- lineages[!is.na(lineages$`week`),] #delete data without week assigned

# Aggregate to weekly counts
lineages<- lineages %>% group_by(week) %>% count(ShortLin)


ggplot(lineages, aes(x = as.Date(`week`), y = `n`, fill = as.factor(`ShortLin`), alluvium = as.factor(`ShortLin`))) + 
  geom_col(width = 5) +
  theme_classic() + 
  scale_x_date(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0))+
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, 
                                   size = 10, color = "black",face ="bold"),
        legend.title = element_blank(),
        plot.title = element_text(size = 14, color = "black",face ="bold"))+
  ggtitle("Samples per week Yale") + 
  ylab("number of samples") +
  scale_fill_manual(values=c(OmiCols)) +
  scale_color_manual(values=c(OmiCols)) +
  coord_fixed(ratio = 0.4, xlim = c(as.Date("2022-01-01"),as.Date("2023-01-20")))

#ggsave(filename="yale_cases.png", units = "cm", scale = 1.2)
