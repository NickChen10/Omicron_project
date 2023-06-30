#Description: Weekly variant frequencies in CT
#Author: Nicholas Chen

#Load packages
library(tidyverse)
library(ggplot2)



##### DATA CLEANING AND FORMATTING #####
path=getwd()

# Read weekly frequency data from CovidTracker.CT
data <- read_csv(paste0(path,"/GitHub/Omicron_Project/Data/SARS-CoV-2 variant frequency in Connecticut.csv"))

data <- data %>%
  mutate(ShortLin = case_when(
    endsWith(`Variant names`, "Alpha") ~ "Alpha",
    endsWith(`Variant names`, "Delta") ~ "Delta",
    startsWith(`Variant names`, "Omicron (BA.1)") ~ "BA.1",
    startsWith(`Variant names`, "Omicron (BA.2)") ~ "BA.2",
    startsWith(`Variant names`, "Omicron (BA.4)") ~ "BA.4",
    endsWith(`Variant names`, "Other") ~ "Others",
    startsWith(`Variant names`, "Omicron (BA.5)") ~ "BA.5",
    startsWith(`Variant names`, "Omicron (XBB)") ~ "XBB",
    endsWith(`Variant names`, "") ~ "Others"
  ))

data <- data %>% group_by(ShortLin,week)%>% summarise(freq=sum(freq))
data$freq[which(data$freq <=0.02)] <-0


##### PLOTTING #####
colors <- c( "Alpha" = "#0047ab",
             "Delta" = "#ff6e40",
             "BA.1" = "#184E27", 
             "BA.2" = "#CE4E50",
             "BA.4" = "#CAC6EF",
             "BA.5" = "#7294D4",
             "XBB" = "#892F2E",
             "Others" = "#515151"
)

ggplot(data,aes(x=`week`,y=freq*100, group = `ShortLin`, color = `ShortLin`)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Weekly variant frequencies in Connecticut",
    y = "Variant frequency (%)",
    x = "") + 
  scale_color_manual(values = colors, breaks = levels(colors)) +
  scale_y_continuous(expand=c(0,0))+
  scale_x_date(expand=c(0,0),date_breaks="3 month", date_labels="%b-%Y") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 15,hjust = -0.2, face = "bold"),
    axis.text.x = element_text(size = 11,angle = 45, vjust = 1, hjust=1, face = "bold"),
    axis.text.y = element_text(size = 12,hjust = 0.5)
  )

#ggsave(filename="frequenciesCT.png", units = "cm")
