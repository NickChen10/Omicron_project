#Description: Vaccination Coverage by Doses in CT 
#Author: Nicholas Chen


#Load required packages
library(tidyverse)
library(ggplot2)
library(viridis)



##### DATA FORMATTING #####
# https://data.cdc.gov/Vaccinations/COVID-19-Vaccination-Trends-in-the-United-States-N/rh2h-3yt2
path=getwd()

# Read in Data and filter to CT  
cdc <- read.csv(paste0(path,"/GitHub/Omicron_Project/Data/COVID-19_Vaccination_Trends_in_the_United_States_National_and_Jurisdictional.csv"))
cdc$Date <- as.Date(cdc$Date,"%m/%d/%Y")

cdc <- cdc %>% filter(Location == "CT" & date_type=="Admin") %>% select(Date,Administered_Dose1_Pop_Pct,Series_Complete_Pop_Pct, Second_Booster_50Plus_Vax_Pct, Additional_Doses_Vax_Pct, Bivalent_Booster_Pop_Pct) %>%
  rename("One Dose" = "Administered_Dose1_Pop_Pct", "Two Doses" = "Series_Complete_Pop_Pct", "Three Doses" = "Additional_Doses_Vax_Pct" ,"Four Doses" = "Second_Booster_50Plus_Vax_Pct", "Bivalent" = "Bivalent_Booster_Pop_Pct") %>% gather("Vaccination Status","Percent Coverage",2:6)

cdc$`Vaccination Status` <-factor(cdc$`Vaccination Status`, levels = c("One Dose","Two Doses","Three Doses","Four Doses", "Bivalent"))


##### PLOTTING #####

ggplot(data = cdc, aes(x=Date, y=`Percent Coverage`,color=`Vaccination Status`,fill=`Vaccination Status`)) + 
  geom_area(position = "identity")+
  geom_line(size=1)+
  scale_x_date(breaks="1 month", expand=c(0,0))+
  geom_hline(yintercept=95,linetype="dashed", size=1)+
  scale_y_continuous(expand = c(0,0), limits = c(0,100))+
  theme_classic()+
  ylab("Population Coverage (%)")+
  xlab("")+
  ggtitle("Vaccination Coverage in Connecticut")+
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 11, angle=60, hjust=1,face = "bold"),
    axis.text.y = element_text(size = 12,hjust = 0.5,face = "bold"),
    axis.title.y = element_text(face="bold"),
    legend.box.background = element_rect(colour = NA),
    legend.title = element_text(face="bold"))+
  scale_color_viridis_d(direction=-1)+
  scale_fill_viridis_d(direction=-1)#+
  #scale_alpha_manual(values = c(0.3,0.4,0.7,0.7,0.8))

ggsave(filename="vax_trends_CT_2.png",width=24,height=10,units = "cm")
