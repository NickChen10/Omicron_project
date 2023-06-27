#Description: Fitting loess regression and logistic growth rates to periods of variant emergence for all of Connecticut.
#Author: Nicholas Chen 

#Load packages


packages = c("dplyr","tidyr","ggplot2","RColorBrewer",
             "plotly","DT","data.table",
             "stats","ggpubr","stringr")

package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

rm(packages, package.check)

library(purrr)
library(readr)
library(DescTools)
library(lubridate)
library(zoo)
library(EpiEstim)
library(ggpubr)
library(httr)
library(jsonlite)
library(cowplot)
library(forcats)




##### DATA CLEANING AND FORMATTING #####

path=getwd()

# Read GISAID data for daily lineage frequencies, format date, remove those with missing lineages or duplicated accession IDs
metadata_CT_all <- read.table(file=paste0(path,"/data/GISAID_metadata.tsv"), sep = "\t", header = TRUE) 
metadata_CT_all$`Collection.date` <- as.Date(metadata_CT_all$`Collection.date`)
metadata_CT_all <- metadata_CT_all %>%
  filter(Lineage != "") %>%
  drop_na(Collection.date)%>% filter(!duplicated(`Accession.ID`))



# Filter badly named lineages
metadata_CT_all <- metadata_CT_all %>% rename(pango_lineage = Lineage) 
metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2.75 (marker override based on Emerging Variants AA substitutions)"] = "BA.2.75"
metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BF.6 (marker override based on Emerging Variants AA substitutions)"] = "BF.6"
metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "XBB (marker override based on Emerging Variants AA substitutions)"] = "XBB"
metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2.13.1 (marker override based on Emerging Variants AA substitutions)"] = "BA.2.13.1"
metadata_CT_all["pango_lineage"][metadata_CT_all["pango_lineage"] == "BA.2" & metadata_CT_all["Collection.date"] > "2022-12-16"] = "XBB.1.5"




# Group by week
metadata_CT_all <- metadata_CT_all %>% group_by(`week` = cut(`Collection.date`, "week")) 
metadata_CT_all$week <- as.Date(metadata_CT_all$week)
metadata_CT_all <- metadata_CT_all[!is.na(metadata_CT_all$`week`),]



# Read Pangolin lineage designation and match to WHO variant designations
names <- read.csv(paste0(path,"/data/Pangolin_variant_names.csv"))
metadata_CT_pango <- merge(metadata_CT_all,names[,c("pango_lineage","who_variants")], by = "pango_lineage", all.x = TRUE)
metadata_CT_pango$who_variants[which(is.na(metadata_CT_pango$who_variants))] <- "Other" #replace NA's with other
sum(is.na(metadata_CT_pango$pango_lineage)) #no missing lineages


# Create dummy column for lineage counting
metadata_CT_pango$`CopyDate` <- as.numeric(metadata_CT_pango$`Collection.date`)

# Remove most recent 3 collection dates (inadequate data)
metadata_CT_pango <- metadata_CT_pango[-c(which(metadata_CT_pango$`Collection.date` >= (max(metadata_CT_pango$`Collection.date`)-2))),] # omit latest 3 days due to reporting bias

# Have all lineages be displayed every day, ensuring no breaks when plotting
metadata_CT_pango <- metadata_CT_pango %>% complete(`Collection.date`, nesting(`pango_lineage`), fill = list(CopyDate = 0))
metadata_CT_pango$`CopyDate` <- ifelse(metadata_CT_pango$`CopyDate` == 0,0,1) # for counting lineages


# Reclassify lineage variant designations
metadata_CT_pango <- metadata_CT_pango %>%
  mutate(ShortLin = case_when(
    endsWith(pango_lineage, "2.12.1") ~ "BA.2", #BA.2.12.1 merged with BA.2 
    startsWith(pango_lineage, "BG") ~ "BA.2",
    startsWith(pango_lineage, "BA.2.75") ~ "BA.2",
    startsWith(pango_lineage, "BA.1.1") ~ "BA.1", #Merging BA.1.1 with BA.1 
    startsWith(pango_lineage, "BA.1") ~ "BA.1",
    startsWith(pango_lineage, "BA.2") ~ "BA.2",
    startsWith(pango_lineage, "BK") ~ "BA.2",
    startsWith(pango_lineage, "BP") ~ "BA.2",
    startsWith(pango_lineage, "BL") ~ "BA.2",
    startsWith(pango_lineage, "BH") ~ "BA.2",
    startsWith(pango_lineage, "BM") ~ "BA.2",
    startsWith(pango_lineage, "CH") ~ "BA.2",
    startsWith(pango_lineage, "BJ") ~ "BA.2",
    startsWith(pango_lineage, "BR") ~ "BA.2",
    startsWith(pango_lineage, "BA.4.6") ~ "BA.4", #Merging BA.4.6 with BA.4
    startsWith(pango_lineage, "BA.4") ~ "BA.4",
    startsWith(pango_lineage, "BA.5") ~ "BA.5",
    startsWith(pango_lineage, "BU") ~ "BA.5",
    startsWith(pango_lineage, "BV") ~ "BA.5",
    startsWith(pango_lineage, "BE") ~ "BA.5",
    startsWith(pango_lineage, "BT") ~ "BA.5",
    startsWith(pango_lineage, "BF.7") ~ "BA.5",
    startsWith(pango_lineage, "BF") ~ "BA.5",
    startsWith(pango_lineage, "BQ.1.1") ~ "BA.5", #BQ merged into BA.5 
    startsWith(pango_lineage, "BQ.2") ~ "BA.5",
    startsWith(pango_lineage, "BQ.3") ~ "BA.5",
    startsWith(pango_lineage, "CN") ~ "BA.5",
    startsWith(pango_lineage, "CD") ~ "BA.5",
    startsWith(pango_lineage, "CE") ~ "BA.5",
    startsWith(pango_lineage, "CF") ~ "BA.5",
    startsWith(pango_lineage, "CL") ~ "BA.5",
    startsWith(pango_lineage, "CG") ~ "BA.5",
    startsWith(pango_lineage, "CK") ~ "BA.5",
    startsWith(pango_lineage, "BZ") ~ "BA.5",
    startsWith(pango_lineage, "XBB.1.5") ~ "XBB.1",
    startsWith(pango_lineage, "XBB.1") ~ "XBB.1", 
    TRUE ~ "Other"
  ))


# Daily lineage proportions
daily_lineages <- metadata_CT_pango %>% group_by(`Collection.date`,`ShortLin`) %>% 
  summarise_at(vars(`CopyDate`),list(n = sum)) %>%
  mutate(freq = round(100*n/sum(n),2))

rm(names, metadata_CT_all)



##### EMERGENCE PERIOD LOESS REGRESSIONS #####


variants <- c("BA.1","BA.2","BA.4","BA.5","XBB.1") #variants of interest

# Note: After plotting I noticed some outlier daily frequencies (single days that are much higher than those around them) 
daily_lineages$freq[which(daily_lineages$Collection.date == "2022-11-23" | daily_lineages$Collection.date == "2022-11-24" & daily_lineages$ShortLin == "BA.5")] <- (daily_lineages$freq[which(daily_lineages$Collection.date == "2022-11-22" & daily_lineages$ShortLin == "BA.5")]+daily_lineages$freq[which(daily_lineages$Collection.date == "2022-11-24" & daily_lineages$ShortLin == "BA.5")])/2
  # Averaging the outlier dates with the two days on either side


# Define date ranges for each variant(first instance of 5% daily frequency to the first instance of the maximum daily frequency)
fun <- function(x){
  interval<- daily_lineages %>% filter(ShortLin == x & freq >= 5)
  max <- min(which(interval$freq == max(interval$freq))) 
  interval <- interval[c(1:max),]
  x <- range(interval$Collection.date)
  
}

for(x in 1:length(variants)){
  assign(paste(variants[x],"interval", sep="_"),fun(variants[x]))
}

# To plot BA.4 and BA.5 in the same pane, set BA.4's emergence periods to match BA.5's
difftime(max(BA.5_interval),min(BA.5_interval)) #96 days 
BA.4_interval <- c(min(BA.4_interval),min(BA.4_interval)+96)

##### LOESS REGRESSION FITTING #####

# Define variants of interest and their associated colors
intervals <- list(BA.1_interval,BA.2_interval,XBB.1_interval)
colors <- c("BA.1" = "#184e27", 
            "BA.2" = "#CE4E50",
            "XBB.1" = "#892F2E",
            "Other" = "grey")
variants <- c("BA.1","BA.2","XBB.1")

# Model BA.1, BA.2, and XBB.1
for(j in 1:length(variants)){
  model_data <- daily_lineages %>% 
    filter(Collection.date >= min(intervals[[j]]) & Collection.date <= max(intervals[[j]])) #Restrict collection date to interval dates
  
  model_data$Collection.date <-as.numeric(model_data$Collection.date - min(model_data$Collection.date)) #Start time series at 0
  model_data$ShortLin <- ifelse(model_data$ShortLin == variants[j], variants[j], "Other") #If not the variant of interest, set to "other"
  model_data$ShortLin <- factor(model_data$ShortLin, levels=c("Other",paste0(variants[j]))) #Factor variants 
  model_data <- model_data %>% group_by(Collection.date, ShortLin) %>% summarize(freq=sum(freq)) #Recalculate daily variant frequencies
  
  p <- ggplot(model_data, aes(x=as.numeric(Collection.date), y=freq, col=ShortLin, fill=ShortLin, group=ShortLin))+ 
    geom_point()+
    geom_smooth(size=2,alpha=0.2)+ #Loess regression 
    scale_color_manual(values=c("#adadad",paste0(colors[j])))+
    scale_fill_manual(values=c("#adadad",paste0(colors[j])))+
    ylab("Daily Frequency")+
    xlab("Days Since 5% Emergence")+
    ggtitle(paste0(variants[j]," Emergence"))+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(limits=c(0,100),expand=c(0,0))+
    theme_bw()+
    theme(legend.position="none")
  
  p_data <- ggplot_build(p)$data[[2]] %>% filter(colour != "#adadad") %>% select(x,y) %>% arrange(x)  %>%  filter(y >= 50) 
    #Pull data from the fitted regression and select the first day the fitted line passes 50% (can be changed) 
  p_min <- min(p_data$x) #First day the 50% threshold is passed 
  p_min <- ifelse(is.infinite(p_min), 0,p_min) #hides the dashed line for plots that don't cross the set threshold
  p_2 <- p+geom_vline(xintercept = p_min, linetype="dashed")
  assign(paste0(variants[j],"_plot"),p_2)
  assign(paste0(variants[j],"_emergence"),min(p_data$x))
}

# Model BA.4 and 5 together (same code as above but modified to allow for two variants)
BA.4_5_plot_data <- daily_lineages %>% 
  filter(Collection.date >= min(BA.5_interval) & Collection.date <= max(BA.5_interval))

BA.4_5_plot_data$Collection.date <-as.numeric(BA.4_5_plot_data$Collection.date - min(BA.4_5_plot_data$Collection.date)) 
BA.4_5_plot_data$ShortLin <- ifelse(BA.4_5_plot_data$ShortLin %in% c("BA.4","BA.5"),BA.4_5_plot_data$ShortLin, "Other")
BA.4_5_plot_data$ShortLin <-factor(BA.4_5_plot_data$ShortLin, levels=c("Other","BA.4","BA.5"))
BA.4_5_plot_data <- BA.4_5_plot_data %>% group_by(Collection.date, ShortLin) %>% summarize(freq=sum(freq))

BA.4.5_plot <- ggplot(BA.4_5_plot_data, aes(x=as.numeric(Collection.date), y=freq, col=ShortLin, fill=ShortLin, group=ShortLin))+ 
  geom_point()+
  geom_smooth(size=2, alpha=0.2)+
  scale_color_manual(values=c("#adadad","#CAC6EF","#7294D4"))+
  scale_fill_manual(values=c("#adadad","#CAC6EF","#7294D4"))+
  ylab("Daily Frequency")+
  xlab("Days Since 5% Emergence")+
  ggtitle("BA.5 Emergence")+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+
  theme_bw()+
  theme(legend.position="none")

BA.4.5_plot_extract <- ggplot_build(BA.4.5_plot)$data[[2]] %>% filter(colour == "#7294D4")%>% select(x,y) %>% arrange(x)  %>%  filter(y >= 50) 
BA.4.5_plot_min <- min(BA.4.5_plot_extract$x) 
BA.4.5_plot_min <- ifelse(is.infinite(BA.4.5_plot_min), 0,BA.4.5_plot_min)
BA.4.5_plot <- BA.4.5_plot+geom_vline(xintercept = BA.4.5_plot_min, linetype="dashed")
BA.4.5_emergence <- BA.4.5_plot_min

rm(BA.4_5_plot_data,BA.4.5_plot_min,BA.4.5_plot_extract,model_data,p,p_2,p_data,j,x,p_min)



# Plot together with emergence period text boxes
arranged_plot<-
  plot_grid(BA.1_plot+geom_label(color="black",fill="white", size=3.5,aes(x=10,y=96,label=paste("Emergence Period (5-50%):",round(BA.1_emergence,digits = 0),"days")))+xlab("")+ylab("Frequency"),
          BA.2_plot+geom_label(color="black",fill="white",size=3.5,aes(x=25,y=96,label=paste("Emergence Period (5-50%):",round(BA.2_emergence,digits = 0),"days")))+xlab("")+ylab(""),
          BA.4.5_plot+geom_label(color="black",fill="white",size=3.5,aes(x=39,y=96,label=paste("Emergence Period (5-50%):",round(BA.4.5_emergence,digits = 0),"days")))+ylab("Frequency"),
          XBB.1_plot+geom_label(color="black",fill="white",size=3.5,aes(x=54,y=96,label=paste("Emergence Period (5-50%):",round(XBB.1_emergence,digits = 0),"days")))+ylab(""), 
          nrow=2,ncol=2)
rm(BA.1_emergence,BA.1_plot,BA.2_emergence,BA.2_plot,BA.4.5_emergence,BA.4.5_plot,XBB.1_emergence,XBB.1_plot,p)



##### LOGISTIC GROWTH RATES ######


# Logistic growth rates of BA.1, BA.2, XBB.1
for(x in 1:length(variants)){
  model_data <- metadata_CT_pango %>% filter(Collection.date >= min(intervals[[x]]) & Collection.date <= max(intervals[[x]]) & CopyDate==1)%>%
    select(Collection.date,ShortLin) #Restrict to interval of interest 
  model_data$ShortLin <- ifelse(model_data$ShortLin == variants[x],model_data$ShortLin, "Other") #Filter to variant of interest 
  model_data$ShortLin <- factor(model_data$ShortLin, levels = c("Other",variants[x])) #Factor 
  model_data$emergence_days<- as.numeric(model_data$Collection.date-min(model_data$Collection.date -1)) #Change date to emergence days
  
  model <- glm(ShortLin ~ emergence_days, data=model_data, family=binomial) #Binomial logistic model 
  assign(paste0(variants[x],"_model"),model)
  
  model_data$ShortLin_num <- as.numeric(model_data$ShortLin) #Make the outcome dichotomous numeric and set the reference variant
  model_data$ShortLin_num <- model_data$ShortLin_num -1
  
  p <- ggplot(model_data, aes(x=emergence_days, y=ShortLin_num))+
    stat_smooth(color=paste0(colors[x]),method="glm",method.args=list(family=binomial),size=2)+
    ylab("Probability of Belonging to Variant Category")+
    xlab("Days Since 5% Emergence")+
    ggtitle(paste0(variants[x]))+
    scale_y_continuous(expand=c(0,0))+
    theme_classic()
  
  assign(paste0(variants[x],"_log_plot"),p) 
}

# Repeat this for BA.4 and 5 together
BA.5_data <- metadata_CT_pango %>% filter(Collection.date >= min(BA.5_interval) & Collection.date <= max(BA.5_interval) & CopyDate==1)%>%
  select(Collection.date,ShortLin)
BA.5_data$ShortLin <- ifelse(BA.5_data$ShortLin == "BA.5",BA.5_data$ShortLin, "Other")
BA.5_data$ShortLin <- factor(BA.5_data$ShortLin, levels = c("Other","BA.5"))
BA.5_data$Collection.date <- as.numeric(BA.5_data$Collection.date)
BA.5_data$emergence_days<- BA.5_data$Collection.date-min(BA.5_data$Collection.date -1)
BA.5_data$ShortLin_num <- as.numeric(BA.5_data$ShortLin)
BA.5_data$ShortLin_num <- BA.5_data$ShortLin_num -1
BA.5_model <- glm(ShortLin ~ emergence_days, data=BA.5_data, family=binomial)

BA.4_data <- metadata_CT_pango %>% filter(Collection.date >= min(BA.4_interval) & Collection.date <= max(BA.4_interval) & CopyDate==1)%>%
  select(Collection.date,ShortLin) 
BA.4_data$ShortLin <- ifelse(BA.4_data$ShortLin == "BA.4",BA.4_data$ShortLin, "Other")
BA.4_data$ShortLin <- factor(BA.4_data$ShortLin, levels = c("Other","BA.4"))
BA.4_data$Collection.date <- as.numeric(BA.4_data$Collection.date)
BA.4_data$emergence_days<- BA.4_data$Collection.date-min(BA.4_data$Collection.date -1)
BA.4_data$ShortLin_num <- as.numeric(BA.4_data$ShortLin)
BA.4_data$ShortLin_num <- BA.4_data$ShortLin_num -1
BA.4_model <- glm(ShortLin ~ emergence_days, data=BA.4_data, family=binomial)


BA.4.5_log_plot <- ggplot()+
  stat_smooth(data = BA.4_data,aes(x=emergence_days,y=ShortLin_num),method="glm",method.args=list(family=binomial),size=2,color="#CAC6EF")+
  stat_smooth(data = BA.5_data,aes(x=emergence_days,y=ShortLin_num),method="glm",method.args=list(family=binomial),size=2,color="#7294D4")+
  ylab("Probability of Belonging to Variant Category")+
  xlab("Days Since 5% Emergence")+
  ggtitle("BA.5")+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()

rm(BA.4_data,BA.5_data)

# Plot logistic growth rates
arranged_log_regressions <- plot_grid(BA.1_log_plot+xlab(""),
                             BA.2_log_plot+ylab("")+xlab(""),
                             BA.4.5_log_plot,
                             XBB.1_log_plot+ylab(""))

# Box plots for the slope coefficients 

box_plot_data<-data.frame(slope=c(BA.1_model$coefficients[2],BA.2_model$coefficients[2],BA.4_model$coefficients[2],BA.5_model$coefficients[2],
                                  XBB.1_model$coefficients[2]),
                          lower=c(confint(BA.1_model)[2,1],confint(BA.2_model)[2,1],confint(BA.4_model)[2,1],confint(BA.5_model)[2,1],
                                  confint(XBB.1_model)[2,1]),
                          upper=c(confint(BA.1_model)[2,2],confint(BA.2_model)[2,2],confint(BA.4_model)[2,2],confint(BA.5_model)[2,2],
                                  confint(XBB.1_model)[2,2]),
                          variant=c("BA.1","BA.2","BA.4","BA.5","XBB.1"))

log_box_plot <- ggplot(box_plot_data, aes(x=forcats::fct_rev(variant), y=slope, fill=variant))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0.3,size=0.78)+
  scale_fill_manual(values=c("#184E27","#CE4E50","#CAC6EF","#7294D4","#892F2E"))+
  ylab("Slope")+
  xlab("Variant Category")+
  theme_bw()+
  ggtitle("")+
  theme(legend.position = "none")+
  coord_flip()+
  theme_classic()


arranged_log_box_plots <- plot_grid(arranged_log_regressions, log_box_plot, rel_widths = c(1,0.6),ncol=2,nrow=1,labels=c("A","B"))


