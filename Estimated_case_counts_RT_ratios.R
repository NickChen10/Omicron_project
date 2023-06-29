#Description: Cumulative case counts by variant
#Author: Nicholas Chen

#Load required packages
library(tidyverse)
library(ggplot2)
library(forcats)

##### DATA FORMATTING #####
path=getwd()

# Read COVIDestim data for Connecticut   
data <- read.csv(file=paste0(path,"/data/SARS-CoV-2 effective reproduction number in Connecticut since November 2021.csv")) 
data$days <- as.Date(data$days)



##### PLOTTING #####
# Filter to variants of interest
variants<- c("Omicron (BA.1)","Omicron (BA.2)","Omicron (BA.4)","Omicron (BA.5)","Omicron (XBB)" )
p1_data<-data %>% filter(variant %in% variants)%>% group_by(variant) %>% summarize(n=sum(I))

# Cumulative estimated cases 
ggplot(p1_data, aes(variant))+
  geom_bar(stat="identity",aes(x=variant,y=n,fill=variant))+ 
  scale_fill_manual(values=c("#184e27","#CE4E50","#CAC6EF","#7294D4","#892F2E"))+
  scale_y_continuous(breaks = seq(0,1600000,by = 200000),limits=c(0,1600000),expand=c(0,0))+
  ylab("Cumulative Estaimted Cases")+
  xlab("Variant")+
  theme_classic()+
  theme(legend.position = "none")
  
  
# RT ratio plots
# Filter to periods of Rt overlap and calculate ratios
p2_data <- data %>% select(days,variant, Rt) %>% filter(!is.na(Rt)) %>% spread(key=variant,value=Rt)%>% 
  mutate(BA1_Delta =  `Omicron (BA.1)`/`Delta`, BA2_BA1 = `Omicron (BA.2)`/`Omicron (BA.1)`, BA4_BA2 = `Omicron (BA.4)`/`Omicron (BA.2)`,
           BA5_BA2 = `Omicron (BA.5)`/`Omicron (BA.2)`, BA5_BA4 = `Omicron (BA.5)`/`Omicron (BA.4)`, XBB_BA5 = `Omicron (XBB)`/`Omicron (BA.5)`)
p2_data <- p2_data %>% select(BA1_Delta,BA2_BA1,BA4_BA2,BA5_BA2,BA5_BA4,XBB_BA5) %>% gather(key="pair",value="ratio",BA1_Delta,BA2_BA1,BA4_BA2,BA5_BA2,BA5_BA4,XBB_BA5)%>%
  filter(!is.na(ratio))

# Restrict data to the first 14 days of overlap (emergence)
p2_data <- p2_data %>% group_by(pair) %>% slice(1:min(n(), 14)) 
  
  
p2_data$pair <- factor(p2_data$pair, levels=c('BA1_Delta','BA2_BA1','BA4_BA2','BA5_BA2','BA5_BA4','XBB_BA5'),labels = c("BA.1 v. Delta", "BA.2 v. BA.1","BA.4 v. BA.2","BA.5 v. BA.2","BA.5 v. BA.4","XBB.1 v. BA.5")) 



ggplot(p2_data, aes(x=fct_rev(pair),y=ratio,color=pair,group=pair))+
  geom_jitter(width = 0.15)+
  scale_color_manual(values=c("#184e27","#CE4E50","#CAC6EF","#7294D4","#7294D4","#892F2E"))+
  geom_hline(yintercept = 1, linetype="dashed",size=1)+
  stat_summary(fun.y = median, fun.min = median, fun.max = median,
              geom = "crossbar", width = 0.5,color="#3d3d3d")+  
  coord_flip(clip = "off")+
  scale_y_continuous(limits = c(0.8,2),breaks=seq(0,2,0.5))+
  ylab("Rt ratio")+
  xlab("")+
  theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(
          size = 10, color = "black",face ="bold"), 
        axis.title=element_text(
          size=12, face ="bold"),
        title=element_text(size=14, face ="bold"),
        plot.margin = margin(0,2,0,0,"cm"),
    plot.background = element_blank())+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(pair) %>% mutate(median=((median(ratio)-1)*100)) %>% filter(!duplicated(median)), 
              aes(label = paste0(round(median,digits=2),"%"), x = pair), y = 2.1, size = 4, fontface = c("bold"),
              check_overlap = TRUE,color="black")
  
