#Description: Cumulative case counts by variant
#Author: Nicholas Chen

#Load required packages
library(tidyverse)
library(ggplot2)
library(forcats)
library(ggforce)

##### DATA FORMATTING #####
path=getwd()

# Read COVIDestim data for Connecticut   
data <- read.csv(file=paste0(path,"/GitHub/Omicron_Project/data/SARS-CoV-2 effective reproduction number in Connecticut since November 2021.csv")) 
data$days <- as.Date(data$days)



##### PLOTTING #####
# Filter to variants of interest
variants<- c("Omicron (BA.1)","Omicron (BA.2)","Omicron (BA.4)","Omicron (BA.5)","Omicron (XBB)" )
p1_data<-data %>% filter(variant %in% variants)%>% group_by(variant) %>% summarize(n=sum(I))
  # Note that I is the variant specific estimated infections (`infections` is not variant specific and loops over)

# Plot 1: Cumulative estimated cases 
ggplot(p1_data, aes(variant))+
  geom_bar(stat="identity",aes(x=variant,y=n,fill=variant))+ 
  geom_hline(yintercept=c(5e+05,1e+06), linetype="dashed",size=0.8,alpha=0.6)+
  scale_fill_manual(values=c("#184e27","#CE4E50","#CAC6EF","#7294D4","#892F2E"))+
  scale_y_continuous(breaks = seq(0,1600000,by = 100000),limits=c(0,1600000),expand=c(0,0),labels = function(x) format(x, scientific = TRUE))+
  ylab("Cumulative Estimated Cases")+
  xlab("")+
  ggtitle("Cumulative Estimated Cases in Connecticut")+
  theme_classic()+
  theme(legend.position = "none",
        plot.title = element_text(size = 12,hjust = -0.2, face = "bold"),
        axis.text.x = element_text(size = 9,angle = 45, vjust = 1, hjust=1, face = "bold"),
        axis.text.y = element_text(size = 11,hjust = 0.5),
        )
#ggsave(filename="barplot.png", width=14,height=14,units = "cm")




# Plot 2: Rt values over time
p2_data <- data %>% filter(variant %in% variants) 

ggplot(p2_data, aes(x = days, y = as.numeric(Rt), group = variant, color = variant)) +
  geom_line(linewidth = 1.2) + 
  geom_ribbon(aes(ymin = as.numeric(rtlowci), 
                  ymax = as.numeric(rtupci), fill = variant), 
              alpha=0.2, 
              linetype="blank",
              color="grey") +
  geom_hline(yintercept = 1) +
  labs(
    title = "Variant-specific effective reproduction number (Rt)",
    y = "Rt value",
    x = ""
  ) + 
  scale_fill_manual(values = c("#184e27","#CE4E50","#CAC6EF","#7294D4","#892F2E")) +
  scale_color_manual(values = c("#184e27","#CE4E50","#CAC6EF","#7294D4","#892F2E")) +
  theme_classic()+
  scale_x_date(expand = c(0,0),limits =c(as.Date("2021-12-02"),as.Date("2023-01-20")))+
  theme(
    plot.title = element_text(size = 15,hjust = 0, face = "bold"),
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 12, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12, hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_text(size = 17, hjust = 0.5)
  ) 

#ggsave(filename="RT_plot.png", width = 20, height = 10, scale = 0.5)



# Plot 3: Variant estimated infections over time
ggplot(p2_data, aes(x = days, y = as.numeric(I), group = variant, color = variant)) +
  geom_line(linewidth = 1.2) + 
  geom_ribbon(aes(ymin = as.numeric(I_025), 
                  ymax = as.numeric(I_975), fill = variant), 
              alpha=0.2, 
              linetype="blank",
              color="grey") +
  geom_hline(yintercept = 1) +
  labs(
    title = "Estimated infections per variant in Connecticut",
    y = "Number of estimated infections",
    x = ""
  ) + 
  scale_fill_manual(values = c("#184e27","#CE4E50","#CAC6EF","#7294D4","#892F2E")) +
  scale_color_manual(values = c("#184e27","#CE4E50","#CAC6EF","#7294D4","#892F2E")) +
  theme_classic()+
  scale_x_date(expand = c(0,0),limits =c(as.Date("2021-12-02"),as.Date("2023-01-20")))+
  scale_y_continuous(expand=c(0,0))+
  theme(
    plot.title = element_text(size = 15,hjust = 0, face = "bold"),
    axis.title.x = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(size = 12, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 12, hjust = 0.5),
    legend.position = "none",
    axis.title.y = element_text(size = 17, hjust = 0.5)
  ) 

ggsave(filename="estimated_infections_plot.png", width = 15, height = 10, scale = 0.5)


# RT ratio plots
# Filter to periods of Rt overlap and calculate ratios
p4_data <- data %>% select(days,variant, Rt) %>% filter(!is.na(Rt)) %>% spread(key=variant,value=Rt)%>% 
  mutate(BA1_Delta =  `Omicron (BA.1)`/`Delta`, BA2_BA1 = `Omicron (BA.2)`/`Omicron (BA.1)`, BA4_BA2 = `Omicron (BA.4)`/`Omicron (BA.2)`,
           BA5_BA2 = `Omicron (BA.5)`/`Omicron (BA.2)`, BA5_BA4 = `Omicron (BA.5)`/`Omicron (BA.4)`, XBB_BA5 = `Omicron (XBB)`/`Omicron (BA.5)`)
p4_data <- p4_data %>% select(BA1_Delta,BA2_BA1,BA4_BA2,BA5_BA2,BA5_BA4,XBB_BA5) %>% gather(key="pair",value="ratio",BA1_Delta,BA2_BA1,BA4_BA2,BA5_BA2,BA5_BA4,XBB_BA5)%>%
  filter(!is.na(ratio))

# Restrict data to the first 14 days of overlap (emergence)
p4_data <- p4_data %>% group_by(pair) %>% slice(1:min(n(), 14)) 
  
  
p4_data$pair <- factor(p4_data$pair, levels=c('BA1_Delta','BA2_BA1','BA4_BA2','BA5_BA2','BA5_BA4','XBB_BA5'),labels = c("BA.1 v. Delta", "BA.2 v. BA.1","BA.4 v. BA.2","BA.5 v. BA.2","BA.5 v. BA.4","XBB.1 v. BA.5")) 


# Plot 4: RT ratios in periods of emergence 
ggplot(p4_data, aes(x=fct_rev(pair),y=ratio,color=pair,group=pair))+
  geom_jitter(width = 0.15)+
  scale_color_manual(values=c("#184e27","#CE4E50","#CAC6EF","#7294D4","#7294D4","#892F2E"))+
  geom_hline(yintercept = 1, linetype="dashed",size=1)+
  stat_summary(fun.y = median, fun.min = median, fun.max = median,
              geom = "crossbar", width = 0.5,color="#3d3d3d")+  
  coord_flip(clip = "off")+
  scale_y_continuous(limits = c(0.8,2),breaks=seq(0,2,0.5))+
  ylab("Rt ratio")+
  xlab("")+
  ggtitle("Rt Ratios of Variant Emergence")+
  theme_classic()+
  theme(legend.position="none",
        axis.text = element_text(size = 10, color = "black",face ="bold"), 
        axis.title=element_text(size=12, face ="bold"),
        title=element_text(size=12, face ="bold"),
        plot.margin = margin(0,2,0,0,"cm"),
        plot.background = element_blank())+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(pair) %>% mutate(median=((median(ratio)-1)*100)) %>% filter(!duplicated(median)), 
              aes(label = paste0(round(median,digits=2),"%"), x = pair), y = 2.1, size = 4, fontface = c("bold"),
              check_overlap = TRUE,color="black")
  
#ggsave(filename="rt_ratios.png", units = "cm")
