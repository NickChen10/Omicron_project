#Description: Using statistical analyses to compare CT values of variants in periods of emergence and corcirculation
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


##### DATA FORMATTING #####

#Read in data
data <- read_excel("C:/Users/nicho/GLab Dropbox/Nicholas Chen/Quantiative_Epi/XBB Project/GLab_SC2_sequencing_data.xlsx", 
                   col_types = c("text", "text", "numeric", 
                                 "text", "date", "text", "text", "date", 
                                 "text", "text", "numeric", "text", 
                                 "text", "text", "text", "text", "text", 
                                 "text", "text", "text", "numeric"))
names(data) <- to_any_case(names(data),case="snake")
data$n1_ge_ml <- as.numeric(data$yale_n_1_ge_m_l)
data[,c("extraction_pcr_date","collection_date")]<-lapply(data[,c("extraction_pcr_date","collection_date")],as.Date)
data <- data %>% filter(filter=="connecticut") %>% select(sample_id,collection_date,yale_n_1_fam,lineage, n1_ge_ml) %>% rename(n1 = yale_n_1_fam) 


#Recode lineages with aggregated lineages
data <- data %>%
  mutate(ShortLin = case_when(
    endsWith(lineage, "2.12.1") ~ "BA.2",
    startsWith(lineage, "Q") ~ "Alpha",
    startsWith(lineage, "B.1.1.7") ~ "Alpha",
    startsWith(lineage, "AY") ~ "Delta",
    startsWith(lineage, "B.1.617") ~ "Delta",
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
    startsWith(lineage, "BF.7") ~ "BA.5",
    startsWith(lineage, "BF") ~ "BA.5", #BF.7 as BA.5
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
    startsWith(lineage, "XBB.1") ~ "XBB.1", #Call all XBBs 1.5
    TRUE ~ "Other"
  ))

data$ShortLin <- factor(data$ShortLin, levels=c("Other","Alpha","Delta","BA.1","BA.2","BA.4","BA.5","XBB.1"))


##### CT VALUES/ GENOME EQUIVALENTS OVER FULL STUDY PERIOD #####

# Using GE/mL 
data<- data %>% filter(!(is.na(n1_ge_ml)) & !is.na(n1))

ggplot(data,aes(x=ShortLin,y=n1_ge_ml, color=ShortLin))+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Variant GE/mL",x=NULL, y = "Ge/mL (N)")+
  scale_color_manual(values=c("#515151","#0047ab","#ff6e40","#184e27","#CE4E50","#CAC6EF","#7294D4","#892F2E"))+
  coord_fixed(ratio = 1/2, clip = 'off', expand = TRUE, ylim=c(3e+03, 4e+10))+
  scale_x_discrete()+
  theme_classic()+
  theme(legend.position="none",axis.text = element_text(
   size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"),
   plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  stat_summary(fun.data = function(x){
    return(c(y = 10.6, label =round((10^(median(x))),0)))
    }, 
    geom = "label", colour = "black",
               size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(`ShortLin`) %>% count(), 
                      aes(label = paste0("n=",n), x = `ShortLin`, y=290),size = 3, fontface = c("italic"),
                      check_overlap = T,color="black")+
  scale_y_log10()



#Using Ct values
ggplot(data,aes(x=ShortLin,y=n1, color=ShortLin))+
  geom_hline(yintercept = 25.2,colour = "lightgrey")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Ct Values by Variant",x=NULL, y = "Ge/mL (N)")+
  scale_color_manual(values=c("#515151","#0047ab","#ff6e40","#184e27","#CE4E50","#CAC6EF","#7294D4","#892F2E"))+
  coord_fixed(ratio = 1/10, clip = 'off', expand = TRUE, ylim = c(40,12))+
  scale_y_reverse(breaks = seq(10,43,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  stat_summary(fun.data = function(x){
    return(c(y = -12, label = -1*round(median(x),1)))
  }, 
  geom = "label", colour = "black",
             size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(`ShortLin`) %>% count(), 
          aes(label = paste0("n=",n), x = `ShortLin`), y = -44, size = 3, fontface = c("italic"),
          check_overlap = TRUE,color="black")



###### CO-CIRCULATION PERIODS ######

# Identify periods of variant overlap to define emergence periods  
data_freq <- data %>% group_by(collection_date,ShortLin) %>% summarize(n=length(ShortLin))%>% mutate(freq=round(100*n/sum(n),2)) %>%
  filter(collection_date >= "2021-12-01" & ShortLin != "Other") %>%mutate(date_num = as.numeric(collection_date))
variants <- unique(data_freq$ShortLin)

p1<-ggplot(data_freq, aes(x=as.Date(collection_date), y=freq, color=ShortLin))+
  geom_smooth(se=F)+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+  
  xlab("Date")+
  theme_classic()

ggplotly(p1)
# Crossover points :
  # 3-13-2022
  # 6-25-2022
  # 12-17-2022
cross_overs <- as.Date(c("2022-03-13","2022-06-25","2022-12-17"))

# Define co-circulation periods as 14 days on either side of the crossovers
bounds<- as.Date(c((cross_overs[1]-14),(cross_overs[1]+14),(cross_overs[2]-14),(cross_overs[2]+14),(cross_overs[3]-14),(cross_overs[3]+14)))



###### CT VALUE STATISTICAL TESTS ######
data_stats <- data %>% filter(ShortLin != "Other")

# Define intervals to test
for(i in c(1,3,5)){
  interval <- data_stats %>% filter(collection_date >= bounds[i] & collection_date <= bounds[i+1])
  assign(paste0("interval_",c(1,1,2,2,3)[i]),interval)
}



# Assess normalcy of each variant in each interval via Shapiro tests
BA.1_1 <- interval_1 %>% filter(ShortLin== "BA.1") 
BA.2_1 <- interval_1 %>% filter(ShortLin== "BA.2")
BA.2_2 <- interval_2 %>% filter(ShortLin== "BA.2")
BA.5_2 <- interval_2 %>% filter(ShortLin== "BA.5")
BA.4_2 <- interval_2 %>% filter(ShortLin== "BA.4")
BA.4_3 <- interval_3 %>% filter(ShortLin== "BA.4")
BA.5_3 <- interval_3 %>% filter(ShortLin== "BA.5")
XBB.1_3 <- interval_3 %>% filter(ShortLin== "XBB.1")

vars<- list(BA.1_1,BA.2_1,BA.2_2,BA.5_2,BA.4_2,BA.4_3,BA.5_3,XBB.1_3)
var_names <-c('BA.1_1','BA.2_1','BA.2_2','BA.5_2','BA.4_2','BA.4_3','BA.5_3','XBB.1_3') 

for(i in 1:length(vars)){
  normalcy<-shapiro.test(vars[[i]]$n1)
  if(normalcy$p.value <=0.05){
    print(paste0(var_names[i]," is not normal"))
  } else{
    print(paste0(var_names[i]," is normal"))
  }
}
  # Some are normal, some aren't 
  # Given that n>30 for all but (BA.4_3), can assume normalcy for the rest

rm(BA.1_1,BA.2_1,BA.2_2,BA.5_2,BA.4_2,BA.4_3,BA.5_3,XBB.1_3, vars,var_names)


# Assess variance via Variance tests and Bartlett tests
interval_1 <- interval_1 %>% filter(ShortLin!= "Delta")
var.test(n1 ~ ShortLin, interval_1) #Interval 1: Equal variance 
bartlett.test(n1 ~ ShortLin, data=interval_2) #Interval 2: Unequal variance 
interval_3 <- interval_3 %>% filter(ShortLin!= "BA.4") #Remove BA.4 in interval 3 due to low counts
var.test(BA.5_3$n1,XBB.1_3$n1) #Interval 3: Equal variance 


# Statistical Tests
t.test(n1~ShortLin,interval_1,var.equal = T) #0.003 (higher Ct in BA.1 than BA.2)
wilcox.test(n1~ShortLin,interval_1) #0.003
  # Good concordance between parametric and non-parametric t-tests 

oneway.test(n1 ~ ShortLin, data=interval_2, var.equal = FALSE) #0.0047
anova<- aov(n1 ~ ShortLin, data=interval_2) 
summary(anova) # 0.044
  # Good concordance between Welch ANOVA and simple ANOVA
TukeyHSD(anova, conf.level=.95) #Significant difference between BA.5 and BA.2

t.test(n1~ShortLin, interval_3, var.equal = T) #0.02 (higher in BA.5 than BA)
wilcox.test(n1~ShortLin, interval_3) #0.03
  # Good concordance between parametric and non-parametric t-tests 



##### PLOTTING #####
#Pull 95% CIs for each mean
p1 <- ggplot(interval_1,aes(x=ShortLin,y=n1, color=ShortLin))+
  geom_hline(yintercept = 25.2,colour = "lightgrey")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Emergence Period 1",x=NULL, y = "Ct (N)")+
  scale_color_manual(values=c("#184E27","#CE4E50"))+
  geom_signif(comparisons = list(c("BA.1","BA.2")),map_signif_level=T,color="black")+ 
  #coord_fixed(ratio = 1/9.5, clip = 'off', expand = TRUE, ylim = c(40,12))+
  coord_cartesian(clip = 'off', expand = TRUE, ylim = c(40,12))+
  scale_y_reverse(breaks = seq(10,43,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  stat_summary(fun.data = function(x){
    return(c(y = -12, label = -1*round(median(x),1)))}, 
  geom = "label", colour = "black",
  size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(`ShortLin`) %>% count(), 
            aes(label = paste0("n=",n), x = `ShortLin`), y = -43, size = 3, fontface = c("italic"),
            check_overlap = TRUE,color="black")

interval_2$ShortLin <- factor(interval_2$ShortLin, levels=c("BA.2","BA.5","BA.4"))
p2 <- ggplot(interval_2,aes(x=ShortLin,y=n1, color=ShortLin))+
  geom_hline(yintercept = 25.2,colour = "lightgrey")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Emergence Period 2",x=NULL, y = NULL)+
  scale_color_manual(values=c("#CE4E50","#7294D4","#CAC6EF"))+  
  geom_signif(comparisons = list(c("BA.5","BA.2")),map_signif_level=T,color="black"
              ,test="t.test")+ 
  #coord_fixed(ratio = 1/6.9, clip = 'off', expand = TRUE, ylim = c(40,12))+
  coord_cartesian(clip = 'off', expand = TRUE, ylim = c(40,12))+
  scale_y_reverse(breaks = seq(10,43,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  stat_summary(fun.data = function(x){
    return(c(y = -12, label = -1*round(median(x),1)))}, 
    geom = "label", colour = "black",
    size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(`ShortLin`) %>% count(), 
            aes(label = paste0("n=",n), x = `ShortLin`), y = -43, size = 3, fontface = c("italic"),
            check_overlap = TRUE,color="black")


p3 <- ggplot(interval_3,aes(x=ShortLin,y=n1, color=ShortLin))+
  geom_hline(yintercept = 25.2,colour = "lightgrey")+
  geom_jitter(alpha = 0.5,size = 0.5, stroke = 2,shape = 21, width = 0.15)+
  geom_boxplot(width=0.5,outlier.shape = NA,colour = "#666666",fill = NA)+
  labs(title="Emergence Period 3",x=NULL, y = NULL)+
  scale_color_manual(values=c("#7294D4","#892F2E"))+
  geom_signif(comparisons = list(c("BA.5","XBB.1")),map_signif_level=T,color="black",test="t.test")+ 
  #coord_fixed(ratio = 1/10, clip = 'off', expand = TRUE, ylim = c(40,12))+
  coord_cartesian(clip = 'off', expand = TRUE, ylim = c(40,12))+
  scale_y_reverse(breaks = seq(10,43,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  stat_summary(fun.data = function(x){
    return(c(y = -12, label = -1*round(median(x),1)))}, 
    geom = "label", colour = "black",
    size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(`ShortLin`) %>% count(), 
            aes(label = paste0("n=",n), x = `ShortLin`), y = -43, size = 3, fontface = c("italic"),
            check_overlap = TRUE, color="black")


(p1 + p2 + p3) #EDIT: For the same size plots with varying sized boxplots 
(p1 + plot_spacer() + p2 + plot_spacer() + p3 + plot_layout(widths=c(1,0.1,1.5,0.1,1))) #EDIT: For the same size boxplots + different size plots 
#ggsave(filename="CT_comparisons.png", units = "cm", width = 35, height = 20)



###### CT VALUE LINEAR REGRESSIONS ######


#Recode lineages including BA.1.1 and BA.2.12.1
data$ShortLin <- ""
lineages <- data %>% 
  mutate(ShortLin = case_when(
    endsWith(lineage, "2.12.1") ~ "BA.2.12.1",
    startsWith(lineage, "Q") ~ "Alpha",
    startsWith(lineage, "B.1.1.7") ~ "Alpha",
    startsWith(lineage, "AY") ~ "Delta",
    startsWith(lineage, "B.1.617") ~ "Delta",
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
    startsWith(lineage, "XBB.1") ~ "XBB.1.5", 
    TRUE ~ "Other"
  ))

lineages$ShortLin <- factor(lineages$ShortLin, levels=c("Other","Alpha","Delta","BA.1","BA.1.1","BA.2","BA.2.12.1","BA.4.6","BA.5","BQ.1.1","XBB.1.5"))




# Group data by week
lineages <- lineages %>% mutate(CopyDate = collection_date)%>% group_by(`week` = cut(`CopyDate`, "week"))
lineages$week <- as.Date(lineages$week)
lineages <- lineages[!is.na(lineages$`week`),] #delete data without week assigned


# Variant counts by week for emergence periods (CLEAN THIS USING MY DATA_FREQ CODE)
lineages_sum <- lineages %>% group_by(week) %>% count(ShortLin)
lineages_count  <- aggregate(lineages_sum$n, by=list(lineages_sum$week), FUN = sum)
colnames(lineages_count) <- c("week","Total")
lineages_sum  <- right_join(lineages_sum, lineages_count, by = "week")
lineages_sum <- lineages_sum %>% mutate(`freq` = ((`n`/`Total`)*100))



# First half of variant scatterplots
variants_1sthalf <- c("BA.1","BA.1.1","BA.2","BA.2.12.1","BA.5")
range <- list()
variant_data <- list()

# Define range of time when the variant is above 10% in the population 
for (i in variants_1sthalf) {
  range[[i]] <- range(lineages_sum$week[which(lineages_sum$ShortLin == i & lineages_sum$freq >= 10)])
  variant_data[[i]] <- lineages %>% filter(ShortLin == i) %>%
    filter(collection_date  >= min(range[[i]]) & 
             collection_date <= max(range[[i]]))
}


variants_graphs_1st_half <- list()


# Plot the first half of the figure 
OmiCols <- c("BA.1" = "#184E27", 
             "BA.1.1"= "#74875B",
             "BA.2" = "#CE4E50",
             "BA.2.12.1" = "#FB7069",
             "BA.4" = "#CAC6EF",
             "BA.4.6" = "#CAC6EF",
             "BA.5" = "#7294D4",
             "BQ.1.1" = "#184294",
             "XBB" = "#892F2E",
             "XBB.1" = "#892F2E",
             "XBB.1.5" = "#892F2E",
             "Others" = "grey"
)

for (i in variants_1sthalf) {
  variants_graphs_1st_half[[i]] <- ggplot(variant_data[[i]], 
                                          aes(x=collection_date, y=n1, 
                                              color=ShortLin, group=ShortLin)) +
    geom_point(size=0.8,alpha=0.5)+
    geom_smooth(method = "lm", se=T, aes(group=ShortLin, fill=ShortLin),size=1, alpha=0.1)+
    scale_y_reverse(breaks=c(15,25,35))+  
    ylab("Ct N")+
    theme_classic()+
    theme(
      axis.text.x=element_text(face = "bold"),
      axis.title.x=element_blank(),
      legend.box.background = element_rect(colour = "black"),
      legend.position = "none")+
    scale_fill_manual(values = c(OmiCols))+
    scale_color_manual(values = c(OmiCols))+
    coord_fixed(xlim= c(as.Date("2022-01-10"),as.Date("2022-07-20")), ylim = c(35,15))
}

# Repeat process for the second half of the figure
variants_2ndhalf <- c("BA.5","BA.4.6","BQ.1.1","XBB.1.5")

range_2 <- list()
variant_data_2 <- list()

for (i in variants_2ndhalf) {
  range_2[[i]] <- range(lineages_sum$week[which(lineages_sum$ShortLin == i & lineages_sum$freq >= 10)])
  variant_data_2[[i]] <- lineages %>% filter(ShortLin == i) %>%
    filter(collection_date  >= min(range_2[[i]]) & 
             collection_date <= max(range_2[[i]]))
}

variants_graphs_2nd_half <- list()

for (i in variants_2ndhalf) {
  
  variants_graphs_2nd_half[[i]] <- ggplot(variant_data_2[[i]] %>% filter(`ShortLin` == i)
                                          , aes(x=collection_date, y=n1, 
                                                color=ShortLin, group=ShortLin)) +
    geom_point(size=0.8,alpha=0.5)+
    geom_smooth(method = "lm", se=T, aes(group=ShortLin, fill=ShortLin),size=1, alpha=0.1)+
    #scale_x_date(breaks="1 months")+
    scale_y_reverse(breaks=c(15,25,35))+    
    #ylab("Ct N")+
    theme_classic()+
    theme(
      axis.text.x=element_text(face = "bold"),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.box.background = element_rect(colour = "black"),
      legend.position = "none")+
    scale_fill_manual(values = c(OmiCols))+
    scale_color_manual(values = c(OmiCols))+
    coord_fixed(xlim= c(as.Date("2022-07-20"),as.Date("2023-01-20")), ylim = c(35,15))
}


# Combine data into one plot
variants_combined <- rbind(variant_data[["BA.1"]],
                           variant_data[["BA.1.1"]],
                           variant_data[["BA.2"]],
                           variant_data[["BA.2.12.1"]],
                           variant_data[["BA.5"]],
                           variant_data_2[["BA.5"]],
                           variant_data_2[["BA.4.6"]],
                           variant_data_2[["BQ.1.1"]],
                           variant_data_2[["XBB.1.5"]]
)

# Combine linear regressions
overview <- ggplot(variants_combined, aes(x=collection_date, y=n1, color=ShortLin, group=ShortLin)) +
  geom_hline(yintercept = 25,colour = "lightgrey")+
  geom_point(size=0.2,alpha=0, show.legend=FALSE)+
  geom_rect(inherit.aes = F,data = data.frame(xmin=bounds[1],xmax=bounds[2],ymin=20,ymax=30), 
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="orange",alpha=0.1)+#EDIT: First emergence period 
  geom_rect(inherit.aes = F,data = data.frame(xmin=bounds[3],xmax=bounds[4],ymin=20,ymax=30), 
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="orange",alpha=0.1)+#EDIT: Second emergence period 
  geom_rect(inherit.aes = F,data = data.frame(xmin=bounds[5],xmax=bounds[6],ymin=20,ymax=30), 
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="orange",alpha=0.1)+#EDIT: First emergence period 
  geom_smooth(method = "lm", se=T, aes(group=ShortLin, fill=ShortLin),size=1.5, alpha=0.1,show.legend=FALSE)+
  scale_y_reverse()+
  coord_cartesian(ylim = c(29, 22))+
  xlab("Collection Date")+
  ylab("Ct N")+
  labs(
    title = "Viral load measured via qPCR Ct nucleocapsid (N) target")+
  theme_classic()+
  theme(
    axis.text.x=element_text(size = 12, face ="bold"),
    axis.text.y=element_text(size = 12, face ="bold"),
    axis.title = element_text(size = 12, face ="bold"),
    plot.title = element_text(size = 15, face ="bold")
  )+
  scale_fill_manual(values = c(OmiCols))+
  scale_color_manual(values = c(OmiCols))
  
overview
#ggsave(filename="CToverview.png", units = "cm", width = 30, height = 10, scale =0.7)


# Plot scatterplots together
(variants_graphs_1st_half[[1]] + plot_spacer() +
   variants_graphs_1st_half[[2]] + variants_graphs_2nd_half[[2]] + 
   variants_graphs_1st_half[[3]]+ variants_graphs_2nd_half[[3]] +
   variants_graphs_1st_half[[4]]+ variants_graphs_2nd_half[[4]]+ 
   variants_graphs_1st_half[[5]]+  variants_graphs_2nd_half[[1]]  +
   plot_layout(ncol = 2))
#ggsave(filename="variantCT.png", units = "cm", width = 30, height = 10, scale =0.9)




###### DATA TABLE GENERATION #####
# Want the median Ct values for each variant in the above plot for the first and last week (using the variants_combined dataset)
table_data <- variants_combined %>% select(ShortLin,n1,collection_date)


  
for(i in 1:length(unique(table_data$ShortLin))){
  x <- table_data %>% filter(ShortLin==unique(table_data$ShortLin)[i])
  start <- x %>% filter(week==min(x$week))
  end <- x %>% filter(week==max(x$week))
  first_week <- summary(start$n1)
  last_week <- summary(end$n1)
  sum<-data.frame(rbind(first_week,last_week))%>% mutate(ShortLin = unique(table_data$ShortLin)[i])
  assign(paste0(unique(table_data$ShortLin)[i],"_iqr"),sum)
  
  model<- lm(n1~collection_date,x)
  coef <- summary(model)$coefficients[2]
  p_val <- summary(model)$coefficients[8]
  model_output<- data.frame(ShortLin = unique(table_data$ShortLin)[i],coefficient=coef, p_val = p_val)
  assign(paste0(unique(table_data$ShortLin)[i],"_linear_model"),model_output)
}



#Issues with BA.5 and BQ.1.1 and BA.2.12.1 - do these manually 
BA.5 <- table_data %>% filter(ShortLin=="BA.5")
start_BA.5 <-BA.5 %>% filter(week==min(BA.5$week))
end_BA.5 <- BA.5 %>% filter(week=="2022-12-05")
first_week_BA.5 <- summary(start_BA.5$n1)
last_week_BA.5 <- summary(end_BA.5$n1)
BA.5_iqr <- data.frame(rbind(first_week_BA.5,last_week_BA.5))%>% mutate(ShortLin="BA.5")


BQ.1.1 <- table_data %>% filter(ShortLin=="BQ.1.1")
start_BQ.1.1 <-BQ.1.1 %>% filter(week==min(BQ.1.1$week))
end_BQ.1.1 <- BQ.1.1 %>% filter(week=="2022-12-26")
first_week_BQ.1.1 <- summary(start_BQ.1.1$n1)
last_week_BQ.1.1 <- summary(end_BQ.1.1$n1)
BQ.1.1_iqr <- data.frame(rbind(first_week_BQ.1.1,last_week_BQ.1.1))%>% mutate(ShortLin="BQ.1.1")


BA.2.12.1 <- table_data %>% filter(ShortLin=="BA.2.12.1")
start_BA.2.12.1 <-BA.2.12.1 %>% filter(week==min(BA.2.12.1$week))
end_BA.2.12.1 <- BA.2.12.1 %>% filter(week=="2022-07-11")
first_week_BA.2.12.1 <- summary(start_BA.2.12.1$n1)
last_week_BA.2.12.1 <- summary(end_BA.2.12.1$n1)
BA.2.12.1_iqr <- data.frame(rbind(first_week_BA.2.12.1,last_week_BA.2.12.1))%>% mutate(ShortLin="BA.2.12.1")



XBB.1.5 <- table_data %>% filter(ShortLin=="XBB.1.5")
start_XBB.1.5 <-XBB.1.5 %>% filter(week=="2022-12-05")
end_XBB.1.5 <- XBB.1.5 %>% filter(week=="2023-01-30")
first_week_XBB.1.5 <- summary(start_XBB.1.5$n1)
last_week_XBB.1.5 <- summary(end_XBB.1.5$n1)
XBB.1.5_iqr <- data.frame(rbind(first_week_XBB.1.5,last_week_XBB.1.5))%>% mutate(ShortLin="XBB.1.5")


export <- rbind(BA.1_iqr,BA.1.1_iqr,BA.2_iqr,BA.2.12.1_iqr,BA.4.6_iqr,BA.5_iqr,BQ.1.1_iqr,XBB.1.5_iqr)%>%
  select(ShortLin,q_1 = X1st.Qu.,med=Median,mean=Mean,q_3=X3rd.Qu.)


export_2 <- rbind(BA.1_linear_model,BA.1.1_linear_model,BA.2_linear_model,BA.2.12.1_linear_model,BA.4.6_linear_model,BA.5_linear_model,BQ.1.1_linear_model, XBB.1.5_linear_model)


#write.csv(export,"median_ct_values.csv")
#write.csv(export_2,"linear_regression_values.csv")

test <-lm(n1~collection_date,XBB.1.5)
range(fitted(test))
test_2 <- ggplot(x, aes(y=n1,x=collection_date))+
  geom_point()  +
  geom_smooth(method="lm")
