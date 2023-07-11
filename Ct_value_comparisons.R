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
library(ggview)
library(patchwork)


##### DATA FORMATTING #####

#Read in data (remove directory at the end)
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
  scale_color_manual(values=c("#74875B","#FB7069"))+
  geom_signif(comparisons = list(c("BA.1","BA.2")),map_signif_level=T,color="black")+ 
  coord_fixed(ratio = 1/9.5, clip = 'off', expand = TRUE, ylim = c(40,12))+
  scale_y_reverse(breaks = seq(10,43,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  stat_summary(fun.data = median.n, geom = "label", colour = "black",
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
  scale_color_manual(values=c("#FB7069","#7294D4","#CAC6EF"))+  
  geom_signif(comparisons = list(c("BA.5","BA.2")),map_signif_level=T,color="black",test="t.test")+ 
  coord_fixed(ratio = 1/6.9, clip = 'off', expand = TRUE, ylim = c(40,12))+
  scale_y_reverse(breaks = seq(10,43,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  stat_summary(fun.data = median.n, geom = "label", colour = "black",
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
  coord_fixed(ratio = 1/10, clip = 'off', expand = TRUE, ylim = c(40,12))+
  scale_y_reverse(breaks = seq(10,43,by = 5))+
  scale_x_discrete()+
  theme_classic()+theme(legend.position="none",axis.text = element_text(
    size = 10, color = "black",face ="bold"), axis.title=element_text(
      size=12, face ="bold"),title=element_text(size=14, face ="bold"),
    plot.margin = margin(0,0,1,0, "cm"),plot.background = element_blank())+
  stat_summary(fun.data = median.n, geom = "label", colour = "black",
               size = 3, fontface = "italic")+
  geom_text(inherit.aes = TRUE, data = . %>% group_by(`ShortLin`) %>% count(), 
            aes(label = paste0("n=",n), x = `ShortLin`), y = -43, size = 3, fontface = c("italic"),
            check_overlap = TRUE, color="black")

plot_grid(p1,p2,p3,nrow=1)



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
OmiCols <- c("BA.1" = "#74875B", 
             "BA.1.1"= "#184E27",
             "BA.2" = "#FB7069",
             "BA.2.12.1" = "#CE4E50",
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
  geom_point(size=0.2,alpha=0.7, color="white",show.legend=FALSE)+
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
ggsave(filename="CToverview.png", units = "cm", width = 30, height = 10, scale =0.7)


# Plot scatterplots together
(variants_graphs_1st_half[[1]] + plot_spacer() +
   variants_graphs_1st_half[[2]] + variants_graphs_2nd_half[[2]] + 
   variants_graphs_1st_half[[3]]+ variants_graphs_2nd_half[[3]] +
   variants_graphs_1st_half[[4]]+ variants_graphs_2nd_half[[4]]+ 
   variants_graphs_1st_half[[5]]+  variants_graphs_2nd_half[[1]]  +
   plot_layout(ncol = 2))
ggsave(filename="variantCT.png", units = "cm", scale = 1.4)

