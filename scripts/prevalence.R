library(ggplot2)
library(lubridate)
library(scales)
library(bdskytools)
library("readxl")
library(ggh4x)
library(tidyverse)
library(gg.gap)
library(zoo)

#read Ne Skygrid skyline
BA.2.2.skygrid <- read_tsv("../results/Ne/BA.2.2.afa.skyline.tsv") %>% 
  filter(Time > "2022") %>% 
  mutate(date = date_decimal(Time,tz = "EST"))

# smooth
BA.2.2.skygrid.z <- zoo(BA.2.2.skygrid$Mean,BA.2.2.skygrid$date)
sq <- seq(as.Date(start(BA.2.2.skygrid.z)), as.Date(end(BA.2.2.skygrid.z), frac = 1), by = "day")
zd <- na.approx(BA.2.2.skygrid.z, x = as.Date, xout = sq)

BA.2.2.skygrid.smooth<- as.data.frame(zd)
names(BA.2.2.skygrid.smooth) <- c("Mean")
date <- rownames(BA.2.2.skygrid.smooth)
rownames(BA.2.2.skygrid.smooth) <- NULL
BA.2.2.skygrid.smooth <- cbind(BA.2.2.skygrid.smooth,date)

BA.2.2.skygrid.z <- zoo(BA.2.2.skygrid$Upper,BA.2.2.skygrid$date)
sq <- seq(as.Date(start(BA.2.2.skygrid.z)), as.Date(end(BA.2.2.skygrid.z), frac = 1), by = "day")
zd <- na.approx(BA.2.2.skygrid.z, x = as.Date, xout = sq)
BA.2.2.skygrid.smooth$Upper <- zd

BA.2.2.skygrid.z <- zoo(BA.2.2.skygrid$Lower,BA.2.2.skygrid$date)
sq <- seq(as.Date(start(BA.2.2.skygrid.z)), as.Date(end(BA.2.2.skygrid.z), frac = 1), by = "day")
zd <- na.approx(BA.2.2.skygrid.z, x = as.Date, xout = sq)
BA.2.2.skygrid.smooth$Lower <- zd

#BDSKY skyline
lf.BA.2.2 <- readLogfile("Beast2/revision/xml_edited/HK_BA_2_2_un_40.fasta.edited.final.v2.0.1.log", burnin=0)

Re_sky.BA.2.2    <- getSkylineSubset(lf.BA.2.2, "reproductiveNumber")
Re_hpd.BA.2.2   <- getMatrixHPD(Re_sky.BA.2.2)
delta_hpd.BA.2.2 <- getHPD(lf.BA.2.2$becomeUninfectiousRate)

plotSkyline(1:16, Re_hpd.BA.2.2, type='step', ylab="R")


origin_max <- mean(lf.BA.2.2$origin_BDSKY_Serial)
treeheight <- mean(lf.BA.2.2$TreeHeight)
#info_span <- 1.667
timegrid.BA.2.2 <- seq(0,origin_max,length.out=96)
Re_gridded.BA.2.2 <- gridSkyline(Re_sky.BA.2.2, lf.BA.2.2$origin, timegrid.BA.2.2)
Re_gridded_hpd.BA.2.2 <- getMatrixHPD(Re_gridded.BA.2.2)

times.BA.2.2 <- 2022.317-timegrid.BA.2.2

df.BA.2.2.Re <- as.data.frame(t(as.data.frame(Re_gridded_hpd.BA.2.2)))

rownames(df.BA.2.2.Re) <- NULL
colnames(df.BA.2.2.Re) <- c("Lower","Value","Upper")

df.BA.2.2.Re$Time <- times.BA.2.2
df.BA.2.2.Re$date <- as.Date(date_decimal(df.BA.2.2.Re$Time,tz = "EST"))

df.BA.2.2.Re.pre <- df.BA.2.2.Re %>% filter(date > "2022-01-05") %>% 
  mutate(Re = Value,
         date = as.character(date)) %>% 
  subset(,c("date","Re"))


smooth_per_day <- function(value_list,date_list){
  zo <- zoo(value_list,date_list)
  sq <- seq(as.Date(start(zo)), as.Date(end(zo), frac = 1), by = "day")
  zd <- na.approx(zo, x = as.Date, xout = sq)
  df.zd <- as.data.frame(zd)
  date <- rownames(df.zd)
  rownames(df.zd) <- NULL
  #return two columns: zd and date
  return (cbind(df.zd,date))
}

df.BA.2.2.Re.pre.smooth <- smooth_per_day(df.BA.2.2.Re.pre$Re,df.BA.2.2.Re.pre$date)
colnames(df.BA.2.2.Re.pre.smooth)<- c("Re","date")


#Prevalence
BA.2.2.estimate <- merge(BA.2.2.skygrid.smooth,df.BA.2.2.Re.pre.smooth,all.x = TRUE,by = "date")

k_values = c(0.05,0.1,0.15,0.2)
generation_times = c(2,3)
incidence.summary <- vector()
prevalence.summary <- vector()
daily.incidence.summary <- vector()

for(k in 1:length(k_values)){
  for(t in 1:length(generation_times)){
    BA.2.2.estimate$prevalence <- BA.2.2.estimate$Mean*BA.2.2.estimate$Re*(1+1/k_values[k])/(generation_times[t]/365)
    BA.2.2.estimate$prevalence_Upper <- as.numeric(BA.2.2.estimate$Upper*BA.2.2.estimate$Re*(1+1/k_values[k])/(generation_times[t]/365))
    BA.2.2.estimate$prevalence_Lower <- as.numeric(BA.2.2.estimate$Lower*BA.2.2.estimate$Re*(1+1/k_values[k])/(generation_times[t]/365))
    BA.2.2.estimate$k <- k_values[k]
    BA.2.2.estimate$generation_time <- generation_times[t]
    
    BA.2.2.estimate.sub <- BA.2.2.estimate %>% 
      filter(date < "2022-04-14")
    
    prevalence.summary <- rbind(prevalence.summary,BA.2.2.estimate.sub %>% subset(,c("date","k","generation_time","prevalence","prevalence_Upper","prevalence_Lower")))
    
    BA.2.2.estimate.sub$Time <- decimal_date(as.Date(BA.2.2.estimate.sub$date))
    
    #total cumulative incidence
    BA.2.2.z <- zoo(BA.2.2.estimate.sub$prevalence,BA.2.2.estimate.sub$Time)
    sq <- seq(start(BA.2.2.z), end(BA.2.2.z, frac = 1), by = 2.72/365)
    zd <- na.approx(BA.2.2.z, xout = sq)
    incidence <- as.data.frame(round(cumsum(zd[,1])))
    incidence <- tibble::rownames_to_column(incidence,"time") %>% rename(.,c('incidence'=`round(cumsum(zd[, 1]))`))
    
    #cumulative incidence upper
    BA.2.2.z <- zoo(BA.2.2.estimate.sub$prevalence_Upper,BA.2.2.estimate.sub$Time)
    sq <- seq(start(BA.2.2.z), end(BA.2.2.z, frac = 1), by = 2.72/365)
    zd <- na.approx(BA.2.2.z, xout = sq)
    incidence_Upper <- as.data.frame(round(cumsum(zd[,1])))
    incidence_Upper <- tibble::rownames_to_column(incidence_Upper,"time") %>% rename(.,c('incidence_Upper'=`round(cumsum(zd[, 1]))`))
    
    #cumulative incidence lower
    BA.2.2.z <- zoo(BA.2.2.estimate.sub$prevalence_Lower,BA.2.2.estimate.sub$Time)
    sq <- seq(start(BA.2.2.z), end(BA.2.2.z, frac = 1), by = 2.72/365)
    zd <- na.approx(BA.2.2.z, xout = sq)
    incidence_Lower <- as.data.frame(round(cumsum(zd[,1])))
    incidence_Lower <- tibble::rownames_to_column(incidence_Lower,"time") %>% rename(.,c('incidence_Lower'=`round(cumsum(zd[, 1]))`))
    
    
    incidence.data <- merge(merge(incidence,incidence_Upper),incidence_Lower)
    incidence.data$k <- k_values[k]
    incidence.data$generation_time <- generation_times[t]
    
    incidence.summary <- rbind(incidence.summary,incidence.data)
    
    #daily incidence
    incidence.data$date <- as.Date(date_decimal(as.numeric(incidence.data$time),tz = "EST"))
    incidence.data.z <- zoo(incidence.data$incidence,incidence.data$date)
    sq <- seq(as.Date(start(incidence.data.z)), as.Date(end(incidence.data.z), frac = 1), by = "day")
    zd <- na.approx(incidence.data.z, x = as.Date, xout = sq)
    
    daily_incidences <- vector()
    for (i in 1:length(zd)) {
      if(i == 1){
        daily_incidence <- as.numeric(zd[,1][i])
      }else{
        daily_incidence <- as.numeric(zd[,1][i]) - as.numeric(zd[,1][i-1])
      }
      daily_incidences <- c(daily_incidences,daily_incidence)
    }
    
    #upper
    incidence.data$date <- as.Date(date_decimal(as.numeric(incidence.data$time),tz = "EST"))
    incidence.data.z <- zoo(incidence.data$incidence_Upper,incidence.data$date)
    sq <- seq(as.Date(start(incidence.data.z)), as.Date(end(incidence.data.z), frac = 1), by = "day")
    zd <- na.approx(incidence.data.z, x = as.Date, xout = sq)
    
    daily_upper_incidences <- vector()
    for (i in 1:length(zd)) {
      if(i == 1){
        daily_upper_incidence <- as.numeric(zd[,1][i])
      }else{
        daily_upper_incidence <- as.numeric(zd[,1][i]) - as.numeric(zd[,1][i-1])
      }
      daily_upper_incidences <- c(daily_upper_incidences,daily_upper_incidence)
    }
    
    #lower
    incidence.data$date <- as.Date(date_decimal(as.numeric(incidence.data$time),tz = "EST"))
    incidence.data.z <- zoo(incidence.data$incidence_Lower,incidence.data$date)
    sq <- seq(as.Date(start(incidence.data.z)), as.Date(end(incidence.data.z), frac = 1), by = "day")
    zd <- na.approx(incidence.data.z, x = as.Date, xout = sq)
    
    daily_lower_incidences <- vector()
    for (i in 1:length(zd)) {
      if(i == 1){
        daily_lower_incidence <- as.numeric(zd[,1][i])
      }else{
        daily_lower_incidence <- as.numeric(zd[,1][i]) - as.numeric(zd[,1][i-1])
      }
      daily_lower_incidences <- c(daily_lower_incidences,daily_lower_incidence)
    }
    
    date <- row.names(as.data.frame(zd))
    df.daily_incidence <- as.data.frame(cbind(date,daily_incidences,daily_upper_incidences,daily_lower_incidences))
    df.daily_incidence$k <- k_values[k]
    df.daily_incidence$generation_time <- generation_times[t]
    
    daily.incidence.summary <- rbind(daily.incidence.summary,df.daily_incidence)
  }
}

#write_tsv(daily.incidence.summary,"../results/prevalence/daily_incidence_summary.tsv")
#write_tsv(incidence.summary,"../results/prevalence/incidence_summary.tsv")
#write_tsv(prevalence.summary,"../results/prevalence/prevalence_summary.tsv")

#epi data
epi_cases <- read_tsv("../data/HK_epi_data.tsv") %>% 
  replace(.,is.na(.),0) %>% 
  mutate(local_cases = `CHP_locally-acquired`+RAT_positive_cases,
         date = as.Date(date),
         cumulative_cases = cumsum(local_cases)) %>% 
  filter(date < as.Date("2022-04-12"))

#daily incidence
for(t in 1:length(generation_times)){
  BA.2.2.daily.incidence <- ggplot()+
    geom_line(data=filter(daily.incidence.summary,generation_time == generation_times[t]),aes(x=as.Date(date),y=as.numeric(daily_incidences),color=as.factor(k)))+
    geom_line(data=epi_cases,aes(x=as.Date(date),y=local_cases),color="black",size=0.8)+
    annotate(geom="text",x=as.Date("2022-01-30"), y=150000, label=bquote(tau~"="~.(generation_times[t])), size=3)+
    geom_ribbon(data=filter(daily.incidence.summary,generation_time == generation_times[t] & k == 0.1),aes(x=as.Date(date), ymin=as.numeric(daily_lower_incidences),ymax=as.numeric(daily_upper_incidences),fill=as.factor(k)), alpha=0.2,fill="#377eb8") +
    labs(x = "Date", y = "Incidence")+
    scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-06"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
    coord_trans(y = "sqrt")+
    scale_colour_brewer(palette = "Set1")+
    scale_fill_brewer(palette = "Set1")+
    scale_y_continuous(labels=c("0.02M","0.05M","0.1M","0.2M"),breaks=c(20000,50000,100000,200000),limits=c(0,220000))+
    #guides(color=guide_legend(nrow=1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
          #legend.direction = "horizontal",
          legend.position="none")
  
  BA.2.2.daily.incidence
  
  ggsave(sprintf("../results/prevalence/BA.2.2_daily_incidence_generation_time_%s.pdf",generation_times[t]),BA.2.2.daily.incidence,width=9.25,height=6,units="cm")
}

#cumulative incidence
for(t in 1:length(generation_times)){
  BA.2.2.incidence <- ggplot()+
    geom_line(data=filter(incidence.summary,generation_time == generation_times[t]),aes(x=as.Date(date_decimal(as.numeric(time),tz = "EST")),y=incidence,color=as.factor(k)))+
    geom_line(data=epi_cases,aes(x=as.Date(date),y=cumulative_cases),color="black",size=0.8)+
    annotate(geom="text",x=as.Date("2022-01-30"), y=3000000, label=bquote(tau~"="~.(generation_times[t])), size=3)+
    geom_ribbon(data=filter(incidence.summary,generation_time == generation_times[t] & k == 0.1),aes(x=as.Date(date_decimal(as.numeric(time),tz = "EST")), ymin=incidence_Lower,ymax=incidence_Upper,fill=as.factor(k)), alpha=0.2,fill="#377eb8") +
    labs(x = "Date", y = "Cumulative incidence")+
    scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-06"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
    coord_trans(y = "sqrt")+
    scale_colour_brewer(palette = "Set1")+
    scale_fill_brewer(palette = "Set1")+
    scale_y_continuous(limits = c(0,5500000),labels=c("0.2M","1M","3M","5M"),breaks=c(200000,1000000,3000000,5000000),expand = c(0.001,0),
                       name = "Cumulative incidence",sec.axis=sec_axis( trans=~./7400000, name="Percentage",
                                                                        breaks=c(200000/7400000,1000000/7400000,3000000/7400000,5000000/7400000),
                                                                        labels=c(round(200000/7400000,3),round(1000000/7400000,3),round(3000000/7400000,3),round(5000000/7400000,3))))+
    #guides(color=guide_legend(nrow=1))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
          #legend.direction = "horizontal",
          legend.position="none")
  
  BA.2.2.incidence
  
  ggsave(sprintf("prevalence/BA.2.2_incidence_generation_time_%s.pdf",generation_times[t]),BA.2.2.incidence,width=10.5,height=6,units="cm")
}

#prevalence
for(t in 1:length(generation_times)){
  BA.2.2.prevalence <- ggplot()+
    geom_line(data=filter(prevalence.summary,generation_time == generation_times[t]),aes(x=as.Date(date),y=prevalence,color=as.factor(k)))+
    #geom_line(data=epi_cases,aes(x=as.Date(date),y=local_cases),color="black")+
    geom_ribbon(data=filter(prevalence.summary,generation_time == generation_times[t] & k == 0.1),aes(x=as.Date(date),ymin=prevalence_Lower,ymax=prevalence_Upper,fill=as.factor(k)), alpha=0.2,fill="#377eb8") +
    labs(x = "Date", y = "Point prevelance")+
    annotate(geom="text",x=as.Date("2022-01-30"), y=200000, label=bquote(tau~"="~.(generation_times[t])), size=3)+
    scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-06"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
    coord_trans(y = "sqrt")+
    scale_colour_brewer(palette = "Set1")+
    scale_fill_brewer(palette = "Set1")+
    scale_y_continuous(labels=c("0.02M","0.1M","0.2M","0.4M","0.6M"),breaks=c(20000,100000,200000,400000,600000),limits=c(0,600000))+
    #scale_y_continuous(labels=label_comma(),breaks=c(0,10,50,100,200,400,800,1000,2000,4000,8000,50000,100000,150000,200000))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
          legend.position="none")
  
  BA.2.2.prevalence
  
  ggsave(sprintf("prevalence/BA.2.2_prevalence_generation_time_%s.pdf",generation_times[t]),BA.2.2.prevalence,width=9.25,height=6,units="cm")
}

