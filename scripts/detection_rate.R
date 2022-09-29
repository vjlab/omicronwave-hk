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

#detection cases
detection.case.input <- read_tsv("../data/HK_epi_data.tsv") %>% 
  filter(date < "2022-03-30") %>% 
  subset(,c("date","PCR_new_tests_smoothed"))

BA.2.2.detection.case <- ggplot()+
  geom_line(data=detection.case.input,aes(x=as.Date(date),y=PCR_new_tests_smoothed),color="black")+
  labs(x = "Date", y = "PCR case detection")+
  scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-06"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
        legend.position="none")

BA.2.2.detection.case
ggsave("../results/relative detection rate/BA.2.2_detection_case.pdf",BA.2.2.detection.case,width=10,height=6,units="cm")

#relative detection rate
epi.data <- read_tsv("../data/HK_epi_data.tsv") %>% 
  filter(date < "2022-03-30") %>% 
  subset(,c("date","CHP_dashboard_new_daily_cs","PCR_new_tests_smoothed")) %>% 
  mutate(positive_rate = CHP_dashboard_new_daily_cs/PCR_new_tests_smoothed,
         date = as.character(date))


BA.2.2.detection.rate <- ggplot()+
  geom_line(data=epi.data,aes(x=as.Date(date),y=positive_rate),color="black")+
  labs(x = "Date", y = "PCR positivity rate")+
  scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-06"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
        legend.position="none")

BA.2.2.detection.rate
ggsave("../results/relative detection rate/BA.2.2_detection_case_rate.pdf",BA.2.2.detection.rate,width=10,height=6,units="cm")

#HK population
population_size <- 7400000
sens <- 1
spec <- 1

detection.input <- merge(BA.2.2.skygrid.smooth,epi.data,all.x = TRUE, by="date") %>% 
  mutate(detection_rate_mean = (population_size/Mean)*(positive_rate-1+spec)/(sens-1+spec)*(1-(1-1/population_size)^PCR_new_tests_smoothed),
         detection_rate_upper = (population_size/Upper)*(positive_rate-1+spec)/(sens-1+spec)*(1-(1-1/population_size)^PCR_new_tests_smoothed),
         detection_rate_lower = (population_size/Lower)*(positive_rate-1+spec)/(sens-1+spec)*(1-(1-1/population_size)^PCR_new_tests_smoothed))


BA.2.2.detection <- ggplot()+
  geom_line(data=detection.input,aes(x=as.Date(date),y=detection_rate_mean),color="black")+
  geom_ribbon(data=detection.input, aes(x=as.Date(date),ymin=detection_rate_lower,ymax=detection_rate_upper), fill="#BABABA", alpha=0.5) +
  labs(x = "Date", y = "Relative case detection rate")+
  scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-06"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
  scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
                labels = trans_format("log2", math_format(2^.x)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
        legend.position="none")

BA.2.2.detection
ggsave("../results/relative detection rate/BA.2.2_relative_detection_rate.pdf",BA.2.2.detection,width=10,height=6,units="cm")

