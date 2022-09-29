library(ggplot2)
library(lubridate)
library(scales)
library(bdskytools)
library("readxl")
library(ggh4x)
library(tidyverse)
library(gg.gap)
library(zoo)

#Skygrid skyline
BA.2.2.skygrid <- read_tsv("../results/Ne/BA.2.2.afa.skyline.tsv") %>% 
  filter(Time > "2022") %>% 
  mutate(date = date_decimal(Time,tz = "EST"))

#plot
BA.2.2.skyline <- ggplot()+
  geom_line(data=BA.2.2.skygrid,aes(x=as.Date(date_decimal(Time,tz = "EST")),y=Mean),color="#3cabd9")+
  geom_ribbon(data=BA.2.2.skygrid, aes(x=as.Date(date_decimal(Time,tz = "EST")),ymin=Lower,ymax=Upper), fill="#3cabd9", alpha=0.5) +
  labs(x = "Date", y = "Ne")+
  scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-06"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  annotation_logticks(size=0.2, short=unit(0.05, "cm"), mid=unit(0.1, "cm"), long=unit(0.15, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
        legend.position="none")

BA.2.2.skyline
ggsave("../results/Ne/BA.2.2_skygrid_skyline.pdf",BA.2.2.skyline,width=10,height=6,units="cm")
