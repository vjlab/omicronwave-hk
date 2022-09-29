library(ggplot2)
library(lubridate)
library(scales)
library(bdskytools)
library("readxl")
library(ggh4x)
library(tidyverse)
library(gg.gap)
library(zoo)

#BDSKY skyline

#data v1
lf.BA.2.2 <- readLogfile("Beast2/HK_BA_2_2_v1.burnin.10000000.log", burnin=0)

Re_sky.BA.2.2    <- getSkylineSubset(lf.BA.2.2, "reproductiveNumber")
Re_hpd.BA.2.2   <- getMatrixHPD(Re_sky.BA.2.2)
delta_hpd.BA.2.2 <- getHPD(lf.BA.2.2$becomeUninfectiousRate)

plotSkyline(1:10, Re_hpd.BA.2.2, type='step', ylab="R")


origin_max <- mean(lf.BA.2.2$origin_BDSKY_Serial)
treeheight <- mean(lf.BA.2.2$TreeHeight)
#info_span <- 1.667
timegrid.BA.2.2 <- seq(0,origin_max,length.out=121)
Re_gridded.BA.2.2 <- gridSkyline(Re_sky.BA.2.2, lf.BA.2.2$origin, timegrid.BA.2.2)
Re_gridded_hpd.BA.2.2 <- getMatrixHPD(Re_gridded.BA.2.2)

times.BA.2.2 <- 2022.317-timegrid.BA.2.2

df.BA.2.2.v1 <- as.data.frame(t(as.data.frame(Re_gridded_hpd.BA.2.2)))

rownames(df.BA.2.2.v1) <- NULL
colnames(df.BA.2.2.v1) <- c("Lower","Value","Upper")

df.BA.2.2.v1$Time <- times.BA.2.2
df.BA.2.2.v1$days <- (treeheight-timegrid.BA.2.2)*365
df.BA.2.2.v1$dataset <- "v1"

#data v2
lf.BA.2.2 <- readLogfile("Beast2/HK_BA_2_2_v2.burnin.0.1.log", burnin=0)

Re_sky.BA.2.2    <- getSkylineSubset(lf.BA.2.2, "reproductiveNumber")
Re_hpd.BA.2.2   <- getMatrixHPD(Re_sky.BA.2.2)
delta_hpd.BA.2.2 <- getHPD(lf.BA.2.2$becomeUninfectiousRate)

plotSkyline(1:10, Re_hpd.BA.2.2, type='step', ylab="R")

origin_max <- mean(lf.BA.2.2$origin_BDSKY_Serial)
treeheight <- mean(lf.BA.2.2$TreeHeight)

timegrid.BA.2.2 <- seq(0,origin_max,length.out=121)
Re_gridded.BA.2.2 <- gridSkyline(Re_sky.BA.2.2, lf.BA.2.2$origin, timegrid.BA.2.2)
Re_gridded_hpd.BA.2.2 <- getMatrixHPD(Re_gridded.BA.2.2)

times.BA.2.2 <- 2022.317-timegrid.BA.2.2

df.BA.2.2.v2 <- as.data.frame(t(as.data.frame(Re_gridded_hpd.BA.2.2)))

rownames(df.BA.2.2.v2) <- NULL
colnames(df.BA.2.2.v2) <- c("Lower","Value","Upper")

df.BA.2.2.v2$Time <- times.BA.2.2
df.BA.2.2.v2$days <- (treeheight-timegrid.BA.2.2)*365
df.BA.2.2.v2$dataset <- "v2"

#data v3
lf.BA.2.2 <- readLogfile("Beast2/HK_BA_2_2_v3.burnin.0.1.log", burnin=0)

Re_sky.BA.2.2    <- getSkylineSubset(lf.BA.2.2, "reproductiveNumber")
Re_hpd.BA.2.2   <- getMatrixHPD(Re_sky.BA.2.2)
delta_hpd.BA.2.2 <- getHPD(lf.BA.2.2$becomeUninfectiousRate)

plotSkyline(1:10, Re_hpd.BA.2.2, type='step', ylab="R")


origin_max <- mean(lf.BA.2.2$origin_BDSKY_Serial)
treeheight <- mean(lf.BA.2.2$TreeHeight)

timegrid.BA.2.2 <- seq(0,treeheight,length.out=121)
Re_gridded.BA.2.2 <- gridSkyline(Re_sky.BA.2.2, lf.BA.2.2$origin, timegrid.BA.2.2)
Re_gridded_hpd.BA.2.2 <- getMatrixHPD(Re_gridded.BA.2.2)

times.BA.2.2 <- 2022.317-timegrid.BA.2.2

df.BA.2.2.v3 <- as.data.frame(t(as.data.frame(Re_gridded_hpd.BA.2.2)))

rownames(df.BA.2.2.v3) <- NULL
colnames(df.BA.2.2.v3) <- c("Lower","Value","Upper")

df.BA.2.2.v3$Time <- times.BA.2.2
df.BA.2.2.v3$days <- (treeheight-timegrid.BA.2.2)*365
df.BA.2.2.v3$dataset <- "v3"

#data v4
lf.BA.2.2 <- readLogfile("Beast2/HK_BA_2_2_v4.burnin.0.1.log", burnin=0)

Re_sky.BA.2.2    <- getSkylineSubset(lf.BA.2.2, "reproductiveNumber")
Re_hpd.BA.2.2   <- getMatrixHPD(Re_sky.BA.2.2)
delta_hpd.BA.2.2 <- getHPD(lf.BA.2.2$becomeUninfectiousRate)

plotSkyline(1:10, Re_hpd.BA.2.2, type='step', ylab="R")

origin_max <- mean(lf.BA.2.2$origin_BDSKY_Serial)
treeheight <- mean(lf.BA.2.2$TreeHeight)

timegrid.BA.2.2 <- seq(0,origin_max,length.out=121)
Re_gridded.BA.2.2 <- gridSkyline(Re_sky.BA.2.2, lf.BA.2.2$origin_BDSKY_Serial, timegrid.BA.2.2)
Re_gridded_hpd.BA.2.2 <- getMatrixHPD(Re_gridded.BA.2.2)

times.BA.2.2 <- 2022.317-timegrid.BA.2.2

df.BA.2.2.v4 <- as.data.frame(t(as.data.frame(Re_gridded_hpd.BA.2.2)))

rownames(df.BA.2.2.v4) <- NULL
colnames(df.BA.2.2.v4) <- c("Lower","Value","Upper")

df.BA.2.2.v4$Time <- times.BA.2.2
df.BA.2.2.v4$days <- (treeheight-timegrid.BA.2.2)*365
df.BA.2.2.v4$dataset <- "v4"

df.BA.2.2 <- rbind(df.BA.2.2.v1,df.BA.2.2.v2,df.BA.2.2.v3,df.BA.2.2.v4)

BA.2.bdsky <- ggplot(df.BA.2.2)+
  geom_line(aes(x=as.Date(date_decimal(Time,tz = "EST")),y=Value,color=dataset),)+
  #geom_ribbon(aes(x=as.Date(date_decimal(Time,tz = "EST")),ymin=Lower,ymax=Upper,fill=dataset), alpha=0.3) +
  scale_fill_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
  scale_color_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
  geom_line(data = Re.dashboard, aes(x=as.Date(date),y=mean),color="#bababa")+
  geom_ribbon(data = Re.dashboard, aes(x=as.Date(date),ymin=lower,ymax=upper), alpha=0.3,fill="#bababa") +
  geom_hline(aes(yintercept=1),linetype = "dashed",color="grey")+
  labs(x = "Date", y = "Re")+
  #scale_fill_manual(values=c("#9970ab","#d7191c","#008837","#313695"))+
  scale_y_continuous(limits = c(0,10.2))+
  scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-02"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
        legend.position=c(0.85,0.6))

BA.2.bdsky

ggsave("../results/Re/BA.2.2_bdsky_skyline.pdf",BA.2.bdsky,width=10,height=6,units="cm")

# calculate mean Re
df.BA.2.2.mean <- df.BA.2.2 %>% 
  mutate(date = as.Date(date_decimal(Time,tz = "EST"))) %>% 
  filter(date > "2022-01-02") %>% 
  group_by(date) %>% 
  mutate(mean_value = mean(Value),
         mean_upper = mean(Upper),
         mean_lower = mean(Lower)) %>% 
  subset(,c("date","mean_value","mean_upper","mean_lower")) %>% 
  distinct()

write_tsv(df.BA.2.2.mean,"BA.2.2_mean_bdsky_skyline.tsv")

BA.2.bdsky.mean <- ggplot()+
  geom_line(data = df.BA.2.2.mean, aes(x=date,y=mean_value),color="#3cabd9")+
  geom_ribbon(data = df.BA.2.2.mean, aes(x=date,ymin=mean_lower,ymax=mean_upper), alpha=0.3,fill="#3cabd9") +
  geom_hline(aes(yintercept=1),linetype = "dashed",color="black")+
  labs(x = "Date", y = "Re")+
  scale_x_date(breaks="2 weeks", expand = c(0.007, 0), limits = c(as.Date("2022-01-01"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
        legend.position=c(0.85,0.6))
BA.2.bdsky.mean

ggsave("../results/Re/BA.2.2_mean_bdsky_skyline.pdf",BA.2.bdsky.mean,width=16.8,height=5,units="cm")
