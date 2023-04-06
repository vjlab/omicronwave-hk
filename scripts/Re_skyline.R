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

#unifrom 20 v2
lf.BA.2.2 <- readLogfile("xml_edited/HK_BA_2_2_un_20_v1.fasta.edited.v2.r3.log", burnin=0.1)

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

df.BA.2.2.un.20.v2 <- as.data.frame(t(as.data.frame(Re_gridded_hpd.BA.2.2)))

rownames(df.BA.2.2.un.20.v2) <- NULL
colnames(df.BA.2.2.un.20.v2) <- c("Lower","Value","Upper")

df.BA.2.2.un.20.v2$Time <- times.BA.2.2
df.BA.2.2.un.20.v2$days <- (treeheight-timegrid.BA.2.2)*365
df.BA.2.2.un.20.v2$dataset <- "un.20"
df.BA.2.2.un.20.v2$sampling_proportion <- median(lf.BA.2.2$samplingProportion_BDSKY_Serial)


#unifrom 40 v2
lf.BA.2.2 <- readLogfile("xml_edited/HK_BA_2_2_un_40.fasta.edited.final.v2.0.1.log", burnin=0)

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

df.BA.2.2.un.40.v2 <- as.data.frame(t(as.data.frame(Re_gridded_hpd.BA.2.2)))

rownames(df.BA.2.2.un.40.v2) <- NULL
colnames(df.BA.2.2.un.40.v2) <- c("Lower","Value","Upper")

df.BA.2.2.un.40.v2$Time <- times.BA.2.2
df.BA.2.2.un.40.v2$days <- (treeheight-timegrid.BA.2.2)*365
df.BA.2.2.un.40.v2$dataset <- "un.40"
df.BA.2.2.un.40.v2$sampling_proportion <- median(lf.BA.2.2$samplingProportion_BDSKY_Serial)


#proportional 
lf.BA.2.2 <- readLogfile("xml_edited/HK_BA_2_2_pro_v2.fasta.edited.final.315.0.1.burnin.log", burnin=0)

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

df.BA.2.2.pro <- as.data.frame(t(as.data.frame(Re_gridded_hpd.BA.2.2)))

rownames(df.BA.2.2.pro) <- NULL
colnames(df.BA.2.2.pro) <- c("Lower","Value","Upper")

df.BA.2.2.pro$Time <- times.BA.2.2
df.BA.2.2.pro$days <- (treeheight-timegrid.BA.2.2)*365
df.BA.2.2.pro$dataset <- "pro"
df.BA.2.2.pro$sampling_proportion <- median(lf.BA.2.2$samplingProportion_BDSKY_Serial)


df.BA.2.2.Re <- rbind(df.BA.2.2.un.20.v2,df.BA.2.2.un.40.v2,df.BA.2.2.pro)


BA.2.bdsky <- ggplot(df.BA.2.2.Re)+
  geom_line(aes(x=as.Date(date_decimal(Time,tz = "EST")),y=Value,color=dataset),)+
  #geom_ribbon(aes(x=as.Date(date_decimal(Time,tz = "EST")),ymin=Lower,ymax=Upper,fill=dataset), alpha=0.3) +
  #scale_fill_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
  scale_color_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
  geom_line(data = Re.dashboard, aes(x=as.Date(date),y=mean),color="#bababa")+
  geom_ribbon(data = Re.dashboard, aes(x=as.Date(date),ymin=lower,ymax=upper), alpha=0.3,fill="#bababa") +
  geom_hline(aes(yintercept=1),linetype = "dashed",color="grey")+
  labs(x = "Date", y = "Re", color = "Dataset")+
  #scale_fill_manual(values=c("#9970ab","#d7191c","#008837","#313695"))+
  scale_y_continuous(limits = c(0,10))+
  scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-02"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
        legend.position=c(0.85,0.6))

BA.2.bdsky

ggsave("figures/BA.2_bdskys_skyline_different_dataset.pdf",BA.2.bdsky,width=14,height=8,units="cm")


#un.20
(BA.2.bdsky.un.20 <- ggplot(df.BA.2.2.un.20.v2)+
    geom_line(aes(x=as.Date(date_decimal(Time,tz = "EST")),y=Value),color="#00addc")+
    geom_ribbon(aes(x=as.Date(date_decimal(Time,tz = "EST")),ymin=Lower,ymax=Upper), alpha=0.3,fill="#00addc") +
    #scale_fill_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
    scale_color_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
    geom_line(data = Re.dashboard, aes(x=as.Date(date),y=mean),color="#bababa")+
    geom_ribbon(data = Re.dashboard, aes(x=as.Date(date),ymin=lower,ymax=upper), alpha=0.3,fill="#bababa") +
    geom_hline(aes(yintercept=1),linetype = "dashed",color="grey")+
    labs(x = "Date", y = "Re", color = "Dataset")+
    #scale_fill_manual(values=c("#9970ab","#d7191c","#008837","#313695"))+
    scale_y_continuous(limits = c(0,10))+
    scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-02"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
          legend.position=c(0.85,0.6)))

ggsave("figures/BA.2_bdskys_skyline_un.20.pdf",width=10,height=6,units="cm")

#un.40
(BA.2.bdsky.un.40 <- ggplot(df.BA.2.2.un.40.v2)+
    geom_line(aes(x=as.Date(date_decimal(Time,tz = "EST")),y=Value),color="#00addc")+
    geom_ribbon(aes(x=as.Date(date_decimal(Time,tz = "EST")),ymin=Lower,ymax=Upper), alpha=0.3,fill="#00addc") +
    #scale_fill_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
    scale_color_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
    geom_line(data = Re.dashboard, aes(x=as.Date(date),y=mean),color="#bababa")+
    geom_ribbon(data = Re.dashboard, aes(x=as.Date(date),ymin=lower,ymax=upper), alpha=0.3,fill="#bababa") +
    geom_hline(aes(yintercept=1),linetype = "dashed",color="grey")+
    labs(x = "Date", y = "Re", color = "Dataset")+
    #scale_fill_manual(values=c("#9970ab","#d7191c","#008837","#313695"))+
    scale_y_continuous(limits = c(0,10))+
    scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-02"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
          legend.position=c(0.85,0.6)))

#ggsave("figures/BA.2_bdskys_skyline_un.40.pdf",width=10,height=6,units="cm")
ggsave("figures/BA.2_bdskys_skyline_un.40_final.pdf",width=16.8,height=5,units="cm")

#proportional
(BA.2.bdsky.pro <- ggplot(df.BA.2.2.pro)+
    geom_line(aes(x=as.Date(date_decimal(Time,tz = "EST")),y=Value),color="#00addc")+
    geom_ribbon(aes(x=as.Date(date_decimal(Time,tz = "EST")),ymin=Lower,ymax=Upper), alpha=0.3,fill="#00addc") +
    #scale_fill_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
    scale_color_manual(values=c("#ed930d","#43a459","#60388b","#e7298a"))+
    geom_line(data = Re.dashboard, aes(x=as.Date(date),y=mean),color="#bababa")+
    geom_ribbon(data = Re.dashboard, aes(x=as.Date(date),ymin=lower,ymax=upper), alpha=0.3,fill="#bababa") +
    geom_hline(aes(yintercept=1),linetype = "dashed",color="grey")+
    labs(x = "Date", y = "Re", color = "Dataset")+
    #scale_fill_manual(values=c("#9970ab","#d7191c","#008837","#313695"))+
    scale_y_continuous(limits = c(0,10))+
    scale_x_date(breaks="2 weeks", limits = c(as.Date("2022-01-02"),as.Date("2022-04-26")),date_labels = "%d-%b",guide = "axis_minor")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9),
          legend.position=c(0.85,0.6)))

ggsave("figures/BA.2_bdskys_skyline_pro.pdf",width=10,height=6,units="cm")
