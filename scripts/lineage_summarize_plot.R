library(ggplot2)
library(tidyverse)
library(readr)
library(lubridate)
library(scales)
library(HDInterval)
library(ggsci)

data.21A <- read.csv("../results/monophyletic clade summary/Supplementary Data 2-Delta_MCC_summary.csv") %>% mutate(Lineage="Delta")
data.21K <- read.csv("../results/monophyletic clade summary/Supplementary Data 1-Omicron_BA.1_MCC_summary.csv") %>% mutate(Lineage="Omicron:BA.1.*")
data.21L <- read.csv("../results/monophyletic clade summary/Supplementary Data 3-Omicron_BA.2_MCC_summary.csv") %>% mutate(Lineage="Omicron:BA.2.*")

combine.data <- rbind(data.21A,data.21K,data.21L) %>% mutate(index = rownames(.))

g.lineage <- ggplot(combine.data) +
  geom_jitter(aes(y = num_samples, x = as.Date(first_sample_date), size = num_samples, fill = Lineage),
              alpha = 0.8,
              shape = 21, 
              #fill = "#69a2d5",
              height = 0.01) +
  scale_fill_manual(values=c("#b72467","#f68b1f","#00addc"))+
  scale_y_continuous("Samples per lineage",
                     expand = c(0.01,0), 
                     limits = c(NA,10000),
                     trans = "log10") +
  scale_x_date("Date of lineage detection",
               date_breaks = "1 month", 
               date_labels = "%b",
               expand = c(0.1,0.1)) +
  theme_classic() +
  theme(legend.position = c(0.9,0.6)) +
  annotation_logticks() +
  labs(size = "No. samples")

g.lineage

# add duration
duration.data <- combine.data %>% filter(num_samples > 1)
detected_duration <- rbind(duration.data[,c("last_sample_date","num_samples","Lineage","index")] %>% rename(,x_date = last_sample_date),
                           duration.data[,c("first_sample_date","num_samples","Lineage","index")] %>% rename(,x_date = first_sample_date)) %>% 
  mutate(Duration = "")

g.duration <- g.lineage+
  geom_line(data=detected_duration, aes(y = num_samples, x = as.Date(x_date),group = index,color=Lineage, alpha = Duration))+
  scale_color_manual(values=c("#b72467","#f68b1f","#00addc"))+
  scale_alpha_manual(values=c(1,0.4))+
  theme(legend.position = "none",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))
g.duration

ggsave("../results/monophyletic clade summary/lineage_samples_duration.pdf",g.duration,width=10,height=6,units="cm")


# delay and duration
dd.input <- combine.data %>% 
  filter(num_samples > 2) %>% 
  mutate(duration = known_duration+undetected_duration)

g.dd <- ggplot(dd.input) +
  geom_jitter(aes(x = undetected_duration, y = duration, fill = Lineage, size = num_samples), 
              #height = 5, 
              #width = 5, 
              shape = 21,
              alpha = 0.7) +
  geom_smooth(method = lm, aes(x = undetected_duration, y = duration),
              se = T, 
              size = 0.7, 
              colour = "black", 
              alpha = 0.3) +
  scale_x_continuous("Detection lag (days)",
                     limits = c(0,90),
                     breaks = seq(0,90, by = 10),
                     expand = c(0, 0), 
                     oob = squish) +
  scale_y_continuous("Lineage duration (days)",
                     limits = c(0,150),
                     breaks = seq(0,150, by = 30),
                     expand = c(0, 0), 
                     oob = squish) +
  theme_classic() +
  theme(legend.position = "none",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8))+
  #theme(aspect.ratio = 1, 
  #      legend.position = c(0.85,0.85)) +
  labs(size = "Samples per lineage", 
       fill = "") +
  scale_fill_manual(values = c("#b72467","#f68b1f","#00addc"))

g.dd

mean(dd.input$undetected_duration)
median(dd.input$undetected_duration)
HDInterval::hdi(dd.input$undetected_duration)
cor.test(dd.input$undetected_duration, dd.input$duration,
         method = "spearman",
         exact = F)

ggsave("../results/monophyletic clade summary/lineage_delay_duration_correlation.pdf",g.dd,width=10,height=6,units="cm")
