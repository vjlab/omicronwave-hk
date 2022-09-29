library("ggpubr")
library("lme4")
library("readr")

metadata <- read_tsv("../data/18_districts_incidence_summary.tsv")

# spearman test
cor.test(metadata$incidence_rate, metadata$income,  method = "spearman")

cor.test(metadata$incidence_rate, metadata$population_density,  method = "spearman")

cor.test(metadata$imported_num, metadata$income,  method = "spearman")

cor.test(metadata$imported_num, metadata$population_density,  method = "spearman")


# plot the relationship between income, population density and incidence rate

g <- ggplot(data = metadata, aes(x=income,y=incidence_rate,fill=population_density,label=District))+
  geom_point(shape=21,size=3)+
  labs(x = "Median monthly household income (HK$)", y = "Incidence Rate (Covid cases/population)")+
  geom_text(hjust=-0.1, vjust=0.5,size=2)+
  scale_y_continuous(limits = c(4500,12500))+
  scale_x_continuous(limits = c(22000,45000))+
  scale_fill_gradientn(values=scales::rescale(c(5000, 10000, 20000, 40000, 60000)), colours=c("cyan", "grey", "red"))+
  #scale_fill_gradient2(low=muted("blue"),high=muted("red"),mid = "white")+
  geom_hline(aes(yintercept=9000),linetype = "dashed",color="grey")+
  geom_vline(aes(xintercept=35000),linetype = "dashed",color="grey")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9))
g

ggsave("../results/18 districts inequality/relationship.pdf",g,width=16,height=9,units="cm")



# plot the relationship between income, population density and imported cases

g1 <- ggplot(data = metadata, aes(x=income,y=imported_num,fill=population_density,label=District))+
  geom_point(shape=21,size=3)+
  labs(x = "Median monthly household income (HK$)", y = "Incidence Rate (Covid cases/population)")+
  geom_text(hjust=-0.1, vjust=0.5,size=2)+
  scale_y_continuous(limits = c(0,170))+
  scale_x_continuous(limits = c(22000,45000))+
  scale_fill_gradientn(values=scales::rescale(c(5000, 10000, 20000, 40000, 60000)), colours=c("cyan", "grey", "red"))+
  #scale_fill_gradient2(low=muted("blue"),high=muted("red"),mid = "white")+
  geom_hline(aes(yintercept=80),linetype = "dashed",color="grey")+
  geom_vline(aes(xintercept=35000),linetype = "dashed",color="grey")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 0.9, vjust = 0.9))
g1

ggsave("../results/18 districts inequality/relationship_imported.pdf",g1,width=16,height=9,units="cm")
