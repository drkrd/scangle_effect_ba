library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)
library(caret)
library(parallel)
library(foreach)
library(doParallel)

mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                            "Broadleaved", "Mixed"))
mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                            "fl2", "AB", "AC", "BC", "fl3")) 



mdlmets.smry <- mdlmets[variable!="RMSE",.(mn=mean(value),
                           sdev=sd(value)), by=c("forest_type", "Forest_attr", "variable", "Metrics")]

mdlmets.smry <- dcast(mdlmets.smry, forest_type~Forest_attr+variable+Metrics, value.var = c("mn", "sdev"))
write.csv(mdlmets.smry, "D:/1_Work/Dropbox/2_Publications/2_paper/results/bp.all.meansdstats.csv")




#Result 2
bplots.ref <- ggplot(data=mdlmets[variable!="RMSE"& 
                                    exp %in% c("fl1", "fl2", "fl3")&
                                    Metrics=="ref"&
                                    Forest_attr!="Basal area"], 
                     aes(x=Forest_attr, y=value, fill=exp))+
  geom_boxplot()+
  theme_few()+
  facet_wrap(forest_type~variable, scales = "free", ncol = 3)+
  labs(x="Forest attribute", y="")+
  scale_fill_manual(values= c("#56B4E9", "#E69F00", "#F0E442"))+
  scale_x_discrete(labels= c("BA", "Vst"))+
  theme(text=element_text(family="serif", size=9*(96/72)),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")

ggsave(bplots.ref, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/bplots.ref.all.final.svg", 
       width=15*1.25, height=23*1.25, units="cm", dpi=640) 


crn.bp.stats <- ciron.mdlmets[,.(mean=mean(value),
                                 sd=sd(value)), by=c("Forest_attr", 
                                                     "variable",
                                                     "exp",
                                                     "Metrics")]


bplots.ref <- ggplot(data=mdlmets[variable!="RMSE"& 
                                    exp %in% c("A", "B", "C", "AB", "AC", "BC")&
                                    Metrics=="ref"&
                                    Forest_attr=="Basal area"], 
                     aes(x=fl, y=value, fill=fl))+
  geom_boxplot()+
  theme_few()+
  facet_wrap(forest_type~variable, scales = "free", ncol = 3)+
  labs(x="Forest attribute", y="")+
  scale_fill_manual(values= c("#56B4E9", "#E69F00"))+
  scale_x_discrete(labels= c("(A+B+C)", "(AB+AC+BC)"))+
  theme(text=element_text(family="serif", size=9*(96/72)),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")







#Result 2 appendix
bplots.ref.ba <- ggplot(data=mdlmets[variable!="RMSE" & 
                                    Forest_attr=="Basal area" &
                                    fl %in% c("fl1", "fl2", "fl3") &
                                    Metrics=="ref"], 
                     aes(fill=fl, y=value, x=exp))+
  geom_boxplot()+
  theme_few()+
  facet_wrap(forest_type~Forest_attr+variable, scales = "free", ncol=3)+
  labs(x="Experiments", y="")+
  scale_fill_manual(values= c("#56B4E9", "#E69F00", "#F0E442"))+
  annotate(xmin=4.5, xmax=8.5, ymin=-Inf, ymax=Inf, geom = "rect", alpha = 0.2 )+
  theme(text=element_text(family="serif", size=9*(96/72)),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")

ggsave(bplots.ref.ba, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/bplots.ref.ba.svg", 
       width=15*1.25, height=23*1.25, units="cm", dpi=640) 


##Result 3
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

bplots.refvoxall.ba <- ggplot(data=mdlmets[variable!="RMSE"& 
                                          exp %in% c("fl1", "fl2", "fl3")&
                                            Forest_attr!="Total volume"], 
                           aes(x=Forest_attr, y=value, fill=Metrics))+
  geom_boxplot()+
  facet_wrap(forest_type~variable, scales = "free", ncol = 3)+
  scale_fill_manual(values=c("red4", "grey"))+
  theme_few() + 
  labs(x="Experiment")+
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text=element_text(family="serif", size=9*(96/72)))+
  scale_x_discrete(labels= c("BA", "Vst", "Vtot"))

ggsave(bplots.refvoxall.ba, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/bplots.refvoxall.final.svg", 
       width=15*1.25, height=23*1.25, units="cm", dpi=640) 



##Result 3 appendix
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

bplots.refvox.ba <- ggplot(data=mdlmets[variable!="RMSE"& 
                                        Forest_attr=="Basal area"], 
                         aes(x=exp, y=value, fill=fl, colour=Metrics))+
  geom_boxplot(outlier.size = 1.5)+
  facet_wrap(forest_type~Forest_attr~variable, scales = "free", ncol = 3)+
  theme_few() + 
  scale_fill_manual(values= c("#56B4E9", "#E69F00", "#F0E442"))+
  scale_colour_manual(values=c("grey60", "black"))+
  labs(x="Experiment")+
  annotate(xmin=4.5, xmax=8.5, ymin=-Inf, ymax=Inf, geom = "rect", alpha = 0.2)+
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text=element_text(family="serif", size=9*(96/72)))

ggsave(bplots.refvox.ba, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/bplots.refvox.ba.svg", 
       width=15*1.25, height=23*1.25, units="cm", dpi=640) 



