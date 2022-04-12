library(data.table)
library(ggplot2)


mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
mdlmets$exp <- factor(mdlmets$exp, levels = c("fl1", "A","B","C",
                                                          "fl2","AB","AC","BC","fl3"))
mdlmets$forest_type <- factor(mdlmets$forest_type, levels = c("Riparian", "Broadleaf", "Coniferous", "Mixed"))

x <- mdlmets[variable=="R2"]
y <- x[value<0]$id
mdlmets <- mdlmets[!id%in%y]

##Result 2
cbPalette <- c("#999999", "#E69F00") , "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

bplots.ref.vst <- ggplot(data=mdlmets[variable!="RMSE"& 
                      Forest_attr=="Stem volume"& Metrics=="ref"], aes(x=exp, y=value, fill=fl))+
  geom_boxplot()+
  facet_wrap(forest_type~Forest_attr+variable, scales = "free", ncol = 3)+
  theme_minimal() + 
  scale_fill_manual(values= c("#56B4E9", "#E69F00", "#F0E442"))+
  labs(x="Experiment")+
  theme(text=element_text(family="serif", size=9*(96/72)),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none")

ggsave(bplots.ref.vst, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/bplots.ref.vst.svg", 
       width=15*1.25, height=23*1.25, units="cm", dpi=320) 






x <-mdlmets[forest_type=="Broadleaved"&
              Forest_attr=="Stem volume"&
              variable=="R2"&
              exp%in%c("fl1","fl2","fl3")&
              Metrics=="ref"]

x$exp <- as.character(x$exp)
df <- data.frame(x$exp, x$value)
df$x.exp <- as.factor(df$x.exp)

xx <- pairwise.var.test(df$x.value, df$x.exp, alternative = "greater")











bp.stats <- mdlmets[variable!="RMSE",.(Mean=mean(value),
                       SD=sd(value)), by=c("forest_type",
                                           "Forest_attr", 
                                           "variable",
                                           "exp",
                                           "Metrics")]


bp.stats2 <- dcast(bp.stats, exp~forest_type+Forest_attr+variable+Metrics, value.var =c("Mean","SD"))

write.csv(bp.stats2, "D:/1_Work/Dropbox/2_Publications/2_paper/results/bp_meanstats.csv")


x <- bp.stats[forest_type=="Broadleaved"]
x2 <- dcast[x, Metrics~Forest_attr, value.var=c("Mean")]

x <- melt(x, measure.vars = c("Mean", "SD"))

ggplot(data=x[variable!="RMSE"], aes(x=exp, y=value))+
  geom_line(aes(linetype=variable.1))+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_minimal()+
  scale_colour_grey(start = 0.4, end = 0.2)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")