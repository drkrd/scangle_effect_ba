bgs.bp.stats.con <- bauges.mdlmets.con[,.(mean=mean(value),
                                      sd=sd(value)), by=c("Forest_attr", 
                                                          "variable",
                                                          "exp",
                                                          "Metrics")]
bgs.bp.stats.con <- melt(bgs.bp.stats.con, id.vars = c(1:4), variable.name = "stat")

brplot.ref.con <- ggplot(data=bgs.bp.stats.con[variable!="RMSE"&Metrics=="ref"], aes(x=exp, y=value, fill=Forest_attr))+
  geom_bar(stat="Identity", position=position_dodge())+
  facet_wrap(stat~variable, scales = "free")+
  theme_minimal()+
  scale_fill_grey(start = 0.7, end = 0.1)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")+
  labs(x="Experiments", y="Value", fill="Forest attribute")
####################################################
bgs.bp.stats.feu <- bauges.mdlmets.feu[,.(mean=mean(value),
                                          sd=sd(value)), by=c("Forest_attr", 
                                                              "variable",
                                                              "exp",
                                                              "Metrics")]
bgs.bp.stats.feu <- melt(bgs.bp.stats.feu, id.vars = c(1:4), variable.name = "stat")

brplot.ref.feu <- ggplot(data=bgs.bp.stats.feu[variable!="RMSE"&Metrics=="vox"], aes(x=exp, y=value, fill=Forest_attr))+
  geom_bar(stat="Identity", position=position_dodge())+
  facet_wrap(stat~variable, scales = "free")+
  theme_minimal()+
  scale_fill_grey(start = 0.7, end = 0.1)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")+
  labs(x="Experiments", y="Value", fill="Forest attribute")

#########################################################################
bgs.bp.stats.mix <- bauges.mdlmets.mix[,.(mean=mean(value),
                                          sd=sd(value)), by=c("Forest_attr", 
                                                              "variable",
                                                              "exp",
                                                              "Metrics")]
bgs.bp.stats.mix <- melt(bgs.bp.stats.mix, id.vars = c(1:4), variable.name = "stat")

brplot.ref.mix <- ggplot(data=bgs.bp.stats.mix[variable!="RMSE"&Metrics=="ref"], aes(x=exp, y=value, fill=Forest_attr))+
  geom_bar(stat="Identity", position=position_dodge())+
  facet_wrap(stat~variable, scales = "free")+
  theme_minimal()+
  scale_fill_grey(start = 0.7, end = 0.1)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")+
  labs(x="Experiments", y="Value", fill="Forest attribute")