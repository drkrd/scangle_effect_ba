





crn.bp.stats.ref <- ciron.mdlmets.ref[variable!="RMSE",.(Mean=mean(value),
                                         SD=sd(value)), by=c("Forest_attr", 
                                                             "variable",
                                                             "exp",
                                                             "Metrics")]
crn.bp.stats.ref <- dcast(crn.bp.stats.ref, exp~variable+Metrics+Forest_attr, value.var =c("Mean","SD"))
write.csv(crn.bp.stats.ref, "D:/1_Work/Dropbox/2_Publications/2_paper/results/crn_bp_meanstats_ref.csv")


crn.bp.stats <- ciron.mdlmets[,.(Mean=mean(value),
                                 SD=sd(value)), by=c("Forest_attr", 
                                                     "variable",
                                                     "exp",
                                                     "Metrics")]


crn.bp.stats <- dcast(crn.bp.stats, exp~Forest_attr+variable+Metrics, value.var =c("Mean","SD"))



write.csv(crn.bp.stats, "D:/1_Work/Dropbox/2_Publications/2_paper/results/crn_bp_meanstats.csv")
crn.bp.stats <- melt(crn.bp.stats, id.vars = c(1:4), variable.name = "stat")

crn.bp.stats1 <- crn.bp.stats[Forest_attr=="Basal area" & variable=="RMSE"]

crn.bp.stats1 <- crn.bp.stats[, value=ifelse(variable =="R2" & stat=="Mean", mean*100 , mean)]

brplot.refvox <- ggplot(data=crn.bp.stats[variable!="RMSE"], 
                     aes(x=exp, y=value, fill=Forest_attr, colour=Metrics))+
  geom_bar(stat="Identity", position=position_dodge2())+
  facet_wrap(stat~variable, scales = "free", ncol=3)+
  theme_minimal()+
  scale_fill_grey(start=0.8, end=0)+

  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")+
  labs(x="Experiments", y="Value", fill="Forest attribute")

ggsave(brplot.refvox, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/barplot_refvox.png", 
       width=24*1.25, height=16*1.25, units="cm", dpi=320) 








old <- -4.91
new <- -4.14

diff <- new-old
pc <- 100*(diff/old)
pc











tests.ref <- ciron.mdlmets[Forest_attr=="Basal area" & variable=="MPE" & exp%in%c("fl1","fl2","fl3") &Metrics=="ref", c("value", "exp") ]
# tests.ref <- tests.ref[exp!="C"]
grp <- as.character(tests.ref$exp)
value <- tests.ref$value
tests.ref <- data.frame(grp, value)
tests.ref$grp <- as.factor(tests.ref$grp)
oneway.test(value~grp, data = tests.ref)
xy <- games_howell_test(tests.ref, value~grp , conf.level = 0.99, detailed = FALSE)
xy$x <- paste0(xy$group1,'-',xy$group2)
setDT(xy)
xy <- xy[order(rank(estimate))]
lvls <- xy$x
xy$x <- factor(xy$x, levels=lvls)
ggplot(xy, aes(x= x, y=estimate)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high))+
  theme(axis.text.x = element_text(angle = 45))+
  geom_hline(yintercept = 0)
sd(xy$estimate)

tests.vox <- ciron.mdlmets[Forest_attr=="Stem volume" & variable=="R2" & fl%in%c("fl1","fl2","fl3") &Metrics=="vox", c("value", "exp") ]
# tests.vox <- tests.vox[exp!="C"]
grp <- as.character(tests.vox$exp)
value <- tests.vox$value
tests.vox <- data.frame(grp, value)
tests.vox$grp <- as.factor(tests.vox$grp)
oneway.test(value~grp, data = tests.vox)
xz <- games_howell_test(tests.vox, value~grp , conf.level = 0.99, detailed = FALSE)
xz$x <- paste0(xz$group1,'-',xz$group2)
setDT(xz)
xz <- xz[order(rank(estimate))]
lvls <- xz$x
xz$x <- factor(xz$x, levels=lvls)
ggplot(xz, aes( x= x, y=estimate)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  theme_minimal()+
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high))+
  theme(axis.text.x = element_text(angle = 90))+
  geom_hline(yintercept = 0)

sd(xz$estimate)


setkey(xy,"x")
setkey(xz,"x")
xyz <- xy[xz]


ggplot(data=xyz, aes(x=estimate, y=i.estimate))+
  geom_label(aes(label=x))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point()




