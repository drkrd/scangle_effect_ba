library(ggplot2)
library(ggsci)
library(data.table)
library(rstatix)





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







mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)


tst <- NULL


tests.ref <- ciron.mdlmets[Forest_attr=="Basal area" & 
                             variable=="MPE"  &
                             Metrics=="ref", c("value", "exp") ]
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
xy$fr_type <- "crn"
xy$fr_at <- "G"
xy$var <- "MPE"
ggplot(xy, aes(x= x, y=estimate)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high))+
  theme(axis.text.x = element_text(angle = 45))+
  geom_hline(yintercept = 0)
sd(xy$estimate)
tst <- rbind(tst, xy)















tst2.vox <- NULL
forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests.vox <- mdlmets[forest_type== forest&
                                    Forest_attr==forattr & 
                             fl %in% c("fl1")&
                                   variable==msr  &
                                   Metrics=="vox", c("value", "exp")]
      # tests.vox <- tests.vox[exp!="C"]
      grp <- as.character(tests.vox$exp)
      value <- tests.vox$value
      tests.vox <- data.frame(grp, value)
      tests.vox$grp <- as.factor(tests.vox$grp)
      oneway.test(value~grp, data = tests.vox)
      xy <- games_howell_test(tests.vox, value~grp , conf.level = 0.99, detailed = FALSE)
      xy$x <- paste0(xy$group1,'-',xy$group2)
      setDT(xy)
      xy <- xy[order(rank(estimate))]
      lvls <- xy$x
      xy$x <- factor(xy$x, levels=lvls)
      xy$fr_type <- forest
      xy$fr_at <- forattr
      xy$var <- msr
      tst2.vox <- rbind(tst2.vox, xy)
    }
  }
}
tst2.vox <- setDT(tst2.vox)
tst2.vox.ns <- tst2.vox[p.adj.signif=="ns"]
tst2.vox <- tst2.vox[,clr:= ifelse(p.adj.signif=="ns", "yes", "no")]


tst2.ref.smry <- tst2.ref[, .(min=min(abs(estimate)),
                              mn=mean(abs(estimate)),
                              max=max(abs(estimate))), by=.(fr_type, fr_at, var)]

tst2.ref.smry$Metrics <- "ref"

tst2.vox.smry <- tst2.vox[, .(min=min(abs(estimate)),
                              mn=mean(abs(estimate)),
                              max=max(abs(estimate))), by=.(fr_type, fr_at, var)]


tst2.vox.smry$Metrics <- "vox"

tst2.smry <- rbind(tst2.ref.smry, tst2.vox.smry)

tst2.smry$var <- factor(tst2.smry$var, levels=c("R2", "rRMSE", "MPE"))

ggplot(data=tst2.smry, aes(x=fr_type, y=mn, colour=Metrics ))+
  geom_point()+
  facet_grid(var~fr_at, scales = "free")+
  geom_errorbar(aes(ymin = min, ymax = max), width=0.3)+
  theme_minimal()
  




tst2.ref <- NULL
tst2.vox <- NULL
forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests.ref <- mdlmets[forest_type== forest&
                             Forest_attr==forattr &
                             variable==msr  &
                             Metrics=="ref", c("value", "exp")]
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
      xy$fr_type <- forest
      xy$fr_at <- forattr
      xy$var <- msr
      tst2.ref <- rbind(tst2.ref, xy)
      tst2.ref.smry <- tst2.ref[, .(min=min(abs(estimate)),
                                    mn=mean(abs(estimate)),
                                    max=max(abs(estimate))), by=.(fr_type, fr_at, var)]
      tst2.ref.smry$Metrics <- "ref"
      
      
      tests.vox <- mdlmets[forest_type== forest&
                             Forest_attr==forattr &
                             variable==msr  &
                             Metrics=="vox", c("value", "exp")]
      # tests.vox <- tests.vox[exp!="C"]
      grp <- as.character(tests.vox$exp)
      value <- tests.vox$value
      tests.vox <- data.frame(grp, value)
      tests.vox$grp <- as.factor(tests.vox$grp)
      oneway.test(value~grp, data = tests.vox)
      xy <- games_howell_test(tests.vox, value~grp , conf.level = 0.99, detailed = FALSE)
      xy$x <- paste0(xy$group1,'-',xy$group2)
      setDT(xy)
      xy <- xy[order(rank(estimate))]
      lvls <- xy$x
      xy$x <- factor(xy$x, levels=lvls)
      xy$fr_type <- forest
      xy$fr_at <- forattr
      xy$var <- msr
      tst2.vox <- rbind(tst2.vox, xy)
      tst2.vox.smry <- tst2.vox[, .(min=min(abs(estimate)),
                                    mn=mean(abs(estimate)),
                                    max=max(abs(estimate))), by=.(fr_type, fr_at, var)]
      tst2.vox.smry$Metrics <- "vox"
      
    }
  }
}
tst2.smry <- rbind(tst2.ref.smry, tst2.vox.smry)


tst2.ref <- setDT(tst2.ref)
tst2.ref.ns <- tst2.ref[p.adj.signif=="ns"]
tst2.ref <- tst2.ref[,clr:= ifelse(p.adj.signif=="ns", "yes", "no")]






tst2.smry$var <- factor(tst2.smry$var, levels=c("R2", "rRMSE", "MPE"))

ggplot(data=tst2.smry, aes(x=fr_type, y=mn, colour=Metrics ))+
  geom_point()+
  facet_grid(var~fr_at, scales = "free")+
  geom_errorbar(aes(ymin = min, ymax = max), width=0.3)+
  theme_minimal()


mdlmets <- setDT(mdlmets)

xx <- mdlmets[forest_type=="Mixed"&
                variable=="R2"&
              Forest_attr=="Total volume"&
              Metrics=="ref"]


  t <- pairwise.var.test(xx$value, xx$exp, p.method = "fdr", "two.sided")









ggsave(x1, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/gh.test.mixed.png", 
       width=24*1.25, height=16*1.25, units="cm", dpi=320) 







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




