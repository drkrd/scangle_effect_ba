library(ggplot2)
library(ggpubr)
library(grid)
library(gtable) 
library(ggsci)
library(data.table)
library(rstatix)
library(ggthemes)
library(RVAideMemoire)


mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]



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
#################################################################################################################

forest.type <- c("Mixed", "Coniferous", "Broadleaf", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
exps <- c("fl1", "A", "B", "C", "fl2", "AB", "AC", "BC", "fl3")
ghtests <- NULL
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests <- mdlmets[forest_type== forest&
                         Forest_attr==forattr & 
                         fl %in% c("fl1", "fl2", "fl3")&
                         Metrics=="ref" &
                         variable==msr  &
                         exp %in% c("fl1", "fl2", "fl3"), c("value", "exp")]
      grp <- as.character(tests$exp)
      value <- tests$value
      tests <- data.frame(grp, value)
      tests$grp <- as.factor(tests$grp)
      t1 <- oneway.test(value~grp, data = tests)
      t2 <- games_howell_test(tests, value~grp , conf.level = 0.95, detailed = FALSE)
      t2$forest_type <- forest
      t2$forest_attr <- forattr
      t2$msr <- msr
      ghtests <- rbind(ghtests, t2)
    }}}


setDT(ghtests)
ghtests1 <- ghtests[forest_attr=="Stem volume", c("group1", "group2", "estimate", "p.adj.signif", "forest_attr", "msr", "forest_type")]


setDT(ghtests)
ghtests <- ghtests[,pair:=paste0(group1,"-",group2)]
ggplot(ghtests[msr=="R2"&forest_attr=="Stem volume"])+
  aes(x=pair, y=estimate,colour=forest_type, shape=forest_type, group=forest_type)+
  geom_point()+
  geom_line()



forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
exps <- c("fl1", "A", "B", "C", "fl2", "AB", "AC", "BC", "fl3")
ghtestsv <- NULL
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests <- mdlmets[forest_type== forest&
                         Forest_attr==forattr & 
                         fl %in% c("fl1", "fl2", "fl3")&
                         Metrics=="vox" &
                         variable==msr  &
                         exp %in% c("A", "B", "C", "AB", "AC", "BC"), c("value", "exp")]
      grp <- as.character(tests$exp)
      value <- tests$value
      tests <- data.frame(grp, value)
      tests$grp <- as.factor(tests$grp)
      t1 <- oneway.test(value~grp, data = tests)
      t2 <- games_howell_test(tests, value~grp , conf.level = 0.95, detailed = FALSE)
      t2$forest_type <- forest
      t2$forest_attr <- forattr
      t2$msr <- msr
      ghtestsv <- rbind(ghtestsv, t2)
    }}}

setDT(ghtestsv)
ghtestsv <- ghtestsv[,pair:=paste0(group1,"-",group2)]

ggplot(ghtestsv[msr=="R2"&forest_attr=="Stem volume"])+
  aes(x=pair, y=estimate,colour=forest_type, shape=forest_type, group=forest_type)+
  geom_point()+
  geom_line()



ghtests$metrics <- "ref"
ghtestsv$metrics <- "vox"
ghtestsall <- rbind(ghtests, ghtestsv)

ghtestsall$pair <- factor(ghtestsall$pair, c("A-B", "A-AB", "A-C", "B-C", "AB-AC", "AB-BC", "AC-BC",
                                          "A-AB", "A-AC", "A-BC", "AB-B", "AC-B", "B-BC",
                                          "AB-C", "AC-C", "BC-C"))

ggplot(ghtestsall[msr=="R2"&forest_attr=="Basal area"])+
  aes(x=pair, y=abs(estimate),colour=metrics, shape=metrics, group=metrics)+
  geom_point()+
  geom_line()+
  facet_grid(forest_type~., scales = "free")





#############################################################################################################################


mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]


mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                            "Broadleaf", "Mixed"))
mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                            "fl2", "AB", "AC", "BC", "fl3")) 

{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Broadleaf" & variable=="R2" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  matrix <- matrix(ncol=6, nrow=6)
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
        matrix[i,j] <-  t2.stk[group1==i&group2==j]$estimate
    }
  }
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Riparian" & variable=="R2" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in 1:6){
    for(j in 1:6){
      if(i<j)
        matrix[i,j] <-  t2.stk[group1==rownames(matrix)[i]&group2==colnames(matrix)[j]]$estimate
    }
  }
  
  
  dt <- as.data.table(matrix)
  dt <- cbind(colnames(dt), dt)
  dt <- melt(dt, measure.vars = c(2,3,4,5,6,7))
  dt$value <- as.numeric(round(dt$value, 2))
  dt$variable <- factor(dt$variable, c("BC", "AC", "AB", "C", "B", "A"))
  dt$V1 <- factor(dt$V1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  mat_plt_rip_bro_r2 <- ggplot(dt)+aes(x=V1, y=variable)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5,linetype = 1)+
    scale_x_discrete(position = "top") +
    scale_colour_manual("white")+
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value), color = "black", size=3)+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}

{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Broadleaf" & variable=="rRMSE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  matrix <- matrix(ncol=6, nrow=6)
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
        matrix[i,j] <-  t2.stk[group1==i&group2==j]$estimate
    }
  }
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Riparian" & variable=="rRMSE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in 1:6){
    for(j in 1:6){
      if(i<j)
        matrix[i,j] <-  t2.stk[group1==rownames(matrix)[i]&group2==colnames(matrix)[j]]$estimate
    }
  }
  
  
  dt <- as.data.table(matrix)
  dt <- cbind(colnames(dt), dt)
  dt <- melt(dt, measure.vars = c(2,3,4,5,6,7))
  dt$value <- as.numeric(round(dt$value, 1))
  dt$variable <- factor(dt$variable, c("BC", "AC", "AB", "C", "B", "A"))
  dt$V1 <- factor(dt$V1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  mat_plt_rip_bro_rmse <- ggplot(dt)+aes(x=V1, y=variable)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5,linetype = 1)+
    scale_x_discrete(position = "top") +
    scale_colour_manual("white")+
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value), color = "black", size=3)+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}

{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Broadleaf" & variable=="MPE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  matrix <- matrix(ncol=6, nrow=6)
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
        matrix[i,j] <-  t2.stk[group1==i&group2==j]$estimate
    }
  }
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Riparian" & variable=="MPE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in 1:6){
    for(j in 1:6){
      if(i<j)
        matrix[i,j] <-  t2.stk[group1==rownames(matrix)[i]&group2==colnames(matrix)[j]]$estimate
    }
  }
  
  
  dt <- as.data.table(matrix)
  dt <- cbind(colnames(dt), dt)
  dt <- melt(dt, measure.vars = c(2,3,4,5,6,7))
  dt$value <- as.numeric(round(dt$value, 1))
  dt$variable <- factor(dt$variable, c("BC", "AC", "AB", "C", "B", "A"))
  dt$V1 <- factor(dt$V1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  mat_plt_rip_bro_mpe <- ggplot(dt)+aes(x=V1, y=variable)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5,linetype = 1)+
    scale_x_discrete(position = "top") +
    scale_colour_manual("white")+
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value), color = "black", size=3)+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}

{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Mixed" & variable=="R2" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  matrix <- matrix(ncol=6, nrow=6)
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
        matrix[i,j] <-  t2.stk[group1==i&group2==j]$estimate
    }
  }
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Coniferous" & variable=="R2" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in 1:6){
    for(j in 1:6){
      if(i<j)
        matrix[i,j] <-  t2.stk[group1==rownames(matrix)[i]&group2==colnames(matrix)[j]]$estimate
    }
  }
  
  
  dt <- as.data.table(matrix)
  dt <- cbind(colnames(dt), dt)
  dt <- melt(dt, measure.vars = c(2,3,4,5,6,7))
  dt$value <- as.numeric(round(dt$value, 2))
  dt$variable <- factor(dt$variable, c("BC", "AC", "AB", "C", "B", "A"))
  dt$V1 <- factor(dt$V1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  mat_plt_con_mix_r2 <- ggplot(dt)+aes(x=V1, y=variable)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5,linetype = 1)+
    scale_x_discrete(position = "top") +
    scale_colour_manual("white")+
    labs(x="Mixed", y="Coniferous")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value), color = "black", size=3)+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}

{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Mixed" & variable=="rRMSE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  matrix <- matrix(ncol=6, nrow=6)
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
        matrix[i,j] <-  t2.stk[group1==i&group2==j]$estimate
    }
  }
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Coniferous" & variable=="rRMSE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in 1:6){
    for(j in 1:6){
      if(i<j)
        matrix[i,j] <-  t2.stk[group1==rownames(matrix)[i]&group2==colnames(matrix)[j]]$estimate
    }
  }
  
  
  dt <- as.data.table(matrix)
  dt <- cbind(colnames(dt), dt)
  dt <- melt(dt, measure.vars = c(2,3,4,5,6,7))
  dt$value <- as.numeric(round(dt$value, 1))
  dt$variable <- factor(dt$variable, c("BC", "AC", "AB", "C", "B", "A"))
  dt$V1 <- factor(dt$V1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  mat_plt_con_mix_rmse <- ggplot(dt)+aes(x=V1, y=variable)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5,linetype = 1)+
    scale_x_discrete(position = "top") +
    scale_colour_manual("white")+
    labs(x="Mixed", y="Coniferous")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value), color = "black", size=3)+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}

{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Mixed" & variable=="MPE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  matrix <- matrix(ncol=6, nrow=6)
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
        matrix[i,j] <-  t2.stk[group1==i&group2==j]$estimate
    }
  }
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Coniferous" & variable=="MPE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- oneway.test(value~grp, data = tests)
  t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = TRUE)
  
  setDT(t2)
  t2 <- t2[,c("group1","group2","estimate", "p.adj.signif")]
  t2.inv <- t2
  names(t2.inv)[1:2] <- c("group2", "group1")
  t2.inv <- t2.inv[, estimate:=estimate*(-1)]
  t2.stk <- rbind(t2, t2.inv)
  
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in 1:6){
    for(j in 1:6){
      if(i<j)
        matrix[i,j] <-  t2.stk[group1==rownames(matrix)[i]&group2==colnames(matrix)[j]]$estimate
    }
  }
  
  
  dt <- as.data.table(matrix)
  dt <- cbind(colnames(dt), dt)
  dt <- melt(dt, measure.vars = c(2,3,4,5,6,7))
  dt$value <- as.numeric(round(dt$value, 1))
  dt$variable <- factor(dt$variable, c("BC", "AC", "AB", "C", "B", "A"))
  dt$V1 <- factor(dt$V1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  mat_plt_con_mix_mpe <- ggplot(dt)+aes(x=V1, y=variable)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5,linetype = 1)+
    scale_x_discrete(position = "top") +
    scale_colour_manual("white")+
    labs(x="Mixed", y="Coniferous")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value), color = "black", size=3)+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}



mat_plts <- ggarrange(mat_plt_rip_bro_r2, mat_plt_rip_bro_rmse, mat_plt_rip_bro_mpe,
                      mat_plt_con_mix_r2, mat_plt_con_mix_rmse, mat_plt_con_mix_mpe)


ggsave(mat_plts, file="D:/1_Work/Dropbox/2_Publications/2_paper/revisions/mat_plots.svg", 
       width=15.9*1.25, height=10.6*1.25, units="cm", dpi=640)

##########################################################################################################################################################

{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Riparian" & variable=="R2" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- pairwise.var.test(tests$value, tests$grp, p.method="bonferroni")
  x <- t1$p.value
  x <- as.data.frame(x)
  x <- cbind(rownames(x), x)
  x <- na.omit(reshape2::melt(x))
  colnames(x) <- c("g1", "g2", "value")
  x1 <- x[,c("g2", "g1", "value")]
  colnames(x1) <- c("g1", "g2", "value")
  x <- rbind(x, x1)
  x <- setDT(x)
  matrix <- matrix(ncol=6, nrow=6)
  vardf <- data.frame()
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
      {
        matrix[i,j] <-  x[g1==i&g2==j]$value
        mdlmets1 <- mdlmets[exp == i]
        mdlmets2 <- mdlmets[exp == j]
        vars <- ifelse(var(mdlmets1$value)>var(mdlmets2$value), 
                       var(mdlmets1$value)/var(mdlmets2$value), 
                       var(mdlmets2$value)/var(mdlmets1$value))
        vardf <- as.data.frame(rbind(vardf, c(i,j,vars), c(j,i,vars)))
      }
    }
  }
  vardf <- unique(vardf)
  colnames(vardf) <- c("g1", "g2", "value")
  setDT(vardf)
  matrix[lower.tri(matrix)] <- NA
  matrix <- melt(matrix) 
  setDT(matrix)
  setkeyv(matrix, c("Var1", "Var2"))
  setkeyv(vardf, c("g1", "g2"))
  
  matrix <- vardf[matrix]
  matrix <- matrix[, id:=paste0(g1, "-", g2)]
  matrix1 <- matrix[id %in% c("A-B", "A-C", "A-AB","A-AC","A-BC",
                            "B-C","B-AB","B-AC","B-BC","C-AB",
                            "C-AC","C-BC","AB-AC","AB-BC","AC-BC",
                            "A-A", "B-B", "C-C", "AB-AB", "AC-AC", "BC-BC")]
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 2))
  matrix1$g2 <- factor(matrix1$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix1$g1 <- factor(matrix1$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Broadleaf" & variable=="R2" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- pairwise.var.test(tests$value, tests$grp, p.method="bonferroni")
  x <- t1$p.value
  x <- as.data.frame(x)
  x <- cbind(rownames(x), x)
  x <- na.omit(reshape2::melt(x))
  colnames(x) <- c("g1", "g2", "value")
  x1 <- x[,c("g2", "g1", "value")]
  colnames(x1) <- c("g1", "g2", "value")
  x <- rbind(x, x1)
  x <- setDT(x)
  matrix <- matrix(ncol=6, nrow=6)
  vardf <- data.frame()
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
      {
        matrix[i,j] <-  x[g1==i&g2==j]$value
        mdlmets1 <- mdlmets[exp == i]
        mdlmets2 <- mdlmets[exp == j]
        vars <- ifelse(var(mdlmets1$value)>var(mdlmets2$value), 
                       var(mdlmets1$value)/var(mdlmets2$value), 
                       var(mdlmets2$value)/var(mdlmets1$value))
        vardf <- as.data.frame(rbind(vardf, c(i,j,vars), c(j,i,vars)))
      }
    }
  }
  vardf <- unique(vardf)
  colnames(vardf) <- c("g2", "g1", "value")
  setDT(vardf)
  matrix[lower.tri(matrix)] <- NA
  matrix <- melt(matrix) 
  setDT(matrix)
  setkeyv(matrix, c("Var1", "Var2"))
  setkeyv(vardf, c("g1", "g2"))
  
  matrix <- vardf[matrix]
  matrix <- matrix[, id:=paste0(g2, "-", g1)]
  matrix2 <- matrix[id %in% c("B-A", "C-A", "AB-A","AC-A","BC-A",
                             "C-B","AB-B","AC-B","BC-B","AB-C",
                             "AC-C","BC-C","AC-AB","BC-AB","BC-AC",
                             "A-A", "B-B", "C-C", "AB-AB", "AC-AC", "BC-BC")]
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 2))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig))+
    scale_colour_manual(values=c("black", "red"))+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}
{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Riparian" & variable=="rRMSE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- pairwise.var.test(tests$value, tests$grp, p.method="bonferroni")
  x <- t1$p.value
  x <- as.data.frame(x)
  x <- cbind(rownames(x), x)
  x <- na.omit(reshape2::melt(x))
  colnames(x) <- c("g1", "g2", "value")
  x1 <- x[,c("g2", "g1", "value")]
  colnames(x1) <- c("g1", "g2", "value")
  x <- rbind(x, x1)
  x <- setDT(x)
  matrix <- matrix(ncol=6, nrow=6)
  vardf <- data.frame()
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
      {
        matrix[i,j] <-  x[g1==i&g2==j]$value
        mdlmets1 <- mdlmets[exp == i]
        mdlmets2 <- mdlmets[exp == j]
        vars <- ifelse(var(mdlmets1$value)>var(mdlmets2$value), 
                       var(mdlmets1$value)/var(mdlmets2$value), 
                       var(mdlmets2$value)/var(mdlmets1$value))
        vardf <- as.data.frame(rbind(vardf, c(i,j,vars), c(j,i,vars)))
      }
    }
  }
  vardf <- unique(vardf)
  colnames(vardf) <- c("g1", "g2", "value")
  setDT(vardf)
  matrix[lower.tri(matrix)] <- NA
  matrix <- melt(matrix) 
  setDT(matrix)
  setkeyv(matrix, c("Var1", "Var2"))
  setkeyv(vardf, c("g1", "g2"))
  
  matrix <- vardf[matrix]
  matrix <- matrix[, id:=paste0(g1, "-", g2)]
  matrix1 <- matrix[id %in% c("A-B", "A-C", "A-AB","A-AC","A-BC",
                              "B-C","B-AB","B-AC","B-BC","C-AB",
                              "C-AC","C-BC","AB-AC","AB-BC","AC-BC",
                              "A-A", "B-B", "C-C", "AB-AB", "AC-AC", "BC-BC")]
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 2))
  matrix1$g2 <- factor(matrix1$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix1$g1 <- factor(matrix1$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Broadleaf" & variable=="rRMSE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- pairwise.var.test(tests$value, tests$grp, p.method="bonferroni")
  x <- t1$p.value
  x <- as.data.frame(x)
  x <- cbind(rownames(x), x)
  x <- na.omit(reshape2::melt(x))
  colnames(x) <- c("g1", "g2", "value")
  x1 <- x[,c("g2", "g1", "value")]
  colnames(x1) <- c("g1", "g2", "value")
  x <- rbind(x, x1)
  x <- setDT(x)
  matrix <- matrix(ncol=6, nrow=6)
  vardf <- data.frame()
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
      {
        matrix[i,j] <-  x[g1==i&g2==j]$value
        mdlmets1 <- mdlmets[exp == i]
        mdlmets2 <- mdlmets[exp == j]
        vars <- ifelse(var(mdlmets1$value)>var(mdlmets2$value), 
                       var(mdlmets1$value)/var(mdlmets2$value), 
                       var(mdlmets2$value)/var(mdlmets1$value))
        vardf <- as.data.frame(rbind(vardf, c(i,j,vars), c(j,i,vars)))
      }
    }
  }
  vardf <- unique(vardf)
  colnames(vardf) <- c("g2", "g1", "value")
  setDT(vardf)
  matrix[lower.tri(matrix)] <- NA
  matrix <- melt(matrix) 
  setDT(matrix)
  setkeyv(matrix, c("Var1", "Var2"))
  setkeyv(vardf, c("g1", "g2"))
  
  matrix <- vardf[matrix]
  matrix <- matrix[, id:=paste0(g2, "-", g1)]
  matrix2 <- matrix[id %in% c("B-A", "C-A", "AB-A","AC-A","BC-A",
                              "C-B","AB-B","AC-B","BC-B","AB-C",
                              "AC-C","BC-C","AC-AB","BC-AB","BC-AC",
                              "A-A", "B-B", "C-C", "AB-AB", "AC-AC", "BC-BC")]
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 2))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig))+
    scale_colour_manual(values=c("black", "red"))+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}
{
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Riparian" & variable=="MPE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- pairwise.var.test(tests$value, tests$grp, p.method="bonferroni")
  x <- t1$p.value
  x <- as.data.frame(x)
  x <- cbind(rownames(x), x)
  x <- na.omit(reshape2::melt(x))
  colnames(x) <- c("g1", "g2", "value")
  x1 <- x[,c("g2", "g1", "value")]
  colnames(x1) <- c("g1", "g2", "value")
  x <- rbind(x, x1)
  x <- setDT(x)
  matrix <- matrix(ncol=6, nrow=6)
  vardf <- data.frame()
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
      {
        matrix[i,j] <-  x[g1==i&g2==j]$value
        mdlmets1 <- mdlmets[exp == i]
        mdlmets2 <- mdlmets[exp == j]
        vars <- ifelse(var(mdlmets1$value)>var(mdlmets2$value), 
                       var(mdlmets1$value)/var(mdlmets2$value), 
                       var(mdlmets2$value)/var(mdlmets1$value))
        vardf <- as.data.frame(rbind(vardf, c(i,j,vars), c(j,i,vars)))
      }
    }
  }
  vardf <- unique(vardf)
  colnames(vardf) <- c("g1", "g2", "value")
  setDT(vardf)
  matrix[lower.tri(matrix)] <- NA
  matrix <- melt(matrix) 
  setDT(matrix)
  setkeyv(matrix, c("Var1", "Var2"))
  setkeyv(vardf, c("g1", "g2"))
  
  matrix <- vardf[matrix]
  matrix <- matrix[, id:=paste0(g1, "-", g2)]
  matrix1 <- matrix[id %in% c("A-B", "A-C", "A-AB","A-AC","A-BC",
                              "B-C","B-AB","B-AC","B-BC","C-AB",
                              "C-AC","C-BC","AB-AC","AB-BC","AC-BC",
                              "A-A", "B-B", "C-C", "AB-AB", "AC-AC", "BC-BC")]
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 2))
  matrix1$g2 <- factor(matrix1$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix1$g1 <- factor(matrix1$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  
  
  mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
  mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
  
  
  mdlmets$forest_type <- factor(mdlmets$forest_type, levels=c("Riparian", "Coniferous", 
                                                              "Broadleaf", "Mixed"))
  mdlmets$exp <- factor(mdlmets$exp, levels=c("fl1","A", "B", "C",
                                              "fl2", "AB", "AC", "BC", "fl3")) 
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Broadleaf" & variable=="MPE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
  grp <- as.character(mdlmets$exp)
  value <- mdlmets$value
  tests <- data.frame(grp, value)
  tests$grp <- as.factor(tests$grp)
  t1 <- pairwise.var.test(tests$value, tests$grp, p.method="bonferroni")
  x <- t1$p.value
  x <- as.data.frame(x)
  x <- cbind(rownames(x), x)
  x <- na.omit(reshape2::melt(x))
  colnames(x) <- c("g1", "g2", "value")
  x1 <- x[,c("g2", "g1", "value")]
  colnames(x1) <- c("g1", "g2", "value")
  x <- rbind(x, x1)
  x <- setDT(x)
  matrix <- matrix(ncol=6, nrow=6)
  vardf <- data.frame()
  rownames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  colnames(matrix) <- c("A", "B", "C", "AB", "AC", "BC")
  for(i in rownames(matrix)){
    for(j in colnames(matrix)){
      if(i==j)
        matrix[i,j] <- NA
      else
      {
        matrix[i,j] <-  x[g1==i&g2==j]$value
        mdlmets1 <- mdlmets[exp == i]
        mdlmets2 <- mdlmets[exp == j]
        vars <- ifelse(var(mdlmets1$value)>var(mdlmets2$value), 
                       var(mdlmets1$value)/var(mdlmets2$value), 
                       var(mdlmets2$value)/var(mdlmets1$value))
        vardf <- as.data.frame(rbind(vardf, c(i,j,vars), c(j,i,vars)))
      }
    }
  }
  vardf <- unique(vardf)
  colnames(vardf) <- c("g2", "g1", "value")
  setDT(vardf)
  matrix[lower.tri(matrix)] <- NA
  matrix <- melt(matrix) 
  setDT(matrix)
  setkeyv(matrix, c("Var1", "Var2"))
  setkeyv(vardf, c("g1", "g2"))
  
  matrix <- vardf[matrix]
  matrix <- matrix[, id:=paste0(g2, "-", g1)]
  matrix2 <- matrix[id %in% c("B-A", "C-A", "AB-A","AC-A","BC-A",
                              "C-B","AB-B","AC-B","BC-B","AB-C",
                              "AC-C","BC-C","AC-AB","BC-AB","BC-AC",
                              "A-A", "B-B", "C-C", "AB-AB", "AC-AC", "BC-BC")]
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 2))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig))+
    scale_colour_manual(values=c("black", "red"))+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}

































mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Riparian" & variable=="R2" & exp %in% c("A", "B") & Metrics=="ref"]
ggplot(mdlmets)+aes(x=exp, y=value)+geom_boxplot()



mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
mdlmets <- mdlmets[exp %in% c("A", "B", "C", "AB", "AC", "BC")]
###########
#paiwise tests between ref and vox
###########
forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
exps <- c("fl1", "A", "B", "C", "fl2", "AB", "AC", "BC", "fl3")
ref.vox.test <- NULL
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      for(ex in exps){
        tests <- mdlmets[forest_type== forest&
                           Forest_attr==forattr & 
                           fl %in% c("fl1", "fl2", "fl3")&
                           variable==msr  &
                           exp==ex, c("value", "Metrics")]
        grp <- as.character(tests$Metrics)
        value <- tests$value
        tests <- data.frame(grp, value)
        tests$grp <- as.factor(tests$grp)
        t1 <- oneway.test(value~grp, data = tests)
        t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = FALSE)
        t3 <- var.test(value~grp, data = tests, alternative = "greater", conf.level = 0.99)
        varp <- t3$p.value
        varF <- t3$statistic
        testF <- var(tests$value[tests$grp=="ref"])/var(tests$value[tests$grp=="vox"])
        t3$conf.int
        yy <- cbind(forest, forattr, msr, ex, t1$p.value, t2, varp, varF, testF)
        ref.vox.test <- rbind(ref.vox.test, yy)
      }}}}



xx <- mdlmets[forest_type== "Mixed"&
                Forest_attr=="Basal area" & 
            fl %in% c("fl1", "fl2", "fl3")&
            variable=="R2"  &
            exp %in% c("A", "B", "C", "AB", "AC", "BC", "fl3"), c("value", "Metrics", "exp")]

###########
#paiwise tests between ref and vox (A+B+C+AB+AC+BC+fl3)
###########
forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
ref.vox.test1 <- NULL
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
        tests <- mdlmets[forest_type== forest&
                           Forest_attr==forattr & 
                           fl %in% c("fl1", "fl2", "fl3")&
                           variable==msr  &
                           exp %in% c("A", "B", "C", "AB", "AC", "BC", "fl3"), c("value", "Metrics")]
        grp <- as.character(tests$Metrics)
        value <- tests$value
        tests <- data.frame(grp, value)
        tests$grp <- as.factor(tests$grp)
        t1 <- oneway.test(value~grp, data = tests)
        t2 <- games_howell_test(tests, value~grp , conf.level = 0.95, detailed = FALSE)
        t3 <- var.test(value~grp, data = tests, alternative = "greater", conf.level = 0.99)
        varp <- t3$p.value
        varF <- t3$statistic
        testF <- var(tests$value[tests$grp=="ref"])/var(tests$value[tests$grp=="vox"])
        t3$conf.int
        m1 <- mean(tests[grp=="ref",]$value)
        m2 <- mean(tests[grp=="vox",]$value)
        v1 <- sd(tests[grp=="ref",]$value)
        v2 <- sd(tests[grp=="vox",]$value)
        pcm <- ((m2-m1)/m1)*100
        pcr <- ((v2-v1)/v1)*100
        yy <- cbind(forest, forattr, msr, t1$p.value, t2, varp, varF, testF, m1, m2, v1, v2, pcm, pcr)
        ref.vox.test1 <- rbind(ref.vox.test1, yy)
      }}}
ref.vox.test1$msr <- factor(ref.vox.test1$msr, levels = c("R2", "rRMSE", "MPE"))


ref.vox.test1$forest <- factor(ref.vox.test1$forest, levels = c("Riparian",
                                                                "Coniferous",
                                                                "Broadleaved",
                                                                "Mixed"))
setDT(ref.vox.test1)

ref.vox.test1 <- ref.vox.test1[,sig:=ifelse(p.adj.signif=="ns", "no", "yes")]
ref.vox.test1 <- ref.vox.test1[,sigvarp:=ifelse(varp<=0.01, "yes", "no")]
ref.vox.test1 <- setDT(ref.vox.test1)

rfvt <- melt(ref.vox.test1, measure.vars = c("estimate", "varF"))

g1 <- ggplot(data=rfvt[forattr=="Stem volume"], aes(x=forattr, y=value, 
                                                     colour=forest, shape=sigvarp))+
  geom_point(size=2)+
  facet_wrap(variable~msr, scales="free" )+
  labs(x = "Forest attribute",
       y = "")+
  scale_colour_manual(values=c("black", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_x_discrete(labels= c("Vst"))+
  theme_few()+
  theme(axis.title.y = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text=element_text(family="serif", size=8*(96/72)))



mdlmets1 <- mdlmets[exp %in% c("A", "B", "C", "AB", "AC", "BC")&Metrics=="ref"]
mdlmets.vst <- mdlmets1[,.(mn=mean(value), sdev=sd(value)), by=c("Forest_attr", "variable", "forest_type", "fl")]
write.csv(mdlmets.vst, "D:/1_Work/Dropbox/2_Publications/2_paper/revisions/table_s3_2.csv")

xx <- dcast(mdlmets.vst, .~fl, value.var = c("mn", "sdev"))





####tests
mdlmets1 <- mdlmets[exp %in% c("A","B","C", "AB", "AC", "BC")&Metrics=="ref"]
xy <- mdlmets1[Forest_attr=="Stem volume"&forest_type=="Coniferous"&variable=="R2"]
ggplot(xy)+aes(x=fl, y=value)+geom_boxplot()
xy <- xy[,.(mn=mean(value)), by=fl]





mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
mdlmets1 <- mdlmets[exp %in% c("A", "B", "C", "AB", "AC", "BC")]
###########
#paiwise tests fl1 and fl2
###########
tests <- mdlmets1[forest_type== forest&
                   Forest_attr==forattr & 
                   fl %in% c("fl1", "fl2")&
                   variable==msr  &
                   Metrics=="ref" , c("value", "fl")]
grp <- as.character(tests$fl)
value <- tests$value
tests <- data.frame(grp, value)
tests$grp <- as.factor(tests$grp)
t1 <- oneway.test(value~grp, data = tests)
t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = FALSE)
t3 <- var.test(value~grp, data = tests, alternative = "greater", conf.level = 0.99)
varp <- t3$p.value
varF <- t3$statistic
testF <- var(tests$value[tests$grp=="fl1"])/var(tests$value[tests$grp=="fl2"])
t3$conf.int




forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
fl1.fl2.test <- NULL
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests <- mdlmets[forest_type== forest&
                          Forest_attr==forattr & 
                          fl %in% c("fl1", "fl2")&
                          variable==msr  &
                          Metrics=="ref" , c("value", "fl")]
      grp <- as.character(tests$fl)
      value <- tests$value
      tests <- data.frame(grp, value)
      tests$grp <- as.factor(tests$grp)
      t1 <- oneway.test(value~grp, data = tests)
      t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = FALSE)
      t3 <- var.test(value~grp, data = tests, alternative = "greater", conf.level = 0.99)
      varp <- t3$p.value
      varF <- t3$statistic
      testF <- var(tests$value[tests$grp=="fl1"])/var(tests$value[tests$grp=="fl2"])
      t3$conf.int
      yy <- cbind(forest, forattr, msr, t1$p.value, t2, varp, varF, testF)
      ref.vox.test <- rbind(ref.vox.test, yy)
    }}}
ref.vox.test1 <- ref.vox.test





ref.vox.test1 <- NULL
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
        tests <- mdlmets[forest_type== forest&
                           Forest_attr==forattr & 
                           fl %in% c("fl1", "fl2", "fl3")&
                           variable==msr , c("value", "Metrics")]
        grp <- as.character(tests$Metrics)
        value <- tests$value
        tests <- data.frame(grp, value)
        tests$grp <- as.factor(tests$grp)
        t1 <- oneway.test(value~grp, data = tests)
        t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = FALSE)
        t3 <- var.test(value~grp, data = tests, alternative = "greater", conf.level = 0.9)
        varp <- t3$p.value
        varF <- t3$statistic
        testF <- var(tests$value[tests$grp=="ref"])/var(tests$value[tests$grp=="vox"])
        t3$conf.int
        yy <- cbind(forest, forattr, msr, ex, t1$p.value, t2, varp, varF, testF)
        ref.vox.test1 <- rbind(ref.vox.test1, yy)
      }}}



for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests <- mdlmets[Forest_attr==forattr &
                         forest_type==forest &
                         variable==msr &
                         Metrics=="ref" &
                         fl %in% c("fl1", "fl2", "fl3"),
                       c("value", "exp")]
        grp <- as.character(tests$exp)
        value <- tests$value
        tests <- data.frame(grp, value)
        tests$grp <- as.factor(tests$grp)
        t1 <- oneway.test(value~grp, data = tests)
        t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = FALSE)
        t3 <- var.test(value~grp, data = tests, alternative = "greater", conf.level = 0.99)
        varp <- t3$p.value
        varF <- t3$statistic
        testF <- var(tests$value[tests$grp=="ref"])/var(tests$value[tests$grp=="vox"])
        t3$conf.int
        yy <- cbind(forest, forattr, msr, ex, t1$p.value, t2, varp, varF, testF)
        ref.vox.test <- rbind(ref.vox.test, yy)
      }}}



ref.vox.test1$ex <- factor(ref.vox.test$ex, levels = c("fl1", "A", "B", "C",
                                                      "fl2", "AB", "AC", "BC", "fl3"))

ref.vox.test1$msr <- factor(ref.vox.test1$msr, levels = c("R2", "rRMSE", "MPE"))


ref.vox.test1$forest <- factor(ref.vox.test1$forest, levels = c("Riparian",
                                                              "Coniferous",
                                                              "Broadleaved",
                                                              "Mixed"))
setDT(ref.vox.test1)

ref.vox.test1 <- ref.vox.test1[,sig:=ifelse(p.adj.signif=="ns", "no", "yes")]
ref.vox.test1 <- ref.vox.test1[,sigvarp:=ifelse(varp<=0.01, "yes", "no")]
ref.vox.test1 <- setDT(ref.vox.test1)
ref.vox.test1 <- ref.vox.test1[,forattr:=ifelse(forattr=="Stem volume", "Vst",
                                                ifelse(forattr=="Total volume", "Vtot","BA"))]



tests <- mdlmets[forest_type=="Riparian"&
                   Forest_attr=="Total volume"& 
                   fl %in% c("fl1", "fl2", "fl3")&
                   variable=="R2"  &
                   exp=="A", c("value", "Metrics")]





rfvt <- melt(ref.vox.test1, measure.vars = c("estimate", "varF"))

g1 <- ggplot(data=rfvt[forattr=="Total volume"], aes(x=forattr, y=value, 
                                                              colour=forest, shape=sigvarp))+
  geom_point(size=2)+
  facet_wrap(variable~msr, scales="free" )+
  labs(x = "Forest attribute",
       y = "")+
  scale_shape_manual(values=c(4))+
  scale_colour_manual(values=c("black", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_x_discrete(labels= c("BA", "Vst"))+
  theme_few()+
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text=element_text(family="serif", size=8*(96/72)))

ggsave(g1, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/ref.meanFtests.svg", 
       width=7.5*1.25, height=*1.25, units="cm", dpi=640) 

g2 <- ggplot(data=ref.vox.test1[forattr!="Total volume"], aes(x=forattr, y=varF, 
                                     colour=forest, shape=sigvarp))+
  geom_point(size=2)+
  facet_wrap(msr~., scales="free" )+
  labs(x = "Experiments",
       y = "")+
  scale_shape_manual(values=c(4))+
  scale_colour_manual(values=c("black", "#0072B2", "#D55E00", "#CC79A7"))+
  scale_x_discrete(labels= c("BA", "Vst"))+
  theme_few()+
  theme(axis.title.y = element_blank(),
        legend.position = "top",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text=element_text(family="serif", size=8*(96/72)))

ggsave(g2, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/ref.vox.ftests.svg", 
       width=15*1.25, height=12*1.25, units="cm", dpi=640) 










tst <- mdlmets[forest_type== forest&
          Forest_attr==forattr & 
          exp %in% c("fl1", "fl2", "fl3")&
          variable=="R2"
          , c("value", "Metrics", "exp")]
















g2 <- ggplot(data=ref.vox.test, aes(x=ex, y=estimate))+
  geom_point()+
  facet_grid(c("msr","forest"), scales="free" )

gt1 = ggplot_gtable(ggplot_build(g1))
gt2 = ggplot_gtable(ggplot_build(g2))
gt1$grobs[grep('strip-t.+1$', gt1$layout$name)] = gt2$grobs[grep('strip-t', gt2$layout$name)]
grid.draw(gt1)

gt.side1 = gtable_filter(gt2, 'strip-r-1')
gt.side2 = gtable_filter(gt2, 'strip-r-2')
gt.side3 = gtable_filter(gt2, 'strip-r-3')

gt1 = gtable_add_cols(gt1, widths=gt.side1$widths[1], pos = -1)
gt1 = gtable_add_grob(gt1, zeroGrob(), t = 1, l = ncol(gt1), b=nrow(gt1))

panel_id <- gt1$layout[grep('panel-.+1$', gt1$layout$name),]
gt1 = gtable_add_grob(gt1, gt.side1, t = panel_id$t[1], l = ncol(gt1))
gt1 = gtable_add_grob(gt1, gt.side2, t = panel_id$t[2], l = ncol(gt1))
gt1 = gtable_add_grob(gt1, gt.side3, t = panel_id$t[3], l = ncol(gt1))

grid.newpage()
grid.draw(gt1)











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









tst2.f.ref <- NULL
n=1
forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests <- mdlmets[forest_type== forest&
                             Forest_attr==forattr & 
                             fl %in% c("fl1", "fl2", "fl3")&
                             variable==msr  &
                             Metrics=="vox", c("value", "exp")]
      grp <- as.character(tests$exp)
      value <- tests$value
      tests <- data.frame(grp, value)
      tests$grp <- as.factor(tests$grp)
      ftst <- pairwise.var.test(tests$value, tests$grp, alternative = "greater")
      ftst <- as.data.table(ftst[["p.value"]])
      ftst <- ftst[c(7,8),c(7,8)]
      
      tst2.f.ref[[n]] <- c(forest,
                           forattr,
                           msr,
                           ftst
                           )
      n=n+1
    }
  }
}
      


yy[c(1,2), c(1,2)]


















tst2.vox <- NULL
forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests.vox <- mdlmets[forest_type==forest&
                             Forest_attr==forattr& 
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
  



#test for differences between the means 
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
tst2.ref <- tst2.ref[x %in% c("A-B", "A-C", "B-C",
                              "AB-AC", "AB-BC", "AC-BC")]

tst2.ref$x <- factor(tst2.ref$x, levels = c("A-B", "A-C", "B-C",
                                            "AB-AC", "AB-BC", "AC-BC"))
tst2.ref$fr_type <- factor(tst2.ref$fr_type, levels = c("Riparian", "Coniferous",
                                                        "Broadleaved", "Mixed"))
tst2.ref$var <- factor(tst2.ref$var, levels = c("R2", "rRMSE","MPE"))

tst2.ref <- tst2.ref[, fl:=ifelse(x %in% c("A-B", "A-C", "B-C"), "fl1", "fl2")]

ref.tst <- ggplot(data=tst2.ref, aes(x=x, y=estimate, colour=fr_at))+
  geom_point(aes(shape=clr))+
  annotate(xmin=3.5, xmax=Inf, ymin=-Inf, ymax=Inf, geom = "rect", alpha = 0.2 )+
  facet_wrap(fr_type~var, scales = "free", ncol=3)+
  theme_few()+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "black"))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text=element_text(family="serif", size=8*(96/72)))

ggsave(ref.tst, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/ref.tst.svg", 
       width=15*1.25, height=23*1.25, units="cm", dpi=640) 





tst2.vox <- setDT(tst2.vox)
tst2.vox.ns <- tst2.vox[p.adj.signif=="ns"]
tst2.vox <- tst2.vox[,clr:= ifelse(p.adj.signif=="ns", "yes", "no")]
tst2.vox <- tst2.vox[x %in% c("A-B", "A-C", "B-C",
                              "AB-AC", "AB-BC", "AC-BC")]
tst2.vox$x <- factor(tst2.vox$x, levels = c("A-B", "A-C", "B-C",
                                            "AB-AC", "AB-BC", "AC-BC"))
tst2.vox$fr_type <- factor(tst2.vox$fr_type, levels = c("Riparian", "Coniferous",
                                                        "Broadleaved", "Mixed"))
tst2.vox$var <- factor(tst2.vox$var, levels = c("R2", "rRMSE","MPE"))

tst2.vox <- tst2.vox[, fl:=ifelse(x %in% c("A-B", "A-C", "B-C"), "fl1", "fl2")]

vox.tst <- ggplot(data=tst2.vox, aes(x=x, y=estimate, colour=fr_at))+
  geom_point(aes(shape=clr))+
  annotate(xmin=3.5, xmax=Inf, ymin=-Inf, ymax=Inf, geom = "rect", alpha = 0.2 )+
  facet_wrap(fr_type~var, scales = "free", ncol=3)+
  theme_few()+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "black"))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text=element_text(family="serif", size=8*(96/72)))




####tests for fl1, fl2, fl3

tst2.ref.fl123 <- NULL
tst2.vox.fl123 <- NULL
forest.type <- c("Mixed", "Coniferous", "Broadleaved", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
for(forest in forest.type){
  for(forattr in forest.attr){
    for(msr in mdl.measure){
      tests.ref <- mdlmets[forest_type== forest&
                             Forest_attr==forattr &
                             variable==msr  &
                             exp %in% c("fl1", "fl2", "fl3")&
                             Metrics=="ref", c("value", "exp")]
      # tests.ref <- tests.ref[exp!="C"]
      grp <- as.character(tests.ref$exp)
      value <- tests.ref$value
      tests.ref <- data.frame(grp, value)
      tests.ref$grp <- as.factor(tests.ref$grp)
      oneway.test(value~grp, data = tests.ref)
      xy <- games_howell_test(tests.ref, value~grp , conf.level = 0.99, detailed = FALSE)
      vtest <- pairwise.var.test(tests$value, tests$grp, alternative = "greater")
      vtest <- as.data.table(vtest[["p.value"]])
      xy$x <- paste0(xy$group1,'-',xy$group2)
      setDT(xy)
      lvls <- xy$x
      xy$x <- factor(xy$x, levels=lvls)
      xy$fr_type <- forest
      xy$fr_at <- forattr
      xy$var <- msr
      xy <- cbind(xy, fsig=rbind(vtest$fl1[7], vtest$fl1[8], vtest$fl2[8]))
      tst2.ref.fl123 <- rbind(tst2.ref.fl123, xy)
      tst2.ref.fl123.smry <- tst2.ref.fl123[, .(min=min(abs(estimate)),
                                    mn=mean(abs(estimate)),
                                    max=max(abs(estimate))), by=.(fr_type, fr_at, var)]
      tst2.ref.fl123.smry$Metrics <- "ref"
      
      
      tests.vox <- mdlmets[forest_type== forest&
                             Forest_attr==forattr &
                             variable==msr  &
                             exp %in% c("fl1", "fl2", "fl3")&
                             Metrics=="ref", c("value", "exp")]
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
      tst2.vox.fl123 <- rbind(tst2.vox.fl123, xy)
      tst2.vox.fl123.smry <- tst2.vox.fl123[, .(min=min(abs(estimate)),
                                    mn=mean(abs(estimate)),
                                    max=max(abs(estimate))), by=.(fr_type, fr_at, var)]
      tst2.vox.fl123.smry$Metrics <- "vox"
      
    }
  }
}




tst2.ref.fl123 <- setDT(tst2.ref.fl123)
tst2.ref.fl123 <- tst2.ref.fl123[, clr:= ifelse(p.adj.signif=="ns", "yes", "no")]
tst2.ref.fl123$x <- factor(tst2.ref.fl123$x, levels = c("fl1-fl2", "fl1-fl3", "fl2-fl3"))
tst2.ref.fl123$fr_type <- factor(tst2.ref.fl123$fr_type, levels = c("Riparian", "Coniferous",
                                                        "Broadleaved", "Mixed"))
tst2.ref.fl123$var <- factor(tst2.ref.fl123$var, levels = c("R2", "rRMSE","MPE"))

vox.tst <- ggplot(data=tst2.ref.fl123, aes(x=x, y=abs(estimate), colour=fr_at, group=fr_at))+
  geom_point(aes(shape=clr), size=2)+
  facet_wrap(fr_type~var, scales = "free", ncol=3)+
  theme_few()+
  scale_colour_manual(values=c("#0072B2", "#D55E00", "black"))+
  theme(axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 0),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        text=element_text(family="serif", size=8*(96/72)))












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




