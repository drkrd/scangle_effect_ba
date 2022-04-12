library(ggplot2)
library(ggpubr)
library(grid)
library(gtable) 
library(ggsci)
library(data.table)
library(rstatix)
library(ggthemes)
library(RVAideMemoire)




#############################################################################################################################


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
################################
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
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 1))
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
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 1))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  mat_plt_rip_bro_r2 <- ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig), size=3)+
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
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 1))
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
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 1))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  mat_plt_rip_bro_rmse <- ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig), size=3)+
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
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 1))
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
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 1))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  mat_plt_rip_bro_mpe <- ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Broadleaf", y="Riparian")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig), size=3)+
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
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Coniferous" & variable=="R2" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
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
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 1))
  matrix1$g2 <- factor(matrix1$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix1$g1 <- factor(matrix1$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  
  
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
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 1))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  mat_plt_con_mix_r2 <- ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Mixed", y="Coniferous")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig), size=3)+
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
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Coniferous" & variable=="rRMSE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
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
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 1))
  matrix1$g2 <- factor(matrix1$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix1$g1 <- factor(matrix1$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  
  
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
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 1))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  mat_plt_con_mix_rmse <- ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Mixed", y="Coniferous")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig), size=3)+
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
  
  mdlmets <- mdlmets[Forest_attr=="Stem volume" & forest_type=="Coniferous" & variable=="MPE" & exp %in% c("A", "B", "C", "AB", "AC", "BC") & Metrics=="ref"]
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
  
  
  
  matrix1$value <- as.numeric(round(as.numeric(matrix1$value), 1))
  matrix1$g2 <- factor(matrix1$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix1$g1 <- factor(matrix1$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  
  
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
  
  matrix2$value <- as.numeric(round(as.numeric(matrix2$value), 1))
  matrix2$g2 <- factor(matrix2$g2, c("BC", "AC", "AB", "C", "B", "A"))
  matrix2$g1 <- factor(matrix2$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  colnames(matrix2) <- c("g1", "g2", "value", "i.value", "id")
  dt <- rbind(matrix1, matrix2)
  dt$g2 <- factor(dt$g2, c("BC", "AC", "AB", "C", "B", "A"))
  dt$g1 <- factor(dt$g1, c("A", "B", "C", "AB", "AC", "BC"))
  
  
  
  
  dt <- dt[, sig:=ifelse(i.value<0.05, "black",
                         ifelse(is.na(i.value), "white", "red"))]
  
  mat_plt_con_mix_mpe <- ggplot(dt)+aes(x=g1, y=g2)+
    geom_tile(fill = NA, color = 'grey', lwd = 0.5, linetype = 1)+
    scale_x_discrete(position = "top") +
    labs(x="Mixed", y="Coniferous")+
    coord_fixed()+
    theme_few()+
    geom_text(aes(label = value, colour=sig), size=3)+
    scale_colour_manual(values=c("black", "red"))+
    theme(text=element_text(family="serif", size=7*(96/72)),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          legend.position = "none")
}


mat_plts <- ggarrange(mat_plt_rip_bro_r2, mat_plt_rip_bro_rmse, mat_plt_rip_bro_mpe,
                      mat_plt_con_mix_r2, mat_plt_con_mix_rmse, mat_plt_con_mix_mpe)
ggsave(mat_plts, file="D:/1_Work/Dropbox/2_Publications/2_paper/revisions/mat_plotsF.svg", 
       width=15.9*1.25, height=10.6*1.25, units="cm", dpi=640)
##########################
##########################


mdlmets <- rbind(ciron.mdlmets, bauges.mdlmets)
mdlmets <- mdlmets[, forest_type:=ifelse(forest_type=="Broadleaved", "Broadleaf", forest_type)]
forest.type <- c("Mixed", "Coniferous", "Broadleaf", "Riparian")
forest.attr <- c("Basal area", "Stem volume", "Total volume")
mdl.measure <- c("R2", "rRMSE", "MPE")
exps <- c("fl1", "fl2", "fl3")
ref.test <- NULL
for(forest in forest.type){
  for(msr in mdl.measure){
    tests <- mdlmets[forest_type==forest&
                       Forest_attr=="Stem volume" & 
                       exp %in% c("fl1", "fl2", "fl3")&
                       variable==msr  &
                       Metrics=="ref", c("value", "exp")]
    grp <- as.character(tests$exp)
    value <- tests$value
    tests <- data.frame(grp, value)
    tests$grp <- as.factor(tests$grp)
    t1 <- oneway.test(value~grp, data = tests)
    t2 <- games_howell_test(tests, value~grp , conf.level = 0.99, detailed = T)
    yy <- cbind(forest, msr, t1$p.value, t2)
    testsx <- tests
    setDT(testsx)
    tests1 <- testsx[, .(mn=mean(value), sdev=sd(value)), by="grp"]
    v1 <- tests1[grp=="fl1"]$mn
    v2 <- tests1[grp=="fl2"]$mn
    v3 <- tests1[grp=="fl3"]$mn
    v4 <- tests1[grp=="fl1"]$sdev
    v5 <- tests1[grp=="fl2"]$sdev
    v6 <- tests1[grp=="fl3"]$sdev
    df <- data.frame(m1=c(v1,v1,v2), m2=c(v2,v3,v3), s1=c(v4,v4,v5), s2=c(v5,v6,v6))
    yy <- cbind(yy, df)
    ref.test <- rbind(ref.test, yy)
  }}
write.csv(ref.test, "D:/1_Work/Dropbox/2_Publications/2_paper/revisions/ref_test1.csv")


