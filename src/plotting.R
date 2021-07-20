
mdlmets_ref1 <- mdlmets_ref1[,for_attr:=ifelse(for_attr=="g", "Basal area",
                                               ifelse(for_attr=="vtig", "Stem volume", "Total volume"))]



mdlmets_ref1_con <- mdlmets_ref1[forest=="con"]
mdlmets_ref1_con <- melt(mdlmets_ref1_con, measure.vars = c("R2", "RMSE", "MAE", "MPE"))
mdlmets_ref1_con$variable <- factor(mdlmets_ref1_con$variable, levels=c("R2", "MPE", "RMSE", "MAE")) 
ggplot(data=mdlmets_ref1_con, aes(x=exp, y=value, group=mets))+
  labs(title = "Coniferes", x="Experiment", y="Value")+
  geom_line(aes(linetype=mets), size=1)+geom_point(aes(shape=mets), size=3)+
  facet_grid(variable~for_attr, scales = "free")+
  theme_base()


mdlmets_ref1_feu <- mdlmets_ref1[forest=="feu"]
mdlmets_ref1_feu <- melt(mdlmets_ref1_feu, measure.vars = c("R2", "RMSE", "MAE", "MPE"))
mdlmets_ref1_feu$variable <- factor(mdlmets_ref1_feu$variable, levels=c("R2", "MPE", "RMSE", "MAE")) 
ggplot(data=mdlmets_ref1_feu, aes(x=exp, y=value, group=mets))+
  labs(title = "Feuillus", x="Experiment", y="Value")+
  geom_line(aes(linetype=mets), size=1)+geom_point(aes(shape=mets), size=3)+
  facet_grid(variable~for_attr, scales = "free")+
  theme_base()


mdlmets_ref1_mix <- mdlmets_ref1[forest=="mix"]
mdlmets_ref1_mix <- melt(mdlmets_ref1_mix, measure.vars = c("R2", "RMSE", "MAE", "MPE"))
mdlmets_ref1_mix$variable <- factor(mdlmets_ref1_mix$variable, levels=c("R2", "MPE", "RMSE", "MAE")) 
ggplot(data=mdlmets_ref1_mix, aes(x=exp, y=value, group=mets))+
  labs(title = "Mixte", x="Experiment", y="Value")+
  geom_line(aes(linetype=mets), size=0.8)+geom_point(aes(shape=mets), size=3)+
  facet_grid(variable~for_attr, scales = "free")+
  theme_base()
###############################################################################################
mdlmets_gfl1 <- rbind(mdlmets_gfl1con, 
                   mdlmets_gfl1vxcon,
                   mdlmets_gfl1acon,
                   mdlmets_gfl1vxbcon,
                   mdlmets_gfl1bcon,
                   mdlmets_gfl1vxacon,
                   mdlmets_gfl1ccon,
                   mdlmets_gfl1vxccon)
mdlmets_gfl1 <- melt(mdlmets_gfl1, id.vars = c("type", "exp", "mets"))
mdlmets_gfl1$fls <- rep("One flight line", nrow(mdlmets_gfl1))


mdlmets_gfl2 <- rbind(mdlmets_gfl2allcon, 
                   mdlmets_gfl2vxallcon,
                   mdlmets_gfl2abcon,
                   mdlmets_gfl2vxabcon,
                   mdlmets_gfl2accon,
                   mdlmets_gfl2vxaccon,
                   mdlmets_gfl2bccon,
                   mdlmets_gfl2vxbccon)
mdlmets_gfl2 <- melt(mdlmets_gfl2, id.vars = c("type", "exp", "mets"))
mdlmets_gfl2$fls <- rep("Two flight lines", nrow(mdlmets_gfl2))

mdlmets_g <- rbind(mdlmets_gfl1, mdlmets_gfl2)
mdlmets_g <- mdlmets_g[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                     ifelse(exp=="all", "allclasses 2 FL",
                                            ifelse(exp=="ab", "mostly ab",
                                                   ifelse(exp=="bc", "mostly bc",
                                                          ifelse(exp=="ac", "mostly ac",
                                                                 ifelse(exp=="ab", "mostly ab", exp))))))]

mdlmets_g <- mdlmets_g[, variable:=ifelse(variable=="r2", "R2",
                                     ifelse(variable=="rmse", "RMSE",
                                            ifelse(variable=="mae", "MAE",
                                                   ifelse(variable=="mpe", "MPE", variable))))]
                                                         

mdlmets_g$fls <- as.factor(mdlmets_g$fls)
mdlmets_g$variable <- factor(mdlmets_g$variable, levels=c("R2", "MPE", "RMSE", "MAE"))


ggplot(data=mdlmets_g, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Coniferes: Basal area", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~fls, scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_con[for_attr=="Basal area" & forest=="con"], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_size_identity()+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))

###############################################################################################################
###############################################################################################################
mdlmets_vtotfl1 <- rbind(mdlmets_vtotfl1con, 
                      mdlmets_vtotfl1vxcon,
                      mdlmets_vtotfl1acon,
                      mdlmets_vtotfl1vxbcon,
                      mdlmets_vtotfl1bcon,
                      mdlmets_vtotfl1vxacon,
                      mdlmets_vtotfl1ccon,
                      mdlmets_vtotfl1vxccon)
mdlmets_vtotfl1 <- melt(mdlmets_vtotfl1, id.vars = c("type", "exp", "mets"))
mdlmets_vtotfl1$fls <- rep("One flight line", nrow(mdlmets_vtotfl1))


mdlmets_vtotfl2 <- rbind(mdlmets_vtotfl2allcon, 
                      mdlmets_vtotfl2vxallcon,
                      mdlmets_vtotfl2abcon,
                      mdlmets_vtotfl2vxabcon,
                      mdlmets_vtotfl2accon,
                      mdlmets_vtotfl2vxaccon,
                      mdlmets_vtotfl2bccon,
                      mdlmets_vtotfl2vxbccon)
mdlmets_vtotfl2 <- melt(mdlmets_vtotfl2, id.vars = c("type", "exp", "mets"))
mdlmets_vtotfl2$fls <- rep("Two flight lines", nrow(mdlmets_vtotfl2))

mdlmets_vtot <- rbind(mdlmets_vtotfl1, mdlmets_vtotfl2)
mdlmets_vtot <- mdlmets_vtot[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                     ifelse(exp=="all", "allclasses 2 FL",
                                            ifelse(exp=="ab", "mostly ab",
                                                   ifelse(exp=="bc", "mostly bc",
                                                          ifelse(exp=="ac", "mostly ac",
                                                                 ifelse(exp=="ab", "mostly ab", exp))))))]

mdlmets_vtot <- mdlmets_vtot[, variable:=ifelse(variable=="r2", "R2",
                                          ifelse(variable=="rmse", "RMSE",
                                                 ifelse(variable=="mae", "MAE",
                                                        ifelse(variable=="mpe", "MPE", variable))))]


mdlmets_vtot$fls <- as.factor(mdlmets_vtot$fls)
mdlmets_vtot$variable <- factor(mdlmets_vtot$variable, levels=c("R2", "MPE", "RMSE", "MAE"))


ggplot(data=mdlmets_vtot, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Coniferes: Total volume", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~fls, scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_con[for_attr=="Total volume" & forest=="con"], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))


#############################################################################################################
mdlmets_vtigfl1 <- rbind(mdlmets_vtigfl1con, 
                         mdlmets_vtigfl1vxcon,
                         mdlmets_vtigfl1acon,
                         mdlmets_vtigfl1vxbcon,
                         mdlmets_vtigfl1bcon,
                         mdlmets_vtigfl1vxacon,
                         mdlmets_vtigfl1ccon,
                         mdlmets_vtigfl1vxccon)
mdlmets_vtigfl1 <- melt(mdlmets_vtigfl1, id.vars = c("type", "exp", "mets"))
mdlmets_vtigfl1$fls <- rep("One flight line", nrow(mdlmets_vtigfl1))


mdlmets_vtigfl2 <- rbind(mdlmets_vtigfl2allcon, 
                         mdlmets_vtigfl2vxallcon,
                         mdlmets_vtigfl2abcon,
                         mdlmets_vtigfl2vxabcon,
                         mdlmets_vtigfl2accon,
                         mdlmets_vtigfl2vxaccon,
                         mdlmets_vtigfl2bccon,
                         mdlmets_vtigfl2vxbccon)
mdlmets_vtigfl2 <- melt(mdlmets_vtigfl2, id.vars = c("type", "exp", "mets"))
mdlmets_vtigfl2$fls <- rep("Two flight lines", nrow(mdlmets_vtigfl2))

mdlmets_vtig <- rbind(mdlmets_vtigfl1, mdlmets_vtigfl2)
mdlmets_vtig <- mdlmets_vtig[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                           ifelse(exp=="all", "allclasses 2 FL",
                                                  ifelse(exp=="ab", "mostly ab",
                                                         ifelse(exp=="bc", "mostly bc",
                                                                ifelse(exp=="ac", "mostly ac",
                                                                       ifelse(exp=="ab", "mostly ab", exp))))))]

mdlmets_vtig <- mdlmets_vtig[, variable:=ifelse(variable=="r2", "R2",
                                                ifelse(variable=="rmse", "RMSE",
                                                       ifelse(variable=="mae", "MAE",
                                                              ifelse(variable=="mpe", "MPE", variable))))]


mdlmets_vtig$fls <- as.factor(mdlmets_vtig$fls)
mdlmets_vtig$variable <- factor(mdlmets_vtig$variable, levels=c("R2", "MPE", "RMSE", "MAE"))


ggplot(data=mdlmets_vtig, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Coniferes: Stem volume", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~fls, scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_con[for_attr=="Stem volume" & forest=="con"], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))

#########################################################################################################
###########################################################################################################
mdlmets_gfl1 <- rbind(mdlmets_gfl1feu, 
                      mdlmets_gfl1vxfeu,
                      mdlmets_gfl1afeu,
                      mdlmets_gfl1vxbfeu,
                      mdlmets_gfl1bfeu,
                      mdlmets_gfl1vxafeu,
                      mdlmets_gfl1cfeu,
                      mdlmets_gfl1vxcfeu)
mdlmets_gfl1 <- melt(mdlmets_gfl1, id.vars = c("type", "exp", "mets"))
mdlmets_gfl1$fls <- rep("One flight line", nrow(mdlmets_gfl1))


mdlmets_gfl2 <- rbind(mdlmets_gfl2allfeu, 
                      mdlmets_gfl2vxallfeu,
                      mdlmets_gfl2abfeu,
                      mdlmets_gfl2vxabfeu,
                      mdlmets_gfl2acfeu,
                      mdlmets_gfl2vxacfeu,
                      mdlmets_gfl2bcfeu,
                      mdlmets_gfl2vxbcfeu)
mdlmets_gfl2 <- melt(mdlmets_gfl2, id.vars = c("type", "exp", "mets"))
mdlmets_gfl2$fls <- rep("Two flight lines", nrow(mdlmets_gfl2))

mdlmets_g <- rbind(mdlmets_gfl1, mdlmets_gfl2)
mdlmets_g <- mdlmets_g[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                     ifelse(exp=="all", "allclasses 2 FL",
                                            ifelse(exp=="ab", "mostly ab",
                                                   ifelse(exp=="bc", "mostly bc",
                                                          ifelse(exp=="ac", "mostly ac",
                                                                 ifelse(exp=="ab", "mostly ab", exp))))))]

mdlmets_g <- mdlmets_g[, variable:=ifelse(variable=="r2", "R2",
                                          ifelse(variable=="rmse", "RMSE",
                                                 ifelse(variable=="mae", "MAE",
                                                        ifelse(variable=="mpe", "MPE", variable))))]


mdlmets_g$fls <- as.factor(mdlmets_g$fls)
mdlmets_g$variable <- factor(mdlmets_g$variable, levels=c("R2", "MPE", "RMSE", "MAE"))


ggplot(data=mdlmets_g, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Feuillus: Basal area", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~fls, scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_feu[for_attr=="Basal area" & forest=="feu"], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))

###############################################################################################################
mdlmets_vtotfl1 <- rbind(mdlmets_vtotfl1feu, 
                         mdlmets_vtotfl1vxfeu,
                         mdlmets_vtotfl1afeu,
                         mdlmets_vtotfl1vxbfeu,
                         mdlmets_vtotfl1bfeu,
                         mdlmets_vtotfl1vxafeu,
                         mdlmets_vtotfl1cfeu,
                         mdlmets_vtotfl1vxcfeu)
mdlmets_vtotfl1 <- melt(mdlmets_vtotfl1, id.vars = c("type", "exp", "mets"))
mdlmets_vtotfl1$fls <- rep("One flight line", nrow(mdlmets_vtotfl1))


mdlmets_vtotfl2 <- rbind(mdlmets_vtotfl2allfeu, 
                         mdlmets_vtotfl2vxallfeu,
                         mdlmets_vtotfl2abfeu,
                         mdlmets_vtotfl2vxabfeu,
                         mdlmets_vtotfl2acfeu,
                         mdlmets_vtotfl2vxacfeu,
                         mdlmets_vtotfl2bcfeu,
                         mdlmets_vtotfl2vxbcfeu)
mdlmets_vtotfl2 <- melt(mdlmets_vtotfl2, id.vars = c("type", "exp", "mets"))
mdlmets_vtotfl2$fls <- rep("Two flight lines", nrow(mdlmets_vtotfl2))

mdlmets_vtot <- rbind(mdlmets_vtotfl1, mdlmets_vtotfl2)
mdlmets_vtot <- mdlmets_vtot[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                           ifelse(exp=="all", "allclasses 2 FL",
                                                  ifelse(exp=="ab", "mostly ab",
                                                         ifelse(exp=="bc", "mostly bc",
                                                                ifelse(exp=="ac", "mostly ac",
                                                                       ifelse(exp=="ab", "mostly ab", exp))))))]

mdlmets_vtot <- mdlmets_vtot[, variable:=ifelse(variable=="r2", "R2",
                                                ifelse(variable=="rmse", "RMSE",
                                                       ifelse(variable=="mae", "MAE",
                                                              ifelse(variable=="mpe", "MPE", variable))))]


mdlmets_vtot$fls <- as.factor(mdlmets_vtot$fls)
mdlmets_vtot$variable <- factor(mdlmets_vtot$variable, levels=c("R2", "MPE", "RMSE", "MAE"))


ggplot(data=mdlmets_vtot, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Feuillus: Total volume", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~fls, scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_feu[for_attr=="Total volume" & forest=="feu"], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))


#############################################################################################################
mdlmets_vtigfl1 <- rbind(mdlmets_vtigfl1feu, 
                         mdlmets_vtigfl1vxfeu,
                         mdlmets_vtigfl1afeu,
                         mdlmets_vtigfl1vxbfeu,
                         mdlmets_vtigfl1bfeu,
                         mdlmets_vtigfl1vxafeu,
                         mdlmets_vtigfl1cfeu,
                         mdlmets_vtigfl1vxcfeu)
mdlmets_vtigfl1 <- melt(mdlmets_vtigfl1, id.vars = c("type", "exp", "mets"))
mdlmets_vtigfl1$fls <- rep("One flight line", nrow(mdlmets_vtigfl1))


mdlmets_vtigfl2 <- rbind(mdlmets_vtigfl2allfeu, 
                         mdlmets_vtigfl2vxallfeu,
                         mdlmets_vtigfl2abfeu,
                         mdlmets_vtigfl2vxabfeu,
                         mdlmets_vtigfl2acfeu,
                         mdlmets_vtigfl2vxacfeu,
                         mdlmets_vtigfl2bcfeu,
                         mdlmets_vtigfl2vxbcfeu)
mdlmets_vtigfl2 <- melt(mdlmets_vtigfl2, id.vars = c("type", "exp", "mets"))
mdlmets_vtigfl2$fls <- rep("Two flight lines", nrow(mdlmets_vtigfl2))

mdlmets_vtig <- rbind(mdlmets_vtigfl1, mdlmets_vtigfl2)
mdlmets_vtig <- mdlmets_vtig[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                           ifelse(exp=="all", "allclasses 2 FL",
                                                  ifelse(exp=="ab", "mostly ab",
                                                         ifelse(exp=="bc", "mostly bc",
                                                                ifelse(exp=="ac", "mostly ac",
                                                                       ifelse(exp=="ab", "mostly ab", exp))))))]

mdlmets_vtig <- mdlmets_vtig[, variable:=ifelse(variable=="r2", "R2",
                                                ifelse(variable=="rmse", "RMSE",
                                                       ifelse(variable=="mae", "MAE",
                                                              ifelse(variable=="mpe", "MPE", variable))))]


mdlmets_vtig$fls <- as.factor(mdlmets_vtig$fls)
mdlmets_vtig$variable <- factor(mdlmets_vtig$variable, levels=c("R2", "MPE", "RMSE", "MAE"))


ggplot(data=mdlmets_vtig, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Feuillus: Stem volume", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~fls, scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_feu[for_attr=="Stem volume" & forest=="feu"], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))
################################################################################################################
################################################################################################################
mdlmets_gfl1 <- rbind(mdlmets_gfl1mix, 
                      mdlmets_gfl1vxmix,
                      mdlmets_gfl1amix,
                      mdlmets_gfl1vxbmix,
                      mdlmets_gfl1bmix,
                      mdlmets_gfl1vxamix,
                      mdlmets_gfl1cmix,
                      mdlmets_gfl1vxcmix)
mdlmets_gfl1 <- melt(mdlmets_gfl1, id.vars = c("type", "exp", "mets"))
mdlmets_gfl1$fls <- rep("One flight line", nrow(mdlmets_gfl1))


mdlmets_gfl2 <- rbind(mdlmets_gfl2allmix, 
                      mdlmets_gfl2vxallmix,
                      mdlmets_gfl2abmix,
                      mdlmets_gfl2vxabmix,
                      mdlmets_gfl2acmix,
                      mdlmets_gfl2vxacmix,
                      mdlmets_gfl2bcmix,
                      mdlmets_gfl2vxbcmix)
mdlmets_gfl2 <- melt(mdlmets_gfl2, id.vars = c("type", "exp", "mets"))
mdlmets_gfl2$fls <- rep("Two flight lines", nrow(mdlmets_gfl2))

mdlmets_g <- rbind(mdlmets_gfl1, mdlmets_gfl2)
mdlmets_g <- mdlmets_g[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                     ifelse(exp=="all", "allclasses 2 FL",
                                            ifelse(exp=="ab", "mostly ab",
                                                   ifelse(exp=="bc", "mostly bc",
                                                          ifelse(exp=="ac", "mostly ac",
                                                                 ifelse(exp=="a", "mostly a", 
                                                                        ifelse(exp=="b", "mostly b",
                                                                               ifelse(exp=="c", "mostly c", exp))))))))]

mdlmets_g <- mdlmets_g[, variable:=ifelse(variable=="r2", "R2",
                                          ifelse(variable=="rmse", "RMSE",
                                                 ifelse(variable=="mae", "MAE",
                                                        ifelse(variable=="mpe", "MPE", variable))))]


mdlmets_g$fls <- as.factor(mdlmets_g$fls)
mdlmets_g$variable <- factor(mdlmets_g$variable, levels=c("R2", "MPE", "RMSE", "MAE"))


ggplot(data=mdlmets_g, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Mixte: Basal area", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~fls, scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_mix[for_attr=="Basal area" & forest=="mix"], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))

###############################################################################################################
mdlmets_vtotfl1 <- rbind(mdlmets_vtotfl1mix, 
                         mdlmets_vtotfl1vxmix,
                         mdlmets_vtotfl1amix,
                         mdlmets_vtotfl1vxbmix,
                         mdlmets_vtotfl1bmix,
                         mdlmets_vtotfl1vxamix,
                         mdlmets_vtotfl1cmix,
                         mdlmets_vtotfl1vxcmix)
mdlmets_vtotfl1 <- melt(mdlmets_vtotfl1, id.vars = c("type", "exp", "mets"))
mdlmets_vtotfl1$fls <- rep("One flight line", nrow(mdlmets_vtotfl1))


mdlmets_vtotfl2 <- rbind(mdlmets_vtotfl2allmix, 
                         mdlmets_vtotfl2vxallmix,
                         mdlmets_vtotfl2abmix,
                         mdlmets_vtotfl2vxabmix,
                         mdlmets_vtotfl2acmix,
                         mdlmets_vtotfl2vxacmix,
                         mdlmets_vtotfl2bcmix,
                         mdlmets_vtotfl2vxbcmix)
mdlmets_vtotfl2 <- melt(mdlmets_vtotfl2, id.vars = c("type", "exp", "mets"))
mdlmets_vtotfl2$fls <- rep("Two flight lines", nrow(mdlmets_vtotfl2))

mdlmets_vtot <- rbind(mdlmets_vtotfl1, mdlmets_vtotfl2)
mdlmets_g <- mdlmets_g[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                     ifelse(exp=="all", "allclasses 2 FL",
                                            ifelse(exp=="ab", "mostly ab",
                                                   ifelse(exp=="bc", "mostly bc",
                                                          ifelse(exp=="ac", "mostly ac",
                                                                 ifelse(exp=="a", "mostly a", 
                                                                        ifelse(exp=="b", "mostly b",
                                                                               ifelse(exp=="c", "mostly c", exp))))))))]

mdlmets_vtot <- mdlmets_vtot[, variable:=ifelse(variable=="r2", "R2",
                                                ifelse(variable=="rmse", "RMSE",
                                                       ifelse(variable=="mae", "MAE",
                                                              ifelse(variable=="mpe", "MPE", variable))))]


mdlmets_vtot$fls <- as.factor(mdlmets_vtot$fls)
mdlmets_vtot$variable <- factor(mdlmets_vtot$variable, levels=c("R2", "MPE", "RMSE", "MAE"))


ggplot(data=mdlmets_vtot, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Mixte: Total volume", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~fls, scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_mix[for_attr=="Total volume" & forest=="mix"], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))


#############################################################################################################
mdlmets_vtigfl1 <- rbind(mdlmets_vtigfl1mix, 
                         mdlmets_vtigfl1vxmix,
                         mdlmets_vtigfl1amix,
                         mdlmets_vtigfl1vxbmix,
                         mdlmets_vtigfl1bmix,
                         mdlmets_vtigfl1vxamix,
                         mdlmets_vtigfl1cmix,
                         mdlmets_vtigfl1vxcmix)
mdlmets_vtigfl1 <- melt(mdlmets_vtigfl1, id.vars = c("type", "exp", "mets"))
mdlmets_vtigfl1$fls <- rep("One flight line", nrow(mdlmets_vtigfl1))


mdlmets_vtigfl2 <- rbind(mdlmets_vtigfl2allmix, 
                         mdlmets_vtigfl2vxallmix,
                         mdlmets_vtigfl2abmix,
                         mdlmets_vtigfl2vxabmix,
                         mdlmets_vtigfl2acmix,
                         mdlmets_vtigfl2vxacmix,
                         mdlmets_vtigfl2bcmix,
                         mdlmets_vtigfl2vxbcmix)
mdlmets_vtigfl2 <- melt(mdlmets_vtigfl2, id.vars = c("type", "exp", "mets"))
mdlmets_vtigfl2$fls <- rep("Two flight lines", nrow(mdlmets_vtigfl2))

mdlmets_vtig <- rbind(mdlmets_vtigfl1, mdlmets_vtigfl2)
mdlmets_g <- mdlmets_g[, exp:=ifelse(exp=="allclasses", "allclasses 1 FL",
                                     ifelse(exp=="all", "allclasses 2 FL",
                                            ifelse(exp=="ab", "mostly ab",
                                                   ifelse(exp=="bc", "mostly bc",
                                                          ifelse(exp=="ac", "mostly ac",
                                                                 ifelse(exp=="a", "mostly a", 
                                                                        ifelse(exp=="b", "mostly b",
                                                                               ifelse(exp=="c", "mostly c", exp))))))))]

mdlmets_vtig <- mdlmets_vtig[, variable:=ifelse(variable=="r2", "R2",
                                                ifelse(variable=="rmse", "RMSE",
                                                       ifelse(variable=="mae", "MAE",
                                                              ifelse(variable=="mpe", "MPE", variable))))]


mdlmets_vtig$fls <- as.factor(mdlmets_vtig$fls)
mdlmets_vtig$variable <- factor(mdlmets_vtig$variable, levels=c("R2", "MPE", "RMSE", "MAE"))

xx$variable <- factor(xx$variable, levels=c("R2", "MPE", "RMSE", "MAE"))



ggplot(data=xx, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "Feuillus: Basal area", x="Experiment", y="Value")+
  theme_base()+
  facet_grid(variable~., scales = "free")+
  scale_fill_grey(start = .4, end = .8)+
  geom_hline(data = mdlmets_ref1_mix[for_attr=="Total volume" & forest=="mix" & id %in% c("vtot_all_old_mix", "vtot_all_vox_mix")], 
             aes(yintercept=value, colour=exp, linetype=mets), size=1)+
  scale_colour_manual(values = c("green4","coral1", "magenta1"))+
  theme(legend.key.size = unit(1, 'cm'))
