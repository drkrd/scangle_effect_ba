ciron.mdlmets.all <- rbind(cironfl1.mdlmets.all,
                           cironfl2.mdlmets.all,
                           cironfl3.mdlmets.all)
ciron.mdlmets.all$exp <- factor(ciron.mdlmets.all$exp,
                                levels = c("One", "Two", "Three"))



func_labeller <- function(variable, value)
{
  c("","Basal area","")
}

bp1 <- ggplot(data=ciron.mdlmets.all[variable!="RMSE"& Metrics=="ref"&Forest_attr=="Basal area"], aes(x=exp, y=value))+
  geom_boxplot(aes(fill=Metrics), outlier.size = 1)+
  facet_wrap(Forest_attr~variable, scales = "free", ncol = 4, 
             labeller = labeller(Forest_attr=as_labeller(func_labeller)))+
  theme_minimal()+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "none")+
  labs(x="Flight line combinations",
       y=NULL)+
  scale_fill_manual(values=c("white", "gray45"))

ggsave(bp1, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/Trial1.png", 
       width=9*1.25, height=5*1.25, units="cm", dpi=320) 




dta <- ciron.mdlmets.all[Metrics=="ref"]

dta <- dta[, .(vr=var(value), mn=mean(value)), by=.(Forest_attr, variable, exp)]

dta <- dcast(dta, Forest_attr+exp~variable, value.var="vr")









bauges.mdlmets.con.all <- rbind(baugesfl1.mdlmets.con.all,
                                baugesfl2.mdlmets.con.all,
                                baugesfl3.mdlmets.con.all)

ggplot(data=bauges.mdlmets.con.all, aes(x=fl, y=value, fill=Metrics))+
  geom_boxplot()+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_base()+
  labs(title="Coniferous",
       x="No. of flight lines",
       y="")+
  scale_fill_manual(values=c("gray10", "gray45"))





bauges.mdlmets.feu.all <- rbind(baugesfl1.mdlmets.feu.all,
                                baugesfl2.mdlmets.feu.all,
                                baugesfl3.mdlmets.feu.all)

ggplot(data=bauges.mdlmets.feu.all[Metrics=="ref"], aes(x=fl, y=value, fill=Metrics))+
  geom_boxplot()+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_base()+
  labs(title="Broadleaved",
       x="No. of flight lines",
       y="")+
  scale_fill_manual(values=c("gray10"))

bauges.mdlmets.mix.all <- rbind(baugesfl1.mdlmets.mix.all,
                                baugesfl2.mdlmets.mix.all,
                                baugesfl3.mdlmets.mix.all)

ggplot(data=bauges.mdlmets.mix.all[Metrics=="ref"], aes(x=fl, y=value, fill=Metrics))+
  geom_boxplot()+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_base()+
  labs(title="Mixed",
       x="No. of flight lines",
       y="")+
  scale_fill_manual(values=c("gray10"))
