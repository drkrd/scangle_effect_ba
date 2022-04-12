##heatmap

##combine ciron and bauges74 information

plotmetsfl1$stratum <- "Riparian"
plotmetsfl1 <- plotmetsfl1[, Plot:= .GRP, by = .(id_placette)]

pmetsfl1.con.all$stratum <- "Coniferous"
pmetsfl1.con.all <- pmetsfl1.con.all[, Plot:= .GRP, by = .(id_placette)]

pmetsfl1.feu.all$stratum <- "Broadleaved"
pmetsfl1.feu.all <- pmetsfl1.feu.all[, Plot:= .GRP, by = .(id_placette)]

pmetsfl1.mix.all$stratum <- "Mixed"
pmetsfl1.mix.all <- pmetsfl1.mix.all[, Plot:= .GRP, by = .(id_placette)]



db <- rbind(plotmetsfl1[,c("id_placette", "cl", "stratum", "Plot")], 
            pmetsfl1.con.all[,c("id_placette", "cl", "stratum", "Plot")],
            pmetsfl1.feu.all[,c("id_placette", "cl", "stratum", "Plot")],
            pmetsfl1.mix.all[,c("id_placette", "cl", "stratum", "Plot")])
#there are 499 flight lines in total (Bauges + Ciron) that were used in the study
db <- db[, .("N"=.N), by=c("stratum", "id_placette", "cl", "Plot")]
db$N <- as.factor(db$N)
db$Plot <- as.factor(db$Plot)

db <- db[, cl:=ifelse(cl=="a", "A",
                      ifelse(cl=="b", "B","C"))]

db$stratum <- factor(db$stratum, levels = c("Riparian", "Coniferous",
                                            "Broadleaved", "Mixed"))

hmap <- ggplot(data=db)+
  aes(x=Plot, y=cl, fill=N)+
  scale_fill_few()+
  geom_tile(colour="white", stat="identity", width=0.9, height=0.9)+
  facet_wrap(stratum~., scales = "free", nrow = 4)+
  theme_base()+
  theme(text=element_text(family="serif", size=9*(96/72)),
        strip.text.x = element_text(size = 11*(96/72)),
        strip.background = element_blank(),
        legend.position = "top",
        panel.spacing = unit(0, "lines"))+
  labs(x= "Plot nos.", y="Class")


ggsave(hmap, file="D:/1_Work/Dropbox/2_Publications/2_paper/revisions/figures/hmap.png", dpi=320) 

