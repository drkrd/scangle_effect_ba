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
library(stringr)

# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
# plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
bauges_db <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db15jul.csv", sep = ",")
colnames(bauges_db)[2] <- "id_placette"
height <- 5

bauges_db <- bauges_db[, newstratum := ifelse(comp_R_G>75 & comp_F_G<25 , "Coniferes",
                                              ifelse(comp_F_G>75 & comp_R_G<25, "Feuillus", "Mixte"))]


bauges_db[, stratum := as.factor(stratum)]
bauges_db[, newstratum := as.factor(newstratum)]
bauges_db[, id_placette := as.factor(id_placette)]
setkey(bauges_db, "id_placette")


# The following flightlines are problematic. There is no voxelisation for them
# 283_74_ONF_un@6.28
# 96_IRSTEA_un@30

plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/flightlines_1/")
opt_independent_files(plots) <- TRUE ####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS



height=5
func_computeall <- function(chunk)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  
  
  id_plac <- sub("\\_n.*", "", basename(tools::file_path_sans_ext(chunk@files)))
  mang <- sub(".*@", "", basename(tools::file_path_sans_ext(chunk@files)))
  if(!is.na(as.numeric(mang)))
  {
    mang <- round(as.numeric(mang),2)
    
  }
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber, ht = height)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber, ht = height)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht= height)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht= height)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
plotmetsfl1 <- catalog_apply(plots, func_computeall)
pmetsfl1 <- plotmetsfl1
pmetsfl1 <- rbindlist(pmetsfl1)

pmetsfl1 <- pmetsfl1[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                       ifelse(meanang>=10&meanang<20, "b",
                                              ifelse(meanang>=20&meanang<30, "c",
                                                     ifelse(meanang>=30&meanang<40,"d","e"))))]
pmetsfl1 <- pmetsfl1[cl!="d" & cl!="e"]
pmetsfl1 <- na.omit(pmetsfl1)
pmetsfl1 <- pmetsfl1[!(id_placette=="96_IRSTEA" & meanang==30.0)]
pmetsfl1 <- pmetsfl1[!(id_placette=="283_74_ONF" & meanang==6.28)]
pmetsfl1 <- pmetsfl1[!pflidr==0]
pmetsfl1.con.all <- pmetsfl1[which(id_placette %in% bauges_db_con$id_placette)]
pmetsfl1.feu.all <- pmetsfl1[which(id_placette %in% bauges_db_feu$id_placette)]
pmetsfl1.mix.all <- pmetsfl1[which(id_placette %in% bauges_db_mix$id_placette)]






bauges_db <- bauges_db[which(id_placette %in% pmetsfl1$id_placette), 
                   c("id_placette", "G75", "volume_total", "volume_tige", "newstratum")]
bauges_db_con <- bauges_db[newstratum=="Coniferes"]
bauges_db_feu <- bauges_db[newstratum=="Feuillus"]
bauges_db_mix <- bauges_db[newstratum=="Mixte"]





pmetsbc <- pmetsfl1.con.all
smry_bc <- pmetsbc[,.N, by=c("id_placette", "cl")]
smry_bc <- dcast(smry_bc, id_placette~cl, value.var ="N" )
smry_bc[is.na(smry_bc)] <- 0
smry_bc <- melt(smry_bc, id.vars = "id_placette")
smry_bc$value <- as.factor(smry_bc$value)
smry_bc$terrain <- rep("Mountainous", nrow(smry_bc))
smry_bc$forest <- rep("Coniferous", nrow(smry_bc))
smry_bc$location <- rep("Bauges", nrow(smry_bc))
p1 <- ggplot(data=smry_bc, aes(x=id_placette, y=variable, fill=value))+
  geom_tile(aes(width=0.9, height=0.9))+
  coord_equal()+
  scale_fill_grey()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none",
        plot.margin=unit(c(1,1,0,1), "cm"))


pmetsbf <- pmetsfl1.feu.all
smry_bf <- pmetsbf[,.N, by=c("id_placette", "cl")]
smry_bf <- dcast(smry_bf, id_placette~cl, value.var ="N" )
smry_bf[is.na(smry_bf)] <- 0
smry_bf <- melt(smry_bf, id.vars = "id_placette")
smry_bf$value <- as.factor(smry_bf$value)
smry_bf$terrain <- rep("Mountainous", nrow(smry_bf))
smry_bf$forest <- rep("Broadleaved", nrow(smry_bf))
smry_bf$location <- rep("Bauges", nrow(smry_bf))
p2 <- ggplot(data=smry_bf, aes(x=id_placette, y=variable, fill=value))+
  geom_tile(aes(width=0.9, height=0.9))+
  coord_equal()+
  scale_fill_grey()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none",
        plot.margin=unit(c(0,1,1,1), "cm"))

pmetsbm <- pmetsfl1.mix.all
smry_bm <- pmetsbm[,.N, by=c("id_placette", "cl")]
smry_bm <- dcast(smry_bm, id_placette~cl, value.var ="N" )
smry_bm[is.na(smry_bm)] <- 0
smry_bm <- melt(smry_bm, id.vars = "id_placette")
smry_bm$value <- as.factor(smry_bm$value)
smry_bm$terrain <- rep("Mountainous", nrow(smry_bm))
smry_bm$forest <- rep("Mixed", nrow(smry_bm))
smry_bm$location <- rep("Bauges", nrow(smry_bm))

p3 <- ggplot(data=smry_bm, aes(x=id_placette, y=variable, fill=value))+
  geom_tile(aes(width=0.9, height=0.9))+
  coord_equal()+
  scale_fill_grey()+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        legend.position="none",
        plot.margin=unit(c(-1,1,-1,1), "cm"))


smry_ciron <- pmetsx[,.N, by=c("id_placette", "cl")]
smry_ciron <- dcast(smry_ciron, id_placette~cl, value.var ="N" )
smry_ciron[is.na(smry_ciron)] <- 0
smry_ciron <- melt(smry_ciron, id.vars = "id_placette")
smry_ciron$value <- as.factor(smry_ciron$value)
smry_ciron$terrain <- rep("Flat", nrow(smry_ciron))
smry_ciron$forest <- rep("Riparian", nrow(smry_ciron))
smry_ciron$location <- rep("Ciron", nrow(smry_ciron))

p4 <- ggplot(data=smry_ciron, aes(x=id_placette, y=variable, fill=value))+
  geom_tile(aes(width=0.9, height=0.9))+
  coord_equal()+
  scale_fill_grey()+
  theme_bw()

legend <- get_legend(p4)

p4 <- p4 + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 legend.position="none",
                 plot.margin=unit(c(-1,1,1,1), "cm"))




grid.arrange(arrangeGrob(p1, p2, nrow=2))


smry <- as.data.table(rbind(smry_ciron, smry_bc, smry_bf, smry_bm))


library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

ggplot(data=smry, aes(x=id_placette, y=variable, fill=value))+
  geom_tile(aes(width=0.9, height=0.9))+
  theme_classic2()+
  labs(x="Plots", y="Classes")+
  theme(axis.text.x = element_blank())+
  scale_fill_grey()+
  facet_wrap(forest ~ ., scales =  'free', nrow=4)
  

                      