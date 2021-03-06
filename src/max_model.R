mets_dbs_max <- data.frame()


for(name in names(aois))
{
  aoi <- lasflightline(aois[[name]])
  flist <- unique(aoi@data$flightlineID)
  meangle_fl <- data.frame()
  
  
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aoi, flightlineID == fl)
    ar <- area_calc(data.frame(x=aoi_subset@data$X, y=aoi_subset@data$Y))
    meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
    meangle_fl <- rbind(meangle_fl,c(meangle,fl, ar))
  }
  
  
  meangle_fl <- meangle_fl[meangle_fl[,3]>650,]
  flid <- meangle_fl[,2][which.max(meangle_fl[,1])]
  aoi_tmp <- lasfilter(aoi, flightlineID == flid)
  
  
  aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  
  
  mets_dbs_max <- as.data.frame(rbind(mets_dbs_max, c(name,
                                                      max(meangle_fl[,1]),
                                                      aoi_meanch,
                                                      aoi_varch,
                                                      aoi_pf,
                                                      aoi_cvlad)))
}


names(mets_dbs_max) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
mets_dbs_max$id_placette <- as.factor(mets_dbs_max$id_placette)
mets_dbs_max$mean_angle <- as.numeric(mets_dbs_max$mean_angle)
mets_dbs_max$meanch <- as.numeric(mets_dbs_max$meanch)
mets_dbs_max$varch <- as.numeric(mets_dbs_max$varch)
mets_dbs_max$pf <- as.numeric(mets_dbs_max$pf)
mets_dbs_max$cvlad <- as.numeric(mets_dbs_max$cvlad)
mets_dbs_max <- right_join(fd_ba_smry, mets_dbs_max, by="id_placette")


model_maxang <- lm(data = mets_dbs_max[-c(21,22),], 
                   formula = log(sum_ba_hec)~log(meanch)+log(varch)+log(pf)+log(cvlad))





#########################################################################################"

library(dplyr)
library(conflicted)
mets_dbs_max <- data.frame()
ar_th <- 0.9*pi*17.5*17.5
for(name in names(aois))
{
  #aoi <- lasflightline(aois[[name]])
  flist <- unique(aois[[name]]@data$flightlineID)
  meangle_fl <- data.frame()
  
  
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aois[[name]], flightlineID == fl)
    ar <- area_calc(data.frame(x=aoi_subset@data$X, y=aoi_subset@data$Y))
    meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
    meangle_fl <- rbind(meangle_fl,c(meangle, fl, ar))
  }
  
  
  meangle_fl <- meangle_fl[meangle_fl[,3]>ar_th,]
  flid <- meangle_fl[,2][which.max(meangle_fl[,1])]
  aoi_tmp <- lasfilter(aois[[name]], flightlineID == flid)
  
  
  aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_pf <- 1-func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  
  
  mets_dbs_max <- as.data.frame(rbind(mets_dbs_max, c(name,
                                                      max(meangle_fl[,1]),
                                                      aoi_meanch,
                                                      aoi_varch,
                                                      aoi_pf,
                                                      aoi_cvlad)))
}


names(mets_dbs_max) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
mets_dbs_max$id_placette <- as.factor(mets_dbs_max$id_placette)
mets_dbs_max$mean_angle <- as.numeric(mets_dbs_max$mean_angle)
mets_dbs_max$meanch <- as.numeric(mets_dbs_max$meanch)
mets_dbs_max$varch <- as.numeric(mets_dbs_max$varch)
mets_dbs_max$pf <- as.numeric(mets_dbs_max$pf)
mets_dbs_max$cvlad <- as.numeric(mets_dbs_max$cvlad)


mets_maxang_mi73 <- right_join(df_mi73[,c(1,48)], 
                          mets_dbs_max, 
                          by=c("Id_plac"="id_placette"))

mets_maxang_mi73 <- mets_maxang_mi73[complete.cases(mets_maxang_mi73),]

mdl_maxang_mi73 <- lm(data = mets_maxang_mi73, 
                 formula = log(G175)~log(meanch)+log(varch)+log(pf)+log(cvlad))

summary(mdl_maxang_mi73)

