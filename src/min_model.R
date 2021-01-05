library(dplyr)
library(conflicted)
mets_dbs_min <- data.frame()
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
  flid <- meangle_fl[,2][which.min(meangle_fl[,1])]
  aoi_tmp <- lasfilter(aois[[name]], flightlineID == flid)
  
  
  aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_pf <- 1-func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  
  
  mets_dbs_min <- as.data.frame(rbind(mets_dbs_min, c(name,
                                                      min(meangle_fl[,1]),
                                                      aoi_meanch,
                                                      aoi_varch,
                                                      aoi_pf,
                                                      aoi_cvlad)))
}


names(mets_dbs_min) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
mets_dbs_min$id_placette <- as.factor(mets_dbs_min$id_placette)
mets_dbs_min$mean_angle <- as.numeric(mets_dbs_min$mean_angle)
mets_dbs_min$meanch <- as.numeric(mets_dbs_min$meanch)
mets_dbs_min$varch <- as.numeric(mets_dbs_min$varch)
mets_dbs_min$pf <- as.numeric(mets_dbs_min$pf)
mets_dbs_min$cvlad <- as.numeric(mets_dbs_min$cvlad)


mets_minang_co73 <- right_join(df_co73[,c(1,48)], 
                      mets_dbs_min, 
                      by=c("Id_plac"="id_placette"))


mets_minang_co73 <- mets_minang_co73[complete.cases(mets_minang_co73),]

mdl_minang_co73 <- lm(data = mets_minang_co73, 
                   formula = log(G175)~log(meanch)+log(varch)+log(pf)+log(cvlad))

summary(mdl_minang_co73)
