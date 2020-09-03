library(dplyr)
library(conflicted)
mets_dbs_min <- data.frame()

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
    meangle_fl <- rbind(meangle_fl,c(meangle, fl, ar))
  }
  
  
  meangle_fl <- meangle_fl[meangle_fl[,3]>650,]
  flid <- meangle_fl[,2][which.min(meangle_fl[,1])]
  aoi_tmp <- lasfilter(aoi, flightlineID == flid)
  
  
  aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
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
mets_dbs_min$meanch <- as.numeric(mets_dbs_min$meanch)
mets_dbs_min$varch <- as.numeric(mets_dbs_min$varch)
mets_dbs_min$pf <- as.numeric(mets_dbs_min$pf)
mets_dbs_min$cvlad <- as.numeric(mets_dbs_min$cvlad)
mets_dbs_min <- right_join(fd_ba_smry, mets_dbs_min, by="id_placette")

model_minang <- lm(data = mets_dbs_min[-c(21,22),], 
                   formula = log(sum_ba_hec)~log(meanch)+log(varch)+log(pf)+log(cvlad))
