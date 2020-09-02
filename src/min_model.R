min_dbs <- data.frame()
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
  flid <- meangle_fl[,2][which.min(meangle_fl[,1])]
  aoi_tmp <- lasfilter(aoi, flightlineID == flid)
  aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  min_dbs <- as.data.frame(rbind(min_dbs, c(name,
                                            min(meangle_fl[,1]),
                                            aoi_meanch,
                                            aoi_varch,
                                            aoi_pf,
                                            aoi_cvlad)))
}
names(min_dbs) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
min_dbs$mean_angle <- as.numeric(min_dbs$mean_angle)
min_dbs$meanch <- log(as.numeric(min_dbs$meanch))
min_dbs$varch <- log(as.numeric(min_dbs$varch))
min_dbs$pf <- log(as.numeric(min_dbs$pf))
min_dbs$cvlad <- log(as.numeric(min_dbs$cvlad))
min_dbs2 <- right_join(fdata, min_dbs, by="id_placette")
min_model <- lm(sum_ba_hec~meanch+varch+pf+cvlad, data=min_dbs2[-c(21,22),])
summary(min_model)
