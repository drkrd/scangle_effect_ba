max_dbs <- data.frame()
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
  max_dbs <- as.data.frame(rbind(max_dbs, c(name,
                                            max(meangle_fl[,1]),
                                            aoi_meanch,
                                            aoi_varch,
                                            aoi_pf,
                                            aoi_cvlad)))
}
names(max_dbs) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
max_dbs$mean_angle <- as.numeric(max_dbs$mean_angle)
max_dbs$meanch <- log(as.numeric(max_dbs$meanch))
max_dbs$varch <- log(as.numeric(max_dbs$varch))
max_dbs$pf <- log(as.numeric(max_dbs$pf))
max_dbs$cvlad <- log(as.numeric(max_dbs$cvlad))
max_dbs2 <- right_join(fdata, max_dbs, by="id_placette")
max_model <- lm(sum_ba_hec~meanch+varch+pf+cvlad, data=max_dbs2[-c(21,22),])
summary(max_model)
