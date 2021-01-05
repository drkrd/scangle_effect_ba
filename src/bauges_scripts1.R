plac74_con <- right_join(plac74_con, dd, by = "Id_plac")
plac74_con <- plac74_con[which(plac74_con$Id_plac %in% names(pc74_con)),]
mod74_con <- lm(log(G175)~log(meanch)+log(varch)+log(pf)+log(cvlad), data = plac74_con)
summary(mod74_con)



lasc <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/17_5m_rad/")


for(name in names(aois))
{
  aoi <- aois[[name]]
  aoi <- retrieve_flightlines(aoi)
  flist <- unique(aoi@data$flightlineID)
  df_fl_ang <- data.frame()
  
  
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aoi, flightlineID == fl)
    ar <- lidR::area(aoi_subset)
    meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
    df_fl_ang <- rbind(df_fl_ang,c(meangle, fl, ar))
  }
  colnames(df_fl_ang) <- c("meangle", "fl", "ar")
  print(df_fl_ang)
  df_fl_ang <- df_fl_ang[which(df_fl_ang[,3]>850),]
  print(df_fl_ang)
  
  for(nfl in df_fl_ang[,2])
  {
    aoi_sub <- lasfilter(aoi, flightlineID == nfl)
    ang <- df_fl_ang[which(df_fl_ang[,2]==nfl),1]
    nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/17_5m_rad/fl",
                 "/", name, "_", toString(round(ang,2)), ".las")
    writeLAS(aoi_sub, nm)
  }
}