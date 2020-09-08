#variation <- numeric(100)
library(lidR)
library(dplyr)
sd_lst <- list(1000)
r2_lst <- list(1000)
high_lst <- list(1000)
for(loop in 1:1000)
{

  
  mets_dbs_rand <- data.frame()
  min_count=0
  
  
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
    x <- sample(1:length(meangle_fl[,1]),1)
    min_fl <- min(meangle_fl[,1])
    max_fl <- max(meangle_fl[,1])
    flid <- meangle_fl[,2][x]
    if(meangle_fl[,1][x]==min_fl) {min_count=min_count+1}
    
    aoi_tmp <- lasfilter(aoi, flightlineID == flid)
    aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    mets_dbs_rand <- as.data.frame(rbind(mets_dbs_rand, c(name,
                                                          meangle_fl[,1][x],
                                                          aoi_meanch,
                                                          aoi_varch,
                                                          aoi_pf,
                                                          aoi_cvlad)))
  }
  
  
  names(mets_dbs_rand) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
  mets_dbs_rand$id_placette <- as.factor(mets_dbs_rand$id_placette)
  mets_dbs_rand$mean_angle <- as.numeric(mets_dbs_rand$mean_angle)
  mets_dbs_rand$meanch <- as.numeric(mets_dbs_rand$meanch)
  mets_dbs_rand$varch <- as.numeric(mets_dbs_rand$varch)
  mets_dbs_rand$pf <- as.numeric(mets_dbs_rand$pf)
  mets_dbs_rand$cvlad <- as.numeric(mets_dbs_rand$cvlad)
  
  
  
  mets_dbs_rand <- right_join(fd_ba_smry, mets_dbs_rand, by="id_placette")
  
  # fd_ba_smry1 <- fd_ba_smry1 %>% 
  #   mutate(Classes=if_else(mean_angle<5,1,
  #                          if_else(mean_angle<10,2,
  #                                  if_else(mean_angle<15,3,
  #                                          if_else(mean_angle<20,4,
  #                                                  if_else(mean_angle<25,5,
  #                                                          if_else(mean_angle<30,6,
  #                                                                  if_else(mean_angle<35,7,
  #                                                                          if_else(mean_angle<40,8,9)))))))))
  
  
  # vars <- length(unique(fd_ba_smry1$Classes))
  
  
  sd_meangle <- sd(mets_dbs_rand$mean_angle)
  model_randang <- lm(data = mets_dbs_rand[-c(21,22),], 
                      formula = log(sum_ba_hec)~log(meanch)+log(varch)+log(pf)+log(cvlad))
  
  r2 <- summary(model_randang)$adj.r.squared
  high_lst[[toString(r2)]] <- mets_dbs_rand$mean_angle

  #print(r2)
  #variation[i] <- vars
  sd_lst[[loop]] <- sd_meangle
  r2_lst[[loop]] <- r2
  print(loop)
}
