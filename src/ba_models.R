variation <- numeric(100)
sd_lst <- numeric(100)
r2_lst <- numeric(100)
high_lst <- list()
for(i in 1:1000)
{
  #print(i)
  rand_dbs <- data.frame()
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
    rand_dbs <- as.data.frame(rbind(rand_dbs, c(name,
                                                meangle_fl[,1][x],
                                                aoi_meanch,
                                                aoi_varch,
                                                aoi_pf,
                                                aoi_cvlad)))
  }
  
  
  names(rand_dbs) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
  rand_dbs$mean_angle <- as.numeric(rand_dbs$mean_angle)
  rand_dbs$meanch <- log(as.numeric(rand_dbs$meanch))
  rand_dbs$varch <- log(as.numeric(rand_dbs$varch))
  rand_dbs$pf <- log(as.numeric(rand_dbs$pf))
  rand_dbs$cvlad <- log(as.numeric(rand_dbs$cvlad))
  
  
  
  fdata1 <- right_join(fdata, rand_dbs, by="id_placette")
  
  # fdata1 <- fdata1 %>% 
  #   mutate(Classes=if_else(mean_angle<5,1,
  #                          if_else(mean_angle<10,2,
  #                                  if_else(mean_angle<15,3,
  #                                          if_else(mean_angle<20,4,
  #                                                  if_else(mean_angle<25,5,
  #                                                          if_else(mean_angle<30,6,
  #                                                                  if_else(mean_angle<35,7,
  #                                                                          if_else(mean_angle<40,8,9)))))))))
  
  
  # vars <- length(unique(fdata1$Classes))
  
  
  sd_meangle <- sd(fdata1$mean_angle)
  model <- lm(sum_ba_hec~meanch+varch+pf+cvlad, data = fdata1[-c(21,21),])
  
  r2 <- summary(model)$adj.r.squared
  high_lst[[toString(r2)]] <- fdata1$mean_angle

  #print(sd_meangle)
  variation[i] <- vars
  sd_lst[i] <- sd_meangle
  r2_lst[i] <- r2
  
  
}
