
library(lidR)
library(dplyr)
library(tidyverse)
library(caret)

sd_lst <- list(1000)
r2_lst <- list(1000)
high_lst <- list(1000)
df_all <- data.frame()
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

  sd_meangle <- sd(mets_dbs_rand$mean_angle)
  
  model_randang_loocv <- train(log(sum_ba_hec)~log(meanch)+log(varch)+log(pf)+log(cvlad),
                         data = mets_dbs_rand,
                         method = "lm",
                         trControl = trainControl(method="LOOCV"))
  
  
  model_randang <- lm(data = mets_dbs_rand, 
                      formula = log(sum_ba_hec)~log(meanch)+log(varch)+log(pf)+log(cvlad))
  


  r2 <- summary(model_randang)$adj.r.squared
  rmser <- rmse(log(mets_dbs_rand$sum_ba_hec), predict(model_randang))
  maer <- mean(abs(log(mets_dbs_rand$sum_ba_hec)-predict(model_randang)))
  r2cv <- model_randang_loocv$results$Rsquared
  rmsercv <- model_randang_loocv$results$RMSE
  maercv <- model_randang_loocv$results$MAE
  
  
  
  df_all <- rbind(df_all, c(r2, exp(rmser), exp(maer), r2cv, 
                            exp(rmsercv), exp(maercv), mets_dbs_rand$mean_angle))
  
  #ang_lst[[toString(r2)]] <- mets_dbs_rand$mean_angle
  
  
  
  

  #print(r2)
  #variation[i] <- vars
  sd_lst[[loop]] <- sd_meangle
  r2_lst[[loop]] <- r2
  print(loop)
}

colnames(df_all11) <- c("r2", "rmser", "maer", "r2cv", "rmsercv", "maercv",
                      1,  2,  3,  4,  5,  6,  7,  8,  9, 9.2,
                      10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                      20, 21, 22, 23, 24, 25, 26, 27, 28, 29)