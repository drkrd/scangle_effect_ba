
library(lidR)
library(dplyr)
library(tidyverse)
library(caret)

list_mets <- list()
list_models <- list()
list_modelscv <- list()
df_all_one <- data.frame()
df_meangle1 <- data.frame()
df_meangle2 <- data.frame()
for(loop in 1:1000)
{

  combinations <- 1
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
    
    l <- length(meangle_fl[,1])
    if(l>combinations) 
    {
      x <- sample((1:l), combinations)
      flid1 <- meangle_fl[,2][x][1]
      flid2 <- meangle_fl[,2][x][2]
    }
   
    else
    {
      x <- sample((1:l), l)
      flid1 <- meangle_fl[,2][x][1]
      flid2 <- 0
    }
 
    
    aoi_tmp <- filter_poi(aoi, flightlineID == flid1 | flightlineID == flid2)

    
    aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    
  
    mets_dbs_rand <- as.data.frame(rbind(mets_dbs_rand, c(name,
                                                          meangle_fl[,1][x][1],
                                                          meangle_fl[,1][x][2],
                                                          aoi_meanch,
                                                          aoi_varch,
                                                          aoi_pf,
                                                          aoi_cvlad)))
    
  

  }
  
  
  names(mets_dbs_rand) <- c("id_placette", "mean_angle1", "mean_angle2", "meanch", "varch", "pf", "cvlad")
  mets_dbs_rand$id_placette <- as.factor(mets_dbs_rand$id_placette)
  mets_dbs_rand$mean_angle1 <- as.numeric(mets_dbs_rand$mean_angle1)
  mets_dbs_rand$mean_angle2 <- as.numeric(mets_dbs_rand$mean_angle2)
  mets_dbs_rand$meanch <- as.numeric(mets_dbs_rand$meanch)
  mets_dbs_rand$varch <- as.numeric(mets_dbs_rand$varch)
  mets_dbs_rand$pf <- as.numeric(mets_dbs_rand$pf)
  mets_dbs_rand$cvlad <- as.numeric(mets_dbs_rand$cvlad)
  
  
  
  mets_dbs_rand <- right_join(fd_ba_smry, mets_dbs_rand, by="id_placette")

  #sd_meangle <- sd(mets_dbs_rand$mean_angle)
  
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
  

  
  
  df_all_one <- as.data.frame(rbind(df_all_one, c(r2, exp(rmser), exp(maer), r2cv, 
                                              exp(rmsercv), exp(maercv),
                                              abs(mets_dbs_rand$mean_angle1-mets_dbs_rand$mean_angle2))))
  length(df_all_one)
  df_meangle1 <- as.data.frame(rbind(df_meangle1, mets_dbs_rand$mean_angle1))
  df_meangle2 <- as.data.frame(rbind(df_meangle2, mets_dbs_rand$mean_angle2))
  list_mets[[loop]] <- as.list(mets_dbs_rand[,c(1,8:11)])
  list_models[[loop]] <- model_randang
  list_modelscv[[loop]] <- model_randang
  

  length(df_meangle1)
  
  
  print(loop)
}

colnames(df_all_one) <- c("r2", "rmser", "maer", "r2cv", "rmsercv", "maercv",
                        "d1",  "d2",  "d3",  "d4",  "d5",  "d6",  "d7",  "d8",  "d9", "d9.2",
                        "d10", "d11", "d12", "13", "d14", "d15", "d16", "d17", "d18", "d19", 
                        "d20", "d21", "d22", "d23", "d24", "d25", "d26", "d27", "d28", "d29")

colnames(df_meangle1) <- c(1,  2,  3,  4,  5,  6,  7,  8,  9, 9.2,
                        10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                        20, 21, 22, 23, 24, 25, 26, 27, 28, 29)

colnames(df_meangle2) <- c(1,  2,  3,  4,  5,  6,  7,  8,  9, 9.2,
                           10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 
                           20, 21, 22, 23, 24, 25, 26, 27, 28, 29)
