
library(lidR)
library(dplyr)
library(tidyverse)
library(caret)
library(ModelMetrics)

list_mets_mi <- list()
list_models_mi <- list()
list_modelscv_mi <- list()
df_all_one_mi <- data.frame()

ar_th <- 0.9*pi*17.5*17.5
#df_meangle <- data.frame()
aois_mi <- aois[which(names(aois) %in% df_mi$Id_plac)]
for(loop in 1:300)
{

  ncombs <- 1
  df_randangs <- data.frame()
  min_miunt=0
  
  for(name in names(aois_mi))
  {

    aoi <- aois[[name]]
    flist <- unique(aoi@data$flightlineID)
    df_fl_ang <- data.frame()
    
    
    for(fl in flist)
    {
      aoi_subset <- lasfilter(aoi, flightlineID == fl)
      ar <- lidR::area(aoi_subset)
      meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
      df_fl_ang <- rbind(df_fl_ang,c(meangle, fl, ar))
    }
    
    
    df_fl_ang <- df_fl_ang[which(df_fl_ang[,3]>ar_th),]
    
    #get the number of flight lines that cover the area entirely 
    l <- length(df_fl_ang[,1])
    #if three flight lines are available, it is possible to make 3c1 (one fl),>> 
    #3c2 (two fl) or 3c3 (three fl) etc. configs. 
    
    #there may be some plots where only one fline exists but the experiment may be>>
    #to test various combinations of flines. In such cases, the available fline data is used>>
    #along with the combinations of flines on other plots.
    
    #also, it is necessary to specify this because there will be an error for nCr, where n<r
    
    #if condition is for n>r where there will be multiple possible combinations
    if(l>ncombs) {
      #generate a combination of numbers from a list (1:l) of numbers
      #these combination of numbers are basically the index values using which we will randomly>>
      #pick flines from the df_fl_ang in each iteration
      x <- sample((1:l), ncombs) 
      flids <- df_fl_ang[,2][x]
      aoi_tmp <- filter_poi(aoi, flightlineID %in% flids)
    }else 
    {
      #else condition is for n<r or n==r where only one possible combinations exists.
      #we will pick all the available flight lines
      x <- sample((1:l), ncombs)
      flids <- df_fl_ang[,2]
      aoi_tmp <- filter_poi(aoi, flightlineID %in% flids)
    }
    
    
    z <- aoi_tmp@data$Z
    rn <- aoi_tmp@data$ReturnNumber
    
    mch <- func_meanch(z, rn)
    vch <- func_varch(z, rn)
    pf <- 1-func_pf(z, rn)
    cvl <- func_cvlad(z, rn)

    df_randangs <- as.data.frame(rbind(df_randangs, c(name,
                                                      df_fl_ang[,1][x],
                                                      mch,
                                                      vch,
                                                      pf,
                                                      cvl)))
  }

  
  names(df_randangs) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
  df_randangs$id_placette <- as.factor(df_randangs$id_placette)
  df_randangs$mean_angle <- as.numeric(df_randangs$mean_angle)
  df_randangs$meanch <- as.numeric(df_randangs$meanch)
  df_randangs$varch <- as.numeric(df_randangs$varch)
  df_randangs$pf <- as.numeric(df_randangs$pf)
  df_randangs$cvlad <- as.numeric(df_randangs$cvlad)
  
  
  
  mets_randangs_mi <- right_join(df_mi[,c(1,48)], 
                            df_randangs, 
                            by=c("Id_plac"="id_placette"))
  #mets_randangs_mi <- mets_randangs_mi[complete.cases(mets_randangs_mi),]
  
  #df_randangs <- df_randangs[-which(df_randangs$pf==0),]

  #sd_meangle <- sd(df_randangs$mean_angle)
  
  mdl_randangs_loocv_mi <- train(log(G175)~log(meanch)+log(varch)+log(pf)+log(cvlad),
                         data = mets_randangs_mi,
                         method = "lm",
                         trControl = trainControl(method="LOOCV"))
  
  
  mdl_randangs_mi <- lm(data = mets_randangs_mi, 
                      formula = log(G175)~log(meanch)+log(varch)+log(pf)+log(cvlad))
  


  r2 <- summary(mdl_randangs_mi)$adj.r.squared
  rmser <- rmse(log(mets_randangs_mi[, 2]), predict(mdl_randangs_mi))
  maer <- mean(abs(log(mets_randangs_mi[, 2])-predict(mdl_randangs_mi)))
  r2cv <- mdl_randangs_loocv_mi$results$Rsquared
  rmsercv <- mdl_randangs_loocv_mi$results$RMSE
  maercv <- mdl_randangs_loocv_mi$results$MAE
  

  
  
  df_all_one_mi <- as.data.frame(rbind(df_all_one_mi, c(r2, exp(rmser), exp(maer), r2cv, 
                                              exp(rmsercv), exp(maercv), 
                                              c(mets_randangs_mi$mean_angle))))
  #length(df_all_one_mi)
  #df_meangle <- as.data.frame(rbind(df_meangle1, df_randangs$mean_angle1))

  list_mets_mi[[loop]] <- list(mets_randangs_mi[,c(1:7)])
  list_models_mi[[loop]] <- mdl_randangs_mi
  list_modelscv_mi[[loop]] <- mdl_randangs_loocv_mi
  
  print(loop)
}


colnames(df_all_one_mi) <- c("r2", "rmser", "maer", "r2cv", "rmsercv", "maercv", names(aois_mi))

