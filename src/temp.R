

aois_ref13 <- list()
for(name in names(aois))
{
  aoi <- lasflightline(aois[[name]])
  flist <- unique(aoi@data$flightlineID)
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aoi, flightlineID == fl)
    meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
    if(all.equal(meangle,df_row13[[name]])==TRUE)
    {
      print(meangle)
      aois_ref13[[name]] <- aoi_subset
    }
  }
}









func_fline <- function(angs_df)
{
  aois_tmp <- list()
  for(name in names(aois))
  {
    aoi <- lasflightline(aois[[name]])
    flist <- unique(aoi@data$flightlineID)
    for(fl in flist)
    {
      aoi_subset <- lasfilter(aoi, flightlineID == fl)
      meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
      if(all.equal(meangle,angs_df[[name]])==TRUE)
      {

        aois_tmp[[name]] <- aoi_subset
      }
    }
  }
  return(aois_tmp)
}









df_all11 <- data.frame()
mets_dbs_checks <- data.frame()



for (name in names(aois_ref))
{
  aoi <- lasflightline(aois[[name]])
  flist <- unique(aoi@data$flightlineID)
  if(length(list_angs[[name]])==1)
  {
    new_angs <- ref_df[[name]]
  }
  else
  {
    new_angs <- setdiff(list_angs[[name]], c(ref_df[[name]]))
  }
  print(new_angs)
  ref_df2 <- ref_df
  for(new_ang in new_angs)
  {
    mets_dbs_checks <- data.frame()
    ref_df2[[name]] <- new_ang
    aois_new <- func_fline(ref_df2)
    for(name_new in names(aois_new))
    {
      aoi_tmp <- aois_new[[name_new]]
      meangle <- abs(mean(aoi_tmp@data$ScanAngleRank))
      aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
      aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
      aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
      aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
      mets_dbs_checks <- as.data.frame(rbind(mets_dbs_checks, c(name_new,
                                                                meangle,
                                                                aoi_meanch,
                                                                aoi_varch,
                                                                aoi_pf,
                                                                aoi_cvlad)))
    }
    names(mets_dbs_checks) <- c("id_placette", "mean_angle", "meanch", "varch", "pf", "cvlad")
    mets_dbs_checks$id_placette <- as.factor(mets_dbs_checks$id_placette)
    mets_dbs_checks$mean_angle <- as.numeric(mets_dbs_checks$mean_angle)
    mets_dbs_checks$meanch <- as.numeric(mets_dbs_checks$meanch)
    mets_dbs_checks$varch <- as.numeric(mets_dbs_checks$varch)
    mets_dbs_checks$pf <- as.numeric(mets_dbs_checks$pf)
    mets_dbs_checks$cvlad <- as.numeric(mets_dbs_checks$cvlad)
    
    
    mets_dbs_checks <- right_join(fd_ba_smry, mets_dbs_checks, by="id_placette")
    model_checksang_loocv <- train(log(sum_ba_hec)~log(meanch)+log(varch)+log(pf)+log(cvlad),
                                 data = mets_dbs_checks[-c(21,22),],
                                 method = "lm",
                                 trControl = trainControl(method="LOOCV"))
    
    
    model_checksang <- lm(data = mets_dbs_checks[-c(21,22),], 
                        formula = log(sum_ba_hec)~log(meanch)+log(varch)+log(pf)+log(cvlad))
    
    r2 <- summary(model_checksang)$adj.r.squared
    rmser <- rmse(log(mets_dbs_checks$sum_ba_hec[-c(21,22)]), predict(model_checksang))
    maer <- mean(abs(log(mets_dbs_checks$sum_ba_hec[-c(21,22)])-predict(model_checksang)))
    r2cv <- model_checksang_loocv$results$Rsquared
    rmsercv <- model_checksang_loocv$results$RMSE
    maercv <- model_checksang_loocv$results$MAE
    
    df_all11 <- rbind(df_all11, c(r2, exp(rmser), exp(maer), r2cv, exp(rmsercv), exp(maercv), mets_dbs_checks$mean_angle))
    print("done")
  }
}



























mets71 <- data.frame()
for(name_new in names(aois_ref71))
{
  aoi_tmp <- aois_ref71[[name_new]]
  meangle <- abs(mean(aoi_tmp@data$ScanAngleRank))
  aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  mets71 <- as.data.frame(rbind(mets71, c(name_new,
                                          meangle,
                                          aoi_meanch,
                                          aoi_varch,
                                          aoi_pf,
                                          aoi_cvlad)))
}


names(mets13) <- c("id_placette", "meangle", "meanch", "varch", "pf", "cvlad")
