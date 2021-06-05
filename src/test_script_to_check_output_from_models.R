library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)

plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/march2021/allpoints/")

opt_independent_files(plots) <- TRUE ####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS

func_computeall <- function(chunk)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  
  
  id_plac <- sub("\\_n.*", "", basename(tools::file_path_sans_ext(chunk@files)))
  mang <- sub(".*@", "", basename(tools::file_path_sans_ext(chunk@files)))
  if(!is.na(as.numeric(mang)))
  {
    mang <- round(as.numeric(mang),2)

  }
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber)
  return(list( 
    id_placette = as.factor(id_plac),
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}

plotmets <- catalog_apply(plots, func_computeall)
plotmets <- rbindlist(plotmets)
setkey(plotmets,"id_placette")
fd_smry <- read.csv("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",")
fd_smry$id_placette <- as.factor(fd_smry$id_placette)
setDT(fd_smry)
setkey(fd_smry, "id_placette")
library(caret)
mets_for_model <- fd_smry[plotmets]
mdl_allpoints <- list()
mdl_allpoints[["G"]] <- train(log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                              data = mets_for_model,
                              method = "lm",
                              trControl = trainControl(method="LOOCV"))


r2 <- c()
se <- c()
for(i in 1:nrow(mets_for_model))
{
  metstr <- mets_for_model[-i]
  metsval <- mets_for_model[i]
  mdl <- lm(log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr), data = metstr)
  s <- summary(mdl)
  ypred <- s[["coefficients"]][1]+
    s[["coefficients"]][2]*log(metsval$meanch)+
    s[["coefficients"]][3]*log(metsval$varch)+
    s[["coefficients"]][4]*log(metsval$pflidr)+
    s[["coefficients"]][5]*log(metsval$cvladlidr)
  se <- c(se, abs(ypred-log(metsval$G_m2_ha)))
}




mdl_allpoints[["vtot"]] <- train(log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                                     data = mets_for_model,
                                     method = "lm",
                                     trControl = trainControl(method="LOOCV"))

y <- train(log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
      data = mets_for_model,
      method = "lm",
      trControl = trainControl(method="LOOCV"))

mdl_allpoints[["vtige"]] <- train(log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                                    data = mets_for_model,
                                    method = "lm",
                                    trControl = trainControl(method="LOOCV"))

z <- train(log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
      data = mets_for_model,
      method = "lm",
      trControl = trainControl(method="LOOCV"))


stats <- lapply(mdl_allpoints, function(x){
  return(list("r2"=x$results$Rsquared, "rmse"=x$results$RMSE, "mae"=x$results$MAE))
})


mets_mdl <- lapply(mdl_allpoints, function(x)
{
  s <- summary(x)
  r2 <- s$adj.r.squared
  rmse <- sqrt(mean((s$residuals)^2))
  mae <- mean(abs(s$residuals))
  return(list("r2"=r2, "rmse"=rmse, "mae"=mae))
})