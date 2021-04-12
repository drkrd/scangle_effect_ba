library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)

typefor <- ""
typepc <- "allpoints"
# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
# plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
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
plotmetsfl2 <- plotmetsfl2[,c("m1","m2"):=tstrsplit(meanang,"_",fixed=T),]
plotmetsfl2 <- plotmetsfl2[,class1:=ifelse(m1>=0&m1<10, "a",
                                           ifelse(m1>=10&m1<20, "b",
                                                  ifelse(m1>=20&m1<30, "c",
                                                         ifelse(m1>=30&m1<40,"d","e"))))]
plotmetsfl2 <- plotmetsfl2[,class2:=ifelse(m2>=0&m2<10, "a",
                                           ifelse(m2>=10&m2<20, "b",
                                                  ifelse(m2>=20&m2<30, "c",
                                                         ifelse(m2>=30&m2<40,"d","e"))))]
setkey(plotmetsfl1,"id_placette")



keycols = c("id_placette", "meanang")
setkeyv(plotmetsfl1, keycols)

plotmetsfl1 <- plotmetsfl1[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                        ifelse(meanang>=10&meanang<20, "b",
                                               ifelse(meanang>=20&meanang<30, "c",
                                                      ifelse(meanang>=30&meanang<40,"d","e"))))]

plotmetsfl1 <- plotmetsfl1[cl!="e"]

plotmetsfl$meanang <- as.factor(plotmetsfl$meanang)
setkey(plotmets, "id")

fd_smry <- read.csv("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",")
fd_smry$id_placette <- as.factor(fd_smry$id_placette)
setDT(fd_smry)
setkey(fd_smry, "id_placette")

mets_for_model <- fd_smry[plotmets]

######################################
#Models for all points################
######################################
library(caret)
mets_for_model <- fd_smry[plotmets]
mdl_allpoints <- list()
mdl_allpoints[["G"]] <- train(log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                                data = mets_for_model,
                                method = "lm",
                                trControl = trainControl(method="LOOCV"))

mdl_allpoints[["V_tige"]] <- train(log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                                    data = mets_for_model,
                                    method = "lm",
                                    trControl = trainControl(method="LOOCV"))

mdl_allpoints[["V_tot"]] <- train(log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                                data = mets_for_model,
                                method = "lm",
                                trControl = trainControl(method="LOOCV"))


#####################################################################
#Several iterations of models with random flight line chosen per plot
#####################################################################



library(parallel)
library(foreach)
library(doParallel)
n <- 2000
set.seed(123)


clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdl_vtig1 <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl1[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fd_smry[mets_for_model]
  
  
  model <- train(log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  

  x <- (list("model"=model,
              "meanang"=c(mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()




mets_mdl_vtig1 <- lapply(mdl_vtig1, function(x)
{
  r2 <- x$model$results$Rsquared
  rmse <- x$model$results$RMSE
  mae <- x$model$results$MAE
  
  return(list("r2"=r2, "rmse"=rmse, "mae"=mae, "flightlines"="any one"))
})
mets_mdl_vtig1 <- rbindlist(mets_mdl_vtig1)



mets_mdl_vtig1 <- reshape2::melt(mets_mdl_vtig1, id.vars=c("flightlines"))


cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(data=mets_mdl_vtig1, aes(x=value, fill=flightlines))+
  geom_histogram(bins=60)+
  facet_wrap(~variable, scales = "free")+
  labs(title = "Stem volume")+
  theme_tufte(base_size = 13)+
  scale_fill_manual(values = cbp1)+
  geom_vline(data=filter(mets_mdl_vtig1, variable=="r2"), aes(xintercept=0.84), colour="red") + 
  geom_vline(data=filter(mets_mdl_vtig1, variable=="rmse"), aes(xintercept=0.16), colour="blue") + 
  geom_vline(data=filter(mets_mdl_vtig1, variable=="mae"), aes(xintercept=0.128), colour="green")
  





  
ggplot(data=mets_mdl_vtig12, aes(x=value, fill=flightlines))+
  geom_boxplot()+
  facet_grid(~variable, scales = "free")+
  theme_tufte(base_size = 13)+
  scale_fill_manual(values = cbp1)


  
  
  
  
  
  
  



