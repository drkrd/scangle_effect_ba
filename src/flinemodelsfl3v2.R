library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggthemes)
library(gstat)
library(caret)
library(parallel)
library(foreach)
library(doParallel)

plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/flightlines_3/")
height=5
####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS
opt_independent_files(plots) <- TRUE
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
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber, ht = height)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber, ht = height)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht= height)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht= height)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
plotmetsfl3 <- catalog_apply(plots, func_computeall)
pmetsfl3 <- plotmetsfl3
pmetsfl3 <- rbindlist(pmetsfl3)
pmetsfl3 <- pmetsfl3[, c("fl1","fl2","fl3"):=tstrsplit(meanang,"_",fixed=T),]
pmetsfl3$fl1 <- as.numeric(pmetsfl3$fl1)
pmetsfl3$fl2 <- as.numeric(pmetsfl3$fl2)
pmetsfl3$fl3 <- as.numeric(pmetsfl3$fl3)

pmetsfl3 <- pmetsfl3[, fl2:= ifelse(is.na(fl2), fl1, fl2),]
pmetsfl3 <- pmetsfl3[, fl3:= ifelse(is.na(fl3), fl2, fl3),]



pmetsfl3 <- pmetsfl3[, cl1:=ifelse(fl1>=0&fl1<10, "a",
                                         ifelse(fl1>=10&fl1<20, "b",
                                                ifelse(fl1>=20&fl1<30, "c",
                                                       ifelse(fl1>=30&fl1<40,"d","e"))))]
pmetsfl3 <- pmetsfl3[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                         ifelse(fl2>=10&fl2<20, "b",
                                                ifelse(fl2>=20&fl2<30, "c",
                                                       ifelse(fl2>=30&fl2<40,"d","e"))))]
pmetsfl3 <- pmetsfl3[, cl3:=ifelse(fl3>=0&fl3<10, "a",
                                         ifelse(fl3>=10&fl3<20, "b",
                                                ifelse(fl3>=20&fl3<30, "c",
                                                       ifelse(fl3>=30&fl3<40,"d","e"))))]



pmetsfl3 <- pmetsfl3[cl1 != "e" & cl2 != "e" & cl3 != "e" & cl1 != "d" & cl2 != "d" & cl3 != "d",]




pmetsfl3 <- pmetsfl3[!fl1 %in% c(37.24, 7.16, 4.60, 16.52) & 
                       !fl2 %in% c(37.24, 7.16, 4.60, 16.52) &
                       !fl2 %in% c(37.24, 7.16, 4.60, 16.52)]




pmetsfl3 <- pmetsfl3[, cl := paste0(sort(.SD), collapse = ""), .SDcols = c("cl1", "cl2", "cl3"),  by = 1:nrow(pmetsfl3)]


##############################################################################################################################
allpcs <- list.files(paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/"),
                     pattern = "*.las",
                     full.names = TRUE)

alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)
names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))
#################################################################
allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/2_Ciron/voxelisation/Results/May/wo_interpolation/voxfiles/flightlines_3/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
allvoxfiles <- allvoxfiles[, c("fl1","fl2","fl3"):=tstrsplit(meanang,"_",fixed=T),]
allvoxfiles$fl1 <- as.numeric(allvoxfiles$fl1)
allvoxfiles$fl2 <- as.numeric(allvoxfiles$fl2)
allvoxfiles$fl3 <- as.numeric(allvoxfiles$fl3)

allvoxfiles <- allvoxfiles[, fl2:= ifelse(is.na(fl2), fl1, fl2),]
allvoxfiles <- allvoxfiles[, fl3:= ifelse(is.na(fl3), fl2, fl3),]

allvoxfiles <- allvoxfiles[, cl1:=ifelse(fl1>=0&fl1<10, "a",
                                         ifelse(fl1>=10&fl1<20, "b",
                                                ifelse(fl1>=20&fl1<30, "c",
                                                       ifelse(fl1>=30&fl1<40,"d","e"))))]
allvoxfiles <- allvoxfiles[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                         ifelse(fl2>=10&fl2<20, "b",
                                                ifelse(fl2>=20&fl2<30, "c",
                                                       ifelse(fl2>=30&fl2<40,"d","e"))))]
allvoxfiles <- allvoxfiles[, cl3:=ifelse(fl3>=0&fl3<10, "a",
                                         ifelse(fl3>=10&fl3<20, "b",
                                                ifelse(fl3>=20&fl3<30, "c",
                                                       ifelse(fl3>=30&fl3<40,"d","e"))))]


allvoxfiles <- allvoxfiles[cl1 != "e" & cl2 != "e" & cl3 != "e" & cl1 != "d" & cl2 != "d" & cl3 != "d",]



allvoxfiles <- allvoxfiles[!fl1 %in% c(37.24, 7.16, 4.60, 16.52) & 
                             !fl2 %in% c(37.24, 7.16, 4.60, 16.52) &
                             !fl2 %in% c(37.24, 7.16, 4.60, 16.52)]



#################################################################

voxall <- rbindlist(apply(allvoxfiles, 1, 
                          func_normvox2, 
                          pth="D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/", 
                          ht=height))
setDT(voxall)

pfcvladvox <- voxall[, .(cvladvox=cv(PADmean, na.rm = TRUE), 
                         sdvfp=sqrt(sum(PADmean*(k1-(sum(k1*PADmean)/sum(PADmean)))^2)/(sum(PADmean)*(length(PADmean[which(PADmean!=0)])-1)/length(PADmean[which(PADmean!=0)]))),
                         pfsumprof=exp(-0.5*sum(PADmean, na.rm = TRUE))), by=.(id_placette, meanang)]

setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(pmetsfl3, c("id_placette", "meanang"))
pmetsfl3 <- pmetsfl3[pfcvladvox]
pmetsfl3 <- pmetsfl3[which(id_placette %in% ciron_db$id_placette)]
pmetsfl3 <- pmetsfl3[!pflidr==0]
#######################################################################

######################################
#############All################
#############################################################################################################################
{
  pmetsfl3 <- unique(pmetsfl3[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl3 <- pmetsfl3[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl3 <- pmetsfl3[, wt := wt/sum(wt), id_placette]
  pmetsfl3$meanang <- as.factor(pmetsfl3$meanang)
  pmetsfl3 <- pmetsfl3[pflidr!=0]
  pmetsfl3.all <- pmetsfl3
  
}


###########################################################################################################
##Generate samples######
########################
dbase <- pmetsfl3.all
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
smplstfl3.all <- foreach(i = 1:10000, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  set.seed(i)
  dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
  return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
  # return(sset)
}
stopCluster(clus)
registerDoSEQ()
smplstfl3.all <- matrix(unlist(smplstfl3.all), nrow = length(unique(dbase$id_placette)))
smplstfl3.all <- unique(as.data.table(t(smplstfl3.all)))




time_log <- data.frame()


start <- Sys.time()
dbase <- pmetsfl3.all
fds <- ciron_db
setkey(fds,"id_placette")
idx.lst <- smplstfl3.all

f1l <- log(g_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(g_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(vtot_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(vtot_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(vtige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(vtige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


##Simulations
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
n <- 5000
cironfl3.all <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  idx <- as.vector(unlist(idx.lst[i]))
  mets_for_model <- dbase[idx]
  dbase$id_placette <- as.factor(dbase$id_placette)
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fds[mets_for_model]
  mets1 <- mets_for_model[,c("g_m2_ha","meanch", "varch", 'pflidr', "cvladlidr")]
  
  # y1 <- log(mets1$G_m2_ha)
  # x1 <- log(mets1$meanch)
  # x2 <- log(mets1$varch)
  # x3 <- log(mets1$pflidr)
  # x4 <- log(mets1$cvladlidr)
  
  # func_linreg <- function(b0 ,b1, b2, b3, b4, sig)
  # {
  #   ypred <- b0+b1*x1+b2*x2+b3*x3+b4*x4
  #   -sum(dnorm(y1, mean = ypred, sd = sig, log=TRUE))
  # }
  # 
  
  m1l <- train(f1l,
               data = mets1,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  G.l.coeff <- m1l$finalModel$coefficients
  G.l.pred <- m1l$pred$pred
  G.l.obs <- m1l$pred$obs
  
  
  mets2 <- mets_for_model[,c("g_m2_ha","meanch", "varch", 'pfsumprof', "cvladvox")]
  m1v <- train(f1v,
               data = mets2,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  G.v.coeff <- m1v$finalModel$coefficients
  G.v.pred <- m1v$pred$pred
  G.v.obs <- m1v$pred$obs
  ##########################################################
  
  ##########################################################
  mets3 <- mets_for_model[,c("vtot_m3_ha", "meanch", "varch", 'pflidr', "cvladlidr")]
  m2l <- train(f2l,
               data = mets3,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtot.l.coeff <- m2l$finalModel$coefficients
  vtot.l.pred <- m2l$pred$pred
  vtot.l.obs <- m2l$pred$obs
  
  mets4 <- mets_for_model[,c("vtot_m3_ha", "meanch", "varch", 'pfsumprof', "cvladvox")]
  m2v <- train(f2v,
               data = mets4,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtot.v.coeff <- m2v$finalModel$coefficients
  vtot.v.pred <- m2v$pred$pred
  vtot.v.obs <- m2v$pred$obs
  ##########################################################
  
  ###########################################################
  mets5 <- mets_for_model[,c("vtige_m3_ha", "meanch", "varch", 'pflidr', "cvladlidr")]
  m3l <- train(f3l,
               data = mets5,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtig.l.coeff <- m3l$finalModel$coefficients
  vtig.l.pred <- m3l$pred$pred
  vtig.l.obs <- m3l$pred$obs
  
  mets6 <- mets_for_model[,c("vtige_m3_ha", "meanch", "varch", 'pfsumprof', "cvladvox")]
  m3v <- train(f3v,
               data = mets6,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtig.v.coeff <- m3v$finalModel$coefficients
  vtig.v.pred <- m3v$pred$pred
  vtig.v.obs <- m3v$pred$obs
  ############################################################
  ############################################################
  
  x <- (list("G.lidr.pred" = G.l.pred,
             "G.lidr.obs" = G.l.obs,
             "G.vox.pred" = G.v.pred,
             "G.vox.obs" = G.v.obs,
             "vtot.lidr.pred" = vtot.l.pred,
             "vtot.lidr.obs" = vtot.l.obs,
             "vtot.vox.pred" = vtot.v.pred,
             "vtot.vox.obs" = vtot.v.obs,
             "vtig.lidr.pred" = vtig.l.pred,
             "vtig.lidr.obs" = vtig.l.obs,
             "vtig.vox.pred" = vtig.v.pred,
             "vtig.vox.obs" = vtig.v.obs,
             "G.l.coeff" = G.l.coeff,
             "G.v.coeff" = G.v.coeff,
             "vtot.l.coeff" = vtot.l.coeff,
             "vtot.v.coeff" = vtot.v.coeff,
             "vtig.l.coeff" = vtig.l.coeff,
             "vtig.v.coeff" = vtig.v.coeff
  ))
  return(x)
}
stopCluster(clus)
registerDoSEQ()

stop <- Sys.time()
time_log <- rbind(time_log, c(n, (stop-start)))




func_mdlmets <- function(obs, pred, for_attr, mettype)
{
  yobs <- exp(obs)
  ypred <- exp(pred)
  see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
  cf <- exp((see^2)/2)
  ypred <- ypred*cf
  SSE <- sum((yobs-ypred)^2)
  SST <- sum((mean(yobs)-yobs)^2)
  R2 <- 1-(SSE/SST)
  aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
  MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
  RMSE <- sqrt(mean((yobs-ypred)^2))
  RMSEpc <- RMSE*100/mean(yobs)
  return(list( "Forest_attr"=for_attr, "Metrics"=mettype, "R2"=R2, "RMSE"=RMSE,"rRMSE"=RMSEpc,"MPE"=MPE))
}


cironfl3.mdlmets.all <- melt(rbindlist(lapply(cironfl3.all, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))

cironfl3.mdlmets.all <- cbind(cironfl3.mdlmets.all, 
                              "exp"=rep("fl3", nrow(cironfl3.mdlmets.all)), 
                              "id"=rep(rep(1:5000, 1, each=6), 4))



ggplot(data=baugesfl3.mdlmets.feu.all[Forest_attr=="Stem volume"], aes( y=value, colour=Metrics))+
  geom_freqpoly()+
  facet_grid(variable~., scales = "free")+
  theme_base()+
  scale_fill_grey()




