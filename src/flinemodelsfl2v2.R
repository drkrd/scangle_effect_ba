library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)
library(caret)
library(parallel)
library(foreach)
library(doParallel)

plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/flightlines_2/")
height=2
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
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber, ht=height)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber, ht=height)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht=height)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht=height)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
plotmetsfl2 <- catalog_apply(plots, func_computeall)
#############################################################################################################
plotmetsfl2 <- rbindlist(plotmetsfl2)
plotmetsfl2 <- plotmetsfl2[, c("fl1","fl2"):=tstrsplit(meanang,"_",fixed=T),]
plotmetsfl2$fl1 <- as.numeric(plotmetsfl2$fl1)
plotmetsfl2$fl2 <- as.numeric(plotmetsfl2$fl2)
plotmetsfl2 <- plotmetsfl2[, fl2:= ifelse(is.na(fl2), fl1, fl2),]

plotmetsfl2 <- plotmetsfl2[, cl1:=ifelse(fl1>=0&fl1<10, "a",
                                           ifelse(fl1>=10&fl1<20, "b",
                                                  ifelse(fl1>=20&fl1<30, "c",
                                                         ifelse(fl1>=30&fl1<40,"d","e"))))]
plotmetsfl2 <- plotmetsfl2[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                           ifelse(fl2>=10&fl2<20, "b",
                                                  ifelse(fl2>=20&fl2<30, "c",
                                                         ifelse(fl2>=30&fl2<40,"d","e"))))]

plotmetsfl2 <- plotmetsfl2[cl1 != "e" & cl2 != "e" & cl1 != "d" & cl2 != "e",]

plotmetsfl2 <- plotmetsfl2[!fl1 %in% c(37.24, 7.16, 4.60, 16.52) & !fl2 %in% c(37.24, 7.16, 4.60, 16.52)]


plotmetsfl2 <- plotmetsfl2[, cl := paste0(sort(.SD), collapse = ""), .SDcols = c("cl1", "cl2"),  by = 1:nrow(plotmetsfl2)]


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

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/2_Ciron/voxelisation/Results/May/wo_interpolation/voxfiles/flightlines_2/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
allvoxfiles <- allvoxfiles[, c("fl1","fl2"):=tstrsplit(meanang,"_",fixed=T),]
allvoxfiles$fl1 <- as.numeric(allvoxfiles$fl1)
allvoxfiles$fl2 <- as.numeric(allvoxfiles$fl2)
allvoxfiles <- allvoxfiles[, cl1:=ifelse(fl1>=0&fl1<10, "a",
                                         ifelse(fl1>=10&fl1<20, "b",
                                                ifelse(fl1>=20&fl1<30, "c",
                                                       ifelse(fl1>=30&fl1<40,"d","e"))))]
allvoxfiles <- allvoxfiles[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                         ifelse(fl2>=10&fl2<20, "b",
                                                ifelse(fl2>=20&fl2<30, "c",
                                                       ifelse(fl2>=30&fl2<40,"d","e"))))]

allvoxfiles <- allvoxfiles[cl1 != "e" & cl2 != "e" & cl1 != "d" & cl2 != "e",]
allvoxfiles <- allvoxfiles[!fl1 %in% c(37.24, 7.16, 4.60, 16.52) & !fl2 %in% c(37.24, 7.16, 4.60, 16.52)]


voxall <- rbindlist(apply(allvoxfiles, 1,
                          func_normvox2, 
                          pth="D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/", 
                          ht=height))
setDT(voxall)

pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), 
                         pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang)]


setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(plotmetsfl2, c("id_placette", "meanang"))
plotmetsfl2 <- plotmetsfl2[pfcvladvox]
#############################################################################################################################
plotmetsfl2 <- unique(plotmetsfl2[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl2 <- plotmetsfl2[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl2 <- plotmetsfl2[, wt := wt/sum(wt), id_placette]
plotmetsfl2$rid <- seq(1:nrow(plotmetsfl2)) 
#############################################################################################################################

plotmetsfl2ab=plotmetsfl2[plotmetsfl2[, .I[cl=="ab"  | all(cl!="ab")], by = id_placette]$V1]
plotmetsfl2ac=plotmetsfl2[plotmetsfl2[, .I[cl=="ac"  | all(cl!="ac")], by = id_placette]$V1]
plotmetsfl2bc=plotmetsfl2[plotmetsfl2[, .I[cl=="bc"  | all(cl!="bc")], by = id_placette]$V1]


l1 <- unique(as.character(plotmetsfl2ab[cl=="ab"]$id_placette))
l2 <- unique(as.character(plotmetsfl2ac[cl=="ac"]$id_placette))
l3 <- unique(as.character(plotmetsfl2bc[cl=="bc"]$id_placette))

cps <- intersect(intersect(l1, l2), l3)

plotmetsfl2ab=plotmetsfl2ab[plotmetsfl2ab[, all(cl != 'ab')| (cl == 'ab' & .BY %in% cps)|!.BY %in% cps, 
                                     by = id_placette]$V1]
plotmetsfl2ab <- unique(plotmetsfl2ab[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl2ab <- plotmetsfl2ab[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl2ab <- plotmetsfl2ab[, wt := wt/sum(wt), id_placette]


plotmetsfl2ac=plotmetsfl2ac[plotmetsfl2ac[, all(cl != 'ac')| (cl == 'ac' & .BY %in% cps)|!.BY %in% cps, 
                                          by = id_placette]$V1]
plotmetsfl2ac <- unique(plotmetsfl2ac[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl2ac <- plotmetsfl2ac[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl2ac <- plotmetsfl2ac[, wt := wt/sum(wt), id_placette]


plotmetsfl2bc=plotmetsfl2bc[plotmetsfl2bc[, all(cl != 'bc')| (cl == 'bc' & .BY %in% cps)|!.BY %in% cps, 
                                          by = id_placette]$V1]
plotmetsfl2bc <- unique(plotmetsfl2bc[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl2bc <- plotmetsfl2bc[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl2bc <- plotmetsfl2bc[, wt := wt/sum(wt), id_placette]
####################################################################################################################

#################################################################################################################################
pmetsfl2.all <- plotmetsfl2

#tabulate per plot the number of pcs belonging to each class
tbl <- with(plotmetsfl2, table(id_placette, cl))

#pick all the plots that have atleast one pc belonging to class a, class b and class c each
tbl2 <- tbl[which(tbl[,2]>0 & tbl[,3]>0 & tbl[,6]>0),]

#subset those pcs and their metrics which satisfy the above condition
pmets.withallcls <- pmetsfl2.all[id_placette %in% rownames(tbl2)]

#subset those pcs and their metrics which donot satisfy the previous condition
pmets.woallcls <- pmetsfl2.all[!id_placette %in% rownames(tbl2)]





#pick only those pcs and their metrics which belong to class a
pmetsfl2.clab <- pmets.withallcls[cl=="ab"]

#combine with the remaining plots
pmetsfl2.clab <- rbind(pmetsfl2.clab, pmets.woallcls)

#compute inverse probabilities
pmetsfl2.clab <- unique(pmetsfl2.clab[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl2.clab <- pmetsfl2.clab[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl2.clab <- pmetsfl2.clab[, wt := wt/sum(wt), id_placette]





#pick only those pcs and their metrics which belong to class b
pmetsfl2.clbc <- pmets.withallcls[cl=="bc"]

#combine with the remaining plots
pmetsfl2.clbc <- rbind(pmetsfl2.clbc, pmets.woallcls)

#compute inverse probabilities
pmetsfl2.clbc <- unique(pmetsfl2.clbc[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl2.clbc <- pmetsfl2.clbc[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl2.clbc <- pmetsfl2.clbc[, wt := wt/sum(wt), id_placette]





#pick only those pcs and their metrics which belong to class b
pmetsfl2.clac <- pmets.withallcls[cl=="ac"]

#combine with the remaining plots
pmetsfl2.clac <- rbind(pmetsfl2.clac, pmets.woallcls)

#compute inverse probabilities
pmetsfl2.clac <- unique(pmetsfl2.clac[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl2.clac <- pmetsfl2.clac[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl2.clac <- pmetsfl2.clac[, wt := wt/sum(wt), id_placette]



###########################################################################################################
###########################################################################################################
##Generate samples######
########################
dbase <- pmetsfl2.all
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
smplst.2all <- foreach(i = 1:10000, .packages=c("dplyr", "data.table", "caret", "sampling")) %dopar% {
  set.seed(i)
  dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
  return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
}
stopCluster(clus)
registerDoSEQ()
smplst.2all <- matrix(unlist(smplst.2all), nrow = 29)
smplst.2all <- unique(as.data.table(t(smplst.2all)))





fd_smry <- fread("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",", drop = "id_placette")
colnames(fd_smry)[1] <- "id_placette"
fd_smry[, id_placette := as.factor(id_placette)]
setkey(fd_smry, "id_placette")



time_log <- data.frame()

start <- Sys.time()
dbase <- pmetsfl2.all
fds <- fd_smry[!id_placette=="14"]
idx.lst <- smplst.2all

f1l <- log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(G_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


##Simulations
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
n <- 1000
ciron.2all <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  idx <- as.vector(unlist(idx.lst[i]))
  mets_for_model <- dbase[idx]
  dbase$id_placette <- as.factor(dbase$id_placette)
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fds[mets_for_model]
  mets1 <- mets_for_model[,c("G_m2_ha","meanch", "varch", 'pflidr', "cvladlidr")]
  
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
  
  
  mets2 <- mets_for_model[,c("G_m2_ha","meanch", "varch", 'pfsumprof', "cvladvox")]
  m1v <- train(f1v,
               data = mets2,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  G.v.coeff <- m1v$finalModel$coefficients
  G.v.pred <- m1v$pred$pred
  G.v.obs <- m1v$pred$obs
  ##########################################################
  
  ##########################################################
  mets3 <- mets_for_model[,c("volume_total_m3_ha", "meanch", "varch", 'pflidr', "cvladlidr")]
  m2l <- train(f2l,
               data = mets3,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtot.l.coeff <- m2l$finalModel$coefficients
  vtot.l.pred <- m2l$pred$pred
  vtot.l.obs <- m2l$pred$obs
  
  mets4 <- mets_for_model[,c("volume_total_m3_ha", "meanch", "varch", 'pfsumprof', "cvladvox")]
  m2v <- train(f2v,
               data = mets4,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtot.v.coeff <- m2v$finalModel$coefficients
  vtot.v.pred <- m2v$pred$pred
  vtot.v.obs <- m2v$pred$obs
  ##########################################################
  
  ###########################################################
  mets5 <- mets_for_model[,c("volume_tige_m3_ha", "meanch", "varch", 'pflidr', "cvladlidr")]
  m3l <- train(f3l,
               data = mets5,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtig.l.coeff <- m3l$finalModel$coefficients
  vtig.l.pred <- m3l$pred$pred
  vtig.l.obs <- m3l$pred$obs
  
  mets6 <- mets_for_model[,c("volume_tige_m3_ha", "meanch", "varch", 'pfsumprof', "cvladvox")]
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
  see <- sqrt(sum((obs-pred)^2)/(length(obs)-4))
  cf <- exp((see^2)/2)
  ypred <- ypred*cf
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
  return(list( "Forest_attr"=for_attr, "Metrics"=mettype, "R2"=R2, "MPE"=MPE, "RMSE"=RMSE,"RMSEpc"=RMSEpc))
}


ciron.mdlmets.2all <- melt(rbindlist(lapply(ciron.2all, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "old")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "old")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "old")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "MPE", "RMSE", "RMSEpc"))

ciron.mdlmets.all1 <- melt(rbindlist(lapply(ciron.all5, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "old")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "old")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "old")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "MPE", "RMSE", "MAE"))

ciron.mdlmets.cla <- melt(rbindlist(lapply(ciron.cla, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "old")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "old")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "old")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "MPE", "RMSE", "MAE"))
ciron.mdlmets.clb <- melt(rbindlist(lapply(ciron.clb, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "old")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "old")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "old")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "MPE", "RMSE", "MAE"))
ciron.mdlmets.clc <- melt(rbindlist(lapply(ciron.clc, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "old")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "old")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "old")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "MPE", "RMSE", "MAE"))

# 
ciron.mdlmets.all <- cbind(ciron.mdlmets.2all, "exp"=rep("2all", nrow(ciron.mdlmets.2all)), "id"=rep(rep(1:1000, 1, each=6), 4))
ciron.mdlmets.all1 <- cbind(ciron.mdlmets.all1, "exp"=rep("all", nrow(ciron.mdlmets.all)), "id"=rep(rep(1:5000, 1, each=6), 4))

ciron.mdlmets.cla <- cbind(ciron.mdlmets.cla, "exp"=rep("mostly A", nrow(ciron.mdlmets.cla)), "id"=rep(rep(1:5000, 1, each=6), 4))
ciron.mdlmets.clb <- cbind(ciron.mdlmets.clb, "exp"=rep("mostly B", nrow(ciron.mdlmets.clb)), "id"=rep(rep(1:5000, 1, each=6), 4))
ciron.mdlmets.clc <- cbind(ciron.mdlmets.clc, "exp"=rep("mostly C", nrow(ciron.mdlmets.clc)), "id"=rep(rep(1:5000, 1, each=6), 4))



# 
ciron.mdlmets <- as.data.table(rbind(ciron.mdlmets.all, ciron.mdlmets.cla, ciron.mdlmets.clb, ciron.mdlmets.clc))

ciron.2mdlmets <- as.data.table(cbind(ciron.mdlmets.all, "fl"=rep("fl2", nrow(ciron.mdlmets.all)))) 

ciron.mets <- rbind(ciron.1mdlmets, ciron.2mdlmets)

ggplot(data=ciron.mets[Forest_attr=="Total volume"], aes(y=value, colour=Metrics))+
  geom_boxplot()+
  facet_grid(variable~fl, scales = "free")+
  theme_base()+
  scale_fill_grey()





























































n <- 5000
set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_gfl2bcpfcvlvx <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl2bc[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fd_smry[mets_for_model]
  
  
  model <- train(log(G_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  
  x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                "RMSE" = model$results$RMSE,
                                "MAE" = model$results$MAE),
             "coeffs" = list(model$finalModel$coefficients),
             "index" = which(plotmetsfl2bc$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()

set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_vtotfl2bcpfcvlvx <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl2bc[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fd_smry[mets_for_model]
  
  
  model <- train(log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  
  x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                "RMSE" = model$results$RMSE,
                                "MAE" = model$results$MAE),
             "coeffs" = list(model$finalModel$coefficients),
             "index" = which(plotmetsfl2bc$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()

set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_vtigfl2bcpfcvlvx <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl2bc[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fd_smry[mets_for_model]
  
  
  model <- train(log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  
  x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                "RMSE" = model$results$RMSE,
                                "MAE" = model$results$MAE),
             "coeffs" = list(model$finalModel$coefficients),
             "index" = which(plotmetsfl2bc$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()








func_extractmdlmets <- function(x, type, exp, mets)
{
  ind <- x[["index"]]
  # dbase <- plotmetsfl1[ind]
  # dbase <- dbase[!id_placette%in%c(21, 22, 26)]
  # ent <- func_entropy(dbase$cl)
  r2 <- x[["modelmets"]]$R2 
  rmse <- x[["modelmets"]]$RMSE
  mae <- x[["modelmets"]]$MAE
  type <- type
  exp <- exp
  mets <- mets
  # names(tbl[which(tbl==max(tbl))])
  return(list("r2"=r2, "rmse"=rmse, "mae"=mae, "type"=type, "exp"=exp, "mets"=mets))
}


mdlmets_vtigfl2bc <- unique(rbindlist(lapply(mdls_vtigfl2bc, 
                                        func_extractmdlmets, 
                                        type="Two flight lines", 
                                        exp="mostly BC", 
                                        mets="old")))

mdlmets_g <- rbind(mdlmets_gfl2, 
                   mdlmets_gfl2pfcvlvx,
                   mdlmets_gfl2ab,
                   mdlmets_gfl2abpfcvlvx,
                   mdlmets_gfl2ac,
                   mdlmets_gfl2acpfcvlvx,
                   mdlmets_gfl2bc,
                   mdlmets_gfl2bcpfcvlvx)
mdlmets_g <- melt(mdlmets_g, id.vars = c("type", "exp", "mets"))
mdlmets_g$for_attr <- rep("Basal area", nrow(mdlmets_g))

mdlmets_vtot <- rbind(mdlmets_vtotfl2, 
                      mdlmets_vtotfl2pfcvlvx,
                      mdlmets_vtotfl2ab,
                      mdlmets_vtotfl2abpfcvlvx,
                      mdlmets_vtotfl2ac,
                      mdlmets_vtotfl2acpfcvlvx,
                      mdlmets_vtotfl2bc,
                      mdlmets_vtotfl2bcpfcvlvx)
mdlmets_vtot <- melt(mdlmets_vtot, id.vars = c("type", "exp", "mets"))
mdlmets_vtot$for_attr <- rep("Total volume", nrow(mdlmets_vtot))

mdlmets_vtig <- rbind(mdlmets_vtigfl2, 
                      mdlmets_vtigfl2pfcvlvx,
                      mdlmets_vtigfl2ab,
                      mdlmets_vtigfl2abpfcvlvx,
                      mdlmets_vtigfl2ac,
                      mdlmets_vtigfl2acpfcvlvx,
                      mdlmets_vtigfl2bc,
                      mdlmets_vtigfl2bcpfcvlvx)
mdlmets_vtig <- melt(mdlmets_vtig, id.vars = c("type", "exp", "mets"))
mdlmets_vtig$for_attr <- rep("Stem volume", nrow(mdlmets_vtig))


mdlmetsfl2 <- rbind(mdlmets_g, mdlmets_vtot, mdlmets_vtig)



mdlmets <- rbind(mdlmetsfl1, mdlmetsfl2)

ggplot(data=mdlmets[for_attr=="Basal area"], aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "")+
  theme_base(base_size = 20)+
  facet_grid(variable~type, scales = "free")+
  geom_hline(data = allmdlmets, aes(yintercept=value, linetype=type, colour=mets), size=1.5)+
  guides(linetype = guide_legend(override.aes = list(size=1)))+
  scale_fill_pander()











