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

# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
# plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/flightlines_1/")
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
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht=2)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht=2)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
plotmetsfl1 <- catalog_apply(plots, func_computeall)
plotmetsfl1 <- rbindlist(plotmetsfl1)
plotmetsfl1 <- plotmetsfl1[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                       ifelse(meanang>=10&meanang<20, "b",
                                              ifelse(meanang>=20&meanang<30, "c",
                                                     ifelse(meanang>=30&meanang<40,"d","e"))))]
plotmetsfl1 <- plotmetsfl1[cl!="e"]
plotmetsfl1 <- plotmetsfl1[!meanang %in% c(37.24, 7.16)]



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

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/2_Ciron/voxelisation/Results/May/wo_interpolation/voxfiles/flightlines_1/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","", basename(file_path_sans_ext(V1)))]
allvoxfiles$meanang <- round(as.numeric(allvoxfiles$meanang))
allvoxfiles[, meanang := sub(".*\\@","", basename(file_path_sans_ext(V1)))]
voxall <- rbindlist(apply(allvoxfiles, 1, 
                          func_normvox2, 
                          pth="D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/", 
                          ht=2))
setDT(voxall)
voxall <- voxall[!meanang %in% c(37.24, 7.16, 49.7, 43.93)]
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang, pfsumvox)]
setkeyv(pfcvladvox, c("id_placette", "meanang"))
pfcvladvox$id_placette <- sub("_.*","", pfcvladvox$id_placette)
pfcvladvox$meanang <- as.numeric(pfcvladvox)
plotmetsfl1$meanang <- as.character(plotmetsfl1$meanang)
setkeyv(plotmetsfl1, c("id_placette", "meanang"))
plotmetsfl1 <- plotmetsfl1[pfcvladvox]
#############################################################################################################################
plotmetsfl1 <- unique(plotmetsfl1[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1 <- plotmetsfl1[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1 <- plotmetsfl1[, wt := wt/sum(wt), id_placette]
plotmetsfl1$meanang <- as.numeric(plotmetsfl1$meanang)
##############################################################################################################################


plotmetsfl1a <- plotmetsfl1[plotmetsfl1[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
plotmetsfl1a <- unique(plotmetsfl1a[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1a <- plotmetsfl1a[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1a <- plotmetsfl1a[, wt := wt/sum(wt), id_placette]


plotmetsfl1b <- plotmetsfl1[plotmetsfl1[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
plotmetsfl1b <- unique(plotmetsfl1b[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1b <- plotmetsfl1b[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1b <- plotmetsfl1b[, wt := wt/sum(wt), id_placette]


plotmetsfl1c <- plotmetsfl1[plotmetsfl1[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]
plotmetsfl1c <- unique(plotmetsfl1c[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1c <- plotmetsfl1c[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1c <- plotmetsfl1c[, wt := wt/sum(wt), id_placette]


#################################################################################################################################
l1 <- unique(plotmetsfl1a[cl=="a"]$id_placette)
l2 <- unique(plotmetsfl1b[cl=="b"]$id_placette)
l3 <- unique(plotmetsfl1c[cl=="c"]$id_placette)

cps <- intersect(intersect(l1, l2), l3)

tbl <- with(plotmetsfl1, table(cl, id_placette))
tblx <- tbl[,which(tbl[1,]>0 & tbl[2,]>0 & tbl[3,]>0)]
pmets <- plotmetsfl1[id_placette %in% colnames(tblx)]
pmets1 <- plotmetsfl1[!id_placette %in% colnames(tblx)]


pmetsfl1.cla <- pmets[cl=="a"]
pmetsfl1.cla <- rbind(pmetsfl1.cla, pmets1)


pmetsfl1.clb <- pmets[cl=="b"]
pmetsfl1.clb <- rbind(pmetsfl1.clb, pmets1)


pmetsfl1.clc <- pmets[cl=="c"]
pmetsfl1.clc <- rbind(pmetsfl1.clc, pmets1)






func_class.ssets <- function(pmets, cls)
{
  #tabulate per plot the number of pcs belonging to each class
  tbl <- with(pmets, table(id_placette, cl))
  
  if(cls=="a")
  {
    #pick all the plots that have atleast one pc belonging to class a
    pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,1]>0),]))]
    
    #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
    pmets.onlycl <- pmets[cl=="a"]
    
    #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
    pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,1]>0),]))]
    
    #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
    pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
  }
  else if(cls=="b")
  {
    #pick all the plots that have atleast one pc belonging to class a
    pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,2]>0),]))]
    
    #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
    pmets.onlycl <- pmets[cl=="b"]
    
    #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
    pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,2]>0),]))]
    
    #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
    pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
  }  
  else
  {
    #pick all the plots that have atleast one pc belonging to class a
    pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,3]>0),]))]
    
    #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
    pmets.onlycl <- pmets[cl=="c"]
    
    #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
    pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,3]>0),]))]
    
    #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
    pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
  }
  return(pmets.cl)
}

pmetsfl1.all <- plotmetsfl1

pmetsfl1.cla <- func_class.ssets(pmetsfl1.all, "a")
pmetsfl1.cla <- unique(pmetsfl1.cla[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl1.cla <- pmetsfl1.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl1.cla <- pmetsfl1.cla[, wt := wt/sum(wt), id_placette]

pmetsfl1.clb <- func_class.ssets(pmetsfl1.all, "b")
pmetsfl1.clb <- unique(pmetsfl1.clb[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl1.clb <- pmetsfl1.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl1.clb <- pmetsfl1.clb[, wt := wt/sum(wt), id_placette]

pmetsfl1.clc <- func_class.ssets(pmetsfl1.all, "c")
pmetsfl1.clc <- unique(pmetsfl1.clc[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl1.clc <- pmetsfl1.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl1.clc <- pmetsfl1.clc[, wt := wt/sum(wt), id_placette]



{
# func_dffilter <- function(df, cls)
# {
#   df1 <- df[id_placette %in% cps]
#   df1 <- df1[df1[, .I[cl==cls  | all(cl!=cls)], by = id_placette]$V1]
#   df2 <- df[!id_placette %in% cps]
#   dfn <- rbind(df1, df2)
#   return(dfn)
# }
# 
# 
# plotmetsfl1a <- func_dffilter(plotmetsfl1, "a")
# plotmetsfl1a <- unique(plotmetsfl1a[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1a <- plotmetsfl1a[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1a <- plotmetsfl1a[, wt := wt/sum(wt), id_placette]
# 
# 
# plotmetsfl1b <- func_dffilter(plotmetsfl1, "b")
# plotmetsfl1b <- unique(plotmetsfl1b[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1b <- plotmetsfl1b[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1b <- plotmetsfl1b[, wt := wt/sum(wt), id_placette]
# 
# plotmetsfl1c <- func_dffilter(plotmetsfl1, "c")
# plotmetsfl1c <- unique(plotmetsfl1c[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1c <- plotmetsfl1c[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1c <- plotmetsfl1c[, wt := wt/sum(wt), id_placette]
# 
# 
# plotmetsfl1a=plotmetsfl1[plotmetsfl1[, all(cl != 'a')| (cl == 'a' & .BY %in% cps)|!.BY %in% cps,
#                                      by = id_placette]$V1]
# 
# plotmetsfl1a <- unique(plotmetsfl1a[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1a <- plotmetsfl1a[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1a <- plotmetsfl1a[, wt := wt/sum(wt), id_placette]
# 
# plotmetsfl1b=plotmetsfl1[plotmetsfl1[, all(cl != 'b')| (cl == 'b' & .BY %in% cps)|!.BY %in% cps,
#                                      by = id_placette]$V1]
# plotmetsfl1b <- unique(plotmetsfl1b[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1b <- plotmetsfl1b[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1b <- plotmetsfl1b[, wt := wt/sum(wt), id_placette]
# 
# plotmetsfl1c=plotmetsfl1[plotmetsfl1[, all(cl != 'c')| (cl == 'c' & .BY %in% cps)|!.BY %in% cps,
#                                      by = id_placette]$V1]
# plotmetsfl1c <- unique(plotmetsfl1c[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1c <- plotmetsfl1c[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1c <- plotmetsfl1c[, wt := wt/sum(wt), id_placette]
}


###########################################################################################################
fd_smry <- fread("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",", drop = "id_placette")
colnames(fd_smry)[1] <- "id_placette"
fd_smry[, id_placette := as.factor(id_placette)]
setkey(fd_smry, "id_placette")
###########################################################################################################
##Generate samples######
########################
dbase <- pmetsfl1.cla
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
smplst.cla <- foreach(i = 1:5000, .packages=c("dplyr", "data.table", "caret", "sampling")) %dopar% {
  set.seed(i)
  dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
  return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
  # return(sset)
  }
stopCluster(clus)
registerDoSEQ()
smplst.cla <- matrix(unlist(smplst.cla), nrow = 30)
smplst.cla <- t(smplst.cla)
smplst.cla <- as.data.table(smplst.cla)
smplst.cla <- unique(smplst.cla)





time_log <- data.frame()

start <- Sys.time()
dbase <- pmetsfl1.cla
fds <- fd_smry
idx.lst <- smplst.cla

f1l <- log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(G_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


##Simulations
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
nme1a <- "ciron"
n <- 5000
ciron.cla <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  idx <- as.vector(unlist(idx.lst[i]))
  mets_for_model <- dbase[idx]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fds[mets_for_model]
  mets1 <- mets_for_model[,c("G_m2_ha","meanch", "varch", 'pflidr', "cvladlidr")]
  
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


{
  # ciron <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  #   idx <- as.vector(unlist(smplst[i]))
  #   mets_for_model <- dbase[idx]
  #   setkey(mets_for_model,"id_placette")
  #   mets_for_model <- fds[mets_for_model]
  #   
  #   func_extractmdlmets <- function(mdl)
  #   {
  #     ypred <- exp(mdl$pred$pred)
  #     yobs <- exp(mdl$pred$obs)
  #     n <- length(yobs)
  #     MPE <- (100/n)*sum((yobs-ypred)/yobs)
  #     RMSE <- sqrt(mean((yobs-ypred)^2))
  #     MAE <- mean(abs(ypred-yobs))
  #     return(list("R2"= mdl$results$Rsquared, "MPE" = MPE, "RMSE" = RMSE,"MAE" = MAE))
  #   }
  #   
  #   #########################################################
  #   model <- train(f1l,
  #                  data = mets_for_model,
  #                  method = "lm",
  #                  trControl = trainControl(method="LOOCV"))
  #   G_lidr_coeff <- list(unlist(model$finalm1$coefficients, recursive = F))
  #   G_lidr_mets <- func_extractmdlmets(model)
  #   
  #   
  #   model <- train(f1v,
  #                  data = mets_for_model,
  #                  method = "lm",
  #                  trControl = trainControl(method="LOOCV"))
  #   G_vox_coeff <- model$finalModel$coefficients
  #   G_vox_mets <- func_extractmdlmets(model)
  #   ##########################################################
  #   
  #   ##########################################################
  #   model <- train(f2l,
  #                  data = mets_for_model,
  #                  method = "lm",
  #                  trControl = trainControl(method="LOOCV"))
  #   Vtot_lidr_coeff <- model$finalModel$coefficients
  #   Vtot_lidr_mets <- func_extractmdlmets(model)
  #   
  #   model <- train(f2v,
  #                  data = mets_for_model,
  #                  method = "lm",
  #                  trControl = trainControl(method="LOOCV"))
  #   Vtot_vox_coeff <- model$finalModel$coefficients
  #   Vtot_vox_mets <- func_extractmdlmets(model)
  #   ##########################################################
  #   
  #   ###########################################################
  #   model <- train(f3l,
  #                  data = mets_for_model,
  #                  method = "lm",
  #                  trControl = trainControl(method="LOOCV"))
  #   Vtig_lidr_coeff <- model$finalModel$coefficients
  #   Vtig_lidr_mets <- func_extractmdlmets(model)
  #   
  #   model <- train(f3v,
  #                  data = mets_for_model,
  #                  method = "lm",
  #                  trControl = trainControl(method="LOOCV"))
  #   Vtig_vox_coeff <- model$finalModel$coefficients
  #   Vtig_vox_mets <- func_extractmdlmets(model)
  #   ############################################################
  #   
  #   
  #   
  #   
  #   x <- (list("G_lidr_mets" = G_lidr_mets,
  #              "G_vox_mets" = G_vox_mets,
  #              "Vtot_lidr_mets" = Vtot_lidr_mets,
  #              "Vtot_vox_mets" = Vtot_vox_mets,
  #              "Vtig_lidr_mets" = Vtig_lidr_mets,
  #              "Vtig_vox_mets" = Vtig_vox_mets,
  #              "G_lidr_coeff" = G_lidr_coeff,
  #              "G_vox_coeff" = G_vox_coeff,
  #              "Vtot_lidr_coeff" = Vtot_lidr_coeff,
  #              "Vtot_vox_coeff" = Vtot_vox_coeff,
  #              "Vtig_lidr_coeff" = Vtig_lidr_coeff,
  #              "Vtig_vox_coeff" = Vtig_vox_coeff
  #   ))
  #   return(x)
  # }
  
  
}


func_mdlmets <- function(obs, pred, for_attr, mettype)
{
  yobs <- exp(obs)
  ypred <- exp(pred)
  see <- sqrt(sum((obs-pred)^2)/(length(obs)-4))
  cf <- exp((see^2)/2)
  ypred <- ypred*cf
  R2 <- cor(pred, obs)^2
  aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
  MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
  RMSE <- sqrt(mean((yobs-ypred)^2))
  MAE <- mean(abs(ypred-yobs))
  return(list( "Forest_attr"=for_attr, "Metrics"=mettype, "R2"=aR2, "MPE"=MPE, "RMSE"=RMSE,"MAE"=MAE))
}


ciron.mdlmets <- lapply(ciron.all, function(x)
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
})
ciron.mdlmets <- rbindlist(ciron.mdlmets)
ciron.mdlmets1 <- melt(ciron.mdlmets, measure.vars = c("R2", "MPE", "RMSE", "MAE"))
ggplot(data=ciron.mdlmets1[Forest_attr=="Total volume"], aes(y=value, fill=Metrics))+
  geom_boxplot()+facet_grid(variable~Forest_attr, scales = "free")










mdlmets <- rbindlist(lapply(ciron, function(x){
  gl <- rbindlist(x[1], idcol = "id")
  gv <- rbindlist(x[2], idcol = "id")
  g <- rbind(gl,gv)
  vtotl <- rbindlist(x[3], idcol = "id")
  vtotv <- rbindlist(x[4], idcol = "id")
  vtot <- rbind(vtotl,vtotv)
  vtigl <- rbindlist(x[5], idcol = "id")
  vtigv <- rbindlist(x[6], idcol = "id")
  vtig <- rbind(vtigl,vtigv)
  df <- rbind(g, vtot, vtig)
  return(df)
}))

mdlmets <- mdlmets[, c("for_attr","type", "mets"):=tstrsplit(id,"_",fixed=T),]
mdlmets <- melt(mdlmets[,-c("id", "mets")], id.vars = c("for_attr", "type"))


ggplot(data=mdlmets[for_attr=="Vtot"], aes(y=value, fill=type))+
  geom_boxplot()+
  labs(title = "")+
  theme_base(base_size = 20)+
  facet_grid(variable~for_attr, scales = "free")+
  geom_hline(data = allmdlmets, aes(yintercept=value, linetype=type, colour=mets), size=1.5)+
  scale_fill_pander()
























nme1b <- paste0("mdlmets_gfl1", n1, n2)
assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
                                      func_extractmdlmets, 
                                      type="One flight line", 
                                      exp="allclasses", 
                                      mets="old"))))
print("BA done")


##Simulations for total volume
set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
nme2a <- paste0("mdls_vtotfl1", n1, n2)
assign(nme2a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fds[mets_for_model]
  
  
  model <- train(f2,
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  ypred <- exp(model$pred$pred)
  yobs <- exp(model$pred$obs)
  n <- length(yobs)
  MPE <- (100/n)*sum((yobs-ypred)/yobs)
  RMSE <- sqrt(mean((yobs-ypred)^2))
  MAE <- mean(abs(ypred-yobs))
  
  
  x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                "RMSE" = RMSE,
                                "MAE" = MAE,
                                "MPE" = MPE),
             "coeffs" = list(model$finalModel$coefficients),
             "index" = which(dbase$meanang %in% mets_for_model$meanang)))
  return(x)
})
stopCluster(clus)
registerDoSEQ()
nme2b <- paste0("mdlmets_vtotfl1", n1, n2)
assign(nme2b, unique(rbindlist(lapply(get(nme2a), 
                                      func_extractmdlmets, 
                                      type="One flight line", 
                                      exp="allclasses", 
                                      mets="old"))))
print("Vtot done")

##Simulations for stem volume
set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
nme3a <- paste0("mdls_vtigfl1", n1, n2)
assign(nme3a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fds[mets_for_model]
  
  
  model <- train(f3,
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  ypred <- exp(model$pred$pred)
  yobs <- exp(model$pred$obs)
  n <- length(yobs)
  MPE <- (100/n)*sum((yobs-ypred)/yobs)
  RMSE <- sqrt(mean((yobs-ypred)^2))
  MAE <- mean(abs(ypred-yobs))
  
  
  x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                "RMSE" = RMSE,
                                "MAE" = MAE,
                                "MPE" = MPE),
             "coeffs" = list(model$finalModel$coefficients),
             "index" = which(dbase$meanang %in% mets_for_model$meanang)))
  return(x)
})
stopCluster(clus)
registerDoSEQ()
nme3b <- paste0("mdlmets_vtigfl1", n1, n2)
assign(nme3b, unique(rbindlist(lapply(get(nme3a), 
                                      func_extractmdlmets, 
                                      type="One flight line", 
                                      exp="allclasses", 
                                      mets="old"))))
print("Vtig done")
stop <- Sys.time()





if(type=="lidr")
{
  f1 <- log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f2 <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f3 <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  start <- Sys.time()
  ##Simulations for basal area
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme1a <- paste0("mdls_gfl1", n1, n2)
  assign(nme1a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f1,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme1b <- paste0("mdlmets_gfl1", n1, n2)
  assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("BA done")
  
  
  ##Simulations for total volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme2a <- paste0("mdls_vtotfl1", n1, n2)
  assign(nme2a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f2,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme2b <- paste0("mdlmets_vtotfl1", n1, n2)
  assign(nme2b, unique(rbindlist(lapply(get(nme2a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("Vtot done")
  
  ##Simulations for stem volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme3a <- paste0("mdls_vtigfl1", n1, n2)
  assign(nme3a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f3,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme3b <- paste0("mdlmets_vtigfl1", n1, n2)
  assign(nme3b, unique(rbindlist(lapply(get(nme3a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("Vtig done")
  stop <- Sys.time()
}else
{
  f1 <- log(G_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f2 <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f3 <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  start <- Sys.time()
  ##Simulations for basal area
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme1a <- paste0("mdls_gfl1vx", n1, n2)
  assign(nme1a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f1,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme1b <- paste0("mdlmets_gfl1vx", n1, n2)
  assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="vox"))))
  print("BA done")
  
  
  ##Simulations for total volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme2a <- paste0("mdls_vtotfl1vx", n1, n2)
  assign(nme2a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f2,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme2b <- paste0("mdlmets_vtotfl1vx", n1, n2)
  assign(nme2b, unique(rbindlist(lapply(get(nme2a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="vox"))))
  print("Vtot done")
  
  
  ##Simulations for stem volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme3a <- paste0("mdls_vtigfl1vx", n1, n2)
  assign(nme3a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f3,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme3b <- paste0("mdlmets_vtigfl1vx", n1, n2)
  assign(nme3b, unique(rbindlist(lapply(get(nme3a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="vox"))))
  print("Vtig done")
  stop <- Sys.time()
}



func_extractmdlmets <- function(model, type, exp, mets)
{
  # ind <- x[["index"]]
  # dbase <- plotmetsfl1[ind]
  # dbase <- dbase[!id_placette%in%c(21, 22, 26)]
  # ent <- func_entropy(dbase$cl)
  ypred <- exp(model$pred$pred)
  yobs <- exp(model$pred$obs)
  n <- length(yobs)
  MPE <- (100/n)*sum((yobs-ypred)/yobs)
  RMSE <- sqrt(mean((yobs-ypred)^2))
  MAE <- mean(abs(ypred-yobs))
  R2 <- model$results$Rsquared
  # names(tbl[which(tbl==max(tbl))])
  return(list("R2"=R2, "MPE"=MPE, "RMSE"=rmse, "MAE"=mae, "type"=type, "exp"=exp, "mets"=mets))
}


mdlmets_gfl1 <- unique(rbindlist(lapply(mdls_vtigfl1, 
                                            func_extractmdlmets, 
                                            type="One flight line", 
                                            exp="all classes 1 FL", 
                                            mets="old")))



mdlmets_g <- rbind(mdlmets_gfl1, 
                   mdlmets_gfl1pfcvlvx,
                   mdlmets_gfl1a,
                   mdlmets_gfl1apfcvlvx,
                   mdlmets_gfl1b,
                   mdlmets_gfl1bpfcvlvx,
                   mdlmets_gfl1c,
                   mdlmets_gfl1cpfcvlvx)
mdlmets_g <- melt(mdlmets_g, id.vars = c("type", "exp", "mets"))
mdlmets_g$for_attr <- rep("Basal area", nrow(mdlmets_g))

mdlmets_vtot <- rbind(mdlmets_vtotfl1, 
                      mdlmets_vtotfl1pfcvlvx,
                      mdlmets_vtotfl1a,
                      mdlmets_vtotfl1apfcvlvx,
                      mdlmets_vtotfl1b,
                      mdlmets_vtotfl1bpfcvlvx,
                      mdlmets_vtotfl1c,
                      mdlmets_vtotfl1cpfcvlvx)
mdlmets_vtot <- melt(mdlmets_vtot, id.vars = c("type", "exp", "mets"))
mdlmets_vtot$for_attr <- rep("Total volume", nrow(mdlmets_vtot))

mdlmets_vtig <- rbind(mdlmets_vtigfl1, 
                      mdlmets_vtigfl1pfcvlvx,
                      mdlmets_vtigfl1a,
                      mdlmets_vtigfl1apfcvlvx,
                      mdlmets_vtigfl1b,
                      mdlmets_vtigfl1bpfcvlvx,
                      mdlmets_vtigfl1c,
                      mdlmets_vtigfl1cpfcvlvx)
mdlmets_vtig <- melt(mdlmets_vtig, id.vars = c("type", "exp", "mets"))
mdlmets_vtig$for_attr <- rep("Stem volume", nrow(mdlmets_vtig))


mdlmets <- rbind(mdlmets_g, mdlmets_vtot, mdlmets_vtig)


ggplot(data=mdlmets, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "")+
  theme_base(base_size = 20)+
  facet_grid(variable~for_attr, scales = "free")+
  geom_hline(data = allmdlmets, aes(yintercept=value, linetype=type, colour=mets), size=1.5)+
  scale_fill_pander()
  
  







start <- Sys.time()
type <- "vox"
dbase <- plotmetsfl1c
fds <- fd_smry
n1 <- "ccc"
n2 <- "xxx"
n <- 1000
func_extractmdlmets <- function(x, type, exp, mets)
{
  ind <- x[["index"]]
  # dbase <- plotmetsfl1[ind]
  # dbase <- dbase[!id_placette%in%c(21, 22, 26)]
  # ent <- func_entropy(dbase$cl)
  r2 <- x[["modelmets"]]$R2 
  rmse <- x[["modelmets"]]$RMSE
  mae <- x[["modelmets"]]$MAE
  mpe <- x[["modelmets"]]$MPE
  type <- type
  exp <- exp
  mets <- mets
  # names(tbl[which(tbl==max(tbl))])
  return(list("r2"=r2, "rmse"=rmse, "mae"=mae, "mpe"=mpe, "type"=type, "exp"=exp, "mets"=mets))
}
if(type=="lidr")
{
  f1 <- log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f2 <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f3 <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  start <- Sys.time()
  ##Simulations for basal area
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme1a <- paste0("mdls_gfl1", n1, n2)
  assign(nme1a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f1,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme1b <- paste0("mdlmets_gfl1", n1, n2)
  assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("BA done")
  
  
  ##Simulations for total volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme2a <- paste0("mdls_vtotfl1", n1, n2)
  assign(nme2a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f2,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme2b <- paste0("mdlmets_vtotfl1", n1, n2)
  assign(nme2b, unique(rbindlist(lapply(get(nme2a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("Vtot done")
  
  ##Simulations for stem volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme3a <- paste0("mdls_vtigfl1", n1, n2)
  assign(nme3a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f3,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme3b <- paste0("mdlmets_vtigfl1", n1, n2)
  assign(nme3b, unique(rbindlist(lapply(get(nme3a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("Vtig done")
  stop <- Sys.time()
}else
{
  f1 <- log(G_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f2 <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f3 <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  start <- Sys.time()
  ##Simulations for basal area
  set.seed(123, kind = "L'Ecuyer-CMRG")
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme1a <- paste0("mdls_gfl1vx", n1, n2)
  assign(nme1a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    return(mets_for_model)
    
    # model <- train(f1,
    #                data = mets_for_model,
    #                method = "lm",
    #                trControl = trainControl(method="LOOCV"))
    # 
    # ypred <- exp(model$pred$pred)
    # yobs <- exp(model$pred$obs)
    # n <- length(yobs)
    # MPE <- (100/n)*sum((yobs-ypred)/yobs)
    # RMSE <- sqrt(mean((yobs-ypred)^2))
    # MAE <- mean(abs(ypred-yobs))
    # 
    # 
    # x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
    #                               "RMSE" = RMSE,
    #                               "MAE" = MAE,
    #                               "MPE" = MPE),
    #            "coeffs" = list(model$finalModel$coefficients),
    #            "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    # return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  # nme1b <- paste0("mdlmets_gfl1vx", n1, n2)
  # assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
  #                                       func_extractmdlmets, 
  #                                       type="One flight line", 
  #                                       exp="allclasses", 
  #                                       mets="vox"))))
  print("BA done")
  
  
  ##Simulations for total volume
  set.seed(123, kind = "L'Ecuyer-CMRG")
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme2a <- paste0("mdls_vtotfl1vx", n1, n2)
  assign(nme2a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    return(mets_for_model)
    
    
  #   model <- train(f2,
  #                  data = mets_for_model,
  #                  method = "lm",
  #                  trControl = trainControl(method="LOOCV"))
  #   
  #   ypred <- exp(model$pred$pred)
  #   yobs <- exp(model$pred$obs)
  #   n <- length(yobs)
  #   MPE <- (100/n)*sum((yobs-ypred)/yobs)
  #   RMSE <- sqrt(mean((yobs-ypred)^2))
  #   MAE <- mean(abs(ypred-yobs))
  #   
  #   
  #   x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
  #                                 "RMSE" = RMSE,
  #                                 "MAE" = MAE,
  #                                 "MPE" = MPE),
  #              "coeffs" = list(model$finalModel$coefficients),
  #              "index" = which(dbase$meanang %in% mets_for_model$meanang)))
  #   return(x)
   })
   stopCluster(clus)
   registerDoSEQ()
   # nme2b <- paste0("mdlmets_vtotfl1vx", n1, n2)
   # assign(nme2b, unique(rbindlist(lapply(get(nme2a),
   #                                       func_extractmdlmets,
   #                                       type="One flight line",
   #                                       exp="allclasses",
   #                                       mets="vox"))))
  print("Vtot done")
  
  
  ##Simulations for stem volume
  set.seed(123, kind = "L'Ecuyer-CMRG")
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme3a <- paste0("mdls_vtigfl1vx", n1, n2)
  assign(nme3a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    return(mets_for_model)
    
    
  #   model <- train(f3,
  #                  data = mets_for_model,
  #                  method = "lm",
  #                  trControl = trainControl(method="LOOCV"))
  #   
  #   ypred <- exp(model$pred$pred)
  #   yobs <- exp(model$pred$obs)
  #   n <- length(yobs)
  #   MPE <- (100/n)*sum((yobs-ypred)/yobs)
  #   RMSE <- sqrt(mean((yobs-ypred)^2))
  #   MAE <- mean(abs(ypred-yobs))
  #   
  #   
  #   x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
  #                                 "RMSE" = RMSE,
  #                                 "MAE" = MAE,
  #                                 "MPE" = MPE),
  #              "coeffs" = list(model$finalModel$coefficients),
  #              "index" = which(dbase$meanang %in% mets_for_model$meanang)))
  #   return(x)
    })
  stopCluster(clus)
  registerDoSEQ()
  # nme3b <- paste0("mdlmets_vtigfl1vx", n1, n2)
  # assign(nme3b, unique(rbindlist(lapply(get(nme3a),
  #                                       func_extractmdlmets,
  #                                       type="One flight line",
  #                                       exp="allclasses",
  #                                       mets="vox"))))
  print("Vtig done")
  stop <- Sys.time()
}


mean(unlist(lapply(mdls_gfl1vxcccxxx, function(x) return(mean(x$meanch)))))


x <- mdls_gfl1vxcccxxx[[1]]
