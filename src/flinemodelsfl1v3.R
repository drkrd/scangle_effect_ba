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
plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/flightlines_1/")

####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS
opt_independent_files(plots) <- TRUE 


height=5
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
plotmetsfl1 <- catalog_apply(plots, func_computeall)
plotmetsfl1 <- rbindlist(plotmetsfl1)
plotmetsfl1 <- plotmetsfl1[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                       ifelse(meanang>=10&meanang<20, "b",
                                              ifelse(meanang>=20&meanang<30, "c",
                                                     ifelse(meanang>=30&meanang<40,"d","e"))))]
plotmetsfl1 <- plotmetsfl1[cl!="d" & cl!="e"]
plotmetsfl1 <- plotmetsfl1[!meanang %in% c(37.24, 7.16, 4.60, 16.52)]



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
allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/2_Ciron/voxelisation/Results/May/wo_interpolation/voxfiles/flightlines_1/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","", basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","", basename(file_path_sans_ext(V1)))]
allvoxfiles$meanang <- round(as.numeric(allvoxfiles$meanang), 2)
#################################################################
voxall <- rbindlist(apply(allvoxfiles, 1, 
                          func_normvox2, 
                          pth="D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/", 
                          ht=height))



voxall$meanang <- as.numeric(voxall$meanang) 
voxall <- voxall[,cl:=ifelse(meanang>=0&meanang<10, "a",
                             ifelse(meanang>=10&meanang<20, "b",
                                    ifelse(meanang>=20&meanang<30, "c",
                                           ifelse(meanang>=30&meanang<40,"d","e"))))]


voxall <- voxall[cl!="d" & cl!="e"]
voxall <- voxall[!meanang %in% c(37.24, 7.16, 4.60, 16.52)]
##################################################################
pfcvladvox <- voxall[, .(cvladvox=cv(PADmean, na.rm = TRUE), 
                         sdvfp=sqrt(sum(PADmean*(k1-(sum(k1*PADmean)/sum(PADmean)))^2)/(sum(PADmean)*(length(PADmean[which(PADmean!=0)])-1)/length(PADmean[which(PADmean!=0)]))),
                         pfsumprof=exp(-0.5*sum(PADmean, na.rm = TRUE))), by=.(id_placette, meanang)]

setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(plotmetsfl1, c("id_placette", "meanang"))
plotmetsfl1 <- plotmetsfl1[pfcvladvox]
#############################################################################################################################
plotmetsfl1 <- unique(plotmetsfl1[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1 <- plotmetsfl1[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1 <- plotmetsfl1[, wt := wt/sum(wt), id_placette]
plotmetsfl1$meanang <- as.numeric(plotmetsfl1$meanang)
plotmetsfl1$rid <- seq(1:nrow(plotmetsfl1)) 
##############################################################################################################################

{
# plotmetsfl1a <- plotmetsfl1[plotmetsfl1[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
# plotmetsfl1a <- unique(plotmetsfl1a[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1a <- plotmetsfl1a[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1a <- plotmetsfl1a[, wt := wt/sum(wt), id_placette]
# 
# 
# plotmetsfl1b <- plotmetsfl1[plotmetsfl1[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
# plotmetsfl1b <- unique(plotmetsfl1b[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1b <- plotmetsfl1b[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1b <- plotmetsfl1b[, wt := wt/sum(wt), id_placette]
# 
# 
# plotmetsfl1c <- plotmetsfl1[plotmetsfl1[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]
# plotmetsfl1c <- unique(plotmetsfl1c[, prob := prop.table(table(cl))[cl], id_placette][])
# plotmetsfl1c <- plotmetsfl1c[, wt := (1/prob)/(1/sum(prob)), id_placette]
# plotmetsfl1c <- plotmetsfl1c[, wt := wt/sum(wt), id_placette]
}

#################################################################################################################################
pmetsfl1.all <- plotmetsfl1

#tabulate per plot the number of pcs belonging to each class
tbl <- with(plotmetsfl1, table(id_placette, cl))

#pick all the plots that have atleast one pc belonging to class a, class b and class c each
tbl2 <- tbl[which(tbl[,1]>0 & tbl[,2]>0 & tbl[,3]>0),]

#subset those pcs and their metrics which satisfy the above condition
pmets.withallcls <- pmetsfl1.all[id_placette %in% rownames(tbl2)]

#subset those pcs and their metrics which donot satisfy the previous condition
pmets.woallcls <- pmetsfl1.all[!id_placette %in% rownames(tbl2)]





#pick only those pcs and their metrics which belong to class a
pmetsfl1.cla <- pmets.withallcls[cl=="a"]

#combine with the remaining plots
pmetsfl1.cla <- rbind(pmetsfl1.cla, pmets.woallcls)

#compute inverse probabilities
pmetsfl1.cla <- unique(pmetsfl1.cla[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl1.cla <- pmetsfl1.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl1.cla <- pmetsfl1.cla[, wt := wt/sum(wt), id_placette]





#pick only those pcs and their metrics which belong to class b
pmetsfl1.clb <- pmets.withallcls[cl=="b"]

#combine with the remaining plots
pmetsfl1.clb <- rbind(pmetsfl1.clb, pmets.woallcls)

#compute inverse probabilities
pmetsfl1.clb <- unique(pmetsfl1.clb[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl1.clb <- pmetsfl1.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl1.clb <- pmetsfl1.clb[, wt := wt/sum(wt), id_placette]





#pick only those pcs and their metrics which belong to class b
pmetsfl1.clc <- pmets.withallcls[cl=="c"]

#combine with the remaining plots
pmetsfl1.clc <- rbind(pmetsfl1.clc, pmets.woallcls)

#compute inverse probabilities
pmetsfl1.clc <- unique(pmetsfl1.clc[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl1.clc <- pmetsfl1.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl1.clc <- pmetsfl1.clc[, wt := wt/sum(wt), id_placette]

{
# func_class.ssets <- function(pmets, cls)
# {
#   #tabulate per plot the number of pcs belonging to each class
#   tbl <- with(pmets, table(id_placette, cl))
#   
#   if(cls=="a")
#   {
#     #pick all the plots that have atleast one pc belonging to class a
#     pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,1]>0),]))]
#     
#     #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
#     pmets.onlycl <- pmets[cl=="a"]
#     
#     #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
#     pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,1]>0),]))]
#     
#     #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
#     pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
#   }
#   else if(cls=="b")
#   {
#     #pick all the plots that have atleast one pc belonging to class a
#     pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,2]>0),]))]
#     
#     #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
#     pmets.onlycl <- pmets[cl=="b"]
#     
#     #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
#     pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,2]>0),]))]
#     
#     #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
#     pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
#   }  
#   else
#   {
#     #pick all the plots that have atleast one pc belonging to class a
#     pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,3]>0),]))]
#     
#     #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
#     pmets.onlycl <- pmets[cl=="c"]
#     
#     #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
#     pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,3]>0),]))]
#     
#     #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
#     pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
#   }
#   return(pmets.cl)
# }
# 
# pmetsfl1.all <- plotmetsfl1
# 
# pmetsfl1.cla <- func_class.ssets(pmetsfl1.all, "a")
# pmetsfl1.cla <- unique(pmetsfl1.cla[, prob := prop.table(table(cl))[cl], id_placette][])
# pmetsfl1.cla <- pmetsfl1.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
# pmetsfl1.cla <- pmetsfl1.cla[, wt := wt/sum(wt), id_placette]
# 
# pmetsfl1.clb <- func_class.ssets(pmetsfl1.all, "b")
# pmetsfl1.clb <- unique(pmetsfl1.clb[, prob := prop.table(table(cl))[cl], id_placette][])
# pmetsfl1.clb <- pmetsfl1.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
# pmetsfl1.clb <- pmetsfl1.clb[, wt := wt/sum(wt), id_placette]
# 
# pmetsfl1.clc <- func_class.ssets(pmetsfl1.all, "c")
# pmetsfl1.clc <- unique(pmetsfl1.clc[, prob := prop.table(table(cl))[cl], id_placette][])
# pmetsfl1.clc <- pmetsfl1.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
# pmetsfl1.clc <- pmetsfl1.clc[, wt := wt/sum(wt), id_placette]
}


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
fd_smry$id_placette <- as.character(fd_smry$id_placette)
setkey(fd_smry, "id_placette")
fd_smry <- ciron_db
fd_smry <- fd_smry[!id_placetten%in%c(14,144)]
###########################################################################################################
##Generate samples######
########################
dbase <- pmetsfl1.all
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
smplst.all <- foreach(i = 1:8500, .packages=c("dplyr", "data.table", "caret", "sampling")) %dopar% {
  set.seed(i)
  dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
  return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
  }
stopCluster(clus)
registerDoSEQ()
smplst.all <- matrix(unlist(smplst.all), nrow = 29)
smplst.all <- unique(as.data.table(t(smplst.all)))

{
# a1 <- c(1,1,1,2,2,2,2,2,3,3,3,4,5,5,5,5,5,5,6,6,7,7,7,7)
# set.seed(9292)
# a2 <- rnorm(24)
# xdt <- as.data.table(cbind(a1, a2))
# idlist <- list()
# for(i in 1:500){
#   set.seed(i)
#   xdt2 <- xdt[, .SD[sample(.N, min(1,.N))], by = a1]
#   idlist[[i]] <- which(xdt$a2%in%xdt2$a2)
#   # return(sset)
# }
# idlist <- matrix(unlist(idlist), nrow = 7)
# idlist <- t(idlist)
# idlist <- as.data.table(idlist)
# idlist <- unique(idlist)
}

fd_smry_bkp <- fds


time_log <- data.frame()

start <- Sys.time()
dbase <- pmetsfl1.all
fds <- fd_smry_bkp[!id_placetten%in%c(14,144)]
idx.lst <- smplst.all

f1l <- log(g_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(g_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(vtot_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(vtot_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(vtige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(vtige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


##Simulations
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
nme1a <- "ciron"
n <- 5
cironfl1.clc <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
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
cironfl1.ref <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  idx <- as.vector(unlist(idx.lst[i]))
  mets_for_model <- dbase[idx]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fds[mets_for_model]
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
  mets1 <- mets_for_model[,c("g_m2_ha","meanch", "varch", 'pflidr', "cvladlidr")]
  m1 <- train(log(g_m2_ha)~log(meanch)+log(varch),
               data = mets1,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  m1.coeff <- m1$finalModel$coefficients
  m1.pred <- m1$pred$pred
  m1.obs <- m1$pred$obs
  
  
  mets2 <- mets_for_model[,c("g_m2_ha","meanch", "varch", 'pflidr', "cvladlidr")]
  m2 <- train(log(g_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
               data = mets2,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  m2.coeff <- m2$finalModel$coefficients
  m2.pred <- m2$pred$pred
  m2.obs <- m2$pred$obs
  ##########################################################
  
  ##########################################################
  mets3 <- mets_for_model[,c("vtot_m3_ha", "meanch", "varch", 'pflidr', "cvladlidr")]
  m3 <- train(log(vtot_m3_ha)~log(meanch)+log(varch),
               data = mets3,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  m3.coeff <- m3$finalModel$coefficients
  m3.pred <- m3$pred$pred
  m3.obs <- m3$pred$obs
  
  mets4 <- mets_for_model[,c("vtot_m3_ha", "meanch", "varch", 'pflidr', "cvladlidr")]
  m4 <- train(log(vtot_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
               data = mets4,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  m4.coeff <- m4$finalModel$coefficients
  m4.pred <- m4$pred$pred
  m4.obs <- m4$pred$obs
  ##########################################################
  
  ###########################################################
  mets5 <- mets_for_model[,c("vtige_m3_ha", "meanch", "varch", 'pflidr', "cvladlidr")]
  m5 <- train(log(vtige_m3_ha)~log(meanch)+log(varch),
              data = mets5,
              method = "lm",
              trControl = trainControl(method="LOOCV"))
  m5.coeff <- m5$finalModel$coefficients
  m5.pred <- m5$pred$pred
  m5.obs <- m5$pred$obs
  
  mets6 <- mets_for_model[,c("vtige_m3_ha", "meanch", "varch", 'pflidr', "cvladlidr")]
  m6 <- train(log(vtige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
              data = mets6,
              method = "lm",
              trControl = trainControl(method="LOOCV"))
  m6.coeff <- m6$finalModel$coefficients
  m6.pred <- m6$pred$pred
  m6.obs <- m6$pred$obs
  ############################################################
  ############################################################
  
  x <- (list("m1.pred" = m1.pred,
             "m1.obs" = m1.obs,
             "m2.pred" = m2.pred,
             "m2.obs" = m2.obs,
             "m3.pred" = m3.pred,
             "m3.obs" = m3.obs,
             "m4.pred" = m4.pred,
             "m4.obs" = m4.obs,
             "m5.pred" = m5.pred,
             "m5.obs" = m5.obs,
             "m6.pred" = m6.pred,
             "m6.obs" = m6.obs
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


ciron.mdlmets.all <- melt(rbindlist(lapply(cironfl1.all, function(x)
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
ciron.mdlmets.cla <- melt(rbindlist(lapply(cironfl1.cla, function(x)
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
ciron.mdlmets.clb <- melt(rbindlist(lapply(cironfl1.clb, function(x)
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
ciron.mdlmets.clc <- melt(rbindlist(lapply(cironfl1.clc, function(x)
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



ciron.mdlmets.ref <- melt(rbindlist(lapply(cironfl1.ref, function(x)
{
  m1.mdlmets <- func_mdlmets(x$m1.obs, x$m1.pred, "G1", "old")
  m2.mdlmets <- func_mdlmets(x$m2.obs, x$m2.pred, "G2", "old")
  m3.mdlmets <- func_mdlmets(x$m5.obs, x$m5.pred, "Vst1", "old")
  m4.mdlmets <- func_mdlmets(x$m6.obs, x$m6.pred, "Vst2", "old")
  m5.mdlmets <- func_mdlmets(x$m3.obs, x$m3.pred, "Vt1", "old")
  m6.mdlmets <- func_mdlmets(x$m4.obs, x$m4.pred, "Vt2", "old")
  

  y <- list(m1.mdlmets,
            m2.mdlmets,
            m3.mdlmets,
            m4.mdlmets,
            m5.mdlmets,
            m6.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))


ciron.mdlmets.all <- cbind(ciron.mdlmets.all, "exp"=rep("fl1", nrow(ciron.mdlmets.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
ciron.mdlmets.cla <- cbind(ciron.mdlmets.cla, "exp"=rep("A", nrow(ciron.mdlmets.cla)), "id"=rep(rep(1:5000, 1, each=6), 4))
ciron.mdlmets.clb <- cbind(ciron.mdlmets.clb, "exp"=rep("B", nrow(ciron.mdlmets.clb)), "id"=rep(rep(1:5000, 1, each=6), 4))
ciron.mdlmets.clc <- cbind(ciron.mdlmets.clc, "exp"=rep("C", nrow(ciron.mdlmets.clc)), "id"=rep(rep(1:5000, 1, each=6), 4))


ciron.mdlmets.ref <- cbind(ciron.mdlmets.ref, "exp"=rep("all", nrow(ciron.mdlmets.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
ciron.mdlmets.ref <- as.data.table(rbind(ciron.mdlmets.ref))
ciron.mdlmets.ref <- melt(ciron.mdlmets.ref)


cironfl1.mdlmets <- as.data.table(rbind(ciron.mdlmets.all, 
                                        ciron.mdlmets.cla, 
                                        ciron.mdlmets.clb,
                                        ciron.mdlmets.clc))
     
cironfl1.mdlmets <- cbind(cironfl1.mdlmets, "fl"=rep("fl1", nrow(cironfl1.mdlmets)))
cironfl2.mdlmets <- cbind(cironfl2.mdlmets, "fl"=rep("fl2", nrow(cironfl2.mdlmets)))
cironfl3.mdlmets.all <- cbind(cironfl3.mdlmets.all, "fl"=rep("fl3", nrow(cironfl3.mdlmets.all)))


ciron.mdlmets <- rbind(cironfl1.mdlmets, cironfl2.mdlmets, cironfl3.mdlmets.all)
ciron.mdlmets$exp <- factor(ciron.mdlmets$exp, levels = c("fl1", "A","B","C",
                                                          "fl2","AC","AB","BC","fl3"))

ciron.mdlmets$forest_type <- rep("Riparian", nrow(ciron.mdlmets))

            
bp_g <- ggplot(data=ciron.mdlmets[variable!="RMSE" & Forest_attr=="Basal area"], aes(fill=Metrics, y=value, x=exp))+
  geom_boxplot()+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_minimal()+
  labs(x="Experiments", y="")+
  scale_fill_grey(start = 1, end = 0.5)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")

bp_vst <- ggplot(data=ciron.mdlmets[variable!="RMSE" & Forest_attr=="Stem volume"], aes(fill=Metrics, y=value, x=exp))+
  geom_boxplot()+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_minimal()+
  labs(x="Experiments", y="")+
  scale_fill_grey(start = 1, end = 0.5)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")
bp_vt <- ggplot(data=ciron.mdlmets[variable!="RMSE" & Forest_attr=="Total volume"], aes(fill=Metrics, y=value, x=exp))+
  geom_boxplot()+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_base()+
  labs(x="Experiments", y="")+
  scale_fill_grey(start = 1, end = 0.5)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")

ggsave(bp_g, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/bp_g.png", 
       width=16*1.25, height=8*1.25, units="cm", dpi=320) 



ggplot(data=ciron.mdlmets[variable!="RMSE" & 
                             Forest_attr=="Basal area" &
                             exp %in% c("fl1", "fl2", "fl3") &
                             Metrics=="ref"], aes(fill=Metrics, y=value, x=exp))+
  geom_boxplot()+
  facet_wrap(variable~forest_type, scales = "free", ncol=1)+
  theme_minimal()+
  labs(x="Experiments", y="")+
  scale_fill_grey(start = 1, end = 0.5)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")







x <- cironfl1.mdlmets[Forest_attr=="Basal area" & variable!="RMSE" ]




gl.coeffs <- NULL
for(i in 1:length(ciron.all))
{
  gl.coeffs <- rbind(gl.coeffs, ciron.all[[i]]$G.l.coeff)
}


colnames(gl.coeffs) <- c("Intercept", "meanch", "varch", "pf", "cvlad")
gl.coeffs <- as.data.table(gl.coeffs)
gl.coeffs1 <- melt(gl.coeffs)
gv.coeffs1 <- gv.coeffs[, .(mean(Intercept),
                            mean(meanch),
                            mean(varch),
                            mean(pf),
                            mean(cvlad))]


ggplot(data=gl.coeffs1[variable!="Intercept"], aes(y=exp(value), colour=variable))+geom_boxplot()

mets_all2 <- mets_all2[, `:=`( g.v = gv.coeffs1$V1 + 
                                 gv.coeffs1$V2*log(meanch) + 
                                 gv.coeffs1$V3*log(varch) + 
                                 gv.coeffs1$V4*log(pflidr) + 
                                 gv.coeffs1$V5*log(cvladlidr))]








