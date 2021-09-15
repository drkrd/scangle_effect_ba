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
bauges_db <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db15jul.csv", sep = ",")
colnames(bauges_db)[2] <- "id_placette"
height <- 5

bauges_db <- bauges_db[, newstratum := ifelse(comp_R_G>75 & comp_F_G<25 , "Coniferes",
                                              ifelse(comp_F_G>75 & comp_R_G<25, "Feuillus", "Mixte"))]
bauges_db <- bauges_db[which(id_placette %in% plot)]
setkey(bauges_db, "id_placette")


plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/flightlines_2/")
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
plotmetsfl2 <- catalog_apply(plots, func_computeall)
pmetsfl2 <- plotmetsfl2
pmetsfl2 <- rbindlist(pmetsfl2)
pmetsfl2 <- pmetsfl2[, c("fl1","fl2"):=tstrsplit(meanang,"_",fixed=T),]
pmetsfl2$fl1 <- as.numeric(pmetsfl2$fl1)
pmetsfl2$fl2 <- as.numeric(pmetsfl2$fl2)
pmetsfl2 <- pmetsfl2[, fl2:= ifelse(is.na(fl2), fl1, fl2),]

pmetsfl2 <- pmetsfl2[, cl1:=ifelse(fl1>=0&fl1<10, "a",
                                   ifelse(fl1>=10&fl1<20, "b",
                                          ifelse(fl1>=20&fl1<30, "c",
                                                 ifelse(fl1>=30&fl1<40,"d","e"))))]
pmetsfl2 <- pmetsfl2[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                   ifelse(fl2>=10&fl2<20, "b",
                                          ifelse(fl2>=20&fl2<30, "c",
                                                 ifelse(fl2>=30&fl2<40,"d","e"))))]

pmetsfl2 <- pmetsfl2[cl1 != "e" & cl2 != "e" & cl1 != "d" & cl2 != "d",]

pmetsfl2 <- pmetsfl2[!(id_placette=="96_IRSTEA" & fl1==30.0)]
pmetsfl2 <- pmetsfl2[!(id_placette=="96_IRSTEA" & fl2==30.0)]
pmetsfl2 <- pmetsfl2[!(id_placette=="283_74_ONF" & fl1==6.28)]
pmetsfl2 <- pmetsfl2[!(id_placette=="283_74_ONF" & fl2==6.28)]



pmetsfl2 <- pmetsfl2[, cl := paste0(sort(.SD), collapse = ""), .SDcols = c("cl1", "cl2"),  by = 1:nrow(pmetsfl2)]


##############################################################################################################################
allpcs <- list.files(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/"),
                     pattern = "*.las",
                     full.names = TRUE)

alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)
names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))
##################################################################################################################################
allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/flightlines_2/"), 
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



allvoxfiles <- allvoxfiles[cl1 != "e" & cl2 != "e" & cl1 != "d" & cl2 != "d",]
allvoxfiles <- allvoxfiles[!(id_placette=="96_IRSTEA_un" & fl1==30.0)]
allvoxfiles <- allvoxfiles[!(id_placette=="96_IRSTEA_un" & fl2==30.0)]
allvoxfiles <- allvoxfiles[!(id_placette=="283_74_ONF_un" & fl1==6.28)]
allvoxfiles <- allvoxfiles[!(id_placette=="283_74_ONF_un" & fl2==6.28)]

#################################################################

voxall <- rbindlist(apply(allvoxfiles, 1, 
                          func_normvox2, 
                          pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", 
                          ht=height))
setDT(voxall)

pfcvladvox <- voxall[, .(cvladvox=cv(PADmean, na.rm = TRUE), 
                         sdvfp=sqrt(sum(PADmean*(k1-(sum(k1*PADmean)/sum(PADmean)))^2)/(sum(PADmean)*(length(PADmean[which(PADmean!=0)])-1)/length(PADmean[which(PADmean!=0)]))),
                         pfsumprof=exp(-0.5*sum(PADmean, na.rm = TRUE))), by=.(id_placette, meanang)]

bauges_db <- bauges_db[which(id_placette %in% pmetsfl2$id_placette), 
                       c("id_placette", "G75", "volume_total", "volume_tige", "newstratum")]
bauges_db_con <- bauges_db[newstratum=="Coniferes"]
bauges_db_feu <- bauges_db[newstratum=="Feuillus"]
bauges_db_mix <- bauges_db[newstratum=="Mixte"]


setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(pmetsfl2, c("id_placette", "meanang"))
pmetsfl2 <- pmetsfl2[pfcvladvox]
pmetsfl2 <- pmetsfl2[which(id_placette %in% bauges_db$id_placette)]
pmetsfl2 <- pmetsfl2[!pflidr==0]


######################################
#############All################
#############################################################################################################################
{
  plotmetsfl2 <- unique(plotmetsfl2[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2 <- plotmetsfl2[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2 <- plotmetsfl2[, wt := wt/sum(wt), id_placette]
  plotmetsfl2$meanang <- as.factor(plotmetsfl2$meanang)
  plotmetsfl2 <- plotmetsfl2[pflidr!=0]
}


######################################
#############Coniferes################
#############################################################################################################################

{
  pmetsfl2.con.all <- pmetsfl2[which(id_placette %in% bauges_db_con$id_placette)]
  pmetsfl2.con.all <- unique(pmetsfl2.con.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl2.con.all <- pmetsfl2.con.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl2.con.all <- pmetsfl2.con.all[, wt := wt/sum(wt), id_placette]
  
  #tabulate per plot the number of pcs belonging to each class
  tbl <- with(pmetsfl1.con.all, table(id_placette, cl))
  
  #pick all the plots that have atleast one pc belonging to class a, class b and class c each
  tbl2 <- tbl[which(tbl[,1]>0 & tbl[,2]>0 & tbl[,3]>0),]
  
  #subset those pcs and their metrics which satisfy the above condition
  pmetsfl1.con.withallcls <- pmetsfl1.con.all[id_placette %in% rownames(tbl2)]
  
  #subset those pcs and their metrics which donot satisfy the previous condition
  pmetsfl1.con.woallcls <- pmetsfl1.con.all[!id_placette %in% rownames(tbl2)]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to clabss a
  pmetsfl1.con.clab <- pmets.withallcls[cl=="ab"]
  
  #combine with the remaining plots
  pmetsfl1.con.clab <- rbind(pmetsfl1.con.clab, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.con.clab <- unique(pmetsfl1.con.clab[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.con.clab <- pmetsfl1.con.clab[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.con.clab <- pmetsfl1.con.clab[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.con.clac <- pmets.withallcls[cl=="ac"]
  
  #combine with the remaining plots
  pmetsfl1.con.clac <- rbind(pmetsfl1.con.clac, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.con.clac <- unique(pmetsfl1.con.clac[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.con.clac <- pmetsfl1.con.clac[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.con.clac <- pmetsfl1.con.clac[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.con.clbc <- pmets.withallcls[cl=="bc"]
  
  #combine with the remaining plots
  pmetsfl1.con.clbc <- rbind(pmetsfl1.con.clbc, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.con.clbc <- unique(pmetsfl1.con.clbc[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.con.clbc <- pmetsfl1.con.clbc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.con.clbc <- pmetsfl1.con.clbc[, wt := wt/sum(wt), id_placette]
  
    #-------------------------------------------------------------------------------------------------------------------------#
  
}


######################################
#############Feuillus#################
###########################################################################################################################

{
  pmetsfl2.feu.all <- pmetsfl2[which(id_placette %in% bauges_db_con$id_placette)]
  pmetsfl2.feu.all <- unique(pmetsfl2.feu.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl2.feu.all <- pmetsfl2.feu.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl2.feu.all <- pmetsfl2.feu.all[, wt := wt/sum(wt), id_placette]
  
  #tabulate per plot the number of pcs belonging to each class
  tbl <- with(pmetsfl1.feu.all, table(id_placette, cl))
  
  #pick all the plots that have atleast one pc belonging to class a, class b and class c each
  tbl2 <- tbl[which(tbl[,1]>0 & tbl[,2]>0 & tbl[,3]>0),]
  
  #subset those pcs and their metrics which satisfy the above condition
  pmetsfl1.feu.withallcls <- pmetsfl1.feu.all[id_placette %in% rownames(tbl2)]
  
  #subset those pcs and their metrics which donot satisfy the previous condition
  pmetsfl1.feu.woallcls <- pmetsfl1.feu.all[!id_placette %in% rownames(tbl2)]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to clabss a
  pmetsfl1.feu.clab <- pmets.withallcls[cl=="ab"]
  
  #combine with the remaining plots
  pmetsfl1.feu.clab <- rbind(pmetsfl1.feu.clab, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.feu.clab <- unique(pmetsfl1.feu.clab[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.feu.clab <- pmetsfl1.feu.clab[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.feu.clab <- pmetsfl1.feu.clab[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.feu.clac <- pmets.withallcls[cl=="ac"]
  
  #combine with the remaining plots
  pmetsfl1.feu.clac <- rbind(pmetsfl1.feu.clac, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.feu.clac <- unique(pmetsfl1.feu.clac[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.feu.clac <- pmetsfl1.feu.clac[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.feu.clac <- pmetsfl1.feu.clac[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.feu.clbc <- pmets.withallcls[cl=="bc"]
  
  #combine with the remaining plots
  pmetsfl1.feu.clbc <- rbind(pmetsfl1.feu.clbc, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.feu.clbc <- unique(pmetsfl1.feu.clbc[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.feu.clbc <- pmetsfl1.feu.clbc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.feu.clbc <- pmetsfl1.feu.clbc[, wt := wt/sum(wt), id_placette]
  
  #-------------------------------------------------------------------------------------------------------------------------#
  
}

######################################
################Mixte#################
###########################################################################################################################
{
  pmetsfl2.mix.all <- pmetsfl2[which(id_placette %in% bauges_db_con$id_placette)]
  pmetsfl2.mix.all <- unique(pmetsfl2.mix.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl2.mix.all <- pmetsfl2.mix.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl2.mix.all <- pmetsfl2.mix.all[, wt := wt/sum(wt), id_placette]
  
  #tabulate per plot the number of pcs belonging to each class
  tbl <- with(pmetsfl1.mix.all, table(id_placette, cl))
  
  #pick all the plots that have atleast one pc belonging to class a, class b and class c each
  tbl2 <- tbl[which(tbl[,1]>0 & tbl[,2]>0 & tbl[,3]>0),]
  
  #subset those pcs and their metrics which satisfy the above condition
  pmetsfl1.mix.withallcls <- pmetsfl1.mix.all[id_placette %in% rownames(tbl2)]
  
  #subset those pcs and their metrics which donot satisfy the previous condition
  pmetsfl1.mix.woallcls <- pmetsfl1.mix.all[!id_placette %in% rownames(tbl2)]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to clabss a
  pmetsfl1.mix.clab <- pmets.withallcls[cl=="ab"]
  
  #combine with the remaining plots
  pmetsfl1.mix.clab <- rbind(pmetsfl1.mix.clab, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.mix.clab <- unique(pmetsfl1.mix.clab[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.mix.clab <- pmetsfl1.mix.clab[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.mix.clab <- pmetsfl1.mix.clab[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.mix.clac <- pmets.withallcls[cl=="ac"]
  
  #combine with the remaining plots
  pmetsfl1.mix.clac <- rbind(pmetsfl1.mix.clac, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.mix.clac <- unique(pmetsfl1.mix.clac[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.mix.clac <- pmetsfl1.mix.clac[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.mix.clac <- pmetsfl1.mix.clac[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.mix.clbc <- pmets.withallcls[cl=="bc"]
  
  #combine with the remaining plots
  pmetsfl1.mix.clbc <- rbind(pmetsfl1.mix.clbc, pmets.woallcls)
  
  #compute inverse probabilities
  pmetsfl1.mix.clbc <- unique(pmetsfl1.mix.clbc[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.mix.clbc <- pmetsfl1.mix.clbc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.mix.clbc <- pmetsfl1.mix.clbc[, wt := wt/sum(wt), id_placette]
  
  #-------------------------------------------------------------------------------------------------------------------------#
  
}




###########################################################################################################
##Generate samples######
########################
dbase <- pmetsfl2.mix.all
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
smplstfl2.mix.all <- foreach(i = 1:10000, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  set.seed(i)
  dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
  return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
  # return(sset)
}
stopCluster(clus)
registerDoSEQ()
smplstfl2.mix.all <- matrix(unlist(smplstfl2.mix.all), nrow = length(unique(dbase$id_placette)))
smplstfl2.mix.all <- unique(as.data.table(t(smplstfl2.mix.all)))
#########################################################################################################
time_log <- data.frame()
#########################################################################################################

start <- Sys.time()
dbase <- pmetsfl2.mix.all
fds <- bauges_db_mix
setkey(fds,"id_placette")
idx.lst <- smplstfl2.mix.all

f1l <- log(G75)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(G75)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


##Simulations
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
n <- 5000
baugesfl2.mix.all1 <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  idx <- as.vector(unlist(idx.lst[i]))
  mets_for_model <- dbase[idx]
  dbase$id_placette <- as.factor(dbase$id_placette)
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fds[mets_for_model]
  mets1 <- mets_for_model[,c("G75","meanch", "varch", 'pflidr', "cvladlidr")]
  
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
  
  
  mets2 <- mets_for_model[,c("G75","meanch", "varch", 'pfsumprof', "cvladvox")]
  m1v <- train(f1v,
               data = mets2,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  G.v.coeff <- m1v$finalModel$coefficients
  G.v.pred <- m1v$pred$pred
  G.v.obs <- m1v$pred$obs
  ##########################################################
  
  ##########################################################
  mets3 <- mets_for_model[,c("volume_total", "meanch", "varch", 'pflidr', "cvladlidr")]
  m2l <- train(f2l,
               data = mets3,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtot.l.coeff <- m2l$finalModel$coefficients
  vtot.l.pred <- m2l$pred$pred
  vtot.l.obs <- m2l$pred$obs
  
  mets4 <- mets_for_model[,c("volume_total", "meanch", "varch", 'pfsumprof', "cvladvox")]
  m2v <- train(f2v,
               data = mets4,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtot.v.coeff <- m2v$finalModel$coefficients
  vtot.v.pred <- m2v$pred$pred
  vtot.v.obs <- m2v$pred$obs
  ##########################################################
  
  ###########################################################
  mets5 <- mets_for_model[,c("volume_tige", "meanch", "varch", 'pflidr', "cvladlidr")]
  m3l <- train(f3l,
               data = mets5,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtig.l.coeff <- m3l$finalModel$coefficients
  vtig.l.pred <- m3l$pred$pred
  vtig.l.obs <- m3l$pred$obs
  
  mets6 <- mets_for_model[,c("volume_tige", "meanch", "varch", 'pfsumprof', "cvladvox")]
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

baugesfl2.mdlmets.con.all <- melt(rbindlist(lapply(baugesfl2.con.all1, function(x)
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
baugesfl2.mdlmets.feu.all <- melt(rbindlist(lapply(baugesfl2.feu.all1, function(x)
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
baugesfl2.mdlmets.mix.all <- melt(rbindlist(lapply(baugesfl2.mix.all1, function(x)
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




baugesfl2.mdlmets.con.all <- cbind(baugesfl2.mdlmets.con.all, "exp"=rep("all", nrow(baugesfl2.mdlmets.con.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl2.mdlmets.con.all <- cbind(baugesfl2.mdlmets.con.all, "fl"=rep("fl2", nrow(baugesfl2.mdlmets.con.all)))


baugesfl2.mdlmets.feu.all <- cbind(baugesfl2.mdlmets.feu.all, "exp"=rep("all", nrow(baugesfl2.mdlmets.feu.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl2.mdlmets.feu.all <- cbind(baugesfl2.mdlmets.feu.all, "fl"=rep("fl2", nrow(baugesfl2.mdlmets.feu.all)))


baugesfl2.mdlmets.mix.all <- cbind(baugesfl2.mdlmets.mix.all, "exp"=rep("all", nrow(baugesfl2.mdlmets.mix.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl2.mdlmets.mix.all <- cbind(baugesfl2.mdlmets.mix.all, "fl"=rep("fl2", nrow(baugesfl2.mdlmets.mix.all)))


ggplot(data=baugesfl2.mdlmets.feu.all[Forest_attr=="Stem volume"], aes( y=value, colour=Metrics))+
  geom_boxplot()+
  facet_grid(variable~., scales = "free")+
  theme_base()+
  scale_fill_grey()
