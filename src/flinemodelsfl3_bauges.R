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

plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/flightlines_3/")
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

pmetsfl3 <- pmetsfl3[!(id_placette=="96_IRSTEA" & fl1==30.0)]
pmetsfl3 <- pmetsfl3[!(id_placette=="96_IRSTEA" & fl2==30.0)]
pmetsfl3 <- pmetsfl3[!(id_placette=="96_IRSTEA" & fl3==30.0)]
pmetsfl3 <- pmetsfl3[!(id_placette=="283_74_ONF" & fl1==6.28)]
pmetsfl3 <- pmetsfl3[!(id_placette=="283_74_ONF" & fl2==6.28)]
pmetsfl3 <- pmetsfl3[!(id_placette=="283_74_ONF" & fl3==6.28)]




pmetsfl3 <- pmetsfl3[, cl := paste0(sort(.SD), collapse = ""), .SDcols = c("cl1", "cl2", "cl3"),  by = 1:nrow(pmetsfl3)]


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
#################################################################

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/flightlines_3/"), 
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
allvoxfiles <- allvoxfiles[!(id_placette=="96_IRSTEA_un" & fl1==30.0)]
allvoxfiles <- allvoxfiles[!(id_placette=="96_IRSTEA_un" & fl2==30.0)]
allvoxfiles <- allvoxfiles[!(id_placette=="96_IRSTEA_un" & fl3==30.0)]
allvoxfiles <- allvoxfiles[!(id_placette=="283_74_ONF_un" & fl1==6.28)]
allvoxfiles <- allvoxfiles[!(id_placette=="283_74_ONF_un" & fl2==6.28)]
allvoxfiles <- allvoxfiles[!(id_placette=="283_74_ONF_un" & fl3==6.28)]

#################################################################

voxall <- rbindlist(apply(allvoxfiles, 1, 
                          func_normvox2, 
                          pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", 
                          ht=height))
setDT(voxall)

pfcvladvox <- voxall[, .(cvladvox=cv(PADmean, na.rm = TRUE), 
                         sdvfp=sqrt(sum(PADmean*(k1-(sum(k1*PADmean)/sum(PADmean)))^2)/(sum(PADmean)*(length(PADmean[which(PADmean!=0)])-1)/length(PADmean[which(PADmean!=0)]))),
                         pfsumprof=exp(-0.5*sum(PADmean, na.rm = TRUE))), by=.(id_placette, meanang)]

bauges_db <- bauges_db[which(id_placette %in% pmetsfl3$id_placette), 
                       c("id_placette", "G75", "volume_total", "volume_tige", "newstratum")]
bauges_db_con <- bauges_db[newstratum=="Coniferes"]
bauges_db_feu <- bauges_db[newstratum=="Feuillus"]
bauges_db_mix <- bauges_db[newstratum=="Mixte"]


setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(pmetsfl3, c("id_placette", "meanang"))
pmetsfl3 <- pmetsfl3[pfcvladvox]
pmetsfl3 <- pmetsfl3[which(id_placette %in% bauges_db$id_placette)]
pmetsfl3 <- pmetsfl3[!pflidr==0]
#######################################################################

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
  pmetsfl3.con.all <- pmetsfl3[which(id_placette %in% bauges_db_con$id_placette)]
  pmetsfl3.con.all <- unique(pmetsfl3.con.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl3.con.all <- pmetsfl3.con.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl3.con.all <- pmetsfl3.con.all[, wt := wt/sum(wt), id_placette]
  pmetsfl3.con.all$meanang <- as.factor(pmetsfl3.con.all$meanang)
  pmetsfl3.con.all <- pmetsfl3.con.all[pflidr!=0]
  ##############################################################################################################################
  
  plotmetsfl2ab_con=plotmetsfl2_con[plotmetsfl2_con[, .I[cl=="ab"  | all(cl!="ab")], by = id_placette]$V1]
  plotmetsfl2ac_con=plotmetsfl2_con[plotmetsfl2_con[, .I[cl=="ac"  | all(cl!="ac")], by = id_placette]$V1]
  plotmetsfl2bc_con=plotmetsfl2_con[plotmetsfl2_con[, .I[cl=="bc"  | all(cl!="bc")], by = id_placette]$V1]
  
  l1_con <- unique(plotmetsfl2ab_con[cl=="ab"]$id_placette)
  l2_con <- unique(plotmetsfl2ac_con[cl=="ac"]$id_placette)
  l3_con <- unique(plotmetsfl2bc_con[cl=="bc"]$id_placette)
  
  cps <- intersect(intersect(l1_con, l2_con), l3_con)
  
  func_dffilter <- function(df, cls)
  {
    df1 <- df[id_placette %in% cps]
    df1 <- df1[df1[, .I[cl==cls  | all(cl!=cls)], by = id_placette]$V1]
    df2 <- df[!id_placette %in% cps]
    dfn <- rbind(df1, df2)
    return(dfn)
  }
  
  
  plotmetsfl2ab_con <- func_dffilter(plotmetsfl2_con, "ab")
  plotmetsfl2ab_con <- unique(plotmetsfl2ab_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2ab_con <- plotmetsfl2ab_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2ab_con <- plotmetsfl2ab_con[, wt := wt/sum(wt), id_placette]
  
  
  plotmetsfl2ac_con <- func_dffilter(plotmetsfl2_con, "ac")
  plotmetsfl2ac_con <- unique(plotmetsfl2ac_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2ac_con <- plotmetsfl2ac_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2ac_con <- plotmetsfl2ac_con[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl2bc_con <- func_dffilter(plotmetsfl2_con, "bc")
  plotmetsfl2bc_con <- unique(plotmetsfl2bc_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2bc_con <- plotmetsfl2bc_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2bc_con <- plotmetsfl2bc_con[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl2ab_con=plotmetsfl2_con[plotmetsfl2_con[, all(cl != 'ab')| (cl == 'ab' & .BY %in% cps)|!.BY %in% cps, 
                                                    by = id_placette]$V1]
  
  plotmetsfl2ab_con <- unique(plotmetsfl2ab_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2ab_con <- plotmetsfl2ab_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2ab_con <- plotmetsfl2ab_con[, wt := wt/sum(wt), id_placette]
  
  
  
  plotmetsfl2ac_con=plotmetsfl2_con[plotmetsfl2_con[, all(cl != 'ac')| (cl == 'ac' & .BY %in% cps)|!.BY %in% cps,
                                                    by = id_placette]$V1]
  
  plotmetsfl2ac_con <- unique(plotmetsfl2ac_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2ac_con <- plotmetsfl2ac_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2ac_con <- plotmetsfl2ac_con[, wt := wt/sum(wt), id_placette]
  
  
  
  plotmetsfl2bc_con=plotmetsfl2_con[plotmetsfl2_con[, all(cl != 'bc')| (cl == 'bc' & .BY %in% cps)|!.BY %in% cps,
                                                    by = id_placette]$V1]
  
  plotmetsfl2bc_con <- unique(plotmetsfl2bc_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2bc_con <- plotmetsfl2bc_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2bc_con <- plotmetsfl2bc_con[, wt := wt/sum(wt), id_placette]
  #-------------------------------------------------------------------------------------------------------------------------#
  
}


######################################
#############Feuillus#################
###########################################################################################################################

{
  pmetsfl3.feu.all <- pmetsfl3[which(id_placette %in% bauges_db_feu$id_placette)]
  pmetsfl3.feu.all <- unique(pmetsfl3.feu.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl3.feu.all <- pmetsfl3.feu.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl3.feu.all <- pmetsfl3.feu.all[, wt := wt/sum(wt), id_placette]
  pmetsfl3.feu.all$meanang <- as.factor(pmetsfl3.feu.all$meanang)
  pmetsfl3.feu.all <- pmetsfl3.feu.all[pflidr!=0]
  ###########################################################################################################################
  
  plotmetsfl2ab_feu=plotmetsfl2_feu[plotmetsfl2_feu[, .I[cl=="ab"  | all(cl!="ab")], by = id_placette]$V1]
  plotmetsfl2ac_feu=plotmetsfl2_feu[plotmetsfl2_feu[, .I[cl=="ac"  | all(cl!="ac")], by = id_placette]$V1]
  plotmetsfl2bc_feu=plotmetsfl2_feu[plotmetsfl2_feu[, .I[cl=="bc"  | all(cl!="bc")], by = id_placette]$V1]
  
  l1_feu <- unique(plotmetsfl2ab_feu[cl=="ab"]$id_placette)
  l2_feu <- unique(plotmetsfl2ac_feu[cl=="ac"]$id_placette)
  l3_feu <- unique(plotmetsfl2bc_feu[cl=="bc"]$id_placette)
  
  cps <- intersect(intersect(l1_feu, l2_feu), l3_feu)
  
  plotmetsfl2ab_feu=plotmetsfl2[plotmetsfl2[, all(cl != 'ab')| (cl == 'ab' & .BY %in% cps)|!.BY %in% cps, 
                                            by = id_placette]$V1]
  plotmetsfl2ab_feu <- unique(plotmetsfl2ab_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2ab_feu <- plotmetsfl2ab_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2ab_feu <- plotmetsfl2ab_feu[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl2ac_feu=plotmetsfl2[plotmetsfl2[, all(cl != 'ac')| (cl == 'ac' & .BY %in% cps)|!.BY %in% cps,
                                            by = id_placette]$V1]
  plotmetsfl2ac_feu <- unique(plotmetsfl2ac_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2ac_feu <- plotmetsfl2ac_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2ac_feu <- plotmetsfl2ac_feu[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl2bc_feu=plotmetsfl2[plotmetsfl2[, all(cl != 'bc')| (cl == 'bc' & .BY %in% cps)|!.BY %in% cps,
                                            by = id_placette]$V1]
  plotmetsfl2bc_feu <- unique(plotmetsfl2bc_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2bc_feu <- plotmetsfl2bc_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2bc_feu <- plotmetsfl2bc_feu[, wt := wt/sum(wt), id_placette]
  #-----------------------------------------------------------------------------------------------------------------------#
  
}

######################################
################Mixte#################
###########################################################################################################################
{
  pmetsfl3.mix.all <- pmetsfl3[which(id_placette %in% bauges_db_mix$id_placette)]
  pmetsfl3.mix.all <- unique(pmetsfl3.mix.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl3.mix.all <- pmetsfl3.mix.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl3.mix.all <- pmetsfl3.mix.all[, wt := wt/sum(wt), id_placette]
  pmetsfl3.mix.all$meanang <- as.factor(pmetsfl3.mix.all$meanang)
  pmetsfl3.mix.all <- pmetsfl3.mix.all[pflidr!=0]
  
  ###########################################################################################################################
  
  plotmetsfl2ab_mix=plotmetsfl2_mix[plotmetsfl2_mix[, .I[cl=="ab"  | all(cl!="ab")], by = id_placette]$V1]
  plotmetsfl2ac_mix=plotmetsfl2_mix[plotmetsfl2_mix[, .I[cl=="ac"  | all(cl!="ac")], by = id_placette]$V1]
  plotmetsfl2bc_mix=plotmetsfl2_mix[plotmetsfl2_mix[, .I[cl=="bc"  | all(cl!="bc")], by = id_placette]$V1]
  
  l1_mix <- unique(plotmetsfl2ab_mix[cl=="ab"]$id_placette)
  l2_mix <- unique(plotmetsfl2ac_mix[cl=="ac"]$id_placette)
  l3_mix <- unique(plotmetsfl2bc_mix[cl=="bc"]$id_placette)
  
  cps <- intersect(intersect(l1_mix, l2_mix), l3_mix)
  
  plotmetsfl2ab_mix=plotmetsfl2[plotmetsfl2[, all(cl != 'ab')| (cl == 'ab' & .BY %in% cps)|!.BY %in% cps, 
                                            by = id_placette]$V1]
  plotmetsfl2ab_mix <- unique(plotmetsfl2ab_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2ab_mix <- plotmetsfl2ab_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2ab_mix <- plotmetsfl2ab_mix[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl2ac_mix=plotmetsfl2[plotmetsfl2[, all(cl != 'ac')| (cl == 'ac' & .BY %in% cps)|!.BY %in% cps,
                                            by = id_placette]$V1]
  plotmetsfl2ac_mix <- unique(plotmetsfl2ac_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2ac_mix <- plotmetsfl2ac_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2ac_mix <- plotmetsfl2ac_mix[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl2bc_mix=plotmetsfl2[plotmetsfl2[, all(cl != 'bc')| (cl == 'bc' & .BY %in% cps)|!.BY %in% cps,
                                            by = id_placette]$V1]
  plotmetsfl2bc_mix <- unique(plotmetsfl2bc_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2bc_mix <- plotmetsfl2bc_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2bc_mix <- plotmetsfl2bc_mix[, wt := wt/sum(wt), id_placette]
  #-----------------------------------------------------------------------------------------------------------------------#
  
}




###########################################################################################################
##Generate samples######
########################
dbase <- pmetsfl3.feu.all
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
smplstfl3.feu.all <- foreach(i = 1:10000, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  set.seed(i)
  dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
  return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
  # return(sset)
}
stopCluster(clus)
registerDoSEQ()
smplstfl3.feu.all <- matrix(unlist(smplstfl3.feu.all), nrow = length(unique(dbase$id_placette)))
smplstfl3.feu.all <- unique(as.data.table(t(smplstfl3.feu.all)))




time_log <- data.frame()


start <- Sys.time()
dbase <- pmetsfl3.feu.all
fds <- bauges_db_feu
setkey(fds,"id_placette")
idx.lst <- smplstfl3.feu.all

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
baugesfl3.feu.all <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
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

baugesfl3.mdlmets.con.all <- melt(rbindlist(lapply(baugesfl3.con.all, function(x)
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
baugesfl3.mdlmets.feu.all <- melt(rbindlist(lapply(baugesfl3.feu.all, function(x)
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
baugesfl3.mdlmets.mix.all <- melt(rbindlist(lapply(baugesfl3.mix.all, function(x)
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


baugesfl3.mdlmets.con.all <- cbind(baugesfl3.mdlmets.con.all, "exp"=rep("all", nrow(baugesfl3.mdlmets.con.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl3.mdlmets.feu.all <- cbind(baugesfl3.mdlmets.feu.all, "exp"=rep("all", nrow(baugesfl3.mdlmets.feu.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl3.mdlmets.mix.all <- cbind(baugesfl3.mdlmets.mix.all, "exp"=rep("all", nrow(baugesfl3.mdlmets.mix.all)), "id"=rep(rep(1:5000, 1, each=6), 4))

ggplot(data=baugesfl3.mdlmets.feu.all[Forest_attr=="Stem volume"], aes( y=value, colour=Metrics))+
  geom_boxplot()+
  facet_grid(variable~., scales = "free")+
  theme_base()+
  scale_fill_grey()







#############################################################################################################################
plotmetsfl3 <- unique(plotmetsfl3[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl3 <- plotmetsfl3[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl3 <- plotmetsfl3[, wt := wt/sum(wt), id_placette]
plotmetsfl3$rid <- seq(1:nrow(plotmetsfl3)) 
#############################################################################################################################

plotmetsfl3ab=plotmetsfl3[plotmetsfl3[, .I[cl=="ab"  | all(cl!="ab")], by = id_placette]$V1]
plotmetsfl3ac=plotmetsfl3[plotmetsfl3[, .I[cl=="ac"  | all(cl!="ac")], by = id_placette]$V1]
plotmetsfl3bc=plotmetsfl3[plotmetsfl3[, .I[cl=="bc"  | all(cl!="bc")], by = id_placette]$V1]


l1 <- unique(as.character(plotmetsfl3ab[cl=="ab"]$id_placette))
l2 <- unique(as.character(plotmetsfl3ac[cl=="ac"]$id_placette))
l3 <- unique(as.character(plotmetsfl3bc[cl=="bc"]$id_placette))

cps <- intersect(intersect(l1, l2), l3)

plotmetsfl3ab=plotmetsfl3ab[plotmetsfl3ab[, all(cl != 'ab')| (cl == 'ab' & .BY %in% cps)|!.BY %in% cps, 
                                          by = id_placette]$V1]
plotmetsfl3ab <- unique(plotmetsfl3ab[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl3ab <- plotmetsfl3ab[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl3ab <- plotmetsfl3ab[, wt := wt/sum(wt), id_placette]


plotmetsfl3ac=plotmetsfl3ac[plotmetsfl3ac[, all(cl != 'ac')| (cl == 'ac' & .BY %in% cps)|!.BY %in% cps, 
                                          by = id_placette]$V1]
plotmetsfl3ac <- unique(plotmetsfl3ac[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl3ac <- plotmetsfl3ac[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl3ac <- plotmetsfl3ac[, wt := wt/sum(wt), id_placette]


plotmetsfl3bc=plotmetsfl3bc[plotmetsfl3bc[, all(cl != 'bc')| (cl == 'bc' & .BY %in% cps)|!.BY %in% cps, 
                                          by = id_placette]$V1]
plotmetsfl3bc <- unique(plotmetsfl3bc[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl3bc <- plotmetsfl3bc[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl3bc <- plotmetsfl3bc[, wt := wt/sum(wt), id_placette]
####################################################################################################################

#################################################################################################################################
pmetsfl3.all <- plotmetsfl3

#tabulate per plot the number of pcs belonging to each class
tbl <- with(plotmetsfl3, table(id_placette, cl))

#pick all the plots that have atleast one pc belonging to class a, class b and class c each
tbl2 <- tbl[which(tbl[,2]>0 & tbl[,3]>0 & tbl[,6]>0),]

#subset those pcs and their metrics which satisfy the above condition
pmets.withallcls <- pmetsfl3.all[id_placette %in% rownames(tbl2)]

#subset those pcs and their metrics which donot satisfy the previous condition
pmets.woallcls <- pmetsfl3.all[!id_placette %in% rownames(tbl2)]





#pick only those pcs and their metrics which belong to class a
pmetsfl3.clab <- pmets.withallcls[cl=="ab"]

#combine with the remaining plots
pmetsfl3.clab <- rbind(pmetsfl3.clab, pmets.woallcls)

#compute inverse probabilities
pmetsfl3.clab <- unique(pmetsfl3.clab[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl3.clab <- pmetsfl3.clab[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl3.clab <- pmetsfl3.clab[, wt := wt/sum(wt), id_placette]





#pick only those pcs and their metrics which belong to class b
pmetsfl3.clbc <- pmets.withallcls[cl=="bc"]

#combine with the remaining plots
pmetsfl3.clbc <- rbind(pmetsfl3.clbc, pmets.woallcls)

#compute inverse probabilities
pmetsfl3.clbc <- unique(pmetsfl3.clbc[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl3.clbc <- pmetsfl3.clbc[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl3.clbc <- pmetsfl3.clbc[, wt := wt/sum(wt), id_placette]





#pick only those pcs and their metrics which belong to class b
pmetsfl3.clac <- pmets.withallcls[cl=="ac"]

#combine with the remaining plots
pmetsfl3.clac <- rbind(pmetsfl3.clac, pmets.woallcls)

#compute inverse probabilities
pmetsfl3.clac <- unique(pmetsfl3.clac[, prob := prop.table(table(cl))[cl], id_placette][])
pmetsfl3.clac <- pmetsfl3.clac[, wt := (1/prob)/(1/sum(prob)), id_placette]
pmetsfl3.clac <- pmetsfl3.clac[, wt := wt/sum(wt), id_placette]



###########################################################################################################
###########################################################################################################
##Generate samples######
########################
dbase <- pmetsfl3.all
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
smplst.3all <- foreach(i = 1:10000, .packages=c("dplyr", "data.table", "caret", "sampling")) %dopar% {
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
  mets_for_model <- plotmetsfl3bc[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
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
             "index" = which(plotmetsfl3bc$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()

set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_vtotfl2bcpfcvlvx <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl3bc[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
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
             "index" = which(plotmetsfl3bc$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()

set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_vtigfl2bcpfcvlvx <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl3bc[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
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
             "index" = which(plotmetsfl3bc$meanang %in% mets_for_model$meanang)))
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











