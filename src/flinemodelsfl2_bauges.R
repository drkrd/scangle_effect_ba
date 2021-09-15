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
  pmetsfl2.con.all$meanang <- as.factor(pmetsfl2.con.all$meanang)
  pmetsfl2.con.all <- pmetsfl2.con.all[pflidr!=0]
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
  pmetsfl2.feu.all <- pmetsfl2[which(id_placette %in% bauges_db_feu$id_placette)]
  pmetsfl2.feu.all <- unique(pmetsfl2.feu.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl2.feu.all <- pmetsfl2.feu.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl2.feu.all <- pmetsfl2.feu.all[, wt := wt/sum(wt), id_placette]
  pmetsfl2.feu.all$meanang <- as.factor(pmetsfl2.feu.all$meanang)
  pmetsfl2.feu.all <- pmetsfl2.feu.all[pflidr!=0]
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
  pmetsfl2.mix.all <- pmetsfl2[which(id_placette %in% bauges_db_mix$id_placette)]
  pmetsfl2.mix.all <- unique(pmetsfl2.mix.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl2.mix.all <- pmetsfl2.mix.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl2.mix.all <- pmetsfl2.mix.all[, wt := wt/sum(wt), id_placette]
  pmetsfl2.mix.all$meanang <- as.factor(pmetsfl2.mix.all$meanang)
  pmetsfl2.mix.all <- pmetsfl2.mix.all[pflidr!=0]
  
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







time_log <- data.frame()


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
baugesfl2.mdlmets.feu.all <- cbind(baugesfl2.mdlmets.feu.all, "exp"=rep("all", nrow(baugesfl2.mdlmets.feu.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl2.mdlmets.mix.all <- cbind(baugesfl2.mdlmets.mix.all, "exp"=rep("all", nrow(baugesfl2.mdlmets.mix.all)), "id"=rep(rep(1:5000, 1, each=6), 4))


ggplot(data=baugesfl2.mdlmets.feu.all[Forest_attr=="Stem volume"], aes( y=value, colour=Metrics))+
  geom_boxplot()+
  facet_grid(variable~., scales = "free")+
  theme_base()+
  scale_fill_grey()



















  start <- Sys.time()
  type <- "vox"
  dbase <- plotmetsfl2ac_feu
  fds <- fd_smry_feu
  n1 <- "ac"
  n2 <- "feu"
  n <- 3500
  func_extractmdlmets <- function(x, type, exp, mets)
  {
    ind <- x[["index"]]
    # dbase <- plotmetsfl2[ind]
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
    f1 <- log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
    f2 <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
    f3 <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
    start <- Sys.time()
    ##Simulations for basal area
    set.seed(123)
    clus <- makeCluster(detectCores() - 1)
    registerDoParallel(clus, cores = detectCores() - 1)
    nme1a <- paste0("mdls_gfl2", n1, n2)
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
    
    
    nme1b <- paste0("mdlmets_gfl2", n1, n2)
    assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
                                          func_extractmdlmets, 
                                          type="One flight line", 
                                          exp=n1, 
                                          mets="old"))))
    print("BA done")
    
    
    ##Simulations for total volume
    set.seed(123)
    clus <- makeCluster(detectCores() - 1)
    registerDoParallel(clus, cores = detectCores() - 1)
    nme2a <- paste0("mdls_vtotfl2", n1, n2)
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
    nme2b <- paste0("mdlmets_vtotfl2", n1, n2)
    assign(nme2b, unique(rbindlist(lapply(get(nme2a), 
                                          func_extractmdlmets, 
                                          type="One flight line", 
                                          exp=n1, 
                                          mets="old"))))
    print("Vtot done")
    
    ##Simulations for stem volume
    set.seed(123)
    clus <- makeCluster(detectCores() - 1)
    registerDoParallel(clus, cores = detectCores() - 1)
    nme3a <- paste0("mdls_vtigfl2", n1, n2)
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
    nme3b <- paste0("mdlmets_vtigfl2", n1, n2)
    assign(nme3b, unique(rbindlist(lapply(get(nme3a), 
                                          func_extractmdlmets, 
                                          type="One flight line", 
                                          exp=n1, 
                                          mets="old"))))
    print("Vtig done")
    stop <- Sys.time()
  }else
  {
    f1 <- log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
    f2 <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
    f3 <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
    start <- Sys.time()
    ##Simulations for basal area
    set.seed(123)
    clus <- makeCluster(detectCores() - 1)
    registerDoParallel(clus, cores = detectCores() - 1)
    nme1a <- paste0("mdls_gfl2vx", n1, n2)
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
    nme1b <- paste0("mdlmets_gfl2vx", n1, n2)
    assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
                                          func_extractmdlmets, 
                                          type="One flight line", 
                                          exp=n1, 
                                          mets="vox"))))
    print("BA done")
    
    
    ##Simulations for total volume
    set.seed(123)
    clus <- makeCluster(detectCores() - 1)
    registerDoParallel(clus, cores = detectCores() - 1)
    nme2a <- paste0("mdls_vtotfl2vx", n1, n2)
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
    nme2b <- paste0("mdlmets_vtotfl2vx", n1, n2)
    assign(nme2b, unique(rbindlist(lapply(get(nme2a), 
                                          func_extractmdlmets, 
                                          type="One flight line", 
                                          exp=n1, 
                                          mets="vox"))))
    print("Vtot done")
    
    
    ##Simulations for stem volume
    set.seed(123)
    clus <- makeCluster(detectCores() - 1)
    registerDoParallel(clus, cores = detectCores() - 1)
    nme3a <- paste0("mdls_vtigfl2vx", n1, n2)
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
    nme3b <- paste0("mdlmets_vtigfl2vx", n1, n2)
    assign(nme3b, unique(rbindlist(lapply(get(nme3a), 
                                          func_extractmdlmets, 
                                          type="One flight line", 
                                          exp=n1, 
                                          mets="vox"))))
    print("Vtig done")
    stop <- Sys.time()
  }




  
  mdlmets_g <- rbind(mdlmets_gfl2allcon, 
                     mdlmets_gfl2vxallcon,
                     mdlmets_gfl2abcon,
                     mdlmets_gfl2vxabcon,
                     mdlmets_gfl2accon,
                     mdlmets_gfl2vxaccon,
                     mdlmets_gfl2bccon,
                     mdlmets_gfl2vxbccon)
  mdlmets_g <- melt(mdlmets_g, id.vars = c("type", "exp", "mets"))
  mdlmets_g$for_attr <- rep("Basal area", nrow(mdlmets_g))
  
  mdlmets_vtot <- rbind(mdlmets_vtotfl2allcon, 
                        mdlmets_vtotfl2vxallcon,
                        mdlmets_vtotfl2abcon,
                        mdlmets_vtotfl2vxabcon,
                        mdlmets_vtotfl2accon,
                        mdlmets_vtotfl2vxaccon,
                        mdlmets_vtotfl2bccon,
                        mdlmets_vtotfl2vxbccon)
  mdlmets_vtot <- melt(mdlmets_vtot, id.vars = c("type", "exp", "mets"))
  mdlmets_vtot$for_attr <- rep("Total volume", nrow(mdlmets_vtot))
  
  mdlmets_vtig <- rbind(mdlmets_vtigfl2allcon, 
                        mdlmets_vtigfl2vxallcon,
                        mdlmets_vtigfl2abcon,
                        mdlmets_vtigfl2vxabcon,
                        mdlmets_vtigfl2accon,
                        mdlmets_vtigfl2vxaccon,
                        mdlmets_vtigfl2bccon,
                        mdlmets_vtigfl2vxbccon)
  mdlmets_vtig <- melt(mdlmets_vtig, id.vars = c("type", "exp", "mets"))
  mdlmets_vtig$for_attr <- rep("Stem volume", nrow(mdlmets_vtig))
  
  mdlmets_all <- rbind(mdlmets_g, mdlmets_vtot, mdlmets_vtig)
  


ggplot(data=mdlmets_all, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "")+
  theme_base(base_size = 20)+
  facet_grid(variable~for_attr, scales = "free")+
  geom_hline(data = allmdlmets, aes(yintercept=value, linetype=type, colour=mets), size=1.5)+
  scale_fill_pander()








