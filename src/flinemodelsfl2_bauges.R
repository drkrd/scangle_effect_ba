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
fd_smry <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db_new.csv", sep = ",")
colnames(fd_smry)[2] <- "id_placette"


fd_smry <- fd_smry[, newstratum := ifelse(G175R/G175>0.75 & G175F/G175<0.25 , "Coniferes", 
                                          ifelse(G175F/G175>0.75 & G175R/G175<0.25, "Feuillus", "Mixte"))]
fd_smry[, stratum := as.factor(stratum)]
fd_smry[, newstratum := as.factor(newstratum)]
fd_smry[, id_placette := as.factor(id_placette)]
setkey(fd_smry, "id_placette")
fd_smry_con <- fd_smry[newstratum=="Coniferes"]
fd_smry_feu <- fd_smry[newstratum=="Feuillus"]
fd_smry_mix <- fd_smry[newstratum=="Mixte"]


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
    mang <- round(as.numeric(mang), 2)
    
  }
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht=6)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht=6)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
plotmetsfl2_orig <- catalog_apply(plots, func_computeall)
#############################################################################################################
plotmetsfl2 <- plotmetsfl2_orig
plotmetsfl2 <- rbindlist(plotmetsfl2)
plotmetsfl2 <- plotmetsfl2[, c("fl2","fl2"):=tstrsplit(meanang,"_",fixed=T),]
plotmetsfl2$fl2 <- as.numeric(plotmetsfl2$fl2)
plotmetsfl2$fl2 <- as.numeric(plotmetsfl2$fl2)
plotmetsfl2 <- plotmetsfl2[, fl2:= ifelse(is.na(fl2), fl2, fl2),]

plotmetsfl2 <- plotmetsfl2[, cl1:=ifelse(fl2>=0&fl2<10, "a",
                                           ifelse(fl2>=10&fl2<20, "b",
                                                  ifelse(fl2>=20&fl2<30, "c",
                                                         ifelse(fl2>=30&fl2<40,"d","e"))))]
plotmetsfl2 <- plotmetsfl2[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                           ifelse(fl2>=10&fl2<20, "b",
                                                  ifelse(fl2>=20&fl2<30, "c",
                                                         ifelse(fl2>=30&fl2<40,"d","e"))))]
plotmetsfl2 <- plotmetsfl2[cl1 != "e" & cl2 != "e",]

plotmetsfl2 <- plotmetsfl2[!(id_placette=="96_IRSTEA" & fl2==30.0)]
plotmetsfl2 <- plotmetsfl2[!(id_placette=="96_IRSTEA" & fl2==30.0)]
plotmetsfl2 <- plotmetsfl2[!(id_placette=="283_74_ONF" & fl2==6.28)]
plotmetsfl2 <- plotmetsfl2[!(id_placette=="283_74_ONF" & fl2==6.28)]



plotmetsfl2 <- plotmetsfl2[, cl := paste0(sort(.SD), collapse = ""), .SDcols = c("cl1", "cl2"),  by = 1:nrow(plotmetsfl2)]


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

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/flightlines_2/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
allvoxfiles <- allvoxfiles[, c("fl2","fl2"):=tstrsplit(meanang,"_",fixed=T),]
allvoxfiles$fl2 <- as.numeric(allvoxfiles$fl2)
allvoxfiles$fl2 <- as.numeric(allvoxfiles$fl2)
allvoxfiles <- allvoxfiles[, cl1:=ifelse(fl2>=0&fl2<10, "a",
                                         ifelse(fl2>=10&fl2<20, "b",
                                                ifelse(fl2>=20&fl2<30, "c",
                                                       ifelse(fl2>=30&fl2<40,"d","e"))))]
allvoxfiles <- allvoxfiles[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                         ifelse(fl2>=10&fl2<20, "b",
                                                ifelse(fl2>=20&fl2<30, "c",
                                                       ifelse(fl2>=30&fl2<40,"d","e"))))]

allvoxfiles <- allvoxfiles[cl1 != "e" & cl2 != "e",]
allvoxfiles <- allvoxfiles[!(id_placette=="96_IRSTEA_un" & fl2==30.0)]
allvoxfiles <- allvoxfiles[!(id_placette=="96_IRSTEA_un" & fl2==30.0)]
allvoxfiles <- allvoxfiles[!(id_placette=="283_74_ONF_un" & fl2==6.28)]
allvoxfiles <- allvoxfiles[!(id_placette=="283_74_ONF_un" & fl2==6.28)]
voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox2, pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", ht=6))
setDT(voxall)
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), 
                         pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), 
                     by=.(id_placette, meanang, pfsumvox)]
setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(plotmetsfl2, c("id_placette", "meanang"))
plotmetsfl2 <- plotmetsfl2[pfcvladvox]

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
  plotmetsfl2_con <- plotmetsfl2[which(id_placette %in% fd_smry_con$id_placette)]
  plotmetsfl2_con <- unique(plotmetsfl2_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2_con <- plotmetsfl2_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2_con <- plotmetsfl2_con[, wt := wt/sum(wt), id_placette]
  plotmetsfl2_con$meanang <- as.factor(plotmetsfl2_con$meanang)
  plotmetsfl2_con <- plotmetsfl2_con[pflidr!=0]
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
  plotmetsfl2_feu <- plotmetsfl2[which(id_placette %in% fd_smry_feu$id_placette)]
  plotmetsfl2_feu <- unique(plotmetsfl2_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2_feu <- plotmetsfl2_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2_feu <- plotmetsfl2_feu[, wt := wt/sum(wt), id_placette]
  plotmetsfl2_feu$meanang <- as.factor(plotmetsfl2_feu$meanang)
  plotmetsfl2_feu <- plotmetsfl2_feu[pflidr!=0]
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
  plotmetsfl2_mix <- plotmetsfl2[which(id_placette %in% fd_smry_mix$id_placette)]
  plotmetsfl2_mix <- unique(plotmetsfl2_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl2_mix <- plotmetsfl2_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl2_mix <- plotmetsfl2_mix[, wt := wt/sum(wt), id_placette]
  plotmetsfl2_mix$meanang <- as.factor(plotmetsfl2_mix$meanang)
  plotmetsfl2_mix <- plotmetsfl2_mix[pflidr!=0]
  
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








