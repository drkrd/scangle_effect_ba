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

# The following flightlines are problematic. There is no voxelisation for them
# 283_74_ONF_un@6.28
# 96_IRSTEA_un@30

plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/flightlines_1/")
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
plotmetsfl1 <- catalog_apply(plots, func_computeall)
plotmetsfl1 <- rbindlist(plotmetsfl1)
plotmetsfl1 <- plotmetsfl1[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                       ifelse(meanang>=10&meanang<20, "b",
                                              ifelse(meanang>=20&meanang<30, "c",
                                                     ifelse(meanang>=30&meanang<40,"d","e"))))]
plotmetsfl1 <- plotmetsfl1[cl!="e"]
plotmetsfl1 <- plotmetsfl1[-c(419,488)]
plotmetsfl1 <- na.omit(plotmetsfl1)
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

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/flightlines_1/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox2, pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", ht=6))
setDT(voxall)
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), 
                         pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), 
                     by=.(id_placette, meanang, pfsumvox)]

pfcvladvox$meanang <- as.numeric(pfcvladvox$meanang)
pfcvladvox <- pfcvladvox[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                     ifelse(meanang>=10&meanang<20, "b",
                                            ifelse(meanang>=20&meanang<30, "c",
                                                   ifelse(meanang>=30&meanang<40,"d","e"))))]
pfcvladvox <- pfcvladvox[cl!="e"]
setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(plotmetsfl1, c("id_placette", "meanang"))
plotmetsfl1 <- plotmetsfl1[pfcvladvox]


######################################
#############All################
#############################################################################################################################
{
  plotmetsfl1 <- unique(plotmetsfl1[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1 <- plotmetsfl1[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1 <- plotmetsfl1[, wt := wt/sum(wt), id_placette]
  plotmetsfl1$meanang <- as.factor(plotmetsfl1$meanang)
  plotmetsfl1 <- plotmetsfl1[pflidr!=0]
}


######################################
#############Coniferes################
#############################################################################################################################

{
  plotmetsfl1_con <- plotmetsfl1[which(id_placette %in% fd_smry_con$id_placette)]
  plotmetsfl1_con <- unique(plotmetsfl1_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1_con <- plotmetsfl1_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1_con <- plotmetsfl1_con[, wt := wt/sum(wt), id_placette]
  plotmetsfl1_con$meanang <- as.factor(plotmetsfl1_con$meanang)
  plotmetsfl1_con <- plotmetsfl1_con[pflidr!=0]
  ##############################################################################################################################
  
  plotmetsfl1a_con=plotmetsfl1_con[plotmetsfl1_con[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
  plotmetsfl1b_con=plotmetsfl1_con[plotmetsfl1_con[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
  plotmetsfl1c_con=plotmetsfl1_con[plotmetsfl1_con[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]
  
  l1_con <- unique(plotmetsfl1a_con[cl=="a"]$id_placette)
  l2_con <- unique(plotmetsfl1b_con[cl=="b"]$id_placette)
  l3_con <- unique(plotmetsfl1c_con[cl=="c"]$id_placette)
  
  cps <- intersect(intersect(l1_con, l2_con), l3_con)
  
  func_dffilter <- function(df, cls)
  {
    df1 <- df[id_placette %in% cps]
    df1 <- df1[df1[, .I[cl==cls  | all(cl!=cls)], by = id_placette]$V1]
    df2 <- df[!id_placette %in% cps]
    dfn <- rbind(df1, df2)
    return(dfn)
  }
  
  
  plotmetsfl1a_con <- func_dffilter(plotmetsfl1_con, "a")
  plotmetsfl1a_con <- unique(plotmetsfl1a_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1a_con <- plotmetsfl1a_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1a_con <- plotmetsfl1a_con[, wt := wt/sum(wt), id_placette]
  
  
  plotmetsfl1b_con <- func_dffilter(plotmetsfl1_con, "b")
  plotmetsfl1b_con <- unique(plotmetsfl1b_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1b_con <- plotmetsfl1b_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1b_con <- plotmetsfl1b_con[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1c_con <- func_dffilter(plotmetsfl1_con, "c")
  plotmetsfl1c_con <- unique(plotmetsfl1c_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1c_con <- plotmetsfl1c_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1c_con <- plotmetsfl1c_con[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1a_con=plotmetsfl1_con[plotmetsfl1_con[, all(cl != 'a')| (cl == 'a' & .BY %in% cps)|!.BY %in% cps, 
                                                   by = id_placette]$V1]

  plotmetsfl1a_con <- unique(plotmetsfl1a_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1a_con <- plotmetsfl1a_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1a_con <- plotmetsfl1a_con[, wt := wt/sum(wt), id_placette]
  
  
  
  plotmetsfl1b_con=plotmetsfl1_con[plotmetsfl1_con[, all(cl != 'b')| (cl == 'b' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  
  plotmetsfl1b_con <- unique(plotmetsfl1b_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1b_con <- plotmetsfl1b_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1b_con <- plotmetsfl1b_con[, wt := wt/sum(wt), id_placette]
  
  
  
  plotmetsfl1c_con=plotmetsfl1_con[plotmetsfl1_con[, all(cl != 'c')| (cl == 'c' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  
  plotmetsfl1c_con <- unique(plotmetsfl1c_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1c_con <- plotmetsfl1c_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1c_con <- plotmetsfl1c_con[, wt := wt/sum(wt), id_placette]
  #-------------------------------------------------------------------------------------------------------------------------#
  
}


######################################
#############Feuillus#################
###########################################################################################################################

{
  plotmetsfl1_feu <- plotmetsfl1[which(id_placette %in% fd_smry_feu$id_placette)]
  plotmetsfl1_feu <- unique(plotmetsfl1_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1_feu <- plotmetsfl1_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1_feu <- plotmetsfl1_feu[, wt := wt/sum(wt), id_placette]
  plotmetsfl1_feu$meanang <- as.factor(plotmetsfl1_feu$meanang)
  plotmetsfl1_feu <- plotmetsfl1_feu[pflidr!=0]
  ###########################################################################################################################
  
  plotmetsfl1a_feu=plotmetsfl1_feu[plotmetsfl1_feu[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
  plotmetsfl1b_feu=plotmetsfl1_feu[plotmetsfl1_feu[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
  plotmetsfl1c_feu=plotmetsfl1_feu[plotmetsfl1_feu[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]
  
  l1_feu <- unique(plotmetsfl1a_feu[cl=="a"]$id_placette)
  l2_feu <- unique(plotmetsfl1b_feu[cl=="b"]$id_placette)
  l3_feu <- unique(plotmetsfl1c_feu[cl=="c"]$id_placette)
  
  cps <- intersect(intersect(l1_feu, l2_feu), l3_feu)
  
  plotmetsfl1a_feu=plotmetsfl1[plotmetsfl1[, all(cl != 'a')| (cl == 'a' & .BY %in% cps)|!.BY %in% cps, 
                                           by = id_placette]$V1]
  plotmetsfl1a_feu <- unique(plotmetsfl1a_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1a_feu <- plotmetsfl1a_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1a_feu <- plotmetsfl1a_feu[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1b_feu=plotmetsfl1[plotmetsfl1[, all(cl != 'b')| (cl == 'b' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  plotmetsfl1b_feu <- unique(plotmetsfl1b_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1b_feu <- plotmetsfl1b_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1b_feu <- plotmetsfl1b_feu[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1c_feu=plotmetsfl1[plotmetsfl1[, all(cl != 'c')| (cl == 'c' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  plotmetsfl1c_feu <- unique(plotmetsfl1c_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1c_feu <- plotmetsfl1c_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1c_feu <- plotmetsfl1c_feu[, wt := wt/sum(wt), id_placette]
  #-----------------------------------------------------------------------------------------------------------------------#
  
}

######################################
################Mixte#################
###########################################################################################################################
{
  plotmetsfl1_mix <- plotmetsfl1[which(id_placette %in% fd_smry_mix$id_placette)]
  plotmetsfl1_mix <- unique(plotmetsfl1_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1_mix <- plotmetsfl1_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1_mix <- plotmetsfl1_mix[, wt := wt/sum(wt), id_placette]
  plotmetsfl1_mix$meanang <- as.factor(plotmetsfl1_mix$meanang)
  plotmetsfl1_mix <- plotmetsfl1_mix[pflidr!=0]
  
  ###########################################################################################################################
  
  plotmetsfl1a_mix=plotmetsfl1_mix[plotmetsfl1_mix[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
  plotmetsfl1b_mix=plotmetsfl1_mix[plotmetsfl1_mix[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
  plotmetsfl1c_mix=plotmetsfl1_mix[plotmetsfl1_mix[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]
  
  l1_mix <- unique(plotmetsfl1a_mix[cl=="a"]$id_placette)
  l2_mix <- unique(plotmetsfl1b_mix[cl=="b"]$id_placette)
  l3_mix <- unique(plotmetsfl1c_mix[cl=="c"]$id_placette)
  
  cps <- intersect(intersect(l1_mix, l2_mix), l3_mix)
  
  plotmetsfl1a_mix=plotmetsfl1[plotmetsfl1[, all(cl != 'a')| (cl == 'a' & .BY %in% cps)|!.BY %in% cps, 
                                           by = id_placette]$V1]
  plotmetsfl1a_mix <- unique(plotmetsfl1a_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1a_mix <- plotmetsfl1a_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1a_mix <- plotmetsfl1a_mix[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1b_mix=plotmetsfl1[plotmetsfl1[, all(cl != 'b')| (cl == 'b' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  plotmetsfl1b_mix <- unique(plotmetsfl1b_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1b_mix <- plotmetsfl1b_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1b_mix <- plotmetsfl1b_mix[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1c_mix=plotmetsfl1[plotmetsfl1[, all(cl != 'c')| (cl == 'c' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  plotmetsfl1c_mix <- unique(plotmetsfl1c_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1c_mix <- plotmetsfl1c_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1c_mix <- plotmetsfl1c_mix[, wt := wt/sum(wt), id_placette]
  #-----------------------------------------------------------------------------------------------------------------------#
  
}



# library(boot)
# func_boot <- function(data, idx, fd)
# {
#   mets_for_model <- data[idx]
# 
#   setkey(mets_for_model,"id_placette")
# 
# 
#   mets_for_model <- fd[mets_for_model]
# 
# 
#   model <- train(log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladlidr),
#                  data = mets_for_model,
#                  method = "lm",
#                  trControl = trainControl(method="LOOCV"))
#   model$results$Rsquared
# }
# 
# boot(plotmetsfl1_con, statistic = func_boot, R=1000, strata = plotmetsfl1_con$id_placette, fd=fd_smry_con)
# 


start <- Sys.time()
type <- "lidr"
dbase <- plotmetsfl1a_feu
fds <- fd_smry_feu
n1 <- "a"
n2 <- "feu"
n <- 3500
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
  f1 <- log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f2 <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f3 <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
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
  f1 <- log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f2 <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f3 <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
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


mdlmets_vtotfl1_con  <- unique(rbindlist(lapply(mdls_vtotfl1_con, 
                                                    func_extractmdlmets, 
                                                    type="One flight line", 
                                                    exp="allclasses", 
                                                    mets="old")))


mdlmets_vtotfl1vx_con  <- unique(rbindlist(lapply(mdls_vtotfl1vx_con, 
                                                  func_extractmdlmets, 
                                                  type="One flight line", 
                                                  exp="allclasses", 
                                                  mets="vox")))



mdlmets_g <- rbind(mdlmets_gfl1mix, 
                   mdlmets_gfl1vxmix)
mdlmets_g <- melt(mdlmets_g, id.vars = c("type", "exp", "mets"))
mdlmets_g$for_attr <- rep("Basal area", nrow(mdlmets_g))

mdlmets_vtot <- rbind(mdlmets_vtotfl1mix, 
                      mdlmets_vtotfl1vxmix)
mdlmets_vtot <- melt(mdlmets_vtot, id.vars = c("type", "exp", "mets"))
mdlmets_vtot$for_attr <- rep("Total volume", nrow(mdlmets_vtot))

mdlmets_vtig <- rbind(mdlmets_vtigfl1mix, 
                      mdlmets_vtigfl1vxmix)
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
  
  




plots <- list()

for(id_plac in unique(voxall$id_placette))
{
  temp <- voxall[id_placette==id_plac]
  temp <- temp[,c("k1", "m", "meanang", "id_placette")]
  plots[[id_plac]] <- ggplot(data=temp, aes(x=k1, y=m, linetype = meanang))+
    theme_minimal()+
    geom_line(size=0.8)+
    coord_flip()+
    labs(title = unique(temp$id_placette))+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
}
plots_con <- plots[names(plots)%in%fd_smry_con$id_placette]
plots_feu <- plots[names(plots)%in%fd_smry_feu$id_placette]
plots_mix <- plots[names(plots)%in%fd_smry_mix$id_placette]


l1 <- c(1,2,3)
l2 <- c(1,4,6)

