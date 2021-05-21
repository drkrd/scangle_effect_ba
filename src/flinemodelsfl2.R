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

typefor <- " "
typepc <- "allpoints"
# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
# plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/flightlines_2/")
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
plotmetsfl2_orig <- catalog_apply(plots, func_computeall)
#############################################################################################################
plotmetsfl2 <- plotmetsfl2_orig
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
plotmetsfl2 <- plotmetsfl2[cl1 != "e" & cl2 != "e",]

plotmetsfl2 <- plotmetsfl2[!fl1 %in% c(37.24, 7.16) & !fl2 %in% c(37.24, 7.16)]


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

allvoxfiles <- allvoxfiles[cl1 != "e" & cl2 != "e",]
allvoxfiles <- allvoxfiles[!fl1 %in% c(37.24, 7.16) & !fl2 %in% c(37.24, 7.16)]
voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox2))
setDT(voxall)
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang, pfsumvox)]
setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(plotmetsfl2, c("id_placette", "meanang"))
plotmetsfl2 <- plotmetsfl2[pfcvladvox]
#############################################################################################################################
plotmetsfl2 <- unique(plotmetsfl2[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl2 <- plotmetsfl2[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl2 <- plotmetsfl2[, wt := wt/sum(wt), id_placette]
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


fd_smry <- fread("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",", drop = "id_placette")
colnames(fd_smry)[1] <- "id_placette"
fd_smry[, id_placette := as.factor(id_placette)]
setkey(fd_smry, "id_placette")



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











