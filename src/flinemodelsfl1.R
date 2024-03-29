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
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber)
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
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox2))
setDT(voxall)
voxall <- voxall[!meanang %in% c(37.24, 7.16, 49.7, 43.93)]
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang, pfsumvox)]
setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(plotmetsfl1, c("id_placette", "meanang"))
plotmetsfl1 <- plotmetsfl1[pfcvladvox]
#############################################################################################################################
plotmetsfl1 <- unique(plotmetsfl1[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1 <- plotmetsfl1[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1 <- plotmetsfl1[, wt := wt/sum(wt), id_placette]
plotmetsfl1$meanang <- as.factor(plotmetsfl1$meanang)
##############################################################################################################################


plotmetsfl1a=plotmetsfl1[plotmetsfl1[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
plotmetsfl1b=plotmetsfl1[plotmetsfl1[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
plotmetsfl1c=plotmetsfl1[plotmetsfl1[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]

l1 <- unique(plotmetsfl1a[cl=="a"]$id_placette)
l2 <- unique(plotmetsfl1b[cl=="b"]$id_placette)
l3 <- unique(plotmetsfl1c[cl=="c"]$id_placette)

cps <- intersect(intersect(l1, l2), l3)



func_dffilter <- function(df, cls)
{
  df1 <- df[id_placette %in% cps]
  df1 <- df1[df1[, .I[cl==cls  | all(cl!=cls)], by = id_placette]$V1]
  df2 <- df[!id_placette %in% cps]
  dfn <- rbind(df1, df2)
  return(dfn)
}


plotmetsfl1a <- func_dffilter(plotmetsfl1, "a")
plotmetsfl1a <- unique(plotmetsfl1a[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1a <- plotmetsfl1a[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1a <- plotmetsfl1a[, wt := wt/sum(wt), id_placette]


plotmetsfl1b <- func_dffilter(plotmetsfl1, "b")
plotmetsfl1b <- unique(plotmetsfl1b[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1b <- plotmetsfl1b[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1b <- plotmetsfl1b[, wt := wt/sum(wt), id_placette]

plotmetsfl1c <- func_dffilter(plotmetsfl1, "c")
plotmetsfl1c <- unique(plotmetsfl1c[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1c <- plotmetsfl1c[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1c <- plotmetsfl1c[, wt := wt/sum(wt), id_placette]


plotmetsfl1a=plotmetsfl1[plotmetsfl1[, all(cl != 'a')| (cl == 'a' & .BY %in% cps)|!.BY %in% cps,
                                     by = id_placette]$V1]

plotmetsfl1a <- unique(plotmetsfl1a[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1a <- plotmetsfl1a[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1a <- plotmetsfl1a[, wt := wt/sum(wt), id_placette]

plotmetsfl1b=plotmetsfl1[plotmetsfl1[, all(cl != 'b')| (cl == 'b' & .BY %in% cps)|!.BY %in% cps,
                                     by = id_placette]$V1]
plotmetsfl1b <- unique(plotmetsfl1b[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1b <- plotmetsfl1b[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1b <- plotmetsfl1b[, wt := wt/sum(wt), id_placette]

plotmetsfl1c=plotmetsfl1[plotmetsfl1[, all(cl != 'c')| (cl == 'c' & .BY %in% cps)|!.BY %in% cps,
                                     by = id_placette]$V1]
plotmetsfl1c <- unique(plotmetsfl1c[, prob := prop.table(table(cl))[cl], id_placette][])
plotmetsfl1c <- plotmetsfl1c[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl1c <- plotmetsfl1c[, wt := wt/sum(wt), id_placette]




fd_smry <- fread("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",", drop = "id_placette")
colnames(fd_smry)[1] <- "id_placette"
fd_smry[, id_placette := as.factor(id_placette)]
setkey(fd_smry, "id_placette")



n <- 5000
##Simulations for basal area
set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_gfl1cpfcvlvx <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl1c[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
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
             "index" = which(plotmetsfl1c$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()

##Simulations for total volume
set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_vtotfl1cpfcvlvx <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl1c[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
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
             "index" = which(plotmetsfl1c$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()

##Simulations for stem volume
set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_vtigfl1cpfcvlvx <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl1c[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
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
             "index" = which(plotmetsfl1c$meanang %in% mets_for_model$meanang)))
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


mdlmets_vtotfl1b <- unique(rbindlist(lapply(mdls_vtigfl1b, 
                                            func_extractmdlmets, 
                                            type="One flight line", 
                                            exp="mostly B", 
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
  
  








