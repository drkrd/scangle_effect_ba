library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)
library(caret)

typefor <- " "
typepc <- "allpoints"
# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
# plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/allpoints/")
opt_independent_files(plots) <- TRUE ####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS

fd_smry <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db_new.csv", sep = ",")
colnames(fd_smry)[2] <- "id_placette"
fd_smry[, id_placette := as.factor(id_placette)]
setkey(fd_smry, "id_placette")
fd_smry_con <- fd_smry[stratum=="Coniferes"]
fd_smry_feu <- fd_smry[stratum=="Feuillus"]
fd_smry_mix <- fd_smry[stratum=="Mixte"]


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

plotmetsall <- catalog_apply(plots, func_computeall)
plotmetsall <- rbindlist(plotmetsall)

##############################################################################################################
plotmetsfl1max_con <- plotmetsfl1[plotmetsfl1[, .I[meanang == max(meanang)], by=id_placette]$V1]
plotmetsfl1min_con <- plotmetsfl1[plotmetsfl1[, .I[meanang == min(meanang)], by=id_placette]$V1]
##############################################################################################################
allpcs <- list.files(paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/"),
                     pattern = "*.las",
                     full.names = TRUE)

alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)
names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/2_Ciron/voxelisation/Results/May/wo_interpolation/voxfiles/allpoints/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox2))
setDT(voxall)
voxall <- voxall[!meanang %in% c(37.24, 7.16, 49.7, 43.93)]
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang, pfsumvox)]
setkeyv(pfcvladvox, c("id_placette"))
setkeyv(plotmetsall, c("id_placette"))
plotmetsall <- plotmetsall[pfcvladvox]
###############################################################################################################
setkey(plotmetsall, "id_placette")
dbase_mets_all <- fd_smry[plotmetsall]
dbase_mets_all_con <- dbase_mets_all[stratum=="Coniferes"]

func_allmodels <- function(db)
{
  mdl_all <- list()
  mdl_all[["G"]] <- train(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                                data = db,
                                method = "lm",
                                trControl = trainControl(method="LOOCV"))
  
  mdl_all[["v_tig"]] <- train(log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                                     data = db,
                                     method = "lm",
                                     trControl = trainControl(method="LOOCV"))
  
  mdl_all[["v_tot"]] <- train(log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                                    data = db,
                                    method = "lm",
                                    trControl = trainControl(method="LOOCV"))
  
  return(list("models" = mdl_all))
}

mdls_all_con <- func_allmodels(dbase_mets_all_con)

allmdlmets <- rbind(mdls_fl1min[[1]]$G$results, 
                    mdls_fl1min[[1]]$v_tot$results, 
                    mdls_fl1min[[1]]$v_tig$results,
                    mdls_fl1minpfcvlvx[[1]]$G$results, 
                    mdls_fl1minpfcvlvx[[1]]$v_tot$results, 
                    mdls_fl1minpfcvlvx[[1]]$v_tig$results,
                    mdls_fl1max[[1]]$G$results, 
                    mdls_fl1max[[1]]$v_tot$results, 
                    mdls_fl1max[[1]]$v_tig$results,
                    mdls_fl1maxpfcvlvx[[1]]$G$results, 
                    mdls_fl1maxpfcvlvx[[1]]$v_tot$results, 
                    mdls_fl1maxpfcvlvx[[1]]$v_tig$results,
                    mdls_all[[1]]$G$results, 
                    mdls_all[[1]]$v_tot$results, 
                    mdls_all[[1]]$v_tig$results,
                    mdls_allpfcvlvx[[1]]$G$results, 
                    mdls_allpfcvlvx[[1]]$v_tot$results, 
                    mdls_allpfcvlvx[[1]]$v_tig$results)

for_attr <- rep(c("Basal area", "Total volume", "Stem volume"), 6)
mets <- rep(c(rep("old", 3), rep("vox", 3)),3)
type <- c(rep("min", 6), rep("max", 6), rep("all", 6))

allmdlmets <- cbind(allmdlmets, type)


allmdlmetst <- rbind(mdls_allt[[1]]$G$results, 
                     mdls_allt[[1]]$v_tot$results, 
                     mdls_allt[[1]]$v_tig$results,
                     mdls_all[[1]]$G$results, 
                     mdls_all[[1]]$v_tot$results, 
                     mdls_all[[1]]$v_tig$results)

allmdlmets <- allmdlmets[-1]

allmdlmets$type <- c(rep("min",3), rep("max",3), rep("all",3))
allmdlmets$met <- c(rep(c("G","vtot","vtig"),3))

names(allmdlmets) <- c("rmse", "r2", "mae", "ref", "met")
allmdlmets <- melt(allmdlmets)
setDT(allmdlmets)
allmdlmets <- allmdlmets[,-1]
names(allmdlmets) <- c("rmse", "r2", "mae", "for_attr", "mets", "type")
allmdlmets <- melt(allmdlmets[], id.vars=c("for_attr", "mets", "type"))








