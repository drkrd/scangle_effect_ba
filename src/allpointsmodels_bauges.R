library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)
library(caret)
library(ggrepel)


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
fd_smry <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db15jul.csv", sep = ",")
colnames(fd_smry)[2] <- "id_placette"
height <- 5


# fd_smry <- fd_smry[, newstratum := ifelse(G175R/G175>0.75 & G175F/G175<0.25 , "Coniferes",
#                                           ifelse(G175F/G175>0.75 & G175R/G175<0.25, "Feuillus", "Mixte"))]

fd_smry <- fd_smry[, newstratum := ifelse(comp_R_G>75 & comp_F_G<25 , "Coniferes",
                                          ifelse(comp_F_G>75 & comp_R_G<25, "Feuillus", "Mixte"))]

setkey(fd_smry, "id_placette")



plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/allpoints_fl/")
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
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber, ht=height)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber, ht=height)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht=height)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht=height)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
plotmetsall_orig <- catalog_apply(plots, func_computeall)
plotmetsflall <- rbindlist(plotmetsall_orig)

##############################################################################################################
allpcs <- list.files(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/"),
                     pattern = "*.las",
                     full.names = TRUE)

alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)
names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/allpoints_fl/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox2, 
                          pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", 
                          ht=height))
setDT(voxall)
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang)]
setkeyv(pfcvladvox, c("id_placette"))
setkeyv(plotmetsflall, c("id_placette"))
plotmetsflall <- plotmetsflall[pfcvladvox]
##############################################################################################################
func_refmodls <- function(db, form)
{
  mdl<- train(form,
              data = db,
              method = "lm",
              trControl = trainControl(method="LOOCV"))
  
  obs <- mdl$pred$obs
  pred <- mdl$pred$pred 
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
  bias <- (sum(yobs-ypred))/length(yobs)
  biaspc <- bias*100/mean(yobs)
  return(list("R2" = aR2, 
              "RMSE" = RMSE,
              "MAE" = MAE,
              "MPE" = MPE))
}

fd_smry <- fd_smry[id_placette %in% plotmetsflall$id_placette]
fd_smry1 <- fd_smry[,c("id_placette", "G75", "volume_tige", "volume_total", "newstratum", "stratum",
                       "G175", "volume_tige_175", "volume_total_175")]
fd_smry_con <- fd_smry1[newstratum=="Coniferes"]
fd_smry_feu <- fd_smry1[newstratum=="Feuillus"]
fd_smry_mix <- fd_smry1[stratum=="Mixte"]

mets.all.all <- fd_smry1[plotmetsflall[id_placette %in% fd_smry1$id_placette]]
mets.all.con <- fd_smry_con[plotmetsflall[id_placette %in% fd_smry_con$id_placette]]
mets.all.feu <- fd_smry_feu[plotmetsflall[id_placette %in% fd_smry_feu$id_placette]]
mets.all.mix <- fd_smry_mix[plotmetsflall[id_placette %in% fd_smry_mix$id_placette]]


f1l <- log(G75)~log(meanch)+log(varch)+log(1-pflidr)+log(cvladlidr)
f1v <- log(G75)~log(meanch)+log(varch)+log(1pfsumprof)+log(cvladvox)
f2l <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


mdl<- train(f1l,
            data = mets.all.feu,
            method = "lm",
            trControl = trainControl(method="LOOCV"))
summary(mdl)

obs <- mdl$pred$obs
pred <- mdl$pred$pred
yobs <- exp(obs)
ypred <- exp(pred)
see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
cf <- exp((see^2)/2)
ypred <- ypred*cf
SSE <- sum((yobs-ypred)^2)
SST <- sum((mean(yobs)-yobs)^2)
R2a <- cor(obs, pred)^2
R2 <- 1-(SSE/SST)
MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
RMSE <- sqrt(mean((yobs-ypred)^2))
RMSEpc <- RMSE*100/mean(yobs)
bias <- (sum(yobs-ypred))/length(yobs)
biaspc <- bias*100/mean(yobs)
xx <- as.data.table(cbind(ypred,yobs,pred,obs))
xx <- cbind(xx, mets.all.con$id_placette)

ggplot(data=xx, aes(x=yobs, y=ypred))+geom_abline()+
  geom_point()+  geom_label_repel(aes(label=V2), size=2, alpha=0.3)+
  geom_smooth(method = 'lm', alpha=0.5)+
  coord_fixed(xlim = c(min(yobs,ypred), 
                       max(yobs,ypred)), 
              ylim = c(min(yobs,ypred),
                       max(yobs,ypred)))+
  theme_few()+
  labs(x = "Observed",
       y = "Predicted")+
  annotate("text", 
           x = -Inf, y = Inf, hjust = 0, vjust = 1,
           label = paste0("R² = ", round(R2, 2), "\n",
                          "RMSE = ", round(RMSE, 2), "\n",
                          "RMSE% = ", round(RMSEpc, 2), "\n",
                          "MPE% = ", round(MPE, 2), "\n",
                          "Bias = ", round(bias, 2), "\n",
                          "Bias% = ", round(biaspc, 2)))





plotmetsflall_con <- plotmetsflall[which(id_placette %in% fd_smry_con$id_placette)]
plotmetsflmax_con <- plotmetsfl1_conx[plotmetsfl1_conx[, .I[meanang == max(meanang)], by=id_placette]$V1]
plotmetsflmin_con <- plotmetsfl1_conx[plotmetsfl1_conx[, .I[meanang == min(meanang)], by=id_placette]$V1]

dbflall_con <- fd_smry_con[plotmetsflall_con] 
dbflmax_con <- fd_smry_con[plotmetsflmax_con] 
dbflmin_con <- fd_smry_con[plotmetsflmin_con] 


mdlmets_ref <- list()
#############################################################################################
plotmetsflall_con <- plotmetsflall[which(id_placette %in% fd_smry_con$id_placette)]
plotmetsflmax_con <- plotmetsfl1_conx[plotmetsfl1_conx[, .I[meanang == max(meanang)], by=id_placette]$V1]
plotmetsflmin_con <- plotmetsfl1_conx[plotmetsfl1_conx[, .I[meanang == min(meanang)], by=id_placette]$V1]

dbflall_con <- fd_smry_con[plotmetsflall_con] 
dbflmax_con <- fd_smry_con[plotmetsflmax_con] 
dbflmin_con <- fd_smry_con[plotmetsflmin_con] 

mdlmets_ref[["g_all_old_con"]] <- func_refmodls(dbflall_con, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_all_vox_con"]] <- func_refmodls(dbflall_con, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["g_min_old_con"]] <- func_refmodls(dbflmin_con, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_min_vox_con"]] <- func_refmodls(dbflmin_con, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["g_max_old_con"]] <- func_refmodls(dbflmax_con, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_max_vox_con"]] <- func_refmodls(dbflmax_con, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))

mdlmets_ref[["vtot_all_old_con"]] <- func_refmodls(dbflall_con, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_all_vox_con"]] <- func_refmodls(dbflall_con, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtot_min_old_con"]] <- func_refmodls(dbflmin_con, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_min_vox_con"]] <- func_refmodls(dbflmin_con, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtot_max_old_con"]] <- func_refmodls(dbflmax_con, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_max_vox_con"]] <- func_refmodls(dbflmax_con, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))

mdlmets_ref[["vtig_all_old_con"]] <- func_refmodls(dbflall_con, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_all_vox_con"]] <- func_refmodls(dbflall_con, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtig_min_old_con"]] <- func_refmodls(dbflmin_con, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_min_vox_con"]] <- func_refmodls(dbflmin_con, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtig_max_old_con"]] <- func_refmodls(dbflmax_con, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_max_vox_con"]] <- func_refmodls(dbflmax_con, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))

###################################################################################################
plotmetsflall_feu <- plotmetsflall[which(id_placette %in% fd_smry_feu$id_placette)]
plotmetsflmax_feu <- plotmetsfl1_feux[plotmetsfl1_feux[, .I[meanang == max(meanang)], by=id_placette]$V1]
plotmetsflmin_feu <- plotmetsfl1_feux[plotmetsfl1_feux[, .I[meanang == min(meanang)], by=id_placette]$V1]

dbflall_feu <- fd_smry_feu[plotmetsflall_feu] 
dbflmax_feu <- fd_smry_feu[plotmetsflmax_feu] 
dbflmin_feu <- fd_smry_feu[plotmetsflmin_feu] 

mdlmets_ref[["g_all_old_feu"]] <- func_refmodls(dbflall_feu, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_all_vox_feu"]] <- func_refmodls(dbflall_feu, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["g_min_old_feu"]] <- func_refmodls(dbflmin_feu, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_min_vox_feu"]] <- func_refmodls(dbflmin_feu, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["g_max_old_feu"]] <- func_refmodls(dbflmax_feu, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_max_vox_feu"]] <- func_refmodls(dbflmax_feu, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))

mdlmets_ref[["vtot_all_old_feu"]] <- func_refmodls(dbflall_feu, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_all_vox_feu"]] <- func_refmodls(dbflall_feu, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtot_min_old_feu"]] <- func_refmodls(dbflmin_feu, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_min_vox_feu"]] <- func_refmodls(dbflmin_feu, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtot_max_old_feu"]] <- func_refmodls(dbflmax_feu, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_max_vox_feu"]] <- func_refmodls(dbflmax_feu, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))

mdlmets_ref[["vtig_all_old_feu"]] <- func_refmodls(dbflall_feu, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_all_vox_feu"]] <- func_refmodls(dbflall_feu, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtig_min_old_feu"]] <- func_refmodls(dbflmin_feu, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_min_vox_feu"]] <- func_refmodls(dbflmin_feu, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtig_max_old_feu"]] <- func_refmodls(dbflmax_feu, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_max_vox_feu"]] <- func_refmodls(dbflmax_feu, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
###################################################################################################
plotmetsflall_mix <- plotmetsflall[which(id_placette %in% fd_smry_mix$id_placette)]
plotmetsflmax_mix <- plotmetsfl1_mixx[plotmetsfl1_mixx[, .I[meanang == max(meanang)], by=id_placette]$V1]
plotmetsflmin_mix <- plotmetsfl1_mixx[plotmetsfl1_mixx[, .I[meanang == min(meanang)], by=id_placette]$V1]

dbflall_mix <- fd_smry_mix[plotmetsflall_mix] 
dbflmax_mix <- fd_smry_mix[plotmetsflmax_mix] 
dbflmin_mix <- fd_smry_mix[plotmetsflmin_mix] 

mdlmets_ref[["g_all_old_mix"]] <- func_refmodls(dbflall_mix, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_all_vox_mix"]] <- func_refmodls(dbflall_mix, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["g_min_old_mix"]] <- func_refmodls(dbflmin_mix, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_min_vox_mix"]] <- func_refmodls(dbflmin_mix, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["g_max_old_mix"]] <- func_refmodls(dbflmax_mix, form =  log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["g_max_vox_mix"]] <- func_refmodls(dbflmax_mix, form =  log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))

mdlmets_ref[["vtot_all_old_mix"]] <- func_refmodls(dbflall_mix, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_all_vox_mix"]] <- func_refmodls(dbflall_mix, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtot_min_old_mix"]] <- func_refmodls(dbflmin_mix, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_min_vox_mix"]] <- func_refmodls(dbflmin_mix, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtot_max_old_mix"]] <- func_refmodls(dbflmax_mix, form =  log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtot_max_vox_mix"]] <- func_refmodls(dbflmax_mix, form =  log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))

mdlmets_ref[["vtig_all_old_mix"]] <- func_refmodls(dbflall_mix, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_all_vox_mix"]] <- func_refmodls(dbflall_mix, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtig_min_old_mix"]] <- func_refmodls(dbflmin_mix, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_min_vox_mix"]] <- func_refmodls(dbflmin_mix, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))
mdlmets_ref[["vtig_max_old_mix"]] <- func_refmodls(dbflmax_mix, form =  log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr))
mdlmets_ref[["vtig_max_vox_mix"]] <- func_refmodls(dbflmax_mix, form =  log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox))


mdlmets_ref1 <- rbindlist(mdlmets_ref, idcol = "id")

mdlmets_ref1 <- mdlmets_ref1[, c("for_attr","exp", "mets","forest"):=tstrsplit(id,"_",fixed=T),]






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








