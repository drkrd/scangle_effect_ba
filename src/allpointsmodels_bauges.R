library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)
library(caret)
library(ggrepel)
library(stringr)



####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
placette <- read.csv(file="D:/1_Work/5_Bauges/T1-Carto/Observatoire/Donnees/table_placette_v20190214_placette.csv", header=TRUE, sep=";", dec=",")
bauges.db <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db15jul.csv", sep = ",")
colnames(bauges.db)[2] <- "id_placette"
height <- 5


# fd_smry <- fd_smry[, newstratum := ifelse(G175R/G175>0.75 & G175F/G175<0.25 , "Coniferes",
#                                           ifelse(G175F/G175>0.75 & G175R/G175<0.25, "Feuillus", "Mixte"))]

bauges.db <- bauges.db[, newstratum := ifelse(comp_R_G>75 & comp_F_G<25 , "Coniferes",
                                              ifelse(comp_F_G>75 & comp_R_G<25, "Feuillus", "Mixte"))]

setkey(bauges.db, "id_placette")



plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/allpoints/")
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
  sdvfplidr <- func_sdvfp(las@data$Z, las@data$ReturnNumber, ht=height)
  
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr,
    sdvfplidr = sdvfplidr))
}
plotmetsall_orig <- catalog_apply(plots, func_computeall)
plotmetsflall <- rbindlist(plotmetsall_orig)
pmetsflall <- plotmetsflall
pmetsflall <- pmetsflall[!id_placette%in%c("283_74_ONF", "96_IRSTEA")]

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

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/allpoints/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox2, 
                          pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", 
                          ht=height))
setDT(voxall)
pfcvladvox <- voxall[, .(cvladvox=cv(PADmean, na.rm = TRUE), 
                         sdvfp=sqrt(sum(PADmean*(k1-(sum(k1*PADmean)/sum(PADmean)))^2)/(sum(PADmean)*(length(PADmean[which(PADmean!=0)])-1)/length(PADmean[which(PADmean!=0)]))),
                         pfsumprof=exp(-0.5*sum(PADmean, na.rm = TRUE))), by=.(id_placette, meanang)]
setkeyv(pfcvladvox, c("id_placette"))
setkeyv(pmetsflall, c("id_placette"))
pmetsflall <- pmetsflall[pfcvladvox]
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
  predicted <- exp(pred)
  see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
  cf <- exp((see^2)/2)
  predicted <- predicted*cf
  SSE <- sum((yobs-predicted)^2)
  SST <- sum((mean(yobs)-yobs)^2)
  R2 <- 1-(SSE/SST)
  aR2 <- 1-((1-R2)*(length(predicted)-1)/(length(predicted)-4-1))
  MPE <- (100/length(predicted))*sum((yobs-predicted)/yobs)
  RMSE <- sqrt(mean((yobs-predicted)^2))
  RMSEpc <- RMSE*100/mean(yobs)
  bias <- (sum(yobs-predicted))/length(yobs)
  biaspc <- bias*100/mean(yobs)
  return(list("R2" = aR2, 
              "RMSE" = RMSE,
              "MAE" = MAE,
              "MPE" = MPE))
}

bauges.db <- bauges.db[id_placette %in% pmetsflall$id_placette]
bauges.db <- bauges.db[,c("id_placette", "G75", "volume_tige", "volume_total", 
                          "comp_F_G","comp_R_G","comp_G","newstratum", "stratum")]

bauges.db2 <- readxl::read_xlsx("data/table_placette - PNR74.xlsx", sheet = "placette" )
bauges.db2 <- as.data.table(bauges.db2)
bauges.db74 <- bauges.db2[A_exclure=="N"]
bauges.db <- bauges.db[id_placette%in%bauges.db74$Id_plac]

bauges.db.con <- bauges.db[newstratum=="Coniferes"]
bauges.db.feu <- bauges.db[newstratum=="Feuillus"]
bauges.db.mix <- bauges.db[newstratum=="Mixte"]

pmetsall.con <- bauges.db.con[pmetsflall[id_placette %in% bauges.db.con$id_placette]]
pmetsall.feu <- bauges.db.feu[pmetsflall[id_placette %in% bauges.db.feu$id_placette]]
pmetsall.mix <- bauges.db.mix[pmetsflall[id_placette %in% bauges.db.mix$id_placette]]


f1l <- log(G75)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(G75)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)








{
mdl.g.ref<- train(f1l,
            data = pmetsall.con,
            method = "lm",
            trControl = trainControl(method="LOOCV"))

mdl.g.vox<- train(f1v,
            data = pmetsall.con,
            method = "lm",
            trControl = trainControl(method="LOOCV"))

mdl.vtige.ref<- train(f2l,
            data = pmetsall.con,
            method = "lm",
            trControl = trainControl(method="LOOCV"))

mdl.vtige.vox<- train(f2v,
            data = pmetsall.con,
            method = "lm",
            trControl = trainControl(method="LOOCV"))

mdl.vtot.ref<- train(f3l,
            data = pmetsall.con,
            method = "lm",
            trControl = trainControl(method="LOOCV"))

mdl.vtot.vox<- train(f3v,
                     data = pmetsall.con,
                     method = "lm",
                     trControl = trainControl(method="LOOCV"))

bgs.mdlmets.con.trad <- rbindlist(list(func_mdlmets(mdl.g.ref$pred$obs, mdl.g.ref$pred$pred, "Basal area", "ref"),
                                       func_mdlmets(mdl.g.vox$pred$obs, mdl.g.vox$pred$pred, "Basal area", "vox"),
                                       func_mdlmets(mdl.vtige.ref$pred$obs, mdl.vtige.ref$pred$pred, "Stem volume", "ref"),
                                       func_mdlmets(mdl.vtige.vox$pred$obs, mdl.vtige.vox$pred$pred, "Stem volume", "vox"),
                                       func_mdlmets(mdl.vtot.ref$pred$obs, mdl.vtot.ref$pred$pred, "Total volume", "ref"),
                                       func_mdlmets(mdl.vtot.vox$pred$obs, mdl.vtot.vox$pred$pred, "Total volume", "vox")))


}
{
  mdl.g.ref<- train(f1l,
                    data = pmetsall.feu,
                    method = "lm",
                    trControl = trainControl(method="LOOCV"))
  
  mdl.g.vox<- train(f1v,
                    data = pmetsall.feu,
                    method = "lm",
                    trControl = trainControl(method="LOOCV"))
  
  mdl.vtige.ref<- train(f2l,
                        data = pmetsall.feu,
                        method = "lm",
                        trControl = trainControl(method="LOOCV"))
  
  mdl.vtige.vox<- train(f2v,
                        data = pmetsall.feu,
                        method = "lm",
                        trControl = trainControl(method="LOOCV"))
  
  mdl.vtot.ref<- train(f3l,
                       data = pmetsall.feu,
                       method = "lm",
                       trControl = trainControl(method="LOOCV"))
  
  mdl.vtot.vox<- train(f3v,
                       data = pmetsall.feu,
                       method = "lm",
                       trControl = trainControl(method="LOOCV"))
  
  bgs.mdlmets.feu.trad <- rbindlist(list(func_mdlmets(mdl.g.ref$pred$obs, mdl.g.ref$pred$pred, "Basal area", "ref"),
                                         func_mdlmets(mdl.g.vox$pred$obs, mdl.g.vox$pred$pred, "Basal area", "vox"),
                                         func_mdlmets(mdl.vtige.ref$pred$obs, mdl.vtige.ref$pred$pred, "Stem volume", "ref"),
                                         func_mdlmets(mdl.vtige.vox$pred$obs, mdl.vtige.vox$pred$pred, "Stem volume", "vox"),
                                         func_mdlmets(mdl.vtot.ref$pred$obs, mdl.vtot.ref$pred$pred, "Total volume", "ref"),
                                         func_mdlmets(mdl.vtot.vox$pred$obs, mdl.vtot.vox$pred$pred, "Total volume", "vox")))
  
  
}
{
  mdl.g.ref<- train(f1l,
                    data = pmetsall.mix,
                    method = "lm",
                    trControl = trainControl(method="LOOCV"))
  
  mdl.g.vox<- train(f1v,
                    data = pmetsall.mix,
                    method = "lm",
                    trControl = trainControl(method="LOOCV"))
  
  mdl.vtige.ref<- train(f2l,
                        data = pmetsall.mix,
                        method = "lm",
                        trControl = trainControl(method="LOOCV"))
  
  mdl.vtige.vox<- train(f2v,
                        data = pmetsall.mix,
                        method = "lm",
                        trControl = trainControl(method="LOOCV"))
  
  mdl.vtot.ref<- train(f3l,
                       data = pmetsall.mix,
                       method = "lm",
                       trControl = trainControl(method="LOOCV"))
  
  mdl.vtot.vox<- train(f3v,
                       data = pmetsall.mix,
                       method = "lm",
                       trControl = trainControl(method="LOOCV"))
  bgs.mdlmets.mix.trad <- rbindlist(list(func_mdlmets(mdl.g.ref$pred$obs, mdl.g.ref$pred$pred, "Basal area", "ref"),
                                         func_mdlmets(mdl.g.vox$pred$obs, mdl.g.vox$pred$pred, "Basal area", "vox"),
                                         func_mdlmets(mdl.vtige.ref$pred$obs, mdl.vtige.ref$pred$pred, "Stem volume", "ref"),
                                         func_mdlmets(mdl.vtige.vox$pred$obs, mdl.vtige.vox$pred$pred, "Stem volume", "vox"),
                                         func_mdlmets(mdl.vtot.ref$pred$obs, mdl.vtot.ref$pred$pred, "Total volume", "ref"),
                                         func_mdlmets(mdl.vtot.vox$pred$obs, mdl.vtot.vox$pred$pred, "Total volume", "vox")))
  
  
}



bgs.mdlmets.con.trad <- cbind(bgs.mdlmets.con.trad, "Coniferous")
bgs.mdlmets.feu.trad <- cbind(bgs.mdlmets.feu.trad, "Broadleaved")
bgs.mdlmets.mix.trad <- cbind(bgs.mdlmets.mix.trad, "Mixed")

bgs.mdlmets.trad <- rbind(bgs.mdlmets.con.trad,
                          bgs.mdlmets.feu.trad,
                          bgs.mdlmets.mix.trad)
bgs.mdlmets.trad <- bgs.mdlmets.trad %>% 
  mutate_if(is.numeric, round, digits=2)


bgs.mdlmets.trad <- as.data.table(t(bgs.mdlmets.trad))
write.csv(bgs.mdlmets.trad, "D:/1_Work/Dropbox/2_Publications/2_paper/results/bgs.mdlmets.csv")





ids <- pmetsall.con

bgs.preds.con <- as.data.table(rbind(cbind(ids$id_placette, exp(mdl.g.ref$pred$obs), func_btcor(mdl.g.ref$pred$obs, mdl.g.ref$pred$pred), "Basal area", "ref"),
                                 cbind(ids$id_placette, exp(mdl.g.vox$pred$obs), func_btcor(mdl.g.vox$pred$obs, mdl.g.vox$pred$pred), "Basal area", "vox"),
                                 cbind(ids$id_placette, exp(mdl.vtige.ref$pred$obs), func_btcor(mdl.vtige.ref$pred$obs, mdl.vtige.ref$pred$pred), "Stem volume", "ref"),
                                 cbind(ids$id_placette, exp(mdl.vtige.vox$pred$obs), func_btcor(mdl.vtige.vox$pred$obs, mdl.vtige.vox$pred$pred), "Stem volume", "vox"),
                                 cbind(ids$id_placette, exp(mdl.vtot.ref$pred$obs), func_btcor(mdl.vtot.ref$pred$obs, mdl.vtot.ref$pred$pred), "Total volume", "ref"),
                                 cbind(ids$id_placette, exp(mdl.vtot.vox$pred$obs), func_btcor(mdl.vtot.vox$pred$obs, mdl.vtot.vox$pred$pred), "Total volume", "vox")))




colnames(bgs.preds.con) <- c("id_placette","observed", "predicted", "Forest_attr", "Metrics")
bgs.preds.con$observed <- as.numeric(bgs.preds.con$observed)
bgs.preds.con$predicted <- as.numeric(bgs.preds.con$predicted)





bgs.preds <- bgs.preds.mix
annotate_npc <- function(label, x, y, ...)
{
  ggplot2::annotation_custom(grid::textGrob(
    x = unit(x, "npc"), y = unit(y, "npc"), label = label, ...))
}


sp3 <- ggplot(data=bgs.preds, aes(x=observed, y=predicted, colour=Metrics))+
  geom_line(aes(group = id_placette), 
            size=0.3*(96/72),
            colour="black",
            linetype="dashed")+
  geom_point(size=1.5)+
  geom_abline()+
  theme_few()+
  facet_wrap(Forest_attr~.)+
  coord_fixed(xlim = c(min(min(bgs.preds$observed),min(bgs.preds$predicted)),max(max(bgs.preds$observed), max(bgs.preds$predicted))),
              ylim = c(min(min(bgs.preds$observed),min(bgs.preds$predicted)),max(max(bgs.preds$observed), max(bgs.preds$predicted))))+
  labs(title="",
       x="Observed",
       y="Predicted")+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")+
  scale_color_jco()
  

ggsave(sp3, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/sp_mix.png", 
       width=24*1.25, height=8*1.25, units="cm", dpi=320) 





obs <- mdl.g.ref$pred$obs
pred <- mdl.g.ref$pred$pred
yobs <- exp(obs)
ypred <-func_btcor(obs, pred)
SSE <- sum((yobs-ypred)^2)
SST <- sum((mean(yobs)-yobs)^2)
R2 <- 1-(SSE/SST)
1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))





















xy <- data.table(id=mets.all.all[!id_placette=="171_IRSTEA"]$id_placette, pred, obs, predicted, yobs)
xy <- xy[!id%in%c("268_74_ONF", "259_74_ONF")]
yobs <- exp(xy$obs)
predicted <- exp(xy$pred)
see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
cf <- exp((see^2)/2)
predicted <- predicted*cf
SSE <- sum((yobs-predicted)^2)
SST <- sum((mean(yobs)-yobs)^2)
R2 <- 1-(SSE/SST)
1-((1-R2)*(length(predicted)-1)/(length(predicted)-4-1))
(100/length(predicted))*sum((yobs-predicted)/yobs)
sqrt(mean((yobs-predicted)^2))
(sum(yobs-predicted))/length(yobs)




RMSEpc <- RMSE*100/mean(yobs)
bias <- (sum(yobs-predicted))/length(yobs)
biaspc <- bias*100/mean(yobs)
xx <- as.data.table(cbind(predicted,yobs,pred,obs))
xx <- cbind(xx, mets.all.con$id_placette)

ggplot(data=xx, aes(x=yobs, y=predicted))+geom_abline()+
  geom_point()+  geom_label_repel(aes(label=V2), size=2, alpha=0.3)+
  geom_smooth(method = 'lm', alpha=0.5)+
  coord_fixed(xlim = c(min(yobs,predicted), 
                       max(yobs,predicted)), 
              ylim = c(min(yobs,predicted),
                       max(yobs,predicted)))+
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








