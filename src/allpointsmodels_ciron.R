library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)
library(caret)
library(stringr)
library(ggsci)



####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
fd_smry <- fread("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",", drop = "id_placette")
colnames(fd_smry)[1] <- "id_placette"
fd_smry$id_placette <- as.character(fd_smry$id_placette)
setkey(fd_smry, "id_placette")


height <- 5

plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/allpoints/")
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
plotmetsflall <- catalog_apply(plots, func_computeall)
plotmetsflall <- rbindlist(plotmetsflall)


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

pmetsflall <- plotmetsflall
allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/2_Ciron/voxelisation/Results/May/wo_interpolation/voxfiles/allpoints/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))


allvoxfiles[, id_placette := sub("\\@.*", "", basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@", "", basename(file_path_sans_ext(V1)))]
voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox2, 
                          pth ="D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/", 
                          ht=height))
setDT(voxall)
pfcvladvox <- voxall[, .(cvladvox=cv(PADmean, na.rm = TRUE), 
                         sdvfp=sqrt(sum(PADmean*(k1-(sum(k1*PADmean)/sum(PADmean)))^2)/(sum(PADmean)*(length(PADmean[which(PADmean!=0)])-1)/length(PADmean[which(PADmean!=0)]))),
                         pfsumprof=exp(-0.5*sum(PADmean, na.rm = TRUE))), by=.(id_placette, meanang)]
setkeyv(pfcvladvox, "id_placette")
setkey(pmetsflall, "id_placette")
pmetsflall <- pmetsflall[pfcvladvox]
mets_all2 <- ciron_db[pmetsflall]
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
  # print(yobs)
  # print(ypred)
  plot(yobs, ypred)
  R2 <- cor(ypred, yobs)^2
  aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
  MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
  RMSE <- sqrt(mean((yobs-ypred)^2))
  MAE <- mean(abs(ypred-yobs))
  return(list("R2" = aR2, 
              "RMSE" = RMSE,
              "MAE" = MAE,
              "MPE" = MPE))
}


refmdls <- list()

colnames(mets_all)[4] <- "vtige_m3_ha"


f1l <- log(g_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(g_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(vtot_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(vtot_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(vtige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(vtige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


# fit <- NULL
# for(row in 1:nrow(mets_all))
# {
#   train <- mets_all[-row]
#   model <- lm(data=train, formula = form.logv)
#   val <- mets_all[row]
#   fit[row] <- predict(model, val)
# }
# obs <- log(mets_all$)


mdl<- train(f3v,
            data = mets_all2[-c(6,7)],
            method = "lm",
            trControl = trainControl(method="LOOCV"))

obs <- mdl$pred$obs
pred <- mdl$pred$pred
yobs <- obs
ypred <- pred


#correction for back-transformation
# see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
# cf <- exp((see^2)/2)
# ypred <- ypred*cf
SST <- sum((mean(yobs)-yobs)^2)
SSE <- sum((yobs-ypred)^2)
SSR <- sum((ypred-mean(yobs))^2)
R2 <- 1-(SSE/SST)
aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
RMSE <- sqrt(mean((yobs-ypred)^2))
RMSEpc <- RMSE*100/mean(yobs)

{
  mdl.g.ref<- train(f1l,
                    data = mets_all2[-c(6,7)],
                    method = "lm",
                    trControl = trainControl(method="LOOCV"))
  
  mdl.g.vox<- train(f1v,
                    data = mets_all2[-c(6,7)],
                    method = "lm",
                    trControl = trainControl(method="LOOCV"))
  
  mdl.vtige.ref<- train(f3l,
                        data = mets_all2[-c(6,7)],
                        method = "lm",
                        trControl = trainControl(method="LOOCV"))
  
  mdl.vtige.vox<- train(f3v,
                        data = mets_all2[-c(6,7)],
                        method = "lm",
                        trControl = trainControl(method="LOOCV"))
  
  mdl.vtot.ref<- train(f2l,
                       data = mets_all2[-c(6,7)],
                       method = "lm",
                       trControl = trainControl(method="LOOCV"))
  
  mdl.vtot.vox<- train(f2v,
                       data = mets_all2[-c(6,7)],
                       method = "lm",
                       trControl = trainControl(method="LOOCV"))
  
  bgs.mdlmets.con.trad <- rbindlist(list(func_mdlmets(mdl.g.ref$pred$obs, mdl.g.ref$pred$pred, "Basal area", "ref"),
                                         func_mdlmets(mdl.g.vox$pred$obs, mdl.g.vox$pred$pred, "Basal area", "vox"),
                                         func_mdlmets(mdl.vtige.ref$pred$obs, mdl.vtige.ref$pred$pred, "Stem volume", "ref"),
                                         func_mdlmets(mdl.vtige.vox$pred$obs, mdl.vtige.vox$pred$pred, "Stem volume", "vox"),
                                         func_mdlmets(mdl.vtot.ref$pred$obs, mdl.vtot.ref$pred$pred, "Total volume", "ref"),
                                         func_mdlmets(mdl.vtot.vox$pred$obs, mdl.vtot.vox$pred$pred, "Total volume", "vox")))
  
  
}

crn.mdlmets.rip.trad <- cbind(bgs.mdlmets.con.trad, "Riparian")
crn.mdlmets.rip.trad <- crn.mdlmets.rip.trad %>% 
  mutate_if(is.numeric, round, digits=2)
crn.mdlmets.rip.trad <- as.data.table(t(crn.mdlmets.rip.trad))
write.csv(crn.mdlmets.rip.trad, "D:/1_Work/Dropbox/2_Publications/2_paper/results/crn.mdlmets.csv")

data <- mets_all2[-c(6,7)]
ggplot(data=data, aes(x=log(pflidr), y=log(pfsumprof)))+geom_point()+
  theme_minimal()+
  geom_abline()+
  coord_fixed(xlim=c(-4,-0.5), ylim = c(-4,-0.5))

ggplot(data=data, aes(x=log(cvladlidr*100), y=log(cvladvox)))+geom_point()+
  theme_minimal()+
  geom_abline()+
  
  coord_fixed(xlim=c(3.75,4.5), ylim = c(3.75,4.5))


ids <- mets_all2[-c(6,7)]

crn.preds.rip <- as.data.table(rbind(cbind(ids$id_placette, exp(mdl.g.ref$pred$obs), func_btcor(mdl.g.ref$pred$obs, mdl.g.ref$pred$pred), "Basal area", "ref"),
                                     cbind(ids$id_placette, exp(mdl.g.vox$pred$obs), func_btcor(mdl.g.vox$pred$obs, mdl.g.vox$pred$pred), "Basal area", "vox"),
                                     cbind(ids$id_placette, exp(mdl.vtige.ref$pred$obs), func_btcor(mdl.vtige.ref$pred$obs, mdl.vtige.ref$pred$pred), "Stem volume", "ref"),
                                     cbind(ids$id_placette, exp(mdl.vtige.vox$pred$obs), func_btcor(mdl.vtige.vox$pred$obs, mdl.vtige.vox$pred$pred), "Stem volume", "vox"),
                                     cbind(ids$id_placette, exp(mdl.vtot.ref$pred$obs), func_btcor(mdl.vtot.ref$pred$obs, mdl.vtot.ref$pred$pred), "Total volume", "ref"),
                                     cbind(ids$id_placette, exp(mdl.vtot.vox$pred$obs), func_btcor(mdl.vtot.vox$pred$obs, mdl.vtot.vox$pred$pred), "Total volume", "vox")))

colnames(crn.preds.rip) <- c("id_placette","observed", "predicted", "Forest_attr", "Metrics")
crn.preds.rip$observed <- as.numeric(crn.preds.rip$observed)
crn.preds.rip$predicted <- as.numeric(crn.preds.rip$predicted)



crn.preds.rip1 <- crn.preds.rip[Forest_attr=="Basal area"]

sp3 <- ggplot(data=crn.preds.rip1, aes(x=observed, y=predicted, colour=Metrics))+
  geom_line(aes(group = id_placette), 
            size=0.3*(96/72),
            colour="black",
            linetype="dashed")+
  geom_point(size=1.5)+
  geom_abline()+
  theme_few()+
  facet_wrap(Forest_attr~.)+
  coord_fixed(xlim = c(min(min(crn.preds.rip1$observed),min(crn.preds.rip1$predicted)),max(max(crn.preds.rip1$observed), max(crn.preds.rip1$predicted))),
              ylim = c(min(min(crn.preds.rip1$observed),min(crn.preds.rip1$predicted)),max(max(crn.preds.rip1$observed), max(crn.preds.rip1$predicted))))+
  labs(title="",
       x="Observed (m2/ha)",
       y="Predicted (m2/ha)")+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")+
  scale_color_jco()

ggsave(sp3, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/sp_rip_g.png", 
       width=8*1.25, height=8*1.25, units="cm", dpi=320) 














