library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)

typefor <- "feuillus"
typepc <- "flightlines"
# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
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
    meanang = as.factor(mang),
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}

plotmets <- catalog_apply(plots, func_computeall)
plotmets <- rbindlist(plotmets)
keycols = c("id_placette", "meanang")
setkeyv(plotmets, keycols)


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






func_normvox2 <- function(x)
{
  txt <- readLines(x[1], n=3)[2:3]
  zmin <- as.numeric(unlist(strsplit(txt[[1]], "\\s+"))[4])
  voxtbl <- fread(x[1], na.strings = "NA" , skip = 5)
  
  
  
  #align DTM with voxel file
  lasnm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/", x[2], "@all.las")
  ls <- readLASheader(lasnm)
  empty_raster <- raster(ncol=round(ls@PHB[["Max X"]]-ls@PHB[["Min X"]]), 
                         nrow=round(ls@PHB[["Max Y"]]-ls@PHB[["Min Y"]]),
                         xmn=ls@PHB[["Min X"]], xmx= ls@PHB[["Max X"]], 
                         ymn=ls@PHB[["Min Y"]], ymx= ls@PHB[["Max Y"]])
  dt <- alldtms[[paste0(x[2],"@all")]]
  dt <- resample(dt, empty_raster, method='bilinear')
  dt_mat <- as.matrix(dt)
  
  
  #normalisation based on the DTM
  voxtbl <- voxtbl[,1:4][,alt := k+zmin+0.5][, dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  voxtbl <- voxtbl[alt>dtm]
  voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
  voxtbl <- voxtbl[, k1:= k1-0.5]
  voxtbl1 <- voxtbl[k1>2]
  ttl <- sum(voxtbl1$PadBVTotal, na.rm = TRUE)
  ttl <- ttl/(pi*15*15)
  pf <- exp(-0.5*ttl)
  voxtbl <- voxtbl[, .(m=mean(PadBVTotal, na.rm = TRUE), s=sum(is.nan(PadBVTotal))/length(PadBVTotal)), by=list(k1)]
  
  # voxtbl <- voxtbl[, .(s1 = sum(!is.nan(y))), by=list(k1)]
  # voxtbl <- voxtbl[, .(s1 = sum(is.nan(PadBVTotal))/length(PadBVTotal)), by=list(k1)]
  # voxtbl <- voxtbl[, .(s2 = sum(is.nan(PadBVTotal))), by=list(k1)]
  # 
  voxtbl <- voxtbl[, c("id_placette","meanang") := .(as.factor(x[2]), as.factor(x[3])),]
  voxtbl <- voxtbl[, pfsumvox:=pf,]
  voxtbl <- voxtbl[1:max(which(m!=0)),]
  voxtbl <- voxtbl[k1>2,]
  
  return(voxtbl)
}


func_normvox3 <- function(x)
{
  txt <- readLines(x[1], n=4)[2:4]
  zmin <- as.numeric(unlist(strsplit(txt[[1]], "\\s+"))[4])
  #zmin <- round(zmin)
  voxtbl <- fread(x[1], na.strings = "NA" , skip = 5)
  
  
  
  #align DTM with voxel file
  xmin <- as.numeric(unlist(strsplit(txt[1], "\\s+"))[2])
  xmax <- as.numeric(unlist(strsplit(txt[2], "\\s+"))[2])
  ymin <- as.numeric(unlist(strsplit(txt[1], "\\s+"))[3])
  ymax <- as.numeric(unlist(strsplit(txt[2], "\\s+"))[3])
  rows <- as.numeric(unlist(strsplit(txt[3], "\\s+"))[3])
  cols <- as.numeric(unlist(strsplit(txt[3], "\\s+"))[2])
  
  empty_raster <- raster(ncol = cols, nrow = rows,
                         xmn = xmin, xmx = xmax, 
                         ymn = ymin, ymx = ymax)
  
  dt <- alldtms[[paste0(x[2],"@all")]]
  dt <- resample(dt, empty_raster, method='bilinear')
  dt <- round(dt)
  dt_mat <- as.matrix(dt)
  
  
  #normalisation based on the DTM
  voxtbl <- voxtbl[, 1:4][, alt := k+zmin+0.5][, dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  voxtbl <- voxtbl[alt>dtm]
  voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
  voxtbl <- voxtbl[, k1:= k1-0.5]
  ttl <- sum(voxtbl1$PadBVTotal, na.rm = TRUE)
  ttl <- ttl/(pi*15*15)
  pf <- exp(-0.5*ttl)
  voxtbl <- voxtbl[, .(m = mean(PadBVTotal, na.rm = TRUE)), by = k1]
  voxtbl <- voxtbl[k1>2]
  
  # voxtbl <- voxtbl[, .(s1 = sum(!is.nan(y))), by=list(k1)]
  # voxtbl <- voxtbl[, .(s1 = sum(is.nan(PadBVTotal))/length(PadBVTotal)), by=list(k1)]
  # voxtbl <- voxtbl[, .(s2 = sum(is.nan(PadBVTotal))), by=list(k1)]
  # 
  voxtbl <- voxtbl[, c("id_placette","meanang") := .(as.factor(x[2]), as.factor(x[3])),]
  voxtbl <- voxtbl[, pfsumvox:=pf,]
  voxtbl <- voxtbl[1:max(which(m!=0)),]
  return(voxtbl)
}



voxall <- rbindlist(apply(allvoxfiles, 1, func_normvox3))
setDT(voxall)
voxall <- voxall[!meanang %in% c(37.24, 7.16, 49.7, 43.93)]
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang, pfsumvox)]
setkeyv(pfcvladvox, "id_placette")
plotmets<- plotmetsfl[pfcvladvox]

for(plot in unique(voxall$id_placette))
{
  tmp <- voxall[,`:=`(pfsumvox=NULL)]
  tmp <- tmp[id_placette==plot]
  tmp <- reshape2::melt(tmp, measure.vars=c("m","s"))
  plt <- ggplot(data=tmp, aes(x=k1, y=value, colour=meanang))+
    geom_line()+facet_grid(~variable)+coord_flip()
  ggsave(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/", typefor ,"/profiles/navoxels/", plot, ".png"),
         plt, device = "png", dpi = "retina")
}


voxall1 <- voxall[which(max(s)<0.8),,by(meanang)] 

voxall1 <-voxall %>% 
  group_by(meanang) %>% 
  filter(voxall[which(voxall1$s<0.8),])

##########################
#Compute cvlad from voxels
##########################

cvladvox <- lapply(voxmergedall, function(x) return(x[k1>2,.(cvs=cv(m, na.rm = TRUE)), by=meanang]))
cvladvox <- lapply(voxmergedall, function(x){
  x <- x[k1>2,.(pad=sum(m, na.rm = TRUE)/(pi*15*15)), by=meanang]
})


cvladvox <- voxmergedall[, .(cvladvox=cv(m, na.rm = TRUE), pfvox=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette,meanang)]

cvladvox <- rbindlist(cvladvox, use.names = TRUE, idcol = TRUE)
names(cvladvox) <- c("id_placette", "meanang", "cvladvox") 
cvladvox$id_placette <- as.factor(cvladvox$id_placette)
setkeyv(cvladvox, keycols)
plotmets<- plotmets[cvladvox]

#-------------------------
#-------------------------


#####################################################
#Compute pf from voxels by inverting beer-lambert law
#####################################################

func_computevoxpf <- function(x, zmin, dt_mat)
{
  id_placette <- x[2]
  meanang <- x[3]
  voxtbl <- fread(x[1], na.strings = "NA" , skip = 5)
  voxtbl <- voxtbl[,1:4][,alt := k+zmin+0.5,][,dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  voxtbl <- voxtbl[alt>dtm]
  voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
  voxtbl <- voxtbl[,k1:= k1-0.5]
  voxtbl <- voxtbl[k1>2]
  ttl <- sum(voxtbl$PadBVTotal, na.rm = TRUE)
  ttl <- ttl/(pi*15*15)
  pf <- exp(-0.5*ttl)
  return(list(
    id_placette = as.factor(id_placette),
    meanang = as.factor(meanang),
    pfvox = pf))
}

voxtblsall <- data.frame()
for(group in unique(allvoxfiles$grp))
{
  voxfiles <- allvoxfiles[allvoxfiles$grp==group,]
  txt <- readLines(voxfiles$V1[[1]], n=3)[2:3]
  zmin <- as.numeric(unlist(strsplit(txt[[1]], "\\s+"))[4])
  lasnm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/", 
                  typefor, 
                  "/15m_rad/allpoints/",
                  group,
                  "@all",
                  ".las")
  ls <- readLASheader(lasnm)
  empty_raster <- raster(ncol=round(ls@PHB[["Max X"]]-ls@PHB[["Min X"]]), 
                         nrow=round(ls@PHB[["Max Y"]]-ls@PHB[["Min Y"]]),
                         xmn=ls@PHB[["Min X"]], xmx= ls@PHB[["Max X"]], 
                         ymn=ls@PHB[["Min Y"]], ymx= ls@PHB[["Max Y"]])
  dt <- alldtms[[paste0(group,"@all")]]
  dt <- resample(dt, empty_raster, method='bilinear')
  dt_mat <- as.matrix(dt)
  
  voxtbls <- apply(voxfiles, 1, func_computevoxpf, zmin=zmin, dt_mat=dt_mat)
  voxtbls <- rbindlist(voxtbls, use.names = TRUE)
  voxtblsall <- rbind(voxtblsall, voxtbls)
  
}

setkeyv(voxtblsall, keycols)
plotmets<- plotmets[voxtblsall]



#####
#####
# library(caret)
# placmes <- placmes_mi74
# setkey(placmes, id_placette)
# mets_for_model <- placmes[plotmets]
# 
# model <- summary(lm(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladvox),
#             data = mets_for_model))
# 
# mdl_mixte_all <- train(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
#                        data = mets_for_model,
#                        method = "lm",
#                        trControl = trainControl(method="LOOCV"))
# 


#####################################################################
#Several iterations of models with random flight line chosen per plot
#####################################################################

library(parallel)
library(foreach)
library(doParallel)
n <- 1000
plotmets <- plotmets[pflidr!=0]
setDT(placette.mes)
placmes <- placette.mes[which(Id_plac%in%plotmets$id_placette),c("Id_plac","G175"),]
names(placmes) <- c("id_placette", "G175")
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdl_lidr <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret"), .combine='rbind') %dopar% {
  mets_for_model <- plotmets[, .SD[sample(.N, min(1,.N))], by = id_placette]
  mets_for_model <- right_join(placmes, mets_for_model, by="id_placette")

  
  
  
  # model <- lm(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
  #             data = mets_for_model)
  # return(summary(model)$adj.r.squared)
  
  
  
  model <- train(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  return(list(summary(model)))

  # return(c(model$results$RMSE,
  #          model$results$MAE,
  #          model$results$Rsquared,
  #          c(summary(model)$coefficients[2:5,4])))
}
stopCluster(clus)


clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdl_pfsumvox <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret"), .combine='rbind') %dopar% {
  mets_for_model <- plotmets[, .SD[sample(.N, min(1,.N))], by = id_placette]
  mets_for_model <- right_join(placmes, mets_for_model, by="id_placette")
  
  
  
  # model <- lm(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
  #             data = mets_for_model)
  # return(summary(model)$adj.r.squared)
  
  
  
  model <- train(log(G175)~log(meanch)+log(varch)+log(pfsumvox)+log(cvladlidr),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  return(list(summary(model)))
  

  # return(c(model$results$RMSE,
  #          model$results$MAE,
  #          model$results$Rsquared,
  #          c(summary(model)$coefficients[2:5,4])))
}
stopCluster(clus)



clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdl_pfsumprof <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret"), .combine='rbind') %dopar% {
  mets_for_model <- plotmets[, .SD[sample(.N, min(1,.N))], by = id_placette]
  mets_for_model <- right_join(placmes, mets_for_model, by="id_placette")
  
  
  
  # model <- lm(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
  #             data = mets_for_model)
  # return(summary(model)$adj.r.squared)
  
  
  
  model <- train(log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladlidr),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  return(list(summary(model)))
  
  # return(c(model$results$RMSE,
  #          model$results$MAE,
  #          model$results$Rsquared,
  #          c(summary(model)$coefficients[2:5,4])))
}
stopCluster(clus)



clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdl_vox <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret"), .combine='rbind') %dopar% {
  mets_for_model <- plotmets[, .SD[sample(.N, min(1,.N))], by = id_placette]
  mets_for_model <- right_join(placmes, mets_for_model, by="id_placette")
  
  
  
  # model <- lm(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
  #             data = mets_for_model)
  # return(summary(model)$adj.r.squared)
  
  
  
  model <- train(log(G175)~log(meanch)+log(varch)+log(pfvox)+log(cvladvox),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  return(list(summary(model)))
  
  # return(c(model$results$RMSE,
  #          model$results$MAE,
  #          model$results$Rsquared,
  #          c(summary(model)$coefficients[2:5,4])))
}
stopCluster(clus)

mdls_fe_pfs <- list("lidr"=mdl_lidr,"pfsumvox"=mdl_pfsumvox,"pfsumprof"=mdl_pfsumprof)


rmse <- lapply(mdls_fe_pfs, function(x)
{
  lapply(x, function(y){
    RSS <- c(crossprod(y$residuals))
    MSE <- RSS / length(y$residuals)
    RMSE <- sqrt(MSE)
    return(RMSE)
  })
})
rmse <- reshape2::melt(rmse)
ggplot(data = rmse, aes(y=value, x=L1))+geom_boxplot()
ggplot(data = rmse, aes(x=value))+geom_histogram()+facet_wrap(~L1, ncol = 1)




arrangedata <- function(x1,x2,x3,x4)
{
  x1 <- setDT(as.data.frame(cbind(x1, rep("lidr",nrow(x1)))))
  x2 <- setDT(as.data.frame(cbind(x2, rep("cvladvox",nrow(x2)))))
  x3 <- setDT(as.data.frame(cbind(x3, rep("pfvox",nrow(x3)))))
  x4 <- setDT(as.data.frame(cbind(x4, rep("vox",nrow(x4)))))
  
  colnames(x1) <- c("rmse", "mae", "r2", "sigmean", "sigvar", "sigpf", "sigcvlad", "type")
  colnames(x2) <- c("rmse", "mae", "r2", "sigmean", "sigvar", "sigpf", "sigcvlad", "type")
  colnames(x3) <- c("rmse", "mae", "r2", "sigmean", "sigvar", "sigpf", "sigcvlad", "type")
  colnames(x4) <- c("rmse", "mae", "r2", "sigmean", "sigvar", "sigpf", "sigcvlad", "type")
  
  
  rmsevals <- rbind(x1[,c(1,8)],x2[,c(1,8)], x3[,c(1,8)], x4[,c(1,8)])
  maevals <- rbind(x1[,c(2,8)],x2[,c(2,8)], x3[,c(2,8)], x4[,c(2,8)])
  r2vals <- rbind(x1[,c(3,8)],x2[,c(3,8)], x3[,c(3,8)], x4[,c(3,8)])
  sigmean <- rbind(x1[,c(4,8)],x2[,c(4,8)], x3[,c(4,8)], x4[,c(4,8)])
  sigvar <- rbind(x1[,c(5,8)],x2[,c(5,8)], x3[,c(5,8)], x4[,c(5,8)])
  sigpf <- rbind(x1[,c(6,8)],x2[,c(6,8)], x3[,c(6,8)], x4[,c(6,8)])
  sigcvlad <- rbind(x1[,c(7,8)],x2[,c(7,8)], x3[,c(7,8)], x4[,c(7,8)])
  
  
  
  return(list(rmse=rmsevals, mae=maevals, r2=r2vals, sigmean=sigmean, sigvar=sigvar, sigpf=sigpf, sigcvlad=sigcvlad))
}


dta <- arrangedata(mdl_lidr, mdl_cvladvox, mdl_pfvox, mdl_vox)



ggplot(data=dta$sigpf, aes(x=type, y=as.numeric(sigpf)))+geom_boxplot()+theme_minimal()

ggplot(data=data_fe$r2, aes(x=as.numeric(r2)))+geom_histogram(bins=60)+facet_wrap(~type, ncol = 1)+theme_minimal()



x <- dta$r2 
x$type <- as.factor(x$type)













profile_plot <- function(profs)
{
  plt <- ggplot(data=profs[profs$k1>2,], aes(x=k1, y=m, colour = meanang, linetype = meanang))+
    scale_colour_brewer(palette="Paired")+
    theme_minimal()+
    geom_line(size=0.8)+
    coord_flip()+
    labs(title = paste0("Vertical profile for ", unique(profs$id_placette)),
         y ="PAD", 
         x = "Height above ground (m)", 
         colour="Mean scan angle", 
         linetype = "Mean scan angle")
  
  
  return(plt)
}

allplots <-  lapply(voxmergedall, profile_plot)


for(name in names(allplots))
{
  nm <- paste0("D:/1_Work/2_Ciron/voxelisation/Results/26jan/w_wt/",name,".png")
  ggsave(nm, allplots[[name]], device="png")
}


















voxtbl <- fread(allvoxfiles$V1[1], na.strings = "NA" , skip = 5)
voxtbl <- voxtbl[,1:4][,alt := k+zmin][,dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
#voxtbl <- na.omit(voxtbl)
voxtbl <- voxtbl[alt>dtm]
voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
voxtbl <- voxtbl[, .(m = mean(PadBVTotal, na.rm = TRUE)), by=list(k1)]
