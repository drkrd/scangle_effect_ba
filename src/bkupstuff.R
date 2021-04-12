library(lidR)
library(data.table)


##here read only normalised points clouds with label format "plotid_n"
plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/mixte/15m_rad/allpoints/")
####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS
opt_independent_files(plots) <- TRUE 

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

plotmets_all_mi <- catalog_apply(plots, func_computeall)
plotmets_all_mi <- rbindlist(plotmets_all_mi)


##########################
#Compute cvlad from voxels
##########################

cvladvox <- lapply(voxmergedall, function(x)
{
  x <- x[k1>2]
  return(x[,.(cvs=as.double(cv(m, na.rm = TRUE))), by=meanang])
})

cvladvox <- reshape2::melt(cvladvox)
names(cvladvox) <- c("meanang","x","cvladvox","id_placette") 
setDT(cvladvox)
cvladvox$id_placette <- as.factor(cvladvox$id_placette)


# cvladvox <- cvladvox %>% 
#   mutate(meanang=round(as.numeric(sub(".*\\_", "", basename(tools::file_path_sans_ext(id)))),2))


cvladvox <- cvladvox[,`:=`(x=NULL)]
# cvladvox <- cvladvox[-124,]

plotmets_all_mi <- right_join(plotmets_all_mi, cvladvox, by=c("id_placette", "meanang"))


#####################################################
#Compute pf from voxels by inverting beer-lambert law
#####################################################

pfvox <- apply(allvoxfiles, 1, function(x){
  id_placette <- x[2]
  meanang <- x[3]
  voxtbl <- fread(x[1], na.strings = "NA" , skip = 5)
  voxtbl <- voxtbl[,1:4][,alt := k+zmin+0.5,][,dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  #voxtbl <- na.omit(voxtbl)
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
})

pfvox <- rbindlist(pfvox)
# pfvox <- pfvox[-124,]
plotmets_all_mi <- right_join(plotmets_all_mi, pfvox, by=c("id_placette", "meanang"))

# plots_metrics <- na.omit(plots_metrics)
# plots_metrics_fe <- plots_metrics_fe[!which(plots_metrics$meanang%in%c(37.83,43.78,49.45,28.53)),]


#####################################################################
#Several iterations of models with random flight line chosen per plot
#####################################################################

fieldmets_mi <- placette.mes[which(placette.mes$Id_plac %in% plotmets_all_mi$id_placette),c("Id_plac","G175")]
colnames(fieldmets_mi) <- c("id_placette", "G175")

library(parallel)
library(foreach)
library(doParallel)

clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdl_all_mi_cvladvox <- foreach(i = 1:1, .packages=c("dplyr", "data.table", "caret"), .combine='rbind') %dopar% {
  mets_for_model <- plotmets_all_mi[, .SD[sample(.N, min(1,.N))], by = id_placette]
  mets_for_model <- right_join(fieldmets_mi, mets_for_model, by="id_placette")
  
  
  
  # model <- lm(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
  #             data = mets_for_model)
  # return(summary(model)$adj.r.squared)
  
  
  
  model <- train(log(G175)~log(meanch)+log(varch)+log(1-pflidr)+log(cvladvox),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  return(c(model$results$RMSE,
           model$results$MAE,
           model$results$Rsquared))
}
stopCluster(clus)











mdl_all_mi <- rbind(mdl_all_mi_lidr, mdl_all_mi_vox, mdl_all_mi_pfvox, mdl_all_mi_cvladvox)
mdl_all_mi <- setDT(as.data.frame(mdl_all_mi))
mdl_all_mi <- cbind(mdl_all_mi, c("lidr", "vox", "pfvox", "cvladvox"))
colnames(mdl_all_mi) <- c("rmse", "mae", "r2", "type")
rownames(mdl_all_mi) <- NULL
mdl_all_mi$type <- as.factor(mdl_all_mi$type)
mdl_all_mi <- melt(mdl_all_mi, id.vars = "type")
mdl_all_mi <- setDT(mdl_all_mi)



mdl_all_fe_cvladvox <- cbind(mdl_all_fe_cvladvox, rep("cvvox",nrow(mdl_all_fe_cvladvox)))
mdl_all_fe_lidr <- cbind(mdl_all_fe_lidr, rep("lidr",nrow(mdl_all_fe_lidr)))
mdl_all_fe_vox <- cbind(mdl_all_fe_vox, rep("pfcvvox",nrow(mdl_all_fe_vox)))
mdl_all_fe_pfvox <- cbind(mdl_all_fe_pfvox, rep("pfvox",nrow(mdl_all_fe_pfvox)))


colnames(mdl_fe_pfvox) <- c("rmse", "mae", "r2", "type")

mdl_fe_lidr <- setDT(as.data.frame(mdl_fe_lidr))
mdl_fe_pfcvvox <- setDT(as.data.frame(mdl_fe_pfcvvox))
mdl_fe_pfvox <- setDT(as.data.frame(mdl_fe_pfvox))
mdl_fe_cvvox <- setDT(as.data.frame(mdl_fe_cvvox))

rmsevals <- rbind(mdl_fe_lidr[,c(1,4)],mdl_fe_cvvox[,c(1,4)], mdl_fe_pfvox[,c(1,4)], mdl_fe_pfcvvox[,c(1,4)])
maevals <- rbind(mdl_fe_lidr[,c(2,4)],mdl_fe_cvvox[,c(2,4)], mdl_fe_pfvox[,c(2,4)], mdl_fe_pfcvvox[,c(2,4)])
r2vals <- rbind(mdl_fe_lidr[,c(3,4)],mdl_fe_cvvox[,c(3,4)], mdl_fe_pfvox[,c(3,4)], mdl_fe_pfcvvox[,c(3,4)])



ggplot(data=r2vals, aes(x=r2))+geom_histogram()+facet_grid(~type)+theme_minimal()
ggplot(data=maevals, aes(y=mae, x=type))+geom_boxplot()+theme_minimal()
plot(plots_metrics$pflidr,plots_metrics$pfvox)










trn <- data.frame()
x <- list.files("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/allpointsallplots/15m_rad/")
for(i in x)
{
  nm <- file_path_sans_ext(basename(i))
  ls <- readLAS(i)
  ls <- grid_terrain(ls, res=0.1, tin())
  crs(ls) <- CRS("+init=epsg:2154")
  ls <- as.data.frame(terrain(ls, opt=c("slope", "aspect"), unit="degrees" ))
  trn <- rbind(trn, c(mean(ls$slope, na.rm=TRUE), var(ls$slope, na.rm=TRUE), 
                      mean(ls$aspect, na.rm=TRUE), var(ls$aspect, na.rm=TRUE)))
}
