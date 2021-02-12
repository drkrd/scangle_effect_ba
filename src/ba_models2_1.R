library(lidR)
library(data.table)


plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/feuillus/15m_rad/flightlines/")
opt_independent_files(plots) <- TRUE ####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS

func_computeall <- function(chunk)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  id_plac <- sub("_n.*", "", basename(tools::file_path_sans_ext(chunk@files)))
  mang <- round(as.numeric(sub(".*\\_", "", basename(tools::file_path_sans_ext(chunk@files)))),2)
  
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

plots_metrics_fe <- catalog_apply(plots, func_computeall)
plots_metrics_fe <- rbindlist(plots_metrics_fe)


##########################
#Compute cvlad from voxels
##########################

cvladvox <- lapply(voxmergedall, function(x)
{
  x <- x[k1>2]
  return(x[,.(cvs=as.double(cv(m, na.rm = TRUE))), by=id])
})

cvladvox <- melt(cvladvox)
names(cvladvox) <- c("id","x","cvladvox","id_placette") 
setDT(cvladvox)


cvladvox <- cvladvox %>% 
  mutate(meanang=round(as.numeric(sub(".*\\_", "", basename(tools::file_path_sans_ext(id)))),2))


cvladvox <- cvladvox[,`:=`(x=NULL, id=NULL)]
cvladvox <- cvladvox[-124,]

plots_metrics_fe <- right_join(plots_metrics_fe, cvladvox, by=c("id_placette", "meanang"))


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
    id_placette = id_placette,
    meanang = meanang,
    pfvox = pf))
})

pfvox <- rbindlist(pfvox)
pfvox <- pfvox[-124,]
plots_metrics_fe <- right_join(plots_metrics_fe, pfvox, by=c("id_placette", "meanang"))

plots_metrics <- na.omit(plots_metrics)
plots_metrics_fe <- plots_metrics_fe[!which(plots_metrics$meanang%in%c(37.83,43.78,49.45,28.53)),]


#####################################################################
#Several iterations of models with random flight line chosen per plot
#####################################################################

library(parallel)
library(foreach)
library(doParallel)

clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdl_fe_pfcvvox <- foreach(i = 1:1000, .packages=c("dplyr", "data.table", "caret"), .combine='rbind') %dopar% {
  mets_for_model <- plots_metrics_fe[, .SD[sample(.N, min(1,.N))], by = id_placette]
  mets_for_model <- right_join(placette_mes_fe, mets_for_model, by="id_placette")
  
  
  
  # model <- lm(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
  #             data = mets_for_model)
  # return(summary(model)$adj.r.squared)
  
  
  
  model <- train(log(G175)~log(meanch)+log(varch)+log(pfvox)+log(cvladvox),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  return(c(model$results$RMSE,
           model$results$MAE,
           model$results$Rsquared))
}
stopCluster(clus)

mdl_fe_cvvox <- cbind(mdl_fe_cvvox, rep("cvvox",nrow(mdl_fe_cvvox)))
mdl_fe_lidr <- cbind(mdl_fe_lidr, rep("lidr",nrow(mdl_fe_lidr)))
mdl_fe_pfcvvox <- cbind(mdl_fe_pfcvvox, rep("pfcvvox",nrow(mdl_fe_pfcvvox)))
mdl_fe_pfvox <- cbind(mdl_fe_pfvox, rep("pfvox",nrow(mdl_fe_pfvox)))


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

       
