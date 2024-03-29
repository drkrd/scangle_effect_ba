library(lidR)
library(data.table)


##here read only normalised points clouds with label format "plotid_n"
plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/coniferes/15m_rad/flightlines/")
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

plotmets_fl_co <- catalog_apply(plots, func_computeall)
plotmets_fl_co <- rbindlist(plotmets_fl_co)


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

plotmets_fl_co <- right_join(plotmets_fl_co, cvladvox, by=c("id_placette", "meanang"))


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
plotmets_fl_co <- right_join(plotmets_fl_co, pfvox, by=c("id_placette", "meanang"))

# plots_metrics <- na.omit(plots_metrics)
# plots_metrics_fe <- plots_metrics_fe[!which(plots_metrics$meanang%in%c(37.83,43.78,49.45,28.53)),]


#####################################################################
#Several iterations of models with random flight line chosen per plot
#####################################################################

library(parallel)
library(foreach)
library(doParallel)

clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdl_co74_cvladvox <- foreach(i = 1:5000, .packages=c("dplyr", "data.table", "caret"), .combine='rbind') %dopar% {
  mets_for_model <- plotmets_fl_co[, .SD[sample(.N, min(1,.N))], by = id_placette]
  mets_for_model <- right_join(placmes_co74, mets_for_model, by="id_placette")
  
  
  
  # model <- lm(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
  #             data = mets_for_model)
  # return(summary(model)$adj.r.squared)
  
  
  
  model <- train(log(G175)~log(meanch)+log(varch)+log(pflidr)+log(cvladvox),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  return(c(model$results$RMSE,
           model$results$MAE,
           model$results$Rsquared))
}
stopCluster(clus)

mdl_co74_cvladvox <- cbind(mdl_co74_cvladvox, rep("cvladvox",nrow(mdl_fe_cvvox)))
mdl_co74_lidr <- cbind(mdl_co74_lidr, rep("lidr",nrow(mdl_fe_lidr)))
mdl_co74_pfcvvox <- cbind(mdl_co74_vox, rep("vox",nrow(mdl_fe_pfcvvox)))
mdl_co74_pfvox <- cbind(mdl_co74_pfvox, rep("pfvox",nrow(mdl_fe_pfvox)))


arrangedata <- function(x1,x2,x3,x4)
{
  x1 <- setDT(as.data.frame(cbind(x1, rep("lidr",nrow(x1)))))
  x2 <- setDT(as.data.frame(cbind(x2, rep("cvladvox",nrow(x2)))))
  x3 <- setDT(as.data.frame(cbind(x3, rep("pfvox",nrow(x3)))))
  x4 <- setDT(as.data.frame(cbind(x4, rep("vox",nrow(x4)))))
  
  colnames(x1) <- c("rmse", "mae", "r2", "type")
  colnames(x2) <- c("rmse", "mae", "r2", "type")
  colnames(x3) <- c("rmse", "mae", "r2", "type")
  colnames(x4) <- c("rmse", "mae", "r2", "type")
  
  
  rmsevals <- rbind(x1[,c(1,4)],x2[,c(1,4)], x3[,c(1,4)], x4[,c(1,4)])
  maevals <- rbind(x1[,c(2,4)],x2[,c(2,4)], x3[,c(2,4)], x4[,c(2,4)])
  r2vals <- rbind(x1[,c(3,4)],x2[,c(3,4)], x3[,c(3,4)], x4[,c(3,4)])
  
  return(list(rmse=rmsevals, mae=maevals, r2=r2vals))
}


dta <- arrangedata(mdl_co74_lidr, mdl_co74_cvladvox, mdl_co74_pfvox, mdl_co74_vox)







colnames(mdl_fe_pfvox) <- c("rmse", "mae", "r2", "type")

mdl_fe_lidr <- setDT(as.data.frame(mdl_fe_lidr))
mdl_fe_pfcvvox <- setDT(as.data.frame(mdl_fe_pfcvvox))
mdl_fe_pfvox <- setDT(as.data.frame(mdl_fe_pfvox))
mdl_fe_cvvox <- setDT(as.data.frame(mdl_fe_cvvox))

rmsevals <- rbind(mdl_fe_lidr[,c(1,4)],mdl_fe_cvvox[,c(1,4)], mdl_fe_pfvox[,c(1,4)], mdl_fe_pfcvvox[,c(1,4)])
maevals <- rbind(mdl_fe_lidr[,c(2,4)],mdl_fe_cvvox[,c(2,4)], mdl_fe_pfvox[,c(2,4)], mdl_fe_pfcvvox[,c(2,4)])
r2vals <- rbind(mdl_fe_lidr[,c(3,4)],mdl_fe_cvvox[,c(3,4)], mdl_fe_pfvox[,c(3,4)], mdl_fe_pfcvvox[,c(3,4)])



ggplot(data=dta$rmse, aes(x=rmse))+geom_histogram()+facet_grid(~type)+theme_minimal()
ggplot(data=maevals, aes(y=mae, x=type))+geom_boxplot()+theme_minimal()
plot(plots_metrics$pflidr,plots_metrics$pfvox)


