library(lidR)
plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/all/15m_rad/")
opt_independent_files(plots) <- TRUE

func_computeall <- function(chunk)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  id_placette <- sub("\\_.*", "", basename(tools::file_path_sans_ext(chunk@files)))
  meanang <- sub(".*\\_", "", basename(tools::file_path_sans_ext(chunk@files)))
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber)
  return(list( 
    id_placette = id_placette,
    meanang = meanang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}

plots_metrics <- catalog_apply(plots.n, func_computeall)
plots_metrics <- rbindlist(plots_metrics)


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
cvladvox <- cvladvox[as.numeric(sub(".*\\_","",id))%%1!=0][,meanang:=sub(".*\\_","",id)][,`:=`(x=NULL, id=NULL)]
cvladvox <- cvladvox[-93]

plots_metrics <- right_join(plots_metrics, cvladvox, by=c("id_placette", "meanang"))


#####################################################
#Compute pf from voxels by inverting beer-lambert law
#####################################################

pfvox <- apply(allvoxfiles, 1, function(x){
  id_placette <- sub("\\_.*", "", basename(tools::file_path_sans_ext(x[1])))
  meanang <- sub(".*\\_", "", basename(tools::file_path_sans_ext(x[1])))
  voxtbl <- fread(x[1], na.strings = "NA" , skip = 5)
  voxtbl <- voxtbl[,1:4][,alt := k+zmin][,dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  #voxtbl <- na.omit(voxtbl)
  voxtbl <- voxtbl[alt>dtm]
  voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
  voxtbl <- voxtbl[k1>2]
  ttl <- sum(voxtbl$PadBVTotal, na.rm = TRUE)
  ttl <- ttl/(pi*17.5*17.5)
  pf <- exp(-0.5*ttl)
  return(list(
    id_placette = id_placette,
    meanang = meanang,
    pfvox = pf))
})

pfvox <- rbindlist(pfvox1)
pfvox <- pfvox[as.numeric(meanang)%%1!=0]
plots_metrics <- right_join(plots_metrics, pfvox1, by=c("id_placette", "meanang"))

plots_metrics <- na.omit(plots_metrics)
plots_metrics <- plots_metrics[!which(plots_metrics$meanang%in%c(37.83,43.78,49.45,28.53)),]


#####################################################################
#Several iterations of models with random flight line chosen per plot
#####################################################################

lls <- vector()
for(i in 1:5000)
{
  mets_for_model <- plots_metrics[, .SD[sample(.N, min(1,.N))], by = id_placette]
  if(isFALSE(list(mets_for_model$meanang)%in%list_meanangs))
  {
    list_meanangs[[i]] <- mets_for_model$meanang
    mets_for_model <- right_join(fd_ba_smry, mets_for_model, by="id_placette")
    model <- lm(data = mets_for_model, 
                formula = log(sum_ba_hec)~log(meanch)+log(varch)+log(pflidr)+log(cvladvox))
    lls <- c(lls, summary(model)$adj.r.squared)
  }
}
beep()


rvals <- cbind(rvals, lllv=lls)
names(rvals) <- c("id","lidr","only_pfvox","pfcvladvox","onlycvlad")
rvals2 <- melt(rvals, id.vars = "id")

rvals <- as.data.frame(rvals) 





library(parallel)
library(foreach)
library(doParallel)

clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)


out <- foreach(i = 1:1000, .packages=c("dplyr", "data.table", "caret"), .combine='c') %dopar% {
  mets_for_model <- plots_metrics[, .SD[sample(.N, min(1,.N))], by = id_placette]
  mets_for_model <- right_join(fd_ba_smry, mets_for_model, by="id_placette")
  
  model <- train(data = mets_for_model, 
                 formula = log(sum_ba_hec)~log(meanch)+log(varch)+log(pflidr)+log(cvladvox),
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  list(model$results$RMSE,
       model$results$Rsquared,
       model$results$MAE)
}
stopCluster(clus)






ggplot(data=rvals2, aes(x=value))+geom_histogram()+facet_grid(~variable)+theme_minimal()
ggplot(data=rvals2, aes(y=value, x=variable))+geom_boxplot()+theme_minimal()
plot(plots_metrics$pflidr,plots_metrics$pfvox)
abline(a=0,b=1)
       
