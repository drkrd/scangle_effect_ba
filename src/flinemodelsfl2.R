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

typefor <- " "
typepc <- "allpoints"
# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
# plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/flightlines_2/")
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
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
plotmetsfl2 <- catalog_apply(plots, func_computeall)
plotmetsfl2 <- rbindlist(plotmetsfl2)
plotmetsfl2 <- plotmetsfl2[, c("m1","m2"):=tstrsplit(meanang,"_",fixed=T),]
plotmetsfl2$m1 <- as.numeric(plotmetsfl2$m1)
plotmetsfl2$m2 <- as.numeric(plotmetsfl2$m2)
plotmetsfl2 <- plotmetsfl2[, m2:= ifelse(is.na(m2), m1, m2),]

plotmetsfl2 <- plotmetsfl2[,class1:=ifelse(m1>=0&m1<10, "a",
                                           ifelse(m1>=10&m1<20, "b",
                                                  ifelse(m1>=20&m1<30, "c",
                                                         ifelse(m1>=30&m1<40,"d","e"))))]
plotmetsfl2 <- plotmetsfl2[,class2:=ifelse(m2>=0&m2<10, "a",
                                           ifelse(m2>=10&m2<20, "b",
                                                  ifelse(m2>=20&m2<30, "c",
                                                         ifelse(m2>=30&m2<40,"d","e"))))]
plotmetsfl2 <- plotmetsfl2[class1 != "e"]
plotmetsfl2 <- plotmetsfl2[class2 != "e"]

plotmetsfl2 <- plotmetsfl2[!m1 %in% c(37.24, 7.16)]
plotmetsfl2 <- plotmetsfl2[!m2 %in% c(37.24, 7.16)]

plotmetsfl2 <- plotmetsfl2[, cid := paste0(sort(.SD), collapse = ""), .SDcols = c("class1", "class2"),  by = 1:nrow(plotmetsfl2)]



plotmetsfl2 <- unique(plotmetsfl2[, prob := prop.table(table(cid))[cid], id_placette][])
plotmetsfl2 <- plotmetsfl2[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl2 <- plotmetsfl2[, wt := wt/sum(wt), id_placette]


fd_smry <- fread("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",", drop = "id_placette")
colnames(fd_smry)[1] <- "id_placette"
fd_smry[, id_placette := as.factor(id_placette)]
setkey(fd_smry, "id_placette")



n <- 10000
set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_gfl2 <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl2[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fd_smry[mets_for_model]
  
  
  model <- train(log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  
  x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                "RMSE" = model$results$RMSE,
                                "MAE" = model$results$MAE),
             "coeffs" = list(model$finalModel$coefficients),
             "index" = which(plotmetsfl2$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()




mdlmets_gfl2 <- lapply(mdls_gfl2, function(x)
{
  
  ind <- x[["index"]]
  dbase <- plotmetsfl2[ind]
  dbase <- dbase[!id_placette%in%c(21, 22, 26)]
  ent <- func_entropy(dbase$cid)
  r2 <- x[["modelmets"]]$R2 
  rmse <- x[["modelmets"]]$RMSE
  mae <- x[["modelmets"]]$MAE
  type <- rep("2 flight lines")
  # names(tbl[which(tbl==max(tbl))])
  return(list("ent"=ent,"r2"=r2, "rmse"=rmse, "mae"=mae, "type"=type))
})
mdlmets_gfl2 <- rbindlist(mdlmets_gfl2)




ggplot(data=mets_mdl_vtig1, aes(x=value, fill=flightlines))+
  geom_histogram(bins=60)+
  facet_wrap(~variable, scales = "free")+
  labs(title = "Stem volume")+
  theme_tufte(base_size = 13)+
  scale_fill_manual(values = cbp1)+
  geom_vline(data=filter(mets_mdl_vtig1, variable=="r2"), aes(xintercept=0.84), colour="red") + 
  geom_vline(data=filter(mets_mdl_vtig1, variable=="rmse"), aes(xintercept=0.16), colour="blue") + 
  geom_vline(data=filter(mets_mdl_vtig1, variable=="mae"), aes(xintercept=0.128), colour="green")










