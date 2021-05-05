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
plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/flightlines_123/")
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
plotmetsfl123 <- catalog_apply(plots, func_computeall)
plotmetsfl123 <- rbindlist(plotmetsfl123)
plotmetsfl123 <- plotmetsfl123[, c("m1","m2","m3"):=tstrsplit(meanang,"_",fixed=T),]
plotmetsfl123$m1 <- as.numeric(plotmetsfl123$m1)
plotmetsfl123$m2 <- as.numeric(plotmetsfl123$m2)
plotmetsfl123$m3 <- as.numeric(plotmetsfl123$m3)

plotmetsfl123 <- plotmetsfl123[, m2:= ifelse(is.na(m2), m1, m2),]
plotmetsfl123 <- plotmetsfl123[, m3:= ifelse(is.na(m3), m2, m3),]
plotmetsfl123 <- plotmetsfl123[, nc:= ifelse(m1!=m2 & m2!=m3, "3", 
                                             ifelse(m1!=m2 & m2==m3, "2", "1"))]


plotmetsfl123 <- plotmetsfl123[,class1:=ifelse(m1>=0&m1<10, "a",
                                               ifelse(m1>=10&m1<20, "b",
                                                      ifelse(m1>=20&m1<30, "c",
                                                             ifelse(m1>=30&m1<40,"d","e"))))]


plotmetsfl123 <- plotmetsfl123[,class2:=ifelse(m2>=0&m2<10, "a",
                                               ifelse(m2>=10&m2<20, "b",
                                                      ifelse(m2>=20&m2<30, "c",
                                                             ifelse(m2>=30&m2<40,"d","e"))))]


plotmetsfl123 <- plotmetsfl123[,class3:=ifelse(m3>=0&m3<10, "a",
                                               ifelse(m3>=10&m3<20, "b",
                                                      ifelse(m3>=20&m3<30, "c",
                                                             ifelse(m3>=30&m3<40,"d","e"))))]


plotmetsfl123 <- plotmetsfl123[class1 != "e"]
plotmetsfl123 <- plotmetsfl123[class2 != "e"]
plotmetsfl123 <- plotmetsfl123[class3 != "e"]


plotmetsfl123 <- plotmetsfl123[!m1 %in% c(37.24, 7.16)]
plotmetsfl123 <- plotmetsfl123[!m2 %in% c(37.24, 7.16)]
plotmetsfl123 <- plotmetsfl123[!m3 %in% c(37.24, 7.16)]


plotmetsfl3 <- plotmetsfl3[, cid := paste0(sort(.SD), collapse = ""), .SDcols = c("class1", "class2", "class3"),  by = 1:nrow(plotmetsfl3)]


plotmetsfl123 <- unique(plotmetsfl123[, prob := prop.table(table(nc))[nc], id_placette][])
plotmetsfl123 <- plotmetsfl123[, wt := (1/prob)/(1/sum(prob)), id_placette]
plotmetsfl123 <- plotmetsfl123[, wt := wt/sum(wt), id_placette]


fd_smry <- fread("D:/1_Work/2_Ciron/Data/var_dendro_ciron.csv", sep = ",", drop = "id_placette")
colnames(fd_smry)[1] <- "id_placette"
fd_smry[, id_placette := as.factor(id_placette)]
setkey(fd_smry, "id_placette")


n <- 10000
set.seed(123)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
mdls_vtotfl123 <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  mets_for_model <- plotmetsfl123[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fd_smry[mets_for_model]
  
  
  model <- train(log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr),
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  
  
  x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                "RMSE" = model$results$RMSE,
                                "MAE" = model$results$MAE),
             "coeffs" = list(model$finalModel$coefficients),
             "index" = which(plotmetsfl123$meanang %in% mets_for_model$meanang)))
  return(x)
}
stopCluster(clus)
registerDoSEQ()


func_entropy = function(cat.vect)
{
  px  = table(cat.vect)/length(cat.vect)
  lpx = log(px, base=2)
  ent = -sum(px*lpx)
  return(ent)
}

mdlmets_vtotfl123 <- lapply(mdls_vtotfl123, function(x)
{
  
  ind <- x[["index"]]
  dbase <- plotmetsfl123[ind]
  #dbase <- dbase[!id_placette%in%c(21, 22, 26)]
  ent <- func_entropy(dbase$nc)
  r2 <- x[["modelmets"]]$R2 
  rmse <- x[["modelmets"]]$RMSE
  mae <- x[["modelmets"]]$MAE
  type <- rep("3 flight lines")
  mx <- names(table(dbase$nc))[which(table(dbase$nc)==max(table(dbase$nc)))]
  if(mx==1) return(list("ent"=ent,"r2"=r2, "rmse"=rmse, "mae"=mae, "type"=type, "max"=mx))
  # names(tbl[which(tbl==max(tbl))])
})
mdlmets_vtotfl123 <- rbindlist(mdlmets_vtotfl123)


ggplot(data = mdlmets_gfl123, aes(x=r2, y=ent))+geom_point()



ggplot(data=mets_mdl_vtig1, aes(x=value, fill=flightlines))+
  geom_histogram(bins=60)+
  facet_wrap(~variable, scales = "free")+
  labs(title = "Stem volume")+
  theme_tufte(base_size = 13)+
  scale_fill_manual(values = cbp1)+
  geom_vline(data=filter(mets_mdl_vtig1, variable=="r2"), aes(xintercept=0.84), colour="red") + 
  geom_vline(data=filter(mets_mdl_vtig1, variable=="rmse"), aes(xintercept=0.16), colour="blue") + 
  geom_vline(data=filter(mets_mdl_vtig1, variable=="mae"), aes(xintercept=0.128), colour="green")










