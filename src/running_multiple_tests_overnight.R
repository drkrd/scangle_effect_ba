l1 <- c("pmetsfl2.mix.clab", 
        "pmetsfl2.mix.clac", 
        "pmetsfl2.mix.clbc",
        "pmetsfl2.feu.all", 
        "pmetsfl2.feu.clab", 
        "pmetsfl2.feu.clac", 
        "pmetsfl2.feu.clbc",
        "pmetsfl2.con.all", 
        "pmetsfl2.con.clab", 
        "pmetsfl2.con.clac", 
        "pmetsfl2.con.clbc")
l2 <- c("smplstfl2.mix.clab", 
        "smplstfl2.mix.clac", 
        "smplstfl2.mix.clbc",
        "smplstfl2.feu.all", 
        "smplstfl2.feu.clab", 
        "smplstfl2.feu.clac", 
        "smplstfl2.feu.clbc",
        "smplstfl2.con.all", 
        "smplstfl2.con.clab", 
        "smplstfl2.con.clac", 
        "smplstfl2.con.clbc")
l3 <- c("bauges.db.mix", 
        "bauges.db.mix", 
        "bauges.db.mix",
        "bauges.db.feu", 
        "bauges.db.feu", 
        "bauges.db.feu", 
        "bauges.db.feu",
        "bauges.db.con", 
        "bauges.db.con", 
        "bauges.db.con", 
        "bauges.db.con")
l4 <- c("baugesfl2.mix.clab1", 
        "baugesfl2.mix.clac1", 
        "baugesfl2.mix.clbc1",
        "baugesfl2.feu.all1", 
        "baugesfl2.feu.clab1", 
        "baugesfl2.feu.clac1", 
        "baugesfl2.feu.clbc1",
        "baugesfl2.con.all1", 
        "baugesfl2.con.clab1", 
        "baugesfl2.con.clac1", 
        "baugesfl2.con.clbc1")


df <- data.frame(l1, l2, l3, l4)




for(i in 1:nrow(df))
{
  dbase <- get(df$l1[i])
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  x <- foreach(i = 1:6500, .packages=c("dplyr", "data.table", "caret", "sampling")) %dopar% {
    set.seed(i)
    dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
    return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
  }
  stopCluster(clus)
  registerDoSEQ()
  assign(df$l2[i], unique(as.data.table(t(matrix(unlist(x), nrow = length(unique(dbase$id_placette)))))))
  
  dbase <- get(df$l1[i])
  fds <- get(df$l3[i])
  idx.lst <- get(df$l2[i])
  
  f1l <- log(G75)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f1v <- log(G75)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f2l <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f2v <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f3l <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f3v <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  
  
  ##Simulations
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  n <- 5000
  assign(df$l4[i],foreach(x = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    idx <- as.vector(unlist(idx.lst[x]))
    mets_for_model <- dbase[idx]
    dbase$id_placette <- as.factor(dbase$id_placette)
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    mets1 <- mets_for_model[,c("G75","meanch", "varch", 'pflidr', "cvladlidr")]
    
    # y1 <- log(mets1$G_m2_ha)
    # x1 <- log(mets1$meanch)
    # x2 <- log(mets1$varch)
    # x3 <- log(mets1$pflidr)
    # x4 <- log(mets1$cvladlidr)
    
    # func_linreg <- function(b0 ,b1, b2, b3, b4, sig)
    # {
    #   ypred <- b0+b1*x1+b2*x2+b3*x3+b4*x4
    #   -sum(dnorm(y1, mean = ypred, sd = sig, log=TRUE))
    # }
    # 
    
    m1l <- train(f1l,
                 data = mets1,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    
    G.l.coeff <- m1l$finalModel$coefficients
    G.l.pred <- m1l$pred$pred
    G.l.obs <- m1l$pred$obs
    
    
    mets2 <- mets_for_model[,c("G75","meanch", "varch", 'pfsumprof', "cvladvox")]
    m1v <- train(f1v,
                 data = mets2,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    G.v.coeff <- m1v$finalModel$coefficients
    G.v.pred <- m1v$pred$pred
    G.v.obs <- m1v$pred$obs
    ##########################################################
    
    ##########################################################
    mets3 <- mets_for_model[,c("volume_total", "meanch", "varch", 'pflidr', "cvladlidr")]
    m2l <- train(f2l,
                 data = mets3,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    vtot.l.coeff <- m2l$finalModel$coefficients
    vtot.l.pred <- m2l$pred$pred
    vtot.l.obs <- m2l$pred$obs
    
    mets4 <- mets_for_model[,c("volume_total", "meanch", "varch", 'pfsumprof', "cvladvox")]
    m2v <- train(f2v,
                 data = mets4,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    vtot.v.coeff <- m2v$finalModel$coefficients
    vtot.v.pred <- m2v$pred$pred
    vtot.v.obs <- m2v$pred$obs
    ##########################################################
    
    ###########################################################
    mets5 <- mets_for_model[,c("volume_tige", "meanch", "varch", 'pflidr', "cvladlidr")]
    m3l <- train(f3l,
                 data = mets5,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    vtig.l.coeff <- m3l$finalModel$coefficients
    vtig.l.pred <- m3l$pred$pred
    vtig.l.obs <- m3l$pred$obs
    
    mets6 <- mets_for_model[,c("volume_tige", "meanch", "varch", 'pfsumprof', "cvladvox")]
    m3v <- train(f3v,
                 data = mets6,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    vtig.v.coeff <- m3v$finalModel$coefficients
    vtig.v.pred <- m3v$pred$pred
    vtig.v.obs <- m3v$pred$obs
    ############################################################
    ############################################################
    
    x <- (list("G.lidr.pred" = G.l.pred,
               "G.lidr.obs" = G.l.obs,
               "G.vox.pred" = G.v.pred,
               "G.vox.obs" = G.v.obs,
               "vtot.lidr.pred" = vtot.l.pred,
               "vtot.lidr.obs" = vtot.l.obs,
               "vtot.vox.pred" = vtot.v.pred,
               "vtot.vox.obs" = vtot.v.obs,
               "vtig.lidr.pred" = vtig.l.pred,
               "vtig.lidr.obs" = vtig.l.obs,
               "vtig.vox.pred" = vtig.v.pred,
               "vtig.vox.obs" = vtig.v.obs,
               "G.l.coeff" = G.l.coeff,
               "G.v.coeff" = G.v.coeff,
               "vtot.l.coeff" = vtot.l.coeff,
               "vtot.v.coeff" = vtot.v.coeff,
               "vtig.l.coeff" = vtig.l.coeff,
               "vtig.v.coeff" = vtig.v.coeff
    ))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  
}







ciron.db <- ciron.db[!id_placetten%in%c("14","144")]

names(ciron.db) <- c("G75", "volume_total", "volume_tige", "id_placette")

l1 <- c("pmetsfl1.all",
        "pmetsfl1.cla", 
        "pmetsfl1.clb", 
        "pmetsfl1.clc",
        "pmetsfl2.all", 
        "pmetsfl2.clab", 
        "pmetsfl2.clac", 
        "pmetsfl2.clbc", 
        "pmetsfl3.all")
l2 <- c("smplstfl1.all", 
        "smplstfl1.cla", 
        "smplstfl1.clb",
        "smplstfl1.clc", 
        "smplstfl2.all",
        "smplstfl2.clab", 
        "smplstfl2.clac", 
        "smplstfl2.clbc",
        "smplstfl3.all")
l3 <- c("ciron.db", 
        "ciron.db", 
        "ciron.db",
        "ciron.db", 
        "ciron.db", 
        "ciron.db", 
        "ciron.db",
        "ciron.db", 
        "ciron.db")
l4 <- c("cironfl1.all",
        "cironfl1.cla", 
        "cironfl1.clb", 
        "cironfl1.clc",
        "cironfl2.all", 
        "cironfl2.clab", 
        "cironfl2.clac", 
        "cironfl2.clbc",
        "cironfl3.all")



df <- data.frame(l1, l2, l3, l4)



time.log <- NULL
for(i in 1:9)
{
  start <- Sys.time()
  dbase <- get(df$l1[i])
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  x <- foreach(ii = 1:10000, .packages=c("dplyr", "data.table", "caret", "sampling")) %dopar% {
    set.seed(ii)
    dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
    return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
  }
  stopCluster(clus)
  registerDoSEQ()
  assign(df$l2[i], unique(as.data.table(t(matrix(unlist(x), nrow = length(unique(dbase$id_placette)))))))
  
  dbase <- get(df$l1[i])
  fds <- get(df$l3[i])
  idx.lst <- get(df$l2[i])
  
  f1l <- log(G75)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f1v <- log(G75)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f2l <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f2v <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f3l <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f3v <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  
  
  ##Simulations
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  n <- 5000
  assign(df$l4[i],foreach(x = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    idx <- as.vector(unlist(idx.lst[x]))
    mets_for_model <- dbase[idx]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    mets1 <- mets_for_model[,c("G75","meanch", "varch", 'pflidr', "cvladlidr")]
    
    # y1 <- log(mets1$G_m2_ha)
    # x1 <- log(mets1$meanch)
    # x2 <- log(mets1$varch)
    # x3 <- log(mets1$pflidr)
    # x4 <- log(mets1$cvladlidr)
    
    # func_linreg <- function(b0 ,b1, b2, b3, b4, sig)
    # {
    #   ypred <- b0+b1*x1+b2*x2+b3*x3+b4*x4
    #   -sum(dnorm(y1, mean = ypred, sd = sig, log=TRUE))
    # }
    # 
    
    m1l <- train(f1l,
                 data = mets1,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    
    G.l.coeff <- m1l$finalModel$coefficients
    G.l.pred <- m1l$pred$pred
    G.l.obs <- m1l$pred$obs
    
    
    mets2 <- mets_for_model[,c("G75","meanch", "varch", 'pfsumprof', "cvladvox")]
    m1v <- train(f1v,
                 data = mets2,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    G.v.coeff <- m1v$finalModel$coefficients
    G.v.pred <- m1v$pred$pred
    G.v.obs <- m1v$pred$obs
    ##########################################################
    
    ##########################################################
    mets3 <- mets_for_model[,c("volume_total", "meanch", "varch", 'pflidr', "cvladlidr")]
    m2l <- train(f2l,
                 data = mets3,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    vtot.l.coeff <- m2l$finalModel$coefficients
    vtot.l.pred <- m2l$pred$pred
    vtot.l.obs <- m2l$pred$obs
    
    mets4 <- mets_for_model[,c("volume_total", "meanch", "varch", 'pfsumprof', "cvladvox")]
    m2v <- train(f2v,
                 data = mets4,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    vtot.v.coeff <- m2v$finalModel$coefficients
    vtot.v.pred <- m2v$pred$pred
    vtot.v.obs <- m2v$pred$obs
    ##########################################################
    
    ###########################################################
    mets5 <- mets_for_model[,c("volume_tige", "meanch", "varch", 'pflidr', "cvladlidr")]
    m3l <- train(f3l,
                 data = mets5,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    vtig.l.coeff <- m3l$finalModel$coefficients
    vtig.l.pred <- m3l$pred$pred
    vtig.l.obs <- m3l$pred$obs
    
    mets6 <- mets_for_model[,c("volume_tige", "meanch", "varch", 'pfsumprof', "cvladvox")]
    m3v <- train(f3v,
                 data = mets6,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
    vtig.v.coeff <- m3v$finalModel$coefficients
    vtig.v.pred <- m3v$pred$pred
    vtig.v.obs <- m3v$pred$obs
    ############################################################
    ############################################################
    
    lst <- (list("G.lidr.pred" = G.l.pred,
               "G.lidr.obs" = G.l.obs,
               "G.vox.pred" = G.v.pred,
               "G.vox.obs" = G.v.obs,
               "vtot.lidr.pred" = vtot.l.pred,
               "vtot.lidr.obs" = vtot.l.obs,
               "vtot.vox.pred" = vtot.v.pred,
               "vtot.vox.obs" = vtot.v.obs,
               "vtig.lidr.pred" = vtig.l.pred,
               "vtig.lidr.obs" = vtig.l.obs,
               "vtig.vox.pred" = vtig.v.pred,
               "vtig.vox.obs" = vtig.v.obs,
               "G.l.coeff" = G.l.coeff,
               "G.v.coeff" = G.v.coeff,
               "vtot.l.coeff" = vtot.l.coeff,
               "vtot.v.coeff" = vtot.v.coeff,
               "vtig.l.coeff" = vtig.l.coeff,
               "vtig.v.coeff" = vtig.v.coeff
    ))
    return(lst)
    if(x%%1000==0)print(x)
  })
  stopCluster(clus)
  registerDoSEQ()
print("Done")  
end <- Sys.time()
time.log <- rbind(rbind(time.log, c(n, (end-start))))
}

