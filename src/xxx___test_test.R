f2l <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
dbase <- pmetsfl1.clc
fds <- fd_smry
idx <- as.vector(unlist(smplst.clc[3709]))


test <- ciron.mdlmets[which(Forest_attr=="Total volume" & exp=="mostly C" & Metrics=="vox" & variable=="R2")] 
test1 <- test[value>0.61]
test2 <- test[value<0.61]
mmean1 <- c()
mvar1 <- c()
mpfv1 <- c()
mcvv1 <- c()
for(id in test1$id)
{
  idx <- as.vector(unlist(smplst.clc[id]))
  dbase <- pmetsfl1.clc
  mets_for_model <- dbase[idx]
  mmean1 <-c(mmean1,var(log(mets_for_model$meanch)))
  mvar1 <- c(mvar1,var(log(mets_for_model$varch)))
  mpfv1 <- c(mpfv1,var(log(mets_for_model$pfsumprof)))
  mcvv1 <- c(mcvv1,var(log(mets_for_model$cvladvox)))
}


mmean2 <- c()
mvar2 <- c()
mpfv2 <- c()
mcvv2 <- c()
for(id in test2$id)
{
  idx <- as.vector(unlist(smplst.clc[id]))
  dbase <- pmetsfl1.clc
  mets_for_model <- dbase[idx]
  mmean2 <-c(mmean2,var(log(mets_for_model$meanch)))
  mvar2 <- c(mvar2,var(log(mets_for_model$varch)))
  mpfv2 <- c(mpfv2,var(log(mets_for_model$pfsumprof)))
  mcvv2 <- c(mcvv2,var(log(mets_for_model$cvladvox)))
}

means <- rbind(cbind(mmean1, rep("high", length(mmean1))), cbind(mmean2, rep("low", length(mmean2))))
means <- as.data.table(means)
means$mmean1 <- as.numeric(means$mmean1)
means <- cbind(means, rep("means", length(means)))

vars <- rbind(cbind(mvar1, rep("high", length(mvar1))), cbind(mvar2, rep("low", length(mvar2))))
vars <- as.data.table(vars)
vars$mvar1 <- as.numeric(vars$mvar1)
vars <- cbind(vars, rep("vars", length(vars)))

pfvs <- rbind(cbind(mpfv1, rep("high", length(mpfv1))), cbind(mpfv2, rep("low", length(mpfv2))))
pfvs <- as.data.table(pfvs)
pfvs$mpfv1 <- as.numeric(pfvs$mpfv1)
pfvs <- cbind(pfvs, rep("pfvs", length(pfvs)))

cvvs <- rbind(cbind(mcvv1, rep("high", length(mcvv1))), cbind(mcvv2, rep("low", length(mcvv2))))
cvvs <- as.data.table(cvvs)
cvvs$mcvv1 <- as.numeric(cvvs$mcvv1)
cvvs <- cbind(cvvs, rep("cvvs", length(cvvs)))

ggplot(data = means, aes(y=mmean1, x=V2, fill=V2))+geom_boxplot()

xyz <- rbind(means, vars, pfvs, cvvs, use.names=F)
colnames(xyz) <- c("val", "v1", "v2")

ggplot(data = xyz, aes(y=val, x=v1, fill=v1))+geom_boxplot()+facet_wrap(~v2, scales = "free")



s1 <- test1[which(value==max(value))]
s2 <- test2[which(value==min(value))]
idx1 <- as.vector(unlist(smplst.clc[s1$id]))
idx2 <- as.vector(unlist(smplst.clc[s2$id]))
dbase <- pmetsfl1.clc
mets_for_model1 <- dbase[idx1]
mets_for_model1 <- fds[mets_for_model1]
mets_for_model2 <- dbase[idx2]
mets_for_model2 <- fds[mets_for_model2]




m1v <- train(f2v,
             data = mets_for_model1,
             method = "lm",
             trControl = trainControl(method="LOOCV"))
m2v <- train(f2v,
             data = mets_for_model2,
             method = "lm",
             trControl = trainControl(method="LOOCV"))

setkey(mets_for_model1, "id_placette")
setkey(mets_for_model2, "id_placette")


mets_12 <- mets_for_model1[mets_for_model2]
x <- ggplot(data=mets_12, aes(x=log(pfsumprof), y=log(i.pfsumprof)))+
  geom_point()+
  geom_text(aes(label=paste0(meanang, " ", id_placette," ", i.meanang)), position = position_nudge(y = -0.05))
print(x)

obs <- m2v$pred$obs
pred <- m2v$pred$pred
yobs <- exp(obs)
ypred <- exp(pred)
see <- sqrt(sum((obs-pred)^2)/(length(obs)-4))
cf <- exp((see^2)/2)
ypred <- ypred*cf
R2 <- cor(ypred, yobs)^2
1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))


mets

table(mets_for_model1$cl)
table(mets_for_model2$cl)





hist(test$value)  
dbase$id_placette <- as.factor(dbase$id_placette)
setkey(mets_for_model,"id_placette")
mets_for_model <- fds[mets_for_model]
mets3 <- mets_for_model[,c("volume_total_m3_ha", "meanch", "varch", 'pfsumprof', "cvladvox")]

m2l <- train(f2v,
             data = mets3,
             method = "lm",
             trControl = trainControl(method="LOOCV"))
vtot.l.coeff <- m2l$finalModel$coefficients
vtot.l.pred <- m2l$pred$pred
vtot.l.obs <- m2l$pred$obsmets3

obs <- vobs
pred <- vpred
yobs <- exp(obs)
ypred <- exp(pred)
see <- sqrt(sum((obs-pred)^2)/(length(obs)-4))
cf <- exp((see^2)/2)
ypred <- ypred*cf
R2 <- cor(ypred, yobs)^2
aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
RMSE <- sqrt(mean((yobs-ypred)^2))
MAE <- mean(abs(ypred-yobs))


ciron.mdlmets.all <- cbind(ciron.mdlmets.all, rep(1:5000, each=6))
vox.all.totv <- ciron.mdlmets.all[Forest_attr=="Total volume" & Metrics=="vox" & variable=="R2"]
hist(vox.all.totv$value)
idx <- as.vector(unlist(smplst.all[4267]))
test_mets <- pmetsfl1.all[idx]
pred <- ciron.all[[4267]]$vtot.vox.pred
obs <- ciron.all[[4267]]$vtot.vox.obs
yobs <- exp(obs)
ypred <- exp(pred)
see <- sqrt(sum((obs-pred)^2)/(length(obs)-4))
cf <- exp((see^2)/2)
ypred <- ypred*cf

test_mets <- as.data.table(cbind(test_mets, cbind(yobs, ypred)))

x1 <- ggplot(data=test_mets, aes(x=log(pflidr), y=log(pfsumprof)))+geom_point()+
  geom_text(aes(label=paste(id_placette," ", meanang)), vjust=1)+geom_abline()


ggplot(data=pmetsfl1.all, aes(x=log(pflidr), y=log(pfsumprof), colour=id_placette))+geom_point()+
  geom_text(aes(label=paste(id_placette," ", meanang)), vjust=1)+geom_abline()


pfs <- pmetsfl1.all[,c("id_placette", "pflidr", "pfsumprof")]
pfs <- pmetsfl1.all[,c("id_placette", "cvladlidr", "cvladvox")]




grid.arrange(x1, x2, ncol=2)
                 
                 


x <- list()
for(i in 1:5000)
{
  
  idx <- as.vector(unlist(smplst.clc[i]))
  x[i]<- diversity(pmetsfl1.clc[idx]$)
}


test_test <- as.data.table(cbind(yobs, ypred, seq(1:30)))

test_test1 <- as.data.table(cbind(vobs1, vpred1, seq(1:30)))


ggplot(data=test_test, aes(x=yobs, y=ypred))+geom_point()+
  geom_text(aes(label=V3),hjust=0, vjust=0)

idx <- as.vector(unlist(smplst.all[978]))
xxxx <- pmetsfl1.all[idx]



plot(vox.clc.totv$value, x)


pmetsfl1.clc <- pmetsfl1.clc[, clid := ifelse(cl=="a", 1, 
                                              ifelse(cl=="b", 2, 3))]


