


fit1v <- lm(f1v, data = mets_allx)
fit4v <- nls(f4v,
             data = mets_allx,
             weights = 1/VARx,
             start=list(a = exp(coef(fit1v)[1]), 
                        b = coef(fit1v)[2], 
                        c = coef(fit1v)[3],
                        d = coef(fit1v)[4],
                        e = coef(fit1v)[5]))
ypred <- predict(fit4v, mets_allx)
yobs <- mets_allx$volume_total_m3_ha
SSE <- sum((yobs-ypred)^2)
SST <- sum((mean(yobs)-yobs)^2)
R2 <- 1-(SSE/SST)
1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
(100/length(ypred))*sum((yobs-ypred)/yobs)
sqrt(mean((yobs-ypred)^2))
(sum(yobs-ypred))/length(yobs)
bias*100/mean(yobs)






pred <- predict(fit1v, mets_allx)
obs <- log(mets_allx$volume_total_m3_ha)
ypred <- exp(pred)
yobs <- exp(obs)
see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
cf <- exp((see^2)/2)
ypred <- ypred*cf
SSE <- sum((yobs-ypred)^2)
SST <- sum((mean(yobs)-yobs)^2)
R2 <- 1-(SSE/SST)
1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
(100/length(ypred))*sum((yobs-ypred)/yobs)
sqrt(mean((yobs-ypred)^2))
(sum(yobs-ypred))/length(yobs)
bias*100/mean(yobs)










f1l <- log(volume_total_m3_ha)~log(meanch)
f1v <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)

f4l <- (volume_total_m3_ha)~a*(meanch)^b
f4v <- (volume_total_m3_ha)~a*(meanch)^b*(varch)^c*(pfsumprof)^d*(sdvfp)^e


fx <- lm(f1l, data=mets_allx)
fy <- nls(f4l,
          data = mets_allx,
          weights = 1/VARx,
          start=list(a = exp(coef(fit1l)[1]), 
                     b = coef(fit1l)[2]))







fit5 <- glm(f5, data = mets_all, family = Gamma(link = "log"))







f1v <- (G_m2_ha)~(meanch)+(varch)+(pfsumprof)+(cvladvox)




f2l <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)

f2l <- (volume_total_m3_ha)~(meanch)+(varch)+(pflidr)+(cvladlidr)
f2v <- (volume_total_m3_ha)~(meanch)+(varch)+(pfsumprof)+(cvladvox)

f1 <- train(f2l, data=mets_all, method="glm" )
f2 <- glm(f2v, data=mets_all, family = poisson())


f3l <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(volume_tige_m3_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)

mets_all <- fd_smry_all[pmetsflall[which(id_placette %in% fd_smry_all$id_placette)]]
refmdls[["Basal area_lidr"]] <- as.data.frame(func_refmodls(mets_all, f1l))
refmdls[["Basal area_vox"]] <- as.data.frame(func_refmodls(mets_all, f1v))
refmdls[["Stem volume_lidr"]] <- as.data.frame(func_refmodls(mets_all, f2l))
refmdls[["Stem volume_vox"]] <- as.data.frame(func_refmodls(mets_all, f2v))
refmdls[["Total volume_lidr"]] <- as.data.frame(func_refmodls(mets_all, f3l))
refmdls[["Total volume_vox"]] <- as.data.frame(func_refmodls(mets_all, f3v))

mdlmets.refmdls <- rbindlist(refmdls, idcol="id")
mdlmets.refmdls <- mdlmets.refmdls[, c("Forest_attr", "Forest_type", "Metrics"):=tstrsplit(id,"_",fixed=T),][,-1]
mdlmets.refmdls <- melt(mdlmets.refmdls, id.vars = c("Forest_attr", "Forest_type", "Metrics"))

ggplot(data=mdlmets.refmdls, aes(x=Forest_type, y=value, fill=Metrics, colour=Forest_attr))+
  geom_bar(stat='identity', position = position_dodge())+facet_wrap(variable~., scales = "free")


mets_all1 <- mets_all 
mets_all1$pfsumprof[-6]

mdl<- train(f1l,
            data = mets_all,
            method = "lm",
            trControl = trainControl(method="LOOCV"))

summary(mdl)

yobs <-  mets_all$G_m2_ha
ypred <- m2$fitted.values
yobs <- exp(obs)
ypred <- exp(pred)
see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
cf <- exp((see^2)/2)
ypred <- ypred*cf
SSE <- sum((yobs-ypred)^2)
SST <- sum((mean(yobs)-yobs)^2)
R2 <- 1-(SSE/SST)
aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
RMSE <- sqrt(mean((yobs-ypred)^2))
RMSEpc <- RMSE*100/mean(yobs)
bias <- (sum(yobs-ypred))/length(yobs)
biaspc <- bias*100/mean(yobs)
xx <- as.data.table(cbind(ypred,yobs,pred,obs))
xx <- cbind(xx, mets_all$id_placette)


ggplot(data=xx, aes(x=yobs, y=ypred))+geom_abline()+
  geom_point()+geom_text(aes(label=V2), hjust = 0, nudge_x = 0.01*(max(ypred)-min(ypred)))+
  coord_fixed(xlim = c(min(yobs,ypred), 
                       max(yobs,ypred)), 
              ylim = c(min(yobs,ypred),
                       max(yobs,ypred)))+
  theme_few()+
  labs(x = "Observed",
       y = "Predicted")+
  annotate("text", 
           x = -Inf, y = Inf, hjust = 0, vjust = 1,
           label = paste0("R² = ", round(R2,2), "\n",
                          "RMSE = ", round(RMSE,2), "\n",
                          "RMSE% = ", round(RMSEpc,2), "\n",
                          "MPE% = ", round(MPE,2), "\n",
                          "Bias = ", round(bias,2), "\n",
                          "Bias% = ", round(biaspc,2)))

mets_mets <- mets_all[mets_all5m]


xyz <- data.frame("R2" =  round(aR2,2),
                  "RMSE" = round(RMSE,2),
                  "RMSEpc" = round(RMSEpc,2),
                  "MPE" =  round(MPE,2))

xyz <- rbind(xyz, c(aR2, RMSE, RMSEpc, MPE))

setDT(xyz)

xyz1 <- melt(xyz, variable.name = "variable")
xyz1 <- cbind(xyz1, "Metrics"=rep(c("old", "vox"), 4))




ggplot(data=pmetsflall, aes(x=log(pflidr), y=log(pfsumprof)))+
  geom_point()+geom_text(aes(label=id_placette))+geom_abline()+ 
  coord_fixed()+
  coord_fixed(xlim = c(min(pflidr, pfsumprof), 
                       max(pflidr, pfsumprof)), 
              ylim = c(min(pflidr, pfsumprof),
                       max(pflidr, pfsumprof)))



# print(yobs)
# print(ypred)
plot(yobs, ypred)
R2 <- cor(ypred, yobs)^2
aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
RMSE <- sqrt(mean((yobs-ypred)^2))
MAE <- mean(abs(ypred-yobs))
aR2
summary(mdl)




#######################################################################################################



SSE <- sum((mets_all5m$volume_total_m3_ha-xx$fitted.values)^2)
SST <- sum((mean(mets_all5m$volume_total_m3_ha)-mets_all5m$volume_total_m3_ha)^2)
R2 <- 1-(SSE/SST)
1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))


