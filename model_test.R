mets_xx <- read.csv("D:/1_Work/mets_ciron2m.csv")

func_mdlmets <- function(pred, obs, cor, ref)
{

  if(cor==T)
  {
    yobs <- exp(obs)
    ypred <- exp(pred)
    #correction for back-transformation
    see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
    cf <- exp((see^2)/2)
    ypred <- ypred*cf
  }
  else{
    yobs=obs
    ypred=pred
  }
  
  SSE <- sum((yobs-ypred)^2)
  SSR <- sum((ypred-mean(yobs))^2)
  SST <- sum((mean(yobs)-yobs)^2)
  R2 <- 1-(SSE/SST)
  R2.adj <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
  RMSE <- sqrt(mean((yobs-ypred)^2))
  MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
  bias <- (sum(yobs-ypred))/length(yobs)
  biaspc <- bias*100/mean(yobs)
  return(data.table("R2"=R2.adj, "RMSE"=RMSE, "MPE"=MPE, "Ref"=ref))
}





combination <- "4"
if(combination=="1"){
  type="lidr.cvlad"
  eqn.log <- log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  eqn.nls <- G_m2_ha~a*(meanch)^b*(varch)^c*(pflidr)^d*(cvladlidr)^e
  
}else if(combination=="2"){
  type="lidr.sdvfp"
  eqn.log <- log(G_m2_ha)~log(meanch)+log(varch)+log(pflidr)+log(sdvfplidr)
  eqn.nls <- G_m2_ha~a*(meanch)^b*(varch)^c*(pflidr)^d*(sdvfplidr)^e
  
}else if(combination=="3"){
  type="vox.cvlad"
  eqn.log <- log(G_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  eqn.nls <- G_m2_ha~a*(meanch)^b*(varch)^c*(pfsumprof)^d*(cvladvox)^e
  
}else if(combination=="4"){
  type="vox.sdvfp"
  eqn.log <- log(G_m2_ha)~log(meanch)+log(varch)+log(pfsumprof)+log(sdvfp)
  eqn.nls <- G_m2_ha~a*(meanch)^b*(varch)^c*(pfsumprof)^d*(sdvfp)^e
}else(print("CHOOSE SOMETHING"))

#creating a copy of the data
mets_allx <- mets_all

#original log-log model
fit.log <- lm(eqn.log,
                 data = mets_allx)
#computation of model metrics
pred <- predict(fit.log)
obs <- log(mets_allx$G_m2_ha)
# obs <- log(mets_allx$G_m2_ha)
# obs <- log(mets_allx$volume_tige_m3_ha)
assign(paste0("log", combination), func_mdlmets(pred, obs, T, paste0("log.", type)))


# 
# fit.nls <- nls(G_m2_ha~a*meanch^b,
#                data = mets_allx,
#                start=list(a = exp(coef(fit.log)[1]),
#                           b = coef(fit.log)[2]))
# 


#nls without weighting
fit.nls <- nls(eqn.nls,
               data = mets_allx,
               start=list(a = exp(coef(fit.log)[1]),
                          b = coef(fit.log)[2],
                          c = coef(fit.log)[3],
                          d = coef(fit.log)[4],
                          e = coef(fit.log)[5]))

yht <- predict(fit.nls)
y <- mets_allx$G_m2_ha
assign(paste0("nls",combination), func_mdlmets(pred=yht, obs=y, F, paste0("nls.", type)))


#predicted values
nls.pred <- predict(fit.nls)
#residuals
nls.res <- mets_allx$G_m2_ha-nls.pred
#standardised residuals for plot
nls.res.std <- nls.res/sd(nls.res)

dfx <- data.frame("id"=mets_allx$id_placette, nls.pred, nls.res.std)
#plot of residuals vs predicted values
ggplot(data=dfx, aes(x=nls.pred, y=nls.res.std))+
  geom_point()+
  geom_text(aes(label=id), nudge_x = 0.6)+
  geom_hline(yintercept = 0)+
  theme_base()

#computation of weights
wt <- 1/lm(abs(nls.res) ~ nls.pred)$fitted.values^2

xx <- data.frame("id"=mets_allx$id_placette, wt)
#weighted nls
fit.wnls <- nls(eqn.nls,
                data = mets_allx,
                weights = wt,
                start=list(a = exp(coef(fit.log)[1]), 
                           b = coef(fit.log)[2], 
                           c = coef(fit.log)[3],
                           d = coef(fit.log)[4],
                           e = coef(fit.log)[5]))




y <- mets_allx$G_m2_ha
# yobs <- mets_allx$G_m2_ha
# yobs <- mets_allx$volume_tige_m3_ha
yht <- predict(fit.wnls)

assign(paste0("wnls",combination), func_mdlmets(pred=yht, obs=y, F, paste0("wnls.", type)))
       





mdlmets3 <- rbind(log1, log2, log3, log4, nls1, nls2, nls3, nls4, wnls1, wnls2, wnls3, wnls4)
mdlmets3 <- melt(mdlmets3, id.vars = "Ref")

mdlmets3 <- mdlmets3[, c("regtype","mettype", "profmet"):=tstrsplit(Ref,".",fixed=T),]

ggplot(data=mdlmets3, aes(x=profmet, y=abs(value), colour=mettype, group=mettype))+
  geom_line(arrow=arrow())+
  geom_point()+
  facet_grid(variable~regtype, scales = "free")+
  theme(axis.text.x = element_text(angle = 45))+
  theme_base()





