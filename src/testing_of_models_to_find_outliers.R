test <- bgs.mdlmetsfl1.feu.all[variable=="R2" & Metrics=="vox" & Forest_attr=="Basal area"]
ggplot(data=test, aes(x=value))+geom_histogram()
testh <- test[value>0.45]
testl <- test[value<0.45]
ggplot(data=testh, aes(x=value))+geom_histogram()
mdlh <- testh[id==977]
mdll <- testl[id==922]
idxh <- as.vector(unlist(smplstfl1.feu.clc[mdlh$id,]))
idxl <- as.vector(unlist(smplstfl1.feu.clc[mdll$id,]))
metsh <-pmetsfl1.feu.clc[idxh,]
metsl <-pmetsfl1.feu.clc[idxl,]
setkeyv(metsh, c("id_placette"))
setkeyv(metsl, c("id_placette"))
metshl <- metsh[metsl]
ggplot(data=metshl, aes(x=log(pfsumprof), y=log(i.pfsumprof)))+
  geom_point()+geom_text(aes(label=paste0(id_placette,i.meanang)))

(testh$id)
