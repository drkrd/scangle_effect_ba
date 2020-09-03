mets_dbs <- data.frame()
for(name in names(aois))
{
  aoi_tmp <- aois[[name]]
  aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
  mets_dbs <- as.data.frame(rbind(mets_dbs, c(name,
                                              aoi_meanch,
                                              aoi_varch,
                                              aoi_pf,
                                              aoi_cvlad)))
}
names(mets_dbs) <- c("id_placette", "meanch", "varch", "pf", "cvlad")
mets_dbs$id_placette <- as.factor(mets_dbs$id_placette)
mets_dbs$meanch <- as.numeric(mets_dbs$meanch)
mets_dbs$varch <- as.numeric(mets_dbs$varch)
mets_dbs$pf <- as.numeric(mets_dbs$pf)
mets_dbs$cvlad <- as.numeric(mets_dbs$cvlad)
mets_dbs <- right_join(fd_ba_smry, mets_dbs, by="id_placette")

model <- lm(data = mets_dbs[-c(21,22),], 
            formula = log(sum_ba_hec)~log(meanch)+log(varch)+log(pf)+log(cvlad))
