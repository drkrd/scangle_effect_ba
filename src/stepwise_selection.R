library(dplyr)
library(conflicted)
library(tidyverse)
library(caret)
library(leaps)
library(MASS)
library(leaps)
df_stdmets_aoi <- data.frame()
df_stdmets <- data.frame()
all_dfs <- list()
all_dfs1 <- data.frame()
ar_th <- 0.9*pi*15*15
for(name in names(aois))
{
  aoi <- aois[[name]]
  flist <- unique(aoi@data$flightlineID)
  meangle_fl <- data.frame()
  
  meanch_orig <- func_meanch(aoi$Z, aoi$ReturnNumber)
  varch_orig <- func_varch(aoi$Z, aoi$ReturnNumber)
  pf_orig <- func_pf(aoi$Z, aoi$ReturnNumber)
  cvlad_orig <- func_cvlad(aoi$Z, aoi$ReturnNumber)
  all_dfs1 <- as.data.frame(rbind(all_dfs1, c(name,
                                              100,
                                              meanch_orig,
                                              varch_orig,
                                              pf_orig,
                                              cvlad_orig)))
                            
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aoi, flightlineID == fl)
    ar <- area_calc(data.frame(x=aoi_subset@data$X, y=aoi_subset@data$Y))
    meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
    meangle_fl <- rbind(meangle_fl,c(meangle, fl, ar))
  }
  
  angs <- c()
  means <- c()
  vars <- c()
  pfs <- c()
  cvlads <- c()
  meangle_fl <- meangle_fl[meangle_fl[,3]>ar_th,]
  df <- data.frame()
  for(nflid in meangle_fl[,2])
  {
    mean_ang <- meangle_fl[,1][which(meangle_fl[,2]==nflid)]
    aoi_fl <- lasfilter(aoi, flightlineID == nflid)
    meanch <- func_meanch(aoi_fl$Z, aoi_fl$ReturnNumber)
    varch <- func_varch(aoi_fl$Z, aoi_fl$ReturnNumber)
    pf <- func_pf(aoi_fl$Z, aoi_fl$ReturnNumber)
    cvlad <- func_cvlad(aoi_fl$Z, aoi_fl$ReturnNumber)
    all_dfs1 <- as.data.frame(rbind(all_dfs1, c(name,
                                                mean_ang,
                                                meanch,
                                                varch,
                                                pf,
                                                cvlad)))
    angs <- c(angs, mean_ang)
    means <- c(means, meanch )
    vars<- c(vars, varch )
    pfs<- c(pfs, pf )
    cvlads <- c(cvlads, cvlad )
  }
  
  df <- as.data.frame(cbind(angs, means, vars, pf, cvlad))
  df <- as.data.frame(rbind(df, c(100,
                                  meanch_orig, 
                                  varch_orig,
                                  pf_orig,
                                  cvlad_orig)))
  all_dfs[name] <- list(df)
  

}
names(all_dfs1) <- c("id_placette", "meanang", "meanch", "varch", "pf", "cvlad")

