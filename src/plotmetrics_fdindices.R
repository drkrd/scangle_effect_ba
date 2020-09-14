##Computation of metrics per flight line
df_df <- data.frame()
for(plotid in names(list_angs))
{
  aoi <- lasflightline(aois[[plotid]])
  flist <- unique(aoi@data$flightlineID)
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aoi, flightlineID == fl)
    meangle_fl <- abs(mean(aoi_subset@data$ScanAngleRank))
    ar <- area_calc(data.frame(x=aoi_subset@data$X, y=aoi_subset@data$Y))
    if(ar>650)
    {
      aoi_meanch <- func_meanch(aoi_subset@data$Z, aoi_subset@data$ReturnNumber)
      aoi_varch <- func_varch(aoi_subset@data$Z, aoi_subset@data$ReturnNumber)
      aoi_pf <- func_pf(aoi_subset@data$Z, aoi_subset@data$ReturnNumber)
      aoi_cvlad <- func_cvlad(aoi_subset@data$Z, aoi_subset@data$ReturnNumber)
      df_df <- as.data.frame(rbind(df_df, c(plotid,
                                            meangle_fl,
                                            aoi_meanch,
                                            aoi_varch,
                                            aoi_pf,
                                            aoi_cvlad)))
    }
  }
}
names(df_df) <- c("id_placette", "meanangle", "meanch", "varch", "gapfr", "cvlad")


##Calculation of the indices for field data
##Reading file without dead trees. There was a problem with R while creating subsets
##Dead tree rows were deleted in excel
fd <- readxl::read_xlsx("D:/1_Work/Dropbox/3_R_codes/Projects/scangle_effect_ba/data/Mesures_placette_Frisbee_all_plots_womort.xlsx",
                        sheet = "Arbres")
ind_df <- data.frame()
for(idp in unique(fd$id_placette))
{
  xx1 <- fd[which(fd$id_placette==idp), c(4,5,10,12)]
  xx1 <- xx1[which(xx1$dbh_in_cm!=0),]
  
  
  
  lst <- xx1 %>% group_by(species) %>% tally()
  sm=0
  for(i in lst$n)
  {
    a=i/sum(lst$n)
    b=log(i/sum(lst$n))
    sm=sm+(a*b)
  }
  sm=-1*sm #shannon
  pilou= sm/log(length(lst$n)) #pilou

  
  
  
  
  w <- owin(c(min(xx1$x), max(xx1$x)), c(min(xx1$y), max(xx1$y)))
  pp1 <- as.ppp(xx1,w)
  pp1_buff <- buffer(pp1, buf.xwid = 2.5, buf.ywid = 2.5)
  nnindices<-nnIndex(pp1_buff,
                     N=4,
                     smark=c("species", "dbh_in_cm", "x", "y"))
  M<-fsasN4(nnindices$nnspecies,match.fun=mingling)
  minglingind <- M$meanI#mingling
  W <- fsasN4(nnangle(nnindices$nndist,
                      nnindices$nnx,
                      nnindices$nny)$nnangle,
              match.fun=uniform.angle,para=72)
  uaind <- W$meanI#Uniform angle 
  ind_df <- as.data.frame(rbind(ind_df, c(idp, 
                                          minglingind, 
                                          uaind, 
                                          length(lst$n), 
                                          sm, 
                                          pilou,
                                          gini(xx1$dbh_in_cm))))
}
colnames(ind_df) <- c("id_placette",
                      "Mingling",
                      "UniformAngle",
                      "NoSpecies",
                      "Shannon",
                      "Pilou",
                      "GiniDBH")



final_df <- right_join(df_df, ind_df, by="id_placette")
