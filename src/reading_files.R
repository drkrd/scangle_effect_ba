library(lidR)
library(tidyr)
library(dplyr)
library(data.table)



fd <- readxl::read_xlsx("D:/1_Work/Dropbox/3_R_codes/Projects/scangle_effect_ba/data/Mesures_placette_Frisbee_all_plots_sansmort.xlsx",
                        sheet = "Arbres")# The original file does not include dead trees. They were deleted before-hand

# fd <- as.data.frame(fd)
# fd_ba_smry <- fd %>%
#   #mutate(Mort=replace(Mort, Mort=="X", "x")) %>% #standardise dead tree entries
#   #filter(is.na(Mort)) %>% #remove dead trees
#   mutate(ba_sqm = pi*(cbh_in_cm/(pi*200))^2) %>% #ba in sqm; D = C/(pi*100) =>  A = pi*(D/2)² => a=pi*(C/(pi*200))²
#   mutate(id_placette=as.factor(id_placette)) %>% 
#   dplyr::select(c(id_placette, cx, cy, ba_sqm)) %>% 
#   group_by(id_placette, cx, cy) %>% 
#   summarise(sum_ba=sum(ba_sqm, na.rm = TRUE)) %>% #sum of basal area per plot
#   mutate(sum_ba_hec=(10000*sum_ba)/706.8583) %>% #extrapolated to 1 hectare
#   ungroup()


##data.table implementation
setDT(fd)
fd_ba_smry <- fd[,ba_sqm := pi*(cbh_in_cm/(pi*200))^2,][
  ,id_placette:=as.factor(id_placette),][
    ,.(id_placette, cx, cy, ba_sqm),][
      ,.(sum_ba=sum(ba_sqm, na.rm = TRUE)), keyby = .(id_placette, cx, cy)][
        ,sum_ba_hec:=(10000*sum_ba)/706.8583,]





plot_loc <- fd_ba_smry[,.(id_placette,cx,cy),]


##read las data into a catalog
lasc <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/17_5m_rad/")

##clipping plots based on plot centre coordinates
aois_un <- list()
for(i in 1:30)
{
  ls<- clip_circle(lasc, 
                   xcenter = fd_ba_smry$cx[i], 
                   ycenter = fd_ba_smry$cy[i], 
                   radius = 15)
   aois_un[[toString(fd_ba_smry$id_placette[i])]] <- ls  
   writeLAS(ls,paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/", toString(fd_ba_smry$id_placette[i]), ".las"))
}




ar_th <- 0.9*pi*15*15
pf74_mi_fls <- list()
for(name in names(pc74_mi))
{
  # aoi <- readLAS(paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/", 
  #                name,
  #                ".las"))
  aoi <- pc74_mi[[name]]
  aoi <- lasflightline(aoi)
  flist <- unique(aoi@data$flightlineID)
  meangle_fl <- data.frame()
  
  
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aoi, flightlineID == fl)
    ar <- area_calc(data.frame(x=aoi_subset@data$X, y=aoi_subset@data$Y))
    meangle <- round(abs(mean(aoi_subset@data$ScanAngleRank)),2)
    if(ar>ar_th)
    {
      nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/mixte/15m_rad/flightlines/",
                   name,"_", meangle)
      writeLAS(aoi_subset,paste0(nm, ".las"))
      pf74_mi_fls[[nm]] <- aoi_subset
    }
  }
}




pc74_co_n <- list()
for(name in names(pc74_co)){
  pc74_co_n[[name]] <- normalize_height(pc74_co[[name]], algorithm = tin())
  # writeLAS(nrm, paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/feuillus/15m_rad/allpoints/",
  #                                  name,"_n.las"))
}



pc74 <- lapply(pc74_co, function(x)
{
  
})
