library(lidR)
library(tidyr)
library(dplyr)
library(data.table)



fd <- readxl::read_xlsx("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/Mesures_placette_Frisbee_all_plots_sansmort.xlsx",
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
      ,.(ba=sum(ba_sqm, na.rm = TRUE)), keyby = .(id_placette, cx, cy)][
        ,sum_ba_hec:=(10000*ba)/706.8583,]

fd_ba_smry <- cbind(fd_ba_smry, seq(1:30))





plot_loc <- fd_ba_smry[,.(id_placette,cx,cy),]


##read las data into a catalog
lasc <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/")
lasc1 <- readLAS("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/14_un@all.las")


lasc1 <- clip_circle(lasc1, 
            xcenter = fd_ba_smry$cx[14], 
            ycenter = fd_ba_smry$cy[14], 
            radius = 9)

writeLAS(lasc ,paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/", "14_n@all.las")) 

##clipping plots based on plot centre coordinates
aois<- list()
for(i in 1:30)
{
  ls<- clip_circle(lasc, 
                   xcenter = fd_ba_smry$cx[i], 
                   ycenter = fd_ba_smry$cy[i], 
                   radius = 15)
   aois[[toString(fd_ba_smry$id_placette[i])]] <- ls  
   # writeLAS(ls,paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/march2021/allpoints/", 
   #                    toString(fd_ba_smry$id_placette[i]), "_n@all.las"))
}



aois <- aois1 
ar_th <- 0.9*pi*15*15
aois_fls <- list()
for(name in names(aois))
{
  # aoi <- readLAS(paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/", 
  #                name,
  #                ".las"))
  aoi <- aois[[name]]
  aoi <- retrieve_flightlines(aoi)
  flist <- unique(aoi@data$flightlineID)
  fls <- c()
  angs <- c()
  
  
  for(fl in flist)
  {
    aoi_subset <- filter_poi(aoi, flightlineID == fl)
    ar <- area_calc(data.frame(x=aoi_subset@data$X, y=aoi_subset@data$Y))
    meangle <- round(abs(mean(aoi_subset@data$ScanAngleRank)),2)
    if(ar>ar_th)
    {
      fls <- c(fls, fl)
      angs <- c(angs, meangle)
      # nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/march2021/flightlines/",
      #              name,"_n@", meangle)
      # writeLAS(aoi_subset,paste0(nm, ".las"))
      # aois_fls[[nm]] <- aoi_subset
    }
    
    # print(fls)
    # nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/march2021/allpointsflonly/",
    #              name,"_allfl")
    # aoi_ss <- lasfilter(aoi, flightlineID %in% fls)
    # aois_allfls[[nm]] <- aoi_ss
    # writeLAS(aoi_ss, paste0(nm, ".las"))
  }
  flstbl <- as.data.table(cbind(fls, angs))
  combos=2
  if(length(flstbl$fls)>=combos)
  {
    allcombos <- combn(fls, combos)
    for(i in 1:ncol(allcombos))
    {
      combos <- allcombos[,i]
      angs <- flstbl$angs[which(flstbl$fls%in%combos)]
      aoi_subset2 <- filter_poi(aoi, flightlineID %in% combos)
      nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/march2021/flightlines2/",
                   name, "_n@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset2, nm)
    }

    
  }
  else{
    aoi_subset2 <- filter_poi(aoi, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/march2021/flightlines2/",
                 name, "_n@", paste(flstbl$angs, collapse = "_"), ".las")
    writeLAS(aoi_subset2, nm)
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
