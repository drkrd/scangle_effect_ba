library(data.table)
library(dplyr)

pc.list <- list.files("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/december2021/allpoints/")

for(name in pc.list)
{
  
  aoi <- readLAS(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/december2021/allpoints/", name))
  name <- tools::file_path_sans_ext(name)
  # ext <- extent(aoi)
  # xc <- ext@xmin+((ext@xmax-ext@xmin)/2)
  # yc <- ext@ymin+((ext@ymax-ext@ymin)/2)
  # aoi <- clip_circle(aoi, xc, yc, 15)
  aoi <- retrieve_flightlines(aoi)
  # writeLAS(aoi,
  #          paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/allpoints/",
  #          name,"_un@all.las"))

  ###########################################################################################
  
  flist <- unique(aoi@data$flightlineID)
  fls <- c()
  angs <- c()
  for(fl in flist)
  {
    aoifl <- filter_poi(aoi, flightlineID == fl)
    ar <- func_areacalc(data.frame(x=aoifl@data$X, y=aoifl@data$Y))
    meangle <- round(abs(mean(aoifl@data$ScanAngleRank)),2)
    if(ar>0.9*pi*15*15)
    {
      fls <- c(fls, fl)
      angs <- c(angs, meangle)
      nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/december2021/flightlines_1/",
                   name,"_un@", meangle, ".las")
      writeLAS(aoifl, nm)
    }
  }
  #############################################################################################
  aoifl2 <- filter_poi(aoi, flightlineID %in% fls)
  nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/december2021/allpoints_fl/",
               name,"_un@allfls", ".las")
  writeLAS(aoifl2, nm)
  ##############################################################################################
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
      nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/december2021/flightlines_2/",
                   name, "_un@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset2, nm)
    }
    
    
  }
  else{
    aoi_subset2 <- filter_poi(aoi, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/december2021/flightlines_2/",
                 name, "_un@", paste(flstbl$angs, collapse = "_"), ".las")
    writeLAS(aoi_subset2, nm)
  }
  ##############################################################################################
  combos=3
  if(length(flstbl$fls)>=combos)
  {
    allcombos <- combn(fls, combos)
    for(i in 1:ncol(allcombos))
    {
      combos <- allcombos[,i]
      angs <- flstbl$angs[which(flstbl$fls%in%combos)]
      aoi_subset3 <- filter_poi(aoi, flightlineID %in% combos)
      nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/december2021/flightlines_3/",
                   name, "_un@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset3, nm)
    }
    
    
  }
  else{
    aoi_subset3 <- filter_poi(aoi, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/december2021/flightlines_3/",
                 name, "_un@", paste(flstbl$angs, collapse = "_"), ".las")
    writeLAS(aoi_subset3, nm)
  }
}