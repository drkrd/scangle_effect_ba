library(data.table)
library(dplyr)
library(lidR)
library(tools)

pc.list <- list.files("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/mar2022/allpoints_30m/", full.names = T)


for(name in pc.list)
{
  id_placette <- sub("\\_n.*", "", basename(tools::file_path_sans_ext(name)))
  name <- file_path_sans_ext(name)
  name <- paste0(name, ".las")
  aoi <- readLAS(name)
  print(area(aoi))
}




for(name in pc.list)
{
  id_placette <- sub("\\_un.*", "", basename(tools::file_path_sans_ext(name)))
  name <- file_path_sans_ext(name)
  name <- paste0(name, ".las")
  aoi <- readLAS(name)
  aoi <- retrieve_flightlines(aoi)
  aoi_n <- normalize_height(aoi, algorithm = tin())
  ext <- extent(aoi)
  xc <- ext@xmin+((ext@xmax-ext@xmin)/2)
  yc <- ext@ymin+((ext@ymax-ext@ymin)/2)
  aoi <- clip_circle(aoi, xc, yc, 15)
  name <- tools::file_path_sans_ext(name)
  # ext <- extent(aoi)
  # xc <- ext@xmin+((ext@xmax-ext@xmin)/2)
  # yc <- ext@ymin+((ext@ymax-ext@ymin)/2)
  # aoi <- clip_circle(aoi, xc, yc, 15)
  aoi <- retrieve_flightlines(aoi)
  # writeLAS(aoi,
  #          paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/allpoints/",
  #          name,"_n@all.las"))

  ###########################################################################################
  
  flist <- unique(aoi@data$flightlineID)
  fls <- c()
  angs <- c()
  for(fl in flist)
  {
    aoifl <- filter_poi(aoi, flightlineID == fl)
    aoifln <- filter_poi(aoi_n, flightlineID == fl)
    ar <- func_areacalc(data.frame(x=aoifl@data$X, y=aoifl@data$Y))
    meangle <- round(abs(mean(aoifl@data$ScanAngleRank)),2)
    if(ar>0.9*pi*15*15)
    {
      fls <- c(fls, fl)
      angs <- c(angs, meangle)
      nm1 <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/mar2022/flightlines_1/",
                   id_placette,"_un@", meangle, ".las")
      writeLAS(aoifl, nm1)
      
      nm2 <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/mar2022/flightlines_1/",
                   id_placette,"_n@", meangle, ".las")
      writeLAS(aoifln, nm2)
      
    }
  }
  #############################################################################################
  aoifl2 <- filter_poi(aoi, flightlineID %in% fls)
  nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/mar2022/allpoints_fl/",
               id_placette,"_un@allfls", ".las")
  writeLAS(aoifl2, nm)
  
  aoifl2n <- filter_poi(aoi_n, flightlineID %in% fls)
  nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/mar2022/allpoints_fl/",
               id_placette,"_n@allfls", ".las")
  writeLAS(aoifl2n, nm)
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
      nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unorm/plots/15m_rad/mar2022/flightlines_2/",
                   id_placette, "_un@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset2, nm)
      
      aoi_subset2n <- filter_poi(aoi, flightlineID %in% combos)
      nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/mar2022/flightlines_2/",
                   id_placette, "_n@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset2n, nm)
    }
    
    
  }
  else{
    aoi_subset2 <- filter_poi(aoi, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/mar2022/flightlines_2/",
                 id_placette, "_un@", paste(flstbl$angs, collapse = "_"), ".las")
    writeLAS(aoi_subset2, nm)
    
    aoi_subset2n <- filter_poi(aoi_n, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/mar2022/flightlines_2/",
                 id_placette, "_n@", paste(flstbl$angs, collapse = "_"), ".las")
    writeLAS(aoi_subset2n, nm)
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
      nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unorm/plots/15m_rad/mar2022/flightlines_3/",
                   id_placette, "_un@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset3, nm)
      
      aoi_subset3n <- filter_poi(aoi_n, flightlineID %in% combos)
      nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/mar2022/flightlines_3/",
                   id_placette, "_n@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset3n, nm)
    }
    
    
  }
  else{
    aoi_subset3 <- filter_poi(aoi, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unorm/plots/15m_rad/mar2022/flightlines_3/",
                 id_placette, "_un@", paste(flstbl$angs, collapse = "_"), ".las")
    writeLAS(aoi_subset3, nm)
    
    aoi_subset3n <- filter_poi(aoi_n, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/mar2022/flightlines_3/",
                 id_placette, "_n@", paste(flstbl$angs, collapse = "_"), ".las")
    writeLAS(aoi_subset3n, nm)
  }
}