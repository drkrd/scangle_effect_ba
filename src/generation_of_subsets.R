library(data.table)
library(lidR)

#read placette centre coordinates
ciron_coords <- setDT(read_xlsx("D:/1_Work/2_Ciron/Data/Field/Mesures_placette_Frisbee_all_plots_aout_2021.xlsx", 
                                sheet = "Info_Plots"))

#read las catalog 
ciron_lascat_norm <- readLAScatalog("Z:/_DATA/Ciron/ULM/FRISBEE_Ciron/01-Lidar/")

# clip plots
ciron_plots <- clip_circle(ciron_lascat_norm, 
                           ciron_coords$x_centre_plot, 
                           ciron_coords$y_centre_plot,
                           radius =  15)


names(ciron_plots) <- c(1:30)

for(name in names(ciron_plots))
{
  nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/august2021/allpoints/",
               name,"_n@all", ".las")
  writeLAS(ciron_plots[[name]], nm)
}


for(name in names(ciron_plots))
{
  aoi <- ciron_plots[[name]]
  ext <- extent(aoi)
  xc <- ext@xmin+((ext@xmax-ext@xmin)/2)
  yc <- ext@ymin+((ext@ymax-ext@ymin)/2)
  aoi <- clip_circle(aoi, xc, yc, 15)
  aoi <- retrieve_flightlines(aoi)
  # writeLAS(aoi,
  #          paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/",
  #          name,"_n@all.las"))
  
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
      nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/august2021/flightlines_1/",
                   name,"_n@", meangle, ".las")
      writeLAS(aoifl, nm)
    }
  }
  #############################################################################################
  aoifl2 <- filter_poi(aoi, flightlineID %in% fls)
  nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/august2021/allpoints_fl/",
               name,"_n@allfls", ".las")
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
      nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/august2021/flightlines_2/",
                   name, "_n@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset2, nm)
    }
    
    
  }
  else{
    aoi_subset2 <- filter_poi(aoi, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/august2021/flightlines_2/",
                 name, "_n@", paste(flstbl$angs, collapse = "_"), ".las")
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
      nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/august2021/flightlines_3/",
                   name, "_n@", paste(angs, collapse = "_"), ".las")
      writeLAS(aoi_subset3, nm)
    }
    
    
  }
  else{
    aoi_subset3 <- filter_poi(aoi, flightlineID %in% flstbl$fls)
    nm <- paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/august2021/flightlines_3/",
                 name, "_n@", paste(flstbl$angs, collapse = "_"), ".las")
    writeLAS(aoi_subset3, nm)
  }
}