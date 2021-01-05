library(lidR)
library(dplyr)
library(conflicted)
mets <- list()

#calculation of area covered pc
#lidr has a simpler option to measure area if the input data is available as a las object using area()
area_calc = function(dfr)
{
  #print(length(dfr$x))
  ch_pts <- chull(dfr$x,dfr$y)
  ch_pts <- c(ch_pts, ch_pts[1])
  dfr <- dfr[ch_pts,]
  dfr <- dfr %>% 
    select(1:2) 
  ch_poly <- Polygon(dfr, hole=F)
  return(ch_poly@area)
}

#aois is a named list of point cloud
for(name in names(aois))
{
  aoi <- lasflightline(aois[[name]])
  flist <- unique(aoi@data$flightlineID)
  meangle_fl <- data.frame()
  
  #computation of mean scan angle and area for the point cloud from each flight line
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aoi, flightlineID == fl)
    ar <- area_calc(data.frame(x=aoi_subset@data$X, y=aoi_subset@data$Y))
    meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
    meangle_fl <- rbind(meangle_fl,c(meangle, fl, ar))
  }
  
  #dropping all flight lines which partially cover the plot
  meangle_fl <- meangle_fl[meangle_fl[,3]>650,]
  
  plot_df <- data.frame()#initialise a dataframe
  for(flid in meangle_fl[,2])
  {
    #filter pointcloud based on flight line
    aoi_tmp <- lasfilter(aoi, flightlineID == flid)
    mean_ang <- meangle_fl[,1][which(meangle_fl[,2]==flid)]
    #compute metrics
    aoi_meanch <- func_meanch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_varch <- func_varch(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_pf <- func_pf(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    aoi_cvlad <- func_cvlad(aoi_tmp@data$Z, aoi_tmp@data$ReturnNumber)
    #add to a dataframe
    plot_df <- as.data.frame(rbind(plot_df, c(flid,
                                              mean_ang,
                                              aoi_meanch, 
                                              aoi_varch, 
                                              aoi_pf, 
                                              aoi_cvlad)))
  }
  names(plot_df) <- c("flid", "angle", "meanch", "varch", "pf", "cvlad")
  
  mets[name] <- list(plot_df)
  
  
}

