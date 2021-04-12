library(lidR)
library(dplyr)
library(future)
library(e1071)
lasc <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/")
rs <- 30 #resolution of the grid
area_filt <- 0.9*rs*rs #90 percent of the cell or more should have point cloud. 

plan(multisession, workers = 8L)
opt_chunk_buffer(lasc) <- 0
opt_chunk_size(lasc)   <- 600
opt_stop_early(lasc) <- TRUE
#plot(lasc, chunk_pattern = TRUE)
opt_select(lasc)       <- "xyzScanAngleRankgpstimeReturnNumber"
opt_filter(lasc) <- "-drop_class 2"


opt <- list(raster_alignment = rs, automerge = TRUE)


func_computeall <- function(chunk, res)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  las <- retrieve_flightlines(las, dt = 30)
  out <- grid_metrics(las, func_scanclasses(X, Y, Z, ScanAngleRank, flightlineID), rs)
  bbox   <- raster::extent(chunk)
  out <- raster::crop(out, bbox)
  return(out)
}

func_scanclasses <- function(x, y, z, sc, flid)
{
  dtbl <- as.data.frame(cbind(x, y, z, sc, flid))
  #dframe <- dframe[dframe$z>2,]
  #added the following because some flight lines had less than 3 points. Area computation was returning an error.
  dtbl <- dtbl %>%
    group_by(flid) %>%
    filter(n() >= 1000) %>% #1000 is an arbitrary number. Generally, a flight line that covers an entire area of 30m x 30m has no. of points far greater than 1000
    ungroup()
  # dtbl[, area_calc(c(.BY, .SD)), by=.(x,y)]
  flist <- unique(dtbl$flid)
  meanlist <- c()
  arlist <- c()
  flist1 <- c()
  fl <- data.frame()
  for (i in flist)
  {
    dtbl1 <- dtbl[which(dtbl$flid == i),]
    
    flist1 <- c(flist1, i)
    
    
    ch_ar <- area_calc(dtbl1)
    arlist <- c(arlist, ch_ar)

    
    mean_val <- (abs(mean(dtbl1$sc)))
    meanlist <- c(meanlist, mean_val)
    
  }
  fl <- as.data.frame(cbind(flist1, arlist, meanlist))
  fl <- fl[which(fl$arlist>0.9*30*30),]
  if(nrow(fl)>0)
  {
    fl <- fl %>%
      mutate(class=ifelse(meanlist>=0&meanlist<10,"a",
                           ifelse(meanlist>=10&meanlist<20,"b",
                                   ifelse(meanlist>=20&meanlist<30,"c",
                                           ifelse(meanlist>=30&meanlist<40,"d","e")))))
    val <- length(unique(fl$class))

    return(as.double(val))
    
  }else{
    return(NA_real_)
  }
}

rast_classes <- catalog_apply(lasc, func_computeall, res = rs, .options = opt)


plot(rast_classes)
