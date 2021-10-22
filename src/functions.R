func_btcor <- function(obs, pred)
{
  yobs <- exp(obs)
  ypred <- exp(pred)
  see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
  cf <- exp((see^2)/2)
  ypred <- ypred*cf
  return(ypred)
}

func_mdlmets <- function(obs, pred, for_attr, mettype)
{
  yobs <- exp(obs)
  ypred <- func_btcor(obs, pred)
  SSE <- sum((yobs-ypred)^2)
  SST <- sum((mean(yobs)-yobs)^2)
  R2 <- 1-(SSE/SST)
  aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
  MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
  RMSE <- sqrt(mean((yobs-ypred)^2))
  RMSEpc <- RMSE*100/mean(yobs)
  return(list( "Forest_attr"=for_attr, "Metrics"=mettype, "R2"=aR2, "RMSE"=RMSE,"rRMSE"=RMSEpc,"MPE"=MPE))
}



func_computeall <- function(chunk)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  
  
  id_plac <- sub("\\_n.*", "", basename(tools::file_path_sans_ext(chunk@files)))
  mang <- sub(".*@", "", basename(tools::file_path_sans_ext(chunk@files)))
  if(!is.na(as.numeric(mang)))
  {
    mang <- round(as.numeric(mang),2)
    
  }
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber, ht = height)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber, ht = height)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht = height)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht = height)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}



myProfilesLAD = function(Z, Zmax, dZ, th)
{
  # creating an empty list
  #print (max(Z))
  #print (Zmax)
  #print ("****")
  list_lad =list()
  
  #creating layer names
  
  z_ini=c(0, Zmax)
  lad_ini=LAD(z_ini, dz=dZ, k=0.5, z0=th)
  
  # initialisation of the list 
  
  for (i in 1:dim(lad_ini)[1])
  {
    list_lad[[i]] = 0
  }
  
  names(list_lad)=lad_ini$z   # adding names
  
  ##### Computation of PAD and converting the result into a list
  
  if (max(Z)>Zmax)
  {
    #print("********")
    #print(max(Z))
    #print("********")
  }
  
  Z=Z[Z<Zmax & Z>0] # to filter outliers and to keep only positive heights
  if(length(Z)>0)
  {
    profil_lad=LAD(Z, dz=dZ, k=0.5, z0=5)
    if (dim(profil_lad)[1]>0)    # test to identify empty PAD profile (no vegetation above 1 m)
    {
      for (i in 1:dim(profil_lad)[1])
      { 
        list_lad[[i]]=profil_lad$lad[i] 
      }
      
    }
  }
  
  #return the result of the function
  return(list_lad)
}
## All following functions compute the values from a dataframe of the las data containing Z and return number only
func_meanch <- function(z, rn, ht)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi))
  {
    return(NULL)
  }
  else
  {
    aoi_f <- aoi[aoi$rn==1,]
    aoi_f_h <- aoi_f[aoi_f$z>ht & aoi_f$z<60,]
    return(mean(aoi_f_h$z))
  }
}

func_varch <- function(z, rn, ht)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi))
  {
    return(NULL)
  } 
  else
  {
    aoi_f <- aoi[aoi$rn==1,]
    aoi_f_h <- aoi_f[aoi_f$z>ht & aoi_f$z<60,]
    return(var(aoi_f_h$z))
  }
}

func_pf <- function(z, rn, ht)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi))
  {
    return(NULL)
  }
  else
  {
    aoi_f <- aoi[aoi$rn==1,]
    num <- aoi_f[aoi_f$z<=ht,]
    return(nrow(num)/nrow(aoi_f))
  }
}

func_cvlad <- function(z, rn, ht)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi)) 
  {
    return(NULL)
  } 
  else
  {
    v <- NULL
    zm = min(ceiling(max(aoi$z)), 60)
    val = as.vector(unlist(myProfilesLAD(aoi$z, zm, dZ = 1, ht)))
    return(sd(val, na.rm = T) / mean(val, na.rm = T ))
  }
}

func_sdvfp <- function(z, rn, ht)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi)) 
  {
    return(NULL)
  } 
  else
  {
    v <- NULL
    zm = min(ceiling(max(aoi$z)), 60)
    val = list(myProfilesLAD(aoi$z, zm, dZ = 1, ht))
    val <- t(rbindlist(val))
    val <- as.data.frame(cbind(rownames(val), val))
    setDT(val)
    colnames(val) <- c("k1", "PADmean")
    val$PADmean <- as.numeric(val$PADmean)
    val$k1 <- as.numeric(val$k1)
    val <- val[, .(sdvfp=sqrt(sum(PADmean*(k1-(sum(k1*PADmean)/sum(PADmean)))^2)/(sum(PADmean)*(length(PADmean[which(PADmean!=0)])-1)/length(PADmean[which(PADmean!=0)]))))]

    return(unlist(val$sdvfp))
  }
}

func_cvladprof <- function(z, rn, ht)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi)) 
  {
    return(NULL)
  } 
  else
  {
    v <- NULL
    zm = min(ceiling(max(aoi$z)), 60)
    val = reshape2::melt(myProfilesLAD(aoi$z, zm, dZ = 1, ht))
    setDT(val)
    val <- val[1:max(which(value!=0))]
    
    return(val)
    }
}

func_areacalc = function(dfr)
{
  #print(length(dfr$x))
  if(nrow(dfr)<4)
  {
    return(0)
  }
  else
  {
    ch_pts <- chull(dfr$x,dfr$y)
    ch_pts <- c(ch_pts, ch_pts[1])
    dfr <- dfr[ch_pts,]
    dfr <- dfr %>% 
      select(1:2) 
    ch_poly <- Polygon(dfr, hole=F)
    return(ch_poly@area)
  }

}


func_normvox2 <- function(x, pth, ht)
{
  x <- as.vector(unlist(x))
  txt <- readLines(x[1], n=8)
  mincorner <- txt[pmatch("#min_corner", txt)]
  zmin <- as.numeric(str_extract_all(mincorner, "(\\d+\\.\\d+)")[[1]][3])
  voxtbl <- fread(x[1], na.strings = "NA" , skip = "i j k")
  
  
  
  #align DTM with voxel file
  lasnm <- paste0(pth, x[2], "@all.las")
  ls <- readLASheader(lasnm)
  empty_raster <- raster(ncol=round(ls@PHB[["Max X"]]-ls@PHB[["Min X"]]), 
                         nrow=round(ls@PHB[["Max Y"]]-ls@PHB[["Min Y"]]),
                         xmn=ls@PHB[["Min X"]], xmx= ls@PHB[["Max X"]], 
                         ymn=ls@PHB[["Min Y"]], ymx= ls@PHB[["Max Y"]])
  dt <- alldtms[[paste0(x[2],"@all")]]
  dt <- resample(dt, empty_raster, method='bilinear')
  dt_mat <- as.matrix(dt)
  
  
  #normalisation based on the DTM
  voxtbl <- voxtbl[,c("i","j","k","PadBVTotal")][, alt := k+zmin+0.5][, dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  voxtbl <- voxtbl[alt>dtm]
  voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
  voxtbl <- voxtbl[, k1:= k1-0.5]
  voxtbl <- voxtbl[k1>ht]
  # ttl <- sum(voxtbl1$PadBVTotal, na.rm = TRUE)
  # ttl <- ttl/(pi*15*15)
  # pf <- exp(-0.5*ttl)
  # voxtbl <- voxtbl[, n:=length(PadBVTotal[!is.na(PadBVTotal)]), by=list(k1)]
  # voxtbl <- voxtbl[n>30]
  voxtbl <- voxtbl[, .(PADmean=mean(PadBVTotal, na.rm = TRUE), 
                       # var=var(PadBVTotal, na.rm = TRUE),
                       # lci=t.test(PadBVTotal, na.rm = TRUE)$conf.int[1],
                       # uci=t.test(PadBVTotal, na.rm = TRUE)$conf.int[2],
                       s=sum(is.nan(PadBVTotal))/length(PadBVTotal)), by=list(k1)]
  
  # voxtbl <- voxtbl[, .(s1 = sum(!is.nan(y))), by=list(k1)]
  # voxtbl <- voxtbl[, .(s1 = sum(is.nan(PadBVTotal))/length(PadBVTotal)), by=list(k1)]
  # voxtbl <- voxtbl[, .(s2 = sum(is.nan(PadBVTotal))), by=list(k1)]
  # 
  voxtbl <- voxtbl[, c("id_placette","meanang") := .(sub("\\_un.*", "", x[2]), x[3]),]
  # voxtbl <- voxtbl[, pfsumvox:=pf,]
  voxtbl <- voxtbl[1:max(which(PADmean!=0))]
  # voxtbl <- voxtbl[k1>ht]
  
  return(voxtbl)
}


func_normvox3 <- function(x)
{
  txt <- readLines(x[1], n=4)[2:4]
  zmin <- as.numeric(unlist(strsplit(txt[[1]], "\\s+"))[4])
  #zmin <- round(zmin)
  voxtbl <- fread(x[1], na.strings = "NA" , skip = 5)
  
  
  
  #align DTM with voxel file
  xmin <- as.numeric(unlist(strsplit(txt[1], "\\s+"))[2])
  xmax <- as.numeric(unlist(strsplit(txt[2], "\\s+"))[2])
  ymin <- as.numeric(unlist(strsplit(txt[1], "\\s+"))[3])
  ymax <- as.numeric(unlist(strsplit(txt[2], "\\s+"))[3])
  rows <- as.numeric(unlist(strsplit(txt[3], "\\s+"))[3])
  cols <- as.numeric(unlist(strsplit(txt[3], "\\s+"))[2])
  
  empty_raster <- raster(ncol = cols, nrow = rows,
                         xmn = xmin, xmx = xmax, 
                         ymn = ymin, ymx = ymax)
  
  dt <- alldtms[[paste0(x[2],"@all")]]
  dt <- resample(dt, empty_raster, method='bilinear')
  dt <- round(dt)
  dt_mat <- as.matrix(dt)
  
  
  #normalisation based on the DTM
  voxtbl <- voxtbl[, 1:4][, alt := k+zmin+0.5][, dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  voxtbl <- voxtbl[alt>dtm]
  voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
  voxtbl <- voxtbl[, k1:= k1-0.5]
  ttl <- sum(voxtbl1$PadBVTotal, na.rm = TRUE)
  ttl <- ttl/(pi*15*15)
  pf <- exp(-0.5*ttl)
  voxtbl <- voxtbl[, .(m = mean(PadBVTotal, na.rm = TRUE)), by = k1]
  voxtbl <- voxtbl[k1>2]
  
  # voxtbl <- voxtbl[, .(s1 = sum(!is.nan(y))), by=list(k1)]
  # voxtbl <- voxtbl[, .(s1 = sum(is.nan(PadBVTotal))/length(PadBVTotal)), by=list(k1)]
  # voxtbl <- voxtbl[, .(s2 = sum(is.nan(PadBVTotal))), by=list(k1)]
  # 
  voxtbl <- voxtbl[, c("id_placette","meanang") := .(as.factor(x[2]), as.factor(x[3])),]
  voxtbl <- voxtbl[, pfsumvox:=pf,]
  voxtbl <- voxtbl[1:max(which(m!=0)),]
  return(voxtbl)
}
