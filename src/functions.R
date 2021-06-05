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
func_meanch <- function(z,rn)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi))
  {
    return(NULL)
  }
  else
  {
    aoi_f <- aoi[aoi$rn==1,]
    aoi_f_h <- aoi_f[aoi_f$z>2 & aoi_f$z<60,]
    return(mean(aoi_f_h$z))
  }
}

func_varch <- function(z, rn)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi))
  {
    return(NULL)
  } 
  else
  {
    aoi_f <- aoi[aoi$rn==1,]
    aoi_f_h <- aoi_f[aoi_f$z>2 & aoi_f$z<60,]
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
  txt <- readLines(x[1], n=3)[2:3]
  zmin <- as.numeric(unlist(strsplit(txt[[1]], "\\s+"))[4])
  voxtbl <- fread(x[1], na.strings = "NA" , skip = 5)
  
  
  
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
  voxtbl <- voxtbl[,1:4][,alt := k+zmin+0.5][, dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  voxtbl <- voxtbl[alt>dtm]
  voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
  voxtbl <- voxtbl[, k1:= k1-0.5]
  voxtbl1 <- voxtbl[k1>2]
  ttl <- sum(voxtbl1$PadBVTotal, na.rm = TRUE)
  ttl <- ttl/(pi*15*15)
  pf <- exp(-0.5*ttl)
  voxtbl <- voxtbl[, .(m=mean(PadBVTotal, na.rm = TRUE), s=sum(is.nan(PadBVTotal))/length(PadBVTotal)), by=list(k1)]
  
  # voxtbl <- voxtbl[, .(s1 = sum(!is.nan(y))), by=list(k1)]
  # voxtbl <- voxtbl[, .(s1 = sum(is.nan(PadBVTotal))/length(PadBVTotal)), by=list(k1)]
  # voxtbl <- voxtbl[, .(s2 = sum(is.nan(PadBVTotal))), by=list(k1)]
  # 
  voxtbl <- voxtbl[, c("id_placette","meanang") := .(sub("\\_un.*", "", x[2]), x[3]),]
  voxtbl <- voxtbl[, pfsumvox:=pf,]
  voxtbl <- voxtbl[1:max(which(m!=0))]
  voxtbl <- voxtbl[k1>ht]
  
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
