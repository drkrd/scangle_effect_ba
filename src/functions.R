myProfilesLAD = function(Z, Zmax, dZ)
{
  # creating an empty list
  #print (max(Z))
  #print (Zmax)
  #print ("****")
  list_lad =list()
  
  #creating layer names
  
  z_ini=c(0, Zmax)
  lad_ini=LAD(z_ini, dz=dZ, k=0.5)
  
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
    profil_lad=LAD(Z, dz=dZ, k=0.5)
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

func_pf <- function(z, rn)
{
  aoi <- as.data.frame(cbind(z, rn))
  if (is.null(aoi))
  {
    return(NULL)
  }
  else
  {
    aoi_f <- aoi[aoi$rn==1,]
    num <- aoi_f[aoi_f$z<=2,]
    return(nrow(num)/nrow(aoi_f))
  }
}

func_cvlad <- function(z,rn)
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
    val = as.vector(unlist(myProfilesLAD(aoi$z, zm, dZ = 1)))
    return(sd(val) / mean(val))
  }
}



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
