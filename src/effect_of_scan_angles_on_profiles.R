library(lidR)
library(dplyr)
library(conflicted)
mets_dbs_min <- data.frame()
ar_th <- 0.9*pi*15*15
pfratio10m <- data.frame()
cvprofs <- list()
for(name in names(aois))
{
  
  aoi <- lasflightline(aois[[name]])
  flist <- unique(aoi@data$flightlineID)
  meangle_fl <- data.frame()
  
  
  for(fl in flist)
  {
    aoi_subset <- lasfilter(aoi, flightlineID == fl)
    ar <- area_calc(data.frame(x=aoi_subset@data$X, y=aoi_subset@data$Y))
    meangle <- abs(mean(aoi_subset@data$ScanAngleRank))
    meangle_fl <- rbind(meangle_fl,c(meangle, fl, ar))
    if(ar>ar_th)
    {
      nm <- paste0(name, "_", round(meangle,2))
      cvprofs[[nm]] <- myProfilesLAD(aoi_subset@data$Z,
                                       Zmax= min(ceiling(max(aoi_subset@data$Z)), 60),
                                       dZ=1)
    }
  }
  
}



func_cvprof <- function(chunk)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  id_plac <- sub("_n.*", "", basename(tools::file_path_sans_ext(chunk@files)))
  mang <- round(as.numeric(sub(".*\\_", "", basename(tools::file_path_sans_ext(chunk@files)))),2)
  cvladlidr <- myProfilesLAD(las@data$Z,
                             Zmax= min(ceiling(max(las@data$Z)), 60),
                             dZ=1)
  l <- length(as.vector(unlist(cvladlidr)))
  k1 <- seq(1,l,1)
  k1 <- k1+1.5
  id_plac <- rep(id_plac,l)
  mang <- rep(mang,l)



  return(list(k1=k1,
              m=as.vector(unlist(cvladlidr)),
              meanang=mang,
              id_placette=id_plac))
}

lidrprofs <- catalog_apply(plots, func_cvprof)
lidrprofs <- do.call(rbind, lapply(lidrprofs, as.data.frame))
lidrprofs$type <- rep("lidr", nrow(lidrprofs))


lidrvoxplots <- list()
for(plac in unique(lidrprofs$id_placette))
{
  tmp <- rbind(lidrprofs[which(lidrprofs$id_placette==plac),], voxprofs[which(voxprofs$id_placette==plac),])
  lidrvoxplots[[plac]] <- ggplot(data=tmp, aes(x=k1, y=m, colour = meanang))+
      scale_colour_brewer(palette="Set1")+
      #theme_minimal()+
      geom_line(size=0.6)+
      coord_flip()+
      facet_grid(cols = vars(type))+
      labs(title = paste0("Vertical profiles for ", plac),
           y ="PAD",
           x = "Height above ground (m)",
           colour="Mean scan angle",
           linetype = "Mean scan angle")
  
  
  
  ggsave(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/", plac, ".png"),
         lidrvoxplots[[plac]], device = "png", dpi = "retina")
}



cvprofs <- list()
for(x in xx){
  nm <- paste0(x[1],"_",x[2])
  cvprofs[[nm]] <- unlist(x[[3]])
}

nms <- sapply(names(cvprofs), function(x) sub("\\_\\d+\\.\\d*","",x))
nms <- unique(nms)


allpads <- list()
padlidr <- list()
for(name in nms)
{
  sset <- cvprofs[sapply(names(cvprofs), function(x)
  {
  startsWith(x,paste0(name,"_"))},
  simplify=TRUE, USE.NAMES = FALSE)]
  
  # sset <- reshape2::melt(sset)
  # names(sset) <- c("m", "k1", "id")
  # sset$id <- as.factor(sset$id)
  # sset$k1 <- as.numeric(sset$k1)
  # padlidr[[name]] <- sset
  # 
  # voxmergedall1[[name]]$type <- 'vox'
  # voxmergedall1[[name]]$id <- as.factor(voxmergedall1[[name]]$id)
  # sset$type <- 'lidr'
  # 
  # sset <- sset[,c(3,2,1,4)]
  # 
  # voxmergedall1[[name]] <- as.data.frame(voxmergedall1[[name]][voxmergedall1[[name]]$k1>2,])
  # 
  # voxmergedall1[[name]]$k1 <- voxmergedall1[[name]]$k1-0.5
  # 
  # allpads[[name]] <- rbind(voxmergedall1[[name]][,c("id","k1","m","type")], sset) 
  # 
  # 
  # 
  # padplots[[name]] <- ggplot(data=allpads[[name]], aes(x=k1, y=m, colour = id, linetype = id))+
  #   scale_colour_brewer(palette="Set1")+
  #   #theme_minimal()+
  #   geom_line(size=0.8)+
  #   coord_flip()+
  #   facet_grid(cols = vars(type))+
  #   labs(title = paste0("Vertical profiles for", " plot ", name),
  #        y ="PAD",
  #        x = "Height above ground (m)",
  #        colour="Mean scan angle",
  #        linetype = "Mean scan angle")
  # 
  # 
  # ggsave(paste0("D:/1_Work/2_Ciron/voxelisation/Results/26jan/w_wt/lidrvsvox/", name, ".png"),
  #        padplots[[name]], device = "png", width = 400, height = 200, units = "mm", dpi = "retina")
  # 
}








