library(lidR)
library(data.table)
library(ggplot2)

typefor <- "feuillus"
typepc <- "allpoints_fl"

plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
opt_independent_files(plots) <- TRUE ####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS

func_cvprof <- function(chunk)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  id_plac <- sub("\\@.*", "", basename(tools::file_path_sans_ext(chunk@files)))
  mang <- sub(".*@", "", basename(tools::file_path_sans_ext(chunk@files)))
  if(!is.na(as.numeric(mang)))
  {
    mang <- round(as.numeric(mang),2)
    
  }
  cvladlidr <- myProfilesLAD(las@data$Z,
                             Zmax= min(ceiling(max(las@data$Z)), 60),
                             dZ=1)
  l <- length(as.vector(unlist(cvladlidr)))
  k1 <- seq(1,l,1)
  k1 <- k1+2.5
  id_plac <- rep(id_plac,l)
  mang <- rep(mang,l)

  return(list(k1=k1,
              m=as.vector(unlist(cvladlidr)),
              meanang=mang,
              id_placette=id_plac))
}

padprofs_lidr <- catalog_apply(plots, func_cvprof)
padprofs_lidr <- do.call(rbind, lapply(padprofs_lidr, function(x) as.data.table(x, keep.rownames=FALSE)))
padprofs_lidr <- padprofs_lidr[, type:="lidr",]
padprofs_lidr <- padprofs_lidr[, meanang := ifelse(!is.na(as.numeric(meanang)),round(.SD,2),meanang), .SDcols="meanang"]
padprofs_lidr <- padprofs_lidr[, c("meanang","id_placette","type"):=lapply(.SD, as.factor), .SDcols=c("meanang","id_placette", "type")]


padprofs_vox <- do.call(rbind, lapply(voxmergedall, function(x) as.data.table(x, keep.rownames=FALSE)))
padprofs_vox <- padprofs_vox[, type:="vox",]
padprofs_vox <- padprofs_vox[k1>2]
padprofs_vox <- padprofs_vox[,c("id","test"):=NULL]
padprofs_vox <- padprofs_vox[, c("meanang","id_placette","type"):=lapply(.SD, as.factor), .SDcols=c("meanang","id_placette", "type")]


lidrvoxplots1 <- list()
for(plac in unique(padprofs_lidr$id_placette))
{
  tmp <- rbind(padprofs_lidr[id_placette==plac], padprofs_vox[id_placette==plac])
  tmp <- tmp[, c("meanang"):=lapply(.SD, factor), .SDcols=c("meanang")]
  lidrvoxplots1[[plac]] <- ggplot(data=tmp, aes(x=k1, y=m, colour = type))+
      scale_colour_brewer(palette="Set1")+
      #theme_minimal()+
      geom_line(size=0.6)+
      coord_flip()+
    #facet_wrap(~type)+
      labs(title = paste0("Vertical profiles for ", plac),
           y ="PAD",
           x = "Height above ground (m)",
           colour="Mean scan angle",
           linetype = "Mean scan angle")
  
  
  
   ggsave(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/feuillus/profiles/allpoints_fl/", plac, ".png"),
          lidrvoxplots1[[plac]], device = "png", dpi = "retina")
}
