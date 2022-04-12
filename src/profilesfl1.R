

height=0
plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/flightlines_1/")
opt_independent_files(plots) <- TRUE ####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS
func_computeall <- function(chunk)
{
  las <- readLAS(chunk)                  # read the chunk
  if (is.empty(las)) return(NULL)        # check if it contains points
  
  
  id_plac <- sub("\\_n.*", "", basename(tools::file_path_sans_ext(chunk@files)))
  mang <- sub(".*@", "", basename(tools::file_path_sans_ext(chunk@files)))
  if(!is.na(as.numeric(mang)))
  {
    mang <- round(as.numeric(mang), 2)
    
  }
  cvladlidr <- func_cvladprof(las@data$Z, las@data$ReturnNumber, ht=height)
  cvladdf <- cbind("id_placette"=rep(id_plac, nrow(cvladlidr)),
                     "meanang"=rep(mang, nrow(cvladlidr)),
                     cvladlidr)

  return(cvladdf)
}
profiles.ref <- catalog_apply(plots, func_computeall)
profiles.ref <- rbindlist(profiles.ref)
colnames(profiles.ref) <- c("id_placette", "meanang", "PAD", "Height")
profiles.ref$Height <- as.numeric(profiles.ref$Height)
profiles.ref$type <- rep("ref", nrow(profiles.ref))
profiles.ref <- profiles.ref[, cl:=ifelse(meanang>=0&meanang<10, "a",
                                          ifelse(meanang>=10&meanang<20, "b",
                                                 ifelse(meanang>=20&meanang<30, "c",
                                                        ifelse(meanang>=30&meanang<40,"d","e"))))]
profiles.ref <- profiles.ref[cl!="d" & cl!="e"]

# The following flight lines are problematic. There is no voxelisation for them because there was a problem
# with the GPS time in the trajectory information, for these flight lines in particular.
# 283_74_ONF_un@6.28
# 96_IRSTEA_un@30
# Dropping the two problematic flight lines mentioned above 
profiles.ref <- profiles.ref[!(id_placette=="96_IRSTEA" & meanang==30.0)]
profiles.ref <- profiles.ref[!(id_placette=="283_74_ONF" & meanang==6.28)]



allpcs <- list.files(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/"),
                    pattern = "*.las",
                    full.names = TRUE)

alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)
names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/flightlines_1/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]





profiles.vox <- rbindlist(apply(allvoxfiles, 1, func_normvox2, 
                                pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", 
                                ht=height))






setDT(profiles.vox)

profiles.vox$meanang <- as.numeric(profiles.vox$meanang)


profiles.vox <- profiles.vox[,c("id_placette", "meanang", "PADmean", "k1")]
colnames(profiles.vox) <- c("id_placette", "meanang", "PAD", "Height")
profiles.vox$type <- rep("vox", nrow(profiles.vox))
profiles.vox <- profiles.vox[, cl:=ifelse(meanang>=0&meanang<10, "a",
                                          ifelse(meanang>=10&meanang<20, "b",
                                                 ifelse(meanang>=20&meanang<30, "c",
                                                        ifelse(meanang>=30&meanang<40,"d","e"))))]
profiles.vox <- profiles.vox[cl!="d" & cl!="e"]

# The following flight lines are problematic. There is no voxelisation for them because there was a problem
# with the GPS time in the trajectory information, for these flight lines in particular.
# 283_74_ONF_un@6.28
# 96_IRSTEA_un@30
# Dropping the two problematic flight lines mentioned above 
profiles.vox <- profiles.vox[!(id_placette=="96_IRSTEA" & meanang==30.0)]
profiles.vox <- profiles.vox[!(id_placette=="283_74_ONF" & meanang==6.28)]

profiles.df <- rbind(profiles.ref, profiles.vox)
profiles.df$id_placette <- as.factor(profiles.df$id_placette)
profiles.df$meanang <- as.factor(profiles.df$meanang)



func_profiles <- function(profs)
{
  plt <- ggplot(data=profs, aes(x=Height, y=PAD, colour = meanang))+
    scale_colour_brewer(palette="Paired")+
    theme_base()+
    geom_line(size=0.8)+
    coord_flip()+
    labs(title = paste0("Vertical profile for ", unique(profs$id_placette)),
         y ="PAD", 
         x = "Height above ground (m)", 
         colour="Mean scan angle", 
         linetype = "Mean scan angle")+
    facet_wrap(~type)+
    theme(text=element_text(family="serif", size=9*(96/72)),
          legend.position = "bottom")
  return(plt)
}


profsfl1 <- NULL
for(plac in unique(profiles.vox$id_placette))
{
  profsfl1[[plac]] <- func_profiles(profiles.vox[id_placette==plac])
}

