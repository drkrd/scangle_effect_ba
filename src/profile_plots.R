


plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/flightlines_2/")
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
profilesfl2.ref <- catalog_apply(plots, func_computeall)
profilesfl2.ref <- rbindlist(profilesfl2.ref)
colnames(profilesfl2.ref) <- c("id_placette", "meanang", "PAD", "Height")
profilesfl2.ref$Height <- as.numeric(profilesfl2.ref$Height)
profilesfl2.ref$type <- rep("ref", nrow(profilesfl2.ref))
profilesfl2.ref <- profilesfl2.ref[, c("fl1","fl2"):=tstrsplit(meanang,"_",fixed=T),]
profilesfl2.ref$fl1 <- as.numeric(profilesfl2.ref$fl1)
profilesfl2.ref$fl2 <- as.numeric(profilesfl2.ref$fl2)
profilesfl2.ref <- profilesfl2.ref[, fl2:= ifelse(is.na(fl2), fl1, fl2),]

profilesfl2.ref <- profilesfl2.ref[, cl1:=ifelse(fl1>=0&fl1<10, "a",
                                   ifelse(fl1>=10&fl1<20, "b",
                                          ifelse(fl1>=20&fl1<30, "c",
                                                 ifelse(fl1>=30&fl1<40,"d","e"))))]
profilesfl2.ref <- profilesfl2.ref[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                   ifelse(fl2>=10&fl2<20, "b",
                                          ifelse(fl2>=20&fl2<30, "c",
                                                 ifelse(fl2>=30&fl2<40,"d","e"))))]

profilesfl2.ref <- profilesfl2.ref[cl1 != "e" & cl2 != "e" & cl1 != "d" & cl2 != "d",]

profilesfl2.ref <- profilesfl2.ref[!(id_placette=="96_IRSTEA" & fl1==30.0)]
profilesfl2.ref <- profilesfl2.ref[!(id_placette=="96_IRSTEA" & fl2==30.0)]
profilesfl2.ref <- profilesfl2.ref[!(id_placette=="283_74_ONF" & fl1==6.28)]
profilesfl2.ref <- profilesfl2.ref[!(id_placette=="283_74_ONF" & fl2==6.28)]

allpcs <- list.files(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/flightlines_1//"),
                     pattern = "*.las",
                     full.names = TRUE)

alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)
names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/flightlines_2//"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
profilesfl2.vox <- rbindlist(apply(allvoxfiles, 1, func_normvox2, 
                                pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", 
                                ht=height))
setDT(profilesfl2.vox)
profilesfl2.vox <- profilesfl2.vox[,c("id_placette", "meanang", "PADmean", "k1")]
colnames(profilesfl2.vox) <- c("id_placette", "meanang", "PAD", "Height")
profilesfl2.vox$type <- rep("vox", nrow(profilesfl2.vox))
profilesfl2.vox <- profilesfl2.vox[, c("fl1","fl2"):=tstrsplit(meanang,"_",fixed=T),]
profilesfl2.vox$fl1 <- as.numeric(profilesfl2.vox$fl1)
profilesfl2.vox$fl2 <- as.numeric(profilesfl2.vox$fl2)
profilesfl2.vox <- profilesfl2.vox[, fl2:= ifelse(is.na(fl2), fl1, fl2),]

profilesfl2.vox <- profilesfl2.vox[, cl1:=ifelse(fl1>=0&fl1<10, "a",
                                                 ifelse(fl1>=10&fl1<20, "b",
                                                        ifelse(fl1>=20&fl1<30, "c",
                                                               ifelse(fl1>=30&fl1<40,"d","e"))))]
profilesfl2.vox <- profilesfl2.vox[, cl2:=ifelse(fl2>=0&fl2<10, "a",
                                                 ifelse(fl2>=10&fl2<20, "b",
                                                        ifelse(fl2>=20&fl2<30, "c",
                                                               ifelse(fl2>=30&fl2<40,"d","e"))))]

profilesfl2.vox <- profilesfl2.vox[cl1 != "e" & cl2 != "e" & cl1 != "d" & cl2 != "d",]

profilesfl2.vox <- profilesfl2.vox[!(id_placette=="96_IRSTEA" & fl1==30.0)]
profilesfl2.vox <- profilesfl2.vox[!(id_placette=="96_IRSTEA" & fl2==30.0)]
profilesfl2.vox <- profilesfl2.vox[!(id_placette=="283_74_ONF" & fl1==6.28)]
profilesfl2.vox <- profilesfl2.vox[!(id_placette=="283_74_ONF" & fl2==6.28)]

profilesfl2.df <- rbind(profilesfl2.ref, profilesfl2.vox)
profilesfl2.df$id_placette <- as.factor(profilesfl2.df$id_placette)
profilesfl2.df$meanang <- as.factor(profilesfl2.df$meanang)



func_profilesfl2 <- function(profs)
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


profsfl2 <- NULL
for(plac in unique(profilesfl2.df$id_placette))
{
  profsfl2[[plac]] <- func_profilesfl2(profilesfl2.df[id_placette==plac])
}

