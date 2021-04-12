library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)

typefor <- "coniferes"
typepc <- "flightlines"
# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
opt_independent_files(plots) <- TRUE ####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS

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
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber)
  return(list( 
    id_placette = as.factor(id_plac),
    meanang = as.factor(mang),
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}

plotmets <- catalog_apply(plots, func_computeall)
plotmets <- rbindlist(plotmets)




if(exists("alldtms"))
{
  if(!identical(sub("@.*","",names(alldtms)), as.character(unique(plotmets$id_placette))))
  {
    
    allpcs <- list.files(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/", typefor, "/15m_rad/allpoints/"),
                         pattern = "*.las",
                         full.names = TRUE)
    
    alldtms <- sapply(allpcs, function(x){
      ls <- readLAS(x)
      dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
      simplify = FALSE,
      USE.NAMES = TRUE)
    names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))
    
  }
}else
{
  allpcs <- list.files(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/", typefor, "/15m_rad/allpoints/"),
                       pattern = "*.las",
                       full.names = TRUE)
  
  alldtms <- sapply(allpcs, function(x){
    ls <- readLAS(x)
    dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
    simplify = FALSE,
    USE.NAMES = TRUE)
  names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))
}

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/", typefor, "/15m_rad/", typepc, "/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))


allvoxfiles[, grp := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
#allvoxfiles[, grp := basename(file_path_sans_ext(V1))]
allvoxfiles[, grp1 := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
#allvoxfiles <- allvoxfiles[allvoxfiles$grp1!="n",]
voxmergedall <- list()
allplots <- list()

func_normvox <- function(var1, zmin, dt_mat)
{
  voxtbl <- fread(var1, na.strings = "NA" , skip = 5)
  voxtbl <- voxtbl[,1:4][,alt := k+zmin+0.5][,dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
  #voxtbl <- na.omit(voxtbl)
  voxtbl <- voxtbl[alt>dtm]
  voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
  voxtbl <- voxtbl[, k1:= k1-0.5]
  voxtbl <- voxtbl[, .(m = mean(PadBVTotal, na.rm = TRUE)), by=list(k1)]
  return(voxtbl)
}

for(group in unique(allvoxfiles$grp))
{
  voxfiles <- allvoxfiles[allvoxfiles$grp==group,]
  txt <- readLines(voxfiles$V1[[1]], n=3)[2:3]
  zmin <- as.numeric(unlist(strsplit(txt[[1]], "\\s+"))[4])
  lasnm <- paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/", 
                  typefor, 
                  "/15m_rad/allpoints/",
                  group,
                  "@all",
                  ".las")
  ls <- readLASheader(lasnm)
  empty_raster <- raster(ncol=round(ls@PHB[["Max X"]]-ls@PHB[["Min X"]]), 
                         nrow=round(ls@PHB[["Max Y"]]-ls@PHB[["Min Y"]]),
                         xmn=ls@PHB[["Min X"]], xmx= ls@PHB[["Max X"]], 
                         ymn=ls@PHB[["Min Y"]], ymx= ls@PHB[["Max Y"]])
  dt <- alldtms[[paste0(group,"@all")]]
  dt <- resample(dt, empty_raster, method='bilinear')
  dt_mat <- as.matrix(dt)
  
  voxtbls <- sapply(voxfiles$V1, func_normvox, zmin=zmin, dt_mat=dt_mat,
                    simplify = FALSE, USE.NAMES = TRUE)
  
  voxmerged <- data.table::rbindlist(voxtbls, idcol = "meanang")
  voxmerged[, meanang := sub(".*@", "", basename(tools::file_path_sans_ext(meanang))),]

  voxmerged[, meanang := as.factor(ifelse(!is.na(as.numeric(meanang)), round(as.numeric(meanang),2), meanang)),]
  voxmerged[, id_placette := group,]
  voxmerged <- voxmerged %>% 
    mutate(test=if_else(lag(m)==0 & m==0,"x", "y")) %>% 
    filter(test!="x")
  voxmergedall[[group]] <- voxmerged
}



##########################
#Compute cvlad from voxels
##########################

cvladvox1 <- lapply(voxmergedall, function(x)
{
  x <- x[k1>2]
  return(x[,.(cvs=as.double(cv(m, na.rm = TRUE))), by=meanang])
})

cvladvox1 <- reshape2::melt(cvladvox1)
names(cvladvox1) <- c("meanang","cvladvox1","id_placette") 
setDT(cvladvox1)
cvladvox1$id_placette <- as.factor(cvladvox1$id_placette)


# cvladvox <- cvladvox %>% 
#   mutate(meanang=round(as.numeric(sub(".*\\_", "", basename(tools::file_path_sans_ext(id)))),2))


cvladvox1 <- cvladvox1[,`:=`(x=NULL)]
# cvladvox <- cvladvox[-124,]

plotmets <- right_join(plotmets, cvladvox1, by=c("id_placette", "meanang"))


















profile_plot <- function(profs)
{
  plt <- ggplot(data=profs[profs$k1>2,], aes(x=k1, y=m, colour = meanang, linetype = meanang))+
    scale_colour_brewer(palette="Paired")+
    theme_minimal()+
    geom_line(size=0.8)+
    coord_flip()+
    labs(title = paste0("Vertical profile for ", unique(profs$id_placette)),
         y ="PAD", 
         x = "Height above ground (m)", 
         colour="Mean scan angle", 
         linetype = "Mean scan angle")
  
  
  return(plt)
}

allplots <-  lapply(voxmergedall, profile_plot)


for(name in names(allplots))
{
  nm <- paste0("D:/1_Work/2_Ciron/voxelisation/Results/26jan/w_wt/",name,".png")
  ggsave(nm, allplots[[name]], device="png")
}


















voxtbl <- fread(allvoxfiles$V1[1], na.strings = "NA" , skip = 5)
voxtbl <- voxtbl[,1:4][,alt := k+zmin][,dtm := dt_mat[cbind(nrow(dt_mat)-j, i+1)]]
#voxtbl <- na.omit(voxtbl)
voxtbl <- voxtbl[alt>dtm]
voxtbl <- voxtbl[, k1:= order(k), by=list(i,j)]
voxtbl <- voxtbl[, .(m = mean(PadBVTotal, na.rm = TRUE)), by=list(k1)]
