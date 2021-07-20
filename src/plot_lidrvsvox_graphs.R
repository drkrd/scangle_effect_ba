library(lidR)
library(data.table)
library(tools)
library(grid)
library(gridExtra)



plots <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/plots/15m_rad/april2021/all_w1fl/")
opt_independent_files(plots) <- TRUE ####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS
height=2

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
  ladlidr <- func_cvladprof(las@data$Z, las@data$ReturnNumber, ht=height)
  ladprof <- as.data.table(unlist(ladlidr))
  ladprof <- cbind(ladprof, 
                   names(ladlidr),
                   rep(mang, length(ladlidr)),
                   rep(id_plac, length(ladlidr)))
  colnames(ladprof) <- c("m", "k1", "meanang", "id_placette")
  ladprof <- ladprof[1:max(which(m!=0))]
  

  return(ladprof)
}
plotmetsfl1 <- catalog_apply(plots, func_computeall)
lidrprof <- rbindlist(plotmetsfl1)
lidrprof$k1 <- as.numeric(lidrprof$k1)

lidrprof[, n1:=sum(m*k1)/sum(m)]





allpcs <- list.files(paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/"),
                     pattern = "*.las",
                     full.names = TRUE)

alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)
names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))

allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/2_Ciron/voxelisation/Results/May/wo_interpolation/voxfiles/all_w1fl/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]
voxall <- rbindlist(apply(allvoxfiles, 1, 
                          func_normvox2, 
                          pth="D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/15m_rad_test/may2021/allpoints/", 
                          ht=height))

pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), 
                         sdvfp=sqrt(sum(m*(k1-(sum(k1*m)/sum(m)))^2)/(sum(m)*(length(m[which(m!=0)])-1)/length(m[which(m!=0)]))),
                         pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang)]

xy <- pfcvladvox[, .SD[sample(.N, min(1,.N))], by = id_placette]


profs <- rbind(cbind(lidrprof[,c("k1", "m", "id_placette", "meanang")], rep("lidr", nrow(lidrprof))),
               cbind(voxall[,c("k1", "m", "id_placette", "meanang")], rep("vox", nrow(voxall))))
profs$id_placette <- as.factor(profs$id_placette)
profs$V2 <- as.factor(profs$V2)
profs$meanang <- as.factor(profs$meanang)

x1 <- ggplot(data=profs[id_placette=="17"], aes(x=k1, y=m, colour=meanang))+
  geom_line()+
  facet_wrap(V2~.)+
  coord_flip()+
  scale_color_discrete(breaks = sort(as.numeric(rownames(meanang))))



x2 <- ggplot(data=lidrprof[id_placette=="6" & k1>2], aes(x=k1, y=m, colour=meanang))+
  geom_line()+
  geom_point(size=1)+
  coord_flip()

grid.arrange(x1,x2, ncol=2)


x1
