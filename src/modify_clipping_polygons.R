library(doFuture)
registerDoFuture()
plan(multiprocess)

#generate DTMs from point cloud

allpcs <- list.files(paste0("D:/1_Work/2_Ciron/Data/ULM/LAS/unnorm/plots/17_5m_rad/"),
                     pattern = "*.las",
                     full.names = TRUE)

alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)


dtbl <- data.table(dtms=names(alldtms), id_placette=file_path_sans_ext(basename(names(alldtms))))


func_polygonedit <- function(dtm)
{
  ext <- extent(dtm)
  xc <- ext@xmin+(ext@xmax-ext@xmin)/2
  yc <- ext@ymin+(ext@ymax-ext@ymin)/2
  dtbl <- data.table("angle"=seq(0,360,0.9), 
                     "xc"=rep(xc,401), 
                     "yc"=rep(yc,401), 
                     "rad"=rep(15,401))
  dtbl[,c("xn","yn"):=.(xc+rad*cos(angle*(pi/180)), yc+rad*sin(angle*(pi/180))),]
  pts <- st_as_sf(dtbl[], coords = c("xn", "yn"), crs=2154)
  y <- apply(dtbl, 1, function(x)
  {
    dtbl1 <- data.table("angle"=rep(x[1],151), 
                        "xc"=rep(xc,151), 
                        "yc"=rep(yc,151),
                        "rad"=seq(0,15,0.1))
    dtbl1[,c("xn","yn"):=.(xc+rad*cos(angle*(pi/180)), yc+rad*sin(angle*(pi/180))),]
    dtbl1[,dist:=sqrt((xn-xc)^2+(yn-yc)^2)]
    pts1 <- st_as_sf(dtbl1, coords = c("xn", "yn"), crs=2154)
    vals <- raster::extract(dtm, pts1)
    dtbl1 <- cbind(dtbl1,vals)
    dtbl1[,dist2 := sqrt((dist-shift(dist))^2+(vals-shift(vals, 1L, type = "lag"))^2)]
    set(dtbl1, 1L, 9L, 0)
    dtbl1[,cumdist := cumsum(dist2)]
    dtbl1 <- cbind(dtbl1, "xn1"=rep(x[5], nrow(dtbl1)), "yn1"=rep(x[6], nrow(dtbl1)))
    dtbl1[,mind:=abs(cumdist-15),]
    return(dtbl1[,.SD[which.min(mind)],])
  })
  y <- rbindlist(y)
  poly1 <- st_cast(st_combine(st_as_sf(y, coords = c("xn", "yn"), crs=2154)), "POLYGON")
  poly2 <- st_cast(st_combine(st_as_sf(y, coords = c("xn1", "yn1"), crs=2154)), "POLYGON")
  return(list(dtm, poly1, poly2))
  
}

polys <- lapply(alldtms[1:5], func_polygonedit)

######################################################################################
library(parallel)
library(foreach)
library(doParallel)
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
polys2 <- foreach(i = alldtms, .final = function(i) setNames(i, names(alldtms)), .packages=c("dplyr", "raster", "sf", "data.table")) %dopar% {
  poly <- func_polygonedit(i)
}
stopCluster(clus)
#######################################################################################



xs <- lapply(seq_along(polys2), function(x){
  nm <- basename(file_path_sans_ext(names(polys2)[[x]])) 
  x1 <- as.data.frame(polys2[[x]][[1]], xy=T)
  ggplot()+geom_raster(data = x1, aes(x=x, y=y, fill=Z))+
    scale_fill_gradientn(colours = terrain.colors(10), na.value="transparent")+
    theme_minimal()+
    geom_sf(data=polys2[[x]][[2]], colour="red", fill=NA)+
    geom_sf(data=polys2[[x]][[3]], colour="black", fill=NA)+
    labs(title=paste0("Plot", nm), x="Longitude", y="Latitude", fill="Altitude")
})

library(gridExtra)

ggsave(
  filename = "plots.pdf", 
  plot = marrangeGrob(xs, nrow=1, ncol=1), 
  width = 15, height = 9
)

plots__woterrain <- sapply(seq_along(polys2), function(x){
  return(st_as_sf(polys2[[x]][[3]]))
})


plots2 <- st_combine(c(plots__woterrain[[1]],
           plots__woterrain[[2]],
           plots__woterrain[[3]],
           plots__woterrain[[4]],
           plots__woterrain[[5]],
           plots__woterrain[[6]],
           plots__woterrain[[7]],
           plots__woterrain[[8]],
           plots__woterrain[[9]],
           plots__woterrain[[10]],
           plots__woterrain[[11]],
           plots__woterrain[[12]],
           plots__woterrain[[13]],
           plots__woterrain[[14]],
           plots__woterrain[[15]],
           plots__woterrain[[16]],
           plots__woterrain[[17]],
           plots__woterrain[[18]],
           plots__woterrain[[19]],
           plots__woterrain[[20]],
           plots__woterrain[[21]],
           plots__woterrain[[22]],
           plots__woterrain[[23]],
           plots__woterrain[[24]],
           plots__woterrain[[25]],
           plots__woterrain[[26]],
           plots__woterrain[[27]],
           plots__woterrain[[28]],
           plots__woterrain[[29]],
           plots__woterrain[[30]]))
