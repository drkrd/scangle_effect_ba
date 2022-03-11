library(sf)
library(basemaps)
library(raster)
library(ggplot2)
library(RStoolbox)
library(ggspatial)
library(data.table)
library(lidR)

ciron.poly <- read_sf("D:/1_Work/2_Ciron/Data/Shape_files/riparian_region.shp")
ciron.plots <- fread("D:/1_Work/2_Ciron/Data/Field/centreplacette2154.csv")
ciron.plots <- ciron.plots[, id_placette:=ifelse(id_placette==9.2,10,
                                                ifelse(id_placette>9.2, id_placette+1, id_placette))]
ciron.plots <- ciron.plots[id_placette!=13.0] 
ciron.plots$forest_type <- "Riparian"
ciron.plots <- st_as_sf(ciron.plots, coords = c("X","Y"))

st_crs(ciron.plots) <- 2154
crs(ciron.poly)
ciron.bbox <- as(extent(ciron.poly), 'SpatialPolygons') 
ciron.bbox <- st_as_sfc(ciron.bbox)
st_crs(ciron.bbox) <- 2154
ciron.bbox <- st_transform(ciron.bbox, 4326)
as.data.frame(ciron.plots)

ciron.bbox <- st_as_sfc(st_bbox(ciron.bbox))




ciron.raster <- basemap_raster(ciron.bbox, 
                               map_service = "osm", 
                               map_type = "topographic",
                               map_res = 0.3)


ciron.raster <- projectRaster(ciron.raster, crs=4326)
ciron.poly <- st_transform(ciron.poly, 4326)
ciron.plots <- st_transform(ciron.plots, 4326)


img.ciron <- 
  ggplot()+
  ggRGB(img=ciron.raster,
        r = 1,
        g = 2,
        b = 3,
        ggLayer = T)+
  geom_sf(data=ciron.poly, fill="blue", legend=TRUE)+
  geom_sf(data=ciron.plots, aes(colour=forest_type, shape=forest_type), size=3, shape=17, legend=TRUE)+
  annotation_scale(location = "br", width_hint = 0.5,
                   pad_x = unit(0.5, "in")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  theme(axis.text.y = element_text(angle=90),
        legend.position = c(0.15, 0.1),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.box="horizontal",
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill='white'),
        text=element_text(family="serif", size=12*(96/72)))+
  scale_color_manual(values="Purple")+
  coord_sf(xlim = c(-0.298, -0.035))+
  guides(color=guide_legend(override.aes=list(fill=NA)))
img.ciron




ggsave(img.ciron, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/ciron.plots.final.png", 
       width=12*1.25, height=12*1.25, units="cm", dpi=320) 












lascat73 <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/02_NUAGE_Z_L93/")
plot(lascat73)


bauges.db <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db15jul.csv")
exclude <- fread("D:/1_Work/__R_codes/Projects/dl_forestry_scangle/data/plots_a_exclure.txt", header = F)
bauges.db <- bauges.db[!Id_plac%in%exclude$V1]
bauges.db <- bauges.db[, newstratum:=ifelse(newstratum=="Coniferes", "Coniferous",
                                            ifelse(newstratum=="Feuillus", "Broadleaf", "Mixed"))]

bauges.db <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db15jul.csv", sep = ",")

#changing name from Id_plac to id_placette.
colnames(bauges.db)[2] <- "id_placette"  

#reclassifying the plots based on the G that includes G for small trees i.e. dbh>7.5cm ; computed by Kamel 
bauges.db <- bauges.db[, newstratum := ifelse(comp_R_G>75 & comp_F_G<25 , "Coniferes",
                                              ifelse(comp_F_G>75 & comp_R_G<25, "Feuillus", "Mixte"))]



bauges.plots <- st_as_sf(bauges.db, coords = c("X","Y"), crs=2154)
bauges.plots <- st_transform(bauges.plots, crs=4326)
x <-c(5.9, 6.5, 6.5, 5.9, 5.9)
y <-c(45.5, 45.5, 46.0, 46.0, 45.5)
xym <- cbind(x, y)
poly <- st_sfc(st_polygon(list(xym)), crs = 4326)


bauges.map <-basemap_raster(poly, 
                     map_service = "osm_stamen", 
                     map_type = "terrain_bg",
                     map_token = "",
                     map_res = 1.1)

bauges.map <- projectRaster(bauges.map, crs=4326)

img.bauges <- ggplot()+
  ggRGB(img=bauges.map,
        r = 1,
        g = 2,
        b = 3,
        ggLayer = T)+
  geom_sf(data=bauges.plots, size=2)+
  annotation_scale(location = "br", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0., "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  theme(axis.text.y = element_text(angle=90),
        legend.position = c(0.85, 0.88),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.box="horizontal",
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill='white'),
        text=element_text(family="serif", size=12*(96/72)))+
  scale_color_manual(values=c("Red"))+
  coord_sf()
img.bauges





ggsave(img.bauges, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/bauges.plots.final.png", 
       width=12*1.25, height=12*1.25, units="cm", dpi=640) 






ggplot(data=bauges.db, aes(x=log(G175)))+geom_histogram()+facet_wrap(newstratum~., scales = "free")















lbggmap(get_stamenmap((bbox = c(left = 5.9, 
                              bottom = 45.5, 
                              right =6.5, 
                              top = 45.9))),
      maptype = "terrain-background", zoom = 18)+
  geom_sf(data=bauges_plots, aes(colour=nwstrtm), inherit.aes = FALSE)

get_map(location=lat_bauges, 
        zoom = "auto", 
        maptype = "terrain-background",
        source = 
        )
bauges_map <- openmap(c(lat2, lon1), c(lat1, lon2), zoom = 10,
                      type = "esri-topo", mergeTiles = TRUE)
plot(bauges_map)


ggplot(data=bauges_map2)


ggmap(bauges_map2)+
  geom_sf(data=bauges_plots, aes(colour=nwstrtm))+
  theme_bw()

  