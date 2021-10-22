library(sf)
library(basemaps)
library(raster)
library(ggplot2)
library(RStoolbox)
library(ggspatial)
library(data.table)

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
                               map_res = 0.5,
                               map_dir = "D:/1_Work/Dropbox/2_Publications/2_paper/basemaps/")


ciron.raster <- projectRaster(ciron.raster, crs=4326)
ciron.poly <- st_transform(ciron.poly, 4326)
ciron.plots <- st_transform(ciron.plots, 4326)


img <- 
  ggplot()+
  ggRGB(img=ciron.raster,
        r = 1,
        g = 2,
        b = 3,
        ggLayer = T)+
  geom_sf(data=ciron.poly, fill="red", alpha=0.4)+
  geom_sf(data=ciron.plots, aes(shape="type"), size=2, shape=17, legend=TRUE)+
  annotation_scale(location = "br", width_hint = 0.5,
                   pad_x = unit(0.5, "in")) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.5, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  theme(axis.text.y = element_text(angle=90),
        legend.title = element_blank(),
        legend.box="horizontal",
        text=element_text(family="serif", size=12*(96/72)))+
  scale_color_manual(values=c("#512E5F"))+
  coord_sf(xlim = c(-0.298, -0.035))
img
ggsave(img, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/ciron_map.png", 
       width=12*1.25, height=12*1.25, units="cm", dpi=320) 






bauges.plots <- st_as_sf(bauges.db, coords = c("X","Y"), crs=2154)
bauges.plots <- st_transform(bauges.plots, crs=4326)
x <-c(5.9, 6.5, 6.5, 5.9, 5.9)
y <-c(45.5, 45.5, 46.0, 46.0, 45.5)
xym <- cbind(x, y)
poly <- st_sfc(st_polygon(list(xym)), crs = 4326)


bauges.map <- basemap_raster(poly, 
                     map_service = "osm", 
                     map_type = "topographic",
                     map_res = 0.5)

bauges.map <- projectRaster(bauges.map, crs=4326)

img <- ggplot()+
  ggRGB(img=bauges.map,
        r = 1,
        g = 2,
        b = 3,
        ggLayer = T)+
  geom_sf(data=bauges.plots, aes(colour=newstratum, shape=newstratum), size=2)+
  annotation_scale(location = "br", width_hint = 0.5) +
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0., "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering)+
  theme(axis.text.y = element_text(angle=90),
        legend.position = c(0.9, 0.9),
        legend.title = element_blank(),
        legend.box="horizontal",
        text=element_text(family="serif", size=12*(96/72)))+
  scale_color_manual(values=c("#512E5F", "#1B4F72", "#CA6F1E"))+
  coord_sf()
img



ggsave(img, file="D:/1_Work/Dropbox/2_Publications/2_paper/results/Trial1.png", 
       width=12*1.25, height=12*1.25, units="cm", dpi=320) 




















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

  