library(lidR)

pc_copy <- pc
pc_copy <- list()

pc_copy <- lapply(pc, function(x)
{
  x <- lidR::retrieve_flightlines(x)
  return(x)
})


##reduce area of the point cloud to an area of a circle of 17.5m radius. 
#15m + 2.5m buffer for edge issues
aois <- lapply(pc_copy, function(x)
{
  ext <- extent(x)
  xc <- ext@xmin+((ext@xmax - ext@xmin)/2)
  yc <- ext@ymin+((ext@ymax - ext@ymin)/2)
  x <- clip_circle(x, xc, yc, 17.5)
  return(x)
})
slist_areas <- lapply(aois, function(x)
{
  flid_areas <- c()
  means <- c()
  meds <- c()
  plot_area <- lidR::area(x)
  flids <- unique(x@data$flightlineID)
  for(flid in flids)
  {
    flid_pc <- filter_poi(x, flightlineID==flid)
    pc_area <- lidR::area(flid_pc)
    pc_mean <- mean(abs(flid_pc@data$ScanAngleRank))
    pc_median <- median(abs(flid_pc@data$ScanAngleRank))
    flid_areas <- c(flid_areas, pc_area)
    means <- c(means, pc_mean)
    meds <- c(meds, pc_median)
  }
  df <- as.data.frame(cbind(flid_areas, means, meds))
  names(df) <- c("Area", "mean_scangle", "med_scangle")
  return(df)
})



df_areas <- data.frame()
for(name in names(list_areas))
{
  x <- as.data.frame(list_areas[[name]])
  x$id_plac <- rep(name, length(x$Area))
  df_areas <- rbind(df_areas, x)
}


df_areas <- df_areas %>% 
  mutate(fac=as.factor(if_else(id_plac%in%plots_731$Id_plac,"73",
                               if_else(id_plac%in%plots_741$Id_plac, "74", "Common"))),
         ar=as.factor(if_else(Area<0.9*pi*17.5*17.5,"0","1")),
         id_plac=as.factor(id_plac)) 

ggplot(data = df_areas, 
       aes(x=id_plac, y=mean_scangle, colour=ar))+
  geom_point()+
  labs(x = "Placettes", y="Mean scan angle (°)")+
  facet_grid(cols=vars(fac), scales = "free")+
  theme(axis.text.x = element_blank())+
  theme(legend.position="none")






list_nflids <- lapply(pc_copy, function(x)
{
  flid_areas <- c()
  counter=0
  plot_area <- lidR::area(x)
  flids <- unique(x@data$flightlineID)
  for(flid in flids)
  {
    flid_pc <- filter_poi(x, flightlineID==flid)
    pc_area <- lidR::area(flid_pc)
    per <- (pc_area/plot_area)*100
    if(per>90) counter=counter+1
  }
  return(counter)
})


list_mets <- lapply(pc_copy, function(x)
{
  meanch <- func_meanch(x@data$Z, x@data$ReturnNumber)
  varch <- func_varch(x@data$Z, x@data$ReturnNumber)
  pf <- func_pf(x@data$Z, x@data$ReturnNumber)
  cvlad <- func_cvlad(x@data$Z, x@data$ReturnNumber)
  
  return(list("meanch" = meanch,
              "varch" = varch,
              "pf" = pf,
              "cvlad" = cvlad))
})


dd  <-  as.data.frame(t(matrix(unlist(list_mets), 
                               nrow=length(unlist(list_mets[1])))))
rownames(dd) <- names(list_mets)
colnames(dd) <- c("meanch", "varch", "pf", "cvlad")


fd_summary <- fd %>% 
  dplyr::select(c(1,4)) %>% 
  mutate(ba_sqm = pi*((Diam/2)^2)) %>% 
  group_by(Id_plac) %>%
  summarise(sum_ba=sum(ba_sqm, na.rm = TRUE)) %>% #sum of basal area per plot
  mutate(sum_ba_hec=(10000*sum_ba)/2827.43) %>% #extrapolated to 1 hectare
  ungroup()

library(dplyr)




pc_copy <- list()

pc_copy <- lapply(pc, function(x)
{
  x <- lidR::retrieve_flightlines(x)
  return(x)
})


list_areas <- lapply(aois, function(x)
{
  flid_areas <- c()
  means <- c()
  meds <- c()
  plot_area <- lidR::area(x)
  flids <- unique(x@data$flightlineID)
  for(flid in flids)
  {
    flid_pc <- filter_poi(x, flightlineID==flid)
    pc_area <- lidR::area(flid_pc)
    pc_mean <- mean(abs(flid_pc@data$ScanAngleRank))
    pc_median <- median(abs(flid_pc@data$ScanAngleRank))
    flid_areas <- c(flid_areas, pc_area)
    means <- c(means, pc_mean)
    meds <- c(meds, pc_median)
  }
  df <- as.data.frame(cbind(flid_areas, means, meds))
  names(df) <- c("Area", "mean_scangle", "med_scangle")
  return(df)
})


list_nflids <- lapply(pc_copy, function(x)
{
  flid_areas <- c()
  counter=0
  plot_area <- lidR::areax)
  flids <- unique(x@data$flightlineID)
  for(flid in flids)
  {
    flid_pc <- filter_poi(x, flightlineID==flid)
    pc_area <- lidR::areaflid_pc)
    per <- (pc_area/plot_area)*100
    if(per>90) counter=counter+1
  }
  return(counter)
})


list_mets <- lapply(pc_copy, function(x)
{
  meanch <- func_meanch(x@data$Z, x@data$ReturnNumber)
  varch <- func_varch(x@data$Z, x@data$ReturnNumber)
  pf <- func_pf(x@data$Z, x@data$ReturnNumber)
  cvlad <- func_cvlad(x@data$Z, x@data$ReturnNumber)
  
  return(list(meanch, varch, pf, cvlad))
})


dd  <-  as.data.frame(t(matrix(unlist(list_mets), nrow=length(unlist(list_mets[1])))))
rownames(dd) <- names(list_mets)

fd <- read.csv2("D:/1_Work/4_Bauges/T1-Carto/Observatoire/Donnees/table_placette_v20190214_arbre.csv")
fd_summary <- fd %>% 
  dplyr::select(c(1,4)) %>% 
  mutate(ba_sqm = pi*((Diam/2)^2)) %>% 
  group_by(Id_plac) %>%
  summarise(sum_ba=sum(ba_sqm, na.rm = TRUE)) %>% #sum of basal area per plot
  mutate(sum_ba_hec=(10000*sum_ba)/2827.43) %>% #extrapolated to 1 hectare
  ungroup()


mets_all <- right_join(mets_dbs, fd_summary, by="Id_plac")

library(dplyr)




plots74 <- readxl::read_xls("D:/1_Work/Dropbox/3_R_codes/Projects/scangle_effect_ba/plots_in_74.xls")
plots73 <- readxl::read_xls("D:/1_Work/Dropbox/3_R_codes/Projects/scangle_effect_ba/plots_in_73.xls")
plots73 <- plots73[-which(plots73$Id_plac %in% plots74$Id_plac),]

df_co73 <- df_co[which(df_co$Id_plac %in% plots73$Id_plac),]
df_fe73 <- df_fe[which(df_fe$Id_plac %in% plots73$Id_plac),]
df_mi73 <- df_mi[which(df_mi$Id_plac %in% plots73$Id_plac),]


df_co74 <- df_co[which(df_co$Id_plac %in% plots74$Id_plac),]
df_fe74 <- df_fe[which(df_fe$Id_plac %in% plots74$Id_plac),]
df_mi74 <- df_mi[which(df_mi$Id_plac %in% plots74$Id_plac),]


