library(lidR)
library(tidyr)
library(dplyr)

fd <- readxl::read_xlsx("D:/1_Work/Dropbox/3_R_codes/Projects/scangle_effect_ba/data/Mesures_placette_Frisbee_all_plots_womort.xlsx",
                        sheet = "Arbres")



fd <- as.data.frame(fd)

fd_ba <- fd[,c(1,2,3,10,11,18,19)]

fd_ba_smry <- fd_ba %>%
  #mutate(Mort=replace(Mort, Mort=="X", "x")) %>% #standardise dead tree entries
  #filter(is.na(Mort)) %>% #remove dead trees
  mutate(dbh_cm = (cbh_in_cm/pi)) %>% #circumference to diameter: c = 2piR = piD
  mutate(ba_sqm = pi*((dbh_cm/200)^2)) %>% #ba in sqm
  mutate(id_placette=as.factor(id_placette)) %>% 
  dplyr::select(c(1,2,3,7,8,9)) %>% 
  group_by(id_placette, cx, cy) %>% 
  summarise(sum_ba=sum(ba_sqm, na.rm = TRUE)) %>% #sum of basal area per plot
  mutate(sum_ba_hec=(10000*sum_ba)/706.8583) %>% #extrapolated to 1 hectare
  ungroup()


plot_loc <- fd_ba_smry[,c(1:3)]


##read las data into a catalog
lasc <- readLAScatalog("D:/1_Work/2_Ciron/Data/ULM/LAS/norm/")

##clipping plots based on plot centre coordinates
aois <- list()
for(i in 1:30)
{
  aois[[toString(fd_ba_smry$id_placette[i])]] <- clip_circle(lasc, 
                                                             xcenter = fd_ba_smry$cx[i], 
                                                             ycenter = fd_ba_smry$cy[i], 
                                                             radius = 15 )
}
