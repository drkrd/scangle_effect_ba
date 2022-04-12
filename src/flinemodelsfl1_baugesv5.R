library(data.table)
library(tools)
library(lidR)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(gstat)
library(caret)
library(parallel)
library(foreach)
library(doParallel)
library(stringr)


# fixing a height threshold of 5 metres. This was based on a few observations on Ciron
height <- 5 

bauges.db <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db15jul.csv", sep = ",")

#changing name from Id_plac to id_placette.
colnames(bauges.db)[2] <- "id_placette"  

#reclassifying the plots based on the G that includes G for small trees i.e. >7.5cm ; computed by Kamel 
bauges.db <- bauges.db[, newstratum := ifelse(comp_R_G>75 & comp_F_G<25 , "Coniferes",
                                              ifelse(comp_F_G>75 & comp_R_G<25, "Feuillus", "Mixte"))]

#setting a key for the database
setkey(bauges.db, "id_placette")

#subset only the necessary columns
bauges.db <- bauges.db[, c("id_placette","G175" ,"G75", "volume_total", "volume_tige", "stratum", "newstratum")]
bauges.db <- bauges.db[, `:=`(volume_total=(volume_total*10000)/(pi*15*15),
                              volume_tige=(volume_tige*10000)/(pi*15*15))]
bauges.db2 <- readxl::read_xlsx("data/table_placette - PNR74.xlsx", sheet = "placette" )
bauges.db2 <- as.data.table(bauges.db2)
bauges.db74 <- bauges.db2[A_exclure=="N"]
bauges.db <- bauges.db[id_placette%in%bauges.db74$Id_plac]
bauges.db <- bauges.db[!id_placette%in%c("96_IRSTEA", "283_74_ONF")]
#subset based on type of forest
bauges.db.con <- bauges.db[newstratum=="Coniferes"]
bauges.db.feu <- bauges.db[newstratum=="Feuillus"]
bauges.db.mix <- bauges.db[newstratum=="Mixte"]



####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
plots <- readLAScatalog("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/15m_rad/may2021/flightlines_1/")
####VERY IMPORTANT WHEN PROCESSING INDIVIDUAL PLOTS AS A CATALOG
opt_independent_files(plots) <- TRUE
height=5


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
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber, ht = height)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber, ht = height)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht = height)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht = height)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
#compute metrics for each point cloud 
plotmetsfl1 <- catalog_apply(plots, func_computeall)
#make a copy of the list of list of metrics per point cloud
pmetsfl1 <- plotmetsfl1
#change the list to a data table
pmetsfl1 <- rbindlist(pmetsfl1)
# add a column 'cl' and classify each point cloud based on the mean scan angle (MSA) with which it was scanned
# 5 class system. MSA is already a part of the point clouds names. 
# While computing the metrics, func_computeall extracts the MSA.
pmetsfl1 <- pmetsfl1[, cl:=ifelse(meanang>=0&meanang<10, "a",
                                  ifelse(meanang>=10&meanang<20, "b",
                                         ifelse(meanang>=20&meanang<30, "c",
                                                ifelse(meanang>=30&meanang<40,"d","e"))))]

#there are a total of 494 flight lines at this point


# classes 'd' and 'e' are not considered in this site.
# High inclinations may be due to the flights making turns. Can't say with certainty
# We are dropping all the flight lines with MSA > 30°
pmetsfl1 <- pmetsfl1[cl!="d" & cl!="e"]

#there are a total of 476 flight lines at this point

# The following flight lines are problematic. There is no voxelisation for them because there was a problem
# with the GPS time in the trajectory information, for these flight lines in particular.
# 283_74_ONF_un@6.28
# 96_IRSTEA_un@30



##############################################################################################################################
######## Generation of DTMs ###########
#######################################

# List all point clouds for which DTMs need to be generated
# Read the point clouds that contain all flight lines. 
# These point clouds were subset from the data with plot centre coordinates
allpcs <- list.files(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/"),
                     pattern = "*.las",
                     full.names = TRUE)

# generate DTMs for each of the point clouds and store in a named list.
# names are plot names
alldtms <- sapply(allpcs, function(x){
  ls <- readLAS(x)
  dtm <- grid_terrain(ls, res = 0.1, algorithm = tin())},
  simplify = FALSE,
  USE.NAMES = TRUE)
names(alldtms) <- basename(file_path_sans_ext(names(alldtms)))
##############################################################################################################################
######## Profile generation from voxels ###########
###########################################################################
#list the vox files 
allvoxfiles <- as.data.table(list.files(paste0("D:/1_Work/5_Bauges/Voxelisation/Results/74/may2021/voxfiles/flightlines_1/"), 
                                        pattern = "*.vox",
                                        full.names = TRUE))
#extract the placette ids and MSA which are in the names of the files
allvoxfiles[, id_placette := sub("\\@.*","",basename(file_path_sans_ext(V1)))]
allvoxfiles[, meanang := sub(".*\\@","",basename(file_path_sans_ext(V1)))]

#convert MSA to numeric
allvoxfiles$meanang <- as.numeric(allvoxfiles$meanang)

#classification of the values for each flight line into the same classification as above to frop unwanted data
allvoxfiles <- allvoxfiles[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                     ifelse(meanang>=10&meanang<20, "b",
                                            ifelse(meanang>=20&meanang<30, "c",
                                                   ifelse(meanang>=30&meanang<40,"d","e"))))]
allvoxfiles <- allvoxfiles[cl!="d" & cl!="e"]

#generating profiles and computing metrics from the same
#returns data in a tabular format
voxall <- rbindlist(apply(allvoxfiles, 1, 
                          func_normvox2, pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", 
                          ht=height))
#convert to a data table
setDT(voxall)
##############################################################################################################################
######## Computation of metrics from profiles ###########
#########################################################

#Computing three metrics from the profile values
pfcvladvox <- voxall[, .(cvladvox=cv(PADmean, na.rm = TRUE), 
                         sdvfp=sqrt(sum(PADmean*(k1-(sum(k1*PADmean)/sum(PADmean)))^2)/(sum(PADmean)*(length(PADmean[which(PADmean!=0)])-1)/length(PADmean[which(PADmean!=0)]))),
                         pfsumprof=exp(-0.5*sum(PADmean, na.rm = TRUE))), by=.(id_placette, meanang)]

#convert MSA to numeric
pfcvladvox$meanang <- as.numeric(pfcvladvox$meanang)


pfcvladvox <- pfcvladvox[id_placette%in%bauges.db$id_placette]
pmetsfl1 <- pmetsfl1[which(id_placette %in% bauges.db$id_placette)]
#only 420 flight lines remain  now
#setting keys to join metrics
setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(pmetsfl1, c("id_placette", "meanang"))
pmetsfl1 <- pmetsfl1[pfcvladvox]


#drop some of the data that have pf of 0 because these values will cause a problem with log models  
pmetsfl1 <- pmetsfl1[!pflidr==0]
#there are 406 flight lines at this point. This is the final number that was used in the study

######################################
#############All################
#############################################################################################################################
{
  pmetsfl1 <- unique(pmetsfl1[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1 <- pmetsfl1[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1 <- pmetsfl1[, wt := wt/sum(wt), id_placette]
}

func_class.ssets <- function(pmets, cls)
{
  #tabulate per plot the number of pcs belonging to each class
  tbl <- with(pmetsfl1, table(id_placette, cl))
  
  if(cls=="a")
  {
    #pick all the plots that have atleast one pc belonging to class a
    pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,1]>0),]))]
    
    #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
    pmets.onlycl <- pmets[cl=="a"]
    
    #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
    pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,1]>0),]))]
    
    #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
    pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
  }
  else if(cls=="b")
  {
    #pick all the plots that have atleast one pc belonging to class a
    pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,2]>0),]))]
    
    #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
    pmets.onlycl <- pmets[cl=="b"]
    
    #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
    pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,2]>0),]))]
    
    #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
    pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
  }  
  else
  {
    #pick all the plots that have atleast one pc belonging to class a
    pmets.withcl <- pmets[which(id_placette %in% rownames(tbl[which(tbl[,3]>0),]))]
    
    #pick all pcs belonging to class a ONLY from each of the plots extracted in the previous step
    pmets.onlycl <- pmets[cl=="c"]
    
    #pick all the remaining plots that don't satisfy the above criteria to make sure all plots are picked at the end
    pmets.wocl <- pmets[!which(id_placette %in% rownames(tbl[which(tbl[,3]>0),]))]
    
    #combine the plots which have only pcs belonging to class a and the plots wo any pcs belonging to class a
    pmets.cl <- rbind(pmets.onlycl, pmets.wocl) 
  }
  return(pmets.cl)
}
######################################
#############Coniferes################
#############################################################################################################################

{
  # select those plots that are classified 'coniferous'
  pmetsfl1.con.all <- pmetsfl1[which(id_placette %in% bauges.db.con$id_placette)]
  # compute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.con.all <- unique(pmetsfl1.con.all[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.con.all <- pmetsfl1.con.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot 
  pmetsfl1.con.all <- pmetsfl1.con.all[, wt := wt/sum(wt), id_placette]


  #we will analyse the effect of having only data acquired such that the MSA belonged to class a at all plots
  #similarly, we will also analyse if the data has MSA belonging to class b at all plots
  #since, it is hard to have all classes at all locations, we need to find plots that satisfy this criteria
  #in the following steps, we will find all the plots that have atleast one fl from each class
  #from these plots, we can then subset datasets for only class A, only class B etc
  #therefore the final results become comparable across the experiments albeit with some noise from the other plots 
  
  
  #tabulate per plot the number of flight lines (fl) belonging to each class
  tbl <- with(pmetsfl1.con.all, table(id_placette, cl))
  
  #pick the plots that have at least one fl belonging to class a, class b and class c each
  tbl2 <- tbl[which(tbl[,1]>0 & tbl[,2]>0 & tbl[,3]>0),]
  
  #subset those plots and their metrics which satisfy the above condition
  pmetsfl1.con.withallcls <- pmetsfl1.con.all[id_placette %in% rownames(tbl2)]
  
  #subset those plots and their metrics which do not satisfy the previous condition
  pmetsfl1.con.woallcls <- pmetsfl1.con.all[!id_placette %in% rownames(tbl2)]
  
  
  
  
  
  #pick only those fls and their metrics which belong to class a
  pmetsfl1.con.cla <- pmetsfl1.con.withallcls[cl=="a"]
  
  #combine with the remaining plots
  pmetsfl1.con.cla <- rbind(pmetsfl1.con.cla, pmetsfl1.con.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.con.cla <- unique(pmetsfl1.con.cla[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.con.cla <- pmetsfl1.con.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot 
  pmetsfl1.con.cla <- pmetsfl1.con.cla[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.con.clb <- pmetsfl1.con.withallcls[cl=="b"]
  
  #combine with the remaining plots
  pmetsfl1.con.clb <- rbind(pmetsfl1.con.clb, pmetsfl1.con.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.con.clb <- unique(pmetsfl1.con.clb[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.con.clb <- pmetsfl1.con.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot
  pmetsfl1.con.clb <- pmetsfl1.con.clb[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.con.clc <- pmetsfl1.con.withallcls[cl=="c"]
  
  #combine with the remaining plots
  pmetsfl1.con.clc <- rbind(pmetsfl1.con.clc, pmetsfl1.con.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.con.clc <- unique(pmetsfl1.con.clc[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.con.clc <- pmetsfl1.con.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot
  pmetsfl1.con.clc <- pmetsfl1.con.clc[, wt := wt/sum(wt), id_placette]
  
  
}

######################################
#############Feuillus#################
###########################################################################################################################

{
  # select those plots that are classified 'coniferous'
  pmetsfl1.feu.all <- pmetsfl1[which(id_placette %in% bauges.db.feu$id_placette)]
  # compute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.feu.all <- unique(pmetsfl1.feu.all[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.feu.all <- pmetsfl1.feu.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot 
  pmetsfl1.feu.all <- pmetsfl1.feu.all[, wt := wt/sum(wt), id_placette]
  
  
  #we will analyse the effect of having only data acquired such that the MSA belonged to class a at all plots
  #similarly, we will also analyse if the data has MSA belonging to class b at all plots
  #since, it is hard to have all classes at all locations, we need to find plots that satisfy this criteria
  #in the following steps, we will find all the plots that have atleast one fl from each class
  #from these plots, we can then subset datasets for only class A, only class B etc
  #therefore the final results become comparable across the experiments albeit with some noise from the other plots 
  
  
  #tabulate per plot the number of flight lines (fl) belonging to each class
  tbl <- with(pmetsfl1.feu.all, table(id_placette, cl))
  
  #pick the plots that have at least one fl belonging to class a, class b and class c each
  tbl2 <- tbl[which(tbl[,1]>0 & tbl[,2]>0 & tbl[,3]>0),]
  
  #subset those plots and their metrics which satisfy the above condition
  pmetsfl1.feu.withallcls <- pmetsfl1.feu.all[id_placette %in% rownames(tbl2)]
  
  #subset those plots and their metrics which do not satisfy the previous condition
  pmetsfl1.feu.woallcls <- pmetsfl1.feu.all[!id_placette %in% rownames(tbl2)]
  
  
  
  
  
  #pick only those fls and their metrics which belong to class a
  pmetsfl1.feu.cla <- pmetsfl1.feu.withallcls[cl=="a"]
  
  #combine with the remaining plots
  pmetsfl1.feu.cla <- rbind(pmetsfl1.feu.cla, pmetsfl1.feu.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.feu.cla <- unique(pmetsfl1.feu.cla[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.feu.cla <- pmetsfl1.feu.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot 
  pmetsfl1.feu.cla <- pmetsfl1.feu.cla[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.feu.clb <- pmetsfl1.feu.withallcls[cl=="b"]
  
  #combine with the remaining plots
  pmetsfl1.feu.clb <- rbind(pmetsfl1.feu.clb, pmetsfl1.feu.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.feu.clb <- unique(pmetsfl1.feu.clb[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.feu.clb <- pmetsfl1.feu.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot
  pmetsfl1.feu.clb <- pmetsfl1.feu.clb[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.feu.clc <- pmetsfl1.feu.withallcls[cl=="c"]
  
  #combine with the remaining plots
  pmetsfl1.feu.clc <- rbind(pmetsfl1.feu.clc, pmetsfl1.feu.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.feu.clc <- unique(pmetsfl1.feu.clc[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.feu.clc <- pmetsfl1.feu.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot
  pmetsfl1.feu.clc <- pmetsfl1.feu.clc[, wt := wt/sum(wt), id_placette]
  
  
}

######################################
################Mixte#################
###########################################################################################################################

{
  # select those plots that are classified 'coniferous'
  pmetsfl1.mix.all <- pmetsfl1[which(id_placette %in% bauges.db.mix$id_placette)]
  # compute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.mix.all <- unique(pmetsfl1.mix.all[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.mix.all <- pmetsfl1.mix.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot 
  pmetsfl1.mix.all <- pmetsfl1.mix.all[, wt := wt/sum(wt), id_placette]
  
  
  #we will analyse the effect of having only data acquired such that the MSA belonged to class a at all plots
  #similarly, we will also analyse if the data has MSA belonging to class b at all plots
  #since, it is hard to have all classes at all locations, we need to find plots that satisfy this criteria
  #in the following steps, we will find all the plots that have atleast one fl from each class
  #from these plots, we can then subset datasets for only class A, only class B etc
  #therefore the final results become comparable across the experiments albeit with some noise from the other plots 
  
  
  #tabulate per plot the number of flight lines (fl) belonging to each class
  tbl <- with(pmetsfl1.mix.all, table(id_placette, cl))
  
  #pick the plots that have at least one fl belonging to class a, class b and class c each
  tbl2 <- tbl[which(tbl[,1]>0 & tbl[,2]>0 & tbl[,3]>0),]
  
  #subset those plots and their metrics which satisfy the above condition
  pmetsfl1.mix.withallcls <- pmetsfl1.mix.all[id_placette %in% rownames(tbl2)]
  
  #subset those plots and their metrics which do not satisfy the previous condition
  pmetsfl1.mix.woallcls <- pmetsfl1.mix.all[!id_placette %in% rownames(tbl2)]
  
  
  
  
  
  #pick only those fls and their metrics which belong to class a
  pmetsfl1.mix.cla <- pmetsfl1.mix.withallcls[cl=="a"]
  
  #combine with the remaining plots
  pmetsfl1.mix.cla <- rbind(pmetsfl1.mix.cla, pmetsfl1.mix.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.mix.cla <- unique(pmetsfl1.mix.cla[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.mix.cla <- pmetsfl1.mix.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot 
  pmetsfl1.mix.cla <- pmetsfl1.mix.cla[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.mix.clb <- pmetsfl1.mix.withallcls[cl=="b"]
  
  #combine with the remaining plots
  pmetsfl1.mix.clb <- rbind(pmetsfl1.mix.clb, pmetsfl1.mix.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.mix.clb <- unique(pmetsfl1.mix.clb[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.mix.clb <- pmetsfl1.mix.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot
  pmetsfl1.mix.clb <- pmetsfl1.mix.clb[, wt := wt/sum(wt), id_placette]
  
  
  
  
  
  #pick only those pcs and their metrics which belong to class b
  pmetsfl1.mix.clc <- pmetsfl1.mix.withallcls[cl=="c"]
  
  #combine with the remaining plots
  pmetsfl1.mix.clc <- rbind(pmetsfl1.mix.clc, pmetsfl1.mix.woallcls)
  
  # recompute the probability of picking a class of data when a flight line is picked, per plot
  pmetsfl1.mix.clc <- unique(pmetsfl1.mix.clc[, prob := prop.table(table(cl))[cl], id_placette][])
  # compute the inverse probability 
  pmetsfl1.mix.clc <- pmetsfl1.mix.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  # compute weights based on the inverse probability
  # these weights ensure that at the end of several iterations all classes invariably picked roughly the same number of times
  # these weights come into action later on when sample() is used to subset one row per plot
  pmetsfl1.mix.clc <- pmetsfl1.mix.clc[, wt := wt/sum(wt), id_placette]
  
  
}

###########################################################################################################
##Generate samples######
########################
dbase <- pmetsfl1.mix.all
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
smplstfl1.mix.all <- foreach(i = 1:10000, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  set.seed(i)
  dbase.sset <- dbase[, .SD[sample(.N, min(1,.N), prob = wt)], by = id_placette]
  return(which(dbase$cvladlidr%in%dbase.sset$cvladlidr & dbase$pflidr%in%dbase.sset$pflidr))
  # return(sset)
}
stopCluster(clus)
registerDoSEQ()
smplstfl1.mix.all <- matrix(unlist(smplstfl1.mix.all), nrow = length(unique(dbase$id_placette)))
smplstfl1.mix.all <- unique(as.data.table(t(smplstfl1.mix.all)))
###########################################################################################################
time_log <- data.frame()
########################################

start <- Sys.time()
dbase <- pmetsfl1.mix.all
fds <- bauges.db.mix
idx.lst <- smplstfl1.mix.all

f1l <- log(G75)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(G75)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


##Simulations
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
n <- 5000
baugesfl1.mix.all2 <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  idx <- as.vector(unlist(idx.lst[i]))
  mets_for_model <- dbase[idx]
  # dbase$id_placette <- as.factor(dbase$id_placette)
  setkey(mets_for_model,"id_placette")
  setkey(fds,"id_placette")
  mets_for_model <- fds[mets_for_model]
  mets_for_model <- mets_for_model[id_placette!="268_74_ONF"]
  mets_for_model <- mets_for_model[id_placette!="137_IRSTEA"]
  
  mets1 <- mets_for_model[,c("G75","meanch", "varch", 'pflidr', "cvladlidr")]
  
  # y1 <- log(mets1$G_m2_ha)
  # x1 <- log(mets1$meanch)
  # x2 <- log(mets1$varch)
  # x3 <- log(mets1$pflidr)
  # x4 <- log(mets1$cvladlidr)
  
  # func_linreg <- function(b0 ,b1, b2, b3, b4, sig)
  # {
  #   ypred <- b0+b1*x1+b2*x2+b3*x3+b4*x4
  #   -sum(dnorm(y1, mean = ypred, sd = sig, log=TRUE))
  # }
  # 
  
  m1l <- train(f1l,
               data = mets1,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  
  G.l.coeff <- m1l$finalModel$coefficients
  G.l.pred <- m1l$pred$pred
  G.l.obs <- m1l$pred$obs
  
  
  mets2 <- mets_for_model[,c("G75","meanch", "varch", 'pfsumprof', "cvladvox")]
  m1v <- train(f1v,
               data = mets2,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  G.v.coeff <- m1v$finalModel$coefficients
  G.v.pred <- m1v$pred$pred
  G.v.obs <- m1v$pred$obs
  ##########################################################
  
  ##########################################################
  mets3 <- mets_for_model[,c("volume_total", "meanch", "varch", 'pflidr', "cvladlidr")]
  m2l <- train(f2l,
               data = mets3,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtot.l.coeff <- m2l$finalModel$coefficients
  vtot.l.pred <- m2l$pred$pred
  vtot.l.obs <- m2l$pred$obs
  
  mets4 <- mets_for_model[,c("volume_total", "meanch", "varch", 'pfsumprof', "cvladvox")]
  m2v <- train(f2v,
               data = mets4,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtot.v.coeff <- m2v$finalModel$coefficients
  vtot.v.pred <- m2v$pred$pred
  vtot.v.obs <- m2v$pred$obs
  ##########################################################
  
  ###########################################################
  mets5 <- mets_for_model[,c("volume_tige", "meanch", "varch", 'pflidr', "cvladlidr")]
  m3l <- train(f3l,
               data = mets5,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtig.l.coeff <- m3l$finalModel$coefficients
  vtig.l.pred <- m3l$pred$pred
  vtig.l.obs <- m3l$pred$obs
  
  mets6 <- mets_for_model[,c("volume_tige", "meanch", "varch", 'pfsumprof', "cvladvox")]
  m3v <- train(f3v,
               data = mets6,
               method = "lm",
               trControl = trainControl(method="LOOCV"))
  vtig.v.coeff <- m3v$finalModel$coefficients
  vtig.v.pred <- m3v$pred$pred
  vtig.v.obs <- m3v$pred$obs
  ############################################################
  ############################################################
  
  x <- (list("G.lidr.pred" = G.l.pred,
             "G.lidr.obs" = G.l.obs,
             "G.vox.pred" = G.v.pred,
             "G.vox.obs" = G.v.obs,
             "vtot.lidr.pred" = vtot.l.pred,
             "vtot.lidr.obs" = vtot.l.obs,
             "vtot.vox.pred" = vtot.v.pred,
             "vtot.vox.obs" = vtot.v.obs,
             "vtig.lidr.pred" = vtig.l.pred,
             "vtig.lidr.obs" = vtig.l.obs,
             "vtig.vox.pred" = vtig.v.pred,
             "vtig.vox.obs" = vtig.v.obs,
             "G.l.coeff" = G.l.coeff,
             "G.v.coeff" = G.v.coeff,
             "vtot.l.coeff" = vtot.l.coeff,
             "vtot.v.coeff" = vtot.v.coeff,
             "vtig.l.coeff" = vtig.l.coeff,
             "vtig.v.coeff" = vtig.v.coeff
  ))
  return(x)
}
stopCluster(clus)
registerDoSEQ()

stop <- Sys.time()
time_log <- rbind(time_log, c(i, (stop-start)))
#####################################################################################

# func_mdlmets <- function(obs, pred, for_attr, mettype)
# {
#   yobs <- exp(obs)
#   ypred <- exp(pred)
#   see <- sqrt(sum((obs-pred)^2)/(length(obs)-5))
#   cf <- exp((see^2)/2)
#   ypred <- ypred*cf
#   SSE <- sum((yobs-ypred)^2)
#   SST <- sum((mean(yobs)-yobs)^2)
#   R2 <- 1-(SSE/SST)
#   aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
#   MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
#   RMSE <- sqrt(mean((yobs-ypred)^2))
#   RMSEpc <- RMSE*100/mean(yobs)
#   return(list( "Forest_attr"=for_attr, "Metrics"=mettype, "R2"=aR2, "RMSE"=RMSE,"rRMSE"=RMSEpc,"MPE"=MPE))
# }

baugesfl1.mdlmets.con.all <- melt(rbindlist(lapply(baugesfl1.con.all, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.con.cla <- melt(rbindlist(lapply(baugesfl1.con.cla, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.con.clb <- melt(rbindlist(lapply(baugesfl1.con.clb, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.con.clc <- melt(rbindlist(lapply(baugesfl1.con.clc, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))

baugesfl1.mdlmets.feu.all <- melt(rbindlist(lapply(baugesfl1.feu.all, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.feu.cla <- melt(rbindlist(lapply(baugesfl1.feu.cla, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.feu.clb <- melt(rbindlist(lapply(baugesfl1.feu.clb, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.feu.clc <- melt(rbindlist(lapply(baugesfl1.feu.clc, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))

baugesfl1.mdlmets.mix.all <- melt(rbindlist(lapply(baugesfl1.mix.all, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.mix.cla <- melt(rbindlist(lapply(baugesfl1.mix.cla, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.mix.clb <- melt(rbindlist(lapply(baugesfl1.mix.clb, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))
baugesfl1.mdlmets.mix.clc <- melt(rbindlist(lapply(baugesfl1.mix.clc, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "ref")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "ref")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "ref")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "RMSE", "rRMSE", "MPE"))


baugesfl1.mdlmets.con.all <- cbind(baugesfl1.mdlmets.con.all, "exp"=rep("fl1", nrow(baugesfl1.mdlmets.con.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.con.all <- cbind(baugesfl1.mdlmets.con.all, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.con.all)))
baugesfl1.mdlmets.con.cla <- cbind(baugesfl1.mdlmets.con.cla, "exp"=rep("A", nrow(baugesfl1.mdlmets.con.cla)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.con.cla <- cbind(baugesfl1.mdlmets.con.cla, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.con.cla)))
baugesfl1.mdlmets.con.clb <- cbind(baugesfl1.mdlmets.con.clb, "exp"=rep("B", nrow(baugesfl1.mdlmets.con.clb)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.con.clb <- cbind(baugesfl1.mdlmets.con.clb, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.con.cla)))
baugesfl1.mdlmets.con.clc <- cbind(baugesfl1.mdlmets.con.clc, "exp"=rep("C", nrow(baugesfl1.mdlmets.con.clc)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.con.clc <- cbind(baugesfl1.mdlmets.con.clc, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.con.cla)))


baugesfl1.mdlmets.feu.all <- cbind(baugesfl1.mdlmets.feu.all, "exp"=rep("fl1", nrow(baugesfl1.mdlmets.feu.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.feu.all <- cbind(baugesfl1.mdlmets.feu.all, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.feu.all)))
baugesfl1.mdlmets.feu.cla <- cbind(baugesfl1.mdlmets.feu.cla, "exp"=rep("A", nrow(baugesfl1.mdlmets.feu.cla)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.feu.cla <- cbind(baugesfl1.mdlmets.feu.cla, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.feu.cla)))
baugesfl1.mdlmets.feu.clb <- cbind(baugesfl1.mdlmets.feu.clb, "exp"=rep("B", nrow(baugesfl1.mdlmets.feu.clb)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.feu.clb <- cbind(baugesfl1.mdlmets.feu.clb, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.feu.cla)))
baugesfl1.mdlmets.feu.clc <- cbind(baugesfl1.mdlmets.feu.clc, "exp"=rep("C", nrow(baugesfl1.mdlmets.feu.clc)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.feu.clc <- cbind(baugesfl1.mdlmets.feu.clc, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.feu.cla)))


baugesfl1.mdlmets.mix.all <- cbind(baugesfl1.mdlmets.mix.all, "exp"=rep("fl1", nrow(baugesfl1.mdlmets.mix.all)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.mix.all <- cbind(baugesfl1.mdlmets.mix.all, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.mix.all)))
baugesfl1.mdlmets.mix.cla <- cbind(baugesfl1.mdlmets.mix.cla, "exp"=rep("A", nrow(baugesfl1.mdlmets.mix.cla)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.mix.cla <- cbind(baugesfl1.mdlmets.mix.cla, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.mix.cla)))
baugesfl1.mdlmets.mix.clb <- cbind(baugesfl1.mdlmets.mix.clb, "exp"=rep("B", nrow(baugesfl1.mdlmets.mix.clb)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.mix.clb <- cbind(baugesfl1.mdlmets.mix.clb, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.mix.cla)))
baugesfl1.mdlmets.mix.clc <- cbind(baugesfl1.mdlmets.mix.clc, "exp"=rep("C", nrow(baugesfl1.mdlmets.mix.clc)), "id"=rep(rep(1:5000, 1, each=6), 4))
baugesfl1.mdlmets.mix.clc <- cbind(baugesfl1.mdlmets.mix.clc, "fl"=rep("fl1", nrow(baugesfl1.mdlmets.mix.cla)))
###################################################################################################
baugesfl1.mdlmets.con <- rbind(baugesfl1.mdlmets.con.all,
                               baugesfl1.mdlmets.con.cla, 
                               baugesfl1.mdlmets.con.clb, 
                               baugesfl1.mdlmets.con.clc)

baugesfl1.mdlmets.feu <- rbind(baugesfl1.mdlmets.feu.all,
                               baugesfl1.mdlmets.feu.cla, 
                               baugesfl1.mdlmets.feu.clb, 
                               baugesfl1.mdlmets.feu.clc)

baugesfl1.mdlmets.mix <- rbind(baugesfl1.mdlmets.mix.all,
                               baugesfl1.mdlmets.mix.cla, 
                               baugesfl1.mdlmets.mix.clb, 
                               baugesfl1.mdlmets.mix.clc)

bauges.mdlmets.con <- rbind(baugesfl1.mdlmets.con, 
                            baugesfl2.mdlmets.con, 
                            baugesfl3.mdlmets.con.all)
bauges.mdlmets.con$exp <- factor(bauges.mdlmets.con$exp, levels=c("fl1","A", "B", "C",
                                                                  "fl2", "AB", "AC", "BC", "fl3")) 
                                                                  
bauges.mdlmets.feu <- rbind(baugesfl1.mdlmets.feu, 
                            baugesfl2.mdlmets.feu, 
                            baugesfl3.mdlmets.feu.all)
bauges.mdlmets.feu$exp <- factor(bauges.mdlmets.feu$exp, levels=c("fl1","A", "B", "C",
                                                                  "fl2", "AB", "AC", "BC", "fl3")) 


bauges.mdlmets.mix <- rbind(baugesfl1.mdlmets.mix, 
                            baugesfl2.mdlmets.mix,
                            baugesfl3.mdlmets.mix.all)
bauges.mdlmets.mix$exp <- factor(bauges.mdlmets.mix$exp, levels=c("fl1","A", "B", "C",
                                                                  "fl2", "AB", "AC", "BC", "fl3")) 
bauges.mdlmets.con$forest_type <- rep("Coniferous", nrow(bauges.mdlmets.con))
bauges.mdlmets.feu$forest_type <- rep("Broadleaved", nrow(bauges.mdlmets.feu))
bauges.mdlmets.mix$forest_type <- rep("Mixed", nrow(bauges.mdlmets.mix))

bauges.mdlmets <- rbind(bauges.mdlmets.con, bauges.mdlmets.feu, bauges.mdlmets.mix)

ggplot(data=bauges.mdlmets.con[variable!="RMSE" & 
                             Forest_attr=="Basal area"], aes(fill=Metrics, y=value, x=exp))+
  geom_boxplot()+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_minimal()+
  labs(x="Experiments", y="")+
  scale_fill_grey(start = 1, end = 0.5)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")




##Result 2
ggplot(data=bauges.mdlmets.mix[variable!="RMSE"& 
                                 Forest_attr=="Basal area" & 
                                 exp %in% c("fl1", "fl2", "fl3")&
                                 Metrics=="ref"], aes(x=exp, y=value, fill=Metrics))+
  geom_boxplot()+
  facet_wrap(Forest_attr~variable, scales = "free")+
  theme_minimal()+
  scale_fill_grey(start = 1, end = 0.5)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")


##Result 3
ggplot(data=bauges.mdlmets[variable!="RMSE" & 
                                 Forest_attr=="Basal area" &
                                 exp %in% c("fl1", "fl2", "fl3") &
                              Metrics=="ref"], aes(fill=Metrics, y=value, x=exp))+
  geom_boxplot()+
  facet_wrap(forest_type~variable, scales = "free")+
  theme_minimal()+
  labs(x="Experiments", y="")+
  scale_fill_grey(start = 1, end = 0.5)+
  theme(text=element_text(family="serif", size=9*(96/72)),
        legend.position = "top")




