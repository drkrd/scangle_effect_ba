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

# l1 <- ifelse(typefor=="conifere","con",
#              ifelse(typefor=="feuillus","feu","mix"))
# l2 <- ifelse(typepc=="flightlines","fl","all")


####IMPORTANT#####
##Here, read only NORMALISED POINT CLOUDS with label format "plotid_n.las"
# plots <- readLAScatalog(paste0("D:/1_Work/5_Bauges/Data/ULM/LAS/norm/plots/", typefor, "/15m_rad/", typepc, "/"))
fd_smry <- fread("D:/1_Work/__R_codes/Projects/scangle_effect_ba/data/bauges_db_new.csv", sep = ",")
colnames(fd_smry)[2] <- "id_placette"


fd_smry <- fd_smry[, newstratum := ifelse((G175R/G175)>0.75 & (G175F/G175)<0.25, "Resineux", 
                                          ifelse((G175F/G175)>0.75 & (G175R/G175)<0.25, "Feuillus", "Mixte"))]
fd_smry[, stratum := as.factor(stratum)]
fd_smry[, newstratum := as.factor(newstratum)]
fd_smry[, id_placette := as.factor(id_placette)]
fd_smry <- fd_smry[which(id_placette %in% plot)]
setkey(fd_smry, "id_placette")


# The following flightlines are problematic. There is no voxelisation for them
# 283_74_ONF_un@6.28
# 96_IRSTEA_un@30

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
  meanch <- func_meanch(las@data$Z, las@data$ReturnNumber)
  varch <- func_varch(las@data$Z, las@data$ReturnNumber)
  pflidr <- func_pf(las@data$Z, las@data$ReturnNumber, ht=5)
  cvladlidr <- func_cvlad(las@data$Z, las@data$ReturnNumber, ht=5)
  return(list( 
    id_placette = id_plac,
    meanang = mang,
    meanch = meanch,
    varch = varch,
    pflidr = pflidr,
    cvladlidr = cvladlidr))
}
plotmetsfl1 <- catalog_apply(plots, func_computeall)
plotmetsfl1 <- rbindlist(plotmetsfl1)
plotmetsfl1 <- plotmetsfl1[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                       ifelse(meanang>=10&meanang<20, "b",
                                              ifelse(meanang>=20&meanang<30, "c",
                                                     ifelse(meanang>=30&meanang<40,"d","e"))))]


plotmetsfl1 <- plotmetsfl1[cl!="d" & cl!="e"]
plotmetsfl1 <- na.omit(plotmetsfl1)
pmetsfl1 <- plotmetsfl1
pmetsfl1 <- pmetsfl1[-which(pmetsfl1$id_placette=="283_74_ONF" & pmetsfl1$meanang==6.28)]



fd_smry <- fd_smry[which(id_placette %in% plotmetsfl1$id_placette), 
                   c("id_placette", "G75", "volume_total", "volume_tige", "newstratum")]
fd_smry_con <- fd_smry[newstratum=="Resineux"]
fd_smry_feu <- fd_smry[newstratum=="Feuillus"]
fd_smry_mix <- fd_smry[newstratum=="Mixte"]
##############################################################################################################################
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


voxall <- rbindlist(apply(allvoxfiles, 1, 
                          func_normvox2, pth ="D:/1_Work/5_Bauges/Data/ULM/LAS/unnorm/plots/15m_rad/may2021/allpoints/", 
                          ht=5))
setDT(voxall)
pfcvladvox <- voxall[, .(cvladvox=cv(m, na.rm = TRUE), 
                         sdvfp=sqrt(sum(m*(k1-(sum(k1*m)/sum(m)))^2)/(sum(m)*(length(m[which(m!=0)])-1)/length(m[which(m!=0)]))),
                         pfsumprof=exp(-0.5*sum(m, na.rm = TRUE))), by=.(id_placette, meanang)]


pfcvladvox$meanang <- as.numeric(pfcvladvox$meanang)
pfcvladvox <- pfcvladvox[,cl:=ifelse(meanang>=0&meanang<10, "a",
                                     ifelse(meanang>=10&meanang<20, "b",
                                            ifelse(meanang>=20&meanang<30, "c",
                                                   ifelse(meanang>=30&meanang<40,"d","e"))))]
pfcvladvox <- pfcvladvox[cl!="d" & cl!="e"]
setkeyv(pfcvladvox, c("id_placette", "meanang"))
setkeyv(pmetsfl1, c("id_placette", "meanang"))
pmetsfl1 <- pmetsfl1[pfcvladvox]
pmetsfl1 <- pmetsfl1[which(id_placette %in% fd_smry$id_placette)]


######################################
#############All################
#############################################################################################################################
{
  plotmetsfl1 <- unique(plotmetsfl1[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1 <- plotmetsfl1[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1 <- plotmetsfl1[, wt := wt/sum(wt), id_placette]
  # plotmetsfl1$meanang <- as.factor(plotmetsfl1$meanang)
  plotmetsfl1 <- plotmetsfl1[pflidr!=0]
  pmetsfl1.all <- plotmetsfl1
}

func_class.ssets <- function(pmets, cls)
{
  #tabulate per plot the number of pcs belonging to each class
  tbl <- with(pmets, table(id_placette, cl))
  
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
  pmetsfl1.con.all <- pmetsfl1[which(id_placette %in% fd_smry_con$id_placette)]

  pmetsfl1.con.cla <- func_class.ssets(pmetsfl1.con.all, "a")
  pmetsfl1.con.cla <- unique(pmetsfl1.con.cla[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.con.cla <- pmetsfl1.con.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.con.cla <- pmetsfl1.con.cla[, wt := wt/sum(wt), id_placette]
  
  pmetsfl1.con.clb <- func_class.ssets(pmetsfl1.con.all, "b")
  pmetsfl1.con.clb <- unique(pmetsfl1.con.clb[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.con.clb <- pmetsfl1.con.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.con.clb <- pmetsfl1.con.clb[, wt := wt/sum(wt), id_placette]
  
  pmetsfl1.con.clc <- func_class.ssets(pmetsfl1.con.all, "c")
  pmetsfl1.con.clc <- unique(pmetsfl1.con.clc[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.con.clc <- pmetsfl1.con.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.con.clc <- pmetsfl1.con.clc[, wt := wt/sum(wt), id_placette]
}

{
  pmetsfl1.con.all <- pmetsfl1[which(id_placette %in% fd_smry_con$id_placette)]
  pmetsfl1.con.all <- pmetsfl1.con.all[pflidr!=0]
  pmetsfl1.con.all <- unique(pmetsfl1.con.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.con.all <- pmetsfl1.con.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.con.all <- pmetsfl1.con.all[, wt := wt/sum(wt), id_placette]
  pmetsfl1.con.all$meanang <- as.factor(pmetsfl1.con.all$meanang)
  ##############################################################################################################################
  
  plotmetsfl1a_con=plotmetsfl1_con[plotmetsfl1_con[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
  plotmetsfl1b_con=plotmetsfl1_con[plotmetsfl1_con[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
  plotmetsfl1c_con=plotmetsfl1_con[plotmetsfl1_con[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]
  
  l1_con <- unique(plotmetsfl1a_con[cl=="a"]$id_placette)
  l2_con <- unique(plotmetsfl1b_con[cl=="b"]$id_placette)
  l3_con <- unique(plotmetsfl1c_con[cl=="c"]$id_placette)
  
  cps <- intersect(intersect(l1_con, l2_con), l3_con)
  
  func_dffilter <- function(df, cls)
  {
    df1 <- df[id_placette %in% cps]
    df1 <- df1[df1[, .I[cl==cls  | all(cl!=cls)], by = id_placette]$V1]
    df2 <- df[!id_placette %in% cps]
    dfn <- rbind(df1, df2)
    return(dfn)
  }
  
  
  plotmetsfl1a_con <- func_dffilter(plotmetsfl1_con, "a")
  plotmetsfl1a_con <- unique(plotmetsfl1a_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1a_con <- plotmetsfl1a_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1a_con <- plotmetsfl1a_con[, wt := wt/sum(wt), id_placette]
  
  
  plotmetsfl1b_con <- func_dffilter(plotmetsfl1_con, "b")
  plotmetsfl1b_con <- unique(plotmetsfl1b_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1b_con <- plotmetsfl1b_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1b_con <- plotmetsfl1b_con[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1c_con <- func_dffilter(plotmetsfl1_con, "c")
  plotmetsfl1c_con <- unique(plotmetsfl1c_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1c_con <- plotmetsfl1c_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1c_con <- plotmetsfl1c_con[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1a_con=plotmetsfl1_con[plotmetsfl1_con[, all(cl != 'a')| (cl == 'a' & .BY %in% cps)|!.BY %in% cps, 
                                                   by = id_placette]$V1]

  plotmetsfl1a_con <- unique(plotmetsfl1a_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1a_con <- plotmetsfl1a_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1a_con <- plotmetsfl1a_con[, wt := wt/sum(wt), id_placette]
  
  
  
  plotmetsfl1b_con=plotmetsfl1_con[plotmetsfl1_con[, all(cl != 'b')| (cl == 'b' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  
  plotmetsfl1b_con <- unique(plotmetsfl1b_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1b_con <- plotmetsfl1b_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1b_con <- plotmetsfl1b_con[, wt := wt/sum(wt), id_placette]
  
  
  
  plotmetsfl1c_con=plotmetsfl1_con[plotmetsfl1_con[, all(cl != 'c')| (cl == 'c' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  
  plotmetsfl1c_con <- unique(plotmetsfl1c_con[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1c_con <- plotmetsfl1c_con[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1c_con <- plotmetsfl1c_con[, wt := wt/sum(wt), id_placette]
  #-------------------------------------------------------------------------------------------------------------------------#
  
}


######################################
#############Feuillus#################
###########################################################################################################################

{
  pmetsfl1.feu.all <- pmetsfl1[which(id_placette %in% fd_smry_feu$id_placette)]
  
  pmetsfl1.feu.cla <- func_class.ssets(pmetsfl1.feu.all, "a")
  pmetsfl1.feu.cla <- unique(pmetsfl1.feu.cla[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.feu.cla <- pmetsfl1.feu.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.feu.cla <- pmetsfl1.feu.cla[, wt := wt/sum(wt), id_placette]
  
  pmetsfl1.feu.clb <- func_class.ssets(pmetsfl1.feu.all, "b")
  pmetsfl1.feu.clb <- unique(pmetsfl1.feu.clb[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.feu.clb <- pmetsfl1.feu.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.feu.clb <- pmetsfl1.feu.clb[, wt := wt/sum(wt), id_placette]
  
  pmetsfl1.feu.clc <- func_class.ssets(pmetsfl1.feu.all, "c")
  pmetsfl1.feu.clc <- unique(pmetsfl1.feu.clc[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.feu.clc <- pmetsfl1.feu.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.feu.clc <- pmetsfl1.feu.clc[, wt := wt/sum(wt), id_placette]
}

{
  pmetsfl1.feu.all <- plotmetsfl1[which(id_placette %in% fd_smry_feu$id_placette)]
  pmetsfl1.feu.all <- pmetsfl1.feu.all[pflidr!=0]
  pmetsfl1.feu.all <- unique(pmetsfl1.feu.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.feu.all <- pmetsfl1.feu.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.feu.all <- pmetsfl1.feu.all[, wt := wt/sum(wt), id_placette]
  pmetsfl1.feu.all$meanang <- as.factor(pmetsfl1.feu.all$meanang)
  ###########################################################################################################################
  
  plotmetsfl1a_feu=plotmetsfl1_feu[plotmetsfl1_feu[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
  plotmetsfl1b_feu=plotmetsfl1_feu[plotmetsfl1_feu[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
  plotmetsfl1c_feu=plotmetsfl1_feu[plotmetsfl1_feu[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]
  
  l1_feu <- unique(plotmetsfl1a_feu[cl=="a"]$id_placette)
  l2_feu <- unique(plotmetsfl1b_feu[cl=="b"]$id_placette)
  l3_feu <- unique(plotmetsfl1c_feu[cl=="c"]$id_placette)
  
  cps <- intersect(intersect(l1_feu, l2_feu), l3_feu)
  
  plotmetsfl1a_feu=plotmetsfl1[plotmetsfl1[, all(cl != 'a')| (cl == 'a' & .BY %in% cps)|!.BY %in% cps, 
                                           by = id_placette]$V1]
  plotmetsfl1a_feu <- unique(plotmetsfl1a_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1a_feu <- plotmetsfl1a_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1a_feu <- plotmetsfl1a_feu[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1b_feu=plotmetsfl1[plotmetsfl1[, all(cl != 'b')| (cl == 'b' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  plotmetsfl1b_feu <- unique(plotmetsfl1b_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1b_feu <- plotmetsfl1b_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1b_feu <- plotmetsfl1b_feu[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1c_feu=plotmetsfl1[plotmetsfl1[, all(cl != 'c')| (cl == 'c' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  plotmetsfl1c_feu <- unique(plotmetsfl1c_feu[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1c_feu <- plotmetsfl1c_feu[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1c_feu <- plotmetsfl1c_feu[, wt := wt/sum(wt), id_placette]
  #-----------------------------------------------------------------------------------------------------------------------#
  
}

######################################
################Mixte#################
###########################################################################################################################

{
  pmetsfl1.mix.all <- pmetsfl1[which(id_placette %in% fd_smry_mix$id_placette)]
  pmetsfl1.mix.all <- pmetsfl1.mix.all[pflidr!=0]
  pmetsfl1.mix.all <- unique(pmetsfl1.mix.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.mix.all <- pmetsfl1.mix.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.mix.all <- pmetsfl1.mix.all[, wt := wt/sum(wt), id_placette]
  
  pmetsfl1.mix.cla <- func_class.ssets(pmetsfl1.mix.all, "a")
  pmetsfl1.mix.cla <- unique(pmetsfl1.mix.cla[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.mix.cla <- pmetsfl1.mix.cla[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.mix.cla <- pmetsfl1.mix.cla[, wt := wt/sum(wt), id_placette]
  
  pmetsfl1.mix.clb <- func_class.ssets(pmetsfl1.mix.all, "b")
  pmetsfl1.mix.clb <- unique(pmetsfl1.mix.clb[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.mix.clb <- pmetsfl1.mix.clb[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.mix.clb <- pmetsfl1.mix.clb[, wt := wt/sum(wt), id_placette]
  
  pmetsfl1.mix.clc <- func_class.ssets(pmetsfl1.mix.all, "c")
  pmetsfl1.mix.clc <- unique(pmetsfl1.mix.clc[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.mix.clc <- pmetsfl1.mix.clc[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.mix.clc <- pmetsfl1.mix.clc[, wt := wt/sum(wt), id_placette]
}

{
  pmetsfl1.mix.all <- plotmetsfl1[which(id_placette %in% fd_smry_mix$id_placette)]
  pmetsfl1.mix.all <- pmetsfl1.mix.all[pflidr!=0]
  pmetsfl1.mix.all <- unique(pmetsfl1.mix.all[, prob := prop.table(table(cl))[cl], id_placette][])
  pmetsfl1.mix.all <- pmetsfl1.mix.all[, wt := (1/prob)/(1/sum(prob)), id_placette]
  pmetsfl1.mix.all <- pmetsfl1.mix.all[, wt := wt/sum(wt), id_placette]
  pmetsfl1.mix.all$meanang <- as.factor(pmetsfl1.mix.all$meanang)
  
  ###########################################################################################################################
  
  plotmetsfl1a_mix=pmetsfl1.mix.all[pmetsfl1.mix.all[, .I[cl=="a"  | all(cl!="a")], by = id_placette]$V1]
  plotmetsfl1b_mix=pmetsfl1.mix.all[pmetsfl1.mix.all[, .I[cl=="b"  | all(cl!="b")], by = id_placette]$V1]
  plotmetsfl1c_mix=pmetsfl1.mix.all[pmetsfl1.mix.all[, .I[cl=="c"  | all(cl!="c")], by = id_placette]$V1]
  
  l1_mix <- unique(plotmetsfl1a_mix[cl=="a"]$id_placette)
  l2_mix <- unique(plotmetsfl1b_mix[cl=="b"]$id_placette)
  l3_mix <- unique(plotmetsfl1c_mix[cl=="c"]$id_placette)
  
  cps <- intersect(intersect(l1_mix, l2_mix), l3_mix)
  
  plotmetsfl1a_mix=plotmetsfl1[plotmetsfl1[, all(cl != 'a')| (cl == 'a' & .BY %in% cps)|!.BY %in% cps, 
                                           by = id_placette]$V1]
  plotmetsfl1a_mix <- unique(plotmetsfl1a_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1a_mix <- plotmetsfl1a_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1a_mix <- plotmetsfl1a_mix[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1b_mix=plotmetsfl1[plotmetsfl1[, all(cl != 'b')| (cl == 'b' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  plotmetsfl1b_mix <- unique(plotmetsfl1b_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1b_mix <- plotmetsfl1b_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1b_mix <- plotmetsfl1b_mix[, wt := wt/sum(wt), id_placette]
  
  plotmetsfl1c_mix=plotmetsfl1[plotmetsfl1[, all(cl != 'c')| (cl == 'c' & .BY %in% cps)|!.BY %in% cps,
                                           by = id_placette]$V1]
  plotmetsfl1c_mix <- unique(plotmetsfl1c_mix[, prob := prop.table(table(cl))[cl], id_placette][])
  plotmetsfl1c_mix <- plotmetsfl1c_mix[, wt := (1/prob)/(1/sum(prob)), id_placette]
  plotmetsfl1c_mix <- plotmetsfl1c_mix[, wt := wt/sum(wt), id_placette]
  #-----------------------------------------------------------------------------------------------------------------------#
  
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
smplstfl1.mix.all <- matrix(unlist(smplstfl1.mix.all), nrow = 47)
smplstfl1.mix.all <- unique(as.data.table(t(smplstfl1.mix.all)))







time_log <- data.frame()




start <- Sys.time()
dbase <- pmetsfl1.mix.all
fds <- fd_smry_mix
idx.lst <- smplstfl1.mix.all

f1l <- log(G75)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f1v <- log(G75)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f2l <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f2v <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
f3l <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
f3v <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)


##Simulations for basal area
clus <- makeCluster(detectCores() - 1)
registerDoParallel(clus, cores = detectCores() - 1)
n <- 500
baugesfl1.mix.all1 <- foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
  idx <- as.vector(unlist(idx.lst[i]))
  mets_for_model <- dbase[idx]
  setkey(mets_for_model,"id_placette")
  mets_for_model <- fds[mets_for_model]
  
  m1l <- train(f1l,
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  G.l.coeff <- m1l$finalm1l$coefficients
  G.l.pred <- m1l$pred$pred
  G.l.obs <- m1l$pred$obs
  
  
  
  m1v <- train(f1v,
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  G.v.coeff <- m1v$finalm1v$coefficients
  G.v.pred <- m1v$pred$pred
  G.v.obs <- m1v$pred$obs
  ##########################################################
  
  ##########################################################
  m2l <- train(f2l,
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  vtot.l.coeff <- m2l$finalm2l$coefficients
  vtot.l.pred <- m2l$pred$pred
  vtot.l.obs <- m2l$pred$obs
  
  m2v <- train(f2v,
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  vtot.v.coeff <- m2v$finalm2v$coefficients
  vtot.v.pred <- m2v$pred$pred
  vtot.v.obs <- m2v$pred$obs
  ##########################################################
  
  ###########################################################
  m3l <- train(f3l,
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  vtig.l.coeff <- m3l$finalm3l$coefficients
  vtig.l.pred <- m3l$pred$pred
  vtig.l.obs <- m3l$pred$obs
  
  m3v <- train(f3v,
                 data = mets_for_model,
                 method = "lm",
                 trControl = trainControl(method="LOOCV"))
  vtig.v.coeff <- m3v$finalm3v$coefficients
  vtig.v.pred <- m3v$pred$pred
  vtig.v.obs <- m3v$pred$obs
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
time_log <- rbind(time_log, c(n, (stop-start)))



func_mdlmets <- function(obs, pred, for_attr, mettype)
{
  yobs <- exp(obs)
  ypred <- exp(pred)
  see <- sqrt(sum((obs-pred)^2)/(length(obs)-4))
  cf <- exp((see^2)/2)
  ypred <- ypred*cf
  R2 <- cor(ypred, yobs)^2
  aR2 <- 1-((1-R2)*(length(ypred)-1)/(length(ypred)-4-1))
  MPE <- (100/length(ypred))*sum((yobs-ypred)/yobs)
  RMSE <- sqrt(mean((yobs-ypred)^2))
  MAE <- mean(abs(ypred-yobs))
  return(list( "Forest_attr"=for_attr, "Metrics"=mettype, "R2"=aR2, "MPE"=MPE, "RMSE"=RMSE,"MAE"=MAE))
}


bgs.mdlmetsfl1.mix.all <- melt(rbindlist(lapply(baugesfl1.mix.all1, function(x)
{
  G.l.mdlmets <- func_mdlmets(x$G.lidr.obs, x$G.lidr.pred, "Basal area", "old")
  G.v.mdlmets <- func_mdlmets(x$G.vox.obs, x$G.vox.pred, "Basal area", "vox")
  vtot.l.mdlmets <- func_mdlmets(x$vtot.lidr.obs, x$vtot.lidr.pred, "Total volume", "old")
  vtot.v.mdlmets <- func_mdlmets(x$vtot.vox.obs, x$vtot.vox.pred, "Total volume", "vox")
  vtig.l.mdlmets <- func_mdlmets(x$vtig.lidr.obs, x$vtig.lidr.pred, "Stem volume", "old")
  vtig.v.mdlmets <- func_mdlmets(x$vtig.vox.obs, x$vtig.vox.pred, "Stem volume", "vox")
  y <- list(G.l.mdlmets,
            G.v.mdlmets,
            vtot.l.mdlmets,
            vtot.v.mdlmets,
            vtig.l.mdlmets,
            vtig.v.mdlmets)
  metdf <- rbindlist(y)
  return(metdf)
})), measure.vars = c("R2", "MPE", "RMSE", "MAE"))


bgs.mdlmetsfl1.mix.all <- cbind(bgs.mdlmetsfl1.mix.all, "exp"=rep("all", nrow(bgs.mdlmetsfl1.mix.all)), "id"=rep(rep(1:500, 1, each=6), 4))

bgs.mdlmetsfl1.mix <- as.data.table(rbind(bgs.mdlmetsfl1.mix.all))

ggplot(data=bgs.mdlmetsfl1.mix, aes(y=value, colour=Metrics))+
  geom_boxplot()+
  facet_grid(variable~Forest_attr, scales = "free")+
  theme_base()+
  scale_fill_grey()



idx <- as.vector(unlist(idx.lst[1]))
mets_for_model <- dbase[idx]
setkey(mets_for_model,"id_placette")
mets_for_model <- fds[mets_for_model]

m1l <- train(f1l,
             data = mets_for_model,
             method = "lm",
             trControl = trainControl(method="LOOCV"))




















# library(boot)
# func_boot <- function(data, idx, fd)
# {
#   mets_for_model <- data[idx]
# 
#   setkey(mets_for_model,"id_placette")
# 
# 
#   mets_for_model <- fd[mets_for_model]
# 
# 
#   model <- train(log(G175)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladlidr),
#                  data = mets_for_model,
#                  method = "lm",
#                  trControl = trainControl(method="LOOCV"))
#   model$results$Rsquared
# }
# 
# boot(plotmetsfl1_con, statistic = func_boot, R=1000, strata = plotmetsfl1_con$id_placette, fd=fd_smry_con)
# 


start <- Sys.time()
type <- "vox"
dbase <- plotmetsfl1_feu
fds <- fd_smry_feu
n1 <- "x"
n2 <- "feu"
n <- 3500
func_extractmdlmets <- function(x, type, exp, mets)
{
  ind <- x[["index"]]
  # dbase <- plotmetsfl1[ind]
  # dbase <- dbase[!id_placette%in%c(21, 22, 26)]
  # ent <- func_entropy(dbase$cl)
  r2 <- x[["modelmets"]]$R2 
  rmse <- x[["modelmets"]]$RMSE
  mae <- x[["modelmets"]]$MAE
  mpe <- x[["modelmets"]]$MPE
  type <- type
  exp <- exp
  mets <- mets
  # names(tbl[which(tbl==max(tbl))])
  return(list("r2"=r2, "rmse"=rmse, "mae"=mae, "mpe"=mpe, "type"=type, "exp"=exp, "mets"=mets))
}
if(type=="lidr")
{
  f1 <- log(G75)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f2 <- log(volume_total)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  f3 <- log(volume_tige)~log(meanch)+log(varch)+log(pflidr)+log(cvladlidr)
  start <- Sys.time()
  ##Simulations for basal area
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme1a <- paste0("mdls_gfl1", n1, n2)
  assign(nme1a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f1,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme1b <- paste0("mdlmets_gfl1", n1, n2)
  assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("BA done")
  
  
  ##Simulations for total volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme2a <- paste0("mdls_vtotfl1", n1, n2)
  assign(nme2a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f2,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme2b <- paste0("mdlmets_vtotfl1", n1, n2)
  assign(nme2b, unique(rbindlist(lapply(get(nme2a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("Vtot done")
  
  ##Simulations for stem volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme3a <- paste0("mdls_vtigfl1", n1, n2)
  assign(nme3a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f3,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme3b <- paste0("mdlmets_vtigfl1", n1, n2)
  assign(nme3b, unique(rbindlist(lapply(get(nme3a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="old"))))
  print("Vtig done")
  stop <- Sys.time()
}else
{
  f1 <- log(G75)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f2 <- log(volume_total)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  f3 <- log(volume_tige)~log(meanch)+log(varch)+log(pfsumprof)+log(cvladvox)
  start <- Sys.time()
  ##Simulations for basal area
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme1a <- paste0("mdls_gfl1vx", n1, n2)
  assign(nme1a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f1,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme1b <- paste0("mdlmets_gfl1vx", n1, n2)
  assign(nme1b, unique(rbindlist(lapply(get(nme1a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="vox"))))
  print("BA done")
  
  
  ##Simulations for total volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme2a <- paste0("mdls_vtotfl1vx", n1, n2)
  assign(nme2a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f2,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme2b <- paste0("mdlmets_vtotfl1vx", n1, n2)
  assign(nme2b, unique(rbindlist(lapply(get(nme2a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="vox"))))
  print("Vtot done")
  
  
  ##Simulations for stem volume
  set.seed(123)
  clus <- makeCluster(detectCores() - 1)
  registerDoParallel(clus, cores = detectCores() - 1)
  nme3a <- paste0("mdls_vtigfl1vx", n1, n2)
  assign(nme3a, foreach(i = 1:n, .packages=c("dplyr", "data.table", "caret")) %dopar% {
    mets_for_model <- dbase[, .SD[sample(.N, min(1,.N), prob=wt) ], by = id_placette]
    setkey(mets_for_model,"id_placette")
    mets_for_model <- fds[mets_for_model]
    
    
    model <- train(f3,
                   data = mets_for_model,
                   method = "lm",
                   trControl = trainControl(method="LOOCV"))
    
    ypred <- exp(model$pred$pred)
    yobs <- exp(model$pred$obs)
    n <- length(yobs)
    MPE <- (100/n)*sum((yobs-ypred)/yobs)
    RMSE <- sqrt(mean((yobs-ypred)^2))
    MAE <- mean(abs(ypred-yobs))
    
    
    x <- (list("modelmets" = list("R2" = model$results$Rsquared, 
                                  "RMSE" = RMSE,
                                  "MAE" = MAE,
                                  "MPE" = MPE),
               "coeffs" = list(model$finalModel$coefficients),
               "index" = which(dbase$meanang %in% mets_for_model$meanang)))
    return(x)
  })
  stopCluster(clus)
  registerDoSEQ()
  nme3b <- paste0("mdlmets_vtigfl1vx", n1, n2)
  assign(nme3b, unique(rbindlist(lapply(get(nme3a), 
                                        func_extractmdlmets, 
                                        type="One flight line", 
                                        exp="allclasses", 
                                        mets="vox"))))
  print("Vtig done")
  stop <- Sys.time()
}


mdlmets_vtotfl1_con  <- unique(rbindlist(lapply(mdls_vtotfl1_con, 
                                                    func_extractmdlmets, 
                                                    type="One flight line", 
                                                    exp="allclasses", 
                                                    mets="old")))


mdlmets_vtotfl1vx_con  <- unique(rbindlist(lapply(mdls_vtotfl1vx_con, 
                                                  func_extractmdlmets, 
                                                  type="One flight line", 
                                                  exp="allclasses", 
                                                  mets="vox")))


mdlmets_gfl1vxcfeu$exp <- rep("mostly c", nrow(mdlmets_gfl1vxcfeu)) 


mdlmets_gx <- rbind(mdlmets_gfl1xfeu, 
                   mdlmets_gfl1vxxfeu)
mdlmets_gx <- melt(mdlmets_gx, id.vars = c("type", "exp", "mets"))
mdlmets_gx$for_attr <- rep("Basal area", nrow(mdlmets_gx))


mdlmets_g <- rbind(mdlmets_gfl1con, 
                   mdlmets_gfl1vxcon,
                   mdlmets_gfl1acon,
                   mdlmets_gfl1vxbcon,
                   mdlmets_gfl1bcon,
                   mdlmets_gfl1vxacon,
                   mdlmets_gfl1ccon,
                   mdlmets_gfl1vxccon)
mdlmets_g <- melt(mdlmets_g, id.vars = c("type", "exp", "mets"))
mdlmets_g$for_attr <- rep("Basal area", nrow(mdlmets_g))

mdlmets_vtot <- rbind(mdlmets_vtotfl1feu, 
                      mdlmets_vtotfl1vxfeu,
                      mdlmets_vtotfl1afeu,
                      mdlmets_vtotfl1vxafeu,
                      mdlmets_vtotfl1bfeu,
                      mdlmets_vtotfl1vxbfeu,
                      mdlmets_vtotfl1cfeu,
                      mdlmets_vtotfl1vxcfeu)
mdlmets_vtot <- melt(mdlmets_vtot, id.vars = c("type", "exp", "mets"))
mdlmets_vtot$for_attr <- rep("Total volume", nrow(mdlmets_vtot))

mdlmets_vtig <- rbind(mdlmets_vtigfl1feu, 
                      mdlmets_vtigfl1vxfeu,
                      mdlmets_vtigfl1afeu,
                      mdlmets_vtigfl1vxafeu,
                      mdlmets_vtigfl1bfeu,
                      mdlmets_vtigfl1vxbfeu,
                      mdlmets_vtigfl1cfeu,
                      mdlmets_vtigfl1vxcfeu)
mdlmets_vtig <- melt(mdlmets_vtig, id.vars = c("type", "exp", "mets"))
mdlmets_vtig$for_attr <- rep("Stem volume", nrow(mdlmets_vtig))


mdlmets_all <- rbind(mdlmets_g, mdlmets_vtot, mdlmets_vtig)


ggplot(data=mdlmets_gx, aes(x=exp, y=value, fill=mets))+
  geom_boxplot()+
  labs(title = "")+
  theme_base(base_size = 20)+
  facet_grid(variable~for_attr, scales = "free")+
  geom_hline(data = allmdlmets, aes(yintercept=value, linetype=type, colour=mets), size=1.5)+
  scale_fill_pander()
  
  




plots <- list()

for(id_plac in unique(voxall$id_placette))
{
  temp <- voxall[id_placette==id_plac]
  temp <- temp[,c("k1", "m", "meanang", "id_placette")]
  plots[[id_plac]] <- ggplot(data=temp, aes(x=k1, y=m, linetype = meanang))+
    theme_minimal()+
    geom_line(size=0.8)+
    coord_flip()+
    labs(title = unique(temp$id_placette))+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
}
plots_con <- plots[names(plots)%in%fd_smry_con$id_placette]
plots_feu <- plots[names(plots)%in%fd_smry_feu$id_placette]
plots_mix <- plots[names(plots)%in%fd_smry_mix$id_placette]


l1 <- c(1,2,3)
l2 <- c(1,4,6)

