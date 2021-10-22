nm <- NULL
lst <- NULL
for(name in names(alldtms)){
  dtm <- alldtms[[name]]
  crs(dtm) <- 2154
  slp <- terrain(dtm, c("slope"), unit="degrees")
  lst <- c(lst,mean(na.omit(as.data.frame(slp)$slope)))
  nm <- c(nm,name)
}
slope.df <- as.data.table(cbind(nm,lst))
slope.df <- slope.df[!nm%in%c("14_un@all", "144_un@all")]
slope.df$lst <- as.numeric(slope.df$lst)
round(min(slope.df$lst),2)
round(mean(slope.df$lst),2)
round(max(slope.df$lst),2)

