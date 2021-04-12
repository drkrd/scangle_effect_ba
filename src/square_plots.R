
##########################################################################
square_scplot <- function(data, xval, yval)
{
  ggplot(data = data)+
    aes_string(x = xval, y = yval)+
    aes(group=id_placette)+
    geom_point()+
    geom_abline()+
    labs(x = xval,
         y = yval)+
    coord_fixed(xlim = c(min(min(data[[xval]]),min(data[[yval]])),max(max(data[[xval]]), max(data[[yval]]))),
                ylim = c(min(min(data[[xval]]),min(data[[yval]])),max(max(data[[xval]]), max(data[[yval]]))))+
    geom_text(aes(label=meanang), size=3, check_overlap = TRUE)+
    theme_bw()
}

square_scplot(plotmets, "pflidr", "pfvox")
###########################################################################


###########################################################################
ggplot(data=mdl_all)+
  aes(x=type, y=value, group=variable)+
  geom_line()+
  labs(x = "metrics combination" )+
  facet_wrap(vars(variable), scales = "free")+
  theme_bw()





cos(th)=15/hyp(da)




