library(reshape2)
df_all_long <- melt(df_all[(df_all$r2<0.6),], id.vars = "r2", measure.vars = c(7:29))

df_all_long$r2 <- as.factor(df_all_long$r2)

df_all_long <- df_all_long[order(df_all_long[,c(1,2)]),]

df_all_long2 <- df_all_long[c(1:1000),]

ggplot(data=df_all_long)+
  aes(x = r2, y= value)+
  geom_boxplot()+
  xlab("Iteration")+
  ylab("Mean scan angle of flightline")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme_classic()


###histogram plots for R²
ggplot(data=df_all)+
  aes(x=r2)+
  geom_histogram()+
  theme_classic2()


gghistogram(
  data = df_all,
  x="r2",
  xlab = "R²",
  title="Histogram of adjusted R² values",
  palette = "jco",
  ggtheme = theme_pubr()
)

