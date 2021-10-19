library(data.table)
library(ggplot2)
library(viridis)

data.raw<-fread ("output/trinity_stats/xn50.out.txt")
data.raw[,`#E`:=as.numeric(gsub("E", "", `Ex`))]

ggplot(data.raw, aes(x=`Ex`, y=`ExN50`))+geom_point(size=2, alpha=0.7, colour="#440154FF")+scale_y_continuous(limits = c(0, 2250))+
  xlab("Top x Percentile of Expression")+theme_bw()
