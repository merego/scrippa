library("reshape")
library("ggplot2")

load("OnlyIBI/RMSLoss.RData")
DMTIBI <- DMT
load("MCSA-IBI/RMSLoss.RData")
DMTMC <- DMT
rm(DMT)

DMTIBI[,6]<-"IBI"
DMTMC[,6]<-"MC"
DMT<-rbind(DMTMC,DMTIBI)
colnames(DMT)[6]<-"Method"
DMTMelt<-melt(DMT)

ggplot(DMTMelt,aes(x=Method,y=value,fill=variable)) + geom_boxplot()
ggplot(DMTMelt,aes(x=variable,y=value,fill=Method)) + geom_boxplot() +
theme(panel.background = element_rect(fill = 'white'), 
      axis.ticks.x = element_line(size=1.2, colour="black"),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(angle=0, vjust=0.5, size=15, face="bold", colour="black"),        
      axis.ticks.y = element_line(size=1.2, colour="black"),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(angle=0, vjust=0.5, size=15, face="bold", colour="black"),        
      plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
      legend.title  = element_blank(),
      legend.text = element_text(lineheight=3, face="bold", color="black", size=15))
