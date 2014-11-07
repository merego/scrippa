library("reshape")
library("ggplot2")

load("OnlyIBI/RMSLoss.RData")
DMTIBI <- DMT
load("MCSA-IBI/RMSLoss.RData")
DMTMC <- DMT
rm(DMT)


DMTIBI[,(ncol(DMTIBI)+1)]<-"IBI"
DMTMC[,(ncol(DMTMC)+1)]<-"MCSA-IBI"
DMT<-rbind(DMTMC,DMTIBI)
colnames(DMT)[(ncol(DMTIBI))]<-"Method"
DMTMelt<-melt(DMT)

if (ncol(DMTIBI)==6) {
  ticsLabels <- c(expression(r["i,i+1"]),expression(r["i,i+2"]),expression(r["i,i+3"]),expression(theta),expression(phi)) 
} else {
  ticsLabels <- c(expression(r["i,i+1"]),expression(r["i,i+2"]),expression(theta),expression(phi))
}

#ggplot(DMTMelt,aes(x=Method,y=value,fill=variable)) + geom_boxplot()
cbbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# To use for fills, add

BoxPlot <- ggplot(DMTMelt,aes(x=variable,y=value,fill=Method)) + geom_boxplot() +
  scale_fill_manual(values=cbbPalette) +
  scale_x_discrete(labels=ticsLabels) +
theme(panel.background = element_rect(fill = 'white'), 
      axis.ticks.x = element_line(size=1.2, colour="black"),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(angle=0, vjust=0.5, size=25, face="bold", colour="black"),        
      axis.ticks.y = element_line(size=1.2, colour="black"),
      axis.title.y = element_blank(),
      axis.text.y  = element_text(angle=0, vjust=0.5, size=25, face="bold", colour="black"),        
      plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
      legend.title  = element_blank(),
      legend.text = element_text(lineheight=3, face="bold", color="black", size=20),
      legend.position = c(0.85,0.9))
print(BoxPlot)
ggsave(BoxPlot,file="RMSLossBoxPlot.eps")
