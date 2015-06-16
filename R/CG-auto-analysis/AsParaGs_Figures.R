library("ggplot2")
library("scales")
thm2<- theme(panel.background = element_rect(fill = 'white'),
             axis.ticks.x = element_line(size=0.8, colour="black"),
             axis.title.x =  element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.text.x  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.ticks.y = element_line(size=0.8, colour="black"),
             axis.text.y  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.title.y = element_text(angle=90, vjust=0.5, size=12, face="bold", colour="black"),
             plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
             legend.title  = element_blank(),
             legend.text = element_text(lineheight=3, face="bold", color="black", size=12),
             strip.text = element_text(lineheight=3, face="bold", color="black", size=12)) # Text for facets header

c5 <- c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")
c10 <- c("#a6cee3",  "#1f78b4",   "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


potEne <- function(k1,k2,theta0r,r0,theta,r) {
 V <- 1/2 * k1 * (r - r0)^2 + 1/2 * k2 * (theta - theta0r)^2
 return(V)
}

theta2r <- function(theta) {
 l <- 3.8
 r <- 2.0 * l * sin(theta/2)
 return(r)
}

deg2rad <- function(deg) {
 rad <- deg/180.0*pi
 return(rad)
}

rad2deg <- function(rad) {
 deg <- rad*180.0/pi
 return(deg)
}

alpha <- function(theta0) {
 alpha <- 3.8 * (cos(theta0/2))
 return(alpha)
}

theta0d <- 88.0 # deg
theta0r <- deg2rad(theta0d) # rad
r0 <- theta2r(theta0r) # A

ReadData <- function(TestIndex) {
 filename <- sprintf("report_%02d.dat",TestIndex)
 tmp<-read.table(filename,header=TRUE)
 filename <- sprintf("AvgLoss%02d.dat",TestIndex)
 los<-read.table(filename)[,1]
 REP <- cbind(tmp,los)
 rm(tmp)
 rm(los)
 return(REP)
}

PlotCorrelations <- function(Surface,TestIndex) {
  #colfunc <- colorRampPalette(c10)  
  qn = quantile(Surface$AvgLoss,c(0.01,0.2,0.5,0.7,0.99))
  qn01<-rescale(qn)
  plt <- ggplot(Surface) + geom_point(aes(kr,ktheta,colour=AvgLoss),size=0.5) +
    geom_abline(intercept = 115, slope=-(alpha(theta0r)^2), size=1.2, color="black") +
    scale_colour_gradientn(colours= topo.colors(20), values=c(qn01)) +
    thm2
  if (TestIndex < 7 || TestIndex == 9) 
    plt <- plt + xlab("kr [kcal mol-1 A-2]") + ylab("ktheta [kcal mol-1 rad-2]")
  else
    plt <- plt + xlab("epsi [kcal mol-1 A-2]") + ylab("ktheta [kcal mol-1 rad-2]")
  
  filename <- sprintf("Correlations_%02d.png",TestIndex)
  png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
  print(plt)
  dev.off() 
  return(plt)
}


TestIndex <- 7
REP <- ReadData(TestIndex)
Surface <- data.frame(REP$Param1,REP$Param2,REP$Param1.1,REP$Param2.1,REP$los)
colnames(Surface)<-c("kr","r0","ktheta","theta0","AvgLoss")
PlotCorrelations(Surface,TestIndex)




