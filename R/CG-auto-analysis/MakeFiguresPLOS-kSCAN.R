refAngle <- 92 # degree
refAngle <- 92/180*pi
refDistance <- 5.33 # A

library("akima")
library("kriging")
library("ggplot2")
xyz<-read.table("k1k2Loss.dat")
bin<-40
# convert to equivalent k-harmonic
xyz[,1] <- 2 * xyz[,1] * 1.3^2 
xyz[,2] <- xyz[,2] * sin(1.6057)^2
kriged <- kriging(xyz[,1],xyz[,2],xyz[,3])
s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = bin), yo = seq(min(xyz[,2]), max(xyz[,2]), length = bin), linear = TRUE,duplicate="mean")

AvgLoss<-read.table("Averagek2Loss.dat")
colnames(AvgLoss)<-c("Avgk1", "Avgk2", "Stdk2", "AvgLoss", "StdLoss")

miy <- min(AvgLoss[,2]-AvgLoss[,3]-1)
may <- max(AvgLoss[,2]+AvgLoss[,3]+1)

cbbPalette <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")

x<-xyz[,1]
y<-xyz[,2]
z<-xyz[,3]
colnames(kriged$map)[3]<-"z"
# Analytical dependence
alpha <- (sqrt(2) * 3.8 * sin(refAngle)/(2*sqrt(1-cos(refAngle))))
minVal <- 20
Tline<-data.frame(xline=s$x, yline=2*minVal - s$y * alpha^2)
gpl1 <- ggplot() + 
        geom_raster(data=kriged$map,aes(x=x,y=y,fill=z),hjust = 0, vjust = 0) +
        scale_fill_gradientn(colours=cbbPalette, limits=c(0,0.2), guide = guide_legend(keyheight = 3)) +        
        #scale_fill_gradient(low = "#F70000", high = "#F7DE00",limits=c(0,0.2), guide = guide_legend(keyheight = 3)) +        
        stat_contour(data=kriged$map,aes(x=x,y=y,z=z)) +
        geom_line(dat=Tline,aes(x=xline,y=yline),  colour="white", size=3.0, linetype="dotted") +
        scale_x_continuous(expand=c(0.01,0), limits=c(min(kriged$map$x),45)) + 
        scale_y_continuous(expand=c(0.01,0.01), limits=range(kriged$map$y)) +
        xlab(expression ( k[r["i,i+1"]] )) + 
        ylab(expression(k[theta])) +        
        theme(
        panel.background = element_rect(fill = 'white'),
        axis.title.x = element_text(face="bold", colour="black", size=30),
        axis.ticks.x = element_line(size=1.2, colour="black"),        
        axis.text.x  = element_text(angle=0, vjust=0.5, size=25, face="bold", colour="black"),
        axis.title.y = element_text(face="bold", colour="black", size=30, angle=0),
        axis.ticks.y = element_line(size=1.2, colour="black"),
        axis.text.y  = element_text(angle=0, vjust=0.5, size=25, face="bold", colour="black"),        
        plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
        legend.title  = element_blank(),
        legend.text = element_text(lineheight=3, face="bold", color="black", size=25))




filename <- "ps/Kr13kThetaLoss_Harmonic_andAnalyticalLine.pdf"
#jpeg(filename,width=2000,height=2000,quality=100)
pdf(filename,width=8.0,height=6.0)
print(gpl1)
dev.off()


