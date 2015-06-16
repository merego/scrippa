library("ggplot2")
library("gridExtra")

#StoredData <- "/home/pmereghetti/data/projects/2014/CGautoTest/FigForPaper/3BEADS-MCMC500-kthetaLOW.RData
#StoredData <- "/home/pmereghetti/data/projects/2014/CGautoTest/FigForPaper/3BEADS-MCMC500-kthetaLOW.RData"
#StoredData <- "/home/pmereghetti/data/projects/2014/CGautoTest/FigForPaper/3BEADS-MCMC2000.RData"
#StoredData <- "/home/pmereghetti/over_ssh/lpgm-pc/2014/CGautoTest/FigForPaper/3BEADS-MCMC10000.RData"
StoredData <- "/home/pmereghetti/data/projects/2014/CGautoTest/FigForPaper/3BEADS-MCMC10000.RData"
#StoredData <- "/home/pmereghetti/data/projects/2014/CGautoTest/FigForPaper/3BEADS-MCMC1000-JSok.RData"
load(StoredData)
thm <- theme(panel.background = element_rect(fill = 'white'),
             panel.border = element_rect(colour = "black", fill=NA, size=1.5),
             axis.ticks.x = element_line(size=1.2, colour="black"),
             axis.title.x = element_text(angle=0, vjust=1., size=15, face="bold", colour="black"),
             axis.text.x  = element_text(angle=0, vjust=1., size=20, face="bold", colour="black"),
             axis.ticks.y = element_line(size=1.2, colour="black"),
             axis.text.y  = element_text(angle=0, vjust=1., size=20, face="bold", colour="black"),
             axis.title.y = element_text(angle=90, vjust=1., size=15, face="bold", colour="black"),
             plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
             panel.grid.major.y = element_line(colour="black",size=0.2),
             panel.grid.major.x = element_line(colour="black",size=0.2))

cbbPalette <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")
# Palette many classes
c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")
kb<-1.098e-3
T<-300
DistribALLSpl <- DATA$DistribALLSpl
AllLossFrame <- DATA$AllLossFrame
potlist <- DATA$potlist
GiuliaDistribs <- DATA$GiuliaDistribs
BestMCDistribs <- DATA$BestMCDistribs
refDistribs <- DATA$refDistribs
ipotindex <- DATA$ipotindex
MatrixAllParams <- DATA$MatrixAllParams
ipotOpti <- DATA$ipotOpti
ipotNparams <- DATA$ipotNparams
  
x<-refDistribs[[1]]$x
y<-refDistribs[[1]]$y
ymax<-which.max(y)
estm<-x[ymax]    
y[y<0.0001]<-0.0001
x.fit <- x[x<7 & x>4]
y.fit <- y[x<7 & x>4]      
BI <- - kb * T * log(y.fit)
BI.fit <- nls(BI ~ 1/2 * k * (x.fit - m)^2 + c, start = list(k = 0.5, m=estm, c=0.01))
kr.fit.ref <- summary(BI.fit)$parameters[1]
r.fit.ref <- summary(BI.fit)$parameters[2]

x<-refDistribs[[2]]$x/180*pi
y<-refDistribs[[2]]$y
estm<-x[ymax]  
ymax<-which.max(y)
y[y<0.0001]<-0.0001
x.fit <- x[x<1.7453 & x>1.3962]
y.fit <- y[x<1.7453 & x>1.3962]            
BI <- - kb * T * log(y.fit)
BI.fit <- nls(BI ~ 1/2 * k * (x.fit - m)^2 + c, start = list(k = 10, m=estm, c=0.01))
ktheta.fit.ref <- summary(BI.fit)$parameters[1]
theta.fit.ref <- summary(BI.fit)$parameters[2]

  ipot<-1
  last <- length(DistribALLSpl[[ipot]]$dists)
  step <- 1
  DistFrame <- data.frame(matrix(0,nrow=last,ncol=8))
  for (i in seq(1,last,by=step) ) {    
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    ymax<-which.max(distrib$y)
    x<-distrib$x
    y<-distrib$y 
    estm<-x[ymax]    
    y[y<0.0001]<-0.0001
    x.fit <- x[x<8 & x>4]
    y.fit <- y[x<8 & x>4]            
    BI <- - kb * T * log(y.fit)
    BI.fit <- nls(BI ~ 1/2 * k * (x.fit - m)^2 + c, start = list(k = 0.5, m=estm, c=0.01))
    #y.fit <- nls(y ~ 1/(s*sqrt(2*pi)) * exp (- (x - m)^2 / (2*s^2) ), start = list(s = 0.5, m=estm))
    ##y.fit <- nls(y ~ 1/(s*sqrt(2*pi)) * exp (- (x - m)^2 / (2*s^2) ), start = list(s = 0.2, m=estm))
    #s<-summary(y.fit)
    k.fit <- summary(BI.fit)$parameters[1]
    eq.fit <- summary(BI.fit)$parameters[2]
    Init <- MatrixAllParams[i,2]
    DistFrame[i,1] <- Init
    DistFrame[i,3] <- eq.fit
    DistFrame[i,4] <- k.fit
    #DistFrame[i,7] <- (5.4450-EqMean)^4 + (0.1795-EqStDev)^2
    dx<-abs(diff(x)[1])
    DistFrame[i,7] <- sum((refDistribs[[1]]$y - y)^2) * dx
  }
  ipot<-2
  #ipot<-4
  for (i in seq(1,last,by=step) ) {    
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    ymax<-which.max(distrib$y)
    x<-distrib$x/180*pi
    y<-distrib$y 
    estm<-x[ymax]    
    y[y<0.0001]<-0.0001
    x.fit <- x[x<1.7453 & x>1.3962]
    y.fit <- y[x<1.7453 & x>1.3962]           
    BI <- - kb * T * log(y.fit)
    BI.fit <- nls(BI ~ 1/2 * k * (x.fit - m)^2 + c, start = list(k = 0.5, m=estm, c=0.01))
    #y.fit <- nls(y ~ 1/(s*sqrt(2*pi)) * exp (- (x - m)^2 / (2*s^2) ), start = list(s = 1.0, m=estm)) # 5
    ##y.fit <- nls(y ~ 1/(s*sqrt(2*pi)) * exp (- (x - m)^2 / (2*s^2) ), start = list(s = 2.0, m=estm))
    #s<-summary(y.fit)
    k.fit <- summary(BI.fit)$parameters[1]
    eq.fit <- summary(BI.fit)$parameters[2]
    Init <-  MatrixAllParams[i,4] * 180 / pi
    #Init <-  180.0-MatrixAllParams[i,11] * 180 / pi
    DistFrame[i,2] <- Init
    DistFrame[i,5] <- eq.fit
    DistFrame[i,6] <- k.fit
    #DistFrame[i,8] <- (91.055-EqMean)^4 + (3.704-EqStDev)^2
    dx<-abs(diff(x)[1])
    DistFrame[i,8] <- sum((refDistribs[[ipot]]$y - y)^2) * dx
  }
  colnames(DistFrame)<-c("InitR13","InitTheta","r13.fit","k13.fit","theta.fit","ktheta.fit","cLossR","cLossA")
  
  AvgLoss<-vector()#rowMeans(AllLossFrame[,c(2:3)]) * range01(1 - abs(AllLossFrame[,2]-AllLossFrame[,3])^2)
  for (i in 1:length(AllLossFrame[,2])) {
  #  if (AllLossFrame[i,2] < 0.25 & AllLossFrame[i,3] < 0.25) {
       AvgLoss[i] <- (AllLossFrame[i,2] + AllLossFrame[i,3] ) * 0.5
       
  #  } else {
  #    AvgLoss[i] <- 1      
  # }
  }
  AvgLoss[1] <- mean(AvgLoss[-1])
  DistFrame$Loss<-AvgLoss
  DistFrame$k13<-MatrixAllParams[,1]
  DistFrame$ktheta<-MatrixAllParams[,3]
  DistFrame$newLoss<-DistFrame$ktheta.fit-ktheta.fit.ref

# kTheta vs kr13 Initial
alpha <- 3.8^2 * (cos(1.57/2))^2 # Vale == to ( 3.8 * sin(1.57) / sqrt(2- 2*cos(1.57)) ) ^2 (Paolo)
ggplot(data=DistFrame) +
  geom_point(aes(ktheta,k13,color=newLoss,size=(1/newLoss)^2)) + 
  geom_abline(intercept = 84, slope=-alpha, size=1.1, color="black") +
  geom_abline(intercept = 84, slope=-alpha, size=1.0, color=c25[3]) + 
  thm + 
  scale_color_gradientn(name="JS",colours=cbbPalette) + ggtitle("Jensen-Shannon Loss") + scale_size_continuous(name="1/JS")
  #DistFrame$ktheta<-MatrixAllParams[-1,10]

# kEqR13<-(kb*T)/(sqrt(2*pi) * DistFrame$EqR13StDev^3)
# kEqTheta<-(kb*T)/(sqrt(2*pi) * (DistFrame$EqThetaStDev/10)^3)
# DistFrame$kEqTheta<-kEqTheta
# DistFrame$kEqR13<-kEqR13
# 
#   DistFrame[1,]$Loss <- 0


# # Theta vs r13 equilib.
# plt <- ggplot(data=DistFrame) + 
#   geom_point(aes(EqThetaMean,EqR13Mean,color=Loss)) + 
#   thm
# filename<-"3Beads10000/EqTheta-vs-Eqr13.jpg"
# jpeg(file=filename, width=800)
# print(plt)
# dev.off()
# 
# # Theta vs r13 initial and Theta vs r13 equilib.
# plt <- ggplot(data=DistFrame) + 
#   geom_point(aes(EqThetaMean,EqR13Mean,color=Loss)) + 
#   geom_point(aes(InitTheta,InitR13,color=Loss)) +
#   thm
# filename<-"3Beads10000/EqTheta-vs-Eqr13-AND-Theta-vs-r13.jpg"
# jpeg(file=filename, width=800)
# print(plt)
# dev.off()





# kTheta vs kr13 Equil.
#p1<-ggplot(data=DistFrame) +
#  geom_point(aes(kEqTheta,kEqR13,color=Loss,size=1/Loss)) +
#  thm

DistFrame$step <- seq(1,(nrow(DistFrame)))
DistFrame$LossR13 <- AllLossFrame[,2]
DistFrame[1,]$LossR13<-0
DistFrame$LossTheta <- AllLossFrame[,3]
DistFrame[1,]$LossTheta <-0
# MSD normalized
DistFrame$cLossR <- DistFrame$cLossR/max(DistFrame$cLossR)
DistFrame$cLossA <- DistFrame$cLossA/max(DistFrame$cLossA)
DistFrame$AvgcLoss <- 0.5 * (DistFrame$cLossR + DistFrame$cLossA)

p1 <- ggplot(DistFrame) + 
  geom_line(aes(step,LossR13,col="JS(R13)"),size=2) + 
  geom_line(aes(step,LossTheta,col="JS(Theta)"),size=2) + 
  geom_line(aes(step,Loss,col="AvgJS"),size=1) + thm + ylab("JS loss") + 
  scale_color_manual(name="JS", values=c("JS(R13)"=c25[1], "JS(Theta)"=c25[2],"AvgJS"="black"))

p2 <- ggplot(DistFrame) +
  geom_point(aes(LossR13,LossTheta),col="black",size=1) + geom_abline(1,size=2) + thm + xlab("JS(r13)") + ylab("JS(theta)")

p3 <- ggplot(DistFrame) + 
  geom_line(aes(step,cLossR,col="MSD(R13)"),size=2) + 
  geom_line(aes(step,cLossA,col="MSD(Theta)"),size=2) + 
  geom_line(aes(step,AvgcLoss,col="AvgMSD"),size=1) + thm + ylab("MSD loss") + 
  scale_color_manual(name="MSD",values=c("MSD(R13)"=c25[1], "MSD(Theta)"=c25[2],"AvgMSD"="black"))
  

p4 <- ggplot(DistFrame) +
  geom_point(aes(cLossR,cLossA),col="black",size=1) + geom_abline(1,size=2) + thm + xlab("MSD(r13)") + ylab("MSD(theta)")

filename<-"3Beads10000/JSLoss-and-MSDLoss.ps"
postscript(file=filename, width=16)
grid.arrange(p1,p2,p3,p4)
dev.off()

# kTheta vs kr13 Initial
alpha <- 3.8^2 * (cos(1.57/2))^2 # Vale == to ( 3.8 * sin(1.57) / sqrt(2- 2*cos(1.57)) ) ^2 (Paolo)
# kthetaeff = ktheta + alpha * k13
# ktheta = ktehtaeff - alpha * k13
# k13 = kthetaeff / alpha - ktheta / alpha
#DistFrame$krfit<- 20.0 + alpha * DistFrame$k13
p1<-ggplot(data=DistFrame) +
  geom_point(aes(ktheta,k13,color=Loss,size=(1/Loss)^2)) + 
  geom_abline(intercept = 70/alpha, slope=-1/alpha, size=2.5, color="black") +
  geom_abline(intercept = 70/alpha, slope=-1/alpha, size=2, color=c25[3]) + 
  thm + 
  scale_color_gradientn(name="JS",colours=cbbPalette) + ggtitle("Jensen-Shannon Loss") + scale_size_continuous(name="1/JS")

p2<-ggplot(data=DistFrame) +
  geom_point(aes(ktheta,k13,color=AvgcLoss,size=(1/AvgcLoss)^2)) + 
  geom_abline(intercept = 70/alpha, slope=-1/alpha, size=2.5, color="black") +
  geom_abline(intercept = 70/alpha, slope=-1/alpha, size=2, color=c25[3]) +
  thm +
  scale_color_gradientn(name="MSD",colours=cbbPalette) + ggtitle("MSD Loss") +  scale_size_continuous(name="1/MSD")

filename<-"3Beads10000/kTheta-vs-kR13-JSLoss-and-MSDLoss.ps"
postscript(file=filename, width=8)
grid.arrange(p1,p2)
dev.off()

# Distributions for frame n
framex <- DistFrame[(DistFrame$Loss<0.23 & DistFrame$InitTheta > 96 & DistFrame$InitR13 < 5.3),]$step[1]
#framex <- DistFrame[(DistFrame$Loss<0.23 & DistFrame$InitTheta < 85 & DistFrame$InitR13 < 5.3),]$step[1]


# Step bassa loss ma fuori dal linea geometrica
p1 <- ggplot(DistFrame) + 
  geom_point(aes(InitTheta,InitR13,color=Loss,size=(1/Loss))) +   geom_point(aes(InitTheta[framex],InitR13[framex]),size=5) + 
  geom_point(aes(InitTheta[684],InitR13[684]),size=5) + 
  geom_line(aes(InitTheta,(2*3.8*sin(InitTheta/180*pi/2))),size=2) + 
  scale_color_gradientn(name="JS",colours=cbbPalette) + 
  ggtitle("Jensen-Shannon Loss") + 
  scale_size_continuous(name="1/JS") + 
  ylim(5.2,5.7) +
  thm

p2 <- ggplot(DistFrame) + 
  geom_point(aes(InitTheta,InitR13,color=LossR13,size=(1/LossR13)^2)) +  
  geom_line(aes(InitTheta,(2*3.8*sin(InitTheta/180*pi/2))),size=2) + 
  scale_color_gradientn(name="JS",colours=cbbPalette) + 
  ggtitle("R13 Jensen-Shannon Loss") + 
  scale_size_continuous(name="1/JS") + 
  ylim(5.2,5.7) +
  thm

p3 <- ggplot(DistFrame) + 
  geom_point(aes(InitTheta,InitR13,color=LossTheta,size=(1/LossTheta)^2)) +  
  geom_line(aes(InitTheta,(2*3.8*sin(InitTheta/180*pi/2))),size=2) + 
  scale_color_gradientn(name="JS",colours=cbbPalette) + 
  ggtitle("Theta Jensen-Shannon Loss") + 
  scale_size_continuous(name="1/JS") + 
  ylim(5.2,5.7) +
  thm

filename<-"3Beads10000/Theta-R13-JSLoss.ps"
postscript(file=filename, paper="a4")
grid.arrange(p1,p2,p3)
dev.off()



DistsR <- data.frame(x=refDistribs[[1]]$x,y=refDistribs[[1]]$y,type=rep("refr13",402))
DistsR <- rbind(DistsR,data.frame(x=DistribALLSpl[[1]]$dists[[framex]]$x,y=DistribALLSpl[[1]]$dists[[framex]]$y,type=rep("R13",402)))
#P <- read.table("../NewData/3Beads/MCMC/OUTPUT-Giulia/r0/Param1.dat")
#dx=diff(P[,1])[1]
#y<-P[,2] / (sum(P[,2]) * dx)
#DistsR <- rbind(DistsR,data.frame(x=P[,1],y=y,type=rep("R13",402)))

DistsA <- data.frame(x=refDistribs[[2]]$x,y=refDistribs[[2]]$y,type=rep("refTheta",402))
DistsA <- rbind(DistsA,data.frame(x=DistribALLSpl[[2]]$dists[[framex]]$x,y=DistribALLSpl[[2]]$dists[[framex]]$y,type=rep("Theta",402)))
#P <- read.table("../NewData/3Beads/MCMC/OUTPUT-Giulia/r0/Param2.dat")
#P[,1]<-P[,1]*180/pi
#dx=diff(P[,1])[1]
#y<-P[,2] / (sum(P[,2]) * dx)
#DistsA <- rbind(DistsA,data.frame(x=P[,1],y=y,type=rep("R13",402)))
title1<-paste("R13 frame",framex,sep="")
title2<-paste("Theta frame",framex,sep="")
p1 <- ggplot(DistsR) + geom_line(aes(x,y,group=type,colour=type),size=2) + thm + xlim(4,7) + xlab("R13") + ggtitle(title1)# + annotate("text",x=5.3,y=1.9,label="r0=5.33")
p2 <- ggplot(DistsA) + geom_line(aes(x,y,group=type,colour=type),size=2) + thm  + xlim(70,120) + xlab("Theta") + ggtitle(title2)# + annotate("text",x=88.9,y=0.1,label="theta0=88.6")

#filename<-"3Beads10000/DistsFrame337.ps"
filename<-"3Beads10000/DistsFrame684.ps"
postscript(file=filename, paper="a4")
grid.arrange(p1,p2)
dev.off()


MCMC10000<-data.frame(DistFrame$step,DistFrame$InitTheta,DistFrame$ktheta,DistFrame$InitR13,DistFrame$k13,DistFrame$LossTheta,DistFrame$LossR13,DistFrame$Loss,DistFrame$cLossA,DistFrame$cLossR,DistFrame$AvgcLoss)
colnames(MCMC10000) <- c("step","theta","ktheta","r13","kr13","JSLossTheta","JSLossR13","JSLossAvg","MSDLossTheta","MSDLossR13","MSDLossAvg")

filename <- "3Beads10000/MCMC10000.dat"
write.table(file=filename,MCMC10000,row.names=FALSE)



# ggplot(data=DistFrame) + geom_point(aes(k13,kEqR13,color=Loss))
# ggplot(data=DistFrame) + geom_point(aes(ktheta,kEqTheta,color=Loss))
# ggplot(data=DistFrame) + geom_point(aes(ktheta,k13,color=Loss))
# ggplot(data=DistFrame) + geom_point(aes(InitR13,EqR13Mean,color=Loss))
# ggplot(data=DistFrame) + geom_point(aes(ktheta,EqThetaMean,color=Loss))

# # k13 initial vs Theta equilb (k alte Theta eq. reach the value set in xml)
# plt <- ggplot(data=DistFrame) + geom_point(aes(k13,EqThetaMean,color=Loss)) +thm
# filename<-"3Beads10000/kr13-vs-EqTheta.jpg"
# jpeg(file=filename, width=800)
# print(plt)
# dev.off()
# # k13 initial vs r13 equilb (k alte Theta eq. reach the value set in xml)
# plt <- ggplot(data=DistFrame) + geom_point(aes(k13,EqThetaMean,color=Loss)) + thm
# filename<-"3Beads10000/kr13-vs-Eqr13.jpg"
# jpeg(file=filename, width=800)
# print(plt)
# dev.off()
# # ktheta iniziale vs ktheta finale no corelazione
# p1<-ggplot(data=DistFrame) + geom_point(aes(ktheta,kEqTheta,color=Loss)) + thm
# # ktheta iniziale vs k13 no correlazione
# p2<-ggplot(data=DistFrame) + geom_point(aes(ktheta,kEqR13,color=Loss)) + thm
# # k13 iniziale vs ktheta finale no corelazione si corr
# p3<-ggplot(data=DistFrame) + geom_point(aes(k13,kEqTheta,color=Loss)) + thm
# # k13 iniziale vs k13 finale no corelazione si corr 
# p4<-ggplot(data=DistFrame) + geom_point(aes(k13,kEqR13,color=Loss)) + thm
# filename<-"3Beads10000/kinit-vs-kEq.jpg"
# jpeg(file=filename, width=800)
# grid.arrange(p1,p2,p3,p4)
# dev.off()


# # ktheta/kr13 initial vs ktheta/kr13 equilib.
# plt <- ggplot(data=DistFrame) + 
#   geom_point(aes(kEqTheta/ktheta,kEqR13/k13,color=Loss)) + 
#   xlim(0,1.0) + 
#   ylim(0.75,1.25) +
#   thm
# filename<-"3Beads10000/kThetaOverKr13-vs-kEqThetaOverkEqr13.jpg"
# jpeg(file=filename, width=800)
# print(plt)
# dev.off()


#refAngle <- 92 # degree
#refAngle <- 92/180*pi
#refDistance <- 5.33 # A
#alpha <- (sqrt(2) * 3.8 * sin(refAngle)/(2*sqrt(1-cos(refAngle))))
#minVal <- 20
#Tline<-data.frame(xline=x, yline=2*minVal - y * alpha^2)
# kriged <- kriging(DistFrame$ktheta,DistFrame$k13,DistFrame$Loss)
# #cbbPalette <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")
# plt <- ggplot() + geom_raster(data=kriged$map,aes(x=x,y=y,fill=pred),hjust = 0, vjust = 0) +
#   scale_fill_gradientn(colours=cbbPalette, limits=c(0,1.0), guide = guide_legend(keyheight = 3)) +
#   stat_contour(data=kriged$map,aes(x=x,y=y,z=pred)) +
#   geom_line(data=DistFrame,aes(ktheta,krfit),size=2,color=cbbPalette[4]) + xlim(0,10) + ylim(1,10) + thm + xlab("ktheta") + ylab("kr13")
# filename<-"3Beads10000/kTheta-kR13-MAP.jpg"
# png(file=filename, width=800)
# print(plt)
# dev.off()