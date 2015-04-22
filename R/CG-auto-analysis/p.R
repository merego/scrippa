
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
  #AvgLoss<-rowMeans(AllLossFrame[-1,c(2,3)])
  AvgLoss<-rowMeans(AllLossFrame[-1,c(2,5)])
  ipot<-1
  last <- length(DistribALLSpl[[ipot]]$dists)
  step <- 1
  DistFrame <- data.frame(matrix(0,nrow=last,ncol=4))
  for (i in seq(1,last,by=step) ) {    
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    ymax<-which.max(distrib$y)
    x<-distrib$x
    y<-distrib$y
    estm<-x[ymax]    
    #y.fit <- nls(y ~ 1/(s*sqrt(2*pi)) * exp (- (x - m)^2 / (2*s^2) ), start = list(s = 0.5, m=estm))
    y.fit <- nls(y ~ 1/(s*sqrt(2*pi)) * exp (- (x - m)^2 / (2*s^2) ), start = list(s = 0.2, m=estm))
    s<-summary(y.fit)
    EqMean<-s$parameters[2]
    EqStDev<-s$parameters[1]
    Init <- MatrixAllParams[i,2]
    DistFrame[i,1] <- Init
    DistFrame[i,3] <- EqMean
    DistFrame[i,4] <- EqStDev    
  }
  #ipot<-2
  ipot<-4
  for (i in seq(1,last,by=step) ) {    
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    ymax<-which.max(distrib$y)
    x<-180.0 - distrib$x
    y<-distrib$y
    estm<-x[ymax]    
    #y.fit <- nls(y ~ 1/(s*sqrt(2*pi)) * exp (- (x - m)^2 / (2*s^2) ), start = list(s = 5.0, m=estm))
    y.fit <- nls(y ~ 1/(s*sqrt(2*pi)) * exp (- (x - m)^2 / (2*s^2) ), start = list(s = 2.0, m=estm))
    s<-summary(y.fit)
    EqMean<-s$parameters[2]
    EqStDev<-s$parameters[1]
    #Init <-  180.0-MatrixAllParams[i,4] * 180 / pi
    Init <-  180.0-MatrixAllParams[i,11] * 180 / pi
    DistFrame[i,2] <- Init
    DistFrame[i,5] <- EqMean
    DistFrame[i,6] <- EqStDev
  }
  colnames(DistFrame)<-c("InitR13","InitTheta","EqR13Mean","EqR13StDev","EqThetaMean","EqThetaStDev")
  
  AvgLoss<-rowMeans(AllLossFrame[-1,c(2:3)])
  DistFrame$Loss<-AvgLoss
  DistFrame$k13<-MatrixAllParams[-1,1]
  #DistFrame$ktheta<-MatrixAllParams[-1,3]
 DistFrame$ktheta<-MatrixAllParams[-1,10]

ggplot(data=DistFrame) + geom_point(aes(EqThetaMean,EqR13Mean,color=Loss))
kEqR13<-(kb*T)/(sqrt(2*pi) * DistFrame$EqR13StDev^3)
kEqTheta<-(kb*T)/(sqrt(2*pi) * (DistFrame$EqThetaStDev)^3)
DistFrame$kEqTheta<-kEqTheta
DistFrame$kEqR13<-kEqR13
ggplot(data=DistFrame) + geom_point(aes(EqThetaStDev,EqR13StDev,color=Loss))
ggplot(data=DistFrame) + geom_point(aes(k13,kEqR13,color=Loss))
ggplot(data=DistFrame) + geom_point(aes(ktheta,kEqTheta,color=Loss))
ggplot(data=DistFrame) + geom_point(aes(ktheta,k13,color=Loss))
ggplot(data=DistFrame) + geom_point(aes(InitR13,EqR13Mean,color=Loss))
ggplot(data=DistFrame) + geom_point(aes(ktheta,EqThetaMean,color=Loss))

#refAngle <- 92 # degree
#refAngle <- 92/180*pi
#refDistance <- 5.33 # A
#alpha <- (sqrt(2) * 3.8 * sin(refAngle)/(2*sqrt(1-cos(refAngle))))
#minVal <- 20
#Tline<-data.frame(xline=x, yline=2*minVal - y * alpha^2)
kriged <- kriging(DistFrame$k13,DistFrame$ktheta,DistFrame$Loss)
cbbPalette <- c("#D55E00", "#56B4E9", "#009E73", "#F0E442")
ggplot() + geom_raster(data=kriged$map,aes(x=y,y=x,fill=pred),hjust = 0, vjust = 0) +
  scale_fill_gradientn(colours=cbbPalette, limits=c(0,1.0), guide = guide_legend(keyheight = 3)) +
  stat_contour(data=kriged$map,aes(x=y,y=x,z=pred))
#  geom_line(dat=Tline,aes(x=xline,y=yline),  colour="white", size=3.0, linetype="dotted")
