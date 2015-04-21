
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
  AvgLoss<-rowMeans(AllLossFrame[-1,c(2:3)])
  ipot<-1
  last <- length(DistribALLSpl[[ipot]]$dists)
  step <- 1
  DistFrame <- data.frame(matrix(0,nrow=last,ncol=4))
  for (i in seq(1,last,by=step) ) {    
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    ymax<-which.max(distrib$y)
    EqR13<-distrib$x[ymax]    
    InitR13 <- MatrixAllParams[i,2]
    DistFrame[i,1] <- InitR13
    DistFrame[i,3] <- EqR13
  }
  ipot<-2
  for (i in seq(1,last,by=step) ) {    
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    ymax<-which.max(distrib$y)
    EqTheta<-180.0-distrib$x[ymax]  
    InitTheta <- 180.0 - MatrixAllParams[i,4] * 180 / pi
    DistFrame[i,2] <- InitTheta
    DistFrame[i,4] <- EqTheta
  }
  colnames(DistFrame)<-c("InitR13","InitTheta","EqR13","EqTheta")
  
  AvgLoss<-rowMeans(AllLossFrame[-1,c(2:3)])
  DistFrame$Loss<-AvgLoss
  DistFrame$k13<-MatrixAllParams[-1,1]
  DistFrame$ktheta<-MatrixAllParams[-1,3]