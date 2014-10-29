AnaAndPlot <- function(ParamAllIterations,ParametersNames,DistNames,TotalNumberOfOptimizedParameters,NumberOfOptimizedPotentials,OptimPotIndex,n_of_runs,optim,ReportTexFID,DistribAllMomenta) {
 

   sf <- paste(Sys.getenv("Rscripts"),"SpecialPlotFunctions.R",sep="")

   source(sf)

   library("ellipse")
   library("akima") # Interpolation 1D/2D
   library("fpc") # Density based clustering (dbscan function)

 
   NumberOfPotentials <- length(ParamAllIterations[[1]])
   # Cast parameters into dataframe for easier operations
   # This contains only parameters for optimized potentials
   ParametersFrame <- as.data.frame(matrix(0, n_of_runs, TotalNumberOfOptimizedParameters))
   ParametersFrameFIT <- as.data.frame(matrix(0, n_of_runs-1, TotalNumberOfOptimizedParameters))
   colnames(ParametersFrame) <- ParametersNames
   colnames(ParametersFrameFIT) <- ParametersNames
   
   # Cast loss function int dataframe
   LossFrame <- as.data.frame(matrix(NA, n_of_runs, NumberOfPotentials))
   colnames(LossFrame) <- DistNames
 
   # Cast Momenta into dataframe
   Momenta <- as.data.frame(matrix(NA, n_of_runs-1, NumberOfPotentials))
   colnames(Momenta) <- DistNames
   
   # Fill the frames
   for (run in 1:n_of_runs) {
     Nstart <- 1
     NstartF <- 1
     for (i in OptimPotIndex) {
       Npara <- length(ParamAllIterations[[run]][[i]]$params) 
       ParametersFrame[run,(Nstart:(Nstart+Npara-1))] <- ParamAllIterations[[run]][[i]]$params 
       if (optim=="MC") {
         if (run<n_of_runs) {
          NparaF <- length(ParamAllIterations[[run]][[i]]$paramsFIT)
          ParametersFrameFIT[run,(NstartF:(NstartF+NparaF-1))] <- ParamAllIterations[[run]][[i]]$paramsFIT 
         }
         NstartF <- NstartF + NparaF
       }
       LossFrame[run,i] <- ParamAllIterations[[run]][[i]]$loss
       Nstart <- Nstart + Npara 
       if (optim=="MC") {
         if (run<n_of_runs) {
          Momenta[run,i] <- DistribAllMomenta[[run]][[i]]$moments
         }
       }
     }
   }

   
if (optim=="MC") {
  cat("\\section{Monte-Carlo simulated annealing}\n",file=ReportTexFID)
  cat("\\input{BestBasedOnAvgLoss_MC}\n",file=ReportTexFID) 
} 
if (optim=="IBI") {
  cat("\\section{Iterative Boltzmann inversion}\n",file=ReportTexFID)
  cat("\\input{BestBasedOnAvgLoss_IBI}\n",file=ReportTexFID)
}
   
   # For each optimized potential plot a pairs of selected parmeters
   #Nstart <- 1
   #for (i in 1:NumberafOptimizedPotentials) {
   #   xyz <- ParametersFrame[,(Nstart:(Nstart+Npara-1))]
   #   LossFrame[run,i] <- ParamAllIterations[[run]][[i]]$loss
   #   Nstart <- Nstart + Npara 
   #}
   
   # Compute Average Loss
   AvgLoss <- rowMeans(LossFrame,na.rm=TRUE)
 
   if (optim=="MC") {
   bin <- 40

   # Paramsurf r13Angle
   if((sum(ParametersNames=="r13_P1")+sum(ParametersNames=="angle_P1"))==2) {
   AvgLossR13Angle <- rowMeans(LossFrame[,c("r13_Loss","angle_Loss")])
   xyz <- cbind(ParametersFrame[2:length(AvgLossR13Angle),c("r13_P1","angle_P1")],AvgLossR13Angle[2:length(AvgLossR13Angle)])
   colnames(xyz)<-c("r13_P1","angle_P1","AvgLossR13Angle")
   s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = bin), yo = seq(min(xyz[,2]), max(xyz[,2]), length = bin), linear = TRUE,duplicate="mean")
   
   # Find largest basin
   zz <- s$z/max(s$z,na.rm=TRUE)
   threshold<-quantile(zz, probs = c(0.1), na.rm = TRUE)
   zzbelow<-zz
   zzbelow[zzbelow<threshold]<-0.0
   zzbelow[zzbelow>threshold]<-1.0
   index<-which(zzbelow==0,arr.ind=TRUE)
   AB<-matrix(0.0,dim(index)[1],2)
   AB[,1]<-s$x[index[,1]]
   AB[,2]<-s$y[index[,2]]
   # Density based clustering
   #clusters<-dbscanCBI(AB,eps=0.5,MinPts=3)
   #LargestPartition <- clusters$clusterlist[[1]]
   #ii<-which(LargestPartition==TRUE)
   #LargestAB <- AB[ii,]

   # Kmeans (automatic selection number of partitions)
   #  Scanning 2 to 5 clusters select best based on Calinski-Harabasz
   clusters<-kmeansruns(AB,krange=2:5,critout=TRUE,runs=10,criterion="ch")
   best<-which.max(clusters$size)
   ii<-which(clusters$cluster==best)
   LargestAB <- AB[ii,]

   filename <- "ps/R13ANGLE_SURF.eps"
   postscript(filename)
   image(s$x, s$y, s$z/max(s$z,na.rm=TRUE), xlab="k r13", ylab="k angle",zlim=c(0.0,0.5))
   contour(s$x, s$y, s$z/max(s$z,na.rm=TRUE), add=TRUE, nlevel=10, lwd=0.3)
   # points(LargestAB)
   lines(ellipse(sqrt(abs(cov(LargestAB))),centre=colMeans(LargestAB)),lwd=3)  

   dev.off()
 cat(
"\\begin{figure}[!hbp]
{\\includegraphics[width=4.0in,angle=-90]{",file=ReportTexFID)
	cat(filename,file=ReportTexFID)
	cat(
"}}
\\caption{k r13 vs k angle vs Loss function}
\\label{fig:r13}
\\end{figure}\n\n",file=ReportTexFID)
   }

   # Paramsurf r14 dihedral
   if((sum(ParametersNames=="r14_P1")+sum(ParametersNames=="dihedral_P1"))==2) {
   AvgLossR14Dih <- rowMeans(LossFrame[,c("r14_Loss","dihedral_Loss")])
   xyz <- cbind(ParametersFrame[2:length(AvgLossR14Dih),c("r14_P1","dihedral_P1")],AvgLossR14Dih[2:length(AvgLossR14Dih)])
   colnames(xyz)<-c("r13_P1","angle_P1","AvgLossR14Dih")
   s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = bin), yo = seq(min(xyz[,2]), max(xyz[,2]), length = bin), linear = TRUE,duplicate="mean")

   # Find largest basin
   zz <- s$z/max(s$z,na.rm=TRUE)
   threshold<-quantile(zz, probs = c(0.1), na.rm = TRUE)
   zzbelow<-zz
   zzbelow[zzbelow<threshold]<-0.0
   zzbelow[zzbelow>threshold]<-1.0
   index<-which(zzbelow==0,arr.ind=TRUE)
   AB<-matrix(0.0,dim(index)[1],2)
   AB[,1]<-s$x[index[,1]]
   AB[,2]<-s$y[index[,2]]
   # Density based clustering
   #clusters<-dbscanCBI(AB,eps=0.5,MinPts=3)
   #LargestPartition <- clusters$clusterlist[[1]]
   #ii<-which(LargestPartition==TRUE)
   #LargestAB <- AB[ii,]

   # Kmeans (automatic selection number of partitions)
   #  Scanning 2 to 5 clusters select best based on Calinski-Harabasz
   clusters<-kmeansruns(AB,krange=2:5,critout=TRUE,runs=10,criterion="ch")
   best<-which.max(clusters$size)
   ii<-which(clusters$cluster==best)
   LargestAB <- AB[ii,]

   
   filename <- "ps/R14Dih_SURF.eps"
   postscript(filename)
   image(s$x, s$y, s$z/max(s$z,na.rm=TRUE), xlab="k r14", ylab="k dihedral",zlim=c(0.0,0.5))
   contour(s$x, s$y, s$z/max(s$z,na.rm=TRUE), add=TRUE, nlevel=10, lwd=0.3)
   # points(LargestAB)
   lines(ellipse(sqrt(abs(cov(LargestAB))),centre=colMeans(LargestAB)),lwd=3)  

   dev.off()
 cat(
"\\begin{figure}[!hbp]
{\\includegraphics[width=4.0in,angle=-90]{",file=ReportTexFID)
	cat(filename,file=ReportTexFID)
	cat(
"}}
\\caption{k r14 vs k dihedral vs Loss function}
\\label{fig:r13}
\\end{figure}\n\n",file=ReportTexFID)

   }
   }

   # Best params based on SqSum and Var
   Vm <- vector()
   SqSum <- vector()
   VarM <- vector()
   VarScaled <- vector()
   for (run in 1:nrow(LossFrame)) {
     SqSum[run] <- sum(LossFrame[run,]^2,na.rm=TRUE)
     if (length(OptimPotIndex)>1) 
       VarM[run] <- sd(LossFrame[run,],na.rm=TRUE)^2
     else
       VarScaled[run] <- 0.0
   }
   SqScaled<-(SqSum-min(SqSum))/(max(SqSum)-min(SqSum))
   if (length(OptimPotIndex)>1) 
     VarScaled<-(VarM-min(VarM))/(max(VarM)-min(VarM))
   Vm <- SqScaled + VarScaled

   

   # Sort 
   #avgSort <- sort(AvgLoss,index.return=TRUE) 
   avgSort <- sort(Vm,index.return=TRUE) 
   SortedParametersFrame <- ParametersFrame[avgSort$ix[1:min(5,length(avgSort$ix))],]
   #s[,10] <- (s[,10]+3.14)*180/3.14
   #s[,8] <- s[,8]*180/3.14
   all_mean <- sapply(SortedParametersFrame,mean)
   all_std <- sapply(SortedParametersFrame,sd)
   SortedParametersFrame <- rbind(SortedParametersFrame,all_mean)
   SortedParametersFrame <- rbind(SortedParametersFrame,all_mean-all_std)
   SortedParametersFrame <- rbind(SortedParametersFrame,all_mean+all_std)
   rownames(SortedParametersFrame)[nrow(SortedParametersFrame)-2]<-c("Mean")
   rownames(SortedParametersFrame)[nrow(SortedParametersFrame)-1]<-c("-StdDev")
   rownames(SortedParametersFrame)[nrow(SortedParametersFrame)]<-c("+StdDev")
   # Save out latex table
   caption <- paste("Best parameters based on average loss function. Params are sorted on average loss function and the top 5 (max) are selected. Based on ",optim," simulations.",sep="")
   tab <- xtable(t(SortedParametersFrame[,]), caption=caption)
   filename <- paste("BestBasedOnAvgLoss_",optim,".tex",sep="")
   print(tab,file=filename,append=F,table.placement = "h", caption.placement="bottom")

   
   Nstart <- 1
   for (i in OptimPotIndex) { 
     Npara <- length(ParamAllIterations[[2]][[i]]$params)
     for (j in 1:Npara) {
       val <- SortedParametersFrame[nrow(SortedParametersFrame)-2,Nstart:Nstart+j-1]
       cat("<pot_parameter type=\"double\">",val,"</pot_parameter>\n")
     }
     for (j in 1:Npara) {
       val <- SortedParametersFrame[nrow(SortedParametersFrame)-1,Nstart:Nstart+j-1]
       cat("<pot_parameter_min type=\"double\">",val,"</pot_parameter_min>\n")
     }
     for (j in 1:Npara) {
       val <- SortedParametersFrame[nrow(SortedParametersFrame),Nstart:Nstart+j-1]
       cat("<pot_parameter_max type=\"double\">",val,"</pot_parameter_max>\n")
     }
     Nstart <- Nstart + Npara
   }
#SKIP#   <pot_parameter type="double">5.5</pot_parameter>
#SKIP#   <pot_parameter type="double">1.1</pot_parameter>
#SKIP#   <pot_parameter_min type="double">1.0</pot_parameter_min>
#SKIP#   <pot_parameter_min type="double">5.2</pot_parameter_min>
#SKIP#   <pot_parameter_min type="double">0.5</pot_parameter_min>
#SKIP#   <pot_parameter_max type="double">20.0</pot_parameter_max>
#SKIP#   <pot_parameter_max type="double">5.6</pot_parameter_max>
#SKIP#   <pot_parameter_max type="double">1.5</pot_parameter_max>
#SKIP#   

   # Plot radar
   
   NumberOfLineToPlot <- 10
   t <- seq(0,2*pi-(2*pi/NumberOfOptimizedPotentials),2*pi/NumberOfOptimizedPotentials)
   rgb.palette <- colorRampPalette(c("red", "green", "blue"),space="rgb")
   colors <- rgb.palette(NumberOfLineToPlot+1)

   filename <- "ps/RadarPlot.eps"
   postscript(filename)
   plot(0,0)
   x <- matrix(0,NumberOfLineToPlot,NumberOfOptimizedPotentials)
   y <- matrix(0,NumberOfLineToPlot,NumberOfOptimizedPotentials)
   for (run in 1:NumberOfLineToPlot) {
    x[run,] <- as.matrix(LossFrame[avgSort$ix[run],OptimPotIndex]*cos(t))
    y[run,] <- as.matrix(LossFrame[avgSort$ix[run],OptimPotIndex]*sin(t))
    polygon(x[run,],y[run,],col=colors[run],density=0)
   }
   legend(-1,0,avgSort$ix[c(1:10)],col=colors,lty=1)
   segment <- max(abs(x),abs(y))
   segment_x <- segment*cos(t)
   segment_y <- segment*sin(t)
   segments(x0=0,y0=0,x1=segment_x,y1=segment_y,col="gray")
   tcircle <- seq(0,2*pi-(2*pi/100),2*pi/100)
   iso <- 0
   for (i in 1:5) {
     iso <- iso + segment/5
     circle <- rep(iso,100)
     tics_x <- circle*cos(tcircle)
     tics_y <- circle*sin(tcircle)
     lines(tics_x,tics_y,lty=2,col="gray")
   }
  #plot.window(xlim=c(min(x),max(x)),ylim=c(min(y),max(y)))
  dev.off()
 cat(
"\\begin{figure}[!hbp]
{\\includegraphics[width=4.0in,angle=-90]{",file=ReportTexFID)
	cat(filename,file=ReportTexFID)
	cat(
"}}
\\caption{Radar Plot Loss function top 10}
\\label{fig:r13}
\\end{figure}\n\n",file=ReportTexFID)



   if (length(OptimPotIndex)>1) { 
   # Pearson correlation among the top 20%
   colsc <- c(rgb(241, 54, 23, maxColorValue=255), 'white', rgb(0, 61, 104, maxColorValue=255))
   colramp <- colorRampPalette(colsc, space='Lab')
   colors <- colramp(100)
   if (optim=="MC") {
     M <- Momenta[,OptimPotIndex]
     PearsonCor <- cor(M)
     colnames(PearsonCor) <- DistNames[OptimPotIndex]
     rownames(PearsonCor) <- DistNames[OptimPotIndex]
     PearsonCor[lower.tri(PearsonCor,diag=T)] = 0.0 # set 0 diagonal and lower tri matrix
     # Add p-values to the lower triangular part of Pearson matrix
     for (i in 1:dim(M)[2]) for (j in i:dim(M)[2]) {
        pvalue <- cor.test(M[,i],M[,j],alternative="two.side")$p.val
        if (pvalue<0.01) {  PearsonCor[j,i] <- 0.01 }
     }
     filename <- "ps/PearsonCorrelation.eps"
     postscript(filename) 
     my.plotcorr(PearsonCor, col=colors[((PearsonCor + 1)/2) * 100], main='Pearson correlations', lower.panel='none', diag='none',mar=0.0 + c(1, 0, 2, 2),outline=TRUE)
     image.plot(legend.only=TRUE,col=colors, zlim=c(-1,1),horizontal=TRUE,legend.mar=6.4,legend.shrink=0.4,legend.width=0.4 )
     dev.off()

# PARAMETERS_CORRELATION     M <- ParametersFrame[avgSort$ix[1:ceiling(length(AvgLoss)*0.1)],]
# PARAMETERS_CORRELATION     MFIT <- ParametersFrameFIT[avgSort$ix[1:ceiling(length(AvgLoss)*0.1)],]
# PARAMETERS_CORRELATION     PearsonCor <- cor(M,MFIT)
# PARAMETERS_CORRELATION     colnames(PearsonCor) <- ParametersNames
# PARAMETERS_CORRELATION     PearsonCor[lower.tri(PearsonCor,diag=T)] = 0.0 # set 0 diagonal and lower tri matrix
# PARAMETERS_CORRELATION     
# PARAMETERS_CORRELATION     # Add p-values to the lower triangular part of Pearson matrix
# PARAMETERS_CORRELATION     for (i in 1:dim(M)[2]) for (j in i:dim(M)[2]) {
# PARAMETERS_CORRELATION        pvalue <- cor.test(M[,i],MFIT[,j],alternative="two.side")$p.val
# PARAMETERS_CORRELATION        if (pvalue<0.05 && is.finite(pvalue)) {  PearsonCor[j,i] <- pvalue }
# PARAMETERS_CORRELATION     }
# PARAMETERS_CORRELATION
# PARAMETERS_CORRELATION     filename <- "ps/PearsonCorrelation.eps"
# PARAMETERS_CORRELATION     postscript(filename)
# PARAMETERS_CORRELATION     my.plotcorr(PearsonCor, col=colors[((PearsonCor + 1)/2) * 100], main='Pearson correlations', lower.panel='none', diag='ellipse',mar=0.0 + c(1, 0, 2, 2),outline=TRUE)
# PARAMETERS_CORRELATION     image.plot(legend.only=TRUE,col=colors, zlim=c(-1,1),horizontal=TRUE,legend.mar=6.4,legend.shrink=0.4,legend.width=0.4 )
# PARAMETERS_CORRELATION     dev.off()
     # Save out latex table
     #tab <- xtable(PearsonCor, caption= "Pearson Correlation")
     #print(tab,file="PearsonCorrelation.tex",append=F,table.placement = "h", caption.placement="bottom")

cat(
"\\begin{figure}[!hbp]
{\\includegraphics[width=4.0in,angle=-90]{",file=ReportTexFID)
	cat(filename,file=ReportTexFID)
	cat(
"}}
\\caption{Pearson linear correlation. Thick outline indicates correlations with a p-value lower than 0.05. P1:epsi or kappa, P2 : r0, P3 : alpha}
\\label{fig:r13}
\\end{figure}\n\n",file=ReportTexFID)

     } # End if length(optimpotindex)>1
   } # End if MC
   colfunc = colorRampPalette(c("black","white"))
  
   DistribPotiSpl <- list() 
   DistribALLSpl <- list() 
      for (i in OptimPotIndex) { 
         # Number of accepted iterations for each potential
         NAccIter <- 0
         for (run in 1:n_of_runs) {
           if (ParamAllIterations[[run]][[i]]$status=="Accepted") {
              NAccIter <- NAccIter + 1
           }
         }
   
         # Set colorscale
         colors <- colfunc(NAccIter)
   

         # Set filename
         filename <- paste("ps/",ParamAllIterations[[1]][[i]]$dist_name,"_",optim,".eps",sep='')
         postscript(filename)
   
 
         cat(
"\\begin{figure}[!hbp]
{\\includegraphics[width=4.0in,angle=-90]{",file=ReportTexFID)
	cat(filename,file=ReportTexFID)
	cat(
"}}
\\caption{(a) Distributions, colobar =  Average Loss function (lower values better fit). For clarity, only accepted distributions are shown. (b) Distribution with associated 
minimum average loss function shown in red, remaining top 4 are shown in green. (c) Loss function, in black accepted iterations. (d) RMS loss function, in black accepted iterations.}
\\label{fig:r13}
\\end{figure}\n\n",file=ReportTexFID)

         # Define Multiplot
         par(mfrow=c(2,2),mai=c(0.4,0.4,0.4,0.4),mgp=c(1.2,0.4,0.0),oma=c(1,1,3,1))
   
         # Plot 1. Distributions colored by average Loss and status == Accepted + Ref distribution
         filename <- paste('INPUT/reference_distributions/',ParamAllIterations[[1]][[i]]$ref_dist_fname,sep="")
         #cat ("Loading :",filename,"\n")
         distrib <- as.data.frame(read.table(filename))
         qunit <- ParamAllIterations[[1]][[i]]$q_unit
         if (ParamAllIterations[[1]][[i]]$q_unit=="deg") {
           xRange <- c(ParamAllIterations[[1]][[i]]$xmin,ParamAllIterations[[1]][[i]]$xmax)*180/pi
         } else {
           xRange <- c(ParamAllIterations[[1]][[i]]$xmin,ParamAllIterations[[1]][[i]]$xmax)
         }
         plot(distrib[,1],distrib[,3]/max(distrib[,3]),type="p",pch=1,xlab="r (A)",ylab="P(r)",main="(a)",xlim=xRange,ylim=c(0,1))
         referenceDistrib <- list(x=distrib[,1],y=distrib[,3]/max(distrib[,3]));
         image.plot(legend.only=TRUE, zlim=range(AvgLoss[2:length(AvgLoss)]), col=colors, args.legend = list(x = "topleft"))
         cindex<-0
         for (run in avgSort$ix) {
           if (ParamAllIterations[[run]][[i]]$status=="Accepted") {
              cindex<-cindex+1
              ParamAllIterations[[run]][[i]]$dist_fname
              filename <- paste('OUTPUT/r',run-1,'/dists/Statistics/PDB/',ParamAllIterations[[1]][[i]]$dist_fname,sep="")
              distrib <- as.data.frame(read.table(filename))
              distrib.spl <- spline(distrib[,1],distrib[,3]/max(distrib[,3]),n=2*length(distrib[,1]))
              DistribPotiSpl[[cindex]] <- distrib.spl
              lines(distrib.spl,col=colors[cindex])
           }
         }
                
         # Plot 2. Distribution with associated min(AvgLoss) + Top five distribution + Ref. distributions 
         # These distributions are ranked by average loss so they may have status Rejected
         filename <- paste('INPUT/reference_distributions/',ParamAllIterations[[1]][[i]]$ref_dist_fname,sep="")
         distrib <- as.data.frame(read.table(filename))
         plot(distrib[,1],distrib[,3]/max(distrib[,3]),type="p",pch=1,xlab="r (A)",ylab="P(r)",main="(b)",xlim=xRange,ylim=c(0,1))
         # min(AvgLoss)
         filename <- paste('OUTPUT/r',(avgSort$ix[1]-1),'/dists/Statistics/PDB/',ParamAllIterations[[avgSort$ix[1]]][[i]]$dist_fname,sep="")
         distrib <- as.data.frame(read.table(filename))
         distrib.spl <- spline(distrib[,1],distrib[,3]/max(distrib[,3]),n=2*length(distrib[,1]))
         lines(distrib.spl,col=2)
         # Remaining Top 5 AvgLoss
         last <- min(5,length(AvgLoss)-1)
         for (ii in 2:last) {
           filename <- paste('OUTPUT/r',(avgSort$ix[ii]-1),'/dists/Statistics/PDB/',ParamAllIterations[[avgSort$ix[1]]][[i]]$dist_fname,sep="")
           distrib <- as.data.frame(read.table(filename))
           distrib.spl <- spline(distrib[,1],distrib[,3]/max(distrib[,3]),n=2*length(distrib[,1]))
           lines(distrib.spl,col="#5aff5e",lty=2)
         }
      
         # Plot 3. Loss function all iterations
         if (optim=="MC") {
            plot(c(2:nrow(LossFrame)),LossFrame[(2:nrow(LossFrame)),i],xlab="Iterations",ylab="Loss",main="(b)",col="#CCCCCC") # all iterations
            accrun<-1
            LossFrameAccepted <- vector()
            AcceptedIterations <- vector()
            for (run in 2:nrow(LossFrame)) {
              if (ParamAllIterations[[run]][[i]]$status=="Accepted") {
                accrun <- accrun + 1
                LossFrameAccepted[accrun] <- LossFrame[run,i]
                AcceptedIterations[accrun] <- run
                points(run,LossFrame[run,i]) # only accepted  (if MC)
              }
            }
         } else if (optim=="IBI") {
            LossFrameAccepted <- vector()
            AcceptedIterations <- vector()
            AcceptedIterations <- c(2:nrow(LossFrame))
            LossFrameAccepted <- LossFrame[(2:nrow(LossFrame)),i]
            plot(c(2:nrow(LossFrame)),LossFrame[(2:nrow(LossFrame)),i],xlab="Iterations",ylab="Loss",main="(b)",col="black") # all iterations
         }
      
         # Plot 4. RMS(loss) function
         RMSLossFrame <- diff(LossFrame[(2:nrow(LossFrame)),i])^2
         if (optim=="MC") {
           plot(c(1:length(RMSLossFrame)),RMSLossFrame,xlab="Iterations",ylab="RMS(Loss)",main="(c)",col="#CCCCCC",ylim=c(0,1)) # all iterations
           for (run in 1:length(RMSLossFrame)) {
             if (ParamAllIterations[[run]][[i]]$status=="Accepted") {
                 points(run,RMSLossFrame[run]) # only accepted (if MC)
     	     }
           }
         } else if (optim=="IBI") {
           plot(c(1:length(RMSLossFrame)),RMSLossFrame,xlab="Iterations",ylab="RMS(Loss)",main="(c)",col="black") # all iterations
         }
      
         title(paste("Interaction :",ParamAllIterations[[1]][[i]]$dist_name), outer = TRUE)
         dev.off()

         param <- list(name=DistNames[[i]], distributions=DistribPotiSpl, acceptedIterations=AcceptedIterations, lossAcc=LossFrameAccepted, ranges=xRange, qunit=qunit, referenceDistrib=referenceDistrib)
         DistribALLSpl[[length(DistribALLSpl)+1]] <- param
         save(DistribALLSpl,file="DistribALLSpl")

      } # end potential loop

#close(ReportTexFID) 
} # end function
