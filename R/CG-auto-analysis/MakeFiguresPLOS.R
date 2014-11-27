# options(echo=TRUE) # if you want see commands in output file
# args <- commandArgs(trailingOnly = TRUE)
# print(args)


lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

morse <- function(x,c0,c1,c2) {
  f <- c0 * ( (1.0 - exp(-c2*(x-c1)) )^2 - 1.0 )
  return(f)
}

Cosine <- function(x,c0,c1) {
  f <- c0 * (1.0 + cos(1.0 * x - c1))
  return(f)
}

HarmonicCosine <- function(x,c0,c1) {
  f <- 0.5 * c0 * (cos(x) - cos(c1))^2
  return(f)
}


Rk <-   1.9858775e-3 # kcal mol^-1 K^-1
Temp <- 300 # K



library("XML") # required for XML parsing
library("gtools") # required for mixedsort
library("lattice") # required for levelplot
library("xtable") # required for latex output
library("fields") # required for plot colorscale
library("corrplot") # 
library("ggplot2")
library("gridExtra")


# # SET SYSTEM PARAMs
# if (length(args)<3) {
#  stop("Usage : MakePLOSFig IBI MC Alpha")
# }
# IBI <- as.logical(args[1])
# MC <- as.logical(args[2])
# Alpha <- as.logical(args[3])
IBI<-TRUE
MC<-FALSE
Alpha<-TRUE


# Distribution in the same order than in Main.XML
#DistributionsType <- list("r13","")
# Load Best iteration MC
if (MC) {
  BestMC <- xmlInternalTreeParse("OUTPUT/BestIteration_Output.xml");
  BestIterationN <- as.numeric(xpathApply(BestMC, "//c_run",xmlValue))+1;
  BestIteration <- paste("OUTPUT/r",BestIterationN,sep="")
  print(BestIteration);
}


# Loop over runs
DistribALLSpl <- list()
files <- list.files(path="OUTPUT", pattern="r*/*_Input.xml", full.names=T, recursive=TRUE)
SortedFiles<-mixedsort(files)
n_of_runs <- length(files)

# Load Reference distributions


print(IBI)
print(MC)
print(Alpha)

MCIBI<-FALSE
if (MC&&IBI) {
  MCIBI<-TRUE
  MC<-FALSE
  IBI<-TRUE
}


if (Alpha) { 
  #potlist<-c("null","r12","r13","r14","theta","phi","vdw")  
  potlist<-c("null","r12","theta","vdw")  
  #potlist<-c("null","r12","vdw")  
} else {
  potlist<-c("null","r12","r13","theta","phi","vdw")  
}

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#  1.  Load distributions and losses
doc0 <- xmlInternalTreeParse(SortedFiles[1]);
src0 <- xpathApply(doc0, "//input/param")
NumberOfPotentials <- xmlSize(src0)
bOptimize <- vector()
refDistribs <- list()
BestMCDistribs <- list()
ipotindex <- vector()
# Load ref distributions and bOptimize vector
for (ipot in 1:NumberOfPotentials) {
  TmpDoc0 <- xmlDoc(src0[[ipot]])
  bOptimize <- append(bOptimize,xpathApply(TmpDoc0,"//bOptimize",xmlValue)[[1]])
  if (bOptimize[ipot] == "true") {
   ipotindex <- append(ipotindex,ipot)
   filename <- paste('INPUT/Param',ipot-1,'.dat',sep="")
   refDistrib <- as.data.frame(read.table(filename))
   if (MC) {
    filename <- paste(BestIteration,'/Param',ipot-1,'.dat',sep="")
    BestMCDistrib <- as.data.frame(read.table(filename))
   } 
   if ((potlist[ipot]=="theta")||(potlist[ipot]=="phi")) {
     refDistrib[,1] <- refDistrib[,1]/pi*180.0
     if (MC) {
      BestMCDistrib[,1] <- BestMCDistrib[,1]/pi*180.0
     }
   }
   refDistrib.spl <- spline(refDistrib[,1],refDistrib[,2]/(sum(refDistrib[,2])*diff(refDistrib[,1])[1]),n=2*length(refDistrib[,1]))
   refDistribs <- lappend(refDistribs,refDistrib.spl)
   if (MC) {
    BestMCDistrib.spl <- spline(BestMCDistrib[,1],BestMCDistrib[,2]/(sum(BestMCDistrib[,2])*diff(BestMCDistrib[,1])[1]),n=2*length(BestMCDistrib[,1]))
    BestMCDistribs <- lappend(BestMCDistribs,BestMCDistrib.spl)
   } 
  }
}
AllLossFrame <- as.data.frame(matrix(NA, n_of_runs, length(ipotindex)))
ipotOpti <- 0
ipotNparams <- vector()
MatrixParams <- matrix(0,(n_of_runs-1),byrow=TRUE)
for (ipot in 1:NumberOfPotentials) {
  distribs.spl <- list()
  FittedPMFs <- list()
  ToBeFittedPMFs <- list()
  losses <- vector()    
  Medians <- vector()    
  acceptedIterations <- vector()
  if (bOptimize[ipot] == "true") {
    ipotOpti <- ipotOpti + 1
    params <- vector()
    for (run in 1:n_of_runs) {
      doc <- xmlInternalTreeParse(SortedFiles[run]);
      src <- xpathApply(doc, "//input/param")
      TmpDoc <- xmlDoc(src[[ipot]])
      #bOptimize <- xpathApply(TmpDoc,"//bOptimize",xmlValue)[[1]]
      #ii <- 0  
      # Print out only if optimized potentials
      #if(bOptimize=="true") {
       if (MC) {
          status <- xpathApply(TmpDoc,"//Accepted",xmlValue)[[1]] 
       }  else {
         status <- "true"
       }
       # Read Distribution
       filename <- paste('OUTPUT/r',run-1,'/Param',ipot-1,'.dat',sep="")
       distrib <- as.data.frame(read.table(filename))
       # Read Fitted PMF (only for IBI)
       if (IBI)  {
        filename <- paste('OUTPUT/r',run-1,'/Param',ipot-1,'.dat.fitted.PMF.dat',sep="")
        Fitted <- as.data.frame(read.table(filename))
        filename <- paste('OUTPUT/r',run-1,'/Param',ipot-1,'.dat.tobefitted.PMF.dat',sep="")
        ToBeFitted <- as.data.frame(read.table(filename))
       }           
       # And compute median 
       forMedian <- which(cumsum(distrib[,2]/sum(distrib[,2]))<0.5)
       Median <- distrib[c(length(forMedian)),1]
       Medians <- append(Medians,Median)

       # And for only accepted distributions (always true if IBI is used)
        if (status=="true") {
            if ((potlist[ipot]=="theta")||(potlist[ipot]=="phi")) {
             distrib[,1] <- distrib[,1]/pi*180.0    
             Fitted[,1] <- Fitted[,1]/pi*180.0    
             ToBeFitted[,1] <- ToBeFitted[,1]/pi*180.0
            }
          
          distrib.spl <- spline(distrib[,1],distrib[,2]/(sum(distrib[,2])*diff(distrib[,1])[1]),n=2*length(distrib[,1]))
          distribs.spl <- lappend(distribs.spl,distrib.spl)
          Fitted.spl <- spline(Fitted[,1],Fitted[,2],n=2*length(Fitted[,1]))
          FittedPMFs <- lappend(FittedPMFs,Fitted.spl)    
          ToBeFittedPMFs <- lappend(ToBeFittedPMFs,ToBeFitted)
          #lineWidth <- 5.0        
          loss <- as.numeric(xpathApply(TmpDoc,"//LossFunction",xmlValue)[[1]])
          losses <- append(losses,loss)   
          acceptedIterations <- append(acceptedIterations,run)          
        }
      # only for MC, these will be used to get the best parameters
      if (MC) {
        AllLossFrame[run,ipot] <- as.numeric(xpathApply(TmpDoc,"//LossFunction",xmlValue)[[1]])
        paramipot <- as.numeric(xpathApply(TmpDoc,"//parameter_try",xmlValue))
        params <- append(params,paramipot)        
      } else {                  
          paramipot <- as.numeric(xpathApply(TmpDoc,"//parameter",xmlValue))
          params <- append(params,paramipot)               
      }
    } # iterations loop
    Nparams <- round(length(params)/(n_of_runs-1)) 
    ipotNparams[ipotOpti] <- Nparams
    if (MC) {
      MatrixParams <- matrix(params,(n_of_runs-1),byrow=TRUE) 
    }   else {
      MatrixParams <- matrix(params,(n_of_runs),byrow=TRUE) }
    if (ipotOpti==1) {
      MatrixAllParams <- MatrixParams }
    else {
      MatrixAllParams <- cbind(MatrixAllParams,MatrixParams)
    }
    DistAndLosses <- list(dists=distribs.spl,losses=losses,Medians=Medians,acceptedIterations=acceptedIterations,params=MatrixParams,FittedPMFs=FittedPMFs,ToBeFittedPMFs=ToBeFittedPMFs)
    DistribALLSpl <- lappend(DistribALLSpl,DistAndLosses)
    
  } # bOptimize if
} # potentials loop

# dists and pot
# DEBUG
Debug <- TRUE
if (Debug) {
 ranges <- c()
 colfunc = colorRampPalette(c("white","black"))
 for (ipot in 1:length(DistribALLSpl)) {
  
   #REFERENCE
   refdist <- as.data.frame(refDistribs[[ipot]])
   refdist[,2] <- refdist[,2]/(sum(refdist[,2])*diff(refdist[,1])[1])
   #refdist <- as.data.frame(list(x=c(0.0,1.0,2.0,3.0,4.0,5.0),y=c(0.0,0.1,0.7,0.1,0.1,0.0)))   
   refdist[,3] <- - Rk * Temp * log(refdist[,2])
   
   
   colnames(refdist) <- c("q","pref","uref")
   
   ggplBASE <- ggplot(data=refdist) +
    theme_bw() +
    theme(#axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.x=element_blank(),
      plot.margin = rep(unit(0,"null"),4),
      panel.margin = unit(0,"null"))
          #labs(x=NULL))
          #axis.ticks.length = unit(0,"null"),
          #axis.ticks.margin = unit(0,"null")) +
 
   ggpl1BASE <- ggplBASE + geom_line(aes(x=q,y=pref),size=2,colour="red") 
   ggpl2BASE <- ggplBASE + geom_line(aes(x=q,y=uref),size=2,colour="red")

   last <- length(DistribALLSpl[[ipot]]$dists)  
   last <- 20
   nplots <- 10
   nblocks <- last/nplots 
   
   plotlist<-list()
   starti <- 0
   endi <- 0
   for (iblock in seq(1,nblocks,by=1) ) {
    starti <- endi + 1
    endi <- starti + nplots - 1
    if (iblock==nblocks) {
      endi <- starti + nplots - 2
    }
    print(starti)
    print(endi)  
    plotlist<-list()
    for (i in seq(starti,endi,by=1) ) {
     ## SET ranges and parameters
     if (Alpha) {
       if (potlist[ipotindex[ipot]]=="r12") {
         ranges[1] <- 5.0
         ranges[2] <- 6.0
         ranges[3] <- -3.0
         ranges[4] <- 3.0         
       } else  if (potlist[ipotindex[ipot]]=="r13") {
         ranges[1] <- 4.0
         ranges[2] <- 6.0
         ranges[3] <- 0.0
         ranges[4] <- 6.0
       } else  if (potlist[ipotindex[ipot]]=="r14") {
         ranges[1] <- 5.0
         ranges[2] <- 7.0
         ranges[3] <- 0.0
         ranges[4] <- 6.0
       } else  if (potlist[ipotindex[ipot]]=="theta") {
         ranges[1] <- 60.0
         ranges[2] <- 120.0
         ranges[3] <- 0.0
         ranges[4] <- 6.0
       } else  if (potlist[ipotindex[ipot]]=="phi") {
         ranges[1] <- 25.0
         ranges[2] <- 75.0
         ranges[3] <- 0.0
         ranges[4] <- 6.0
       }  } else {
         if (potlist[ipotindex[ipot]]=="r12") {
           ranges[1] <- 4.0
           ranges[2] <- 8.0
         } else  if (potlist[ipotindex[ipot]]=="r13") {
           ranges[1] <- 4.0
           ranges[2] <- 8.0    
         } else  if (potlist[ipotindex[ipot]]=="theta") {
           ranges[1] <- 60.0
           ranges[2] <- 120.0
         } else  if (potlist[ipotindex[ipot]]=="phi") {
           ranges[1] <- 20.0
           ranges[2] <- 120.0
         }
       }
     dist <- as.data.frame(DistribALLSpl[[ipot]]$dists[[i]])
     dist[,2] <- dist[,2]/(sum(dist[,2])*diff(dist[,1])[1])
     dist[,3] <- - Rk * Temp * log(dist[,2])
     Fitted <- as.data.frame(DistribALLSpl[[ipot]]$FittedPMFs[i])
     #ind1<- (dist[,1]>ranges[1]) & (dist[,1]<ranges[2])
     #ind2<- (Fitted$x>ranges[1]) & (Fitted$x<ranges[2])
     #shift <- (min(dist[ind1,3],na.rm=TRUE) - min(Fitted$y[ind2],na.rm=TRUE))
     #Fitted$y <- Fitted$y+shift
     
     ToBeFitted <- as.data.frame(DistribALLSpl[[ipot]]$ToBeFittedPMFs[i])
        
     if (i==1) {
       Uprev <- refdist[,3]
     } else {
       dprev <- as.data.frame(DistribALLSpl[[ipot]]$dists[[i-1]])
       dprev[,2] <- dprev[,2]/(sum(dprev[,2])*diff(dprev[,1])[1])
       Uprev <- -Rk * Temp * log(dprev[,2])
     }
     
     Unext <- as.data.frame(dist[,1])
     Unext[,2] <- Uprev + Rk * Temp * log(dist[,2]/refdist[,2])
     
     colnames(dist) <- c("q","p","u")
     colnames(Unext) <- c("q","unext")
     iterval <- i-1
     ggpl1 <- ggpl1BASE + geom_line(data=dist,aes(x=q,y=p),size=2) +
       scale_x_continuous(limits=c(ranges[1],ranges[2])) + 
       scale_y_continuous(name="p") 
     xlabel<-ranges[2]-ranges[1]/100*10
     ylabel<-ranges[4]-ranges[2]/100*5
     ggpl2 <- ggpl2BASE + geom_line(data=dist,aes(x=q,y=u),size=2) +        
       geom_line(data=ToBeFitted,aes(x=V1,y=V2),lty=1,colour="yellow",size=2) + 
       geom_line(data=Fitted,aes(x=x,y=y),lty=2,colour="green",size=2) +        
       scale_x_continuous(limits=c(ranges[1],ranges[2])) + 
       scale_y_continuous(name="u",limits=c(ranges[3],ranges[4])) +
       annotate("text",x=xlabel,y=ylabel,label=iterval) 
     
     plotlist <-lappend(plotlist,ggpl1)
     plotlist <-lappend(plotlist,ggpl2)
   }
   do.call(grid.arrange,c(plotlist,ncol=4))
   filename<-paste("FigForPaper/Debug/plot",ipot,"block",iblock,".eps",sep="")
   ggsave(filename,do.call(arrangeGrob,c(plotlist,ncol=4)))
  } # END BLOCKS LOOP
 } # End potentials loop
} #END IF DEBUG


#  2. Figures about distributions
colfunc = colorRampPalette(c("white","black"))
reflist<-list()
for (ipot in 1:length(DistribALLSpl)) {
  #par(fig=c(0.0,1.0,0.0,1.0),mar=c(8,8,1,2),mgp=c(1.0,1.4,0.0),oma=c(0.1,0.1,0.1,0.1),new=FALSE)
  NaccIter <- length(DistribALLSpl[[ipot]]$losses)
  loss <- DistribALLSpl[[ipot]]$losses
  SortLoss <- sort(loss,index.return=TRUE)
  #if (IBI) {
  #  SortLoss$ix<-c(length(loss):1)
  #}
  colors <- colfunc(NaccIter)
  #lineWidth<-5.0/(xx^0.15)
  lineWidth<-seq(1.0,7.0,length.out = NaccIter )
  ranges <- c()
  if (Alpha) {
   if (potlist[ipotindex[ipot]]=="r12") {
    ranges[1] <- 4.0
    ranges[2] <- 8.0
   } else if (potlist[ipotindex[ipot]]=="r13") {
    ranges[1] <- 4.0
    ranges[2] <- 8.0
   } else if (potlist[ipotindex[ipot]]=="r14") {
    ranges[1] <- 4.0
    ranges[2] <- 8.0
   } else if (potlist[ipotindex[ipot]]=="theta") {
    ranges[1] <- 60.0
    ranges[2] <- 120.0
   } else if (potlist[ipotindex[ipot]]=="phi") {
    ranges[1] <- 20.0
    ranges[2] <- 120.0
   }  
  } else {
   if (potlist[ipotindex[ipot]]=="r12") {
      ranges[1] <- 4.0
      ranges[2] <- 8.0
   } else if (potlist[ipotindex[ipot]]=="r13") {
      ranges[1] <- 4.0
      ranges[2] <- 8.0    
   } else if (potlist[ipotindex[ipot]]=="theta") {
      ranges[1] <- 60.0
      ranges[2] <- 120.0
   } else if (potlist[ipotindex[ipot]]=="phi") {
      ranges[1] <- 20.0
      ranges[2] <- 120.0
   }         
  }  
  
  filename <- paste("FigForPaper/Param",ipotindex[ipot]-1,".eps",sep="")
  #tiff(filename, width = 1200, height = 1200, units = 'px', compression = c("none") )
  postscript(filename,height = 5, width=10)
  par(mar=c(7,5,1,1),mgp=c(5,2,0),oma=c(0.0,0.0,0.0,0.0),new=FALSE)
  refDistrib <- refDistribs[[ipot]]

  last <- length(DistribALLSpl[[ipot]]$dists)
  step <- round(length(DistribALLSpl[[ipot]]$dists) / 11)
  # Find max y
  maxyv<-vector()
  maxyv[1]<-max(refDistrib$y)
  for (i in seq(1,last,by=step) ) {
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    maxyv[i+1]<-max(distrib$y)
  }
  if (MC) {
    maxy<-max(maxyv,na.rm=TRUE)
  } else {
    if (Alpha)
      maxy<-median(maxyv,na.rm=TRUE) + 0.28*median(maxyv,na.rm=TRUE)
    else
      maxy<-max(maxyv,na.rm=TRUE)
  }
  if (MCIBI) {
    plot(refDistrib$x,refDistrib$y,type="l",lty=2,col=cbbPalette[7],lwd=5.0, xlab="", ylab="", xlim=ranges, ylim=c(0,maxy), cex.axis=3.4, cex.lab=3.5, frame=TRUE, axes=FALSE) # MCSA-IBI
    axis(side = 1, tick = TRUE, font=2, cex.axis=3.0, cex.lab=2.2) 
    #axis(side = 2, tick = TRUE, at=c(0.5,1.0), font=2, cex.axis=3.0, cex.lab=2.2)  
  } else {
    plot(refDistrib$x,refDistrib$y,type="l",lty=2,col=cbbPalette[7],lwd=5.0, xlab="", ylab="", xlim=ranges, ylim=c(0,maxy), cex.axis=3.4, cex.lab=3.5, frame=TRUE, axes=FALSE)
    #axis(side = 2, tick = TRUE, at=c(0.5,1.0), font=2, cex.axis=3.0, cex.lab=2.2)  
  }
  for (i in seq(1,last,by=step) ) {
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    #lineWidth <- lineWidth - 3.0/NaccIter
    #lines(distrib$x,distrib$y,col=colors[SortLoss$ix[i]],lwd=lineWidth[SortLoss$ix[i]])
    lines(distrib$x,distrib$y,col=colors[i],lwd=lineWidth[i])
  } # End iteration loop for distributions
  # Plot in red the reference one
  lines(refDistrib$x,refDistrib$y,type="l",lty=2,col=cbbPalette[7],lwd=7.0)
  # Plot in green the last one
  if (MC)  {   
    distrib <- BestMCDistribs[[ipot]]
    lines(distrib$x,distrib$y,type="l",lty=4,col=cbbPalette[4],lwd=7.0)     
  } else {  
    distrib <- DistribALLSpl[[ipot]]$dists[[last-1]]
    lines(distrib$x,distrib$y,type="l",lty=4,col=cbbPalette[4],lwd=7.0)         
  }
  
  dev.off()
}

# 3. Figures about Losses correlations (Lossfunction correlations) 
#    It also generate the matrix mt.scaled to plot the scaled loss function (section 5. Figures about Loss function )
#    Only for IBI and OnlyIBI 
if (!MC) {
  m<-matrix(0,length(DistribALLSpl),length(DistribALLSpl[[1]]$losses))
  for (ipot in 1:length(DistribALLSpl)) { 
    m[ipot,]<- DistribALLSpl[[ipot]]$losses
    
  }
  
  mt <- t(m)
  #mt.scaled<-scale(mt)
  mt.scaled<-mt
  
  
  # Shift to positve value only
  abmin <- min(mt.scaled)
  mt.scaled<- mt.scaled + abs(abmin)
  
  # Shift by 1 each for better visualization
  for (ipot in 1:length(DistribALLSpl)) 
    mt.scaled[,ipot] <- mt.scaled[,ipot] + 1.0*ipot
  
  minx <- 1
  maxx <- length(DistribALLSpl[[1]]$losses)
  miny <- min(mt.scaled)
  maxy <- max(mt.scaled)
  
  if (ipot==5) {
    names<-c("","","","","")
  } else if (ipot==4)  {
    names<-c("","","","")
  }  else if (ipot==2)  {
    names<-c("","")
  }  else if (ipot==1)  {
    names<-c("")
  }
  colnames(mt.scaled) <- names     
#NOTUSED#  Spearman <- cor(mt.scaled[,],method="spearman")    
#NOTUSED#  #Spearman <- abs(Spearman)
#NOTUSED#  #Spearman[lower.tri(Spearman,diag=T)] = 0.0 # set 0 diagonal and lower tri matrix
#NOTUSED#  pmat <- matrix(0,dim(Spearman)[1],dim(Spearman)[2])
#NOTUSED#  # Add p-values to the lower triangular part of Pearson matrix
#NOTUSED#  for (i in 1:dim(mt.scaled)[2]) for (j in i:dim(mt.scaled)[2]) {
#NOTUSED#    pvalue <- cor.test(mt.scaled[,i],mt.scaled[,j],alternative="two.side",method="spearman")$p.val
#NOTUSED#    pmat[i, j] <- pmat[j, i] <- pvalue
#NOTUSED#    #if (pvalue<0.01) {  Spearman[j,i] <- 0.01 }
#NOTUSED#  }
#NOTUSED#  filename <- "FigForPaper/SpearmanCorrelation.eps"
#NOTUSED#  postscript(filename) 
#NOTUSED#  #my.plotcorr(Spearman, col=colors[((Spearman + 1)/2) * 100], main='Spearman correlations', lower.panel='none', diag='none',mar=0.0 + c(1, 0, 2, 2),outline=TRUE)
#NOTUSED#  #image.plot(legend.only=TRUE,col=colors, zlim=c(-1,1),horizontal=TRUE,legend.mar=6.4,legend.shrink=0.4,legend.width=0.4 )
#NOTUSED#  #corrplot(Spearman)  
#NOTUSED#  # Reverse significant with insignificant, just for visualizing reasons
#NOTUSED#  pmm <- pmat
#NOTUSED#  pmm[pmat>0.01] <- 0.0
#NOTUSED#  pmm[pmat<=0.01] <- 1.0
#NOTUSED#  corrplot(Spearman,p.mat=pmm,sig.level =0.01,method="color",type="lower",diag=F,cl.cex=1.4,pch="*")
#NOTUSED#  dev.off()  
#NOTUSED#  caption <- paste("Correlations among loss functions",sep="")
#NOTUSED#  tab <- xtable((Spearman[,]), caption=caption)
#NOTUSED#  filename <- "FigForPaper/SpearmanCorrelation.tex"
#NOTUSED#  print(tab,file=filename,append=F,table.placement = "h", caption.placement="bottom")  
}

# # 4. Figures about correlations (Distribution correlations)
# #    Only for IBI and OnlyIBI 
# #if (!MC) {
#   mtp<-matrix(0,length(DistribALLSpl),length(DistribALLSpl[[1]]$Medians))
#   for (ipot in 1:length(DistribALLSpl)) 
#     mtp[ipot,]<- DistribALLSpl[[ipot]]$Medians
#   
#   mtp <- t(mtp)
# 
#   if (ipot==5) {
#     names<-c("","","","","")
#   } else  {
#     names<-c("","","","")
#   }
#   colnames(mtp) <- names     
#   Spearman <- cor(mtp[,],method="spearman")    
#   pmat <- matrix(0,dim(Spearman)[1],dim(Spearman)[2])
#   # Add p-values to the lower triangular part of Pearson matrix
#   for (i in 1:dim(mtp)[2]) for (j in i:dim(mtp)[2]) {
#     pvalue <- cor.test(mtp[,i],mtp[,j],alternative="two.side",method="spearman")$p.val
#     pmat[i, j] <- pmat[j, i] <- pvalue
#     #if (pvalue<0.01) {  Spearman[j,i] <- 0.01 }
#   }
#   filename <- "FigForPaper/DistributionCorrelation.eps"
#   postscript(filename) 
#   pmm <- pmat
#   pmm[pmat>0.01] <- 0.0
#   pmm[pmat<=0.01] <- 1.0
#   corrplot(Spearman,p.mat=pmm,sig.level =0.01,method="color",type="lower",diag=F,cl.cex=1.4,pch="*")
#   dev.off()  
#   caption <- paste("Correlations among distributions",sep="")
#   tab <- xtable((Spearman[,]), caption=caption)
#   filename <- "FigForPaper/DistributionCorrelation.tex"
#   print(tab,file=filename,append=F,table.placement = "h", caption.placement="bottom")  
# #}



# 5. Figures about Loss function
filename <- paste("FigForPaper/LossFunction.eps")
postscript(filename, height=5, width=10)
par(mar=c(7,8,2,1),mgp=c(4,1.5,0),oma=c(0.0,0.0,0.0,0.0),new=FALSE)

# Find x-y ranges Form MC only
if (MC) {
  minx<-1
  maxx<-0
  miny<-1000.0
  maxy<-0.0
  for (ipot in 1:length(DistribALLSpl)) {
    loss <- DistribALLSpl[[ipot]]$losses
    acceptedIterations <- DistribALLSpl[[ipot]]$acceptedIterations
    if (max(acceptedIterations,na.rm=TRUE)>maxx)
      maxx <- max(acceptedIterations,na.rm=TRUE)
    
    if (max(loss,na.rm=TRUE)>maxy)
      maxy <- max(loss,na.rm=TRUE) 
    
    if (min(loss,na.rm=TRUE)<miny)
      miny <- min(loss,na.rm=TRUE)    
  }  
}
# Do the first plot
if (MC) { 
  if (Alpha)
    maxy <- maxy + 4.0
  else
    maxy <- maxy + 3.5
  plot(1, type="n", axes=FALSE, frame=TRUE, xlab='', ylab='', xlim=c(minx,maxx), ylim=c(miny,maxy), cex.axis=3.0, cex.lab=3.0) # MC  
}else{
  if (Alpha) {
    plot(1, type="n", axes=FALSE, frame=TRUE, xlab='', ylab='',  xlim=c(minx,maxx), ylim=c(miny,maxy+5.0), cex.axis=3.0, cex.lab=3.0) # IBI and OnlyIBI
  } else {
    plot(1, type="n", axes=FALSE, frame=TRUE, xlab='', ylab='',  xlim=c(minx,maxx), ylim=c(miny+1,maxy+5.0), cex.axis=3.0, cex.lab=3.0) # IBI and OnlyIBI
  }
}
axis(side = 1, tick = TRUE, font=2, cex.axis=3.0, cex.lab=2.2) 
axis(side = 2, tick = FALSE, at=seq(miny,maxy), font=2, cex.axis=3.0, cex.lab=2.2, labels=FALSE)  
mtext("Iteration", side=1, line=4, cex=3.5)
#mtext("[a.u.]", side=2, line=1.0, cex=3.5)
# Add overalys
for (ipot in 1:length(DistribALLSpl)) {
  
  if (MC) 
    loss <- DistribALLSpl[[ipot]]$losses # MC
  else
    loss <- mt.scaled[,ipot]+abs(miny) # IBI and OnlyIBI

  
  Acc <- DistribALLSpl[[ipot]]$acceptedIterations
  # Set colors and point type
  if (ipot==1) {
    color<-"steelblue2"
    color<-cbbPalette[1]
    ptype<-ipot+20
  } else if (ipot==2) {
    color<-"#d50000"
    color<-cbbPalette[2]
    ptype<-ipot+20
  } else if (ipot==3) {
    if (length(DistribALLSpl)==5) {
      color<-"#00d254"
      color<-cbbPalette[4]
      ptype<-ipot+20
    } else {
      color<-"#e56cff"
      color<-cbbPalette[3]
      ptype<-ipot+21
    }
  } else if (ipot==4)  {
    if (length(DistribALLSpl)==5) {
      color<-"#e56cff"
      color<-cbbPalette[3]
      ptype<-ipot+20
    } else {
      color<-"#fffd43"
      ptype<-ipot+21
    }
  } else if (ipot==5) {
    color<-"#fffd43"
    color<-cbbPalette[5]
    ptype<-ipot+20
  }
  #points(Acc,loss+1*ipot-1, pch=ptype, bg=color, cex=2.0)
  lines(Acc,loss+1*ipot-1, col="black", lwd=6.5)
  lines(Acc,loss+1*ipot-1, col=color, lwd=5.0)  
  #symbols(Acc,loss,circles=rep(0.5,length(Acc)),inches=1/8,ann=F,bg=color,add=TRUE)
} # end ipot loop
dev.off()


# 6.  Best Parameters selection (For montecarlo only)
if (MC) {
  # Best params based on SqSum and Var
  Vm <- vector()
  SqSum <- vector()
  VarM <- vector()
  VarScaled <- vector()
  for (run in 1:nrow(AllLossFrame)) {
    SqSum[run] <- sum(AllLossFrame[run,]^2,na.rm=TRUE)
    #if (length(OptimPotIndex)>1) 
      VarM[run] <- sd(AllLossFrame[run,],na.rm=TRUE)^2
    #else
    #  VarScaled[run] <- 0.0
  }
  SqScaled<-(SqSum-min(SqSum))/(max(SqSum)-min(SqSum))
  #if (length(OptimPotIndex)>1) 
    VarScaled<-(VarM-min(VarM))/(max(VarM)-min(VarM))
  Vm <- SqScaled + VarScaled
  
  # Sort 
  #avgSort <- sort(AvgLoss,index.return=TRUE) 
  avgSort <- sort(Vm,index.return=TRUE) 
  SortedParametersFrame<-as.data.frame(MatrixAllParams[avgSort$ix[1:min(10,length(avgSort$ix))]-1,])
  #s[,10] <- (s[,10]+3.14)*180/3.14
  #s[,8] <- s[,8]*180/3.14
  all_mean <- sapply(SortedParametersFrame,median)
  all_quantiles <- sapply(SortedParametersFrame,quantile)
  
  SortedParametersFrame <- rbind(SortedParametersFrame,all_quantiles[3,])
  SortedParametersFrame <- rbind(SortedParametersFrame,all_quantiles[2,])
  SortedParametersFrame <- rbind(SortedParametersFrame,all_quantiles[4,])
  rownames(SortedParametersFrame)[nrow(SortedParametersFrame)-2]<-c("Median")
  rownames(SortedParametersFrame)[nrow(SortedParametersFrame)-1]<-c("25pct")
  rownames(SortedParametersFrame)[nrow(SortedParametersFrame)]<-c("75pct")
  # Save out latex table
  caption <- "Best parameters based on average loss function. Params are sorted on average loss function and the top 5 (max) are selected. Based on MC-SA simulations."
  tab <- xtable(t(SortedParametersFrame[,]), caption=caption)
  filename <- "BestBasedOnAvgLoss_MCSA.tex"
  print(tab,file=filename,append=F,table.placement = "h", caption.placement="bottom")

  Nstart <- 1
  for (i in 1:ipotOpti) { 
    Nparams <- ipotNparams[i]
    for (j in 1:Nparams) {
      val <- SortedParametersFrame[nrow(SortedParametersFrame)-2,Nstart:Nstart+j-1]
      cat("<parameter>",val,"</parameter>\n")
    }
    for (j in 1:Nparams) {
      val <- SortedParametersFrame[nrow(SortedParametersFrame)-1,Nstart:Nstart+j-1]
      cat("<parameter_min>",val,"</parameter_min>\n")
    }
    for (j in 1:Nparams) {
      val <- SortedParametersFrame[nrow(SortedParametersFrame),Nstart:Nstart+j-1]
      cat("<parameter_max>",val,"</parameter_max>\n")
    }
    Nstart <- Nstart + Nparams
  }
}


# 7. Boxplot of RMS(Loss)
# Use here ggplot2
doit<-FALSE
if(doit) {
if(!MC) {
 filename <- paste("FigForPaper/RMSLossFunction.eps")
 postscript(filename, height=5, width=10)
 par(mar=c(4,4,1,1),mgp=c(5,2,0),oma=c(1.0,1.0,1.0,1.0),new=FALSE)
 DMT <- as.data.frame(diff(mt))
 if (ipot==5) {
   colnames(DMT) <- c("r13","r14","r15","Angle","Dihedral")
 } else {
   colnames(DMT) <- c("r13","r14","Angle","Dihedral")
 }
 boxplot(DMT)
 save(DMT,file="RMSLoss.RData")
 dev.off()
}
}

# AcceptedIterations<-vector()
# # Load distributions and loss functions
# for (run in 1:n_of_runs) {
#   # Parse XML
#   doc <- xmlInternalTreeParse(SortedFiles[run]);
#   src <- xpathApply(doc, "//input/param")
#   NumberOfPotentials <- xmlSize(src)
#   TotalNumberOfOptimizedParameters <- 0
#   NumberOfOptimizedPotentials <- 0
#   ParametersNames <- list()
#   OptimPotIndex <- vector()
#   
#   for (i in 1:NumberOfPotentials) {
#     TmpDoc <- xmlDoc(src[[i]])
#     bOptimize <- xpathApply(TmpDoc,"//bOptimize",xmlValue)[[1]]
#     
#     ii <- 0  
#     if(bOptimize=="true") {
#       ii <- ii + 1
#       OptimPotIndex[ii] <- i
#       status <- xpathApply(TmpDoc,"//Accepted",xmlValue)[[1]] 
#       AcceptedIterations[run] <-
#       filename <- paste('OUTPUT/r',run-1,'/Param',i-1,'.dat',sep="")
#       distrib <- as.data.frame(read.table(filename))
#       distrib.spl <- spline(distrib[,1],distrib[,3]/max(distrib[,3]),n=2*length(distrib[,1]))
#       DistribPotiSpl[[cindex]] <- distrib.spl
#       #DistribALLSpl[[ipot]]$referenceDistrib
#       
# #       param <- list(name=DistNames[[i]], distributions=DistribPotiSpl, acceptedIterations=AcceptedIterations, lossAcc=LossFrameAccepted, ranges=xRange, qunit=qunit, referenceDistrib=referenceDistrib)
# #       DistribALLSpl[[length(DistribALLSpl)+1]] <- param
# #       save(DistribALLSpl,file="DistribALLSpl")
#     }
#   }
# }
