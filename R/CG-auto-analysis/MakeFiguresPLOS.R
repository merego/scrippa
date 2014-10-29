lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

library("XML") # required for XML parsing
library("gtools") # required for mixedsort
library("lattice") # required for levelplot
library("xtable") # required for latex output
library("fields") # required for plot colorscale
library("corrplot") # 

# Distribution in the same order than in Main.XML
#DistributionsType <- list("r13","")

# Loop over runs
DistribALLSpl <- list()
files <- list.files(path="OUTPUT", pattern="r*/*_Input.xml", full.names=T, recursive=TRUE)
SortedFiles<-mixedsort(files)
n_of_runs <- length(files)

# Load Reference distributions

# SET SYSTEM PARAMs
IBI <- FALSE
MC <- TRUE
Alpha<-TRUE 

if ((MC&IBI)) {
  cat ("choose MC, IBI or OnlyIBI !\n") 
  exit()
}

if (Alpha) {
  degpot <- 4
} else {
  degpot <- 3
}


# Load distributions and losses
doc0 <- xmlInternalTreeParse(SortedFiles[1]);
src0 <- xpathApply(doc0, "//input/param")
NumberOfPotentials <- xmlSize(src0)
bOptimize <- vector()
refDistribs <- list()
ipotindex <- vector()
# Load ref distributions and bOptimize vector
for (ipot in 1:NumberOfPotentials) {
  TmpDoc0 <- xmlDoc(src0[[ipot]])
  bOptimize <- append(bOptimize,xpathApply(TmpDoc0,"//bOptimize",xmlValue)[[1]])
  if (bOptimize[ipot] == "true") {
   ipotindex <- append(ipotindex,ipot)
   filename <- paste('INPUT/Param',ipot-1,'.dat',sep="")
   refDistrib <- as.data.frame(read.table(filename))
   if (ipot>degpot)
     refDistrib[,1] <- refDistrib[,1]/pi*180.0
   refDistrib.spl <- spline(refDistrib[,1],refDistrib[,2]/max(refDistrib[,2]),n=2*length(refDistrib[,1]))
   refDistribs <- lappend(refDistribs,refDistrib.spl)
  }
}
distribs.spl<-list()
for (ipot in 1:NumberOfPotentials) {
  distribs.spl <- list()
  losses <- vector()
  acceptedIterations <- vector()
  if (bOptimize[ipot] == "true") {
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
       # And only accepted distributions (always true if IBI is used)
        if (status=="true") {
          filename <- paste('OUTPUT/r',run-1,'/Param',ipot-1,'.dat',sep="")
          distrib <- as.data.frame(read.table(filename))
          if (ipot>degpot)
            distrib[,1] <- distrib[,1]/pi*180.0          
          distrib.spl <- spline(distrib[,1],distrib[,2]/max(distrib[,2]),n=2*length(distrib[,1]))
          distribs.spl <- lappend(distribs.spl,distrib.spl)
          #lineWidth <- 5.0        
          loss <- as.numeric(xpathApply(TmpDoc,"//LossFunction",xmlValue)[[1]])
          losses <- append(losses,loss)   
          acceptedIterations <- append(acceptedIterations,run)          
        }      
    } # iterations loop
    DistAndLosses <- list(dists=distribs.spl,losses=losses,acceptedIterations=acceptedIterations)
    DistribALLSpl <- lappend(DistribALLSpl,DistAndLosses)
  } # bOptimize if
} # potentials loop


# Figures about distributions
colfunc = colorRampPalette(c("white","black"))
for (ipot in 1:length(DistribALLSpl)) {
  #par(fig=c(0.0,1.0,0.0,1.0),mar=c(8,8,1,2),mgp=c(1.0,1.4,0.0),oma=c(0.1,0.1,0.1,0.1),new=FALSE)
  NaccIter <- length(DistribALLSpl[[ipot]]$losses)
  loss <- DistribALLSpl[[ipot]]$losses
  colors <- colfunc(NaccIter)
  ranges <- c()
  if (ipot==1) {
    ranges[1] <- 4.5
    ranges[2] <- 6.5
  } else if (ipot==2)  {
    ranges[1] <- 4.0
    ranges[2] <- 6.5
  } else if (ipot==3)  {
    ranges[1] <- 5.0
    ranges[2] <- 7.5
  } else if (ipot==4)  {
    ranges[1] <- 75.0
    ranges[2] <- 110.0
  } else if (ipot==5)  {
    ranges[1] <- 20.0
    ranges[2] <- 85.0
  }  
  filename <- paste("FigForPaper/Param",ipotindex[ipot]-1,".eps",sep="")
  #tiff(filename, width = 1200, height = 1200, units = 'px', compression = c("none") )
  postscript(filename,height = 5, width=10)
  par(mar=c(7,5,1,1),mgp=c(5,2,0),oma=c(0.0,0.0,0.0,0.0),new=FALSE)
  refDistrib <- refDistribs[[ipot]]
  if (IBI) {
    plot(refDistrib$x,refDistrib$y,type="l",lty=2,col=554,lwd=5.0, xlab="", ylab="", xlim=ranges, ylim=c(0,1), cex.axis=3.4, cex.lab=3.5, frame=TRUE, axes=FALSE) # OnlyIBI
    axis(side = 1, tick = TRUE, font=2, cex.axis=2.5, cex.lab=2.2) 
    axis(side = 2, tick = TRUE, at=c(0.5,1.0), font=2, cex.axis=2.5, cex.lab=2.2)  
  } else {
    plot(refDistrib$x,refDistrib$y,type="l",lty=2,col=554,lwd=5.0, xlab="", ylab="", xlim=ranges, ylim=c(0,1), cex.axis=3.4, cex.lab=3.5, frame=TRUE, axes=FALSE)
    axis(side = 2, tick = TRUE, at=c(0.5,1.0), font=2, cex.axis=2.5, cex.lab=2.2)  
  }
  lineWidth <- 5.0
  for (i in 1:length(DistribALLSpl[[ipot]]$dists)) {
    distrib <- DistribALLSpl[[ipot]]$dists[[i]]
    lineWidth <- lineWidth - 3.0/NaccIter
    lines(distrib$x,distrib$y,col=colors[i],lwd=lineWidth)
  } # End iteration loop for distributions
  lines(refDistrib$x,refDistrib$y,type="l",lty=2,col=554,lwd=5.0)
  
  dev.off()
}

# Figures about correlations (Lossfunction correlations)
# Only for IBI and OnlyIBI 
if (!MC) {
  m<-matrix(0,length(DistribALLSpl),length(DistribALLSpl[[1]]$losses))
  for (ipot in 1:length(DistribALLSpl)) 
    m[ipot,]<- DistribALLSpl[[ipot]]$losses
  
  mt <- t(m)
  mt.scaled<-scale(mt)
  
  
  # Shift to positve value only
  abmin <- min(mt.scaled)
  mt.scaled<- mt.scaled + abs(abmin)
  
  # Shift by 1 each for better visualization
  for (ipot in 1:length(DistribALLSpl)) 
    mt.scaled[,ipot] <- mt.scaled[,ipot] + 3*ipot
  
  minx <- 2
  maxx <- length(DistribALLSpl[[1]]$losses)
  miny <- min(mt.scaled)
  maxy <- max(mt.scaled)
  
  if (ipot==5) {
    names<-c("","","","","")
  } else  {
    names<-c("","","","")
  }
  colnames(mt.scaled) <- names     
  Spearman <- cor(mt.scaled,method="spearman")     
  #Spearman[lower.tri(Spearman,diag=T)] = 0.0 # set 0 diagonal and lower tri matrix
  pmat <- matrix(0,dim(Spearman)[1],dim(Spearman)[2])
  # Add p-values to the lower triangular part of Pearson matrix
  for (i in 1:dim(mt.scaled)[2]) for (j in i:dim(mt.scaled)[2]) {
    pvalue <- cor.test(mt.scaled[,i],mt.scaled[,j],alternative="two.side",method="spearman")$p.val
    pmat[i, j] <- pmat[j, i] <- pvalue
    #if (pvalue<0.01) {  Spearman[j,i] <- 0.01 }
  }
  filename <- "FigForPaper/SpearmanCorrelation.eps"
  postscript(filename) 
  #my.plotcorr(Spearman, col=colors[((Spearman + 1)/2) * 100], main='Spearman correlations', lower.panel='none', diag='none',mar=0.0 + c(1, 0, 2, 2),outline=TRUE)
  #image.plot(legend.only=TRUE,col=colors, zlim=c(-1,1),horizontal=TRUE,legend.mar=6.4,legend.shrink=0.4,legend.width=0.4 )
  #corrplot(Spearman)  
  # Reverse significant with insignificant, just for visualizing reasons
  pmm <- pmat
  pmm[pmat>0.01] <- 0.0
  pmm[pmat<=0.01] <- 1.0
  corrplot(Spearman,p.mat=pmm,sig.level =0.01,method="color",type="lower",diag=F,cl.cex=1.4,pch="*")
  dev.off()  
  caption <- paste("Correlations among loss functions",sep="")
  tab <- xtable((Spearman[,]), caption=caption)
  filename <- "FigForPaper/SpearmanCorrelation.tex"
  print(tab,file=filename,append=F,table.placement = "h", caption.placement="bottom")  
}



# Figures about Loss function
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
  maxy <- maxy + 4.0
  plot(1, type="n", axes=FALSE, frame=TRUE, xlab='', ylab='', xlim=c(minx,maxx), ylim=c(miny,maxy), cex.axis=3.0, cex.lab=3.0) # MC  
}else{
  plot(1, type="n", axes=FALSE, frame=TRUE, xlab='', ylab='',  xlim=c(minx,maxx), ylim=c(miny+3,maxy+5), cex.axis=3.0, cex.lab=3.0) # IBI and OnlyIBI
}
axis(side = 1, tick = TRUE, font=2, cex.axis=2.5, cex.lab=2.2) 
axis(side = 2, tick = FALSE, at=seq(miny,maxy), font=2, cex.axis=2.5, cex.lab=2.2, labels=FALSE)  
mtext("Iteration", side=1, line=4, cex=3.5)
mtext("[a.u.]", side=2, line=1.0, cex=3.5)
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
    ptype<-ipot+20
  } else if (ipot==2) {
    color<-"#d50000"
    ptype<-ipot+20
  } else if (ipot==3) {
    if (length(DistribALLSpl)==5) {
      color<-"#00d254"
      ptype<-ipot+20
    } else {
      color<-"#e56cff"
      ptype<-ipot+21
    }
  } else if (ipot==4)  {
    if (length(DistribALLSpl)==5) {
      color<-"#e56cff"
      ptype<-ipot+20
    } else {
      color<-"#fffd43"
      ptype<-ipot+21
    }
  } else if (ipot==5) {
    color<-"#fffd43"
    ptype<-ipot+20
  }
  points(Acc,loss+1*ipot-1, pch=ptype, bg=color, cex=2.0)
  #symbols(Acc,loss,circles=rep(0.5,length(Acc)),inches=1/8,ann=F,bg=color,add=TRUE)
} # end ipot loop
dev.off()

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
