#library("fields") # required for plot colorscale

oI<-read.table("type.dat") #oI = 0 : MC, oI=1 : MC+IBI, oI= 2 : OnlyIBI, 

library("fields")
library("ellipse")
library("xtable")

sf <- paste(Sys.getenv("Rscripts"),"SpecialPlotFunctions.R",sep="")
source(sf)

load("DistribALLSpl")

colfunc = colorRampPalette(c("black","white"))

names <- vector()
# For each potential
   # mpg (labels,tics-mark-label,tics-mark)
   # mar (bottom,left,top,right) in lines
   # oma (bottom,left,top,right) in lines
   #par(mfrow=c(3,2),mar=c(4,4,4,2),mgp=c(2.5,1.2,0.0),oma=c(0.1,0.1,0.1,0.1))
   letter <- c("a","b","c","d","e")
for (poti in 1:length(DistribALLSpl)) {

   #par(fig=c(0.0,1.0,0.0,1.0),mar=c(8,8,1,2),mgp=c(1.0,1.4,0.0),oma=c(0.1,0.1,0.1,0.1),new=FALSE)
   NaccIter <- length(DistribALLSpl[[poti]]$acceptedIterations)
   loss <- DistribALLSpl[[poti]]$lossAcc
   colors <- colfunc(NaccIter)


   name <- DistribALLSpl[[poti]]$name
   names[poti] <- name
 
   cat(name,"\n")   
   filename <- paste("FigForPaper/",name,".eps",sep="")
   #tiff(filename, width = 1200, height = 1200, units = 'px', compression = c("none") )
   postscript(filename,height = 5, width=10)
   par(mar=c(7,5,1,1),mgp=c(5,2,0),oma=c(0.0,0.0,0.0,0.0),new=FALSE)
   
   # refernce distrib
   referenceDistrib <- DistribALLSpl[[poti]]$referenceDistrib
   xLab <- paste("r [",DistribALLSpl[[poti]]$qunit,"]",sep="")
   ranges <- DistribALLSpl[[poti]]$ranges
   let <- paste("(",letter[poti],")",sep="")
   if (name=="r13_Loss") {
     ranges[1] <- 4.5
     ranges[2] <- 6.5
   } else if (name=="r14_Loss")  {
     ranges[1] <- 4.0
     ranges[2] <- 6.5
   } else if (name=="r15_Loss")  {
     ranges[1] <- 5.0
     ranges[2] <- 7.5
   } else if (name=="angle_Loss")  {
     ranges[1] <- 75.0
     ranges[2] <- 110.0
   } else if (name=="dihedral_Loss")  {
     ranges[1] <- 20.0
     ranges[2] <- 85.0
   }
    if (oI[1]==0||oI[1]==2) {
   plot(referenceDistrib$x,referenceDistrib$y,type="l",lty=2,col=554,lwd=5.0, xlab="", ylab="", xlim=ranges, ylim=c(0,1), cex.axis=3.4, cex.lab=3.5, frame=TRUE, axes=FALSE)
    axis(side = 2, tick = TRUE, at=c(0.5,1.0), font=2, cex.axis=2.5, cex.lab=2.2)  
} else {
  plot(referenceDistrib$x,referenceDistrib$y,type="l",lty=2,col=554,lwd=5.0, xlab=xLab, ylab="", xlim=ranges, ylim=c(0,1), cex.axis=3.4, cex.lab=3.5, frame=TRUE, axes=FALSE) # OnlyIBI
   axis(side = 1, tick = TRUE, font=2, cex.axis=2.5, cex.lab=2.2) 
   axis(side = 2, tick = TRUE, at=c(0.5,1.0), font=2, cex.axis=2.5, cex.lab=2.2)  
 }

  # box("figure",lty="dashed", col="blue")
  # box("inner", lty="dotted", col="green")
   # For each accepted distributions
   lineWidth <- 5.0
   for (i in 1:length(DistribALLSpl[[poti]]$distributions)) {
         distrib <- DistribALLSpl[[poti]]$distributions[[i]]
         lineWidth <- lineWidth - 3.0/NaccIter
         lines(distrib$x,distrib$y,col=colors[i],lwd=lineWidth)
   } # End iteration loop for distributions
   
}

   dev.off()
   #dev.new()
   


   library("corrplot")
   # Only for IBI and OnlyIBI 
   if (oI[1]!=0) {
     m<-matrix(0,length(DistribALLSpl),length(DistribALLSpl[[1]]$lossAcc))
     for (poti in 1:length(DistribALLSpl)) 
       m[poti,]<- DistribALLSpl[[poti]]$lossAcc

       mt <- t(m)
       mt.scaled<-scale(mt)
     
     # Shift to positve value only
     abmin <- min(mt.scaled)
     mt.scaled<- mt.scaled + abs(abmin)
     
     # Shift by 1 each for better visualization
     for (poti in 1:length(DistribALLSpl)) 
      mt.scaled[,poti] <- mt.scaled[,poti] + 3*poti

     minx <- 2
     maxx <- length(DistribALLSpl[[1]]$lossAcc)
     miny <- min(mt.scaled)
     maxy <- max(mt.scaled)
     
     if (length(names)==5) {
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


filename <- paste("FigForPaper/LossFunction.eps")
postscript(filename, height=5, width=10)

par(mar=c(7,8,2,1),mgp=c(4,1.5,0),oma=c(0.0,0.0,0.0,0.0),new=FALSE)

   # Find x-y ranges Form MC only
   if (oI[1]==0) {
    minx<-1
    maxx<-0
    miny<-1000.0
    maxy<-0.0
    for (poti in 1:length(DistribALLSpl)) {
     loss <- DistribALLSpl[[poti]]$lossAcc
     Acc <- DistribALLSpl[[poti]]$acceptedIterations
    if (max(Acc,na.rm=TRUE)>maxx)
	maxx <- max(Acc,na.rm=TRUE)

    if (max(loss,na.rm=TRUE)>maxy)
	maxy <- max(loss,na.rm=TRUE) 

    if (min(loss,na.rm=TRUE)<miny)
	miny <- min(loss,na.rm=TRUE)

    }
   }
   

   # Do the first plot
   if (oI!=0) { 
    plot(1, type="n", axes=FALSE, frame=TRUE, xlab='', ylab='',  xlim=c(minx,maxx), ylim=c(miny+3,maxy+5), cex.axis=3.0, cex.lab=3.0) # IBI and OnlyIBI
   }else{
    maxy <- maxy + 4.0
    plot(1, type="n", axes=FALSE, frame=TRUE, xlab='', ylab='', xlim=c(minx,maxx), ylim=c(miny,maxy), cex.axis=3.0, cex.lab=3.0) # MC
   }
  
   axis(side = 1, tick = TRUE, font=2, cex.axis=2.5, cex.lab=2.2) 
   #axis(side = 2, tick = TRUE, at=seq(miny+0.1,maxy), font=2, cex.axis=2.5, cex.lab=2.2, labels=FALSE)  
   axis(side = 2, tick = FALSE, at=seq(miny,maxy), font=2, cex.axis=2.5, cex.lab=2.2, labels=FALSE)  
   mtext("Iteration", side=1, line=4, cex=3.5)
   mtext("[a.u.]", side=2, line=1.0, cex=3.5)


# Add overalys
for (poti in 1:length(DistribALLSpl)) {
  
  if (oI[1]!=0) {
    loss <- mt.scaled[,poti]+abs(miny) # IBI and OnlyIBI
  }
  else
    loss <- DistribALLSpl[[poti]]$lossAcc # MC
  
  Acc <- DistribALLSpl[[poti]]$acceptedIterations
  
  # Set colors and point type
  if (poti==1) {
    color<-"steelblue2"
    ptype<-poti+20
  } else if (poti==2) {
    color<-"#d50000"
    ptype<-poti+20
  } else if (poti==3) {
    if (length(DistribALLSpl)==5) {
      color<-"#00d254"
      ptype<-poti+20
    } else {
      color<-"#e56cff"
      ptype<-poti+21
    }
  } else if (poti==4)  {
    if (length(DistribALLSpl)==5) {
      color<-"#e56cff"
      ptype<-poti+20
    } else {
      color<-"#fffd43"
      ptype<-poti+21
    }
  } else if (poti==5) {
    color<-"#fffd43"
    ptype<-poti+20
  }
  
  points(Acc,loss+1*poti-1, pch=ptype, bg=color, cex=2.0)
  #symbols(Acc,loss,circles=rep(0.5,length(Acc)),inches=1/8,ann=F,bg=color,add=TRUE)
} # end poti loop
   dev.off()




