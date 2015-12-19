lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}
# x: the vector
# n: the number of samples
# centered: if FALSE, then average current sample and previous (n-1) samples
#           if TRUE, then average symmetrically in past and future. (If n is even, use one more sample from future.)
movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

c4 <- brewer.pal(4,"Set1") # color blind safe, bw safe, printer friendly
cf <- brewer.pal(10,"Spectral")

# function to find medoid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}

# Normalize distribution
pnorm <- function(p,type) {
  smoothing <- TRUE
  if (smoothing) {
    p.ks <- locpoly(p[,1],p[,2],bandwidth = 0.05) 
    x<-p.ks$x
    y<-p.ks$y  
  } else {
    x<-p[,1]     
    y<-p[,2]
  }
  x <- ConvertUnit(x,type)
  dx <- diff(x)[1]
  y[is.na(y)]<-0
  Z <- sum(y)*dx
  pn <- matrix(0,nrow=length(x),ncol=2)
  pn[,1] <- x
  pn[,2] <- y/Z
  return(pn)
}
# Convert Units
ConvertUnit <- function(x,type) {
  switch(type,
         r13 = x,
         r14 = x,
         r15 = x,
         theta = rad2deg(x),
         phi = rad2deg(x))
}
# Conversion deg to rad
deg2rad <- function(deg) {
  rad <- deg/180.0*pi
  return(rad)
}

# Conversion rad to deg
rad2deg <- function(rad) {
  deg <- rad*180.0/pi
  return(deg)
}
# Set Labels
SetLabels <- function(type) {
  switch(type,
         r13 = expression ( paste(r["i,i+2"], " [", ring(A), "]" ) ),
         r14 = expression ( paste(r["i,i+3"], " [", ring(A), "]" ) ),
         r15 = expression ( paste(r["i,i+4"], " [", ring(A), "]" ) ),
         theta = expression ( paste(theta," [Deg]",sep=" ") ),
         phi =  expression ( paste(phi," [Deg]",sep=" ") )
  )
}
# General theme for PNG (ok for 600 dpi resolution 6.5x3.25in)
thm2<- theme(panel.background = element_rect(fill = 'white'),
             panel.border = element_rect(colour = "black", fill=NA, size=1),
             axis.ticks.x = element_line(size=0.8, colour="black"),
             axis.title.x =  element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.text.x  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.ticks.y = element_line(size=0.8, colour="black"),
             axis.text.y  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.title.y = element_text(angle=90, vjust=0.5, size=12, face="bold", colour="black"),
             plot.title = element_text(lineheight=3, face="bold", color="black", size=12),
             legend.title  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             legend.text = element_text(lineheight=3, face="bold", color="black", size=12),
             strip.text = element_text(lineheight=3, face="bold", color="black", size=12)) # Text for facets header


library("cluster")
library("bio3d")
library("ggplot2")
library("miscTools")
library("fields")
library("MASS")
library("RColorBrewer")
require(gridExtra)
library("KernSmooth")

plotphipsi <- FALSE
plotdists <- FALSE
plotrmsd <- TRUE
clustering <- FALSE
halpha<-FALSE

# Load Trajectory
his <- read.pdb("../Data/LongSimulation/ana/HISo.pdb",multi=TRUE)
ns_step <- 0.2 # number of ns per frame

if (packageVersion("bio3d")>="2.1.0")
  his$xyz.models<-his$xyz

Nsteps <- nrow(his$xyz.models[])
x <- c(1:Nsteps) * ns_step




if (plotphipsi) {
  
  
   Allphi <- vector()
   Allpsi <- vector()
   AllphiFilename <- ("../Data/LongSimulation/ana/AllAtomTraj/PDBs/Allphi.dat")
   AllpsiFilename <- ("../Data/LongSimulation/ana/AllAtomTraj/PDBs/Allpsi.dat")
   recompute <- FALSE
   if (recompute) {
     for (index in 1:5000) {
       filename=paste("../Data/LongSimulation/ana/AllAtomTraj/PDBs/h",sprintf("%03d",index),".pdb",sep="")
       if (file.exists(filename)) {
        cat(filename,"\n")
        pdbAAi <- read.pdb(filename,multi=TRUE)
        tori <- torsion.pdb(pdbAAi)
        tori$phi[is.na(tori$phi)]<-0
        tori$psi[is.na(tori$psi)]<-0
        Allphi<-append(Allphi,tori$phi)
        Allpsi<-append(Allpsi,tori$psi)
       } 
     }
     save(Allphi,file=AllphiFilename)
     save(Allpsi,file=AllpsiFilename)
   } else {
     load(AllphiFilename)
     load(AllpsiFilename)
   }
   ngridpoints <- 80
   f1<-kde2d(Allphi,Allpsi,lims=c(-100,0,-100,0),n=ngridpoints)
   PhiPsiFrame <- data.frame()
   PhiPsiFrame<- expand.grid(x=seq(-100,0,,ngridpoints),y=seq(-100,0,,ngridpoints))
   PhiPsiV<-as.vector(as.matrix(f1$z))
   PhiPsiFrame$z <- PhiPsiV/(sum(PhiPsiV)*max(diff(PhiPsiFrame$x))*max(diff(PhiPsiFrame$y)))
   ggpltBASE <- ggplot(data=PhiPsiFrame) +  thm2    
   
   
#    runindexes <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,43,45,46,47,48,49,50,54,62,67,68,73,76,92,94,101,109,134,140,159,160,291,423,428,442)
#    frameindexes <- seq(1,1571,by=50)
#    runs <- rep(runindexes,each=(length(frameindexes)*20)) # 20 is the number of torsions, i.e. length(Allphi)
#    eachrun<-rep(frameindexes,each=20)
#    frames<-rep(eachrun,length(runindexes))
#    Allphi <- vector()
#    Allpsi <- vector()
#    for (run in runindexes) {      
#     for (frame in frameindexes) {
#      fik <- paste("../../MCSA-PhiPsiAnalysis/AA/run_",run,"/h",sprintf("%03d",frame),".pdb",sep="")
#      pdbAAi <- read.pdb(fik,multi=TRUE)
#      tori <- torsion.pdb(pdbAAi)
#      #tori$phi[is.na(tori$phi)] <- -50
#      #tori$psi[is.na(tori$psi)] <- -50
#      Allphi <- append(Allphi,tori$phi)
#      Allpsi <- append(Allpsi,tori$psi)  
#     }
#    }
#    DF<-data.frame(runs=runs,frames=frames,xvar=Allphi,yvar=Allpsi)
#    cbPalette <- c("#D55E00", "#56B4E9", "#F0E442")
#    ggplt <- ggpltBASE + geom_point(data=DF,aes(x=xvar,y=yvar,colour=runs),na.rm=TRUE) + 
#      scale_colour_gradientn(colours = cbPalette)
    
   
   ggplts <- ggpltBASE +
             stat_contour(aes(x=x,y=y,z=z,fill=..level..),
                          geom="polygon",
                          size=1.0,
                          breaks = c(0.0001, 0.00015, 0.0002, 0.00025, 0.0003, 0.0005, 0.0006)) +
             scale_fill_gradientn(name="Average\nLoss", colours= cf) +
             xlab(expression(paste(phi," (Deg)",sep=""))) +
             ylab(expression(paste(psi," (Deg)",sep=""))) +     
             scale_x_continuous(expand=c(0.01,0), limits=c(-180,0)) +
             scale_y_continuous(expand=c(0.01,0.01), limits=c(-180,0)) +   
             theme(legend.position="none")
   
   filename <- "pdf/phipsi-LongSimulation.pdf"
   ggsave(file=filename,plot=ggplts,device=pdf,width=3.5, height=3.25,units="in")  
}


# RMSD initial structure
if (plotrmsd) {
  rmsd1 <- rmsd(his$xyz.models[3000,],his$xyz.models,fit=TRUE)
  # Moving average
  rmsdvAVG <- movingAverage(rmsd1,20,centered=TRUE)
  rmsdFrame<-data.frame(x=x,y=rmsd1,avg=rmsdvAVG)
  plt <- ggplot(data=rmsdFrame,aes(x,y)) +    
    geom_line(aes(x,avg),col="black",size=0.5) +
    xlab("Timestep (ns)") + 
    ylab("RMSD (\uc5)") + 
    scale_x_continuous(expand=c(0.01,0.01)) +
    scale_y_continuous(expand=c(0.01,0.01),limits=c(2.0,4)) +
    thm2
  filename <- "pdf/rmsd-LongSimulation.pdf"
  ggsave(file=filename,plot=plt,device=pdf,width=3.25, height=3.25,units="in")  
}

# for (index in 1:10:Nsteps) {
#   fitted <- fit.xyz(his$xyz.models[1,],his$xyz.models[index,])  
#   #oriented <- as.vector(orient.pdb(fitted))
#   mo <- matrix(fitted,nrow=20,byrow = TRUE)  
#   if (index==1) 
#     plot(mo[1,2],mo[1,3])    
#   
#   #plot(mo[1,2],mo[1,3],xlim=c(-4,4),ylim=c(-4,4))
#   
#   points(mo[1,2],mo[1,3],col="red")
#   points(mo[2,2],mo[2,3],col="blue")
#   points(mo[3,2],mo[3,3],col="green")
#   points(mo[3,2],mo[3,3],col="yellow")
# }

# Clustering
if (clustering) {
 # RMSD Matrix
 recomputeClusters<-TRUE
 if (recomputeClusters) {
   vRMSD <- vector()
   for (index in 1:Nsteps) {
     vRMSD <- c(vRMSD,rmsd(his$xyz.models[index,],his$xyz.models[c(index:Nsteps),],fit=TRUE)) 
     print(paste("Step",index,sep=" "))
   }
   RMSDMatrix <- symMatrix(vRMSD,nrow = Nsteps)
   # Clustering
   RMSDDist <- as.dist(RMSDMatrix)
 
   filename <- "RMSDDist.Rdata"
   save(RMSDDist,file=filename)
 } else {
   load("RMSDDist.Rdata")
 }
 
 if (recomputeClusters) {
  # H clust
   fit <- hclust(RMSDDist,method="ward")
   groups <- cutree(fit, k=5)
   #plot(groups)
 
   # Partitioning Around Medoids
   asw <- numeric()
   for (k in 3:10)
     asw[[k]] <- pam(RMSDDist, k) $ silinfo $ avg.width
   k.best <- which.max(asw)
   cat("silhouette-optimal number of clusters:", k.best, "\n")
   k <- k.best
   fitpam<-pam(RMSDDist,k=k)
 
   # Write out medoids and full clusters
   lista<-c(1:length(his$xyz.models[,1]))
   for (i in 1:length(fitpam$medoids)) {
     index <- fitpam$medoids[i]
     filename <- paste("Medoid_",i,"_frame_",index,".pdb",sep="")
     write.pdb(file=filename,xyz=his$xyz.models[index,])
     clusterindexes <- lista[fitpam$clustering==i]
     filename <- paste("Cluster_",i,".pdb",sep="")
     write.pdb(file=filename,xyz=his$xyz.models[clusterindexes,])  
   }
 }
}


# Distributions
if (plotdists) {
  dists <- list()
  AllDists <- list()
  # Load Distributions 
  dists$ang<-pnorm(read.table("../Data/LongSimulation/ana/Distributions/Param3.dat"),"theta")
  dists$dih<-pnorm(read.table("../Data/LongSimulation/ana/Distributions/Param4.dat"),"phi")
  dists$HB13<-pnorm(read.table("../Data/LongSimulation/ana/Distributions/Param1.dat"),"r13")
  dists$HB14<-pnorm(read.table("../Data/LongSimulation/ana/Distributions/Param2.dat"),"r14")
  if(halpha) {
   dists$HB15<-pnorm(read.table("../Data/LongSimulation/ana/Distributions/Param15.dat"),"r15")
  }
  AllDists<-lappend(AllDists,dists)
  # REferencedists
  dists$ang<-pnorm(read.table("RefDists/Param3.dat"),"theta")
  dists$dih<-pnorm(read.table("RefDists/Param4.dat"),"phi")
  dists$HB13<-pnorm(read.table("RefDists/Param1.dat"),"r13")
  dists$HB14<-pnorm(read.table("RefDists/Param2.dat"),"r14")
  if (halpha) {
   dists$HB15<-pnorm(read.table("RefDists/Param15.dat"),"r15")
  }
  AllDists<-lappend(AllDists,dists)
    
  d<-as.data.frame(AllDists)
  
  ggpltBASE <- ggplot(data=d) + 
    scale_y_continuous(breaks=c(0.0,0.5,1.0)) + 
    thm2 + theme(axis.text.y  = element_blank(),
                 axis.ticks.y = element_blank())
  
  xlab <- SetLabels("theta")
  pang <- ggpltBASE + 
    geom_line(aes(x=ang.1.1, y = ang.2.1),col="black",size=1) + 
    geom_line(aes(x=ang.1, y = ang.2),col=c4[1],size=1) +     
    xlim(75,110) +
    xlab(xlab) + 
    theme(plot.margin = unit(c(0.5,0.4,0,0), "line"),
          axis.text.y  = element_blank(),
          axis.title.y  = element_blank())
  
  xlab <- SetLabels("phi")
  pdih <- ggpltBASE +
    geom_line(aes(x=dih.1.1, y = dih.2.1), col="black",size=1) + 
    geom_line(aes(x=dih.1, y = dih.2), col=c4[1],size=1,) +     
    xlim(30,120) +
    xlab(xlab) + 
    theme(plot.margin = unit(c(0.5,0.4,0,0), "line"),
          axis.text.y  = element_blank(),
          axis.title.y  = element_blank())
  
  xlab <- SetLabels("r13")
  pHB13 <- ggpltBASE + 
    geom_line(aes(x=HB13.1.1, y = HB13.2.1),col="black",size=1) +
    geom_line(aes(x=HB13.1, y = HB13.2),col=c4[1],size=1) +    
    xlim(4,7) +  xlab(xlab)  +
    theme(plot.margin = unit(c(0.0,0.4,0,0), "line"),
          axis.text.y  = element_blank(),
          axis.title.y  = element_blank())
  
  xlab <- SetLabels("r14")
  pHB14 <- ggpltBASE + 
    geom_line(aes(x=HB14.1.1, y = HB14.2.1),col="black",size=1) +
    geom_line(aes(x=HB14.1, y = HB14.2),col=c4[1],size=1) +    
    xlim(4,8.5) +  xlab(xlab) +
    theme(plot.margin = unit(c(0.0,0.4,0,0), "line"),
          axis.text.y  = element_blank(),
          axis.title.y  = element_blank())
#  
#   if (halpha) { 
#   xlab <- SetLabels("r15")
#   pHB15 <- ggpltBASE + 
#     geom_line(aes(x=HB15.1, y = HB15.2),col=c4[1],size=1) +
#     geom_line(aes(x=HB15.1, y = HB15.2.1),col="black",size=1) +
#     xlim(5,8) +  xlab(xlab) +   ylab("")
#   }

  filename <- "pdf/Distributions-LongSimulation.pdf"
  pdf(filename,width=3.25,height=3.25)
  if (halpha) {
   grid.arrange(pang,pdih,pHB13,pHB14,pHB15)
  } else {
   grid.arrange(pang,pdih,pHB13,pHB14)
  }
  dev.off()
}
