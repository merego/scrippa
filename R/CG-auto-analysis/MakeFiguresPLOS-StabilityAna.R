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


# function to find medoid in cluster i
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}

library("cluster")
library("bio3d")
library("ggplot2")
library("miscTools")
library("fields")
library("MASS")
require(gridExtra)

plotphipsi <- TRUE
plotdists <- FALSE
plotrmsd <- FALSE

# Load Trajectory
his <- read.pdb("HISo.pdb",multi=TRUE)
ns_step <- 1 # number of ns per frame

if (packageVersion("bio3d")=="2.1.0")
  his$xyz.models<-his$xyz

Nsteps <- nrow(his$xyz.models[])
x <- c(1:Nsteps) * ns_step




if (plotphipsi) {
  
  
   Allphi <- vector()
   Allpsi <- vector()
   AllphiFilename <- ("AllAtomTraj/PDBs/Allphi.dat")
   AllpsiFilename <- ("AllAtomTraj/PDBs/Allpsi.dat")
   reload <- FALSE
   if (reload) {
     for (index in 1:2000) {
       filename=paste("AllAtomTraj/PDBs/h",sprintf("%03d",index),".pdb",sep="")
       cat(filename,"\n")
       pdbAAi <- read.pdb(filename,multi=TRUE)
       tori <- torsion.pdb(pdbAAi)
       tori$phi[is.na(tori$phi)]<-0
       tori$psi[is.na(tori$psi)]<-0
       Allphi<-append(Allphi,tori$phi)
       Allpsi<-append(Allpsi,tori$psi)
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
   PhiPsiFrame$z <- PhiPsiV/max(PhiPsiV)
   ggpltBASE <- ggplot(data=PhiPsiFrame) +         
     theme(panel.background = element_rect(fill = 'white',colour="black",size=1.0), 
           panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           axis.title.x = element_text(face="bold", colour="black", size=30), 
           axis.ticks.x = element_line(size=1.2, colour="black"), 
           axis.text.x  = element_text(angle=0, vjust=0.5, size=25, face="bold", colour="black"), 
           axis.title.y = element_text(face="bold", colour="black", size=30, angle=90), 
           axis.ticks.y = element_line(size=1.2, colour="black"), 
           axis.text.y  = element_text(angle=0, vjust=0.5, size=25, face="bold", colour="black"), 
           plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
           legend.title  = element_blank(), 
           legend.text = element_text(lineheight=3, face="bold", color="black", size=25),
           legend.key.height = unit(2,"cm"))
    
   
   
   runindexes <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,43,45,46,47,48,49,50,54,62,67,68,73,76,92,94,101,109,134,140,159,160,291,423,428,442)
   frameindexes <- seq(1,1571,by=50)
   runs <- rep(runindexes,each=(length(frameindexes)*20)) # 20 is the number of torsions, i.e. length(Allphi)
   eachrun<-rep(frameindexes,each=20)
   frames<-rep(eachrun,length(runindexes))
   Allphi <- vector()
   Allpsi <- vector()
   for (run in runindexes) {      
    for (frame in frameindexes) {
     fik <- paste("../../MCSA-PhiPsiAnalysis/AA/run_",run,"/h",sprintf("%03d",frame),".pdb",sep="")
     pdbAAi <- read.pdb(fik,multi=TRUE)
     tori <- torsion.pdb(pdbAAi)
     #tori$phi[is.na(tori$phi)] <- -50
     #tori$psi[is.na(tori$psi)] <- -50
     Allphi <- append(Allphi,tori$phi)
     Allpsi <- append(Allpsi,tori$psi)  
    }
   }
   DF<-data.frame(runs=runs,frames=frames,xvar=Allphi,yvar=Allpsi)
   cbPalette <- c("#D55E00", "#56B4E9", "#F0E442")
   ggplt <- ggpltBASE + geom_point(data=DF,aes(x=xvar,y=yvar,colour=runs),na.rm=TRUE) + 
     scale_colour_gradientn(colours = cbPalette)
    
   
   ggplts <- ggplt + 
     xlab(expression(phi)) +
     ylab(expression(psi)) +
     scale_x_continuous(expand=c(0.01,0), limits=c(-180,0)) +
     scale_y_continuous(expand=c(0.01,0.01), limits=c(-180,0)) +
     stat_contour(aes(x=x,y=y,z=z,fill=..level..),
                  size=1.0,
                  breaks = c(0.05, 0.1, 0.2, 0.4, 0.6, 1.0)) +
     scale_fill_continuous(guide = 'none')
   
   print(ggplts)
   
   filename <- "imgs/phipsi.pdf"
   pdf(filename,width=8.0,height=6.0)
   print(ggplts)
   dev.off()
}


# RMSD initial structure
if (plotrmsd) {
  rmsd1 <- rmsd(his$xyz.models[1,],his$xyz.models,fit=TRUE)
  # Moving average
  rmsdvAVG <- movingAverage(rmsd1,20,centered=TRUE)
  rmsdFrame<-data.frame(x=x,y=rmsd1,avg=rmsdvAVG)
  ggplt <- ggplot(data=rmsdFrame,aes(x,y)) +
    geom_line(col="gray") +
    geom_line(aes(x,avg),col="black",size=2) +
    xlab("Timestep (ns)") + 
    ylab("RMSD (\uc5)") +
    theme(panel.background = element_rect(fill = 'white',colour="black",size=1.0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_text(face="bold", colour="black", size=30), 
          axis.ticks.x = element_line(size=1.2, colour="black"), 
          axis.text.x  = element_text(angle=0, vjust=0.5, size=25, face="bold", colour="black"), 
          axis.title.y = element_text(face="bold", colour="black", size=30, angle=90), 
          axis.ticks.y = element_line(size=1.2, colour="black"), 
          axis.text.y  = element_text(angle=0, vjust=0.5, size=25, face="bold", colour="black"), 
          plot.title = element_text(lineheight=3, face="bold", color="black", size=30), 
          legend.title  = element_blank(), 
          legend.text = element_text(lineheight=3, face="bold", color="black", size=30))
  
  filename <- "imgs/rmsd.pdf"
  pdf(filename,width=8.0,height=6.0)
  print(ggplt)
  dev.off()
  
  hist(rmsd1, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
  lines(density(rmsd1), col="gray", lwd=3)
  
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

# # RMSD Matrix
# recomputeClusters=FALSE
# if (recomputeClusters) {
#   vRMSD <- vector()
#   for (index in 1:Nsteps) {
#     vRMSD <- c(vRMSD,rmsd(his$xyz.models[index,],his$xyz.models[c(index:Nsteps),],fit=TRUE)) 
#     print(paste("Step",index,sep=" "))
#   }
#   RMSDMatrix <- symMatrix(vRMSD,nrow = Nsteps)
#   # Clustering
#   RMSDDist <- as.dist(RMSDMatrix)
# 
#   filename <- "RMSDDist.Rdata"
#   save(RMSDDist,file=filename)
# } else {
#   load("RMSDDist.Rdata")
# }
# 
# if (recomputeClusters) {
#  # H clust
#   fit <- hclust(RMSDDist,method="ward")
#   groups <- cutree(fit, k=5)
#   #plot(groups)
# 
#   # Partitioning Around Medoids
#   asw <- numeric()
#   for (k in 3:10)
#     asw[[k]] <- pam(RMSDDist, k) $ silinfo $ avg.width
#   k.best <- which.max(asw)
#   cat("silhouette-optimal number of clusters:", k.best, "\n")
#   k <- k.best
#   fitpam<-pam(RMSDDist,k=k)
# 
#   # Write out medoids and full clusters
#   lista<-c(1:length(his$xyz.models[,1]))
#   for (i in 1:length(fitpam$medoids)) {
#     index <- fitpam$medoids[i]
#     filename <- paste("Medoid_",i,"_frame_",index,".pdb",sep="")
#     write.pdb(file=filename,xyz=his$xyz.models[index,])
#     clusterindexes <- lista[fitpam$clustering==i]
#     filename <- paste("Cluster_",i,".pdb",sep="")
#     write.pdb(file=filename,xyz=his$xyz.models[clusterindexes,])  
#   }
# }
# dev.off()
# 
# 
# Plot phi-theta correlations
#if (plotphitheta) {
#  PhiTheta<-read.table("Phi-Psi_Map/Results_20140924_105055/Statistics/PDB/Phi_Teta(+).csv")
#  PhiThetaList<-list()
#  PhiThetaList$x <- seq(-180,180,,nrow(PhiTheta)-1)
#  PhiThetaList$y <- seq(-180,180,,ncol(PhiTheta)-1)
#  PhiThetaList$z <- as.matrix(PhiTheta)
#  grid.list<- list( x= seq( -180,180,,800), y=  seq( -180,180,,800))
#  #isurf<-interp.surface.grid(PhiThetaList,grid.list)
#  #image(isurf,xlim=c(70,110),ylim=c(20,80),xlab="Theta",ylab="Phi")
#  
#  PhiThetaFrame<-data.frame()
#  PhiThetaFrame<-expand.grid(x=seq(-180,180,,nrow(PhiTheta)),y=seq(-180,180,,ncol(PhiTheta)))
#  PhiThetaV <- as.vector(as.matrix(PhiTheta))
#  PhiThetaFrame$z<-PhiThetaV/max(PhiThetaV)
#  
#  ggpltBASE <- ggplot(data=PhiThetaFrame) + 
#    theme(panel.background = element_rect(fill = 'white'), 
#          axis.title.x = element_text(face="bold", colour="black", size=20), 
#          axis.ticks.x = element_line(size=1.2, colour="black"), 
#          axis.text.x  = element_text(angle=0, vjust=0.5, size=15, face="bold", colour="black"), 
#          axis.title.y = element_text(face="bold", colour="black", size=20, angle=90), 
#          axis.ticks.y = element_line(size=1.2, colour="black"), 
#          axis.text.y  = element_text(angle=0, vjust=0.5, size=15, face="bold", colour="black"), 
#          plot.title = element_text(lineheight=3, face="bold", color="black", size=30), 
#          legend.title  = element_blank(), 
#          legend.text = element_text(lineheight=3, face="bold", color="black", size=15))
#  
#  ggplt <- ggpltBASE + 
#    geom_tile(aes(x=x,y=y,fill=z)) + 
#    stat_contour(aes(x=x,y=y,z=z),size=0.5) + 
#    scale_fill_gradient(low="green", high="red") +  xlim(-180,180) + ylim(-180,180)
#  print(ggplt)
#}

# Distributions
if (plotdists) {
  dists <- list()
  AllDists <- list()
  # Load Distributions 
  dists$ang<-read.table("Distributions/dists/Statistics/PDB/Bond_angle_CA-CA-CA.dat")
  dists$dih<-read.table("Distributions/dists/Statistics/PDB/Dihedral_angle_CA-CA-CA-CA.dat")
  dists$HB13<-read.table("Distributions/dists/Statistics/PDB/Atom_Distance_2_CA-CA.dat")
  dists$HB14<-read.table("Distributions/dists/Statistics/PDB/Atom_Distance_3_CA-CA.dat")
  dists$HB15<-read.table("Distributions/dists/Statistics/PDB/Atom_Distance_4_CA-CA.dat")
  AllDists<-lappend(AllDists,dists)
  # REferencedists
  dists$ang<-read.table("ReferenceDists/bondang_NPDB.dat")
  dists$dih<-read.table("ReferenceDists/dihang_NPDB.dat")
  dists$HB13<-read.table("ReferenceDists/r13_NPDB.dat")
  dists$HB14<-read.table("ReferenceDists/r14_NPDB.dat")
  dists$HB15<-read.table("ReferenceDists/r15_NPDB.dat")
  AllDists<-lappend(AllDists,dists)
  
  NClusters <- length(AllDists)
  
  d<-as.data.frame(AllDists)
  
  ggpltBASE <- ggplot(data=d) + 
    scale_y_continuous(breaks=c(0.0,0.5,1.0)) +
    theme(panel.background = element_rect(fill = 'white'), 
          axis.title.x = element_text(face="bold", colour="black", size=30), 
          axis.ticks.x = element_line(size=1.2, colour="black"), 
          axis.text.x  = element_text(angle=0, vjust=0.5, size=20, face="bold", colour="black"), 
          axis.title.y = element_text(face="bold", colour="black", size=30, angle=90), 
          axis.ticks.y = element_line(size=1.2, colour="black"), 
          axis.text.y  = element_text(angle=0, vjust=0.5, size=20, face="bold", colour="black"), 
          plot.title = element_text(lineheight=3, face="bold", color="black", size=30), 
          legend.title  = element_blank(), 
          legend.text = element_text(lineheight=3, face="bold", color="black", size=30))
  
  pang <- ggpltBASE + 
    geom_line(aes(x=ang.V1, y = ang.V2/max(ang.V2)),col="black",size=2) + 
    geom_line(aes(x=ang.V1, y = ang.V2.1/max(ang.V2.1)),col="red",size=2,linetype="dashed") + 
    xlim(75,110) +
    xlab(expression(theta)) + 
    ylab("")  
  
  pdih <- ggpltBASE +
    geom_line(aes(x=dih.V1, y = dih.V2/max(dih.V2)), col="black",size=2,) + 
    geom_line(aes(x=dih.V1, y = dih.V2.1/max(dih.V2.1)), col="red",size=2,linetype="dashed") + 
    xlim(30,80) +
    xlab(expression(phi)) + 
    ylab("")
  
  pHB13 <- ggpltBASE + 
    geom_line(aes(x=HB13.V1, y = HB13.V2/max(HB13.V2)),col="black",size=2) +
    geom_line(aes(x=HB13.V1, y = HB13.V2.1/max(HB13.V2.1)),col="red",size=2,linetype="dashed") +
    xlim(4,7) +  xlab(expression ( r["i,i+1"] )) +  ylab("")
  
  pHB14 <- ggpltBASE + 
    geom_line(aes(x=HB14.V1, y = HB14.V2/max(HB14.V2)),col="black",size=2) +
    geom_line(aes(x=HB14.V1, y = HB14.V2.1/max(HB14.V2.1)),col="red",size=2,linetype="dashed") +
    xlim(4,6.5) +  xlab(expression ( r["i,i+2"] )) +   ylab("")
  
  pHB15 <- ggpltBASE + 
    geom_line(aes(x=HB15.V1, y = HB15.V2/max(HB15.V2)),col="black",size=2) +
    geom_line(aes(x=HB15.V1, y = HB15.V2.1/max(HB15.V2.1)),col="red",size=2,linetype="dashed") +
    xlim(5,8) +  xlab(expression ( r["i,i+3"] )) +   ylab("")
  
  filename <- "imgs/Distributions.pdf"
  pdf(filename,width=8.0,height=6.0)
  grid.arrange(pang,pdih,pHB13,pHB14,pHB15)
  dev.off()
}
