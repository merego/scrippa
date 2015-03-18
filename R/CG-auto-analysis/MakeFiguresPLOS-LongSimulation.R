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

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbbPaletteL <- c("#b8b8b8", "#e6bc5d", "#a7d1e9", "#6b9e90", "#f0ea95", "#5691b2", "#d59664", "#cca7bc")

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
plotdists <- TRUE
plotrmsd <- TRUE


BaseDirAlpha="../../Alpha/AlphaWithOptiParam-LongSimulation/ana/"
BaseDir310="../../310/310WithOptiParam-LongSimulation/ana/"
# Load Trajectory
filename<-paste(BaseDirAlpha,"HISo.pdb",sep="")
hisA <- read.pdb(filename,multi=TRUE)
filename<-paste(BaseDir310,"HISo.pdb",sep="")
his3 <- read.pdb(filename,multi=TRUE)
ns_step <- 1 # number of ns per frame (same for alpha and 310)

if (packageVersion("bio3d")>="2.1.0") {
  hisA$xyz.models<-hisA$xyz
  his3$xyz.models<-his3$xyz
}

Nsteps <- nrow(hisA$xyz.models[]) # (same for alpha and 310)
x <- c(1:Nsteps) * ns_step  # (same for alpha and 310)



if (plotphipsi) {
  
  recompute3 <- FALSE
  recomputeA <- FALSE
  
   AllphiA <- vector()
   AllpsiA <- vector()
   Allphi3 <- vector()
   Allpsi3 <- vector()
   AllphiFilenameA <- (paste(BaseDirAlpha,"AlphaAllphiMD.dat",sep=""))
   AllpsiFilenameA <- (paste(BaseDirAlpha,"AlphaAllpsiMD.dat",sep=""))
   AllphiFilename3 <- (paste(BaseDir310,"310AllphiMD.dat",sep=""))
   AllpsiFilename3 <- (paste(BaseDir310,"310AllpsiMD.dat",sep=""))
   # alpha MD run
   if (recomputeA) {
     for (index in 1:2000) {
       filename=paste(BaseDirAlpha,"AllAtomTraj/PDBs/h",sprintf("%03d",index),".pdb",sep="")
       if (file.exists(filename)) {
        cat(filename,"\n")
        pdbAAi <- read.pdb(filename,multi=TRUE)
        tori <- torsion.pdb(pdbAAi)
        tori$phi[is.na(tori$phi)]<-0
        tori$psi[is.na(tori$psi)]<-0
        AllphiA<-append(AllphiA,tori$phi)
        AllpsiA<-append(AllpsiA,tori$psi)
       } 
     }
     save(AllphiA,file=AllphiFilenameA)
     save(AllpsiA,file=AllpsiFilenameA)
   } else {
     load(AllphiFilenameA)
     load(AllpsiFilenameA)
   }
   # 310 MD run   
   if (recompute3) {
     for (index in 1:2000) {
       filename=paste(BaseDir310,"AllAtomTraj/PDBs/h",sprintf("%03d",index),".pdb",sep="")
       if (file.exists(filename)) {
         cat(filename,"\n")
         pdbAAi <- read.pdb(filename,multi=TRUE)
         tori <- torsion.pdb(pdbAAi)
         tori$phi[is.na(tori$phi)]<-0
         tori$psi[is.na(tori$psi)]<-0
         #theta <- angle.xyz(his3$xyz.models[index,])
         #psi <- torsion.xyz(his3$xyz.models[index,])
         Allphi3<-append(Allphi3,tori$phi)
         Allpsi3<-append(Allpsi3,tori$psi)
       } 
     }
     save(Allphi3,file=AllphiFilename3)
     save(Allpsi3,file=AllpsiFilename3)
   } else {
     load(AllphiFilename3)
     load(AllpsiFilename3)
   } 
   
   ngridpoints <- 80
   f1A<-kde2d(AllphiA,AllpsiA,lims=c(-100,0,-100,0),n=ngridpoints)
   f13<-kde2d(Allphi3,Allpsi3,lims=c(-100,0,-100,0),n=ngridpoints)
   PhiPsiFrame <- data.frame()
   PhiPsiFrame<- expand.grid(x=seq(-100,0,,ngridpoints),y=seq(-100,0,,ngridpoints))
   PhiPsiV<-as.vector(as.matrix(f1A$z))
   PhiPsiFrame$zA <- PhiPsiV/max(PhiPsiV)
   PhiPsiV<-as.vector(as.matrix(f13$z))
   PhiPsiFrame$z3 <- PhiPsiV/max(PhiPsiV)   
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
    
   
   # MONTE-CARLO RUNS
   AllphiFilenameA_MC <- (paste(BaseDirAlpha,"AlphaAllphiMC.dat",sep=""))
   AllpsiFilenameA_MC <- (paste(BaseDirAlpha,"AlphaAllpsiMC.dat",sep=""))
   AllphiFilename3_MC <- (paste(BaseDir310,"310AllphiMC.dat",sep=""))
   AllpsiFilename3_MC <- (paste(BaseDir310,"310AllpsiMC.dat",sep=""))   
   runindexes <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,35,36,37,38,39,40,41,43,45,46,47,48,49,50,54,62,67,68,73,76,92,94,101,109,134,140,159,160,291,423,428,442)
   frameindexes <- seq(1,1751,by=200)   
   runs <- rep(runindexes,each=(length(frameindexes)*20)) # 20 is the number of torsions for alpha
   eachrun<-rep(frameindexes,each=20)
   frames<-rep(eachrun,length(runindexes))
   AllphiA <- vector()
   AllpsiA <- vector()
   if (recomputeA) {
    for (run in runindexes) {      
     for (frame in frameindexes) {
      fik <- paste(BaseDirAlpha,"../../MCSA-PhiPsiAnalysis/AA/run_",run,"/h",sprintf("%03d",frame),".pdb",sep="")    
      if (file.exists(fik)) {
       pdbAAi <- read.pdb(fik,multi=TRUE)
       tori <- torsion.pdb(pdbAAi)
       length(tori$phi)
       length(tori$psi)
       #tori$phi[is.na(tori$phi)] <- -50
       #tori$psi[is.na(tori$psi)] <- -50
       AllphiA <- append(AllphiA,tori$phi)
       AllpsiA <- append(AllpsiA,tori$psi)  
      }
     }     
    }
    save(AllphiA,file=AllphiFilenameA_MC)
    save(AllpsiA,file=AllpsiFilenameA_MC)
   } else {
     load(AllphiFilenameA_MC)
     load(AllpsiFilenameA_MC)     
   }     
   DFA<-data.frame(runs=runs,frames=frames,xvar=AllphiA,yvar=AllpsiA)   
   # 310
   runindexes <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,36,37,38,39,40,43,44,45,46,47,48,49,55,58,59,61,62,67,68,69,70,75,77,86,88,89,90,107,110,127,130,150,155,177,214,227,302,329,379,429)
   #runindexes <- c(1,2,3)
   frameindexes <- seq(1,1751,by=200)   
   runs <- rep(runindexes,each=(length(frameindexes)*17)) # 17 is the number of torsions for 310
   eachrun<-rep(frameindexes,each=17)
   frames<-rep(eachrun,length(runindexes))
   Allphi3 <- vector()
   Allpsi3 <- vector()
   if (recompute3) {
    for (run in runindexes) {      
     for (frame in frameindexes) {
       fik <- paste(BaseDir310,"../../MCSA-PhiPsiAnalysis/AA/run_",run,"/h",sprintf("%03d",frame),".pdb",sep="")
       if (file.exists(fik)) {
         pdbAAi <- read.pdb(fik,multi=TRUE)
         tori <- torsion.pdb(pdbAAi)
         #tori$phi[is.na(tori$phi)] <- -50
         #tori$psi[is.na(tori$psi)] <- -50
         Allphi3 <- append(Allphi3,tori$phi)
         Allpsi3 <- append(Allpsi3,tori$psi)  
       } else {
         cat(fik,"\n")
       }       
     }     
    }
    save(Allphi3,file=AllphiFilename3_MC)
    save(Allpsi3,file=AllpsiFilename3_MC)
    } else {
      load(AllphiFilename3_MC)
      load(AllpsiFilename3_MC)
    }   
   
   DF3<-data.frame(runs=runs,frames=frames,xvar=Allphi3,yvar=Allpsi3)   
   #cbPalette1 <- c("#D55E00", "#56B4E9", "#F0E442")
   #cbPalette2 <- c("#D55E00", "#655555", "#F0E442")
   ggplt <- ggpltBASE +   
     geom_point(data=DF3,aes(x=xvar,y=yvar),colour=cbbPalette[7],alpha=1/5,na.rm=TRUE) +
     geom_point(data=DFA,aes(x=xvar,y=yvar),colour=cbbPalette[3],alpha=1/20,na.rm=TRUE)      
     
     #scale_colour_gradientn(colours = cbPalette2)
   
    
   
   ggplts <- ggplt + 
     xlab(  expression ( paste(phi, " [deg]" ) ) ) + 
     ylab(  expression ( paste(psi, " [deg]" ) ) ) + 
     scale_x_continuous(expand=c(0.01,0), limits=c(-180,0)) +
     scale_y_continuous(expand=c(0.01,0.01), limits=c(-180,0)) +
     stat_contour(aes(x=x,y=y,z=zA,fill=..level..),colour=cbbPalette[3],
                  size=1.0,
                  breaks = c(0.05, 0.1, 0.2, 0.4, 0.6, 1.0)) +
     stat_contour(aes(x=x,y=y,z=z3,fill=..level..),colour=cbbPalette[7],
                  size=1.0,
                  breaks = c(0.1, 0.2, 0.4, 0.6, 1.0)) +
     scale_fill_continuous(guide = 'none')
   
   print(ggplts)
   
   filename <- "phipsi.pdf"
   pdf(filename,width=8.0,height=6.0)
   print(ggplts)
   dev.off()
}


# RMSD from initial structure
if (plotrmsd) {
  rmsdA <- rmsd(hisA$xyz.models[1,],hisA$xyz.models,fit=TRUE)
  rmsd3 <- rmsd(his3$xyz.models[1,],his3$xyz.models,fit=TRUE)
  # Moving average
  rmsdvAVGA <- movingAverage(rmsdA,20,centered=TRUE)
  rmsdvAVG3 <- movingAverage(rmsd3,20,centered=TRUE)
  rmsdFrame<-data.frame(x=x,yA=rmsdA,avgA=rmsdvAVGA,y3=rmsd3,avg3=rmsdvAVG3)  
  ggplt <- ggplot(data=rmsdFrame,aes(x,yA)) +
    geom_line(col=cbbPaletteL[3]) +
    geom_line(aes(x,avgA),col=cbbPalette[3],size=2) +
    geom_line(aes(x,y3),col=cbbPaletteL[7],size=0.4) + 
    geom_line(aes(x,avg3),col=cbbPalette[7],size=2) +
    xlab("Timestep [ns]") + 
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
  
  filename <- "rmsd.pdf"
  pdf(filename,width=8.0,height=6.0)
  print(ggplt)
  dev.off()
  
  #hist(rmsd1, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
  #lines(density(rmsd1), col="gray", lwd=3)
  
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

# # RMSD Matrix Only for 310
# recomputeClusters<-FALSE
# if (recomputeClusters) {
#   indexes <- seq(1,Nsteps,by=100)   
#   vRMSD <- vector()
#   for (index in indexes) {
#     vRMSD <- c(vRMSD,rmsd(his3$xyz.models[index,],his3$xyz.models[indexes,],fit=TRUE)) 
#     print(paste("Step",index,sep=" "))
#   }
#   RMSDMatrix <- matrix(vRMSD,nrow = length(indexes))
#   # Clustering
#   RMSDDist <- as.dist(RMSDMatrix)
# 
#   filename <- "Clustering/RMSDDist.Rdata"
#   save(RMSDDist,file=filename)
# } else {
#   load("Clustering/RMSDDist.Rdata")
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
#   for (k in 3:5)
#     asw[[k]] <- pam(RMSDDist, k) $ silinfo $ avg.width
#   k.best <- which.max(asw)
#   cat("silhouette-optimal number of clusters:", k.best, "\n")
#   k <- k.best
#   fitpam<-pam(RMSDDist,k=k)
# 
#   # Write out medoids and full clusters
#   #lista<-c(1:length(his3$xyz.models[,1]))
#   for (i in 1:length(fitpam$medoids)) {
#     index <- indexes[fitpam$medoids[i]]
#     filename <- paste("Clustering/Medoid_",i,"_frame_",index,".pdb",sep="")
#     write.pdb(file=filename,xyz=his3$xyz.models[index,])
#     clusterindexes <- indexes[fitpam$clustering==i]
#     filename <- paste("Clustering/Cluster_",i,".pdb",sep="")
#     write.pdb(file=filename,xyz=his3$xyz.models[clusterindexes,])  
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
  
  # Load Alpha Distributions 
  dists <- list()
  AllDists <- list()
  dists$ang<-read.table(paste(BaseDirAlpha,"Distributions/dists/Statistics/PDB/Bond_angle_CA-CA-CA.dat",sep=""))
  dists$dih<-read.table(paste(BaseDirAlpha,"Distributions/dists/Statistics/PDB/Dihedral_angle_CA-CA-CA-CA.dat",sep=""))
  dists$HB13<-read.table(paste(BaseDirAlpha,"Distributions/dists/Statistics/PDB/Atom_Distance_2_CA-CA.dat",sep=""))
  dists$HB14<-read.table(paste(BaseDirAlpha,"Distributions/dists/Statistics/PDB/Atom_Distance_3_CA-CA.dat",sep=""))
  dists$HB15<-read.table(paste(BaseDirAlpha,"Distributions/dists/Statistics/PDB/Atom_Distance_4_CA-CA.dat",sep="")) 
  AllDists<-lappend(AllDists,dists)
  # REferencedists
  dists$ang<-read.table(paste(BaseDirAlpha,"ReferenceDists/bondang_NPDB.dat",sep=""))
  dists$dih<-read.table(paste(BaseDirAlpha,"ReferenceDists/dihang_NPDB.dat",sep=""))
  dists$HB13<-read.table(paste(BaseDirAlpha,"ReferenceDists/r13_NPDB.dat",sep=""))
  dists$HB14<-read.table(paste(BaseDirAlpha,"ReferenceDists/r14_NPDB.dat",sep=""))
  dists$HB15<-read.table(paste(BaseDirAlpha,"ReferenceDists/r15_NPDB.dat",sep=""))
  AllDists<-lappend(AllDists,dists)
  dA<-as.data.frame(AllDists)
  # Load 310 Distributions 
  dists <- list()
  AllDists <- list()
  dists$ang<-read.table(paste(BaseDir310,"Distributions/dists/Statistics/PDB/Bond_angle_CA-CA-CA.dat",sep=""))
  dists$dih<-read.table(paste(BaseDir310,"Distributions/dists/Statistics/PDB/Dihedral_angle_CA-CA-CA-CA.dat",sep=""))
  dists$HB13<-read.table(paste(BaseDir310,"Distributions/dists/Statistics/PDB/Atom_Distance_2_CA-CA.dat",sep=""))
  dists$HB14<-read.table(paste(BaseDir310,"Distributions/dists/Statistics/PDB/Atom_Distance_3_CA-CA.dat",sep=""))
  AllDists<-lappend(AllDists,dists)
  # REferencedists
  ang<-read.table(paste(BaseDir310,"ReferenceDists/bondang_NPDB.dat",sep=""))
  dists$ang <- spline(ang[,1]*180/pi,ang[,2],n=length(dists$ang[,1]))
  dih<-read.table(paste(BaseDir310,"ReferenceDists/dihang_NPDB.dat",sep=""))
  dists$dih <- spline(dih[,1]*180/pi,dih[,2],n=length(dists$dih[,1]))
  HB13<-read.table(paste(BaseDir310,"ReferenceDists/r13_NPDB.dat",sep=""))
  dists$HB13 <- spline(HB13[,1],HB13[,2],n=length(dists$HB13[,1]))
  HB14<-read.table(paste(BaseDir310,"ReferenceDists/r14_NPDB.dat",sep=""))  
  dists$HB14 <- spline(HB14[,1],HB14[,2],n=length(dists$HB14[,1]))
  AllDists<-lappend(AllDists,dists)  
  d3<-as.data.frame(AllDists)
  
  
  ggpltBASE <- ggplot(data=dA) + 
    scale_y_continuous(breaks=c(0.0,0.5,1.0)) +
    theme(panel.background = element_rect(fill = 'white', colour="black", size=1),
          axis.title.x = element_text(face="bold", colour="black", size=30), 
          axis.ticks.x = element_line(size=1.2, colour="black"), 
          axis.text.x  = element_text(angle=0, vjust=0.5, size=20, face="bold", colour="black"), 
          axis.title.y = element_text(face="bold", colour="black", size=30, angle=90), 
          axis.ticks.y = element_line(size=1.2, colour="black"), 
          axis.text.y  = element_text(angle=0, vjust=0.5, size=20, face="bold", colour="black"), 
          panel.grid = element_blank(),
          plot.title = element_text(lineheight=3, face="bold", color="black", size=30), 
          legend.title  = element_blank(), 
          legend.text = element_text(lineheight=3, face="bold", color="black", size=30))
  
  pang <- ggpltBASE + 
    geom_line(data=dA,aes(x=ang.V1, y = ang.V2/max(ang.V2)),col=cbbPalette[3],size=2) + 
    geom_line(data=dA,aes(x=ang.V1, y = ang.V2.1/max(ang.V2.1)),col=cbbPaletteL[3],size=2,linetype="dashed") + 
    geom_line(data=d3,aes(x=ang.V1, y = ang.V2/max(ang.V2)),col=cbbPalette[7],size=2) + 
    geom_line(data=d3,aes(x=ang.x, y = ang.y/max(ang.y)),col=cbbPaletteL[7],size=2,linetype="dashed") + 
    xlim(75,110) +
    xlab(  expression ( paste(theta, " [deg]" ) ) ) + 
    ylab("")   
  
  pdih <- ggpltBASE +
    geom_line(data=dA,aes(x=dih.V1, y = dih.V2/max(dih.V2)),col=cbbPalette[3],size=2) + 
    geom_line(data=dA,aes(x=dih.V1, y = dih.V2.1/max(dih.V2.1)),col=cbbPaletteL[3],size=2,linetype="dashed") + 
    geom_line(data=d3,aes(x=dih.V1, y = dih.V2/max(dih.V2)),col=cbbPalette[7],size=2) + 
    geom_line(data=d3,aes(x=dih.x, y = dih.y/max(dih.y)),col=cbbPaletteL[7],size=2,linetype="dashed") + 
    xlim(30,100) +
    xlab(  expression ( paste(phi, " [deg]" ) ) ) + 
    ylab("")
    
    mmm <- expression ( paste(r["i,i+2"], " [", ring(A), "]" ) )
  pHB13 <- ggpltBASE + 
    geom_line(data=dA,aes(x=HB13.V1, y = HB13.V2/max(HB13.V2)),col=cbbPalette[3],size=2) + 
    geom_line(data=dA,aes(x=HB13.V1, y = HB13.V2.1/max(HB13.V2.1)),col=cbbPaletteL[3],size=2,linetype="dashed") + 
    geom_line(data=d3,aes(x=HB13.V1, y = HB13.V2/max(HB13.V2)),col=cbbPalette[7],size=2) + 
    geom_line(data=d3,aes(x=HB13.x, y = HB13.y/max(HB13.y)),col=cbbPaletteL[7],size=2,linetype="dashed") + 
    xlim(4.5,6.5) +  xlab( expression ( paste(r["i,i+2"], " [", ring(A), "]" ) ) ) +  ylab("")
  
  pHB14 <- ggpltBASE + 
    geom_line(data=dA,aes(x=HB14.V1, y = HB14.V2/max(HB14.V2)),col=cbbPalette[3],size=2) + 
    geom_line(data=dA,aes(x=HB14.V1, y = HB14.V2.1/max(HB14.V2.1)),col=cbbPaletteL[3],size=2,linetype="dashed") + 
    geom_line(data=d3,aes(x=HB14.V1, y = HB14.V2/max(HB14.V2)),col=cbbPalette[7],size=2) + 
    geom_line(data=d3,aes(x=HB14.x, y = HB14.y/max(HB14.y)),col=cbbPaletteL[7],size=2,linetype="dashed") + 
    xlim(4,7.0) +  xlab(  expression ( paste(r["i,i+3"], " [", ring(A), "]" ) ) ) +   ylab("")
 
   pHB15 <- ggpltBASE + 
    geom_line(aes(x=HB15.V1, y = HB15.V2/max(HB15.V2)),col=cbbPalette[3],size=2) +
    geom_line(aes(x=HB15.V1, y = HB15.V2.1/max(HB15.V2.1)),col=cbbPaletteL[3],size=2,linetype="dashed") +
    xlim(5,7.5) +  xlab(  expression ( paste(r["i,i+4"], " [", ring(A), "]" ) ) ) +   ylab("")


  filename <- "Distributions.pdf"
  pdf(filename,width=8.0,height=6.0)
  grid.arrange(pang,pdih,pHB13,pHB14,pHB15)

  dev.off()
}
