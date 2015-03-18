 options(echo=TRUE) # if you want see commands in output file
 args <- commandArgs(trailingOnly = TRUE)
 print(args)

library("XML") # required for XML parsing
library("gtools") # required for mixedsort
library("lattice") # required for levelplot
library("xtable") # required for latex output
library("fields") # required for plot colorscale
#library("corrplot") # 
library("ggplot2")
library("gridExtra")
library("reshape")

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

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
 c10 <- c("dodgerblue2",
          "#E31A1C", # red
          "green4",
          "#6A3D9A", # purple
          "#FF7F00", # orange
          "black",
          "gold1",
          "skyblue2",
          "#FB9A99", # lt pink
          "gray70"
          )

LoadData <- function() {




   # Distribution in the same order than in Main.XML
   #DistributionsType <- list("r13","")
   # Load Best iteration MC
   BestIterationN <- NA
   if (MC || MCIBI) {
     BestMC <- xmlInternalTreeParse("OUTPUT/BestIteration_Output.xml");
     BestIterationN <- as.numeric(xpathApply(BestMC, "//c_run",xmlValue))+1;
     BestIteration <- paste("OUTPUT/r",BestIterationN,sep="")     
     print(BestIterationN);
   }
   
   
   # Loop over runs
   DistribALLSpl <- list()
   files <- list.files(path="OUTPUT", pattern="r*/*_Input.xml", full.names=T, recursive=TRUE)
   SortedFiles<-mixedsort(files)
   n_of_runs <- length(files)
   
   
   # Load r15 only for 310  helices
   R15Dists310 <- list()
   for (run in 1:n_of_runs) {
     filename <- paste('OUTPUT/r',run-1,'/RECOMPUTED/R15.dat',sep="")
     distrib <- as.data.frame(read.table(filename))
     R15Dist310.spl <- spline(distrib[,1],distrib[,2]/(sum(distrib[,2])*diff(distrib[,1])[1]),n=2*length(distrib[,1]))
     R15Dists310 <- lappend(R15Dists310,R15Dist310.spl)
   }
   filename <- "../R15310/R15Ref.dat"
   distrib <-  as.data.frame(read.table(filename))
   R15Ref310 <- spline(distrib[,1],distrib[,2]/(sum(distrib[,2])*diff(distrib[,1])[1]),n=2*length(distrib[,1]))
   filename <- "../R15310/R15Giulia.dat"
   distrib <-  as.data.frame(read.table(filename))
   R15Giulia310 <- spline(distrib[,1],distrib[,2]/(sum(distrib[,2])*diff(distrib[,1])[1]),n=2*length(distrib[,1]))
 
   
   
   if (Alpha) { 
     potlist<-c("null","r12","r13","r14","theta","phi","vdw")  
   } else {
     potlist<-c("null","r12","r13","theta","phi","vdw")  
   }
   
   
   #  1.  Load distributions and losses
   doc0 <- xmlInternalTreeParse(SortedFiles[1]);
   src0 <- xpathApply(doc0, "//input/param")
   NumberOfPotentials <- xmlSize(src0)
   bOptimize <- vector()
   refDistribs <- list()
   BestMCDistribs <- list()  
   ipotindex <- vector()
   # Giulia distributions
   GiuliaDistribs <- list()
   # Load ref distributions and bOptimize vector
   for (ipot in 1:NumberOfPotentials) {
     TmpDoc0 <- xmlDoc(src0[[ipot]])
     bOptimize <- append(bOptimize,xpathApply(TmpDoc0,"//bOptimize",xmlValue)[[1]])
     if (bOptimize[ipot] == "true") {
      ipotindex <- append(ipotindex,ipot)
      filename <- paste('INPUT/Param',ipot-1,'.dat',sep="")
      refDistrib <- as.data.frame(read.table(filename))
      filename <- paste('GiuliaDistribs/Param',ipot-1,'.dat',sep="")
      GiuliaDistrib <- as.data.frame(read.table(filename))
      if (MC) {
       filename <- paste(BestIteration,'/Param',ipot-1,'.dat',sep="")
       BestMCDistrib <- as.data.frame(read.table(filename))      
      } 
      if ((potlist[ipot]=="theta")||(potlist[ipot]=="phi")) {
        refDistrib[,1] <- refDistrib[,1]/pi*180.0
        GiuliaDistrib[,1] <- GiuliaDistrib[,1]/pi*180.0
        if (MC) {
         BestMCDistrib[,1] <- BestMCDistrib[,1]/pi*180.0
        }
      }
      refDistrib.spl <- spline(refDistrib[,1],refDistrib[,2],n=2*length(refDistrib[,1]))
      refDistrib.spl$y <- refDistrib.spl$y/(sum(refDistrib.spl$y)*diff(refDistrib.spl$x)[1])
      GiuliaDistrib.spl <- spline(GiuliaDistrib[,1],GiuliaDistrib[,2],n=2*length(GiuliaDistrib[,1]))
      GiuliaDistrib.spl$y <- GiuliaDistrib.spl$y/(sum(GiuliaDistrib.spl$y)*diff(GiuliaDistrib.spl$x)[1])
      #refDistrib.spl$y <- refDistrib.spl$y/(max(refDistrib.spl$y))
      refDistribs <- lappend(refDistribs,refDistrib.spl)
      GiuliaDistribs <- lappend(GiuliaDistribs,GiuliaDistrib.spl)
      if (MC) {
       BestMCDistrib.spl <- spline(BestMCDistrib[,1],BestMCDistrib[,2]/(sum(BestMCDistrib[,2])*diff(BestMCDistrib[,1])[1]),n=2*length(BestMCDistrib[,1]))
       #BestMCDistrib.spl <- spline(BestMCDistrib[,1],BestMCDistrib[,2]/(max(BestMCDistrib[,2])),n=2*length(BestMCDistrib[,1]))
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
                if(IBI) {
                  Fitted[,1] <- Fitted[,1]/pi*180.0    
                  ToBeFitted[,1] <- ToBeFitted[,1]/pi*180.0
                }
               }
             distrib[distrib[,1]<0,2]<-0 # remove <0 region
             distrib.spl <- spline(distrib[,1],distrib[,2],n=2*length(distrib[,1]))
             distrib.spl$y <- distrib.spl$y / (sum(distrib.spl$y)*diff(distrib.spl$x)[1])
             #distrib.spl$y <- distrib.spl$y / (max(distrib.spl$y))
             distribs.spl <- lappend(distribs.spl,distrib.spl)
             if (IBI) {
               Fitted.spl <- spline(Fitted[,1],Fitted[,2],n=2*length(Fitted[,1]))
               FittedPMFs <- lappend(FittedPMFs,Fitted.spl)    
               ToBeFittedPMFs <- lappend(ToBeFittedPMFs,ToBeFitted)
             }
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
   
   
   DATA <- list(AllLossFrame=AllLossFrame,DistribALLSpl=DistribALLSpl,potlist=potlist,GiuliaDistribs=GiuliaDistribs,BestMCDistribs=BestMCDistribs,ipotindex=ipotindex,refDistribs=refDistribs,MatrixAllParams=MatrixAllParams,ipotOpti=ipotOpti,ipotNparams=ipotNparams,BestFrame=BestIterationN,R15Dists310=R15Dists310,R15Ref=R15Ref310,R15Giu=R15Giulia310)
   return(DATA)

} # END load data function
 
 

Bestone <- function(DATA,ipot,MC,IBI) {
  MCIBI<-FALSE
  if (MC&&IBI) {
    MCIBI<-TRUE
    MC<-FALSE
    IBI<-TRUE
  }  
  if (MC) { cat("MCSA\n")}  
  if (IBI) { cat("Only IBI\n")}
  if (MCIBI) { cat("MCSA-IBI\n")}
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
  last <- length(DistribALLSpl[[ipot]]$dists)
  step <- round(length(DistribALLSpl[[ipot]]$dists) / 11)
  cat("Number of iterations :",last,"\n")
  if (MC)  {
    BestDistrib <- data.frame(x=BestMCDistribs[[ipot]]$x,y=BestMCDistribs[[ipot]]$y)    
    BestFrame <- DATA$BestFrame
  } else if (MCIBI) {      
    BestFrame <- last - 1
    BestDistrib <- data.frame(x= DistribALLSpl[[ipot]]$dists[[BestFrame]]$x,y= DistribALLSpl[[ipot]]$dists[[BestFrame]]$y)                     
  } else {      
      BestFrame <- last
      BestDistrib <- data.frame(x= DistribALLSpl[[ipot]]$dists[[BestFrame]]$x,y= DistribALLSpl[[ipot]]$dists[[BestFrame]]$y)                     
  }  
  BestOne <- list(BestDistrib=BestDistrib, BestFrame=BestFrame)
  return(BestOne)
}

#  2. Figures about distributions
plotDistributions <- function(DATA) {
  
  DistribALLSpl <- DATA$DistribALLSpl
  AllLossFrame <- DATA$AllLossFrame
  potlist <- DATA$potlist
  GiuliaDistribs <- DATA$GiuliaDistribs  
  refDistribs <- DATA$refDistribs
  ipotindex <- DATA$ipotindex
  MatrixAllParams <- DATA$MatrixAllParams
  ipotOpti <- DATA$ipotOpti
  ipotNparams <- DATA$ipotNparams
  R15Dists310 <- DATA$R15Dists310
  R15Ref310 <- DATA$R15Ref
  R15Giulia310 <- DATA$R15Giu   

  thm <- theme(panel.background = element_rect(fill = 'white'),
               panel.border = element_rect(colour = "black", fill=NA, size=1.5),
               axis.ticks.x = element_line(size=1.2, colour="black"),
               axis.title.x = element_blank(),
               axis.text.x  = element_text(angle=0, vjust=1., size=30, face="bold", colour="black"),
               axis.ticks.y = element_blank(),
               axis.text.y  = element_blank(),
               axis.title.y = element_blank(),
               plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
               legend.title  = element_blank(),
               legend.text = element_blank(),
               legend.position = "none")
  
  for (ipot in 1:length(DistribALLSpl)) {
    NaccIter <- length(DistribALLSpl[[ipot]]$losses)
    ranges <- c()
    if (Alpha) {
      if (potlist[ipotindex[ipot]]=="r12") {
        ranges[1] <- 4.5
        ranges[2] <- 6.5
      } else if (potlist[ipotindex[ipot]]=="r13") {
        ranges[1] <- 4.0
        ranges[2] <- 6.5
      } else if (potlist[ipotindex[ipot]]=="r14") {
        ranges[1] <- 5.0
        ranges[2] <- 8.0
      } else if (potlist[ipotindex[ipot]]=="theta") {
        ranges[1] <- 70.0
        ranges[2] <- 110.0
      } else if (potlist[ipotindex[ipot]]=="phi") {
        ranges[1] <- 20.0
        ranges[2] <- 80.0
      }  
    } else {
      if (potlist[ipotindex[ipot]]=="r12") {
        ranges[1] <- 4.5
        ranges[2] <- 6.5
      } else if (potlist[ipotindex[ipot]]=="r13") {
        ranges[1] <- 4.0
        ranges[2] <- 7.0    
      } else if (potlist[ipotindex[ipot]]=="theta") {
        ranges[1] <- 70.0
        ranges[2] <- 110.0
      } else if (potlist[ipotindex[ipot]]=="phi") {
        ranges[1] <- 20.0
        ranges[2] <- 120.0
      }         
    }  
    
    last <- length(DistribALLSpl[[ipot]]$dists)
    step <- round(length(DistribALLSpl[[ipot]]$dists) / 11)
    # Reshape into data frame
    DistFrame <- data.frame()  
    for (i in seq(1,last,by=step) ) {    
      distrib <- DistribALLSpl[[ipot]]$dists[[i]]
      mm <- matrix(0,nrow=length(distrib$x),ncol=3)
      mm[,1] <- distrib$x
      mm[,2] <- distrib$y
      mm[,3] <- rep(i,length(distrib$x))
      DistFrame <- rbind(DistFrame, mm)
    }
    colnames(DistFrame)<-c("x","y","run")
    
    # Referece
    refDistrib <- data.frame(x=refDistribs[[ipot]]$x,y=refDistribs[[ipot]]$y)
    # Gulia
    GiuliaDistrib <- data.frame(x=GiuliaDistribs[[ipot]]$x,y=GiuliaDistribs[[ipot]]$y)
    # Bestone
    BestOne <- Bestone(DATA,ipot,MC,IBI)
    BestDistrib <- BestOne$BestDistrib
    BestFrame <- BestOne$BestFrame    
    
    
    if (MCIBI) {
      thm <- thm 
    } else {
      thm <- thm + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }
    
    plt <- ggplot(data=DistFrame) +   
      geom_line(aes(x,y,group=run,alpha=run),size=3,color=c25[ipot]) + 
      geom_line(data=refDistrib,aes(x,y),size=3,color="black") +
      geom_line(data=GiuliaDistrib,aes(x,y),size=3,color=c10[7],linetype=2) +
      geom_line(data=BestDistrib,aes(x,y),size=3.5,color="black") +
      geom_line(data=BestDistrib,aes(x,y),size=2.5,color="#ff5eff") +
      xlim(ranges) + 
      thm
    
    filename <- paste("Distributions/",Type,"/",IBIn,MCn,"-Param",ipotindex[ipot]-1,".pdf",sep="")    
    #png(filename,width=2000,bg="transparent")
    pdf(filename,width=8.0,height=4.0)
    print(plt)
    dev.off()
  }
 
# PLOT DISTRIBUTIONS FOR R15 310
    DistFrameR15 <- data.frame()  
    for (i in seq(1,last,by=step) ) {    
      distrib <- R15Dists310[[i]]
      mm <- matrix(0,nrow=length(distrib$x),ncol=3)
      mm[,1] <- distrib$x
      mm[,2] <- distrib$y
      mm[,3] <- rep(i,length(distrib$x))
      DistFrameR15 <- rbind(DistFrameR15,mm)
    }
    colnames(DistFrameR15)<-c("x","y","run")
    
    # Referece
    refDistrib <- data.frame(x=R15Ref310$x,y=R15Ref310$y)
    # Gulia
    GiuliaDistrib <- data.frame(x=R15Giulia310$x,y=R15Giulia310$y)
    
    if (MCIBI) {
      thm <- thm 
    } else {
      thm <- thm + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }
    ranges[1] <- 5.0
    ranges[2] <- 10.0
    plt <- ggplot(data=DistFrameR15) +   
      geom_line(aes(x,y,group=run,alpha=run),size=3,color=c25[3]) + 
      geom_line(data=refDistrib,aes(x,y),size=3,color="black") +
      geom_line(data=GiuliaDistrib,aes(x,y),size=3,color=c10[7],linetype=2) +
      xlim(ranges) + 
      thm
    
    filename <- paste("Distributions/310/",IBIn,MCn,"-R15.pdf",sep="")    
    #png(filename,width=2000,bg="transparent")
    pdf(filename,width=8.0,height=4.0)
    print(plt)
    dev.off()

} # end plotDistributions

#  5. Figures about Loss function
plotLoss <- function(DATA) {
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
  
  thm <- theme(panel.background = element_rect(fill = 'white'),
             panel.border = element_rect(colour = "black", fill=NA, size=1.5),
             axis.ticks.x = element_line(size=1.2, colour="black"),
             axis.title.x = element_blank(),
             axis.text.x  = element_text(angle=0, vjust=1., size=30, face="bold", colour="black"),
             axis.ticks.y = element_line(size=1.2, colour="black"),
             axis.text.y  = element_text(angle=0, vjust=0.0, size=30, face="bold", colour="black"),
             axis.title.y = element_blank(),
             plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
             legend.title  = element_blank(),
             legend.text = element_blank(),
             legend.position = "none",
             panel.grid.major.y = element_line(colour="black",size=0.8),
             panel.grid.major.x = element_blank())
  minx<-1
  maxx<-0
  miny<-0
  maxy<-5
  LossFrame <- data.frame()   
  if (MC) {   AvgLoss<-rowMeans(AllLossFrame,na.rm = TRUE) 
  }  else { 
    AvgLoss <- matrix(0,nrow=1,ncol=length(DistribALLSpl[[1]]$losses))
    for (ipot in 1:length(DistribALLSpl)) {
      AvgLoss <- AvgLoss + as.vector(DistribALLSpl[[ipot]]$losses)      
    }
    AvgLoss <- AvgLoss / 5.0
  }
  
  Iter<-c(1:length(AvgLoss))
  mm <- matrix(0,nrow=length(AvgLoss)-1,ncol=2)
  mm[,1] <- t(Iter[-1])
  mm[,2] <- t(AvgLoss[-1])
  AvgLossFrame <- as.data.frame(mm)
  colnames(AvgLossFrame) <- c("x","y")
  for (ipot in 1:length(DistribALLSpl)) {
    loss <- DistribALLSpl[[ipot]]$losses
    acceptedIterations <- DistribALLSpl[[ipot]]$acceptedIterations
    if (max(acceptedIterations,na.rm=TRUE)>maxx)
      maxx <- max(acceptedIterations,na.rm=TRUE)  
    
    Acc <- DistribALLSpl[[ipot]]$acceptedIterations
    loss <- DistribALLSpl[[ipot]]$losses
    mm <- matrix(0,nrow=length(Acc),ncol=3)
    mm[,1] <- Acc
    mm[,2] <- loss + 1*ipot-1
    mm[,3] <- rep(ipot,length(Acc))
    LossFrame <- rbind(LossFrame, mm)    
  }  
  colnames(LossFrame)<-c("x","y","ipot")

  plt <- ggplot(data=LossFrame,aes(x,y)) + 
  geom_line(color="black",size=3) + 
  geom_line(aes(color=factor(ipot)),size=2) + 
  scale_color_manual(values = c25) +
  facet_grid(ipot ~ .) +
  thm    
  plt <- ggplot(data=LossFrame,aes(x,y,group=ipot)) + 
  geom_line(color="black",size=3) + 
  geom_line(aes(color=factor(ipot),size=2)) + 
  scale_color_manual(values = c25) +
  thm
  filename <- paste("Distributions/",Type,"/",IBIn,MCn,"-Loss.pdf",sep="")    
  #png(filename,width=2000,bg="transparent")
  pdf(filename,width=8.0,height=4.0)
  print(plt)
  dev.off()

  plt <- ggplot(data=AvgLossFrame,aes(x,y)) + geom_line(size=3) + thm
  filename <- paste("Distributions/",Type,"/",IBIn,MCn,"-AVGLoss.pdf",sep="")    
  #png(filename,width=2000,bg="transparent")
  pdf(filename,width=8.0,height=4.0)
  print(plt)
  dev.off()
  
  #return(plt)
} # end plotLoss


plotRadar <- function(systems) {
  sys <- systems[[1]]
  StoredData <- paste("/home/pmereghetti/data/projects/2014/CGautoTest/FigForPaper/",sys[1],"-",sys[2],"-",sys[3],"-DATA.RData",sep="")
  load(StoredData)
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
  
  thm <- theme(panel.background = element_rect(fill = 'white'),
               panel.border = element_rect(colour = "black", fill=NA, size=1.5),
               axis.ticks.x = element_line(size=1.2, colour="black"),
               axis.title.x = element_blank(),
               axis.text.x  = element_text(angle=0, vjust=1., size=30, face="bold", colour="black"),
               axis.ticks.y = element_line(size=1.2, colour="black"),
               axis.text.y  = element_text(angle=0, vjust=0., size=30, face="bold", colour="black"),
               axis.title.y = element_blank(),
               plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
               legend.title  = element_blank(),               
               legend.text = element_text(angle=0, vjust=1., size=20, face="bold", colour="black"),
               plot.title =  element_text(angle=0, vjust=0.0, size=20, face="bold", colour="black"),
               plot.margin= unit(c(0, 1, -0.5, 1),'lines'))
  
  isys <- 0
  ipotOpti <- DATA$ipotOpti
  mm <- matrix(0,nrow=ipotOpti,ncol=4) 
  for (sys in systems) {     
    StoredData <- paste("/home/pmereghetti/data/projects/2014/CGautoTest/FigForPaper/",sys[1],"-",sys[2],"-",sys[3],"-DATA.RData",sep="")
    load(StoredData)
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
    
    isys <- isys + 1    
    
    if (isys == 2) {      
      cat("MC")
      BestOne <- Bestone(DATA,ipot,sys[2],sys[1])    
      BestFrame <- BestOne$BestFrame-1    # loss function is -1
      loss <- AllLossFrame[BestFrame,]
      mm[,isys] <- t(loss[ipotindex])
    } else {
      cat("MCIBI or IBI")

      for (ipot in c(1:length(DistribALLSpl))) {        
        BestOne <- Bestone(DATA,ipot,sys[2],sys[1])    
        BestFrame <- BestOne$BestFrame-1 # loss function is - 1
        mm[ipot,isys] <- DistribALLSpl[[ipot]]$losses[BestFrame]
      }
    }    
  }
  
  LossFrame <- as.data.frame(mm)    
  LossFrame[,4]<-potlist[ipotindex]
  if (sys[3]) {    
    labb<- c(expression(r["i,i+2"]), expression(r["i,i+3"]), expression(r["i,i+4"]), expression(theta), expression(phi))
  } else {     
    labb<- c(expression(r["i,i+2"]), expression(r["i,i+3"]),  expression(theta), expression(phi))     
  }
  colnames(LossFrame) <- c("IBI","MCSA","MCSA-IBI","FF")
  LossFrameM <- melt(LossFrame)
  plt <- ggplot(LossFrameM,aes(x=FF,y=sqrt(value),fill=variable)) + 
    geom_bar(stat="identity",position="dodge") + 
    scale_fill_manual(values = cbbPalette) + 
    scale_x_discrete("F",labels=labb) +    
    scale_y_continuous(breaks=seq(0, 1.0, 0.5)) +     
    thm 
  #filename <- paste("Distributions/",Type,"/RadarPlot.pdf",sep="")    
  #png(filename,width=2000,bg="transparent")
  #pdf(filename,width=8.0,height=4.0)  
  #print(plt)
  #dev.off()
  return(plt)
}

#######################################################



# # # SET SYSTEM PARAMs
 if (length(args)<4) {
  stop("Usage : MakePLOSFig IBI MC Alpha reloadXML")
 }
 IBI <- as.logical(args[1])
 MC <- as.logical(args[2])
 Alpha <- as.logical(args[3])
 reloadXML <- as.logical(args[4])
#  IBI<-TRUE
#  MC<-FALSE
#  Alpha<-TRUE
#  reloadXML <- FALSE

IBIn <- ""
MCn <- ""
Type <- "310"
if (IBI) IBIn <- "IBI"
if (MC) MCn <- "MC"
if (Alpha) Type <- "Alpha"

StoredData <- paste("/home/pmereghetti/data/projects/2014/CGautoTest/FigForPaper/",IBI,"-",MC,"-",Alpha,"-DATA.RData",sep="")

MCIBI<-FALSE
if (MC&&IBI) {
  MCIBI<-TRUE
  MC<-FALSE
  IBI<-TRUE
}

library("gtable")
cat(StoredData,"\n")
if (reloadXML) {
   DATA <- LoadData()
   save(DATA,file=StoredData)
} else {
   load(StoredData)
   # Distributions
   if (!Alpha) {
     c25 <- c25[-3]     
   }
   plotDistributions(DATA)
   
   # Loss Function
   plotLoss(DATA) 
   
#    # Radar plot
    systems <- list(A=c(TRUE,FALSE,TRUE,FALSE),B=c(FALSE,TRUE,TRUE,FALSE),C=c(TRUE,TRUE,TRUE,FALSE)) #alpha   IBI, MC, MCIBI
    p1<- plotRadar(systems)
    systems <- list(A=c(TRUE,FALSE,FALSE,FALSE),B=c(FALSE,TRUE,FALSE,FALSE),C=c(TRUE,TRUE,FALSE,FALSE)) #310   IBI, MC, MCIBI   
    p2<- plotRadar(systems) 
    
    # Extract legend from first plot    
    legend = gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box")     

    filename <- paste("Distributions/RadarPlot.pdf",sep="")    
    pdf(filename,width=8.0,height=4.0)  
    grid.arrange(arrangeGrob(p1 + theme(legend.position="none") + ggtitle("Alpha"),
                             p2 + theme(legend.position="none") + ggtitle("310") ),
                 legend,widths=c(0.8, 0.2),ncol=2)
    dev.off()
    
}




# # 6.  Best Parameters selection (For montecarlo only)
# if (MC) {
#   # Best params based on SqSum and Var
#   Vm <- vector()
#   SqSum <- vector()
#   VarM <- vector()
#   VarScaled <- vector()
#   for (run in 1:nrow(AllLossFrame)) {
#     SqSum[run] <- sum(AllLossFrame[run,]^2,na.rm=TRUE)
#     #if (length(OptimPotIndex)>1) 
#       VarM[run] <- sd(AllLossFrame[run,],na.rm=TRUE)^2
#     #else
#     #  VarScaled[run] <- 0.0
#   }
#   SqScaled<-(SqSum-min(SqSum))/(max(SqSum)-min(SqSum))
#   #if (length(OptimPotIndex)>1) 
#     VarScaled<-(VarM-min(VarM))/(max(VarM)-min(VarM))
#   Vm <- SqScaled + VarScaled
#   
#   # Sort 
#   #avgSort <- sort(AvgLoss,index.return=TRUE) 
#   avgSort <- sort(Vm,index.return=TRUE) 
#   SortedParametersFrame<-as.data.frame(MatrixAllParams[avgSort$ix[1:min(10,length(avgSort$ix))]-1,])
#   #s[,10] <- (s[,10]+3.14)*180/3.14
#   #s[,8] <- s[,8]*180/3.14
#   all_mean <- sapply(SortedParametersFrame,median)
#   all_quantiles <- sapply(SortedParametersFrame,quantile)
#   
#   SortedParametersFrame <- rbind(SortedParametersFrame,all_quantiles[3,])
#   SortedParametersFrame <- rbind(SortedParametersFrame,all_quantiles[2,])
#   SortedParametersFrame <- rbind(SortedParametersFrame,all_quantiles[4,])
#   rownames(SortedParametersFrame)[nrow(SortedParametersFrame)-2]<-c("Median")
#   rownames(SortedParametersFrame)[nrow(SortedParametersFrame)-1]<-c("25pct")
#   rownames(SortedParametersFrame)[nrow(SortedParametersFrame)]<-c("75pct")
#   # Save out latex table
#   caption <- "Best parameters based on average loss function. Params are sorted on average loss function and the top 5 (max) are selected. Based on MC-SA simulations."
#   tab <- xtable(t(SortedParametersFrame[,]), caption=caption)
#   filename <- "BestBasedOnAvgLoss_MCSA.tex"
#   print(tab,file=filename,append=F,table.placement = "h", caption.placement="bottom")
# 
#   Nstart <- 1
#   for (i in 1:ipotOpti) { 
#     Nparams <- ipotNparams[i]
#     for (j in 1:Nparams) {
#       val <- SortedParametersFrame[nrow(SortedParametersFrame)-2,Nstart:Nstart+j-1]
#       cat("<parameter>",val,"</parameter>\n")
#     }
#     for (j in 1:Nparams) {
#       val <- SortedParametersFrame[nrow(SortedParametersFrame)-1,Nstart:Nstart+j-1]
#       cat("<parameter_min>",val,"</parameter_min>\n")
#     }
#     for (j in 1:Nparams) {
#       val <- SortedParametersFrame[nrow(SortedParametersFrame),Nstart:Nstart+j-1]
#       cat("<parameter_max>",val,"</parameter_max>\n")
#     }
#     Nstart <- Nstart + Nparams
#   }
# }
# 
# 

