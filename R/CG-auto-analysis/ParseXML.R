library("XML") # required for XML parsing
library("gtools") # required for mixedsort
library("lattice") # required for levelplot
library("xtable") # required for latex output
library("fields") # required for plot colorscale

# Loop over runs
ParamAllIterations <- list()
ParamAllIterationsMC <- list()
ParamAllIterationsIBI <- list()
DistribAllMomenta <- list()
files <- list.files(path="OUTPUT", pattern="r*/*_Input.xml", full.names=T, recursive=TRUE)
SortedFiles<-mixedsort(files)
n_of_runs <- length(files)
n_of_runs_MC <- 0
n_of_runs_IBI <- 0
for (run in 1:n_of_runs) {

  # Parse XML
  doc <- xmlInternalTreeParse(SortedFiles[run]);
  src <- xpathApply(doc, "//input/param")
  NumberOfPotentials <- xmlSize(src)
  ParamIteration_i <- list()
  Momenta_i <- list()

  # Check which optimization method is used
  optim <- xpathApply(doc,"//optim",xmlValue)[[1]]
  
  TotalNumberOfOptimizedParameters <- 0
  NumberOfOptimizedPotentials <- 0
  ParametersNames <- list()
  DistNames <- list()
  OptimPotIndex <- vector()
  ii <- 0
  for (i in 1:NumberOfPotentials) {
   
     TmpDoc <- xmlDoc(src[[i]])
     dist_name <- xpathApply(TmpDoc,"//dist_name",xmlValue)[[1]]

     
     if (length(xpathApply(TmpDoc,"//pot_name",xmlValue)) != 0)
        pot_name <- xpathApply(TmpDoc,"//pot_name",xmlValue)[[1]]
     
     bOptimize <- xpathApply(TmpDoc,"//bOptimize",xmlValue)[[1]]
     DistNames[i] <- (paste(dist_name,"_Loss",sep=""))

     # Store parameters only for optimized potentials 
     if(bOptimize=="true") {
        ii <- ii + 1
        OptimPotIndex[ii] <- i
        dist_fname <- xpathApply(TmpDoc,"//dist_fname",xmlValue)[[1]]
        ref_dist_fname <- xpathApply(TmpDoc,"//ref_dist_fname",xmlValue)[[1]]
        LossFunction <- as.numeric(xpathApply(TmpDoc,"//LossFunction",xmlValue)[[1]])
        xmin <- as.numeric(xpathApply(TmpDoc,"//xmin",xmlValue)[[1]])
        xmax <- as.numeric(xpathApply(TmpDoc,"//xmax",xmlValue)[[1]])
        q_unit <- xpathApply(TmpDoc,"//q_unit",xmlValue)[[1]]
        if (optim=="MC") { 
          status <- xpathApply(TmpDoc,"//status",xmlValue)[[1]] 
        } 
        else if (optim=="IBI") {     status <- "Accepted"     }
	else { status <- "NA" }
        NumberOfOptimizedPotentials <- NumberOfOptimizedPotentials + 1
        if (optim=="MC") { 
          if (length(xpathApply(TmpDoc,"//pot_parameter_try",xmlValue))==0) {
               pot_parameters <- NA
          } else {   pot_parameters <- as.numeric(xpathApply(TmpDoc,"//pot_parameter_try",xmlValue))}
        } else if (optim=="IBI") {
          pot_parameters <- as.numeric(xpathApply(TmpDoc,"//pot_parameter",xmlValue))
        }
        if (length(xpathApply(TmpDoc,"//pot_parameter_FITTING",xmlValue))!=0)
          potparFIT <- as.numeric(xpathApply(TmpDoc,"//pot_parameter_FITTING",xmlValue))
        else
          potparFIT <- NA
        param <- list(pname=pot_name, dist_name=dist_name, dist_fname=dist_fname, ref_dist_fname=ref_dist_fname, loss=LossFunction, params=pot_parameters, xmin=xmin, xmax=xmax, status=status, q_unit=q_unit, paramsFIT=potparFIT)
        ParamIteration_i[[i]] <- param
        if (length(ParametersNames)==0) {
          ParametersNames <- c( paste(dist_name,"_P",1:length(param$params),sep="") )
	} else { 
          ParametersNames <- c( ParametersNames, paste(dist_name,"_P",1:length(param$params),sep="") )
        }
        TotalNumberOfOptimizedParameters <- TotalNumberOfOptimizedParameters + length(param$params)

        # Load distribution and store IQR 
        if (run<n_of_runs) {
          filename=paste("OUTPUT/r",run,"/dists/Statistics/PDB/",dist_fname,sep="")
          cat(filename,"\n")
          Draw <- read.table(filename)
          distrib <- Draw[,2]
          forMedian <- which(cumsum(distrib/sum(distrib))<0.5)
          forQR25 <- which(cumsum(distrib/sum(distrib))<0.25)
          forQR75 <- which(cumsum(distrib/sum(distrib))<0.75)
          QR25 <- Draw[c(length(forQR25)),1]
          QR75 <- Draw[c(length(forQR75)),1]
          Median <- Draw[c(length(forMedian)),1]
          IQR <- QR75-QR25
          moments <- IQR
          momenta <- list(pname=pot_name,dist_fname=dist_fname,moments=moments)
          Momenta_i[[i]] <- momenta
        }

     } else {
        dist_fname <- NA 
        pot_name <- NA
        LossFunction <- NA 
        pot_parameters <- NA
        param <- list(pname=pot_name, dist_name=dist_name, dist_fname=dist_fname, loss=LossFunction, params=pot_parameters)
        ParamIteration_i[[i]] <- param
        momenta <- list(pname=pot_name,dist_fname=dist_fname,moments=NA)
        Momenta_i[[i]] <- momenta
     }
    
  }
  #ParamAllIterations[[run]] <- ParamIteration_i 
  if (optim=="MC") {
     n_of_runs_MC <- n_of_runs_MC + 1
     ParamAllIterationsMC[[n_of_runs_MC]] <- ParamIteration_i 
     if (run<n_of_runs) {
       DistribAllMomenta[[n_of_runs_MC]] <- Momenta_i
     }
  } else if (optim=="IBI") {
     n_of_runs_IBI <- n_of_runs_IBI + 1
     ParamAllIterationsIBI[[n_of_runs_IBI]] <- ParamIteration_i 
  }
}
sf <- paste(Sys.getenv("Rscripts"),"AnaAndPlot.R",sep="")
source(sf)

#ParamAllIterations <- ParamAllIterationsMC
#n_of_runs_MC <- n_of_runs
#optim <- "MC"

# Open report file
ReportTexFID <- file("Report.tex","w+")
# Report header
   cat(
"\\documentclass[journal=jpcbfk,manuscript=article]{achemso}
\\usepackage{graphicx}
\\usepackage{pslatex}
\\usepackage{array}
\\usepackage{tabularx}
\\usepackage{indentfirst}
\\usepackage{wasysym}
\\usepackage{lastpage}
\\usepackage{longtable}
\\usepackage{subfigure}
\\addtolength{\\topmargin}{-0.87in}
\\addtolength{\\textheight}{1.0in}
\\newcommand{\\squeezeup}{\\vspace{-5.5mm}}\n
\\title{CG-autoparam Report}\n
\\begin{document}\n",file=ReportTexFID)

if (length(ParamAllIterationsMC)!=0) {
  AnaAndPlot(ParamAllIterationsMC,ParametersNames,DistNames,TotalNumberOfOptimizedParameters,NumberOfOptimizedPotentials,OptimPotIndex,n_of_runs_MC,optim="MC",ReportTexFID,DistribAllMomenta)
}

if (length(ParamAllIterationsIBI)!=0) {
  AnaAndPlot(ParamAllIterationsIBI,ParametersNames,DistNames,TotalNumberOfOptimizedParameters,NumberOfOptimizedPotentials,OptimPotIndex,n_of_runs_IBI,optim="IBI",ReportTexFID,DistribAllMomenta)
}

cat("\\end{document}",file=ReportTexFID)
close(ReportTexFID)

######ParamAllIterations <- ParamAllIterationsMC
######n_of_runs <- n_of_runs_MC
######optim <- "MC"
######
######NumberOfPotentials <- length(ParamAllIterations[[1]])
####### Cast parameters into dataframe for easier operations
####### This contains only parameters for optimized potentials
######ParametersFrame <- as.data.frame(matrix(0, n_of_runs, TotalNumberOfOptimizedParameters))
######colnames(ParametersFrame) <- ParametersNames
######
####### Cast loss function int dataframe
######LossFrame <- as.data.frame(matrix(NA, n_of_runs, NumberOfPotentials))
######colnames(LossFrame) <- DistNames
######
####### Fill the frames
######for (run in 1:n_of_runs) {
######  Nstart <- 1
######  for (i in OptimPotIndex) {
######    Npara <- length(ParamAllIterations[[run]][[i]]$params) 
######    ParametersFrame[run,(Nstart:(Nstart+Npara-1))] <- ParamAllIterations[[run]][[i]]$params 
######    LossFrame[run,i] <- ParamAllIterations[[run]][[i]]$loss
######    Nstart <- Nstart + Npara 
######  }
######}
######
######
####### For each optimized potential plot a pairs of selected parmeters
#######Nstart <- 1
#######for (i in 1:NumberafOptimizedPotentials) {
#######   xyz <- ParametersFrame[,(Nstart:(Nstart+Npara-1))]
#######   LossFrame[run,i] <- ParamAllIterations[[run]][[i]]$loss
#######   Nstart <- Nstart + Npara 
#######}
######
####### Compute Average Loss
######AvgLoss <- rowMeans(LossFrame,na.rm=TRUE)
######
######sink("Report.tex")
######
######cat("
######\\documentclass{article}
######\\begin{document}
######")
######
######
####### Best params based on average loss function
######avgSort <- sort(AvgLoss,index.return=TRUE)
######SortedParametersFrame <- ParametersFrame[avgSort$ix,]
#######s[,10] <- (s[,10]+3.14)*180/3.14
#######s[,8] <- s[,8]*180/3.14
######all_mean <- sapply(SortedParametersFrame,mean)
######all_std <- sapply(SortedParametersFrame,sd)
######SortedParametersFrame <- rbind(SortedParametersFrame,all_mean)
######SortedParametersFrame <- rbind(SortedParametersFrame,all_std)
######rownames(SortedParametersFrame)[nrow(SortedParametersFrame)-1]<-c("Mean")
######rownames(SortedParametersFrame)[nrow(SortedParametersFrame)]<-c("StdDev")
####### Save out latex table
######tab <- xtable(SortedParametersFrame[1:10,], caption= "Best parameters based on average loss function. Params are sorted on average loss function and the top 10 are selected")
#######print(tab,file="BestBasedOnAvgLoss.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(tab),by=1))
######print(tab,append=T,table.placement = "h", caption.placement="bottom")
######
######
######
####### Pearson correlation among the top 20%
######PearsonCor <- cor(ParametersFrame[1:ceiling(length(AvgLoss)*0.2),])
######PearsonCor[lower.tri(PearsonCor,diag=T)] = 0.0 # set 0 diagonal and lower tri matrix
######
######
######
####### Save out latex table
######tab <- xtable(PearsonCor, caption= "Pearson Correlation")
#######print(tab,file="BestBasedOnAvgLoss.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(tab),by=1))
######print(tab,append=T,table.placement = "h", caption.placement="bottom")
######cat("
######\\end{document}
######")
######
######sink()
######
####### Split MC plot from IBI plot
####### based on n_of_runs_MC and n_of_runs_IBI
####### MC runs : MCrunStart to MCrunEND
####### IBI runs : IBIrunStart to IBIrunEND
######
######
######if (optim=="MC") {
######   for (i in OptimPotIndex) { 
######      # Number of accepted iterations for each potential
######      NAccIter <- 0
######      for (run in 1:n_of_runs) {
######        if (ParamAllIterations[[run]][[i]]$status=="Accepted") {
######           NAccIter <- NAccIter + 1
######        }
######      }
######
######      # Set colorscale
######      colors <- colfunc(NAccIter)
###### 
######      filename <- paste("ps/",ParamAllIterations[[1]][[i]]$dist_name,"_Monte-Carlo.eps",sep='')
######      postscript(filename)
######
######      # Define Multiplot
######      par(mfrow=c(2,2),mai=c(0.4,0.4,0.4,0.4),mgp=c(1.2,0.4,0.0),oma=c(1,1,3,1))
######
######      # Plot 1. Distributions colored by average Loss and status == Accepted + Ref distribution
######      filename <- paste('INPUT/reference_distributions/',ParamAllIterations[[1]][[i]]$ref_dist_fname,sep="")
######      cat ("Loading :",filename,"\n")
######      distrib <- as.data.frame(read.table(filename))
######      if (ParamAllIterations[[1]][[i]]$q_unit=="deg") {
######        xRange <- c(ParamAllIterations[[1]][[i]]$xmin,ParamAllIterations[[1]][[i]]$xmax)*180/pi
######      } else {
######        xRange <- c(ParamAllIterations[[1]][[i]]$xmin,ParamAllIterations[[1]][[i]]$xmax)
######      }
######      plot(distrib[,1],distrib[,3]/max(distrib[,3]),type="p",pch=1,xlab="r (A)",ylab="P(r)",main="(a)",xlim=xRange)
######      image.plot(legend.only=TRUE, zlim=range(AvgLoss[2:length(AvgLoss)]), col=colors, args.legend = list(x = "topleft"))
######      for (run in avgSort$ix) {
######        if (ParamAllIterations[[run]][[i]]$status=="Accepted") {
######           ParamAllIterations[[run]][[i]]$dist_fname
######           filename <- paste('OUTPUT/r',run-1,'/dists/Statistics/PDB/',ParamAllIterations[[1]][[i]]$dist_fname,sep="")
######           distrib <- as.data.frame(read.table(filename))
######           distrib.spl <- spline(distrib[,1],distrib[,3]/max(distrib[,3]),n=2*length(distrib[,1]))
######           lines(distrib.spl,col=colors[run])
######        }
######      }
######    
######      # Plot 2. Distribution with associated min(AvgLoss) + Top five distribution + Ref. distributions 
######      # These distributions are ranked by average loss so they may have status Rejected
######      filename <- paste('INPUT/reference_distributions/',ParamAllIterations[[1]][[i]]$ref_dist_fname,sep="")
######      distrib <- as.data.frame(read.table(filename))
######      plot(distrib[,1],distrib[,3]/max(distrib[,3]),type="p",pch=1,xlab="r (A)",ylab="P(r)",main="(b)",xlim=xRange)
######      # min(AvgLoss)
######      filename <- paste('OUTPUT/r',(avgSort$ix[1]-1),'/dists/Statistics/PDB/',ParamAllIterations[[avgSort$ix[1]]][[i]]$dist_fname,sep="")
######      distrib <- as.data.frame(read.table(filename))
######      distrib.spl <- spline(distrib[,1],distrib[,3]/max(distrib[,3]),n=2*length(distrib[,1]))
######      lines(distrib.spl,col=2)
######      # Remaining Top 5 AvgLoss
######      last <- min(5,length(AvgLoss)-1)
######      for (ii in 2:last) {
######        filename <- paste('OUTPUT/r',(avgSort$ix[ii]-1),'/dists/Statistics/PDB/',ParamAllIterations[[avgSort$ix[1]]][[i]]$dist_fname,sep="")
######        distrib <- as.data.frame(read.table(filename))
######        distrib.spl <- spline(distrib[,1],distrib[,3]/max(distrib[,3]),n=2*length(distrib[,1]))
######        lines(distrib.spl,col="#5aff5e",lty=2)
######      }
######        
######
######   
######      # Plot 3. Loss function all iterations
######      plot(c(2:nrow(LossFrame)),LossFrame[(2:nrow(LossFrame)),i],xlab="Iterations",ylab="Loss",main="(b)",col="#CCCCCC") # all iterations
######      for (run in 2:nrow(LossFrame)) {
######        if (ParamAllIterations[[run]][[i]]$status=="Accepted") {
######          points(run,LossFrame[run,i]) # only accepted  (if MC)
######        }
######      }
######     
######   
######      # Plot 4. RMS(loss) function
######      RMSLossFrame <- diff(LossFrame[(2:nrow(LossFrame)),1])^2
######      plot(c(1:length(RMSLossFrame)),RMSLossFrame,xlab="Iterations",ylab="RMS(Loss)",main="(c)",col="#CCCCCC") # all iterations
######      for (run in 1:length(RMSLossFrame)) {
######        if (ParamAllIterations[[run]][[i]]$status=="Accepted") {
######            points(run,RMSLossFrame[run]) # only accepted (if MC)
######	}
######      }
######   
######      title(paste("Interaction :",ParamAllIterations[[1]][[i]]$dist_name), outer = TRUE)
######      dev.off()
######   }
######}
######
