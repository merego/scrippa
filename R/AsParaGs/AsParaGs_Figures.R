library("ggplot2")
library("scales")
library("gridExtra")
library("KernSmooth")
library("MASS")


# Aggiornare la tabella se i test sono cambiati.
attributes <- data.frame(read.table("attributes",header=TRUE))

# General theme for PNG (ok for 600 dpi resolution 6.5x3.25in)
thm2<- theme(panel.background = element_rect(fill = 'white'),
             panel.border = element_rect(colour = "black", fill=NA, size=1),
             axis.ticks.x = element_line(size=0.8, colour="black"),
             axis.title.x =  element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.text.x  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.ticks.y = element_line(size=0.8, colour="black"),
             axis.text.y  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.title.y = element_text(angle=90, vjust=0.5, size=12, face="bold", colour="black"),
             plot.title = element_text(lineheight=3, face="bold", color="black", size=30),
             legend.title  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             legend.text = element_text(lineheight=3, face="bold", color="black", size=12),
             strip.text = element_text(lineheight=3, face="bold", color="black", size=12)) # Text for facets header

# Color palette
c5 <- c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")
c10 <- c("#a6cee3",  "#1f78b4",   "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
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

# Append to list
lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
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

# Conversion theta to r
theta2r <- function(theta) {
 l <- 3.8
 r <- 2.0 * l * sin(theta/2)
 return(r)
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

# Theoretical slope for kr/ktheta correations
alpha <- function(theta0) {
 alpha <- 3.8 * (cos(theta0/2))
 return(alpha)
}

# Reference values
theta0d <- 88.0 # deg
theta0r <- deg2rad(theta0d) # rad
r0 <- theta2r(theta0r) # A

# Read reportXX.dat 
# XX == TestIndex
AsParaGs.version <- "new"
ReadData <- function(TestIndex) {
 filename <- sprintf("report%02d.dat",TestIndex)
 REP <- read.table(filename,header=TRUE)
 # Nella vecchia versione di AsParaGs (before 30/06/2015)
 # la loss dell'iterazione t era riferita ai parametri contenuti dell'XML dell'iterazione t-1
 # anche nel report.dat quindi sono sfasati e qui li correggo.
 if (AsParaGs.version=="old") {
   Loss <- REP$AvgLoss[-1]
   Last <- nrow(REP) 
   REP <- REP[-Last,]
   REP$AvgLoss <- Loss
 }   
 # If IBI add WRMS
 if (TestIndex %in% attributes[attributes$SIM=="IBI",1]) { 
   df <- LoadWRMS(TestIndex)
   df.mat <- matrix(df[,1],50,2)
   WRMS.avg <- rowMeans(df.mat)
   REP$WRMS <- WRMS.avg
 }
 
 return(REP)
}

kkLoss <- function(REP,TestIndex) {
  
  # Tests >= 10 have r13,14,theta,phi (<10 only r13,theta)
  if (TestIndex %in% attributes[attributes$r14==TRUE,1]) {
    k1i <- which(colnames(REP)=="Param1") # k13
    k2i <- which(colnames(REP)=="Param1.2") # ktheta
  } else {
    k1i <- which(colnames(REP)=="Param1") # k13
    k2i <- which(colnames(REP)=="Param1.1") # ktheta
  }  
  REP.sub <- cbind(REP[,c(k1i,k2i)],REP$AvgLoss)
  colnames(REP.sub) <- c("k1","k2","AvgLoss")
  return(REP.sub)
}


# Main correlation plot 
PlotCorrelations <- function(REP,TestIndex,fitting=FALSE) {
  
  REP.sub <- kkLoss(REP,TestIndex)    
  
  # Fit intercept of theoretical line
  best.loss <- sort(REP.sub$AvgLoss,index.return=TRUE)
  dfmin <-  REP.sub[best.loss$ix[1:20],] # Top 10
  if (fitting) {
    x <- dfmin[,1]
    y <- dfmin[,2]
    initslope <- (alpha(theta0r))^2 # Theoretical slope
    fit <- nls(y ~ a + b * x, algorithm = "port", start=c(a=120, b=-initslope), upper=c(a=150, b=-(initslope-0.1)), lower=c(a=1, b=-(initslope+0.1)) ) # Fit with 
    keff <- summary(fit)$parameters[1]
    slope <- initslope    
  } else {
    keff <- 115
    slope <- initslope
  }
  
  # Fit intercept of theoretical line (free-fit)
  best.loss <- sort(REP.sub$AvgLoss,index.return=TRUE)
  dfmin <-  REP.sub[best.loss$ix[1:20],] # Top 10
  if (fitting) {    
    x <- dfmin[,1]
    y <- dfmin[,2]    
    # Convert to equivalent harmonic k
    if (TestIndex %in% attributes$Morse) { 
      x <- 2 * x * 1.3^2
    }
    fit.free <- lm(y ~ x)       
  }
  

  # Get some quantiles to create non-linear color scale
  qn = quantile(REP.sub$AvgLoss,c(0.01,0.2,0.5,0.7,0.99))
  qn01<-rescale(qn)
  
  # plot Correlations
  plt <- ggplot(REP.sub) + geom_point(aes(k1,k2,colour=AvgLoss),size=0.5) +
    geom_abline(intercept = keff, slope=-(initslope), size=0.8, color="black") +
    geom_point(data=dfmin,aes(k1,k2),colour="black",size=1) +
    scale_colour_gradientn(name="Average\nLoss", colours= topo.colors(20), values=c(qn01)) +
    scale_x_continuous(expand=c(0.01,0.01)) + # remove white spaces left right
    scale_y_continuous(expand=c(0.01,0.01)) + # remove white spaces bottom top
    thm2
  
  # Add proper xlab, ylab
  ylabel <- expression ( paste ( k[theta], " [ kcal ", mol^"-1" , rad^"-1", "]" , sep = " ")   )
  if (TestIndex %in% attributes$Morse) {    
    xlabel <- expression ( paste ( epsilon[r["i,i+2"]], " [ kcal ", mol^"-1" , Å^"-1", "]" , sep = " ")   )
  } else {
    xlabel <- expression ( paste ( k[r["i,i+2"]], " [ kcal ", mol^"-1" , Å^"-1", "]" , sep = " ")   )    
  }
  plt <- plt + xlab(xlabel) + ylab(ylabel)   
  
  # plot Correlation Lines
  plt1 <- ggplot(REP.sub) + geom_abline(intercept = keff, slope=-(initslope), size=0.8, color="black") +
    geom_point(data=dfmin,aes(k1,k2),colour="black",size=1) +    
    scale_x_continuous(expand=c(0.01,0.01)) + # remove white spaces left right
    scale_y_continuous(expand=c(0.01,0.01)) + # remove white spaces bottom top
    thm2
  
  # Add proper xlab, ylab
  ylabel <- expression ( paste ( k[theta], " [ kcal ", mol^"-1" , rad^"-1", "]" , sep = " ")   )
  if (TestIndex %in% attributes[attributes$Morse==TRUE,1]) {    
    xlabel <- expression ( paste ( epsilon[r["i,i+2"]], " [ kcal ", mol^"-1" , Å^"-1", "]" , sep = " ")   )
  } else {
    xlabel <- expression ( paste ( k[r["i,i+2"]], " [ kcal ", mol^"-1" , Å^"-1", "]" , sep = " ")   )    
  }
  plt1 <- plt1 + xlab(xlabel) + ylab(ylabel)   
  
    
  return(list(Corr=plt,CorrLinesFit=fit.free,R2=summary(fit.free)$adj.r.squared))
}

PlotCorrelationLines <- function(fits) {
  rse<-sqrt(deviance(fit))/(sqrt(8))  
}

# Load dists 
LoadDist <- function(idx,TestIndex) {
  # Position on lpgm-pc is different than t800
  hostname<-Sys.info()["nodename"]
  
  if (hostname=="t800") {
    dir <- sprintf("../test%02d/OUTPUT/r",TestIndex)
  } else  {
    dir <- sprintf("../Data/test%02d/OUTPUT/r",TestIndex)  
  }
  
  df <- data.frame()
  filename <- paste(dir,idx,'/Param1.dat',sep="")
  df <- LoadCheck(filename,df,type="r13",ptype="Dist")
  filename <- paste(dir,idx,'/Param2.dat',sep="")
  df <- LoadCheck(filename,df,type="r14",ptype="Dist")  
  filename <- paste(dir,idx,'/Param3.dat',sep="")
  df <- LoadCheck(filename,df,type="theta",ptype="Dist")  
  filename <- paste(dir,idx,'/Param4.dat',sep="")
  df <- LoadCheck(filename,df,type="phi",ptype="Dist")      
  filename <- paste(dir,idx,'/Param15.dat',sep="")
  df <- LoadCheck(filename,df,type="phi",ptype="Dist")      
  
  
#   filename <- paste(dir,idx,'/Param1.dat',sep="")
#   p1 <- as.data.frame(read.table(filename))
#   filename <- paste(dir,idx,'/Param2.dat',sep="")
#   p2 <- as.data.frame(read.table(filename))
#   filename <- paste(dir,idx,'/Param3.dat',sep="")
#   p3 <- as.data.frame(read.table(filename))
#   filename <- paste(dir,idx,'/Param4.dat',sep="")
#   p4 <- as.data.frame(read.table(filename))
#   filename <- paste(dir,idx,'/Param15.dat',sep="")
#   p15 <- as.data.frame(read.table(filename))  
#   p1.norm <- pnorm(p1,"r13")
#   p2.norm <- pnorm(p2,"r14")
#   p3.norm <- pnorm(p3,"theta")
#   p4.norm <- pnorm(p4,"phi")
#   p15.norm <- pnorm(p15,"r15")
#   df <- data.frame(rbind(p1.norm,p2.norm,p3.norm,p4.norm,p15.norm))
#   df$id <- c(rep("r13",nrow(p1.norm)),rep("r14",nrow(p2.norm)),rep("theta",nrow(p3.norm)),rep("phi",nrow(p4.norm)),rep("r15",nrow(p15.norm)))
  colnames(df) <- c("x","y","id")
  return(df)
}

# Check if file exists load and append (column) to dataframe
LoadCheck <- function(filename,df,type="r13",ptype="PMF") {
  
  if (file.exists(filename)) {
    p <- as.data.frame(read.table(filename))
    if (ptype=="PMF") {
      p[,1] <- ConvertUnit(p[,1],type)
    } else if (ptype=="Dist") {
      p <- pnorm(p,type)
    } else {
      p <- p
    }
    
    tmp <- data.frame(p,rep(type,nrow(p)))

    if (length(df)==0) {
      df <- tmp
    }  else  {
      df <- data.frame(rbind(df,tmp))
    }    
  }  
  return(df)
}
# Load Fitted PMF
LoadPMF <- function(idx,TestIndex) {
  # Position on lpgm-pc is different than t800
  hostname<-Sys.info()["nodename"]
  
  
  if (hostname=="t800") {
    dir <- sprintf("../test%02d/OUTPUT/r",TestIndex)
  } else {
    dir <- sprintf("../Data/test%02d/OUTPUT/r",TestIndex)  
  }
  
  df <- data.frame()
  filename <- paste(dir,idx,'/Param1.dat.fitted.PMF.dat',sep="")
  df <- LoadCheck(filename,df,type="r13",ptype="PMF")
  filename <- paste(dir,idx,'/Param2.dat.fitted.PMF.dat',sep="")
  df <- LoadCheck(filename,df,type="r14",ptype="PMF")  
  filename <- paste(dir,idx,'/Param3.dat.fitted.PMF.dat',sep="")
  df <- LoadCheck(filename,df,type="theta",ptype="PMF")  
  filename <- paste(dir,idx,'/Param4.dat.fitted.PMF.dat',sep="")
  df <- LoadCheck(filename,df,type="phi",ptype="PMF")      
  df.fitted <- df
  
  df <- data.frame()
  filename <- paste(dir,idx,'/Param1.dat.tobefitted.PMF.dat',sep="")
  df <- LoadCheck(filename,df,type="r13",ptype="PMF")
  filename <- paste(dir,idx,'/Param2.dat.tobefitted.PMF.dat',sep="")
  df <- LoadCheck(filename,df,type="r14",ptype="PMF")  
  filename <- paste(dir,idx,'/Param3.dat.tobefitted.PMF.dat',sep="")
  df <- LoadCheck(filename,df,type="theta",ptype="PMF")  
  filename <- paste(dir,idx,'/Param4.dat.tobefitted.PMF.dat',sep="")
  df <- LoadCheck(filename,df,type="phi",ptype="PMF") 
  df.tobefitted <- df
  
  df <- list(df.fitted=df.fitted,df.tobefitted=df.tobefitted)
  return(df)
}


# Load reference dists
LoadRefDist <- function() {
  filename <- "RefDists/Param1.dat"  
  p1 <- read.table(filename)
  # R14
  filename <- "RefDists/Param2.dat"    
  p2 <- read.table(filename)
  # Theta
  filename <- "RefDists/Param3.dat"    
  p3 <- read.table(filename)
  # Phi
  filename <- "RefDists/Param4.dat"    
  p4 <- read.table(filename)
  # R15
  filename <- "RefDists/Param15.dat"      
  p15 <- read.table(filename)  
  p1.norm <- pnorm(p1,"r13")
  p2.norm <- pnorm(p2,"r14")
  p3.norm <- pnorm(p3,"theta")
  p4.norm <- pnorm(p4,"phi")
  p15.norm <- pnorm(p15,"r15")
  df <- data.frame(rbind(p1.norm,p2.norm,p3.norm,p4.norm,p15.norm))
  df$id <- c(rep("r13",nrow(p1.norm)),rep("r14",nrow(p2.norm)),rep("theta",nrow(p3.norm)),rep("phi",nrow(p4.norm)),rep("r15",nrow(p15.norm)))
  colnames(df) <- c("x","y","id")
  return(df)
}

# Load dists Giulia
LoadDistGiulia <- function() {
  filename <- "DistsGiulia/r13_HOOVER.dat"  
  p1 <- read.table(filename)
  # R14
  filename <- "DistsGiulia/r14_HOOVER.dat"     
  p2 <- read.table(filename)
  # Theta
  filename <- "DistsGiulia/bondang_HOOVER.deg.dat"  
  p3 <- read.table(filename)
  # Phi
  filename <- "DistsGiulia/dihang_HOOVER.deg.dat"    
  p4 <- read.table(filename)
  # R15
  filename <- "DistsGiulia/r15_HOOVER.dat"  
  p15 <- read.table(filename)
  p1.norm <- pnorm(p1,"r13")
  p2.norm <- pnorm(p2,"r14")
  p3.norm <- pnorm(p3,"theta")
  p4.norm <- pnorm(p4,"phi")
  p15.norm <- pnorm(p15,"r15")
  df <- data.frame(rbind(p1.norm,p2.norm,p3.norm,p4.norm,p15.norm))
  df$id <- c(rep("r13",nrow(p1.norm)),rep("r14",nrow(p2.norm)),rep("theta",nrow(p3.norm)),rep("phi",nrow(p4.norm)),rep("r15",nrow(p15.norm)))
  colnames(df) <- c("x","y","id")
  return(df)
}

# Load reference dists
LoadWRMS <- function(TestIndex) {
  
  # Position on lpgm-pc is different than t800
  hostname<-Sys.info()["nodename"]
  
  
  if (hostname=="t800") {
    dir <- sprintf("../test%02d",TestIndex)
  } else {
    dir <- sprintf("../Data/test%02d",TestIndex)  
  }
  
  df <- data.frame()
  filename <- paste(dir,'/Param1.wrms.dat',sep="")
  df <- LoadCheck(filename,df,type="r13",ptype="wrms")
  filename <- paste(dir,'/Param2.wrms.dat',sep="")
  df <- LoadCheck(filename,df,type="r14",ptype="wrms")  
  filename <- paste(dir,'/Param3.wrms.dat',sep="")
  df <- LoadCheck(filename,df,type="theta",ptype="wrms")  
  filename <- paste(dir,'/Param4.wrms.dat',sep="")
  df <- LoadCheck(filename,df,type="phi",ptype="wrms")
  
  colnames(df) <- c("wrms","id")
  return(df)
}

# Set specific plot ranges for each term
SetRanges <- function(type) {
  switch(type,
         r13 = c(4.5,6.5),
         r14 = c(4.5,7.0),
         r15 = c(6.0,10.0),
         theta = c(70,110),
         phi = c(-180,180)) 
}

# Set specific plot ranges for each term
SetRangesPMF <- function(type) {
  switch(type,
         r13 = c(5.0,6.0),
         r14 = c(4.5,7.0),
         r15 = c(6.0,10.0),
         theta = c(80,100),
         phi = c(-180,180)) 
}

# Set specific plot ranges for each term
SetRangesGiuliaCheck <- function(type) {
  switch(type,
         r13 = c(3,8),#r13 = c(4.5,6.5),
         r14 = c(2,12),#r14 = c(4.5,7.0),
         r15 = c(6.0,10.0),#r15 = c(6.0,10.0),
         theta = c(60,170),#theta = c(70,110),
         phi = c(-180,180)) #phi = c(30,100)) 
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

# Set Labels
SetLabels <- function(type) {
  switch(type,
         r13 = expression ( paste(r["i,i+2"], " [", ring(A), "]" ) ),
         r14 = expression ( paste(r["i,i+3"], " [", ring(A), "]" ) ),
         r15 = expression ( paste(r["i,i+4"], " [", ring(A), "]" ) ),
         theta = paste(expression(theta),"(Deg)",sep=" "),
         phi = paste(expression(phi),"(Deg)",sep=" "))
}

# Plot dists
PlotDists <- function(df,df.best,df.ref,df.Giu,type,TestIndex) {
    
  r <- SetRanges(type)  
  xlab <- SetLabels(type)
  
  
  # Get some quantiles to create non-linear color scale
  qn = quantile(df$AvgLoss,c(0.01,0.2,0.5,0.7,0.99))
  qn01<-rescale(qn)
  
  # Get max y within restricted range
  xx <- df[df$id==type,]
  maxy <- max( xx[ xx[,1] < r[2] & xx[,1] > r[1], 2] )
  
  plt <- ggplot(data=df[df$id==type,]) + 
    geom_line(aes(x,y,group=run,colour=AvgLoss,alpha=1/AvgLoss),size=0.6) + 
    geom_line(data=df.Giu[df.Giu$id==type,],aes(x,y),size=1.5,color="#FFFF00",linetype="dashed") + 
    geom_line(data=df.ref[df.ref$id==type,],aes(x,y),size=1.5) + 
    geom_line(data=df.best[df.best$id==type,],aes(x,y),size=1.8,color="black") + 
    geom_line(data=df.best[df.best$id==type,],aes(x,y),size=1.5,color="#FF0033") + 
    scale_colour_gradientn(name="Average\nLoss", colours= gray.colors(20), values=c(qn01)) +
    scale_x_continuous(expand=c(0.01,0.01),limits=c(r[1],r[2])) + # remove white spaces left right 
    scale_y_continuous(limits=c(0,maxy)) +
    scale_alpha_continuous(guide=FALSE) +   
    xlab(xlab) +
    ylab("") +    
    thm2 +
    theme(axis.ticks.y = element_blank(),
          axis.text.y  = element_blank())
  #if (fig.format=="PNG") {
    filename <- sprintf("png/Dist_Test%02d_%s.png",TestIndex,type)
    png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
    print(plt)
    dev.off()
  #} else {
    filename <- sprintf("ps/Dist_Test%02d_%s.ps",TestIndex,type)
    cairo_ps(file = filename, width = 6.5, height=3.25, pointsize = 12)
    print(plt)
    dev.off()
  #}

  
#   # DEBUG giulia distr
#   r <- SetRangesGiuliaCheck(type)  
#   plt2 <- ggplot(data=df[df$id==type,]) +     
#     geom_line(data=df.Giu[df.Giu$id==type,],aes(x,y),size=1.5,color="#FF0099",linetype="dashed") + 
#     geom_line(data=df.ref[df.ref$id==type,],aes(x,y),size=1.5) +         
#     scale_x_continuous(expand=c(0.01,0.01),limits=c(r[1],r[2])) + # remove white spaces left right     
#     scale_alpha_continuous(guide=FALSE) +   
#     xlab(xlab) +
#     ylab("") +    
#     thm2 +
#     theme(axis.ticks.y = element_blank(),
#           axis.text.y  = element_blank())
#    filename <- paste("png/Dist_GiuliaCheck",type,".png",sep="")
#    png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
#    print(plt2)
#    dev.off()
  
  return(arrangeGrob(plt))
}
  

# Plot dists
PlotPMFs <- function(df.fitted,df.tobefitted,type,TestIndex) {
  
  r <- SetRangesPMF(type)  
  xlab <- SetLabels(type)
  
  
  # Get some quantiles to create non-linear color scale
  qn = quantile(df.fitted$AvgLoss,c(0.01,0.2,0.5,0.7,0.99))
  qn01<-rescale(qn)
  
  # Get max y within restricted range
  xy <- df.tobefitted[df.tobefitted$id==type,]
  maxy <- max( xy[ xy[,1] < r[2] & xy[,1] > r[1], ]$Run )

  
  # Remove 0 values (can be different in different columns)
  dRun <- df.tobefitted$Run # Current run 
  dSum <- df.tobefitted$Sum # Sum
  dRef <- df.tobefitted$Ref # Ref
  dPrev <- df.tobefitted$Prev # Prev
  dRun[dRun==0]<-NA
  dSum[dSum==0]<-NA
  dRef[dRef==0]<-NA
  dPrev[dPrev==0]<-NA
  df.Run <- data.frame(x=df.tobefitted$x,Run=dRun,id=df.tobefitted$id,run=df.tobefitted$run)
  df.Sum <- data.frame(x=df.tobefitted$x,Sum=dSum,id=df.tobefitted$id,run=df.tobefitted$run)
  df.Ref <- data.frame(x=df.tobefitted$x,Ref=dRef,id=df.tobefitted$id,run=df.tobefitted$run)
  df.Prev <- data.frame(x=df.tobefitted$x,Prev=dPrev,id=df.tobefitted$id,run=df.tobefitted$run)
  
  plt <- list()
  #for (i in c(1:max(df.fitted$run))) {
  for (i in c(1:20)) {
    df1 <- df.fitted[df.fitted$run==i,]
    df2.Run <- df.Run[df.Run$run==i,] #df.tobefitted[df.tobefitted$run==i,]  
    df2.Sum <- df.Sum[df.Sum$run==i,] 
    df2.Ref <- df.Ref[df.Ref$run==i,] 
    df2.Prev <- df.Prev[df.Prev$run==i,] 
    # geom_line(data=df2[df2$id==type,],aes(x,Prev),size=1.0,color=c5[5],linetype=3) +     
    plt[[i]] <- ggplot(data=na.omit(df1[df1$id==type,])) +
      geom_line(aes(x,y),size=1.0,color="black") +
      geom_line(data=na.omit(df2.Run[df2.Run$id==type,]),aes(x,Run),size=1.0,color=c5[1]) +
      geom_line(data=na.omit(df2.Sum[df2.Sum$id==type,]),aes(x,Sum),size=1.0,color="orange") +
      geom_line(data=na.omit(df2.Ref[df2.Ref$id==type,]),aes(x,Ref),size=0.5,color=c5[4],linetype=2) +
      geom_line(data=na.omit(df2.Prev[df2.Prev$id==type,]),aes(x,Prev),size=0.5,color=c5[5],linetype=3) +
      scale_x_continuous(expand=c(0.01,0.01),limits=c(r[1],r[2])) + # remove white spaces left right
      scale_y_continuous(limits=c(0,maxy)) +
      xlab(xlab) +
      ylab("") +
      thm2 +
      theme(axis.ticks.y = element_blank(),
            axis.text.y  = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm")) #+ facet_grid(run~.)
  }
  #if (fig.format=="PNG") {
  filename <- sprintf("png/PMF_Test%02d_%s.png",TestIndex,type)
  png(file = filename, width = 6.5, height=6.5, units = 'in', type = "cairo", res = 600)
  #print(plt)
  do.call("grid.arrange", c(plt, ncol=5))
  dev.off()
  #} else {
#   filename <- sprintf("ps/PMF_Test%02d_%s.ps",TestIndex,type)
#   cairo_ps(file = filename, width = 6.5, height=6.5, pointsize = 12)
#   print(plt)
#   dev.off()
  #}
  
  
  #   # DEBUG giulia distr
  #   r <- SetRangesGiuliaCheck(type)  
  #   plt2 <- ggplot(data=df[df$id==type,]) +     
  #     geom_line(data=df.Giu[df.Giu$id==type,],aes(x,y),size=1.5,color="#FF0099",linetype="dashed") + 
  #     geom_line(data=df.ref[df.ref$id==type,],aes(x,y),size=1.5) +         
  #     scale_x_continuous(expand=c(0.01,0.01),limits=c(r[1],r[2])) + # remove white spaces left right     
  #     scale_alpha_continuous(guide=FALSE) +   
  #     xlab(xlab) +
  #     ylab("") +    
  #     thm2 +
  #     theme(axis.ticks.y = element_blank(),
  #           axis.text.y  = element_blank())
  #    filename <- paste("png/Dist_GiuliaCheck",type,".png",sep="")
  #    png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
  #    print(plt2)
  #    dev.off()
  
  #return(arrangeGrob(plt))
}



LoadPMF.wrap <- function(REP,TestIndex,reload=FALSE) {  
  REP <- ReadData(TestIndex)
  if (reload) {      
    ridx <- c(1:(nrow(REP)-1))
    GlobalIdx <- 0
    for (idx in ridx) {
      PMFs <- LoadPMF(idx,TestIndex)
      tmp <- PMFs[[1]] # Fitted      
      tmp1 <- PMFs[[2]] # ToBeFitted
      tmp$run <- rep(idx,nrow(tmp)) # Add run information
      tmp$AvgLoss <- rep(REP[REP$Run==idx,]$AvgLoss,nrow(tmp)) # Add average loss      
      tmp1$run <- rep(idx,nrow(tmp1)) # Add run information
      tmp1$AvgLoss <- rep(REP[REP$Run==idx,]$AvgLoss,nrow(tmp1)) # Add average loss
      if (GlobalIdx == 0) {    
        df <- tmp
        df1 <- tmp1
      } else {
        df <- rbind(df,tmp)
        df1 <- rbind(df1,tmp1)
      }
      GlobalIdx <- GlobalIdx + 1
    }
    colnames(df) <- c("x","y","id","run","AvgLoss")
    colnames(df1) <- c("x","Run","Prev","Diff","Sum","Ref","id","run","AvgLoss")
    df.fitted <- df
    df.tobefitted <- df1
  
  } else {
    
#     # Load df  
#     load(file="df.ref.Rdata")
#     load(file="df.Giu.Rdata")
#     filename <- sprintf("dfTest%02d.best.Rdata",TestIndex)
#     load(file=filename)
#     filename <- sprintf("dfTest%02d.Rdata",TestIndex)
#     load(file=filename)
  }
  PMFs <- list(df.fitted,df.tobefitted)
  return (PMFs)
}

LoadDists.wrap <- function(REP,TestIndex,reload=FALSE) {  
  REP <- ReadData(TestIndex)
  if (reload) {  
    N <- min((nrow(REP)-2),100) # Include up to 100 interations
    ridx <- sample(max(REP$Run),N)-1 # 0 to N included 
    GlobalIdx <- 0
    for (idx in ridx) {
      tmp <- LoadDist(idx,TestIndex)
      tmp$run <- rep(idx,nrow(tmp)) # Add run information
      tmp$AvgLoss <- rep(REP[REP$Run==idx,]$AvgLoss,nrow(tmp)) # Add average loss
      if (GlobalIdx == 0) {    
        df <- tmp
      } else {
        df <- rbind(df,tmp)
      }
      GlobalIdx <- GlobalIdx + 1
    }
    colnames(df) <- c("x","y","id","run","AvgLoss")
    
    # Best Iteration distributions
    bestr <- REP[which.min(REP$AvgLoss),]$Run
    df.best <- LoadDist(bestr,TestIndex)
    
    
    # Ref distributions for 310
    df.ref <- LoadRefDist()  
    
    # Distributions for 310 Giulia
    df.Giu <- LoadDistGiulia()  
    
    # Store dfs  
    
    save(file="df.ref.Rdata",df.ref)
    save(file="df.Giu.Rdata",df.Giu)
    filename <- sprintf("dfTest%02d.best.Rdata",TestIndex)
    save(file=filename,df.best)
    filename <- sprintf("dfTest%02d.Rdata",TestIndex)
    save(file=filename,df.best)  
    
  } else {
    
    # Load df  
    load(file="df.ref.Rdata")
    load(file="df.Giu.Rdata")
    filename <- sprintf("dfTest%02d.best.Rdata",TestIndex)
    load(file=filename)
    filename <- sprintf("dfTest%02d.Rdata",TestIndex)
    load(file=filename)
  }
  Dists <- list(df,df.best,df.ref,df.Giu)
  return (Dists)
}

PlotLoss <- function(REP,TestIndex) {  
  plt<- ggplot(REP) + 
    geom_line(aes(Run,AvgLoss),size=1.0) +
    xlab("Iteration") +
    ylab("Average Loss") +    
    thm2
  filename <- sprintf("png/Loss_Test%02d.png",TestIndex)
  png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
  print(plt)
  dev.off()
  #} else {
  filename <- sprintf("ps/Loss_Test%02d.ps",TestIndex)
  cairo_ps(file = filename, width = 6.5, height=3.25, pointsize = 12)
  print(plt)
  dev.off()
  
}

## END FUNCTIONS ##


# # Correlation plots for SI (loop from Test 2 to 14)
# fits <- list() # Fitted correlation linear models
# R2s <- list() # Fitted correlation linear models R2
# for (TestIndex in seq(2,14)) {
#   if (TestIndex %in% attributes[attributes$SIM=="MC",1]) { 
#     REP <- ReadData(TestIndex)
#     plots <- PlotCorrelations(REP,TestIndex,fitting=TRUE)    
#     plt <- plots$Corr    
#     fits <- lappend(fits,plots$CorrLinesFit)
#     R2s <- lappend(R2s,plots$R2)
#     filename <- sprintf("png/Correlations_%02d_SI.png",TestIndex)
#     png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
#     print(plt)
#     dev.off()
#     filename <- sprintf("ps/Correlations_%02d_SI.ps",TestIndex)
#     cairo_ps(file = filename, width = 6.5, height=3.25, pointsize = 12)
#     print(plt)
#     dev.off()    
#   }
# }
# 
# i <- 1
# for (TestIndex in seq(2,14)) {
#   fit <- fits[[i]]
#   coef <- coefficients(fit)
#   x <- data.frame(c(1:20))
#   y <- coef[1] + coef[2] * x    
#   
#   if (i==1) {
#     df <- data.frame(x=x,y=y,TestIndex=rep(TestIndex,nrow(x)))  
#   } else {
#     df <- rbind(df, data.frame(x=x,y=y,TestIndex=rep(TestIndex,nrow(x)))  )
#   }  
#   i <- i + 1
# }
# colnames(df) <- c("x","y","TestIndex")
# ggplot(df) + 
#   geom_line(aes(x,y,group=factor(TestIndex),colour=factor(TestIndex))) +
#   scale_color_manual(values=c25) +
#   thm2

# 
# # #############################################################################################
# # Correlation plot Main text. Only for Test = 9
# TestIndex <- 9
# REP <- ReadData(TestIndex)
# plots <- PlotCorrelations(REP,TestIndex,fitting=TRUE)
# plt <- plots$Corr
# filename <- sprintf("png/Correlations_%02d_Main.png",TestIndex)
# png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
# print(plt)
# dev.off()
# filename <- sprintf("ps/Correlations_%02d_Main.ps",TestIndex)
# cairo_ps(file = filename, width = 6.5, height=3.25, pointsize = 12)
# print(plt)
# dev.off()


# # Check ergodicity of MCMC
# TestIndex <- 14
# REP <- ReadData(TestIndex)
# lambda<-0.001
# beta<-1/lambda
# rate<-beta
# AvgLoss.h <- hist(REP$AvgLoss,25)
# AvgLoss.fit.gamma <- fitdistr(REP$AvgLoss,"gamma")
# shape <- AvgLoss.fit.gamma$estimate[1]
# rate <- AvgLoss.fit.gamma$estimate[2]
# dfloss <- data.frame(x.simP=AvgLoss.h$mids, y.simP=AvgLoss.h$density, y.theoP=dgamma(AvgLoss.h$mids,shape=shape,rate=rate))
# plt <- ggplot(dfloss) + geom_point(aes(x=x.simP,y=y.simP),size=3) + 
#   geom_line(aes(x=x.simP,y=y.theoP),size=2,col="red") + 
#   xlab("Loss") + ylab("P(Loss)")  + thm2
# #if (fig.format=="PNG") {
#   filename <- sprintf("png/ErgodicityTest_%02d.png",TestIndex)
#   png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
#   print(plt)
#   dev.off()
# #} else {
#   filename <- sprintf("ps/ErgodicityTest_%02d.ps",TestIndex)
#   cairo_ps(file = filename, width = 6.5, height=3.25, pointsize = 12)
#   print(plt)
#   dev.off()
#}
# 
# 
############################################################################################
# Distribution plot Main text (Test = 13 or Test = 14 I should chose which one)
# Load distributions from a sample of N runs from Test = 13/14
# Set reload TRUE if you want to sample another sample,
# however, it you need to have the raw simulations data in ../test13/ or ../test14/
# 
TestIndex <- 31
REP <- ReadData(TestIndex)
# Load Dists
# Set reload true if you want to read from raw data
# otherwise it will use data frame already stored in ana/df.XXX
# For tests 2 to 14 raw data is only present on t800 so reload, 
# if you are not on t800, set it to false.
reload=TRUE
if (TestIndex<14) {
  reload=FALSE
}
Dists <- LoadDists.wrap(REP,TestIndex,reload=reload)
df <- Dists[[1]]
df.best <- Dists[[2]]
df.ref <- Dists[[3]]
df.Giu <- Dists[[4]]


# Plot distributios
if (TestIndex %in% attributes[attributes$r13==TRUE,1]) {
  PlotDists(df,df.best,df.ref,df.Giu,"r13",TestIndex)
} 
if (TestIndex %in% attributes[attributes$r14==TRUE,1]) {
  PlotDists(df,df.best,df.ref,df.Giu,"r14",TestIndex)
}
if (TestIndex %in% attributes[attributes$theta==TRUE,1]) {
  PlotDists(df,df.best,df.ref,df.Giu,"theta",TestIndex)
}
if (TestIndex %in% attributes[attributes$phi==TRUE,1]) {
  PlotDists(df,df.best,df.ref,df.Giu,"phi",TestIndex)
}

# These are only for IBI 
# Plot PMFs
if (TestIndex %in% attributes[attributes$SIM=="IBI",1]) {
  
  # Load PMFs
  PMFs <- LoadPMF.wrap(REP,TestIndex,reload=TRUE)
  df.fitted <- PMFs[[1]]
  df.tobefitted <- PMFs[[2]]
  
  if (TestIndex %in% attributes[attributes$r13==TRUE,1]) {
    PlotPMFs(df.fitted,df.tobefitted,"r13",TestIndex)  
  }
  if (TestIndex %in% attributes[attributes$r14==TRUE,1]) {
    PlotPMFs(df.fitted,df.tobefitted,"r14",TestIndex)  
  }
  if (TestIndex %in% attributes[attributes$theta==TRUE,1]) {
    PlotPMFs(df.fitted,df.tobefitted,"theta",TestIndex)  
  }
  if (TestIndex %in% attributes[attributes$phi==TRUE,1]) {
    PlotPMFs(df.fitted,df.tobefitted,"phi",TestIndex)  
  }
}

# Plot Loss
REP <- ReadData(TestIndex)
PlotLoss(REP,TestIndex)

