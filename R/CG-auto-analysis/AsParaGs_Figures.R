library("ggplot2")
library("scales")

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
             legend.title  = element_blank(),
             legend.text = element_text(lineheight=3, face="bold", color="black", size=12),
             strip.text = element_text(lineheight=3, face="bold", color="black", size=12)) # Text for facets header

# Color palette
c5 <- c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")
c10 <- c("#a6cee3",  "#1f78b4",   "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Normalize distribution
pnorm <- function(p,type) {
  #p.ks <- locpoly(p[,1],p[,2],gridsize=600,bandwidth = 0.05) 
  #x<-p.ks$x
  #y<-p.ks$y  
  x<-p[,1] 
  x <- ConvertUnit(x,type)
  y<-p[,2]
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

# Read report_XX.dat and  AvgLossXX.dat
# XX == TestIndex
ReadData <- function(TestIndex) {
 filename <- sprintf("report_%02d.dat",TestIndex)
 tmp<-read.table(filename,header=TRUE)
 filename <- sprintf("AvgLoss%02d.dat",TestIndex)
 los<-read.table(filename)[,1]
 REP <- cbind(tmp,los)
 rm(tmp)
 rm(los)
 return(REP)
}

# Main correlation plot 
PlotCorrelations <- function(Surface,TestIndex,fitting=FALSE) {
  
  # Fit intercept of theoretical line
  best.loss <- sort(Surface$AvgLoss,index.return=TRUE)
  dfmin <- Surface[best.loss$ix[1:10],] # Top 10
  if (fitting) {
    y <- dfmin[,3]
    x <- dfmin[,1]
    initslope <- (alpha(theta0r))^2 # Theoretical slope
    fit <- nls(y ~ a + b * x, algorithm = "port", start=c(a=120, b=-initslope), upper=c(a=150, b=-(initslope-0.1)), lower=c(a=100, b=-(initslope+0.1)) ) # Fit with 
    keff <- summary(fit)$parameters[1]
    slope <- initslope
  } else {
    keff <- 115
    slope <- initslope
  }
  
  # Get some quantiles to create non-linear color scale
  qn = quantile(Surface$AvgLoss,c(0.01,0.2,0.5,0.7,0.99))
  qn01<-rescale(qn)
  
  # plot
  plt <- ggplot(Surface) + geom_point(aes(kr,ktheta,colour=AvgLoss),size=0.5) +
    geom_abline(intercept = keff, slope=-(initslope), size=0.8, color="black") +
    geom_point(data=dfmin,aes(kr,ktheta),colour="black",size=1) +
    scale_colour_gradientn(colours= topo.colors(20), values=c(qn01)) +
    scale_x_continuous(expand=c(0.01,0.01)) + # remove white spaces left right
    scale_y_continuous(expand=c(0.01,0.01)) + # remove white spaces bottom top
    thm2
  
  # Add proper xlab, ylab
  if (TestIndex < 7 || TestIndex == 9) {
    plt <- plt + xlab("kr [kcal mol-1 A-2]") + ylab("ktheta [kcal mol-1 rad-2]")
  } else {
    plt <- plt + xlab("epsi [kcal mol-1 A-2]") + ylab("ktheta [kcal mol-1 rad-2]")
  }
  return(plt)
}

# Load dists 
LoadDist <- function(idx,TestIndex) {
  dir <- sprintf("../test%02d/OUTPUT/r",TestIndex)
  filename <- paste(dir,idx,'/Param1.dat',sep="")
  p1 <- as.data.frame(read.table(filename))
  filename <- paste(dir,idx,'/Param2.dat',sep="")
  p2 <- as.data.frame(read.table(filename))
  filename <- paste(dir,idx,'/Param3.dat',sep="")
  p3 <- as.data.frame(read.table(filename))
  filename <- paste(dir,idx,'/Param4.dat',sep="")
  p4 <- as.data.frame(read.table(filename))
  filename <- paste(dir,idx,'/Param15.dat',sep="")
  p15 <- as.data.frame(read.table(filename))
  df <- data.frame(rbind(pnorm(p1,"r13"),pnorm(p2,"r14"),pnorm(p3,"theta"),pnorm(p4,"phi"),pnorm(p15,"r15")))
  df$id <- c(rep("r13",nrow(p1)),rep("r14",nrow(p2)),rep("theta",nrow(p3)),rep("phi",nrow(p4)),rep("r15",nrow(p15)))
  colnames(df) <- c("x","y","id")
  return(df)
}


# Load reference dists
LoadRefDist <- function(idx) {
  filename <- paste("../../TestCorrelations-1/INPUT/Param1.dat",sep='')    
  p1 <- read.table(filename)
  # R14
  filename <- paste("../../TestCorrelations-1/INPUT/Param2.dat",sep='')    
  p2 <- read.table(filename)
  # Theta
  filename <- paste("../../TestCorrelations-1/INPUT/Param3.dat",sep='')    
  p3 <- read.table(filename)
  # Phi
  filename <- paste("../../TestCorrelations-1/INPUT/Param4.dat",sep='')    
  p4 <- read.table(filename)
  # R15
  filename <- paste("../../TestCorrelations-1/INPUT/Param15.dat",sep='')    
  p5 <- read.table(filename)
  df <- data.frame(rbind(pnorm(p1,"r13"),pnorm(p2,"r14"),pnorm(p3,"theta"),pnorm(p4,"phi"),pnorm(p15,"r15")))
  df$id <- c(rep("r13",nrow(p1)),rep("r14",nrow(p2)),rep("theta",nrow(p3)),rep("phi",nrow(p4)),rep("r15",nrow(p15)))
  colnames(df) <- c("x","y","id")
  return(df)
}


# Set specific plot ranges for each term
SetRanges <- function(type) {
  switch(type,
         r13 = c(4.5,6.5),
         r14 = c(4.0,6.5),
         r15 = c(5.0,8.0),
         theta = c(70,110),
         phi = c(20,80)) 
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
         r14 = expression ( paste(r["i,i+2"], " [", ring(A), "]" ) ),
         r15 = expression ( paste(r["i,i+2"], " [", ring(A), "]" ) ),
         theta = paste(expression(theta),"(Deg)",sep=" "),
         phi = paste(expression(phi),"(Deg)",sep=" "))
}

# Plot dists
PlotDists <- function(df,type) {
  r <- SetRanges(type)  
  xlab <- SetLabels(type)
  ggplot(data=df[df$id==type,]) + 
    geom_line(aes(x,y,group=run,colour=AvgLoss),size=1) + 
    geom_line(data=df.ref[df.ref$id==type,],aes(x,y),size=2) + 
    geom_line(data=df.best[df.best$id==type,],aes(x,y),size=2,color="black") + 
    geom_line(data=df.best[df.best$id==type,],aes(x,y),size=1.8,color="#ff5eff") + 
    xlim(r[1],r[2]) + 
    xlab(xlab) +
    ylab("") +
    thm2
}
  

##########

# Correlation plots for SI
for (TestIndex in seq(2,10)) {
  REP <- ReadData(TestIndex)
  Surface <- data.frame(REP$Param1,REP$Param2,REP$Param1.1,REP$Param2.1,REP$los)
  colnames(Surface)<-c("kr","r0","ktheta","theta0","AvgLoss")
  plt <- PlotCorrelations(Surface,TestIndex,fitting=TRUE)
  filename <- sprintf("Correlations_%02d_SI.png",TestIndex)
  png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
  print(plt)
  dev.off() 
}

# Correlation plot Main text (Test = 9)
TestIndex <- 9
REP <- ReadData(TestIndex)
Surface <- data.frame(REP$Param1,REP$Param2,REP$Param1.1,REP$Param2.1,REP$los)
colnames(Surface)<-c("kr","r0","ktheta","theta0","AvgLoss")
plt <- PlotCorrelations(Surface,TestIndex,fitting=TRUE)
filename <- sprintf("Correlations_%02d_Main.png",TestIndex)
png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
print(plt)
dev.off() 


# Load N randomly sampled distributions
ridx <- sample(10,8)-1 # 0 to N included 
GlobalIdx <- 0
for (idx in ridx) {
  tmp <- LoadDist(idx,TestIndex)
  tmp$run <- rep(idx,nrow(tmp)) # Add run information
  tmp$AvgLoss <- rep(Surface$AvgLoss[idx+1],nrow(tmp)) # Add average loss
  if (GlobalIdx == 0) {    
    df <- tmp
  } else {
    df <- rbind(df,tmp)
  }
  GlobalIdx <- GlobalIdx + 1
}
colnames(df) <- c("x","y","id","run","AvgLoss")

# Best Iteration distributions
bestr <- which.min(Surface$AvgLoss)
bestr <- 5
TestIndex <- 9
df.best <- LoadDist(bestr,TestIndex)


# Ref distributions for 310
#df.ref <- LoadRefDist()
df.ref <- df.best

# Plot distributions
PlotDists(df,"r13")
