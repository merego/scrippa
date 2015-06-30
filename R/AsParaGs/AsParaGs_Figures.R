library("ggplot2")
library("scales")
library("gridExtra")
library("KernSmooth")


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
AsParaGs.version <- "old"
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
    fit <- nls(y ~ a + b * x, algorithm = "port", start=c(a=120, b=-initslope), upper=c(a=150, b=-(initslope-0.1)), lower=c(a=20, b=-(initslope+0.1)) ) # Fit with 
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
    scale_colour_gradientn(name="Average\nLoss", colours= topo.colors(20), values=c(qn01)) +
    scale_x_continuous(expand=c(0.01,0.01)) + # remove white spaces left right
    scale_y_continuous(expand=c(0.01,0.01)) + # remove white spaces bottom top
    thm2
  
  # Add proper xlab, ylab
  ylabel <- expression ( paste ( k[theta], " [ kcal ", mol^"-1" , rad^"-1", "]" , sep = " ")   )
  if (TestIndex == 7 || TestIndex == 8 || TestIndex == 12 || TestIndex == 13) {    
    xlabel <- expression ( paste ( epsilon[r["i,i+2"]], " [ kcal ", mol^"-1" , Å^"-1", "]" , sep = " ")   )
  } else {
    xlabel <- expression ( paste ( k[r["i,i+2"]], " [ kcal ", mol^"-1" , Å^"-1", "]" , sep = " ")   )    
  }
  plt <- plt + xlab(xlabel) + ylab(ylabel)    
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
PlotDists <- function(df,df.best,df.ref,df.Giu,type) {
    
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
    filename <- paste("png/Dist_",type,".png",sep="")
    png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
    print(plt)
    dev.off()
  #} else {
    filename <- paste("ps/Dist_",type,".ps",sep="")
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
  

## END FUNCTIONS ##


# Correlation plots for SI
for (TestIndex in seq(2,13)) {
 REP <- ReadData(TestIndex)
 Surface <- data.frame(REP$Param1,REP$Param2,REP$Param1.1,REP$Param2.1,REP$AvgLoss)
 colnames(Surface)<-c("kr","r0","ktheta","theta0","AvgLoss")
 plt <- PlotCorrelations(Surface,TestIndex,fitting=TRUE)
 #if (fig.format=="PNG") {
   filename <- sprintf("png/Correlations_%02d_SI.png",TestIndex)
   png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
   print(plt)
   dev.off()
 #} else {
   filename <- sprintf("ps/Correlations_%02d_SI.ps",TestIndex)
   cairo_ps(file = filename, width = 6.5, height=3.25, pointsize = 12)
   print(plt)
   dev.off()
 #}
}

#############################################################################################
# Correlation plot Main text (Test = 9)
TestIndex <- 9
REP <- ReadData(TestIndex)
Surface <- data.frame(REP$Param1,REP$Param2,REP$Param1.1,REP$Param2.1,REP$AvgLoss)
colnames(Surface)<-c("kr","r0","ktheta","theta0","AvgLoss")
plt <- PlotCorrelations(Surface,TestIndex,fitting=TRUE)
#if (fig.format=="PNG") {
  filename <- sprintf("png/Correlations_%02d_Main.png",TestIndex)
  png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
  print(plt)
  dev.off()
#} else  {
  filename <- sprintf("ps/Correlations_%02d_Main.ps",TestIndex)
  cairo_ps(file = filename, width = 6.5, height=3.25, pointsize = 12)
  print(plt)
  dev.off()
#}

# Check ergodicity of MCMC
TestIndex <- 13
REP <- ReadData(TestIndex)
lambda<-0.08
beta<-1/lambda
fith<-hist(REP$AvgLoss,25)
dx <- diff(fith$mids)[1]
theoP <- exp(-beta * fith$mids)/dx
dfloss <- data.frame(x.simP=fith$mids, y.simP=fith$density, y.theoP=theoP)
plt <- ggplot(dfloss) + geom_line(aes(x=x.simP,y=y.simP),size=2) + 
  geom_line(aes(x=x.simP,y=y.theoP),size=2,col="red") + 
  xlab("Loss") + ylab("P(Loss)")  + thm2
#if (fig.format=="PNG") {
  filename <- sprintf("png/ErgodicityTest_%02d.png",TestIndex)
  png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
  print(plt)
  dev.off()
#} else {
  filename <- sprintf("ps/ErgodicityTest_%02d.ps",TestIndex)
  cairo_ps(file = filename, width = 6.5, height=3.25, pointsize = 12)
  print(plt)
  dev.off()
#}


############################################################################################
# Distribution plot Main text (Test = 13)
# Load distributions from a sample of N runs from Test = 13
# Set reload TRUE if you want to sample another sample,
# however, it you need to have the raw simulations data in ../test13/ 
# 
reload <- FALSE
TestIndex <- 13
REP <- ReadData(TestIndex)
if (reload) {
  ridx <- sample(max(REP$Run),100)-1 # 0 to N included 
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
  save(file="dfTest13.Rdata",df)
  save(file="df.ref.Rdata",df.ref)
  save(file="df.Giu.Rdata",df.Giu)
  save(file="dfTest13.best.Rdata",df.best)
  
} else {
  
  # Load df
  load(file="dfTest13.Rdata")  
  load(file="df.ref.Rdata")
  load(file="df.Giu.Rdata")
  load(file="dfTest13.best.Rdata")
}



# Plot distributions
p13 <- PlotDists(df,df.best,df.ref,df.Giu,"r13")
p14 <- PlotDists(df,df.best,df.ref,df.Giu,"r14")
p15 <- PlotDists(df,df.best,df.ref,df.Giu,"r15")
ptheta <- PlotDists(df,df.best,df.ref,df.Giu,"theta")
pphi <- PlotDists(df,df.best,df.ref,df.Giu,"phi")










