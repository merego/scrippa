#!/usr/bin/R

range01 <- function(x){(x-min(x))/(max(x)-min(x))}


matframe <- as.data.frame(read.table('meat.dat'))
#matframe <- matframe[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
#names <- c("r13 epsi","r13 r0","r13 alpha","r14 epsi","r14 r0","r14 alpha","r15 epsi","r15 r0","r15 alpha","KL")
names <- c("r13 k","r13 r0","r13 loss","r14 k","r14 r0","r14 loss","r15 k","r15 r0","r15 loss","angle k","angle r0","angle loss","dih k","dih r0","dih loss","Average Loss")
colnames(matframe) <- names

standardised <- as.data.frame(scale(t(matframe[,1:ncol(matframe)])))
standardised <- t(standardised)

#cat("Doing PLSR and PLSDA\n")
#y <- matframe[,ncol(matframe)]
#nPLScomp = 3
#library('pls')
#plsrmodel <- plsr(y ~ . , data=data.frame(standardised), ncomp=nPLScomp, validation="LOO")


library('akima')
#par(mfrow=c(3,2),mai=c(0.2,0.2,0.2,0.2),mgp=c(0.5,1.0,0.0))

bin <- 25

# r13
filename <- paste("../imgs/","r13SURF.eps",sep='')
postscript(filename)
xyz <- matframe[,c(1,2,3)]
r13sort <- sort(xyz[,3],index.return=TRUE)
s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = bin), yo = seq(min(xyz[,2]), max(xyz[,2]), length = bin), linear = TRUE)
s$z[is.na(s$z)] <- 1.0
image(s$x, s$y, s$z/max(s$z,na.rm=TRUE), xlab=names[10], ylab=names[11],zlim=c(0.0,0.5))
contour(s$x, s$y, s$z/max(s$z,na.rm=TRUE), add=TRUE, nlevel=10, lwd=0.3)
dev.off()
# r14
filename <- paste("../imgs/","r14SURF.eps",sep='')
postscript(filename)
xyz <- matframe[,c(4,5,6)]
r14sort <- sort(xyz[,3],index.return=TRUE)
s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = bin), yo = seq(min(xyz[,2]), max(xyz[,2]), length = bin), linear = TRUE)
s$z[is.na(s$z)] <- 1.0
image(s$x, s$y, s$z/max(s$z,na.rm=TRUE), xlab=names[4], ylab=names[5])
contour(s$x, s$y, s$z/max(s$z,na.rm=TRUE), add=TRUE, nlevel=5, lwd=0.3)
dev.off()
# r15
filename <- paste("../imgs/","r15SURF.eps",sep='')
postscript(filename)
xyz <- matframe[,c(7,8,9)]
r15sort <- sort(xyz[,3],index.return=TRUE)
s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = bin), yo = seq(min(xyz[,2]), max(xyz[,2]), length = bin), linear = TRUE)
s$z[is.na(s$z)] <- 1.0
image(s$x, s$y, s$z/max(s$z,na.rm=TRUE), xlab=names[7], ylab=names[8])
contour(s$x, s$y, s$z/max(s$z,na.rm=TRUE), add=TRUE, nlevel=5, lwd=0.3)
dev.off()
# angle
filename <- paste("../imgs/","angleSURF.eps",sep='')
postscript(filename)
xyz <- matframe[,c(10,11,12)]
Anglesort <- sort(xyz[,3],index.return=TRUE)
xyz[,2] <- xyz[,2]*180/3.14
s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = bin), yo = seq(min(xyz[,2]), max(xyz[,2]), length = bin), linear = TRUE)
s$z[is.na(s$z)] <- 1.0
image(s$x, s$y, s$z/max(s$z,na.rm=TRUE), xlab=names[10], ylab=names[11],zlim=c(0.0,0.5))
contour(s$x, s$y, s$z/max(s$z,na.rm=TRUE), add=TRUE, nlevel=10, lwd=0.3)
dev.off()
# dihedral
filename <- paste("../imgs/","dihedralSURF.eps",sep='')
postscript(filename)
xyz <- matframe[,c(13,14,15)]
Dihsort <- sort(xyz[,3],index.return=TRUE)
xyz[,2] <- (xyz[,2]+3.14)*180/3.14
s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = bin), yo = seq(min(xyz[,2]), max(xyz[,2]), length = bin), linear = TRUE)
s$z[is.na(s$z)] <- 1.0
image(s$x, s$y, s$z/max(s$z,na.rm=TRUE), xlab=names[13], ylab=names[14])
contour(s$x, s$y, s$z/max(s$z,na.rm=TRUE), add=TRUE, nlevel=5, lwd=0.3)
dev.off()


# Best params based on average loss function
avgSort<-sort(matframe[,ncol(matframe)],index.return=TRUE)
avgSort$ix 
# Save out latex table
library(xtable)
s <- matframe[avgSort$ix[1:10],c(1,2,4,5,7,8,10,11,13,14)]
s[,10] <- (s[,10]+3.14)*180/3.14
s[,8] <- s[,8]*180/3.14
all_mean <- sapply(s,mean)
all_std <- sapply(s,sd)
s <- rbind(s,all_mean)
s <- rbind(s,all_std)
rownames(s)[nrow(s)-1]<-c("Mean")
rownames(s)[nrow(s)]<-c("StdDev")
tab <- xtable(s, caption= "Best parameters based on average loss function. Params are sorted on average loss function and the top 10 are selected", align=c("c","c","c","c","c","c","c","c","c","c","c"))
#print(tab,file="BestBasedOnAvgLoss.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(tab),by=1))
print(tab,file="BestBasedOnAvgLoss.tex",append=T,table.placement = "h", caption.placement="bottom")

# Best params based on common sorted data
common <- Reduce(intersect,list(r13sort$ix[1:100],r14sort$ix[1:100],r15sort$ix[1:100],Anglesort$ix[1:100],Dihsort$ix[1:100]))
avgSort<-sort(matframe[,ncol(matframe)],index.return=TRUE)
avgSort$ix 
# Save out latex table
s <- matframe[common,c(1,2,4,5,7,8,10,11,13,14)]
s[,10] <- (s[,10]+3.14)*180/3.14
s[,8] <- s[,8]*180/3.14
all_mean <- sapply(s,mean)
all_std <- sapply(s,sd)
s <- rbind(s,all_mean)
s <- rbind(s,all_std)
rownames(s)[nrow(s)-1]<-c("Mean")
rownames(s)[nrow(s)]<-c("StdDev")
tab <- xtable(s, caption= "Best parameters based on best common loss function. Each potential is sorted based on its own loss function, 
then the iterations which are common among the top 50s are selected. The first 10 are shown here.", align=c("c","c","c","c","c","c","c","c","c","c","c"))
#print(tab,file="BestBasedOnAvgLoss.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(tab),by=1))
print(tab,file="BestBasedOnBestCommonLoss.tex",append=T,table.placement = "h", caption.placement="bottom")

lossRaw <- matframe[,c(3,6,9,12,15)]
loss <- apply(lossRaw,2,range01) 

SumSqLoss<-rowSums(loss^2)
# Sort more similar to each other (closer to average)
Sphericity <- rowSums((loss-rowMeans(loss))^2)
TopS <- sort(Sphericity,index.return=TRUE)
Npc <- TopS$ix[1:ceiling(nrow(loss)*0.2)]
IsotropicSSL <- SumSqLoss[Npc]
# Sort them by lowest SumSqLoss
IsotropicSSLSorted <- sort(IsotropicSSL,index.return=TRUE)
TopTen <- Npc[IsotropicSSLSorted$ix[1:10]]
# Sort sum of square of loss 
#SumSqLossSort<-sort(SumSqLoss,index.return=TRUE)
# Among the top 10%  select the one which are more similar each other (closer to average value)
#SumSqLossSort10pc <- SumSqLossSort$ix[1:ceiling(nrow(loss)*0.1)]
#TopS <- sort(rowSums(loss[SumSqLossSort10pc,]-rowMeans(loss[SumSqLossSort10pc,])^2),index.return=TRUE)
#Top  <- SumSqLossSort10pc[TopS$ix]


# Best based on mi spider area and more isotropic
# Convert to polar coordinates
Nloss <- ncol(loss)
theta <- seq(0, (2*pi - 2*pi/Nloss), 2*pi/(Nloss))
# And back to cartesian coordinates
x <- as.matrix(loss[,]) %*% diag(as.vector(cos(theta)))
y <- as.matrix(loss[,]) %*% diag(as.vector(sin(theta)))
#plot(t(x[TopTen,]),t(y[TopTen,]),col=1:10)
#for (i in 1:10) {
# polygon(x[TopTen[i],],y[TopTen[i],],col=i,density=0.0)
#}
library("splancs")
Area <- vector();
for (i in 1:nrow(loss[Npc,])) {
  m <- cbind(x[i,],y[i,])
  Area <- rbind(Area,areapl(m))
}
AreaSort <- sort(Area,index.return=TRUE)
TopTen <- AreaSort$ix[1:10]
plot(t(x[TopTen,]),t(y[TopTen,]),col=1:10)
for (i in 1:10) {
 polygon(x[TopTen[i],],y[TopTen[i],],col=i,density=0.0)
}
# Save out latex table
s <- matframe[TopTen,c(1,2,4,5,7,8,10,11,13,14)]
s[,10] <- (s[,10]+3.14)*180/3.14
s[,8] <- s[,8]*180/3.14
all_mean <- sapply(s,mean)
all_std <- sapply(s,sd)
s <- rbind(s,all_mean)
s <- rbind(s,all_std)
rownames(s)[nrow(s)-1]<-c("Mean")
rownames(s)[nrow(s)]<-c("StdDev")
tab <- xtable(s, caption= "Best parameters ...", align=c("c","c","c","c","c","c","c","c","c","c","c"))
#print(tab,file="BestBasedOnAvgLoss.tex",append=T,table.placement = "h", caption.placement="bottom", hline.after=seq(from=-1,to=nrow(tab),by=1))
print(tab,file="BestBasedOnMinCircArea.tex",append=T,table.placement = "h", caption.placement="bottom")




# Radar plot
filename <- paste("../imgs/","RadarPlot.eps",sep='')
postscript(filename)
radarchart(as.data.frame(loss[TopTen,]), axistype=2, seg=3, plty=1,maxmin=FALSE,plwd=2.0)
dev.off()

# Super plot
#filename <- paste("./imgs/","ParametersKDplots.ps",sep='')
#postscript(filename);
#N <- length(matframe)-1
#Nplots <- N*(N-1)/2
#Nxplots <- ceiling(sqrt(Nplots))
#Nyplots <- ceiling(Nplots/Nxplots)
#par(mfrow=c(Nxplots,Nyplots),mai=c(0.2,0.2,0.2,0.2),mgp=c(0.5,1.0,0.0))
##require(KernSmooth)
#library('akima')
#for (i in 1:(N-1)) {
# for (j in (i+1):N) {
#  xyz <- matframe[,c(i,j,(length(matframe)))]
#  s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = 50), yo = seq(min(xyz[,2]), max(xyz[,2]), length = 50), linear = FALSE)
#  s <- interp(xyz[,1],xyz[,2],xyz[,3],  xo = seq(min(xyz[,1]), max(xyz[,1]), length = 50), yo = seq(min(xyz[,2]), max(xyz[,2]), length = 50), linear = TRUE)
#  
# # Kernl smoothed 2D-density 
# #z <- bkde2D(xy, .5)
# #persp(z$fhat)
# #surface3d(z$x1,z$x2,z$fhat*10,col="red")
#
# image(s$x, s$y, s$z, xlab=names[i], ylab=names[j])
# 
# contour(s$x, s$y, s$z, add=TRUE, nlevel=5, lwd=0.1)
# }
#}
#dev.off()


