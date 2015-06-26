library("ggplot2")
library("gridExtra")
thm2<- theme(panel.background = element_rect(fill = 'white'),
             axis.ticks.x = element_line(size=0.8, colour="black"),
             axis.title.x =  element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.text.x  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.ticks.y = element_line(size=0.8, colour="black"),
             axis.text.y  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.title.y = element_text(angle=90, vjust=0.5, size=12, face="bold", colour="black"),
             plot.title = element_text(lineheight=3, face="bold", color="black", size=12),
             legend.title  = element_blank(),
             legend.text = element_text(lineheight=3, face="bold", color="black", size=12),
             strip.text = element_text(lineheight=3, face="bold", color="black", size=12)) # Text for facets header

c5 <- c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")

pnorm <- function(p) {
  #p.ks <- locpoly(p[,1],p[,2],gridsize=600,bandwidth = 0.05) 
  #x<-p.ks$x
  #y<-p.ks$y  
  x<-p[,1]
  y<-p[,2]
  dx <- diff(x)[1]
  y[is.na(y)]<-0
  Z <- sum(y)*dx
  pn <- matrix(0,nrow=length(x),ncol=2)
  pn[,1] <- x
  pn[,2] <- y/Z
  return(pn)
}

potEne <- function(k1,k2,theta0r,r0,theta,r) {
 V <- 1/2 * k1 * (r - r0)^2 + 1/2 * k2 * (theta - theta0r)^2
 return(V)
}

theta2r <- function(theta) {
 l <- 3.8
 r <- 2.0 * l * sin(theta/2)
 return(r)
}

deg2rad <- function(deg) {
 rad <- deg/180.0*pi
 return(rad)
}

rad2deg <- function(rad) {
 deg <- rad*180.0/pi
 return(deg)
}

alpha <- function(theta0) {
 alpha <- 3.8 * (cos(theta0/2))
 return(alpha)
}



theta0d <- 88.0 # deg
theta0r <- deg2rad(theta0d) # rad
r0 <- theta2r(theta0r) # A

theta0r_max <- deg2rad(100.0)
rv <- vector()
tv <- vector()
drv <- vector()
dtv <- vector()
drv_Linear <- vector()
for (t in seq(theta0r,theta0r_max,length.out=100)) {
 r <- theta2r(t)  
 rv <- append(rv,r)
 tv <- append(tv,t)
 drv <- append(drv,(r - r0))
 drv_Linear <- append(drv_Linear, ( alpha(theta0r) * (t - theta0r) ) )
 dtv <- append(dtv,(t - theta0r))
}

df <- data.frame(drv=drv,dtv=dtv,drv_Linear=drv_Linear)

plt1 <- ggplot(df) + 
 geom_point(aes(rad2deg(theta0r+dtv),(r0+drv_Linear),colour="A"), size=2.5, col="gray80") + 
 geom_point(aes(rad2deg(theta0r+dtv),(r0+drv_Linear),colour="A"), size=2, col=c5[1]) + 
 geom_line(aes(rad2deg(theta0r+dtv),(r0+drv),colour="B"),col=c5[4], size=1.5) + 
 xlab("Theta [deg]") +
 ylab("r [A]") +
 thm2
plt2 <- ggplot(df) +
  geom_point(aes(rad2deg(dtv),(drv_Linear - drv))) +
  xlab("Delta Theta [deg]") +
  ylab("(r_linear - r) [A]") +
  thm2
filename <- paste("rtheta.png",sep='')
png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
grid.arrange(plt1,plt2,ncol=2)
dev.off()

theta0r_max <- deg2rad(120.0)
gi<-0
TN<-5
KN<-40
Surfaces <- data.frame(matrix(0.0,nrow=(TN*KN*KN),ncol=6))
for (theta in seq(theta0r,theta0r_max,length.out=TN) ) {
 theta <- 1.536
 dt <- theta - theta0r
 r <- theta2r(theta)
 dr <- r - r0
 for ( k1i in seq(1,20,length.out=KN)) {
   for ( k2i in seq(1,20,length.out=KN)) {
     k1 <- k1i + runif(1)*1
     k2 <- k2i + runif(1)*1
     gi <- gi + 1
     V <- potEne(k1,k2,theta0r,r0,theta,r)
     dV <- V -6.71e-7
     Surfaces[gi,] <- c(k1,k2,V,dt,dr,dV)
   }
 }
}

colnames(Surfaces)<-c("k1","k2","V","dt","dr","dV")
colfunc <- colorRampPalette(c5)
plt <- ggplot(Surfaces) + geom_point(aes(k1,k2,colour=abs(dV))) +
  geom_abline(intercept = 115, slope=-(alpha(theta0r)^2), size=1.5, color="black") + 
  scale_colour_gradientn(colours=colfunc(3)) +
  xlab("kr [kcal mol-1 A-2]") + ylab("ktheta [kcal mol-1 rad-2]") +
  thm2
filename <- paste("k1k2.png",sep='')
png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
print(plt)
dev.off()


load("refDistribs.Rdata")
PrefTheta.tmp <- refDistribs[[2]]
t <- PrefTheta.tmp$x
dt <- diff(t)[1]
Pt <- PrefTheta.tmp$y
PrefTheta <- data.frame(t=t,Pt=Pt)
Maxt <- t[which.max(Pt)]

Proba8 <- signif((sum(Pt[t<(Maxt-8)]) + sum(Pt[t>(Maxt+8)])) * dt*100, digits=4)
Proba5 <- signif((sum(Pt[t<(Maxt-5)]) + sum(Pt[t>(Maxt+5)])) * dt*100, digits=4)

plt <- ggplot(PrefTheta) + geom_line(aes(t,Pt),size=1.5) + 
  geom_vline(xintercept=Maxt,col="red") + 
  geom_ribbon(data=subset(PrefTheta,t<(Maxt-5)), aes(x=t,ymax=Pt), ymin=0, fill=c5[4]) + 
  geom_ribbon(data=subset(PrefTheta,t>(Maxt+5)), aes(x=t,ymax=Pt), ymin=0, fill=c5[4]) +
  geom_ribbon(data=subset(PrefTheta,t<(Maxt-8)), aes(x=t,ymax=Pt), ymin=0, fill=c5[5]) + 
  geom_ribbon(data=subset(PrefTheta,t>(Maxt+8)), aes(x=t,ymax=Pt), ymin=0, fill=c5[5]) +
  xlim(75,125) +
  annotate("text", label=paste(Proba8,"%",sep=''), x=104, y=0.02) + annotate("text", label=paste(Proba5,"5%",sep=''), x=101, y=0.04) +
  thm2
filename <- paste("Ptheta.png",sep='')
png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
print(plt)
dev.off()


filename="df.Rdata"
load(filename) 
dfa<-df[c(100:625),]

colfunc <- colorRampPalette(c5)
plt <- ggplot(dfa) +
  geom_point(aes(kr13,ktheta,color=newLoss)) + 
  scale_colour_gradientn(colours=colfunc(3)) +
  geom_abline(intercept = 350, slope=-(alpha(theta0r)^2), size=1.5, color="black") + 
  xlab("kr [kcal mol-1 A-2]") + ylab("ktheta [kcal mol-1 rad-2]") +
  thm2
filename <- paste("k1k2_SIM.png",sep='')
png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
print(plt)
dev.off()

# Read Best one
bestr <- which.min(df$newLoss)
#bestr <- 79
filename <- paste("../../NewData/FullScan/r",bestr,"/RECOMPUTED/Param1.dat",sep="")
p1 <- read.table(filename)
filename <- paste("../../NewData/FullScan/r",bestr,"/RECOMPUTED/Param3.dat",sep="")
p3 <- read.table(filename)

# Ref distributions for 310
# R13
filename <- paste("../../NewData/310/TestCorrelations-1/INPUT/Param1.dat",sep='')    
p1.ref <- read.table(filename)
# Theta
filename <- paste("../../NewData/310/TestCorrelations-1/INPUT/Param3.dat",sep='')    
p3.ref <- read.table(filename)


distr<-data.frame(pnorm(p1),pnorm(p1.ref),pnorm(p3),pnorm(p3.ref))

colnames(distr)<-c("p1.x","p1.y","p1.ref.x","p1.ref.y","p3.x","p3.y","p3.ref.x","p3.ref.y")

plt1 <- ggplot(distr) + 
  geom_line(aes(p1.x,p1.y)) +  
  geom_line(aes(p1.ref.x,p1.ref.y),col="red") +
  xlim(4,6) + 
  xlab("r13") + ylab("p(r13)") +
  ggtitle("red=ref,black=sim") +
  thm2

plt2<- ggplot(distr) + 
  geom_line(aes(rad2deg(p3.x),p3.y)) +  
  geom_line(aes(rad2deg(p3.ref.x),p3.ref.y),col="red") +
  xlim(80,100) + 
  xlab("theta") + ylab("p(theta)") +
  ggtitle("red=ref,black=sim") +
  thm2

filename <- paste("best_SIM.png",sep='')
png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
grid.arrange(plt1,plt2,ncol=2) 
dev.off()


report<-read.table("../../NewData/310/TestCorrelations-1/OUTPUT/report.txt",header=TRUE)
AvgLoss<-read.table("../../NewData/310/TestCorrelations-1/OUTPUT/Loss.dat")
df <- data.frame(as.numeric(report[,3]),as.numeric(report[,7]),AvgLoss[,1])

colnames(df) <- c("kr13","ktheta","Loss")


colfunc <- colorRampPalette(c5)
df.sorted <- sort(df$Loss,index.return=TRUE)
dfmin <- df[df.sorted$ix[1:10],] # Top 10
y <- dfmin[,2]
x <- dfmin[,1]
initslope <- (alpha(theta0r))^2
fit <- nls(y ~ a + b * x, algorithm = "port", start=c(a=120, b=-initslope), upper=c(a=150, b=-(initslope-0.1)), lower=c(a=100, b=-(initslope+0.1)) ) # Fit with 
keff <- summary(fit)$parameters[1]
plt <- ggplot(df) +
  geom_point(aes(kr13,ktheta,color=Loss),size=0.6) + 
  scale_colour_gradientn(colours=colfunc(3)) +
  geom_abline(intercept = keff, slope=-initslope, size=1.5, color="black") + 
  geom_point(data=dfmin,aes(kr13,ktheta),color="black") +
  xlab("kr [kcal mol-1 A-2]") + ylab("ktheta [kcal mol-1 rad-2]") +
  thm2
filename <- paste("k1k2_SIM2.png",sep='')
png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
print(plt)
dev.off()

# Three dists among the top 10
#bestr <- which.min(df$Loss)
thm3<- theme(panel.background = element_rect(fill = 'white'),
             plot.title = element_blank(),
             axis.ticks.x = element_line(size=0.8, colour="black"),
             axis.title.x =  element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.text.x  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.ticks.y = element_line(size=0.8, colour="black"),
             axis.text.y  = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black"),
             axis.title.y = element_text(angle=0, vjust=0.5, size=12, face="bold", colour="black")) 
bestr<-8970
filename <- paste("../../NewData/310/TestCorrelations-1/OUTPUT/r",bestr,"/Param1.dat",sep="")
p1 <- read.table(filename)
filename <- paste("../../NewData/310/TestCorrelations-1/OUTPUT/r",bestr,"/Param3.dat",sep="")
p3 <- read.table(filename)
distr<-data.frame(pnorm(p1),pnorm(p1.ref),pnorm(p3),pnorm(p3.ref))
colnames(distr)<-c("p1.x","p1.y","p1.ref.x","p1.ref.y","p3.x","p3.y","p3.ref.x","p3.ref.y")
plt11 <- ggplot(distr) + 
  geom_line(aes(p1.x,p1.y),size=1) +  
  geom_line(aes(p1.ref.x,p1.ref.y),size=1,col=c5[1]) +
  xlim(4.5,6) + 
  xlab("r13") + ylab("A") +
  ggtitle("red=ref,black=sim") +
  thm3
plt21<- ggplot(distr) + 
  geom_line(aes(rad2deg(p3.x),p3.y),size=1) +  
  geom_line(aes(rad2deg(p3.ref.x),p3.ref.y),size=1,col=c5[1]) +
  xlim(80,100) + 
  xlab("theta") + ylab("") +
  ggtitle("red=ref,black=sim") +
  thm3
bestr<-8788
filename <- paste("../../NewData/310/TestCorrelations-1/OUTPUT/r",bestr,"/Param1.dat",sep="")
p1 <- read.table(filename)
filename <- paste("../../NewData/310/TestCorrelations-1/OUTPUT/r",bestr,"/Param3.dat",sep="")
p3 <- read.table(filename)
distr<-data.frame(pnorm(p1),pnorm(p1.ref),pnorm(p3),pnorm(p3.ref))
colnames(distr)<-c("p1.x","p1.y","p1.ref.x","p1.ref.y","p3.x","p3.y","p3.ref.x","p3.ref.y")
plt12 <- ggplot(distr) + 
  geom_line(aes(p1.x,p1.y),size=1) +  
  geom_line(aes(p1.ref.x,p1.ref.y),size=1,col=c5[1]) +
  xlim(4.5,6) + 
  xlab("r13") + ylab("B") +
  ggtitle("red=ref,black=sim") +
  thm3
plt22<- ggplot(distr) + 
  geom_line(aes(rad2deg(p3.x),p3.y),size=1) +  
  geom_line(aes(rad2deg(p3.ref.x),p3.ref.y),size=1,col=c5[1]) +
  xlim(80,100) + 
  xlab("theta") + ylab("") +
  ggtitle("red=ref,black=sim") +
  thm3 
bestr<-8559
filename <- paste("../../NewData/310/TestCorrelations-1/OUTPUT/r",bestr,"/Param1.dat",sep="")
p1 <- read.table(filename)
filename <- paste("../../NewData/310/TestCorrelations-1/OUTPUT/r",bestr,"/Param3.dat",sep="")
p3 <- read.table(filename)
distr<-data.frame(pnorm(p1),pnorm(p1.ref),pnorm(p3),pnorm(p3.ref))
colnames(distr)<-c("p1.x","p1.y","p1.ref.x","p1.ref.y","p3.x","p3.y","p3.ref.x","p3.ref.y")
plt13 <- ggplot(distr) + 
  geom_line(aes(p1.x,p1.y),size=1) +  
  geom_line(aes(p1.ref.x,p1.ref.y),size=1,col=c5[1]) +
  xlim(4.5,6) + 
  xlab("r13") + ylab("C") +
  ggtitle("red=ref,black=sim") +
  thm3
plt23<- ggplot(distr) + 
  geom_line(aes(rad2deg(p3.x),p3.y),size=1) +  
  geom_line(aes(rad2deg(p3.ref.x),p3.ref.y),size=1,col=c5[1]) +
  xlim(80,100) + 
  xlab("theta") + ylab("") +
  ggtitle("red=ref,black=sim") +
  thm3
filename <- paste("three_BEST.png",sep='')
png(file = filename, width = 6.0, height=6.5, units = 'in', type = "cairo", res = 600)
grid.arrange(plt11,plt21,plt12,plt22,plt13,plt23,ncol=2) 
dev.off()
