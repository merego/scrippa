library("ggplot2")
thm2<- theme(panel.background = element_rect(fill = 'white'),
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

c5 <- c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")


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
