library(reshape)
library(ggplot2)
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
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Load data
r13<-read.table("r13.dat")
r14<-read.table("r14.dat")
r15<-read.table("r15.dat")
theta<-read.table("theta.dat")
phi<-read.table("phi.dat")

# Convert to data frame
istep <- 20
df <- data.frame(r13)
df$r14<-r14[,2]
df$r15<-r15[,2]
df$theta<-theta[,2]
df$phi<-phi[,2]
df$step <- df[,1]*istep
colnames(df) <- c("frame","r13","r14","r15","theta","phi","step")


# Compute Acf
acf<-acf(df$r13,lag.max=35)
df.acf<-data.frame(acf$lag*istep)
df.acf$r13 <- acf$acf
acf<-acf(df$r14,lag.max=35)
df.acf$r14 <- acf$acf
acf<-acf(df$r15,lag.max=35)
df.acf$r15 <- acf$acf
acf<-acf(df$theta,lag.max=35)
df.acf$theta <- acf$acf
acf<-acf(df$phi,lag.max=35)
df.acf$phi <- acf$acf
colnames(df.acf)[1] <- "step"


fit<-nls(r14~exp(a*step), data=df.acf, start = list(a = 0)) # r14 exponential decay
df.acf.fit <- data.frame(df.acf$step)
df.acf.fit$fit <- predict(fit)
colnames(df.acf.fit) <- c("step","r14.fit")

df.acf.m <- melt.data.frame(df.acf,id.vars="step")
plt <- ggplot(data=df.acf.m) +
  geom_line(aes(x=step,y=value,group=variable,color=variable),size=2) +
  geom_line(data=df.acf.fit,aes(x=step,y=r14.fit),color="red",linetype="dashed",size=1) +
  xlab("step (dump frequency)") +
  ylab("acf") +
  scale_color_manual(values=cbbPalette) +
  ggtitle("310 (17aa), 100000 step, start from folded 310, Giulia FIELD") +
  annotate("text",x=300,y=0.5,label="Dashed red exp fitting on r14. tau ~ 134 step",size=2) +
  thm2

# Save plot to file
filename <- paste("ACF.png",sep='')
png(file = filename, width = 6.5, height=3.25, units = 'in', type = "cairo", res = 600)
print(plt)
dev.off()

#fit<-nls(r13 ~ exp(step * x), data=df, start = list(b = 0))