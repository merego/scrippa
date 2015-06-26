library(reshape)
library(ggplot2)
library(gridExtra)
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


avgloss <- vector()
for (i in seq(1:6)) {
  filename <- paste("r_",i,"/avgloss.dat",sep="")
  tmp <- read.table(filename)
  avgloss <- append(avgloss,tmp[,2])
}

si <- rep(c(1:10),6) # Simulation index
length <- rep(seq(25,200,35),each=10)# Length x 1000

df <- data.frame(avgloss,si,length)

p1 <- ggplot(df,aes(si,avgloss)) + geom_line(aes(group=factor(length),color=factor(length)),size=2) + thm2
p2 <- ggplot(df,aes(x=factor(length),avgloss)) +
  geom_boxplot(aes(fill=factor(length))) + 
  xlab("Simulation steps (x1000)") +
  thm2
grid.arrange(p1,p2,nrow=2)

# Test for unequal variances
bartlett.test(avgloss~length,data=df)
df.2<-df[df$length>95,]
bartlett.test(avgloss~length,data=df.2)