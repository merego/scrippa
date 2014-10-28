lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
}

OddRatio <- function (mat,conf.level=0.95,lnout=TRUE) {
# Odd ratios
# Matrix should be col 1 success, col 2 faliure 
  odd1 <- mat[1,1]
  results <- list(OR=-1,conf.lev.l=-1,conf.lev.u=-1)
  odd1 <- mat[1,1]/mat[1,2]
  odd2 <- mat[2,1]/mat[2,2]
  if(odd1<=1e-8)
    odd1 <- 1e-4
  if(odd2<=1e-8)
    odd2 <- 1e-4
  OR <- odd1/odd2
  logOR <- log(OR)
  results$OR <- logOR  
  z <- qnorm(1-((1-conf.level)/2))
  conf.lev.l <- logOR - z*sqrt(1/(mat[1,1]+1e-1) + 1/(mat[1,2]+1e-1) + 1/(mat[2,1]+1e-1) + 1/(mat[2,2]+1e-1))          
  conf.lev.u <- logOR + z*sqrt(1/(mat[1,1]+1e-1) + 1/(mat[1,2]+1e-1) + 1/(mat[2,1]+1e-1) + 1/(mat[2,2]+1e-1))            
  if (!lnout) {
    conf.lev.l <- exp(conf.lev.l)
    conf.lev.u <- exp(conf.lev.u)      
    results$OR <- OR  
  }
  
  results$conf.lev.l <- conf.lev.l
  results$conf.lev.u <- conf.lev.u
    
  return (results)
}


Cycles <- 3
aalist<-list("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
CountsTop <- list()
CountsWorse <- list()



cyc <- 1
ORCycle <- list()
while (cyc <= Cycles) {
  fnameTop <- paste("../Cycle_",cyc,"/Docking/Top100WebLogo.txt",sep='')
  fnameWorse <- paste("../Cycle_",cyc,"/Docking/Worse100WebLogo.txt",sep='')
  CountTop <- read.table(fnameTop)
  CountTop <- CountTop[,c(2:(length(aalist)+1))]  
  colnames(CountTop)<-aalist
  CountWorse <- read.table(fnameWorse)
  CountWorse <- CountWorse[,c(2:(length(aalist)+1))]  
  colnames(CountWorse)<-aalist
  CountsTop <- lappend(CountsTop,CountTop)
  CountsWorse <- lappend(CountsWorse,CountWorse)
  
  ContTable <- matrix(0,nrow=2,ncol=2)
  rownames(ContTable) <- c("Top","Worse")
  colnames(ContTable) <- c("Present","Absent")
  
  N <- sum(CountTop[1,])
  
  ORPositions <- list()
  for (position in 1:nrow(CountTop)) {
    ORaa <- list()
    for (aaindex in 1:ncol(CountTop)) {
      ContTable[1,1] <- CountTop[position,aaindex]
      ContTable[1,2] <- N - ContTable[1,1]
      ContTable[2,1] <- CountWorse[position,aaindex]
      ContTable[2,2] <- N - ContTable[2,1]
      
      OddR <- OddRatio(ContTable,conf.level=0.95,lnout=TRUE)
      
      if (is.infinite(OddR$OR)) {
        print(ContTable)
      }
      ORaa <- lappend(ORaa,OddR)
    }
    ORPositions <- lappend(ORPositions,ORaa)
  }    
  ORCycle <- lappend(ORCycle,ORPositions)
  
  cyc <- cyc + 1  
}

# Reshape as matrix to plot
# Different plot for each position
rgb.palette <- colorRampPalette(c("red", "green", "blue"), space = "rgb")
colors<-rgb.palette(length(aalist))
par(mfrow=c(5,4),mar=c(0.0,2.5,5,0),mgp=c(1,0.4,0))
cyclev<-c(1:Cycles)
library("ggplot2")
for (position in 1:nrow(CountTop)) {
  FirstPlot<-TRUE
  ORMatrix <- matrix(0,nrow=Cycles,ncol=length(aalist))
  ORMatrixCIL <- matrix(0,nrow=Cycles,ncol=length(aalist)) # Lower confidence interval
  ORMatrixCIU <- matrix(0,nrow=Cycles,ncol=length(aalist)) # Upper confidence interval
  for (cyc in 1:Cycles) {
    for (aaindex in 1:length(aalist)) {
      ORMatrix[cyc,aaindex] <- ORCycle[[cyc]][[position]][[aaindex]]$OR
      ORMatrixCIL[cyc,aaindex] <- ORCycle[[cyc]][[position]][[aaindex]]$conf.lev.l
      ORMatrixCIU[cyc,aaindex] <- ORCycle[[cyc]][[position]][[aaindex]]$conf.lev.u
    }
  }
  #if (FirstPlot) {
  #  plot(cyclev,ORMatrix[,1],'l',lwd=2.0,col=1,xlim=c(1,Cycles),ylim=range(ORMatrix))
  #  FirstPlot<-FALSE
  #}
  gplotM<-data.frame("Cycle"=c(1:Cycles))
  plist <- list()
  for (aaindex in 2:length(aalist)) {
    gplotM$OR <- ORMatrix[,aaindex]
    gplotM$ORCIL <- ORMatrixCIL[,aaindex]
    gplotM$ORCIU <- ORMatrixCIU[,aaindex]
    pl <- ggplot(gplotM, aes(Cycle)) + 
      geom_line(aes(y=OR), colour="blue", size=2) + 
      geom_ribbon(aes(ymin=ORCIL, ymax=ORCIU), alpha=0.2) + 
      ggtitle(as.character(aalist[aaindex])) +
      scale_x_discrete(labels=c(1:Cycles)) +
      theme(axis.title.x = element_text(face="bold", colour="black", size=20),
            axis.text.x  = element_text(angle=0, vjust=0.5, size=30, colour="black"),
            axis.text.y  = element_text(angle=90, vjust=0.5, size=30, colour="black"),
            axis.title.y = element_blank(),            
            plot.title = element_text(lineheight=3, face="bold", color="black", size=30))
    
    plist <- lappend(plist,pl)
    #plot(cyclev,ORMatrix[,aaindex],'l',lwd=2.0,col=1,xlim=c(1,Cycles),ylab="",xlab="",main=aalist[aaindex])        
    #lines(cyclev,ORMatrix[,aaindex],lwd=2.0,col=1)
  }
  filename <- paste("OddRationPosition_",position,".jpg",sep="")
  jpeg(filename,width=2000,height=2000,quality=100)
  multiplot(plotlist=plist,cols=4)
  dev.off()
}
#title("Log odd ratios", line = -1, outer = TRUE)

