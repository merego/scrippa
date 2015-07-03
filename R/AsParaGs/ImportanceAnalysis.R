###################################################
# Feature selection 
# 
library("caret")


## FUNCTIONS ##

# Read reportXX.dat 
# XX == TestIndex
ReadData <- function(TestIndex) {
  filename <- sprintf("report%02d.dat",TestIndex)
  REP <- read.table(filename,header=TRUE)
  return(REP)
}

## END FUNCTIONS ##

# Load report
TestIndex <- 13
REP <- ReadData(TestIndex)

# Select a random subsample of dimension (max) 1000.
N <- min(nrow(REP),1000)
ridx <- sample(max(REP$Run),N) 
Col.Names <- colnames(REP)

# Get the indexes of columns which contains parameters
Col.index <- grep("Param",Col.Names) 

# Parameter matrix (indep. variables)
XX <- data.frame(REP[ridx,Col.index])

# Average loss (dependent variable)
YY <- REP$AvgLoss[ridx]

# Find out correspondence paramter potential term.
# The vector pterm indicates the at which potential term belong 
# a given parameter.
pterm <- vector()
pot.index <- grep("Potential",Col.Names) 
res.index <- grep("Result",Col.Names) 
for (iterm in c(1:length(pot.index)) ) {
  tmp <- c((pot.index[iterm]+1):(res.index[iterm]-1))
  pterm <- append(pterm,rep(iterm,length(tmp)))
}

# Remove constant variance columns
# These are removed form the parameter matrix XX and from
# the pterm vector
XX.clean <- XX[,apply(XX,2,var) != 0]
pterm.clean <- pterm[apply(XX,2,var) != 0]
Constant.Variance.Names <- colnames(XX)[apply(XX,2,var)==0]
cat("These parameters have constant variance : ", Constant.Variance.Names, "\n")
df <- data.frame(XX.clean,YY)



# Set resampling method 
control <- trainControl(method="repeatedCV", number=10, repeats = 3, returnResamp = "all")

# Fitting 
fit.svmRadial <- train(YY ~., 
                       data=df, 
                       method="svmRadial", 
                       tuneLength=5,  
                       preProcess="scale", 
                       trControl=control, 
                       importance=TRUE)

# Check predition quality 
YY.pred <- predict(fit.svmRadial)
cor(YY,YY.pred)
lmfit.sum <- summary(lm(YY.pred ~ YY))
Adj.R2 <- lmfit.sum$adj.r.squared 
Fstat <- lmfit.sum$fstatistic

# Variable importance.
# VarImp function is not available for svm so what it does here is just univariate
# variable importance based on R^2 of fitted "loess" models on each outcome~predictor.
# eg. for variable 1 : lfit<-loess(YY~XX[,1]), summary(lm(YY~predict(lfit))), and look at the  R^2
# To perform a feature selection using svm see below. 
vimp <- varImp(fit.svmRadial)$importance
dfimp <- data.frame(vimp,pterm.clean)
vimp.sorted <- sort(dfimp$Overall,decreasing=TRUE,index.return=TRUE)
SortedParameters <- rownames(vimp)[vimp.sorted$ix]
SortedPotentialTerms <- unique(dfimp$pterm.clean[vimp.sorted$ix])
write.table(file="SortedParameters.dat",SortedParameters)
write.table(file="SortedPotentialTerms.dat",SortedPotentialTerms,col.names=FALSE,row.names=FALSE)

# Recursive feature elimination
# To perform a feature selection using svm you can use the recursive feature elimination function.
# It takes longer, but it seems to be equivalent to the varImp univariate method above.
# svmrfe<-rfe(XX.clean, YY, sizes=c(1,2), rfeControl = rfeControl(functions = caretFuncs,method="repeatedCV", number=10, repeats = 3), method = "svmRadial",  tuneLength=5)
# xyplot(svmrfe)