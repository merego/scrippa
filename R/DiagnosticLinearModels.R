# Diagnostic for multivariate linear models
# 

library(car)
library(mvoutlier)
library(lmtest)
library(tseries)


# Sample data
set.seed(123)
N  <- 100
X1 <- rnorm(N, 175, 7)
X2 <- rnorm(N,  30, 8)
X3 <- 0.3*X1 - 0.2*X2 + rnorm(N, 0, 5)
XX <- data.frame(X1,X2,X3) # row observation, col independent variables
Y  <- 0.5*X1 - 0.3*X2 - 0.4*X3 + 10 + rnorm(N, 0, 5) # row observation, col=1 dep. variable
dfRegr <- data.frame(XX, Y)
xyMat <- data.matrix(dfRegr)

# Nonparametric multivariate outlier detection with package 
aqRes <- aq.plot(xyMat)
outliers <- which(aqRes$outliers)

# Leverage values
cooksDst <- cooks.distance(fit)
summary(cooksDst)

# Influential values
inflRes <- influence.measures(fit)
summary(inflRes)

# Normality assumptions
Estnd <- rstandard(fit)
Estud <- rstudent(fit)

# Histogram studentized residuals
hist(Estud, main="Histogram studentized residals", breaks="FD", freq=FALSE)
curve(dnorm(x, mean=0, sd=1), col="red", lwd=2, add=TRUE)

# QQplot
qqPlot(Estud, distribution="norm", pch=20, main="QQ-Plot studentized residuals")
qqline(Estud, col="red", lwd=2)

# Shapiro normality test
shapiro.test(Estud)

# Independence and homoscedasticy
# Residuals plot
spreadLevelPlot(fit, pch=20)

# Autocorrelation
durbinWatsonTest(fit)

# Breusch-Pagan-Test homoscedasticy
bptest(fit)

# Similar to Breusch-Pagan-Test 
ncvTest(fit)

# Non linearity test
resettest(fit,power=2)
resettest(fit,power=3)
resettest(fit,power=4)


# Correlation
(Rx <- cor(XX))