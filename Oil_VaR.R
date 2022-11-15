#### Energy Economics VaR for Oil Markets #####

### Title: How Risky are the Oil Markets ? Risk Measures revisited with fat-tailed distributions ###

## Job is to find the 1-day ahead VaR for the oil markets ####
## using structural breaks and different distributions

## Load the following libraries
library(forecast)
library(rugarch)
library(fGarch)
library(PearsonDS)
library(Johnson)
library(SuppDists)   # Johnson's SU distribution
library(strucchange) # to identify structural breaks
library(esback) ## For conducting ES tests
library(GAS)  ## For backtesting
library(VaRES) # For estimating the ES 
library(PerformanceAnalytics) # For estimating the ES
library(xts)## Plot in ugarchroll works best with xts objects
library(cvar) ## Package for ES calculations
library(Dowd) ## Package for test of independence of VaR
library(tseries) ## For conducting basic tests

### We will use both Brent and WTI Oil data series ###

Brent <- read.csv("~//Desktop//Air_Desktop//Research//Energy_Related//Brent.csv")
WTI <- read.csv("~//Desktop//Air_Desktop//Research//Energy_Related//WTI.csv")  
WTI_Date <- as.Date.character(WTI$Date,format = "%d/%m/%Y")
Brent_Date   <- as.Date.character(Brent$Date,format = "%d/%m/%Y")
Brent   <- data.frame(Brent_Date,Brent$Price)     
WTI <- data.frame(WTI_Date,WTI$Price)
colnames(Brent) <- c("Date","Price")
colnames(WTI) <- c("Date","Price")

#### Structural Breaks for WTI and Brent #####

WTI_breaks <- breakpoints(WTI$Price~1)
Brent_breaks   <- breakpoints(Brent$Price~1)
WTI_samples <- WTI_breaks$breakpoints
Brent_samples <- Brent_breaks$breakpoints

Returns_Brent   <- log(Brent$Price[-1]/Brent$Price[1:(nrow(Brent)-1)])
Returns_WTI <- log(WTI$Price[-1]/WTI$Price[1:(nrow(WTI)-1)])
Des_Returns_Brent <- basicStats(Returns_Brent)
Des_Returns_WTI <- basicStats(Returns_WTI)

Holdout_Returns_Brent   <- Returns_Brent[(length(Returns_Brent)-999):length(Returns_Brent)] 
Holdout_Returns_WTI <- Returns_WTI[(length(Returns_WTI)-999):length(Returns_WTI)]

#### Convert into xts objects

WTI_xts <- xts(Returns_WTI,order.by = WTI_Date[-1])
Brent_xts   <- xts(Returns_Brent,order.by = Brent_Date[-1])
names(WTI_xts)    <- "Returns"
names(Brent_xts)   <- "Returns"

### Subsamples XTS  ###

WTI_sub1    <- WTI_xts[1:WTI_samples[1]]
WTI_sub2    <- WTI_xts[WTI_samples[1]:WTI_samples[2]]
WTI_sub3    <- WTI_xts[WTI_samples[2]:WTI_samples[3]]
WTI_sub4    <- WTI_xts[WTI_samples[3]:WTI_samples[4]]
WTI_sub5    <- WTI_xts[WTI_samples[4]:nrow(WTI_xts)]

Brent_sub1    <- Brent_xts[1:Brent_samples[1]]
Brent_sub2    <- Brent_xts[Brent_samples[1]:Brent_samples[2]]
Brent_sub3    <- Brent_xts[Brent_samples[2]:Brent_samples[3]]
Brent_sub4    <- Brent_xts[Brent_samples[3]:Brent_samples[4]]
Brent_sub5    <- Brent_xts[Brent_samples[4]:nrow(Brent_xts)]

###Dummy variables for structural breaks ####

WTI_Dummy1 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy2 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy3 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy4 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy5 <- array(data=0,dim=length(WTI_xts))

WTI_Dummy1[1:WTI_samples[1]]  <- 1
WTI_Dummy2[(WTI_samples[1]+1):WTI_samples[2]] <- 1
WTI_Dummy3[(WTI_samples[2]+1):WTI_samples[3]] <- 1
WTI_Dummy4[(WTI_samples[3]+1):WTI_samples[4]] <- 1
WTI_Dummy5[(WTI_samples[4]+1):length(WTI_xts)] <- 1

Brent_Dummy1 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy2 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy3 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy4 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy5 <- array(data=0,dim=length(Brent_xts))

Brent_Dummy1[1:Brent_samples[1]]  <- 1
Brent_Dummy2[(Brent_samples[1]+1):Brent_samples[2]] <- 1
Brent_Dummy3[(Brent_samples[2]+1):Brent_samples[3]] <- 1
Brent_Dummy4[(Brent_samples[3]+1):Brent_samples[4]] <- 1
Brent_Dummy5[(Brent_samples[4]+1):length(Brent_xts)] <- 1

WTI_Dummy <- cbind(WTI_Dummy1,WTI_Dummy2,WTI_Dummy3,WTI_Dummy4)
Brent_Dummy <- cbind(Brent_Dummy1,Brent_Dummy2,Brent_Dummy3,Brent_Dummy4)

### Descriptive Statistics and Basic Tests ###

Descr_WTI <-  cbind(basicStats(WTI_sub1),basicStats(WTI_sub2),basicStats(WTI_sub3),basicStats(WTI_sub4),basicStats(WTI_sub5))
Descr_Brent   <- cbind(basicStats(Brent_sub1),basicStats(Brent_sub2),basicStats(Brent_sub3),basicStats(Brent_sub4),basicStats(Brent_sub5))
ADF_WTI <- cbind(adf.test(WTI_sub1),adf.test(WTI_sub2),adf.test(WTI_sub3),adf.test(WTI_sub4),adf.test(WTI_sub5))
ADF_Brent   <- cbind(adf.test(Brent_sub1),adf.test(Brent_sub2),adf.test(Brent_sub3),adf.test(Brent_sub4),adf.test(Brent_sub5))
pp_WTI <- cbind(pp.test(WTI_sub1),pp.test(WTI_sub2),pp.test(WTI_sub3),pp.test(WTI_sub4),pp.test(WTI_sub5))
pp_Brent   <- cbind(pp.test(Brent_sub1),pp.test(Brent_sub2),pp.test(Brent_sub3),pp.test(Brent_sub4),pp.test(Brent_sub5))
jb_WTI <- cbind(jarque.bera.test(WTI_sub1),jarque.bera.test(WTI_sub2),jarque.bera.test(WTI_sub3),jarque.bera.test(WTI_sub4),jarque.bera.test(WTI_sub5))
jb_Brent   <- cbind(jarque.bera.test(Brent_sub1),jarque.bera.test(Brent_sub2),jarque.bera.test(Brent_sub3),jarque.bera.test(Brent_sub4),jarque.bera.test(Brent_sub5))

LQ_WTI5 <- cbind(Box.test(WTI_sub1^2,lag=5,type="Lj"),Box.test(WTI_sub2^2,lag=5,type="Lj"),Box.test(WTI_sub3^2,lag=5,type="Lj"),Box.test(WTI_sub4^2,lag=5,type="Lj"),Box.test(WTI_sub5^2,lag=5,type="Lj"))
LQ_Brent5   <- cbind(Box.test(Brent_sub1^2,lag=5,type="Lj"),Box.test(Brent_sub2^2,lag=5,type="Lj"),Box.test(Brent_sub3^2,lag=5,type="Lj"),Box.test(Brent_sub4^2,lag=5,type="Lj"),Box.test(Brent_sub5^2,lag=5,type="Lj"))

LQ_WTI10 <- cbind(Box.test(WTI_sub1^2,lag=10,type="Lj"),Box.test(WTI_sub2^2,lag=10,type="Lj"),Box.test(WTI_sub3^2,lag=10,type="Lj"),Box.test(WTI_sub4^2,lag=10,type="Lj"),Box.test(WTI_sub5^2,lag=10,type="Lj"))
LQ_Brent10   <- cbind(Box.test(Brent_sub1^2,lag=10,type="Lj"),Box.test(Brent_sub2^2,lag=10,type="Lj"),Box.test(Brent_sub3^2,lag=10,type="Lj"),Box.test(Brent_sub4^2,lag=10,type="Lj"),Box.test(Brent_sub5^2,lag=10,type="Lj"))

#### Plot the Returns and Oil Prices with Breakpoints indicated

par(mfrow=c(1,2))
plot(WTI,type="l",col="blue",main="Structural Breaks for WTI Oil",ylab="Price/Barrel (USD)",xaxt="n",xlab="")
abline(v=WTI$Date[WTI_samples[1]],col="brown")
abline(v=WTI$Date[WTI_samples[2]],col="brown")
abline(v=WTI$Date[WTI_samples[3]],col="brown")
abline(v=WTI$Date[WTI_samples[4]],col="brown")
Date_formatted <- format(WTI$Date[c(1,WTI_samples,nrow(WTI))],"%b-%y")
axis(1,at=c(min(WTI$Date),WTI$Date[WTI_samples],max(WTI$Date)),labels=Date_formatted,las = 2)

plot(Brent,type="l",col="blue",main="Structural Breaks for Brent Oil",ylab="Price/Barrel (USD)",xlab="",xaxt="n")
abline(v=Brent$Date[Brent_samples[1]],col="brown")
abline(v=Brent$Date[Brent_samples[2]],col="brown")
abline(v=Brent$Date[Brent_samples[3]],col="brown")
abline(v=Brent$Date[Brent_samples[4]],col="brown")
Date_formatted <- format(Brent$Date[c(1,Brent_samples,nrow(Brent))],"%b-%y")
axis(1,at=c(min(Brent$Date),Brent$Date[Brent_samples],max(Brent$Date)),labels=Date_formatted,las = 2)

### Plotting the Density and Histogram for Returns

par(mfrow=c(1,2))
densityPlot(as.timeSeries(WTI_xts),hist = TRUE,title = FALSE,grid=FALSE)
title(main="Density Plot for WTI Oil",xlab="Returns",ylab="Density/Frequency",col.main="blue",cex.main=0.1)
legend(x=0.06,y=25,legend = c("Normal","Empirical"),bty="n",cex=0.7,text.col = c("gray","brown"))

densityPlot(as.timeSeries(WTI_xts),hist = TRUE,title = FALSE,grid=FALSE)
title(main="Density Plot for Brent Oil",xlab="Returns",ylab="Density/Frequency",col.main="blue",cex.main=0.1)
legend(x=0.06,y=25,legend = c("Normal","Empirical"),bty="n",cex=0.7,text.col = c("gray","brown"),col=c("gray","brown"))

### Garch Models

WTI_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model1 <- ugarchfit(WTI_Spec1,WTI_xts,out.sample = 1000)
WTI_Model2 <- ugarchfit(WTI_Spec2,WTI_xts,out.sample = 1000)
WTI_Model3 <- ugarchfit(WTI_Spec3,WTI_xts,out.sample = 1000)
WTI_Model4 <- ugarchfit(WTI_Spec4,WTI_xts,out.sample = 1000)

Brent_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model1 <- ugarchfit(Brent_Spec1,Brent_xts,out.sample = 1000)
Brent_Model2 <- ugarchfit(Brent_Spec2,Brent_xts,out.sample = 1000)
Brent_Model3 <- ugarchfit(Brent_Spec3,Brent_xts,out.sample = 1000)
Brent_Model4 <- ugarchfit(Brent_Spec4,Brent_xts,out.sample = 1000)

#### egarch models

WTI_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model5 <- ugarchfit(WTI_Spec5,WTI_xts,out.sample = 1000)
WTI_Model6 <- ugarchfit(WTI_Spec6,WTI_xts,out.sample = 1000)
WTI_Model7 <- ugarchfit(WTI_Spec7,WTI_xts,out.sample = 1000)
WTI_Model8 <- ugarchfit(WTI_Spec8,WTI_xts,out.sample = 1000)

Brent_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model5 <- ugarchfit(Brent_Spec5,Brent_xts,out.sample = 1000)
Brent_Model6 <- ugarchfit(Brent_Spec6,Brent_xts,out.sample = 1000)
Brent_Model7 <- ugarchfit(Brent_Spec7,Brent_xts,out.sample = 1000)
Brent_Model8 <- ugarchfit(Brent_Spec8,Brent_xts,out.sample = 1000)

### gjrgarch models

WTI_Spec9   <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec10  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec11  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec12  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model9  <- ugarchfit(WTI_Spec9,WTI_xts,out.sample = 1000)
WTI_Model10 <- ugarchfit(WTI_Spec10,WTI_xts,out.sample = 1000)
WTI_Model11 <- ugarchfit(WTI_Spec11,WTI_xts,out.sample = 1000)
WTI_Model12 <- ugarchfit(WTI_Spec12,WTI_xts,out.sample = 1000)

Brent_Spec9  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec10 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec11 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec12 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model9  <- ugarchfit(Brent_Spec9,Brent_xts,out.sample = 1000)
Brent_Model10 <- ugarchfit(Brent_Spec10,Brent_xts,out.sample = 1000)
Brent_Model11 <- ugarchfit(Brent_Spec11,Brent_xts,out.sample = 1000)
Brent_Model12 <- ugarchfit(Brent_Spec12,Brent_xts,out.sample = 1000)

### Store the loglikelihood values, AIC and BIC ####

Llhood  <- matrix(data=0,nrow=12,ncol=2)
colnames(Llhood) <- c("Brent","WTI")
rownames(Llhood) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")
  
Llhood[1,1]  <-  Brent_Model1@fit$LLH
Llhood[2,1]  <-  Brent_Model2@fit$LLH
Llhood[3,1]  <-  Brent_Model3@fit$LLH
Llhood[4,1]  <-  Brent_Model4@fit$LLH
Llhood[5,1]  <-  Brent_Model5@fit$LLH
Llhood[6,1]  <-  Brent_Model6@fit$LLH
Llhood[7,1]  <-  Brent_Model7@fit$LLH
Llhood[8,1]  <-  Brent_Model8@fit$LLH
Llhood[9,1]  <-  Brent_Model9@fit$LLH
Llhood[10,1] <- Brent_Model10@fit$LLH
Llhood[11,1] <- Brent_Model11@fit$LLH
Llhood[12,1] <- Brent_Model12@fit$LLH

Llhood[1,2]  <-  WTI_Model1@fit$LLH
Llhood[2,2]  <-  WTI_Model2@fit$LLH
Llhood[3,2]  <-  WTI_Model3@fit$LLH
Llhood[4,2]  <-  WTI_Model4@fit$LLH
Llhood[5,2]  <-  WTI_Model5@fit$LLH
Llhood[6,2]  <-  WTI_Model6@fit$LLH
Llhood[7,2]  <-  WTI_Model7@fit$LLH
Llhood[8,2]  <-  WTI_Model8@fit$LLH
Llhood[9,2]  <-  WTI_Model9@fit$LLH
Llhood[10,2] <- WTI_Model10@fit$LLH
Llhood[11,2] <- WTI_Model11@fit$LLH
Llhood[12,2] <- WTI_Model12@fit$LLH

BIC  <- matrix(data=0,nrow=12,ncol=2)
colnames(BIC) <- c("Brent","WTI")
rownames(BIC) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")

BIC[1,1]  <-  infocriteria(Brent_Model1)[2]
BIC[2,1]  <-  infocriteria(Brent_Model2)[2]
BIC[3,1]  <-  infocriteria(Brent_Model3)[2]
BIC[4,1]  <-  infocriteria(Brent_Model4)[2]
BIC[5,1]  <-  infocriteria(Brent_Model5)[2]
BIC[6,1]  <-  infocriteria(Brent_Model6)[2]
BIC[7,1]  <-  infocriteria(Brent_Model7)[2]
BIC[8,1]  <-  infocriteria(Brent_Model8)[2]
BIC[9,1]  <-  infocriteria(Brent_Model9)[2]
BIC[10,1] <- infocriteria(Brent_Model10)[2]
BIC[11,1] <- infocriteria(Brent_Model11)[2]
BIC[12,1] <- infocriteria(Brent_Model12)[2]

BIC[1,2]  <-  infocriteria(WTI_Model1)[2]
BIC[2,2]  <-  infocriteria(WTI_Model2)[2]
BIC[3,2]  <-  infocriteria(WTI_Model3)[2]
BIC[4,2]  <-  infocriteria(WTI_Model4)[2]
BIC[5,2]  <-  infocriteria(WTI_Model5)[2]
BIC[6,2]  <-  infocriteria(WTI_Model6)[2]
BIC[7,2]  <-  infocriteria(WTI_Model7)[2]
BIC[8,2]  <-  infocriteria(WTI_Model8)[2]
BIC[9,2]  <-  infocriteria(WTI_Model9)[2]
BIC[10,2] <- infocriteria(WTI_Model10)[2]
BIC[11,2] <- infocriteria(WTI_Model11)[2]
BIC[12,2] <- infocriteria(WTI_Model12)[2]

AIC  <- matrix(data=0,nrow=12,ncol=2)
colnames(AIC) <- c("Brent","WTI")
rownames(AIC) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")

AIC[1,1]  <-  infocriteria(Brent_Model1)[1]
AIC[2,1]  <-  infocriteria(Brent_Model2)[1]
AIC[3,1]  <-  infocriteria(Brent_Model3)[1]
AIC[4,1]  <-  infocriteria(Brent_Model4)[1]
AIC[5,1]  <-  infocriteria(Brent_Model5)[1]
AIC[6,1]  <-  infocriteria(Brent_Model6)[1]
AIC[7,1]  <-  infocriteria(Brent_Model7)[1]
AIC[8,1]  <-  infocriteria(Brent_Model8)[1]
AIC[9,1]  <-  infocriteria(Brent_Model9)[1]
AIC[10,1] <- infocriteria(Brent_Model10)[1]
AIC[11,1] <- infocriteria(Brent_Model11)[1]
AIC[12,1] <- infocriteria(Brent_Model12)[1]

AIC[1,2]  <-  infocriteria(WTI_Model1)[1]
AIC[2,2]  <-  infocriteria(WTI_Model2)[1]
AIC[3,2]  <-  infocriteria(WTI_Model3)[1]
AIC[4,2]  <-  infocriteria(WTI_Model4)[1]
AIC[5,2]  <-  infocriteria(WTI_Model5)[1]
AIC[6,2]  <-  infocriteria(WTI_Model6)[1]
AIC[7,2]  <-  infocriteria(WTI_Model7)[1]
AIC[8,2]  <-  infocriteria(WTI_Model8)[1]
AIC[9,2]  <-  infocriteria(WTI_Model9)[1]
AIC[10,2] <- infocriteria(WTI_Model10)[1]
AIC[11,2] <- infocriteria(WTI_Model11)[1]
AIC[12,2] <- infocriteria(WTI_Model12)[1]

## Conclusion: eGarch with SST and JSU are the best models #####
## Model5 for SST and Model7 for JSU, Model8 for Normal required to calculate Pearson 

### Store the coefficients of the models

WTI_Coef5   <- round(coef(WTI_Model5),5)
WTI_Coef6   <- round(coef(WTI_Model6),5)
WTI_Coef7   <- round(coef(WTI_Model7),5)
WTI_Coef8   <- round(coef(WTI_Model8),5)
WTI_Coef_Matrix   <- list(WTI_Coef5,WTI_Coef6,WTI_Coef7,WTI_Coef8)
names(WTI_Coef_Matrix) <- c("SST","SGED","JSU","Normal")

WTI_Robust_t5   <- round(WTI_Model5@fit$robust.tval,2)
WTI_Robust_t6   <- round(WTI_Model6@fit$robust.tval,2)
WTI_Robust_t7   <- round(WTI_Model7@fit$robust.tval,2)
WTI_Robust_t8   <- round(WTI_Model8@fit$robust.tval,2)

WTI_Robust_t   <- list(WTI_Robust_t5,WTI_Robust_t6,WTI_Robust_t7,WTI_Robust_t8)
names(WTI_Robust_t) <- c("SST","SGED","JSU","Normal")

Brent_Coef5   <- round(coef(Brent_Model5),5)
Brent_Coef6   <- round(coef(Brent_Model6),5)
Brent_Coef7   <- round(coef(Brent_Model7),5)
Brent_Coef8   <- round(coef(Brent_Model8),5)
Brent_Coef_Matrix   <- list(Brent_Coef5,Brent_Coef6,Brent_Coef7,Brent_Coef8)
names(Brent_Coef_Matrix) <- c("SST","SGED","JSU","Normal")

Brent_Robust_t5   <- round(Brent_Model5@fit$robust.tval,2)
Brent_Robust_t6   <- round(Brent_Model6@fit$robust.tval,2)
Brent_Robust_t7   <- round(Brent_Model7@fit$robust.tval,2)
Brent_Robust_t8   <- round(Brent_Model8@fit$robust.tval,2)

Brent_Robust_t   <- list(Brent_Robust_t5,Brent_Robust_t6,Brent_Robust_t7,Brent_Robust_t8)
names(Brent_Robust_t) <- c("SST","SGED","JSU","Normal")

### Ljung-box test statistic

WTI5_LQ1  <- Box.test(WTI_Model5@fit$residuals/WTI_Model5@fit$sigma,lag=1,type="Lj")
WTI5_LQ5  <- Box.test(WTI_Model5@fit$residuals/WTI_Model5@fit$sigma,lag=5,type="Lj")
WTI5_LQ10 <- Box.test(WTI_Model5@fit$residuals/WTI_Model5@fit$sigma,lag=10,type="Lj")

WTI6_LQ1  <- Box.test(WTI_Model6@fit$residuals/WTI_Model6@fit$sigma,lag=1,type="Lj")
WTI6_LQ5  <- Box.test(WTI_Model6@fit$residuals/WTI_Model6@fit$sigma,lag=5,type="Lj")
WTI6_LQ10 <- Box.test(WTI_Model6@fit$residuals/WTI_Model6@fit$sigma,lag=10,type="Lj")

WTI7_LQ1  <- Box.test(WTI_Model7@fit$residuals/WTI_Model7@fit$sigma,lag=1,type="Lj")
WTI7_LQ5  <- Box.test(WTI_Model7@fit$residuals/WTI_Model7@fit$sigma,lag=5,type="Lj")
WTI7_LQ10 <- Box.test(WTI_Model7@fit$residuals/WTI_Model7@fit$sigma,lag=10,type="Lj")

WTI8_LQ1  <- Box.test(WTI_Model8@fit$residuals/WTI_Model8@fit$sigma,lag=1,type="Lj")
WTI8_LQ5  <- Box.test(WTI_Model8@fit$residuals/WTI_Model8@fit$sigma,lag=5,type="Lj")
WTI8_LQ10 <- Box.test(WTI_Model8@fit$residuals/WTI_Model8@fit$sigma,lag=10,type="Lj")

WTI5_LQ <- c(WTI5_LQ1$statistic,WTI5_LQ5$statistic,WTI5_LQ10$statistic)
WTI6_LQ <- c(WTI6_LQ1$statistic,WTI6_LQ5$statistic,WTI6_LQ10$statistic)
WTI7_LQ <- c(WTI7_LQ1$statistic,WTI7_LQ5$statistic,WTI7_LQ10$statistic)
WTI8_LQ <- c(WTI8_LQ1$statistic,WTI8_LQ5$statistic,WTI8_LQ10$statistic)

Brent5_LQ1  <- Box.test(Brent_Model5@fit$residuals/Brent_Model5@fit$sigma,lag=1,type="Lj")
Brent5_LQ5  <- Box.test(Brent_Model5@fit$residuals/Brent_Model5@fit$sigma,lag=5,type="Lj")
Brent5_LQ10 <- Box.test(Brent_Model5@fit$residuals/Brent_Model5@fit$sigma,lag=10,type="Lj")

Brent6_LQ1  <- Box.test(Brent_Model6@fit$residuals/Brent_Model6@fit$sigma,lag=1,type="Lj")
Brent6_LQ5  <- Box.test(Brent_Model6@fit$residuals/Brent_Model6@fit$sigma,lag=5,type="Lj")
Brent6_LQ10 <- Box.test(Brent_Model6@fit$residuals/Brent_Model6@fit$sigma,lag=10,type="Lj")

Brent7_LQ1  <- Box.test(Brent_Model7@fit$residuals/Brent_Model7@fit$sigma,lag=1,type="Lj")
Brent7_LQ5  <- Box.test(Brent_Model7@fit$residuals/Brent_Model7@fit$sigma,lag=5,type="Lj")
Brent7_LQ10 <- Box.test(Brent_Model7@fit$residuals/Brent_Model7@fit$sigma,lag=10,type="Lj")

Brent8_LQ1  <- Box.test(Brent_Model8@fit$residuals/Brent_Model8@fit$sigma,lag=1,type="Lj")
Brent8_LQ5  <- Box.test(Brent_Model8@fit$residuals/Brent_Model8@fit$sigma,lag=5,type="Lj")
Brent8_LQ10 <- Box.test(Brent_Model8@fit$residuals/Brent_Model8@fit$sigma,lag=10,type="Lj")

Brent5_LQ <- c(Brent5_LQ1$statistic,Brent5_LQ5$statistic,Brent5_LQ10$statistic)
Brent6_LQ <- c(Brent6_LQ1$statistic,Brent6_LQ5$statistic,Brent6_LQ10$statistic)
Brent7_LQ <- c(Brent7_LQ1$statistic,Brent7_LQ5$statistic,Brent7_LQ10$statistic)
Brent8_LQ <- c(Brent8_LQ1$statistic,Brent8_LQ5$statistic,Brent8_LQ10$statistic)

WTI_res_stat <-   cbind(WTI5_LQ,WTI6_LQ,WTI7_LQ,WTI8_LQ)
rownames(WTI_res_stat) <- c("LQ1","LQ5","LQ10")
colnames(WTI_res_stat) <- c("SST","SGED","JSU","Normal")

Brent_res_stat <-   cbind(Brent5_LQ,Brent6_LQ,Brent7_LQ,Brent8_LQ)
rownames(Brent_res_stat) <- c("LQ1","LQ5","LQ10")
colnames(Brent_res_stat) <- c("SST","SGED","JSU","Normal")

#### 1-day ahead Rolling Forecasts ####
## Full sample for roll i.e. last 1000 points for forecast and rest for model

WTI_Roll5 <-  ugarchroll(WTI_Spec5,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll6 <-  ugarchroll(WTI_Spec6,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll7 <-  ugarchroll(WTI_Spec7,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll8 <-  ugarchroll(WTI_Spec8,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll5 <-  ugarchroll(Brent_Spec5,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll6 <-  ugarchroll(Brent_Spec6,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll7 <-  ugarchroll(Brent_Spec7,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll8 <-  ugarchroll(Brent_Spec8,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


#### Store the VaR for the out sample using rolling window forecasts

WTI_VaR_SST   <- WTI_Roll5@forecast$VaR
WTI_VaR_SGED  <- WTI_Roll6@forecast$VaR
WTI_VaR_JSU   <- WTI_Roll7@forecast$VaR
WTI_VaR_Norm  <- WTI_Roll8@forecast$VaR

Brent_VaR_SST   <- Brent_Roll5@forecast$VaR
Brent_VaR_SGED  <- Brent_Roll6@forecast$VaR
Brent_VaR_JSU   <- Brent_Roll7@forecast$VaR
Brent_VaR_Norm  <- Brent_Roll8@forecast$VaR

## Significance Levels

signific_levels <- c(0.01,0.05,0.10)

### Backtest the results

WTI_Test_SST_1  <- BacktestVaR(WTI_VaR_SST[,4],WTI_VaR_SST[,1],signific_levels[1])
WTI_Test_SST_5  <- BacktestVaR(WTI_VaR_SST[,4],WTI_VaR_SST[,2],signific_levels[2])
WTI_Test_SST_10 <- BacktestVaR(WTI_VaR_SST[,4],WTI_VaR_SST[,3],signific_levels[3])

WTI_Test_SGED_1  <- BacktestVaR(WTI_VaR_SGED[,4],WTI_VaR_SGED[,1],signific_levels[1])
WTI_Test_SGED_5  <- BacktestVaR(WTI_VaR_SGED[,4],WTI_VaR_SGED[,2],signific_levels[2])
WTI_Test_SGED_10 <- BacktestVaR(WTI_VaR_SGED[,4],WTI_VaR_SGED[,3],signific_levels[3])

WTI_Test_JSU_1  <- BacktestVaR(WTI_VaR_JSU[,4],WTI_VaR_JSU[,1],signific_levels[1])
WTI_Test_JSU_5  <- BacktestVaR(WTI_VaR_JSU[,4],WTI_VaR_JSU[,2],signific_levels[2])
WTI_Test_JSU_10 <- BacktestVaR(WTI_VaR_JSU[,4],WTI_VaR_JSU[,3],signific_levels[3])

WTI_Test_Norm_1  <- BacktestVaR(WTI_VaR_Norm[,4],WTI_VaR_Norm[,1],signific_levels[1])
WTI_Test_Norm_5  <- BacktestVaR(WTI_VaR_Norm[,4],WTI_VaR_Norm[,2],signific_levels[2])
WTI_Test_Norm_10 <- BacktestVaR(WTI_VaR_Norm[,4],WTI_VaR_Norm[,3],signific_levels[3])

Brent_Test_SST_1  <- BacktestVaR(Brent_VaR_SST[,4],Brent_VaR_SST[,1],signific_levels[1])
Brent_Test_SST_5  <- BacktestVaR(Brent_VaR_SST[,4],Brent_VaR_SST[,2],signific_levels[2])
Brent_Test_SST_10 <- BacktestVaR(Brent_VaR_SST[,4],Brent_VaR_SST[,3],signific_levels[3])

Brent_Test_SGED_1  <- BacktestVaR(Brent_VaR_SGED[,4],Brent_VaR_SGED[,1],signific_levels[1])
Brent_Test_SGED_5  <- BacktestVaR(Brent_VaR_SGED[,4],Brent_VaR_SGED[,2],signific_levels[2])
Brent_Test_SGED_10 <- BacktestVaR(Brent_VaR_SGED[,4],Brent_VaR_SGED[,3],signific_levels[3])

Brent_Test_JSU_1  <- BacktestVaR(Brent_VaR_JSU[,4],Brent_VaR_JSU[,1],signific_levels[1])
Brent_Test_JSU_5  <- BacktestVaR(Brent_VaR_JSU[,4],Brent_VaR_JSU[,2],signific_levels[2])
Brent_Test_JSU_10 <- BacktestVaR(Brent_VaR_JSU[,4],Brent_VaR_JSU[,3],signific_levels[3])

Brent_Test_Norm_1  <- BacktestVaR(Brent_VaR_Norm[,4],Brent_VaR_Norm[,1],signific_levels[1])
Brent_Test_Norm_5  <- BacktestVaR(Brent_VaR_Norm[,4],Brent_VaR_Norm[,2],signific_levels[2])
Brent_Test_Norm_10 <- BacktestVaR(Brent_VaR_Norm[,4],Brent_VaR_Norm[,3],signific_levels[3])

#### Saving the test stats and pvalues for 1%, 5% and 10% VaR

## 1%
WTI_Backtest_teststat1 <- matrix(c(WTI_Test_SST_1$LRuc[1],WTI_Test_SST_1$LRcc[1],WTI_Test_SST_1$DQ$stat,
                            WTI_Test_SGED_1$LRuc[1],WTI_Test_SGED_1$LRcc[1],WTI_Test_SGED_1$DQ$stat,
                            WTI_Test_JSU_1$LRuc[1],WTI_Test_JSU_1$LRcc[1],WTI_Test_JSU_1$DQ$stat,
                            WTI_Test_Norm_1$LRuc[1],WTI_Test_Norm_1$LRcc[1],WTI_Test_Norm_1$DQ$stat),nrow=4,byrow = TRUE)
rownames(WTI_Backtest_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_Backtest_teststat1)  <- c("Kupiec","Christoffersen","DQ")

Brent_Backtest_teststat1  <- matrix(c(Brent_Test_SST_1$LRuc[1],Brent_Test_SST_1$LRcc[1],Brent_Test_SST_1$DQ$stat,
                            Brent_Test_SGED_1$LRuc[1],Brent_Test_SGED_1$LRcc[1],Brent_Test_SGED_1$DQ$stat,
                            Brent_Test_JSU_1$LRuc[1],Brent_Test_JSU_1$LRcc[1],Brent_Test_JSU_1$DQ$stat,
                            Brent_Test_Norm_1$LRuc[1],Brent_Test_Norm_1$LRcc[1],Brent_Test_Norm_1$DQ$stat),nrow=4,byrow = TRUE)
rownames(Brent_Backtest_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_Backtest_teststat1)  <- c("Kupiec","Christoffersen","DQ")

WTI_Backtest_pvalues1 <- matrix(c(WTI_Test_SST_1$LRuc[2],WTI_Test_SST_1$LRcc[2],WTI_Test_SST_1$DQ$pvalue,
                             WTI_Test_SGED_1$LRuc[2],WTI_Test_SGED_1$LRcc[2],WTI_Test_SGED_1$DQ$pvalue,
                             WTI_Test_JSU_1$LRuc[2],WTI_Test_JSU_1$LRcc[2],WTI_Test_JSU_1$DQ$pvalue,
                             WTI_Test_Norm_1$LRuc[2],WTI_Test_Norm_1$LRcc[2],WTI_Test_Norm_1$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(WTI_Backtest_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_Backtest_pvalues1)  <- c("Kupiec","Christoffersen","DQ")

Brent_Backtest_pvalues1  <- matrix(c(Brent_Test_SST_1$LRuc[2],Brent_Test_SST_1$LRcc[2],Brent_Test_SST_1$DQ$pvalue,
                            Brent_Test_SGED_1$LRuc[2],Brent_Test_SGED_1$LRcc[2],Brent_Test_SGED_1$DQ$pvalue,
                            Brent_Test_JSU_1$LRuc[2],Brent_Test_JSU_1$LRcc[2],Brent_Test_JSU_1$DQ$pvalue,
                            Brent_Test_Norm_1$LRuc[2],Brent_Test_Norm_1$LRcc[2],Brent_Test_Norm_1$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Brent_Backtest_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_Backtest_pvalues1)  <- c("Kupiec","Christoffersen","DQ")

### 5%
WTI_Backtest_teststat5 <- matrix(c(WTI_Test_SST_5$LRuc[1],WTI_Test_SST_5$LRcc[1],WTI_Test_SST_5$DQ$stat,
                                     WTI_Test_SGED_5$LRuc[1],WTI_Test_SGED_5$LRcc[1],WTI_Test_SGED_5$DQ$stat,
                                     WTI_Test_JSU_5$LRuc[1],WTI_Test_JSU_5$LRcc[1],WTI_Test_JSU_5$DQ$stat,
                                     WTI_Test_Norm_5$LRuc[1],WTI_Test_Norm_5$LRcc[1],WTI_Test_Norm_5$DQ$stat),nrow=4,byrow = TRUE)
rownames(WTI_Backtest_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_Backtest_teststat5)  <- c("Kupiec","Christoffersen","DQ")

Brent_Backtest_teststat5  <- matrix(c(Brent_Test_SST_5$LRuc[1],Brent_Test_SST_5$LRcc[1],Brent_Test_SST_5$DQ$stat,
                                    Brent_Test_SGED_5$LRuc[1],Brent_Test_SGED_5$LRcc[1],Brent_Test_SGED_5$DQ$stat,
                                    Brent_Test_JSU_5$LRuc[1],Brent_Test_JSU_5$LRcc[1],Brent_Test_JSU_5$DQ$stat,
                                    Brent_Test_Norm_5$LRuc[1],Brent_Test_Norm_5$LRcc[1],Brent_Test_Norm_5$DQ$stat),nrow=4,byrow = TRUE)
rownames(Brent_Backtest_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_Backtest_teststat5)  <- c("Kupiec","Christoffersen","DQ")

WTI_Backtest_pvalues5 <- matrix(c(WTI_Test_SST_5$LRuc[2],WTI_Test_SST_5$LRcc[2],WTI_Test_SST_5$DQ$pvalue,
                                    WTI_Test_SGED_5$LRuc[2],WTI_Test_SGED_5$LRcc[2],WTI_Test_SGED_5$DQ$pvalue,
                                    WTI_Test_JSU_5$LRuc[2],WTI_Test_JSU_5$LRcc[2],WTI_Test_JSU_5$DQ$pvalue,
                                    WTI_Test_Norm_5$LRuc[2],WTI_Test_Norm_5$LRcc[2],WTI_Test_Norm_5$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(WTI_Backtest_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_Backtest_pvalues5)  <- c("Kupiec","Christoffersen","DQ")

Brent_Backtest_pvalues5  <- matrix(c(Brent_Test_SST_5$LRuc[2],Brent_Test_SST_5$LRcc[2],Brent_Test_SST_5$DQ$pvalue,
                                   Brent_Test_SGED_5$LRuc[2],Brent_Test_SGED_5$LRcc[2],Brent_Test_SGED_5$DQ$pvalue,
                                   Brent_Test_JSU_5$LRuc[2],Brent_Test_JSU_5$LRcc[2],Brent_Test_JSU_5$DQ$pvalue,
                                   Brent_Test_Norm_5$LRuc[2],Brent_Test_Norm_5$LRcc[2],Brent_Test_Norm_5$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Brent_Backtest_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_Backtest_pvalues5)  <- c("Kupiec","Christoffersen","DQ")

### 10%

WTI_Backtest_teststat10 <- matrix(c(WTI_Test_SST_10$LRuc[1],WTI_Test_SST_10$LRcc[1],WTI_Test_SST_10$DQ$stat,
                                     WTI_Test_SGED_10$LRuc[1],WTI_Test_SGED_10$LRcc[1],WTI_Test_SGED_10$DQ$stat,
                                     WTI_Test_JSU_10$LRuc[1],WTI_Test_JSU_10$LRcc[1],WTI_Test_JSU_10$DQ$stat,
                                     WTI_Test_Norm_10$LRuc[1],WTI_Test_Norm_10$LRcc[1],WTI_Test_Norm_10$DQ$stat),nrow=4,byrow = TRUE)
rownames(WTI_Backtest_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_Backtest_teststat10)  <- c("Kupiec","Christoffersen","DQ")

Brent_Backtest_teststat10  <- matrix(c(Brent_Test_SST_10$LRuc[1],Brent_Test_SST_10$LRcc[1],Brent_Test_SST_10$DQ$stat,
                                    Brent_Test_SGED_10$LRuc[1],Brent_Test_SGED_10$LRcc[1],Brent_Test_SGED_10$DQ$stat,
                                    Brent_Test_JSU_10$LRuc[1],Brent_Test_JSU_10$LRcc[1],Brent_Test_JSU_10$DQ$stat,
                                    Brent_Test_Norm_10$LRuc[1],Brent_Test_Norm_10$LRcc[1],Brent_Test_Norm_10$DQ$stat),nrow=4,byrow = TRUE)
rownames(Brent_Backtest_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_Backtest_teststat10)  <- c("Kupiec","Christoffersen","DQ")

WTI_Backtest_pvalues10 <- matrix(c(WTI_Test_SST_10$LRuc[2],WTI_Test_SST_10$LRcc[2],WTI_Test_SST_10$DQ$pvalue,
                                    WTI_Test_SGED_10$LRuc[2],WTI_Test_SGED_10$LRcc[2],WTI_Test_SGED_10$DQ$pvalue,
                                    WTI_Test_JSU_10$LRuc[2],WTI_Test_JSU_10$LRcc[2],WTI_Test_JSU_10$DQ$pvalue,
                                    WTI_Test_Norm_10$LRuc[2],WTI_Test_Norm_10$LRcc[2],WTI_Test_Norm_10$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(WTI_Backtest_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_Backtest_pvalues10)  <- c("Kupiec","Christoffersen","DQ")

Brent_Backtest_pvalues10  <- matrix(c(Brent_Test_SST_10$LRuc[2],Brent_Test_SST_10$LRcc[2],Brent_Test_SST_10$DQ$pvalue,
                                   Brent_Test_SGED_10$LRuc[2],Brent_Test_SGED_10$LRcc[2],Brent_Test_SGED_10$DQ$pvalue,
                                   Brent_Test_JSU_10$LRuc[2],Brent_Test_JSU_10$LRcc[2],Brent_Test_JSU_10$DQ$pvalue,
                                   Brent_Test_Norm_10$LRuc[2],Brent_Test_Norm_10$LRcc[2],Brent_Test_Norm_10$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Brent_Backtest_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_Backtest_pvalues10)  <- c("Kupiec","Christoffersen","DQ")

WTI_Backtest_pvalues1  <- round(WTI_Backtest_pvalues1,4)
WTI_Backtest_pvalues5  <- round(WTI_Backtest_pvalues5,4)
WTI_Backtest_pvalues10 <- round(WTI_Backtest_pvalues10,4)

Brent_Backtest_pvalues1  <- round(Brent_Backtest_pvalues1,4)
Brent_Backtest_pvalues5  <- round(Brent_Backtest_pvalues5,4)
Brent_Backtest_pvalues10 <- round(Brent_Backtest_pvalues10,4)

WTI_Backtest_teststat1  <- round(WTI_Backtest_teststat1,4)
WTI_Backtest_teststat5  <- round(WTI_Backtest_teststat5,4)
WTI_Backtest_teststat10 <- round(WTI_Backtest_teststat10,4)

Brent_Backtest_teststat1  <- round(Brent_Backtest_teststat1,4)
Brent_Backtest_teststat5  <- round(Brent_Backtest_teststat5,4)
Brent_Backtest_teststat10 <- round(Brent_Backtest_teststat10,4)

### Using VaRTest and VaRDuration Test individually ###

## VaR Test in R

WTI_Norm_VaR_Test1  <- VaRTest(0.01,WTI_VaR_Norm[,4],WTI_VaR_Norm[,1])
WTI_Norm_VaR_Test5  <- VaRTest(0.05,WTI_VaR_Norm[,4],WTI_VaR_Norm[,2])
WTI_Norm_VaR_Test10 <- VaRTest(0.1,WTI_VaR_Norm[,4],WTI_VaR_Norm[,3])

WTI_JSU_VaR_Test1  <- VaRTest(0.01,WTI_VaR_JSU[,4],WTI_VaR_JSU[,1])
WTI_JSU_VaR_Test5  <- VaRTest(0.05,WTI_VaR_JSU[,4],WTI_VaR_JSU[,2])
WTI_JSU_VaR_Test10 <- VaRTest(0.1,WTI_VaR_JSU[,4],WTI_VaR_JSU[,3])

WTI_SGED_VaR_Test1  <- VaRTest(0.01,WTI_VaR_SGED[,4],WTI_VaR_SGED[,1])
WTI_SGED_VaR_Test5  <- VaRTest(0.05,WTI_VaR_SGED[,4],WTI_VaR_SGED[,2])
WTI_SGED_VaR_Test10 <- VaRTest(0.1,WTI_VaR_SGED[,4],WTI_VaR_SGED[,3])

WTI_SST_VaR_Test1  <- VaRTest(0.01,WTI_VaR_SST[,4],WTI_VaR_SST[,1])
WTI_SST_VaR_Test5  <- VaRTest(0.05,WTI_VaR_SST[,4],WTI_VaR_SST[,2])
WTI_SST_VaR_Test10 <- VaRTest(0.1,WTI_VaR_SST[,4],WTI_VaR_SST[,3])

Brent_Norm_VaR_Test1  <- VaRTest(0.01,Brent_VaR_Norm[,4],Brent_VaR_Norm[,1])
Brent_Norm_VaR_Test5  <- VaRTest(0.05,Brent_VaR_Norm[,4],Brent_VaR_Norm[,2])
Brent_Norm_VaR_Test10 <- VaRTest(0.1,Brent_VaR_Norm[,4],Brent_VaR_Norm[,3])

Brent_JSU_VaR_Test1  <- VaRTest(0.01,Brent_VaR_JSU[,4],Brent_VaR_JSU[,1])
Brent_JSU_VaR_Test5  <- VaRTest(0.05,Brent_VaR_JSU[,4],Brent_VaR_JSU[,2])
Brent_JSU_VaR_Test10 <- VaRTest(0.1,Brent_VaR_JSU[,4],Brent_VaR_JSU[,3])

Brent_SGED_VaR_Test1  <- VaRTest(0.01,Brent_VaR_SGED[,4],Brent_VaR_SGED[,1])
Brent_SGED_VaR_Test5  <- VaRTest(0.05,Brent_VaR_SGED[,4],Brent_VaR_SGED[,2])
Brent_SGED_VaR_Test10 <- VaRTest(0.1,Brent_VaR_SGED[,4],Brent_VaR_SGED[,3])

Brent_SST_VaR_Test1  <- VaRTest(0.01,Brent_VaR_SST[,4],Brent_VaR_SST[,1])
Brent_SST_VaR_Test5  <- VaRTest(0.05,Brent_VaR_SST[,4],Brent_VaR_SST[,2])
Brent_SST_VaR_Test10 <- VaRTest(0.1,Brent_VaR_SST[,4],Brent_VaR_SST[,3])

## Store the pvalues at a place

WTI_VaR_Test_pvalues1 <- matrix(c(WTI_SST_VaR_Test1$uc.LRp,WTI_SST_VaR_Test1$cc.LRp,
                                    WTI_SGED_VaR_Test1$uc.LRp,WTI_SGED_VaR_Test1$cc.LRp,
                                    WTI_JSU_VaR_Test1$uc.LRp,WTI_JSU_VaR_Test1$cc.LRp,
                                    WTI_Norm_VaR_Test1$uc.LRp,WTI_Norm_VaR_Test1$cc.LRp),ncol=2,byrow = TRUE)

WTI_VaR_Test_pvalues5 <- matrix(c(WTI_SST_VaR_Test5$uc.LRp,WTI_SST_VaR_Test5$cc.LRp,
                                    WTI_SGED_VaR_Test5$uc.LRp,WTI_SGED_VaR_Test5$cc.LRp,
                                    WTI_JSU_VaR_Test5$uc.LRp,WTI_JSU_VaR_Test5$cc.LRp,
                                    WTI_Norm_VaR_Test5$uc.LRp,WTI_Norm_VaR_Test5$cc.LRp),ncol=2,byrow = TRUE)

WTI_VaR_Test_pvalues10 <- matrix(c(WTI_SST_VaR_Test10$uc.LRp,WTI_SST_VaR_Test10$cc.LRp,
                                    WTI_SGED_VaR_Test10$uc.LRp,WTI_SGED_VaR_Test10$cc.LRp,
                                    WTI_JSU_VaR_Test10$uc.LRp,WTI_JSU_VaR_Test10$cc.LRp,
                                    WTI_Norm_VaR_Test10$uc.LRp,WTI_Norm_VaR_Test10$cc.LRp),ncol=2,byrow = TRUE)

Brent_VaR_Test_pvalues1 <- matrix(c(Brent_SST_VaR_Test1$uc.LRp,Brent_SST_VaR_Test1$cc.LRp,
                                    Brent_SGED_VaR_Test1$uc.LRp,Brent_SGED_VaR_Test1$cc.LRp,
                                    Brent_JSU_VaR_Test1$uc.LRp,Brent_JSU_VaR_Test1$cc.LRp,
                                    Brent_Norm_VaR_Test1$uc.LRp,Brent_Norm_VaR_Test1$cc.LRp),ncol=2,byrow = TRUE)

Brent_VaR_Test_pvalues5 <- matrix(c(Brent_SST_VaR_Test5$uc.LRp,Brent_SST_VaR_Test5$cc.LRp,
                                    Brent_SGED_VaR_Test5$uc.LRp,Brent_SGED_VaR_Test5$cc.LRp,
                                    Brent_JSU_VaR_Test5$uc.LRp,Brent_JSU_VaR_Test5$cc.LRp,
                                    Brent_Norm_VaR_Test5$uc.LRp,Brent_Norm_VaR_Test5$cc.LRp),ncol=2,byrow = TRUE)

Brent_VaR_Test_pvalues10 <- matrix(c(Brent_SST_VaR_Test10$uc.LRp,Brent_SST_VaR_Test10$cc.LRp,
                                     Brent_SGED_VaR_Test10$uc.LRp,Brent_SGED_VaR_Test10$cc.LRp,
                                     Brent_JSU_VaR_Test10$uc.LRp,Brent_JSU_VaR_Test10$cc.LRp,
                                     Brent_Norm_VaR_Test10$uc.LRp,Brent_Norm_VaR_Test10$cc.LRp),ncol=2,byrow = TRUE)

## Storing the test stats
WTI_VaR_Test_teststat1 <- matrix(c(WTI_SST_VaR_Test1$uc.LRstat,WTI_SST_VaR_Test1$cc.LRstat,
                                    WTI_SGED_VaR_Test1$uc.LRstat,WTI_SGED_VaR_Test1$cc.LRstat,
                                    WTI_JSU_VaR_Test1$uc.LRstat,WTI_JSU_VaR_Test1$cc.LRstat,
                                    WTI_Norm_VaR_Test1$uc.LRstat,WTI_Norm_VaR_Test1$cc.LRstat),ncol=2,byrow = TRUE)

WTI_VaR_Test_teststat5 <- matrix(c(WTI_SST_VaR_Test5$uc.LRstat,WTI_SST_VaR_Test5$cc.LRstat,
                                    WTI_SGED_VaR_Test5$uc.LRstat,WTI_SGED_VaR_Test5$cc.LRstat,
                                    WTI_JSU_VaR_Test5$uc.LRstat,WTI_JSU_VaR_Test5$cc.LRstat,
                                    WTI_Norm_VaR_Test5$uc.LRstat,WTI_Norm_VaR_Test5$cc.LRstat),ncol=2,byrow = TRUE)

WTI_VaR_Test_teststat10 <- matrix(c(WTI_SST_VaR_Test10$uc.LRstat,WTI_SST_VaR_Test10$cc.LRstat,
                                     WTI_SGED_VaR_Test10$uc.LRstat,WTI_SGED_VaR_Test10$cc.LRstat,
                                     WTI_JSU_VaR_Test10$uc.LRstat,WTI_JSU_VaR_Test10$cc.LRstat,
                                     WTI_Norm_VaR_Test10$uc.LRstat,WTI_Norm_VaR_Test10$cc.LRstat),ncol=2,byrow = TRUE)

Brent_VaR_Test_teststat1 <- matrix(c(Brent_SST_VaR_Test1$uc.LRstat,Brent_SST_VaR_Test1$cc.LRstat,
                                  Brent_SGED_VaR_Test1$uc.LRstat,Brent_SGED_VaR_Test1$cc.LRstat,
                                  Brent_JSU_VaR_Test1$uc.LRstat,Brent_JSU_VaR_Test1$cc.LRstat,
                                  Brent_Norm_VaR_Test1$uc.LRstat,Brent_Norm_VaR_Test1$cc.LRstat),ncol=2,byrow = TRUE)

Brent_VaR_Test_teststat5 <- matrix(c(Brent_SST_VaR_Test5$uc.LRstat,Brent_SST_VaR_Test5$cc.LRstat,
                                  Brent_SGED_VaR_Test5$uc.LRstat,Brent_SGED_VaR_Test5$cc.LRstat,
                                  Brent_JSU_VaR_Test5$uc.LRstat,Brent_JSU_VaR_Test5$cc.LRstat,
                                  Brent_Norm_VaR_Test5$uc.LRstat,Brent_Norm_VaR_Test5$cc.LRstat),ncol=2,byrow = TRUE)

Brent_VaR_Test_teststat10 <- matrix(c(Brent_SST_VaR_Test10$uc.LRstat,Brent_SST_VaR_Test10$cc.LRstat,
                                   Brent_SGED_VaR_Test10$uc.LRstat,Brent_SGED_VaR_Test10$cc.LRstat,
                                   Brent_JSU_VaR_Test10$uc.LRstat,Brent_JSU_VaR_Test10$cc.LRstat,
                                   Brent_Norm_VaR_Test10$uc.LRstat,Brent_Norm_VaR_Test10$cc.LRstat),ncol=2,byrow = TRUE)
## Round them up to 4 decimals

WTI_VaR_Test_pvalues1  <- round(WTI_VaR_Test_pvalues1,4)
WTI_VaR_Test_pvalues5  <- round(WTI_VaR_Test_pvalues5,4)
WTI_VaR_Test_pvalues10 <- round(WTI_VaR_Test_pvalues10,4)

Brent_VaR_Test_pvalues1  <- round(Brent_VaR_Test_pvalues1,4)
Brent_VaR_Test_pvalues5  <- round(Brent_VaR_Test_pvalues5,4)
Brent_VaR_Test_pvalues10 <- round(Brent_VaR_Test_pvalues10,4)

WTI_VaR_Test_teststat1  <- round(WTI_VaR_Test_teststat1,4)
WTI_VaR_Test_teststat5  <- round(WTI_VaR_Test_teststat5,4)
WTI_VaR_Test_teststat10 <- round(WTI_VaR_Test_teststat10,4)

Brent_VaR_Test_teststat1  <- round(Brent_VaR_Test_teststat1,4)
Brent_VaR_Test_teststat5  <- round(Brent_VaR_Test_teststat5,4)
Brent_VaR_Test_teststat10 <- round(Brent_VaR_Test_teststat10,4)

rownames(WTI_VaR_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Test_pvalues1)  <- c("Kupiec","Christoffersen")
rownames(WTI_VaR_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Test_pvalues5)  <- c("Kupiec","Christoffersen")
rownames(WTI_VaR_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Test_pvalues10)  <- c("Kupiec","Christoffersen")

rownames(Brent_VaR_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Test_pvalues1)  <- c("Kupiec","Christoffersen")
rownames(Brent_VaR_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Test_pvalues5)  <- c("Kupiec","Christoffersen")
rownames(Brent_VaR_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Test_pvalues10)  <- c("Kupiec","Christoffersen")

rownames(WTI_VaR_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Test_teststat1)  <- c("Kupiec","Christoffersen")
rownames(WTI_VaR_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Test_teststat5)  <- c("Kupiec","Christoffersen")
rownames(WTI_VaR_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Test_teststat10)  <- c("Kupiec","Christoffersen")

rownames(Brent_VaR_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Test_teststat1)  <- c("Kupiec","Christoffersen")
rownames(Brent_VaR_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Test_teststat5)  <- c("Kupiec","Christoffersen")
rownames(Brent_VaR_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Test_teststat10)  <- c("Kupiec","Christoffersen")


## VaR Duration test

WTI_Norm_VaR_Dur_Test1  <- VaRDurTest(0.01,WTI_VaR_Norm[,4],WTI_VaR_Norm[,1])
WTI_Norm_VaR_Dur_Test5  <- VaRDurTest(0.05,WTI_VaR_Norm[,4],WTI_VaR_Norm[,2])
WTI_Norm_VaR_Dur_Test10 <- VaRDurTest(0.1,WTI_VaR_Norm[,4],WTI_VaR_Norm[,3])

WTI_JSU_VaR_Dur_Test1  <- VaRDurTest(0.01,WTI_VaR_JSU[,4],WTI_VaR_JSU[,1])
WTI_JSU_VaR_Dur_Test5  <- VaRDurTest(0.05,WTI_VaR_JSU[,4],WTI_VaR_JSU[,2])
WTI_JSU_VaR_Dur_Test10 <- VaRDurTest(0.1,WTI_VaR_JSU[,4],WTI_VaR_JSU[,3])

WTI_SGED_VaR_Dur_Test1  <- VaRDurTest(0.01,WTI_VaR_SGED[,4],WTI_VaR_SGED[,1])
WTI_SGED_VaR_Dur_Test5  <- VaRDurTest(0.05,WTI_VaR_SGED[,4],WTI_VaR_SGED[,2])
WTI_SGED_VaR_Dur_Test10 <- VaRDurTest(0.1,WTI_VaR_SGED[,4],WTI_VaR_SGED[,3])

WTI_SST_VaR_Dur_Test1  <- VaRDurTest(0.01,WTI_VaR_SST[,4],WTI_VaR_SST[,1])
WTI_SST_VaR_Dur_Test5  <- VaRDurTest(0.05,WTI_VaR_SST[,4],WTI_VaR_SST[,2])
WTI_SST_VaR_Dur_Test10 <- VaRDurTest(0.1,WTI_VaR_SST[,4],WTI_VaR_SST[,3])

Brent_Norm_VaR_Dur_Test1  <- VaRDurTest(0.01,Brent_VaR_Norm[,4],Brent_VaR_Norm[,1])
Brent_Norm_VaR_Dur_Test5  <- VaRDurTest(0.05,Brent_VaR_Norm[,4],Brent_VaR_Norm[,2])
Brent_Norm_VaR_Dur_Test10 <- VaRDurTest(0.1,Brent_VaR_Norm[,4],Brent_VaR_Norm[,3])

Brent_JSU_VaR_Dur_Test1  <- VaRDurTest(0.01,Brent_VaR_JSU[,4],Brent_VaR_JSU[,1])
Brent_JSU_VaR_Dur_Test5  <- VaRDurTest(0.05,Brent_VaR_JSU[,4],Brent_VaR_JSU[,2])
Brent_JSU_VaR_Dur_Test10 <- VaRDurTest(0.1,Brent_VaR_JSU[,4],Brent_VaR_JSU[,3])

Brent_SGED_VaR_Dur_Test1  <- VaRDurTest(0.01,Brent_VaR_SGED[,4],Brent_VaR_SGED[,1])
Brent_SGED_VaR_Dur_Test5  <- VaRDurTest(0.05,Brent_VaR_SGED[,4],Brent_VaR_SGED[,2])
Brent_SGED_VaR_Dur_Test10 <- VaRDurTest(0.1,Brent_VaR_SGED[,4],Brent_VaR_SGED[,3])

Brent_SST_VaR_Dur_Test1  <- VaRDurTest(0.01,Brent_VaR_SST[,4],Brent_VaR_SST[,1])
Brent_SST_VaR_Dur_Test5  <- VaRDurTest(0.05,Brent_VaR_SST[,4],Brent_VaR_SST[,2])
Brent_SST_VaR_Dur_Test10 <- VaRDurTest(0.1,Brent_VaR_SST[,4],Brent_VaR_SST[,3])

### Store the test results and p=values
WTI_VaR_Dur_Test_pvalues1 <- matrix(c(WTI_SST_VaR_Dur_Test1$LRp,WTI_SGED_VaR_Dur_Test1$LRp,
                                        WTI_JSU_VaR_Dur_Test1$LRp,WTI_Norm_VaR_Dur_Test1$LRp),nrow=4,byrow = TRUE)
WTI_VaR_Dur_Test_pvalues5 <- matrix(c(WTI_SST_VaR_Dur_Test5$LRp,WTI_SGED_VaR_Dur_Test5$LRp,
                                        WTI_JSU_VaR_Dur_Test5$LRp,WTI_Norm_VaR_Dur_Test5$LRp),nrow=4,byrow = TRUE)
WTI_VaR_Dur_Test_pvalues10 <- matrix(c(WTI_SST_VaR_Dur_Test10$LRp,WTI_SGED_VaR_Dur_Test10$LRp,
                                        WTI_JSU_VaR_Dur_Test10$LRp,WTI_Norm_VaR_Dur_Test10$LRp),nrow=4,byrow = TRUE)

Brent_VaR_Dur_Test_pvalues1 <- matrix(c(Brent_SST_VaR_Dur_Test1$LRp,Brent_SGED_VaR_Dur_Test1$LRp,
                                        Brent_JSU_VaR_Dur_Test1$LRp,Brent_Norm_VaR_Dur_Test1$LRp),nrow=4,byrow = TRUE)
Brent_VaR_Dur_Test_pvalues5 <- matrix(c(Brent_SST_VaR_Dur_Test5$LRp,Brent_SGED_VaR_Dur_Test5$LRp,
                                        Brent_JSU_VaR_Dur_Test5$LRp,Brent_Norm_VaR_Dur_Test5$LRp),nrow=4,byrow = TRUE)
Brent_VaR_Dur_Test_pvalues10 <- matrix(c(Brent_SST_VaR_Dur_Test10$LRp,Brent_SGED_VaR_Dur_Test10$LRp,
                                         Brent_JSU_VaR_Dur_Test10$LRp,Brent_Norm_VaR_Dur_Test10$LRp),nrow=4,byrow = TRUE)

WTI_VaR_Dur_Test_teststat1 <- matrix(c(WTI_SST_VaR_Dur_Test1$uLL,WTI_SST_VaR_Dur_Test1$rLL,
                                         WTI_SGED_VaR_Dur_Test1$uLL,WTI_SGED_VaR_Dur_Test1$rLL,
                                        WTI_JSU_VaR_Dur_Test1$uLL,WTI_JSU_VaR_Dur_Test1$rLL,
                                        WTI_Norm_VaR_Dur_Test1$uLL,WTI_Norm_VaR_Dur_Test1$rLL),nrow=4,byrow = TRUE)
WTI_VaR_Dur_Test_teststat5 <- matrix(c(WTI_SST_VaR_Dur_Test5$uLL,WTI_SST_VaR_Dur_Test5$rLL,
                                         WTI_SGED_VaR_Dur_Test5$uLL,WTI_SGED_VaR_Dur_Test5$rLL,
                                         WTI_JSU_VaR_Dur_Test5$uLL,WTI_JSU_VaR_Dur_Test5$rLL,
                                         WTI_Norm_VaR_Dur_Test5$uLL,WTI_Norm_VaR_Dur_Test5$rLL),nrow=4,byrow = TRUE)
WTI_VaR_Dur_Test_teststat10 <- matrix(c(WTI_SST_VaR_Dur_Test10$uLL,WTI_SST_VaR_Dur_Test10$rLL,
                                         WTI_SGED_VaR_Dur_Test10$uLL,WTI_SGED_VaR_Dur_Test10$rLL,
                                         WTI_JSU_VaR_Dur_Test10$uLL,WTI_JSU_VaR_Dur_Test10$rLL,
                                         WTI_Norm_VaR_Dur_Test10$uLL,WTI_Norm_VaR_Dur_Test10$rLL),nrow=4,byrow = TRUE)

Brent_VaR_Dur_Test_teststat1 <- matrix(c(Brent_SST_VaR_Dur_Test1$uLL,Brent_SST_VaR_Dur_Test1$rLL,
                                         Brent_SGED_VaR_Dur_Test1$uLL,Brent_SGED_VaR_Dur_Test1$rLL,
                                         Brent_JSU_VaR_Dur_Test1$uLL,Brent_JSU_VaR_Dur_Test1$rLL,
                                         Brent_Norm_VaR_Dur_Test1$uLL,Brent_Norm_VaR_Dur_Test1$rLL),nrow=4,byrow = TRUE)
Brent_VaR_Dur_Test_teststat5 <- matrix(c(Brent_SST_VaR_Dur_Test5$uLL,Brent_SST_VaR_Dur_Test5$rLL,
                                         Brent_SGED_VaR_Dur_Test5$uLL,Brent_SGED_VaR_Dur_Test5$rLL,
                                         Brent_JSU_VaR_Dur_Test5$uLL,Brent_JSU_VaR_Dur_Test5$rLL,
                                         Brent_Norm_VaR_Dur_Test5$uLL,Brent_Norm_VaR_Dur_Test5$rLL),nrow=4,byrow = TRUE)
Brent_VaR_Dur_Test_teststat10 <- matrix(c(Brent_SST_VaR_Dur_Test10$uLL,Brent_SST_VaR_Dur_Test10$rLL,
                                          Brent_SGED_VaR_Dur_Test10$uLL,Brent_SGED_VaR_Dur_Test10$rLL,
                                          Brent_JSU_VaR_Dur_Test10$uLL,Brent_JSU_VaR_Dur_Test10$rLL,
                                          Brent_Norm_VaR_Dur_Test10$uLL,Brent_Norm_VaR_Dur_Test10$rLL),nrow=4,byrow = TRUE)

## Round them up to 4 decimals

WTI_VaR_Dur_Test_pvalues1  <- round(WTI_VaR_Dur_Test_pvalues1,4)
WTI_VaR_Dur_Test_pvalues5  <- round(WTI_VaR_Dur_Test_pvalues5,4)
WTI_VaR_Dur_Test_pvalues10 <- round(WTI_VaR_Dur_Test_pvalues10,4)

Brent_VaR_Dur_Test_pvalues1  <- round(Brent_VaR_Dur_Test_pvalues1,4)
Brent_VaR_Dur_Test_pvalues5  <- round(Brent_VaR_Dur_Test_pvalues5,4)
Brent_VaR_Dur_Test_pvalues10 <- round(Brent_VaR_Dur_Test_pvalues10,4)

WTI_VaR_Dur_Test_teststat1  <- round(WTI_VaR_Dur_Test_teststat1,4)
WTI_VaR_Dur_Test_teststat5  <- round(WTI_VaR_Dur_Test_teststat5,4)
WTI_VaR_Dur_Test_teststat10 <- round(WTI_VaR_Dur_Test_teststat10,4)

Brent_VaR_Dur_Test_teststat1  <- round(Brent_VaR_Dur_Test_teststat1,4)
Brent_VaR_Dur_Test_teststat5  <- round(Brent_VaR_Dur_Test_teststat5,4)
Brent_VaR_Dur_Test_teststat10 <- round(Brent_VaR_Dur_Test_teststat10,4)

rownames(WTI_VaR_Dur_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
rownames(WTI_VaR_Dur_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
rownames(WTI_VaR_Dur_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")

rownames(Brent_VaR_Dur_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
rownames(Brent_VaR_Dur_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
rownames(Brent_VaR_Dur_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")

rownames(WTI_VaR_Dur_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Dur_Test_teststat1)  <- c("Unrestricted","Restricted")
rownames(WTI_VaR_Dur_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Dur_Test_teststat5)  <- c("Unrestricted","Restricted")
rownames(WTI_VaR_Dur_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(WTI_VaR_Dur_Test_teststat10)  <- c("Unrestricted","Restricted")

rownames(Brent_VaR_Dur_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Dur_Test_teststat1)  <- c("Unrestricted","Restricted")
rownames(Brent_VaR_Dur_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Dur_Test_teststat5)  <- c("Unrestricted","Restricted")
rownames(Brent_VaR_Dur_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Brent_VaR_Dur_Test_teststat10)  <- c("Unrestricted","Restricted")


## QL ratio and FZL ratio

QL_SGED_1  <- round(Brent_Test_SGED_1$Loss$Loss/Brent_Test_Norm_1$Loss$Loss,4)
QL_SGED_5  <- round(Brent_Test_SGED_5$Loss$Loss/Brent_Test_Norm_5$Loss$Loss,4)
QL_SGED_10 <- round(Brent_Test_SGED_10$Loss$Loss/Brent_Test_Norm_10$Loss$Loss,4)

QL_SST_1  <- round(Brent_Test_SST_1$Loss$Loss/Brent_Test_Norm_1$Loss$Loss,4)
QL_SST_5  <- round(Brent_Test_SST_5$Loss$Loss/Brent_Test_Norm_5$Loss$Loss,4)
QL_SST_10 <- round(Brent_Test_SST_10$Loss$Loss/Brent_Test_Norm_10$Loss$Loss,4)

QL_JSU_1  <- round(Brent_Test_JSU_1$Loss$Loss/Brent_Test_Norm_1$Loss$Loss,4)
QL_JSU_5  <- round(Brent_Test_JSU_5$Loss$Loss/Brent_Test_Norm_5$Loss$Loss,4)
QL_JSU_10 <- round(Brent_Test_JSU_10$Loss$Loss/Brent_Test_Norm_10$Loss$Loss,4)

QL <- matrix(c(QL_SST_1,QL_SST_5,QL_SST_10,QL_SGED_1,QL_SGED_5,QL_SGED_10,QL_JSU_1,QL_JSU_5,QL_JSU_10),nrow=3,byrow=TRUE)
rownames(QL) <- c("SST","SGED","JSU")
colnames(QL) <- c("1%","5%","10%")

AE_SST_1  <- round(Brent_Test_SST_1$AE,4)
AE_SST_5  <- round(Brent_Test_SST_5$AE,4)
AE_SST_10 <- round(Brent_Test_SST_10$AE,4)

AE_SGED_1  <- round(Brent_Test_SGED_1$AE,4)
AE_SGED_5  <- round(Brent_Test_SGED_5$AE,4)
AE_SGED_10 <- round(Brent_Test_SGED_10$AE,4)

AE_JSU_1  <- round(Brent_Test_JSU_1$AE,4)
AE_JSU_5  <- round(Brent_Test_JSU_5$AE,4)
AE_JSU_10 <- round(Brent_Test_JSU_10$AE,4)

AE_Norm_1  <- round(Brent_Test_Norm_1$AE,4)
AE_Norm_5  <- round(Brent_Test_Norm_5$AE,4)
AE_Norm_10 <- round(Brent_Test_Norm_10$AE,4)

AE <- matrix(c(AE_SST_1,AE_SST_5,AE_SST_10,AE_SGED_1,AE_SGED_5,AE_SGED_10,AE_JSU_1,AE_JSU_5,AE_JSU_10,AE_Norm_1,AE_Norm_5,AE_Norm_10),nrow=4,byrow=TRUE)
rownames(AE) <- c("SST","SGED","JSU","Norm")
colnames(AE) <- c("1%","5%","10%")

### Check the summary stats of the standardized residuals

WTI_Res8 <- WTI_Model8@fit$residuals/WTI_Model8@fit$sigma
Brent_Res8   <- Brent_Model8@fit$residuals/Brent_Model8@fit$sigma

WTI_Descr8 <- basicStats(WTI_Res8)
Brent_Descr8 <- basicStats(Brent_Res8)

### Fitting distributions on the standardized residuals 

WTI_Dist8 <- pearsonFitML(WTI_Res8)
Brent_Dist8   <- pearsonFitML(Brent_Res8)
WTI_JSU_param <- JohnsonFit(WTI_Res8,moment ="find")
Brent_JSU_param <- JohnsonFit(Brent_Res8,moment ="find")

### Store the parameters for different distributions  

WTI_P4_Param  <- WTI_Dist8[2:5]
Brent_P4_Param    <- Brent_Dist8[2:5]

### Checking the fit using KS tests

WTI_test1 <- ks.test(WTI_Res8,rpearson(length(WTI_Res8),params = WTI_Dist8))
Brent_test1 <- ks.test(Brent_Res8,rpearson(length(Brent_Res8),params = Brent_Dist8))


### Calculate the quantiles

WTI_Exceed_Pearson  <- array(data=0,dim=length(signific_levels))
Brent_Exceed_Pearson    <- array(data=0,dim=length(signific_levels))

for(j in 1:3)
{
  Z_Pearson_WTI <- qpearson(signific_levels[j],WTI_Dist8)
  Z_Pearson_Brent   <- qpearson(signific_levels[j],Brent_Dist8)
  
  mu_WTI   <- array(data=0,dim=length(Returns_WTI)-1000)
  VaR_WTI   <- array(data=0,dim=length(Returns_WTI)-1001)
  mu_Brent   <- array(data=0,dim=length(Returns_Brent)-1000)
  VaR_Brent   <- array(data=0,dim=length(Returns_Brent)-1001)
  
  for(i in 2:(length(Returns_WTI)-1000))
  {
    mu_WTI[i]    <- WTI_Model8@fit$coef[1]+WTI_Model8@fit$coef[2]*Returns_WTI[i-1]
    VaR_WTI[i-1] <- mu_WTI[i]+WTI_Model8@fit$sigma[i]*Z_Pearson_WTI
  }

  for(i in 2:(length(Returns_Brent)-1000))
  {
    mu_Brent[i] <- Brent_Model8@fit$coef[1]+Brent_Model8@fit$coef[2]*Returns_Brent[i-1]
    VaR_Brent[i-1] <- mu_Brent[i]+Brent_Model8@fit$sigma[i]*Z_Pearson_Brent
  }
  
  WTI_Violations <- VaRTest(signific_levels[j],Returns_WTI[1:(length(VaR_WTI))],VaR_WTI)
  Brent_Violations   <- VaRTest(signific_levels[j],Returns_Brent[1:(length(VaR_Brent))],VaR_Brent)
  WTI_Exceed_Pearson[j] <- WTI_Violations$actual.exceed
  Brent_Exceed_Pearson[j] <- Brent_Violations$actual.exceed
  
}


#### Out of sample VaR forecast  #####

#### Using Hold out sample

Holdout_WTI_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Brent_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_WTI_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Brent_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_WTI_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Brent_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_WTI_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Brent_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_WTI_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Brent_Pearson) <- c("1%","5%","10%")
rownames(Holdout_WTI_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Brent_Pearson_pvalues) <- c("1%","5%","10%")

## Store the volatility and mu forecast
#signific_levels <- c(0.1,0.05,0.01)

WTI_Holdout_Mu <- as.data.frame(WTI_Roll8)[,'Mu']
WTI_Holdout_Sigma <- as.data.frame(WTI_Roll8)[,'Sigma']
Brent_Holdout_Mu <- as.data.frame(Brent_Roll8)[,'Mu']
Brent_Holdout_Sigma <- as.data.frame(Brent_Roll8)[,'Sigma']

for(j in 1:3)
{
 
  Z_Pearson_WTI <- qpearson(signific_levels[j],WTI_Dist8)
  Z_Pearson_Brent   <- qpearson(signific_levels[j],Brent_Dist8)

  WTI_Holdout_VaR <- WTI_Holdout_Mu + WTI_Holdout_Sigma*Z_Pearson_WTI
  Brent_Holdout_VaR   <- Brent_Holdout_Mu + Brent_Holdout_Sigma*Z_Pearson_Brent
  WTI_Violations <- BacktestVaR(Holdout_Returns_WTI[1:(length(Holdout_Returns_WTI)-1)],WTI_Holdout_VaR[-1],signific_levels[j])
  Brent_Violations   <- BacktestVaR(Holdout_Returns_Brent[1:(length(Holdout_Returns_Brent)-1)],Brent_Holdout_VaR[-1],signific_levels[j])
  Holdout_WTI_Pearson[j,] <- c(WTI_Violations$LRuc[1],WTI_Violations$LRcc[1],WTI_Violations$DQ$stat,WTI_Violations$AE,WTI_Violations$Loss$Loss)
  Holdout_Brent_Pearson[j,]   <- c(Brent_Violations$LRuc[1],Brent_Violations$LRcc[1],Brent_Violations$DQ$stat,Brent_Violations$AE,Brent_Violations$Loss$Loss)
  Holdout_WTI_Pearson_pvalues[j,] <- c(WTI_Violations$LRuc[2],WTI_Violations$LRcc[2],WTI_Violations$DQ$pvalue)
  Holdout_Brent_Pearson_pvalues[j,] <- c(Brent_Violations$LRuc[2],Brent_Violations$LRcc[2],Brent_Violations$DQ$pvalue)
}


#### Subsample Analysis ####

# ### JSU specification

WTI_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

WTI_Model_Sample_JSU1 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub1,out.sample = 250)
WTI_Model_Sample_JSU2 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub2,out.sample = 250)
WTI_Model_Sample_JSU3 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub3,out.sample = 250)
WTI_Model_Sample_JSU4 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub4,out.sample = 250)
WTI_Model_Sample_JSU5 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub5,out.sample = 250)

Brent_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

Brent_Model_Sample_JSU1 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub1,out.sample = 250)
Brent_Model_Sample_JSU2 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub2,out.sample = 250)
Brent_Model_Sample_JSU3 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub3,out.sample = 250)
Brent_Model_Sample_JSU4 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub4,out.sample = 250)
Brent_Model_Sample_JSU5 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub5,out.sample = 250)


## Rollover Analysis

WTI_Roll_JSU_Sample1 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample2 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample3 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample4 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample5 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_JSU_Sample1 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample2 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample3 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample4 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample5 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### SGED specification

WTI_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

WTI_Model_Sample_SGED1 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub1,out.sample = 250)
WTI_Model_Sample_SGED2 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub2,out.sample = 250)
WTI_Model_Sample_SGED3 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub3,out.sample = 250)
WTI_Model_Sample_SGED4 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub4,out.sample = 250)
WTI_Model_Sample_SGED5 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub5,out.sample = 250)

Brent_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

Brent_Model_Sample_SGED1 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub1,out.sample = 250)
Brent_Model_Sample_SGED2 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub2,out.sample = 250)
Brent_Model_Sample_SGED3 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub3,out.sample = 250)
Brent_Model_Sample_SGED4 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub4,out.sample = 250)
Brent_Model_Sample_SGED5 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub5,out.sample = 250)


## Rollover Analysis

WTI_Roll_SGED_Sample1 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample2 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample3 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample4 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample5 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_SGED_Sample1 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample2 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample3 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample4 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample5 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### SST Specification
# 
WTI_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
# 
WTI_Model_Sample_SST1 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub1,out.sample = 250)
WTI_Model_Sample_SST2 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub2,out.sample = 250)
WTI_Model_Sample_SST3 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub3,out.sample = 250)
WTI_Model_Sample_SST4 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub4,out.sample = 250)
WTI_Model_Sample_SST5 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub5,out.sample = 250)

Brent_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")

Brent_Model_Sample_SST1 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub1,out.sample = 250)
Brent_Model_Sample_SST2 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub2,out.sample = 250)
Brent_Model_Sample_SST3 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub3,out.sample = 250)
Brent_Model_Sample_SST4 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub4,out.sample = 250)
Brent_Model_Sample_SST5 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub5,out.sample = 250)


## Rollover Analysis

WTI_Roll_SST_Sample1 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample2 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample3 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample4 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample5 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_SST_Sample1 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample2 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample3 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample4 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample5 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### Normal Specification

WTI_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model_Sample_Normal1 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)
WTI_Model_Sample_Normal2 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)
WTI_Model_Sample_Normal3 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)
WTI_Model_Sample_Normal4 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)
WTI_Model_Sample_Normal5 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)

Brent_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model_Sample_Normal1 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)
Brent_Model_Sample_Normal2 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)
Brent_Model_Sample_Normal3 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)
Brent_Model_Sample_Normal4 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)
Brent_Model_Sample_Normal5 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)

WTI_Roll_Normal_Sample1 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample2 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample3 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample4 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample5 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_Normal_Sample1 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample2 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample3 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample4 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample5 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### Sample Residuals Statistics

WTI_Res_Sample1 <- c(normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample2 <- c(normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample3 <- c(normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample4 <- c(normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample5 <- c(normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)

Brent_Res_Sample1 <- c(normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample2 <- c(normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample3 <- c(normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample4 <- c(normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample5 <- c(normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)

WTI_Res_Sample_Test <- rbind(WTI_Res_Sample1,WTI_Res_Sample2,WTI_Res_Sample3,WTI_Res_Sample4,WTI_Res_Sample5)
Brent_Res_Sample_Test   <- rbind(Brent_Res_Sample1,Brent_Res_Sample2,Brent_Res_Sample3,Brent_Res_Sample4,Brent_Res_Sample5)
rownames(WTI_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")
colnames(WTI_Res_Sample_Test) <- c("JB","SW","KS") 
rownames(Brent_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")
colnames(Brent_Res_Sample_Test) <- c("JB","SW","KS") 

### Pearson Specification

WTI_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$sigma)

WTI_Param_Pearson_Subsample <- rbind(WTI_Param_Pearson_Sample1,WTI_Param_Pearson_Sample2,WTI_Param_Pearson_Sample3,WTI_Param_Pearson_Sample4,WTI_Param_Pearson_Sample5)
rownames(WTI_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")

Brent_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$sigma)

Brent_Param_Pearson_Subsample <- rbind(Brent_Param_Pearson_Sample1,Brent_Param_Pearson_Sample2,Brent_Param_Pearson_Sample3,Brent_Param_Pearson_Sample4,Brent_Param_Pearson_Sample5)
rownames(Brent_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")

### johnson SU parameter


WTI_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$sigma,moment = "find")

WTI_Param_Johnson_Subsample <- rbind(WTI_Param_Johnson_Sample1,WTI_Param_Johnson_Sample2,WTI_Param_Johnson_Sample3,WTI_Param_Johnson_Sample4,WTI_Param_Johnson_Sample5)
rownames(WTI_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")

Brent_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$sigma,moment = "find")

Brent_Param_Johnson_Subsample <- rbind(Brent_Param_Johnson_Sample1,Brent_Param_Johnson_Sample2,Brent_Param_Johnson_Sample3,Brent_Param_Johnson_Sample4,Brent_Param_Johnson_Sample5)
rownames(Brent_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")


#### Backtesting and Storing the results

WTI_Test_SST1_Sample1 <- BacktestVaR(WTI_Roll_SST_Sample1@forecast$VaR[,4],WTI_Roll_SST_Sample1@forecast$VaR[,1],0.01)
WTI_Test_SST5_Sample1 <- BacktestVaR(WTI_Roll_SST_Sample1@forecast$VaR[,4],WTI_Roll_SST_Sample1@forecast$VaR[,2],0.05)
WTI_Test_SST10_Sample1 <- BacktestVaR(WTI_Roll_SST_Sample1@forecast$VaR[,4],WTI_Roll_SST_Sample1@forecast$VaR[,3],0.1)

Brent_Test_SST1_Sample1 <- BacktestVaR(Brent_Roll_SST_Sample1@forecast$VaR[,4],Brent_Roll_SST_Sample1@forecast$VaR[,1],0.01)
Brent_Test_SST5_Sample1 <- BacktestVaR(Brent_Roll_SST_Sample1@forecast$VaR[,4],Brent_Roll_SST_Sample1@forecast$VaR[,2],0.05)
Brent_Test_SST10_Sample1 <- BacktestVaR(Brent_Roll_SST_Sample1@forecast$VaR[,4],Brent_Roll_SST_Sample1@forecast$VaR[,3],0.1)

WTI_Test_SGED1_Sample1 <- BacktestVaR(WTI_Roll_SGED_Sample1@forecast$VaR[,4],WTI_Roll_SGED_Sample1@forecast$VaR[,1],0.01)
WTI_Test_SGED5_Sample1 <- BacktestVaR(WTI_Roll_SGED_Sample1@forecast$VaR[,4],WTI_Roll_SGED_Sample1@forecast$VaR[,2],0.05)
WTI_Test_SGED10_Sample1 <- BacktestVaR(WTI_Roll_SGED_Sample1@forecast$VaR[,4],WTI_Roll_SGED_Sample1@forecast$VaR[,3],0.1)

Brent_Test_SGED1_Sample1 <- BacktestVaR(Brent_Roll_SGED_Sample1@forecast$VaR[,4],Brent_Roll_SGED_Sample1@forecast$VaR[,1],0.01)
Brent_Test_SGED5_Sample1 <- BacktestVaR(Brent_Roll_SGED_Sample1@forecast$VaR[,4],Brent_Roll_SGED_Sample1@forecast$VaR[,2],0.05)
Brent_Test_SGED10_Sample1 <- BacktestVaR(Brent_Roll_SGED_Sample1@forecast$VaR[,4],Brent_Roll_SGED_Sample1@forecast$VaR[,3],0.1)

WTI_Test_JSU1_Sample1 <- BacktestVaR(WTI_Roll_JSU_Sample1@forecast$VaR[,4],WTI_Roll_JSU_Sample1@forecast$VaR[,1],0.01)
WTI_Test_JSU5_Sample1 <- BacktestVaR(WTI_Roll_JSU_Sample1@forecast$VaR[,4],WTI_Roll_JSU_Sample1@forecast$VaR[,2],0.05)
WTI_Test_JSU10_Sample1 <- BacktestVaR(WTI_Roll_JSU_Sample1@forecast$VaR[,4],WTI_Roll_JSU_Sample1@forecast$VaR[,3],0.1)

Brent_Test_JSU1_Sample1 <- BacktestVaR(Brent_Roll_JSU_Sample1@forecast$VaR[,4],Brent_Roll_JSU_Sample1@forecast$VaR[,1],0.01)
Brent_Test_JSU5_Sample1 <- BacktestVaR(Brent_Roll_JSU_Sample1@forecast$VaR[,4],Brent_Roll_JSU_Sample1@forecast$VaR[,2],0.05)
Brent_Test_JSU10_Sample1 <- BacktestVaR(Brent_Roll_JSU_Sample1@forecast$VaR[,4],Brent_Roll_JSU_Sample1@forecast$VaR[,3],0.1)

WTI_Test_Normal1_Sample1 <- BacktestVaR(WTI_Roll_Normal_Sample1@forecast$VaR[,4],WTI_Roll_Normal_Sample1@forecast$VaR[,1],0.01)
WTI_Test_Normal5_Sample1 <- BacktestVaR(WTI_Roll_Normal_Sample1@forecast$VaR[,4],WTI_Roll_Normal_Sample1@forecast$VaR[,2],0.05)
WTI_Test_Normal10_Sample1 <- BacktestVaR(WTI_Roll_Normal_Sample1@forecast$VaR[,4],WTI_Roll_Normal_Sample1@forecast$VaR[,3],0.1)

Brent_Test_Normal1_Sample1 <- BacktestVaR(Brent_Roll_Normal_Sample1@forecast$VaR[,4],Brent_Roll_Normal_Sample1@forecast$VaR[,1],0.01)
Brent_Test_Normal5_Sample1 <- BacktestVaR(Brent_Roll_Normal_Sample1@forecast$VaR[,4],Brent_Roll_Normal_Sample1@forecast$VaR[,2],0.05)
Brent_Test_Normal10_Sample1 <- BacktestVaR(Brent_Roll_Normal_Sample1@forecast$VaR[,4],Brent_Roll_Normal_Sample1@forecast$VaR[,3],0.1)

##### Sample 2


WTI_Test_SST1_Sample2 <- BacktestVaR(WTI_Roll_SST_Sample2@forecast$VaR[,4],WTI_Roll_SST_Sample2@forecast$VaR[,1],0.01)
WTI_Test_SST5_Sample2 <- BacktestVaR(WTI_Roll_SST_Sample2@forecast$VaR[,4],WTI_Roll_SST_Sample2@forecast$VaR[,2],0.05)
WTI_Test_SST10_Sample2 <- BacktestVaR(WTI_Roll_SST_Sample2@forecast$VaR[,4],WTI_Roll_SST_Sample2@forecast$VaR[,3],0.1)

Brent_Test_SST1_Sample2 <- BacktestVaR(Brent_Roll_SST_Sample2@forecast$VaR[,4],Brent_Roll_SST_Sample2@forecast$VaR[,1],0.01)
Brent_Test_SST5_Sample2 <- BacktestVaR(Brent_Roll_SST_Sample2@forecast$VaR[,4],Brent_Roll_SST_Sample2@forecast$VaR[,2],0.05)
Brent_Test_SST10_Sample2 <- BacktestVaR(Brent_Roll_SST_Sample2@forecast$VaR[,4],Brent_Roll_SST_Sample2@forecast$VaR[,3],0.1)

WTI_Test_SGED1_Sample2 <- BacktestVaR(WTI_Roll_SGED_Sample2@forecast$VaR[,4],WTI_Roll_SGED_Sample2@forecast$VaR[,1],0.01)
WTI_Test_SGED5_Sample2 <- BacktestVaR(WTI_Roll_SGED_Sample2@forecast$VaR[,4],WTI_Roll_SGED_Sample2@forecast$VaR[,2],0.05)
WTI_Test_SGED10_Sample2 <- BacktestVaR(WTI_Roll_SGED_Sample2@forecast$VaR[,4],WTI_Roll_SGED_Sample2@forecast$VaR[,3],0.1)

Brent_Test_SGED1_Sample2 <- BacktestVaR(Brent_Roll_SGED_Sample2@forecast$VaR[,4],Brent_Roll_SGED_Sample2@forecast$VaR[,1],0.01)
Brent_Test_SGED5_Sample2 <- BacktestVaR(Brent_Roll_SGED_Sample2@forecast$VaR[,4],Brent_Roll_SGED_Sample2@forecast$VaR[,2],0.05)
Brent_Test_SGED10_Sample2 <- BacktestVaR(Brent_Roll_SGED_Sample2@forecast$VaR[,4],Brent_Roll_SGED_Sample2@forecast$VaR[,3],0.1)

WTI_Test_JSU1_Sample2 <- BacktestVaR(WTI_Roll_JSU_Sample2@forecast$VaR[,4],WTI_Roll_JSU_Sample2@forecast$VaR[,1],0.01)
WTI_Test_JSU5_Sample2 <- BacktestVaR(WTI_Roll_JSU_Sample2@forecast$VaR[,4],WTI_Roll_JSU_Sample2@forecast$VaR[,2],0.05)
WTI_Test_JSU10_Sample2 <- BacktestVaR(WTI_Roll_JSU_Sample2@forecast$VaR[,4],WTI_Roll_JSU_Sample2@forecast$VaR[,3],0.1)

Brent_Test_JSU1_Sample2 <- BacktestVaR(Brent_Roll_JSU_Sample2@forecast$VaR[,4],Brent_Roll_JSU_Sample2@forecast$VaR[,1],0.01)
Brent_Test_JSU5_Sample2 <- BacktestVaR(Brent_Roll_JSU_Sample2@forecast$VaR[,4],Brent_Roll_JSU_Sample2@forecast$VaR[,2],0.05)
Brent_Test_JSU10_Sample2 <- BacktestVaR(Brent_Roll_JSU_Sample2@forecast$VaR[,4],Brent_Roll_JSU_Sample2@forecast$VaR[,3],0.1)

WTI_Test_Normal1_Sample2 <- BacktestVaR(WTI_Roll_Normal_Sample2@forecast$VaR[,4],WTI_Roll_Normal_Sample2@forecast$VaR[,1],0.01)
WTI_Test_Normal5_Sample2 <- BacktestVaR(WTI_Roll_Normal_Sample2@forecast$VaR[,4],WTI_Roll_Normal_Sample2@forecast$VaR[,2],0.05)
WTI_Test_Normal10_Sample2 <- BacktestVaR(WTI_Roll_Normal_Sample2@forecast$VaR[,4],WTI_Roll_Normal_Sample2@forecast$VaR[,3],0.1)

Brent_Test_Normal1_Sample2 <- BacktestVaR(Brent_Roll_Normal_Sample2@forecast$VaR[,4],Brent_Roll_Normal_Sample2@forecast$VaR[,1],0.01)
Brent_Test_Normal5_Sample2 <- BacktestVaR(Brent_Roll_Normal_Sample2@forecast$VaR[,4],Brent_Roll_Normal_Sample2@forecast$VaR[,2],0.05)
Brent_Test_Normal10_Sample2 <- BacktestVaR(Brent_Roll_Normal_Sample2@forecast$VaR[,4],Brent_Roll_Normal_Sample2@forecast$VaR[,3],0.1)

### Sample3


WTI_Test_SST1_Sample3 <- BacktestVaR(WTI_Roll_SST_Sample3@forecast$VaR[,4],WTI_Roll_SST_Sample3@forecast$VaR[,1],0.01)
WTI_Test_SST5_Sample3 <- BacktestVaR(WTI_Roll_SST_Sample3@forecast$VaR[,4],WTI_Roll_SST_Sample3@forecast$VaR[,2],0.05)
WTI_Test_SST10_Sample3 <- BacktestVaR(WTI_Roll_SST_Sample3@forecast$VaR[,4],WTI_Roll_SST_Sample3@forecast$VaR[,3],0.1)

Brent_Test_SST1_Sample3 <- BacktestVaR(Brent_Roll_SST_Sample3@forecast$VaR[,4],Brent_Roll_SST_Sample3@forecast$VaR[,1],0.01)
Brent_Test_SST5_Sample3 <- BacktestVaR(Brent_Roll_SST_Sample3@forecast$VaR[,4],Brent_Roll_SST_Sample3@forecast$VaR[,2],0.05)
Brent_Test_SST10_Sample3 <- BacktestVaR(Brent_Roll_SST_Sample3@forecast$VaR[,4],Brent_Roll_SST_Sample3@forecast$VaR[,3],0.1)

WTI_Test_SGED1_Sample3 <- BacktestVaR(WTI_Roll_SGED_Sample3@forecast$VaR[,4],WTI_Roll_SGED_Sample3@forecast$VaR[,1],0.01)
WTI_Test_SGED5_Sample3 <- BacktestVaR(WTI_Roll_SGED_Sample3@forecast$VaR[,4],WTI_Roll_SGED_Sample3@forecast$VaR[,2],0.05)
WTI_Test_SGED10_Sample3 <- BacktestVaR(WTI_Roll_SGED_Sample3@forecast$VaR[,4],WTI_Roll_SGED_Sample3@forecast$VaR[,3],0.1)

Brent_Test_SGED1_Sample3 <- BacktestVaR(Brent_Roll_SGED_Sample3@forecast$VaR[,4],Brent_Roll_SGED_Sample3@forecast$VaR[,1],0.01)
Brent_Test_SGED5_Sample3 <- BacktestVaR(Brent_Roll_SGED_Sample3@forecast$VaR[,4],Brent_Roll_SGED_Sample3@forecast$VaR[,2],0.05)
Brent_Test_SGED10_Sample3 <- BacktestVaR(Brent_Roll_SGED_Sample3@forecast$VaR[,4],Brent_Roll_SGED_Sample3@forecast$VaR[,3],0.1)

WTI_Test_JSU1_Sample3 <- BacktestVaR(WTI_Roll_JSU_Sample3@forecast$VaR[,4],WTI_Roll_JSU_Sample3@forecast$VaR[,1],0.01)
WTI_Test_JSU5_Sample3 <- BacktestVaR(WTI_Roll_JSU_Sample3@forecast$VaR[,4],WTI_Roll_JSU_Sample3@forecast$VaR[,2],0.05)
WTI_Test_JSU10_Sample3 <- BacktestVaR(WTI_Roll_JSU_Sample3@forecast$VaR[,4],WTI_Roll_JSU_Sample3@forecast$VaR[,3],0.1)

Brent_Test_JSU1_Sample3 <- BacktestVaR(Brent_Roll_JSU_Sample3@forecast$VaR[,4],Brent_Roll_JSU_Sample3@forecast$VaR[,1],0.01)
Brent_Test_JSU5_Sample3 <- BacktestVaR(Brent_Roll_JSU_Sample3@forecast$VaR[,4],Brent_Roll_JSU_Sample3@forecast$VaR[,2],0.05)
Brent_Test_JSU10_Sample3 <- BacktestVaR(Brent_Roll_JSU_Sample3@forecast$VaR[,4],Brent_Roll_JSU_Sample3@forecast$VaR[,3],0.1)

WTI_Test_Normal1_Sample3 <- BacktestVaR(WTI_Roll_Normal_Sample3@forecast$VaR[,4],WTI_Roll_Normal_Sample3@forecast$VaR[,1],0.01)
WTI_Test_Normal5_Sample3 <- BacktestVaR(WTI_Roll_Normal_Sample3@forecast$VaR[,4],WTI_Roll_Normal_Sample3@forecast$VaR[,2],0.05)
WTI_Test_Normal10_Sample3 <- BacktestVaR(WTI_Roll_Normal_Sample3@forecast$VaR[,4],WTI_Roll_Normal_Sample3@forecast$VaR[,3],0.1)

Brent_Test_Normal1_Sample3 <- BacktestVaR(Brent_Roll_Normal_Sample3@forecast$VaR[,4],Brent_Roll_Normal_Sample3@forecast$VaR[,1],0.01)
Brent_Test_Normal5_Sample3 <- BacktestVaR(Brent_Roll_Normal_Sample3@forecast$VaR[,4],Brent_Roll_Normal_Sample3@forecast$VaR[,2],0.05)
Brent_Test_Normal10_Sample3 <- BacktestVaR(Brent_Roll_Normal_Sample3@forecast$VaR[,4],Brent_Roll_Normal_Sample3@forecast$VaR[,3],0.1)

### Sample 4


WTI_Test_SST1_Sample4 <- BacktestVaR(WTI_Roll_SST_Sample4@forecast$VaR[,4],WTI_Roll_SST_Sample4@forecast$VaR[,1],0.01)
WTI_Test_SST5_Sample4 <- BacktestVaR(WTI_Roll_SST_Sample4@forecast$VaR[,4],WTI_Roll_SST_Sample4@forecast$VaR[,2],0.05)
WTI_Test_SST10_Sample4 <- BacktestVaR(WTI_Roll_SST_Sample4@forecast$VaR[,4],WTI_Roll_SST_Sample4@forecast$VaR[,3],0.1)

Brent_Test_SST1_Sample4 <- BacktestVaR(Brent_Roll_SST_Sample4@forecast$VaR[,4],Brent_Roll_SST_Sample4@forecast$VaR[,1],0.01)
Brent_Test_SST5_Sample4 <- BacktestVaR(Brent_Roll_SST_Sample4@forecast$VaR[,4],Brent_Roll_SST_Sample4@forecast$VaR[,2],0.05)
Brent_Test_SST10_Sample4 <- BacktestVaR(Brent_Roll_SST_Sample4@forecast$VaR[,4],Brent_Roll_SST_Sample4@forecast$VaR[,3],0.1)

WTI_Test_SGED1_Sample4 <- BacktestVaR(WTI_Roll_SGED_Sample4@forecast$VaR[,4],WTI_Roll_SGED_Sample4@forecast$VaR[,1],0.01)
WTI_Test_SGED5_Sample4 <- BacktestVaR(WTI_Roll_SGED_Sample4@forecast$VaR[,4],WTI_Roll_SGED_Sample4@forecast$VaR[,2],0.05)
WTI_Test_SGED10_Sample4 <- BacktestVaR(WTI_Roll_SGED_Sample4@forecast$VaR[,4],WTI_Roll_SGED_Sample4@forecast$VaR[,3],0.1)

Brent_Test_SGED1_Sample4 <- BacktestVaR(Brent_Roll_SGED_Sample4@forecast$VaR[,4],Brent_Roll_SGED_Sample4@forecast$VaR[,1],0.01)
Brent_Test_SGED5_Sample4 <- BacktestVaR(Brent_Roll_SGED_Sample4@forecast$VaR[,4],Brent_Roll_SGED_Sample4@forecast$VaR[,2],0.05)
Brent_Test_SGED10_Sample4 <- BacktestVaR(Brent_Roll_SGED_Sample4@forecast$VaR[,4],Brent_Roll_SGED_Sample4@forecast$VaR[,3],0.1)

WTI_Test_JSU1_Sample4 <- BacktestVaR(WTI_Roll_JSU_Sample4@forecast$VaR[,4],WTI_Roll_JSU_Sample4@forecast$VaR[,1],0.01)
WTI_Test_JSU5_Sample4 <- BacktestVaR(WTI_Roll_JSU_Sample4@forecast$VaR[,4],WTI_Roll_JSU_Sample4@forecast$VaR[,2],0.05)
WTI_Test_JSU10_Sample4 <- BacktestVaR(WTI_Roll_JSU_Sample4@forecast$VaR[,4],WTI_Roll_JSU_Sample4@forecast$VaR[,3],0.1)

Brent_Test_JSU1_Sample4 <- BacktestVaR(Brent_Roll_JSU_Sample4@forecast$VaR[,4],Brent_Roll_JSU_Sample4@forecast$VaR[,1],0.01)
Brent_Test_JSU5_Sample4 <- BacktestVaR(Brent_Roll_JSU_Sample4@forecast$VaR[,4],Brent_Roll_JSU_Sample4@forecast$VaR[,2],0.05)
Brent_Test_JSU10_Sample4 <- BacktestVaR(Brent_Roll_JSU_Sample4@forecast$VaR[,4],Brent_Roll_JSU_Sample4@forecast$VaR[,3],0.1)

WTI_Test_Normal1_Sample4 <- BacktestVaR(WTI_Roll_Normal_Sample4@forecast$VaR[,4],WTI_Roll_Normal_Sample4@forecast$VaR[,1],0.01)
WTI_Test_Normal5_Sample4 <- BacktestVaR(WTI_Roll_Normal_Sample4@forecast$VaR[,4],WTI_Roll_Normal_Sample4@forecast$VaR[,2],0.05)
WTI_Test_Normal10_Sample4 <- BacktestVaR(WTI_Roll_Normal_Sample4@forecast$VaR[,4],WTI_Roll_Normal_Sample4@forecast$VaR[,3],0.1)

Brent_Test_Normal1_Sample4 <- BacktestVaR(Brent_Roll_Normal_Sample4@forecast$VaR[,4],Brent_Roll_Normal_Sample4@forecast$VaR[,1],0.01)
Brent_Test_Normal5_Sample4 <- BacktestVaR(Brent_Roll_Normal_Sample4@forecast$VaR[,4],Brent_Roll_Normal_Sample4@forecast$VaR[,2],0.05)
Brent_Test_Normal10_Sample4 <- BacktestVaR(Brent_Roll_Normal_Sample4@forecast$VaR[,4],Brent_Roll_Normal_Sample4@forecast$VaR[,3],0.1)

## Sample 5


WTI_Test_SST1_Sample5 <- BacktestVaR(WTI_Roll_SST_Sample5@forecast$VaR[,4],WTI_Roll_SST_Sample5@forecast$VaR[,1],0.01)
WTI_Test_SST5_Sample5 <- BacktestVaR(WTI_Roll_SST_Sample5@forecast$VaR[,4],WTI_Roll_SST_Sample5@forecast$VaR[,2],0.05)
WTI_Test_SST10_Sample5 <- BacktestVaR(WTI_Roll_SST_Sample5@forecast$VaR[,4],WTI_Roll_SST_Sample5@forecast$VaR[,3],0.1)

Brent_Test_SST1_Sample5 <- BacktestVaR(Brent_Roll_SST_Sample5@forecast$VaR[,4],Brent_Roll_SST_Sample5@forecast$VaR[,1],0.01)
Brent_Test_SST5_Sample5 <- BacktestVaR(Brent_Roll_SST_Sample5@forecast$VaR[,4],Brent_Roll_SST_Sample5@forecast$VaR[,2],0.05)
Brent_Test_SST10_Sample5 <- BacktestVaR(Brent_Roll_SST_Sample5@forecast$VaR[,4],Brent_Roll_SST_Sample5@forecast$VaR[,3],0.1)

WTI_Test_SGED1_Sample5 <- BacktestVaR(WTI_Roll_SGED_Sample5@forecast$VaR[,4],WTI_Roll_SGED_Sample5@forecast$VaR[,1],0.01)
WTI_Test_SGED5_Sample5 <- BacktestVaR(WTI_Roll_SGED_Sample5@forecast$VaR[,4],WTI_Roll_SGED_Sample5@forecast$VaR[,2],0.05)
WTI_Test_SGED10_Sample5 <- BacktestVaR(WTI_Roll_SGED_Sample5@forecast$VaR[,4],WTI_Roll_SGED_Sample5@forecast$VaR[,3],0.1)

Brent_Test_SGED1_Sample5 <- BacktestVaR(Brent_Roll_SGED_Sample5@forecast$VaR[,4],Brent_Roll_SGED_Sample5@forecast$VaR[,1],0.01)
Brent_Test_SGED5_Sample5 <- BacktestVaR(Brent_Roll_SGED_Sample5@forecast$VaR[,4],Brent_Roll_SGED_Sample5@forecast$VaR[,2],0.05)
Brent_Test_SGED10_Sample5 <- BacktestVaR(Brent_Roll_SGED_Sample5@forecast$VaR[,4],Brent_Roll_SGED_Sample5@forecast$VaR[,3],0.1)

WTI_Test_JSU1_Sample5 <- BacktestVaR(WTI_Roll_JSU_Sample5@forecast$VaR[,4],WTI_Roll_JSU_Sample5@forecast$VaR[,1],0.01)
WTI_Test_JSU5_Sample5 <- BacktestVaR(WTI_Roll_JSU_Sample5@forecast$VaR[,4],WTI_Roll_JSU_Sample5@forecast$VaR[,2],0.05)
WTI_Test_JSU10_Sample5 <- BacktestVaR(WTI_Roll_JSU_Sample5@forecast$VaR[,4],WTI_Roll_JSU_Sample5@forecast$VaR[,3],0.1)

Brent_Test_JSU1_Sample5 <- BacktestVaR(Brent_Roll_JSU_Sample5@forecast$VaR[,4],Brent_Roll_JSU_Sample5@forecast$VaR[,1],0.01)
Brent_Test_JSU5_Sample5 <- BacktestVaR(Brent_Roll_JSU_Sample5@forecast$VaR[,4],Brent_Roll_JSU_Sample5@forecast$VaR[,2],0.05)
Brent_Test_JSU10_Sample5 <- BacktestVaR(Brent_Roll_JSU_Sample5@forecast$VaR[,4],Brent_Roll_JSU_Sample5@forecast$VaR[,3],0.1)

WTI_Test_Normal1_Sample5 <- BacktestVaR(WTI_Roll_Normal_Sample5@forecast$VaR[,4],WTI_Roll_Normal_Sample5@forecast$VaR[,1],0.01)
WTI_Test_Normal5_Sample5 <- BacktestVaR(WTI_Roll_Normal_Sample5@forecast$VaR[,4],WTI_Roll_Normal_Sample5@forecast$VaR[,2],0.05)
WTI_Test_Normal10_Sample5 <- BacktestVaR(WTI_Roll_Normal_Sample5@forecast$VaR[,4],WTI_Roll_Normal_Sample5@forecast$VaR[,3],0.1)

Brent_Test_Normal1_Sample5 <- BacktestVaR(Brent_Roll_Normal_Sample5@forecast$VaR[,4],Brent_Roll_Normal_Sample5@forecast$VaR[,1],0.01)
Brent_Test_Normal5_Sample5 <- BacktestVaR(Brent_Roll_Normal_Sample5@forecast$VaR[,4],Brent_Roll_Normal_Sample5@forecast$VaR[,2],0.05)
Brent_Test_Normal10_Sample5 <- BacktestVaR(Brent_Roll_Normal_Sample5@forecast$VaR[,4],Brent_Roll_Normal_Sample5@forecast$VaR[,3],0.1)


WTI_Sample1 <- matrix(c(WTI_Test_Normal1_Sample1$LRuc[1],WTI_Test_SST1_Sample1$LRuc[1],WTI_Test_SGED1_Sample1$LRuc[1],WTI_Test_JSU1_Sample1$LRuc[1],
  WTI_Test_Normal1_Sample1$LRcc[1],WTI_Test_SST1_Sample1$LRcc[1],WTI_Test_SGED1_Sample1$LRcc[1],WTI_Test_JSU1_Sample1$LRcc[1],
  WTI_Test_Normal1_Sample1$DQ$stat,WTI_Test_SST1_Sample1$DQ$stat,WTI_Test_SGED1_Sample1$DQ$stat,WTI_Test_JSU1_Sample1$DQ$stat,
  WTI_Test_Normal1_Sample1$AE,WTI_Test_SST1_Sample1$AE,WTI_Test_SGED1_Sample1$AE,WTI_Test_JSU1_Sample1$AE,
  (WTI_Test_Normal1_Sample1$Loss$Loss/WTI_Test_Normal1_Sample1$Loss$Loss),(WTI_Test_SST1_Sample1$Loss$Loss/WTI_Test_Normal1_Sample1$Loss$Loss),
  (WTI_Test_SGED1_Sample1$Loss$Loss/WTI_Test_Normal1_Sample1$Loss$Loss),(WTI_Test_JSU1_Sample1$Loss$Loss/WTI_Test_Normal1_Sample1$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(WTI_Sample1) <- c("Normal","SST","SGED","JSU")    
rownames(WTI_Sample1) <- c("UC","CC","DQ","AE","QL")    

WTI_Sample2 <- matrix(c(WTI_Test_Normal1_Sample2$LRuc[1],WTI_Test_SST1_Sample2$LRuc[1],WTI_Test_SGED1_Sample2$LRuc[1],WTI_Test_JSU1_Sample2$LRuc[1],
                          WTI_Test_Normal1_Sample2$LRcc[1],WTI_Test_SST1_Sample2$LRcc[1],WTI_Test_SGED1_Sample2$LRcc[1],WTI_Test_JSU1_Sample2$LRcc[1],
                          WTI_Test_Normal1_Sample2$DQ$stat,WTI_Test_SST1_Sample2$DQ$stat,WTI_Test_SGED1_Sample2$DQ$stat,WTI_Test_JSU1_Sample2$DQ$stat,
                          WTI_Test_Normal1_Sample2$AE,WTI_Test_SST1_Sample2$AE,WTI_Test_SGED1_Sample2$AE,WTI_Test_JSU1_Sample2$AE,
                          (WTI_Test_Normal1_Sample2$Loss$Loss/WTI_Test_Normal1_Sample2$Loss$Loss),(WTI_Test_SST1_Sample2$Loss$Loss/WTI_Test_Normal1_Sample2$Loss$Loss),
                          (WTI_Test_SGED1_Sample2$Loss$Loss/WTI_Test_Normal1_Sample2$Loss$Loss),(WTI_Test_JSU1_Sample2$Loss$Loss/WTI_Test_Normal1_Sample2$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(WTI_Sample2) <- c("Normal","SST","SGED","JSU")    
rownames(WTI_Sample2) <- c("UC","CC","DQ","AE","QL")    

WTI_Sample3 <- matrix(c(WTI_Test_Normal1_Sample3$LRuc[1],WTI_Test_SST1_Sample3$LRuc[1],WTI_Test_SGED1_Sample3$LRuc[1],WTI_Test_JSU1_Sample3$LRuc[1],
                          WTI_Test_Normal1_Sample3$LRcc[1],WTI_Test_SST1_Sample3$LRcc[1],WTI_Test_SGED1_Sample3$LRcc[1],WTI_Test_JSU1_Sample3$LRcc[1],
                          WTI_Test_Normal1_Sample3$DQ$stat,WTI_Test_SST1_Sample3$DQ$stat,WTI_Test_SGED1_Sample3$DQ$stat,WTI_Test_JSU1_Sample3$DQ$stat,
                          WTI_Test_Normal1_Sample3$AE,WTI_Test_SST1_Sample3$AE,WTI_Test_SGED1_Sample3$AE,WTI_Test_JSU1_Sample3$AE,
                          (WTI_Test_Normal1_Sample3$Loss$Loss/WTI_Test_Normal1_Sample3$Loss$Loss),(WTI_Test_SST1_Sample3$Loss$Loss/WTI_Test_Normal1_Sample3$Loss$Loss),
                          (WTI_Test_SGED1_Sample3$Loss$Loss/WTI_Test_Normal1_Sample3$Loss$Loss),(WTI_Test_JSU1_Sample3$Loss$Loss/WTI_Test_Normal1_Sample3$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(WTI_Sample3) <- c("Normal","SST","SGED","JSU")    
rownames(WTI_Sample3) <- c("UC","CC","DQ","AE","QL")    

WTI_Sample4 <- matrix(c(WTI_Test_Normal1_Sample4$LRuc[1],WTI_Test_SST1_Sample4$LRuc[1],WTI_Test_SGED1_Sample4$LRuc[1],WTI_Test_JSU1_Sample4$LRuc[1],
                          WTI_Test_Normal1_Sample4$LRcc[1],WTI_Test_SST1_Sample4$LRcc[1],WTI_Test_SGED1_Sample4$LRcc[1],WTI_Test_JSU1_Sample4$LRcc[1],
                          WTI_Test_Normal1_Sample4$DQ$stat,WTI_Test_SST1_Sample4$DQ$stat,WTI_Test_SGED1_Sample4$DQ$stat,WTI_Test_JSU1_Sample4$DQ$stat,
                          WTI_Test_Normal1_Sample4$AE,WTI_Test_SST1_Sample4$AE,WTI_Test_SGED1_Sample4$AE,WTI_Test_JSU1_Sample4$AE,
                          (WTI_Test_Normal1_Sample4$Loss$Loss/WTI_Test_Normal1_Sample4$Loss$Loss),(WTI_Test_SST1_Sample4$Loss$Loss/WTI_Test_Normal1_Sample4$Loss$Loss),
                          (WTI_Test_SGED1_Sample4$Loss$Loss/WTI_Test_Normal1_Sample4$Loss$Loss),(WTI_Test_JSU1_Sample4$Loss$Loss/WTI_Test_Normal1_Sample4$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(WTI_Sample4) <- c("Normal","SST","SGED","JSU")    
rownames(WTI_Sample4) <- c("UC","CC","DQ","AE","QL")    

WTI_Sample5 <- matrix(c(WTI_Test_Normal1_Sample5$LRuc[1],WTI_Test_SST1_Sample5$LRuc[1],WTI_Test_SGED1_Sample5$LRuc[1],WTI_Test_JSU1_Sample5$LRuc[1],
                          WTI_Test_Normal1_Sample5$LRcc[1],WTI_Test_SST1_Sample5$LRcc[1],WTI_Test_SGED1_Sample5$LRcc[1],WTI_Test_JSU1_Sample5$LRcc[1],
                          WTI_Test_Normal1_Sample5$DQ$stat,WTI_Test_SST1_Sample5$DQ$stat,WTI_Test_SGED1_Sample5$DQ$stat,WTI_Test_JSU1_Sample5$DQ$stat,
                          WTI_Test_Normal1_Sample5$AE,WTI_Test_SST1_Sample5$AE,WTI_Test_SGED1_Sample5$AE,WTI_Test_JSU1_Sample5$AE,
                          (WTI_Test_Normal1_Sample5$Loss$Loss/WTI_Test_Normal1_Sample5$Loss$Loss),(WTI_Test_SST1_Sample5$Loss$Loss/WTI_Test_Normal1_Sample5$Loss$Loss),
                          (WTI_Test_SGED1_Sample5$Loss$Loss/WTI_Test_Normal1_Sample5$Loss$Loss),(WTI_Test_JSU1_Sample5$Loss$Loss/WTI_Test_Normal1_Sample5$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(WTI_Sample5) <- c("Normal","SST","SGED","JSU")    
rownames(WTI_Sample5) <- c("UC","CC","DQ","AE","QL")    


Brent_Sample1 <- matrix(c(Brent_Test_Normal1_Sample1$LRuc[1],Brent_Test_SST1_Sample1$LRuc[1],Brent_Test_SGED1_Sample1$LRuc[1],Brent_Test_JSU1_Sample1$LRuc[1],
                          Brent_Test_Normal1_Sample1$LRcc[1],Brent_Test_SST1_Sample1$LRcc[1],Brent_Test_SGED1_Sample1$LRcc[1],Brent_Test_JSU1_Sample1$LRcc[1],
                          Brent_Test_Normal1_Sample1$DQ$stat,Brent_Test_SST1_Sample1$DQ$stat,Brent_Test_SGED1_Sample1$DQ$stat,Brent_Test_JSU1_Sample1$DQ$stat,
                          Brent_Test_Normal1_Sample1$AE,Brent_Test_SST1_Sample1$AE,Brent_Test_SGED1_Sample1$AE,Brent_Test_JSU1_Sample1$AE,
                          (Brent_Test_Normal1_Sample1$Loss$Loss/Brent_Test_Normal1_Sample1$Loss$Loss),(Brent_Test_SST1_Sample1$Loss$Loss/Brent_Test_Normal1_Sample1$Loss$Loss),
                          (Brent_Test_SGED1_Sample1$Loss$Loss/Brent_Test_Normal1_Sample1$Loss$Loss),(Brent_Test_JSU1_Sample1$Loss$Loss/Brent_Test_Normal1_Sample1$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Brent_Sample1) <- c("Normal","SST","SGED","JSU")    
rownames(Brent_Sample1) <- c("UC","CC","DQ","AE","QL")    

Brent_Sample2 <- matrix(c(Brent_Test_Normal1_Sample2$LRuc[1],Brent_Test_SST1_Sample2$LRuc[1],Brent_Test_SGED1_Sample2$LRuc[1],Brent_Test_JSU1_Sample2$LRuc[1],
                          Brent_Test_Normal1_Sample2$LRcc[1],Brent_Test_SST1_Sample2$LRcc[1],Brent_Test_SGED1_Sample2$LRcc[1],Brent_Test_JSU1_Sample2$LRcc[1],
                          Brent_Test_Normal1_Sample2$DQ$stat,Brent_Test_SST1_Sample2$DQ$stat,Brent_Test_SGED1_Sample2$DQ$stat,Brent_Test_JSU1_Sample2$DQ$stat,
                          Brent_Test_Normal1_Sample2$AE,Brent_Test_SST1_Sample2$AE,Brent_Test_SGED1_Sample2$AE,Brent_Test_JSU1_Sample2$AE,
                          (Brent_Test_Normal1_Sample2$Loss$Loss/Brent_Test_Normal1_Sample2$Loss$Loss),(Brent_Test_SST1_Sample2$Loss$Loss/Brent_Test_Normal1_Sample2$Loss$Loss),
                          (Brent_Test_SGED1_Sample2$Loss$Loss/Brent_Test_Normal1_Sample2$Loss$Loss),(Brent_Test_JSU1_Sample2$Loss$Loss/Brent_Test_Normal1_Sample2$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Brent_Sample2) <- c("Normal","SST","SGED","JSU")    
rownames(Brent_Sample2) <- c("UC","CC","DQ","AE","QL")    

Brent_Sample3 <- matrix(c(Brent_Test_Normal1_Sample3$LRuc[1],Brent_Test_SST1_Sample3$LRuc[1],Brent_Test_SGED1_Sample3$LRuc[1],Brent_Test_JSU1_Sample3$LRuc[1],
                          Brent_Test_Normal1_Sample3$LRcc[1],Brent_Test_SST1_Sample3$LRcc[1],Brent_Test_SGED1_Sample3$LRcc[1],Brent_Test_JSU1_Sample3$LRcc[1],
                          Brent_Test_Normal1_Sample3$DQ$stat,Brent_Test_SST1_Sample3$DQ$stat,Brent_Test_SGED1_Sample3$DQ$stat,Brent_Test_JSU1_Sample3$DQ$stat,
                          Brent_Test_Normal1_Sample3$AE,Brent_Test_SST1_Sample3$AE,Brent_Test_SGED1_Sample3$AE,Brent_Test_JSU1_Sample3$AE,
                          (Brent_Test_Normal1_Sample3$Loss$Loss/Brent_Test_Normal1_Sample3$Loss$Loss),(Brent_Test_SST1_Sample3$Loss$Loss/Brent_Test_Normal1_Sample3$Loss$Loss),
                          (Brent_Test_SGED1_Sample3$Loss$Loss/Brent_Test_Normal1_Sample3$Loss$Loss),(Brent_Test_JSU1_Sample3$Loss$Loss/Brent_Test_Normal1_Sample3$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Brent_Sample3) <- c("Normal","SST","SGED","JSU")    
rownames(Brent_Sample3) <- c("UC","CC","DQ","AE","QL")    

Brent_Sample4 <- matrix(c(Brent_Test_Normal1_Sample4$LRuc[1],Brent_Test_SST1_Sample4$LRuc[1],Brent_Test_SGED1_Sample4$LRuc[1],Brent_Test_JSU1_Sample4$LRuc[1],
                          Brent_Test_Normal1_Sample4$LRcc[1],Brent_Test_SST1_Sample4$LRcc[1],Brent_Test_SGED1_Sample4$LRcc[1],Brent_Test_JSU1_Sample4$LRcc[1],
                          Brent_Test_Normal1_Sample4$DQ$stat,Brent_Test_SST1_Sample4$DQ$stat,Brent_Test_SGED1_Sample4$DQ$stat,Brent_Test_JSU1_Sample4$DQ$stat,
                          Brent_Test_Normal1_Sample4$AE,Brent_Test_SST1_Sample4$AE,Brent_Test_SGED1_Sample4$AE,Brent_Test_JSU1_Sample4$AE,
                          (Brent_Test_Normal1_Sample4$Loss$Loss/Brent_Test_Normal1_Sample4$Loss$Loss),(Brent_Test_SST1_Sample4$Loss$Loss/Brent_Test_Normal1_Sample4$Loss$Loss),
                          (Brent_Test_SGED1_Sample4$Loss$Loss/Brent_Test_Normal1_Sample4$Loss$Loss),(Brent_Test_JSU1_Sample4$Loss$Loss/Brent_Test_Normal1_Sample4$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Brent_Sample4) <- c("Normal","SST","SGED","JSU")    
rownames(Brent_Sample4) <- c("UC","CC","DQ","AE","QL")    

Brent_Sample5 <- matrix(c(Brent_Test_Normal1_Sample5$LRuc[1],Brent_Test_SST1_Sample5$LRuc[1],Brent_Test_SGED1_Sample5$LRuc[1],Brent_Test_JSU1_Sample5$LRuc[1],
                          Brent_Test_Normal1_Sample5$LRcc[1],Brent_Test_SST1_Sample5$LRcc[1],Brent_Test_SGED1_Sample5$LRcc[1],Brent_Test_JSU1_Sample5$LRcc[1],
                          Brent_Test_Normal1_Sample5$DQ$stat,Brent_Test_SST1_Sample5$DQ$stat,Brent_Test_SGED1_Sample5$DQ$stat,Brent_Test_JSU1_Sample5$DQ$stat,
                          Brent_Test_Normal1_Sample5$AE,Brent_Test_SST1_Sample5$AE,Brent_Test_SGED1_Sample5$AE,Brent_Test_JSU1_Sample5$AE,
                          (Brent_Test_Normal1_Sample5$Loss$Loss/Brent_Test_Normal1_Sample5$Loss$Loss),(Brent_Test_SST1_Sample5$Loss$Loss/Brent_Test_Normal1_Sample5$Loss$Loss),
                          (Brent_Test_SGED1_Sample5$Loss$Loss/Brent_Test_Normal1_Sample5$Loss$Loss),(Brent_Test_JSU1_Sample5$Loss$Loss/Brent_Test_Normal1_Sample5$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Brent_Sample5) <- c("Normal","SST","SGED","JSU")    
rownames(Brent_Sample5) <- c("UC","CC","DQ","AE","QL")   

### Storing the p-values for subsamples 

# 1 percent
WTI_pvalues_1pc_Sample1 <- matrix(c(WTI_Test_Normal1_Sample1$LRuc[2],WTI_Test_SST1_Sample1$LRuc[2],WTI_Test_SGED1_Sample1$LRuc[2],WTI_Test_JSU1_Sample1$LRuc[2],
                          WTI_Test_Normal1_Sample1$LRcc[2],WTI_Test_SST1_Sample1$LRcc[2],WTI_Test_SGED1_Sample1$LRcc[2],WTI_Test_JSU1_Sample1$LRcc[2],
                          WTI_Test_Normal1_Sample1$DQ$pvalue,WTI_Test_SST1_Sample1$DQ$pvalue,WTI_Test_SGED1_Sample1$DQ$pvalue,WTI_Test_JSU1_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_1pc_Sample2 <- matrix(c(WTI_Test_Normal1_Sample2$LRuc[2],WTI_Test_SST1_Sample2$LRuc[2],WTI_Test_SGED1_Sample2$LRuc[2],WTI_Test_JSU1_Sample2$LRuc[2],
                                      WTI_Test_Normal1_Sample2$LRcc[2],WTI_Test_SST1_Sample2$LRcc[2],WTI_Test_SGED1_Sample2$LRcc[2],WTI_Test_JSU1_Sample2$LRcc[2],
                                      WTI_Test_Normal1_Sample2$DQ$pvalue,WTI_Test_SST1_Sample2$DQ$pvalue,WTI_Test_SGED1_Sample2$DQ$pvalue,WTI_Test_JSU1_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_1pc_Sample3 <- matrix(c(WTI_Test_Normal1_Sample3$LRuc[2],WTI_Test_SST1_Sample3$LRuc[2],WTI_Test_SGED1_Sample3$LRuc[2],WTI_Test_JSU1_Sample3$LRuc[2],
                                      WTI_Test_Normal1_Sample3$LRcc[2],WTI_Test_SST1_Sample3$LRcc[2],WTI_Test_SGED1_Sample3$LRcc[2],WTI_Test_JSU1_Sample3$LRcc[2],
                                      WTI_Test_Normal1_Sample3$DQ$pvalue,WTI_Test_SST1_Sample3$DQ$pvalue,WTI_Test_SGED1_Sample3$DQ$pvalue,WTI_Test_JSU1_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_1pc_Sample4 <- matrix(c(WTI_Test_Normal1_Sample4$LRuc[2],WTI_Test_SST1_Sample4$LRuc[2],WTI_Test_SGED1_Sample4$LRuc[2],WTI_Test_JSU1_Sample4$LRuc[2],
                                      WTI_Test_Normal1_Sample4$LRcc[2],WTI_Test_SST1_Sample4$LRcc[2],WTI_Test_SGED1_Sample4$LRcc[2],WTI_Test_JSU1_Sample4$LRcc[2],
                                      WTI_Test_Normal1_Sample4$DQ$pvalue,WTI_Test_SST1_Sample4$DQ$pvalue,WTI_Test_SGED1_Sample4$DQ$pvalue,WTI_Test_JSU1_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_1pc_Sample5 <- matrix(c(WTI_Test_Normal1_Sample5$LRuc[2],WTI_Test_SST1_Sample5$LRuc[2],WTI_Test_SGED1_Sample5$LRuc[2],WTI_Test_JSU1_Sample5$LRuc[2],
                                      WTI_Test_Normal1_Sample5$LRcc[2],WTI_Test_SST1_Sample5$LRcc[2],WTI_Test_SGED1_Sample5$LRcc[2],WTI_Test_JSU1_Sample5$LRcc[2],
                                      WTI_Test_Normal1_Sample5$DQ$pvalue,WTI_Test_SST1_Sample5$DQ$pvalue,WTI_Test_SGED1_Sample5$DQ$pvalue,WTI_Test_JSU1_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)


Brent_pvalues_1pc_Sample1 <- matrix(c(Brent_Test_Normal1_Sample1$LRuc[2],Brent_Test_SST1_Sample1$LRuc[2],Brent_Test_SGED1_Sample1$LRuc[2],Brent_Test_JSU1_Sample1$LRuc[2],
                                      Brent_Test_Normal1_Sample1$LRcc[2],Brent_Test_SST1_Sample1$LRcc[2],Brent_Test_SGED1_Sample1$LRcc[2],Brent_Test_JSU1_Sample1$LRcc[2],
                                      Brent_Test_Normal1_Sample1$DQ$pvalue,Brent_Test_SST1_Sample1$DQ$pvalue,Brent_Test_SGED1_Sample1$DQ$pvalue,Brent_Test_JSU1_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_1pc_Sample2 <- matrix(c(Brent_Test_Normal1_Sample2$LRuc[2],Brent_Test_SST1_Sample2$LRuc[2],Brent_Test_SGED1_Sample2$LRuc[2],Brent_Test_JSU1_Sample2$LRuc[2],
                                      Brent_Test_Normal1_Sample2$LRcc[2],Brent_Test_SST1_Sample2$LRcc[2],Brent_Test_SGED1_Sample2$LRcc[2],Brent_Test_JSU1_Sample2$LRcc[2],
                                      Brent_Test_Normal1_Sample2$DQ$pvalue,Brent_Test_SST1_Sample2$DQ$pvalue,Brent_Test_SGED1_Sample2$DQ$pvalue,Brent_Test_JSU1_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_1pc_Sample3 <- matrix(c(Brent_Test_Normal1_Sample3$LRuc[2],Brent_Test_SST1_Sample3$LRuc[2],Brent_Test_SGED1_Sample3$LRuc[2],Brent_Test_JSU1_Sample3$LRuc[2],
                                      Brent_Test_Normal1_Sample3$LRcc[2],Brent_Test_SST1_Sample3$LRcc[2],Brent_Test_SGED1_Sample3$LRcc[2],Brent_Test_JSU1_Sample3$LRcc[2],
                                      Brent_Test_Normal1_Sample3$DQ$pvalue,Brent_Test_SST1_Sample3$DQ$pvalue,Brent_Test_SGED1_Sample3$DQ$pvalue,Brent_Test_JSU1_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_1pc_Sample4 <- matrix(c(Brent_Test_Normal1_Sample4$LRuc[2],Brent_Test_SST1_Sample4$LRuc[2],Brent_Test_SGED1_Sample4$LRuc[2],Brent_Test_JSU1_Sample4$LRuc[2],
                                      Brent_Test_Normal1_Sample4$LRcc[2],Brent_Test_SST1_Sample4$LRcc[2],Brent_Test_SGED1_Sample4$LRcc[2],Brent_Test_JSU1_Sample4$LRcc[2],
                                      Brent_Test_Normal1_Sample4$DQ$pvalue,Brent_Test_SST1_Sample4$DQ$pvalue,Brent_Test_SGED1_Sample4$DQ$pvalue,Brent_Test_JSU1_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_1pc_Sample5 <- matrix(c(Brent_Test_Normal1_Sample5$LRuc[2],Brent_Test_SST1_Sample5$LRuc[2],Brent_Test_SGED1_Sample5$LRuc[2],Brent_Test_JSU1_Sample5$LRuc[2],
                                      Brent_Test_Normal1_Sample5$LRcc[2],Brent_Test_SST1_Sample5$LRcc[2],Brent_Test_SGED1_Sample5$LRcc[2],Brent_Test_JSU1_Sample5$LRcc[2],
                                      Brent_Test_Normal1_Sample5$DQ$pvalue,Brent_Test_SST1_Sample5$DQ$pvalue,Brent_Test_SGED1_Sample5$DQ$pvalue,Brent_Test_JSU1_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)

# 5 percent

WTI_pvalues_5pc_Sample1 <- matrix(c(WTI_Test_Normal5_Sample1$LRuc[2],WTI_Test_SST5_Sample1$LRuc[2],WTI_Test_SGED5_Sample1$LRuc[2],WTI_Test_JSU5_Sample1$LRuc[2],
                                      WTI_Test_Normal5_Sample1$LRcc[2],WTI_Test_SST5_Sample1$LRcc[2],WTI_Test_SGED5_Sample1$LRcc[2],WTI_Test_JSU5_Sample1$LRcc[2],
                                      WTI_Test_Normal5_Sample1$DQ$pvalue,WTI_Test_SST5_Sample1$DQ$pvalue,WTI_Test_SGED5_Sample1$DQ$pvalue,WTI_Test_JSU5_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_5pc_Sample2 <- matrix(c(WTI_Test_Normal5_Sample2$LRuc[2],WTI_Test_SST5_Sample2$LRuc[2],WTI_Test_SGED5_Sample2$LRuc[2],WTI_Test_JSU5_Sample2$LRuc[2],
                                      WTI_Test_Normal5_Sample2$LRcc[2],WTI_Test_SST5_Sample2$LRcc[2],WTI_Test_SGED5_Sample2$LRcc[2],WTI_Test_JSU5_Sample2$LRcc[2],
                                      WTI_Test_Normal5_Sample2$DQ$pvalue,WTI_Test_SST5_Sample2$DQ$pvalue,WTI_Test_SGED5_Sample2$DQ$pvalue,WTI_Test_JSU5_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_5pc_Sample3 <- matrix(c(WTI_Test_Normal5_Sample3$LRuc[2],WTI_Test_SST5_Sample3$LRuc[2],WTI_Test_SGED5_Sample3$LRuc[2],WTI_Test_JSU5_Sample3$LRuc[2],
                                      WTI_Test_Normal5_Sample3$LRcc[2],WTI_Test_SST5_Sample3$LRcc[2],WTI_Test_SGED5_Sample3$LRcc[2],WTI_Test_JSU5_Sample3$LRcc[2],
                                      WTI_Test_Normal5_Sample3$DQ$pvalue,WTI_Test_SST5_Sample3$DQ$pvalue,WTI_Test_SGED5_Sample3$DQ$pvalue,WTI_Test_JSU5_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_5pc_Sample4 <- matrix(c(WTI_Test_Normal5_Sample4$LRuc[2],WTI_Test_SST5_Sample4$LRuc[2],WTI_Test_SGED5_Sample4$LRuc[2],WTI_Test_JSU5_Sample4$LRuc[2],
                                      WTI_Test_Normal5_Sample4$LRcc[2],WTI_Test_SST5_Sample4$LRcc[2],WTI_Test_SGED5_Sample4$LRcc[2],WTI_Test_JSU5_Sample4$LRcc[2],
                                      WTI_Test_Normal5_Sample4$DQ$pvalue,WTI_Test_SST5_Sample4$DQ$pvalue,WTI_Test_SGED5_Sample4$DQ$pvalue,WTI_Test_JSU5_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_5pc_Sample5 <- matrix(c(WTI_Test_Normal5_Sample5$LRuc[2],WTI_Test_SST5_Sample5$LRuc[2],WTI_Test_SGED5_Sample5$LRuc[2],WTI_Test_JSU5_Sample5$LRuc[2],
                                      WTI_Test_Normal5_Sample5$LRcc[2],WTI_Test_SST5_Sample5$LRcc[2],WTI_Test_SGED5_Sample5$LRcc[2],WTI_Test_JSU5_Sample5$LRcc[2],
                                      WTI_Test_Normal5_Sample5$DQ$pvalue,WTI_Test_SST5_Sample5$DQ$pvalue,WTI_Test_SGED5_Sample5$DQ$pvalue,WTI_Test_JSU5_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)


Brent_pvalues_5pc_Sample1 <- matrix(c(Brent_Test_Normal5_Sample1$LRuc[2],Brent_Test_SST5_Sample1$LRuc[2],Brent_Test_SGED5_Sample1$LRuc[2],Brent_Test_JSU5_Sample1$LRuc[2],
                                    Brent_Test_Normal5_Sample1$LRcc[2],Brent_Test_SST5_Sample1$LRcc[2],Brent_Test_SGED5_Sample1$LRcc[2],Brent_Test_JSU5_Sample1$LRcc[2],
                                    Brent_Test_Normal5_Sample1$DQ$pvalue,Brent_Test_SST5_Sample1$DQ$pvalue,Brent_Test_SGED5_Sample1$DQ$pvalue,Brent_Test_JSU5_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_5pc_Sample2 <- matrix(c(Brent_Test_Normal5_Sample2$LRuc[2],Brent_Test_SST5_Sample2$LRuc[2],Brent_Test_SGED5_Sample2$LRuc[2],Brent_Test_JSU5_Sample2$LRuc[2],
                                    Brent_Test_Normal5_Sample2$LRcc[2],Brent_Test_SST5_Sample2$LRcc[2],Brent_Test_SGED5_Sample2$LRcc[2],Brent_Test_JSU5_Sample2$LRcc[2],
                                    Brent_Test_Normal5_Sample2$DQ$pvalue,Brent_Test_SST5_Sample2$DQ$pvalue,Brent_Test_SGED5_Sample2$DQ$pvalue,Brent_Test_JSU5_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_5pc_Sample3 <- matrix(c(Brent_Test_Normal5_Sample3$LRuc[2],Brent_Test_SST5_Sample3$LRuc[2],Brent_Test_SGED5_Sample3$LRuc[2],Brent_Test_JSU5_Sample3$LRuc[2],
                                    Brent_Test_Normal5_Sample3$LRcc[2],Brent_Test_SST5_Sample3$LRcc[2],Brent_Test_SGED5_Sample3$LRcc[2],Brent_Test_JSU5_Sample3$LRcc[2],
                                    Brent_Test_Normal5_Sample3$DQ$pvalue,Brent_Test_SST5_Sample3$DQ$pvalue,Brent_Test_SGED5_Sample3$DQ$pvalue,Brent_Test_JSU5_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_5pc_Sample4 <- matrix(c(Brent_Test_Normal5_Sample4$LRuc[2],Brent_Test_SST5_Sample4$LRuc[2],Brent_Test_SGED5_Sample4$LRuc[2],Brent_Test_JSU5_Sample4$LRuc[2],
                                    Brent_Test_Normal5_Sample4$LRcc[2],Brent_Test_SST5_Sample4$LRcc[2],Brent_Test_SGED5_Sample4$LRcc[2],Brent_Test_JSU5_Sample4$LRcc[2],
                                    Brent_Test_Normal5_Sample4$DQ$pvalue,Brent_Test_SST5_Sample4$DQ$pvalue,Brent_Test_SGED5_Sample4$DQ$pvalue,Brent_Test_JSU5_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_5pc_Sample5 <- matrix(c(Brent_Test_Normal5_Sample5$LRuc[2],Brent_Test_SST5_Sample5$LRuc[2],Brent_Test_SGED5_Sample5$LRuc[2],Brent_Test_JSU5_Sample5$LRuc[2],
                                    Brent_Test_Normal5_Sample5$LRcc[2],Brent_Test_SST5_Sample5$LRcc[2],Brent_Test_SGED5_Sample5$LRcc[2],Brent_Test_JSU5_Sample5$LRcc[2],
                                    Brent_Test_Normal5_Sample5$DQ$pvalue,Brent_Test_SST5_Sample5$DQ$pvalue,Brent_Test_SGED5_Sample5$DQ$pvalue,Brent_Test_JSU5_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)

# 10 percent

WTI_pvalues_10pc_Sample1 <- matrix(c(WTI_Test_Normal10_Sample1$LRuc[2],WTI_Test_SST10_Sample1$LRuc[2],WTI_Test_SGED10_Sample1$LRuc[2],WTI_Test_JSU10_Sample1$LRuc[2],
                                      WTI_Test_Normal10_Sample1$LRcc[2],WTI_Test_SST10_Sample1$LRcc[2],WTI_Test_SGED10_Sample1$LRcc[2],WTI_Test_JSU10_Sample1$LRcc[2],
                                      WTI_Test_Normal10_Sample1$DQ$pvalue,WTI_Test_SST10_Sample1$DQ$pvalue,WTI_Test_SGED10_Sample1$DQ$pvalue,WTI_Test_JSU10_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_10pc_Sample2 <- matrix(c(WTI_Test_Normal10_Sample2$LRuc[2],WTI_Test_SST10_Sample2$LRuc[2],WTI_Test_SGED10_Sample2$LRuc[2],WTI_Test_JSU10_Sample2$LRuc[2],
                                      WTI_Test_Normal10_Sample2$LRcc[2],WTI_Test_SST10_Sample2$LRcc[2],WTI_Test_SGED10_Sample2$LRcc[2],WTI_Test_JSU10_Sample2$LRcc[2],
                                      WTI_Test_Normal10_Sample2$DQ$pvalue,WTI_Test_SST10_Sample2$DQ$pvalue,WTI_Test_SGED10_Sample2$DQ$pvalue,WTI_Test_JSU10_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_10pc_Sample3 <- matrix(c(WTI_Test_Normal10_Sample3$LRuc[2],WTI_Test_SST10_Sample3$LRuc[2],WTI_Test_SGED10_Sample3$LRuc[2],WTI_Test_JSU10_Sample3$LRuc[2],
                                      WTI_Test_Normal10_Sample3$LRcc[2],WTI_Test_SST10_Sample3$LRcc[2],WTI_Test_SGED10_Sample3$LRcc[2],WTI_Test_JSU10_Sample3$LRcc[2],
                                      WTI_Test_Normal10_Sample3$DQ$pvalue,WTI_Test_SST10_Sample3$DQ$pvalue,WTI_Test_SGED10_Sample3$DQ$pvalue,WTI_Test_JSU10_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_10pc_Sample4 <- matrix(c(WTI_Test_Normal10_Sample4$LRuc[2],WTI_Test_SST10_Sample4$LRuc[2],WTI_Test_SGED10_Sample4$LRuc[2],WTI_Test_JSU10_Sample4$LRuc[2],
                                      WTI_Test_Normal10_Sample4$LRcc[2],WTI_Test_SST10_Sample4$LRcc[2],WTI_Test_SGED10_Sample4$LRcc[2],WTI_Test_JSU10_Sample4$LRcc[2],
                                      WTI_Test_Normal10_Sample4$DQ$pvalue,WTI_Test_SST10_Sample4$DQ$pvalue,WTI_Test_SGED10_Sample4$DQ$pvalue,WTI_Test_JSU10_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

WTI_pvalues_10pc_Sample5 <- matrix(c(WTI_Test_Normal10_Sample5$LRuc[2],WTI_Test_SST10_Sample5$LRuc[2],WTI_Test_SGED10_Sample5$LRuc[2],WTI_Test_JSU10_Sample5$LRuc[2],
                                      WTI_Test_Normal10_Sample5$LRcc[2],WTI_Test_SST10_Sample5$LRcc[2],WTI_Test_SGED10_Sample5$LRcc[2],WTI_Test_JSU10_Sample5$LRcc[2],
                                      WTI_Test_Normal10_Sample5$DQ$pvalue,WTI_Test_SST10_Sample5$DQ$pvalue,WTI_Test_SGED10_Sample5$DQ$pvalue,WTI_Test_JSU10_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)


Brent_pvalues_10pc_Sample1 <- matrix(c(Brent_Test_Normal10_Sample1$LRuc[2],Brent_Test_SST10_Sample1$LRuc[2],Brent_Test_SGED10_Sample1$LRuc[2],Brent_Test_JSU10_Sample1$LRuc[2],
                                    Brent_Test_Normal10_Sample1$LRcc[2],Brent_Test_SST10_Sample1$LRcc[2],Brent_Test_SGED10_Sample1$LRcc[2],Brent_Test_JSU10_Sample1$LRcc[2],
                                    Brent_Test_Normal10_Sample1$DQ$pvalue,Brent_Test_SST10_Sample1$DQ$pvalue,Brent_Test_SGED10_Sample1$DQ$pvalue,Brent_Test_JSU10_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_10pc_Sample2 <- matrix(c(Brent_Test_Normal10_Sample2$LRuc[2],Brent_Test_SST10_Sample2$LRuc[2],Brent_Test_SGED10_Sample2$LRuc[2],Brent_Test_JSU10_Sample2$LRuc[2],
                                    Brent_Test_Normal10_Sample2$LRcc[2],Brent_Test_SST10_Sample2$LRcc[2],Brent_Test_SGED10_Sample2$LRcc[2],Brent_Test_JSU10_Sample2$LRcc[2],
                                    Brent_Test_Normal10_Sample2$DQ$pvalue,Brent_Test_SST10_Sample2$DQ$pvalue,Brent_Test_SGED10_Sample2$DQ$pvalue,Brent_Test_JSU10_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_10pc_Sample3 <- matrix(c(Brent_Test_Normal10_Sample3$LRuc[2],Brent_Test_SST10_Sample3$LRuc[2],Brent_Test_SGED10_Sample3$LRuc[2],Brent_Test_JSU10_Sample3$LRuc[2],
                                    Brent_Test_Normal10_Sample3$LRcc[2],Brent_Test_SST10_Sample3$LRcc[2],Brent_Test_SGED10_Sample3$LRcc[2],Brent_Test_JSU10_Sample3$LRcc[2],
                                    Brent_Test_Normal10_Sample3$DQ$pvalue,Brent_Test_SST10_Sample3$DQ$pvalue,Brent_Test_SGED10_Sample3$DQ$pvalue,Brent_Test_JSU10_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_10pc_Sample4 <- matrix(c(Brent_Test_Normal10_Sample4$LRuc[2],Brent_Test_SST10_Sample4$LRuc[2],Brent_Test_SGED10_Sample4$LRuc[2],Brent_Test_JSU10_Sample4$LRuc[2],
                                    Brent_Test_Normal10_Sample4$LRcc[2],Brent_Test_SST10_Sample4$LRcc[2],Brent_Test_SGED10_Sample4$LRcc[2],Brent_Test_JSU10_Sample4$LRcc[2],
                                    Brent_Test_Normal10_Sample4$DQ$pvalue,Brent_Test_SST10_Sample4$DQ$pvalue,Brent_Test_SGED10_Sample4$DQ$pvalue,Brent_Test_JSU10_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Brent_pvalues_10pc_Sample5 <- matrix(c(Brent_Test_Normal10_Sample5$LRuc[2],Brent_Test_SST10_Sample5$LRuc[2],Brent_Test_SGED10_Sample5$LRuc[2],Brent_Test_JSU10_Sample5$LRuc[2],
                                    Brent_Test_Normal10_Sample5$LRcc[2],Brent_Test_SST10_Sample5$LRcc[2],Brent_Test_SGED10_Sample5$LRcc[2],Brent_Test_JSU10_Sample5$LRcc[2],
                                    Brent_Test_Normal10_Sample5$DQ$pvalue,Brent_Test_SST10_Sample5$DQ$pvalue,Brent_Test_SGED10_Sample5$DQ$pvalue,Brent_Test_JSU10_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)


WTI_pvalues_1pc   <- rbind(WTI_pvalues_1pc_Sample1,WTI_pvalues_1pc_Sample2,WTI_pvalues_1pc_Sample3,WTI_pvalues_1pc_Sample4,WTI_pvalues_1pc_Sample5)
WTI_pvalues_5pc   <- rbind(WTI_pvalues_5pc_Sample1,WTI_pvalues_5pc_Sample2,WTI_pvalues_5pc_Sample3,WTI_pvalues_5pc_Sample4,WTI_pvalues_5pc_Sample5)
WTI_pvalues_10pc  <- rbind(WTI_pvalues_10pc_Sample1,WTI_pvalues_10pc_Sample2,WTI_pvalues_10pc_Sample3,WTI_pvalues_10pc_Sample4,WTI_pvalues_10pc_Sample5)

Brent_pvalues_1pc   <- rbind(Brent_pvalues_1pc_Sample1,Brent_pvalues_1pc_Sample2,Brent_pvalues_1pc_Sample3,Brent_pvalues_1pc_Sample4,Brent_pvalues_1pc_Sample5)
Brent_pvalues_5pc   <- rbind(Brent_pvalues_5pc_Sample1,Brent_pvalues_5pc_Sample2,Brent_pvalues_5pc_Sample3,Brent_pvalues_5pc_Sample4,Brent_pvalues_5pc_Sample5)
Brent_pvalues_10pc  <- rbind(Brent_pvalues_10pc_Sample1,Brent_pvalues_10pc_Sample2,Brent_pvalues_10pc_Sample3,Brent_pvalues_10pc_Sample4,Brent_pvalues_10pc_Sample5)

rownames(WTI_pvalues_1pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(WTI_pvalues_1pc) <- c("Normal","SST","SGED","JSU")

rownames(WTI_pvalues_5pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(WTI_pvalues_5pc) <- c("Normal","SST","SGED","JSU")

rownames(WTI_pvalues_10pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(WTI_pvalues_10pc) <- c("Normal","SST","SGED","JSU")

rownames(Brent_pvalues_1pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Brent_pvalues_1pc) <- c("Normal","SST","SGED","JSU")

rownames(Brent_pvalues_5pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Brent_pvalues_5pc) <- c("Normal","SST","SGED","JSU")

rownames(Brent_pvalues_10pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Brent_pvalues_10pc) <- c("Normal","SST","SGED","JSU")

#### Repeating for Pearson Distribution  #####

#signific_levels <- c(0.1,0.05,0.01)
## Store the volatility and mu forecast for Sample 1

Holdout_Sample1_WTI_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample1_Brent_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample1_WTI_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample1_Brent_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample1_WTI_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample1_Brent_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample1_WTI_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample1_Brent_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample1_WTI_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample1_Brent_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample1_WTI_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample1_Brent_Pearson_pvalues) <- c("1%","5%","10%")

Sample1_WTI_Holdout_Mu <- as.data.frame(WTI_Roll_Normal_Sample1)[,'Mu']
Sample1_WTI_Holdout_Sigma <- as.data.frame(WTI_Roll_Normal_Sample1)[,'Sigma']
Sample1_Brent_Holdout_Mu <- as.data.frame(Brent_Roll_Normal_Sample1)[,'Mu']
Sample1_Brent_Holdout_Sigma <- as.data.frame(Brent_Roll_Normal_Sample1)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample1_WTI <- qpearson(signific_levels[j],WTI_Param_Pearson_Sample1)
  Z_Pearson_Sample1_Brent   <- qpearson(signific_levels[j],Brent_Param_Pearson_Sample1)
  
  Sample1_WTI_Holdout_VaR <- Sample1_WTI_Holdout_Mu + Sample1_WTI_Holdout_Sigma*Z_Pearson_Sample1_WTI
  Sample1_Brent_Holdout_VaR   <- Sample1_Brent_Holdout_Mu + Sample1_Brent_Holdout_Sigma*Z_Pearson_Sample1_Brent
  Sample1_WTI_Violations <- BacktestVaR(Returns_WTI[(WTI_samples[1]-250+1):(WTI_samples[1]-1)],Sample1_WTI_Holdout_VaR[-1],signific_levels[j])
  Sample1_Brent_Violations   <- BacktestVaR(Returns_Brent[(Brent_samples[1]-250+1):(Brent_samples[1]-1)],Sample1_Brent_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample1_WTI_Pearson[j,] <- c(Sample1_WTI_Violations$LRuc[1],Sample1_WTI_Violations$LRcc[1],Sample1_WTI_Violations$DQ$stat,Sample1_WTI_Violations$AE,Sample1_WTI_Violations$Loss$Loss)
  Holdout_Sample1_Brent_Pearson[j,]   <- c(Sample1_Brent_Violations$LRuc[1],Sample1_Brent_Violations$LRcc[1],Sample1_Brent_Violations$DQ$stat,Sample1_Brent_Violations$AE,Sample1_Brent_Violations$Loss$Loss)
  Holdout_Sample1_WTI_Pearson_pvalues[j,] <- c(Sample1_WTI_Violations$LRuc[2],Sample1_WTI_Violations$LRcc[2],Sample1_WTI_Violations$DQ$pvalue)
  Holdout_Sample1_Brent_Pearson_pvalues[j,] <- c(Sample1_Brent_Violations$LRuc[2],Sample1_Brent_Violations$LRcc[2],Sample1_Brent_Violations$DQ$pvalue)
}


## Store the volatility and mu forecast for Sample 2

Holdout_Sample2_WTI_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample2_Brent_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample2_WTI_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample2_Brent_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample2_WTI_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample2_Brent_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample2_WTI_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample2_Brent_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample2_WTI_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample2_Brent_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample2_WTI_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample2_Brent_Pearson_pvalues) <- c("1%","5%","10%")

Sample2_WTI_Holdout_Mu <- as.data.frame(WTI_Roll_Normal_Sample2)[,'Mu']
Sample2_WTI_Holdout_Sigma <- as.data.frame(WTI_Roll_Normal_Sample2)[,'Sigma']
Sample2_Brent_Holdout_Mu <- as.data.frame(Brent_Roll_Normal_Sample2)[,'Mu']
Sample2_Brent_Holdout_Sigma <- as.data.frame(Brent_Roll_Normal_Sample2)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample2_WTI <- qpearson(signific_levels[j],WTI_Param_Pearson_Sample2)
  Z_Pearson_Sample2_Brent   <- qpearson(signific_levels[j],Brent_Param_Pearson_Sample2)
  
  Sample2_WTI_Holdout_VaR <- Sample2_WTI_Holdout_Mu + Sample2_WTI_Holdout_Sigma*Z_Pearson_Sample2_WTI
  Sample2_Brent_Holdout_VaR   <- Sample2_Brent_Holdout_Mu + Sample2_Brent_Holdout_Sigma*Z_Pearson_Sample2_Brent
  Sample2_WTI_Violations <- BacktestVaR(Returns_WTI[(WTI_samples[2]-250+1):(WTI_samples[2]-1)],Sample2_WTI_Holdout_VaR[-1],signific_levels[j])
  Sample2_Brent_Violations   <- BacktestVaR(Returns_Brent[(Brent_samples[2]-250+1):(Brent_samples[2]-1)],Sample2_Brent_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample2_WTI_Pearson[j,] <- c(Sample2_WTI_Violations$LRuc[1],Sample2_WTI_Violations$LRcc[1],Sample2_WTI_Violations$DQ$stat,Sample2_WTI_Violations$AE,Sample2_WTI_Violations$Loss$Loss)
  Holdout_Sample2_Brent_Pearson[j,]   <- c(Sample2_Brent_Violations$LRuc[1],Sample2_Brent_Violations$LRcc[1],Sample2_Brent_Violations$DQ$stat,Sample2_Brent_Violations$AE,Sample2_Brent_Violations$Loss$Loss)
  Holdout_Sample2_WTI_Pearson_pvalues[j,] <- c(Sample2_WTI_Violations$LRuc[2],Sample2_WTI_Violations$LRcc[2],Sample2_WTI_Violations$DQ$pvalue)
  Holdout_Sample2_Brent_Pearson_pvalues[j,] <- c(Sample2_Brent_Violations$LRuc[2],Sample2_Brent_Violations$LRcc[2],Sample2_Brent_Violations$DQ$pvalue)
}

## Store the volatility and mu forecast for Sample 3

Holdout_Sample3_WTI_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample3_Brent_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample3_WTI_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample3_Brent_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample3_WTI_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample3_Brent_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample3_WTI_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample3_Brent_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample3_WTI_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample3_Brent_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample3_WTI_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample3_Brent_Pearson_pvalues) <- c("1%","5%","10%")

Sample3_WTI_Holdout_Mu <- as.data.frame(WTI_Roll_Normal_Sample3)[,'Mu']
Sample3_WTI_Holdout_Sigma <- as.data.frame(WTI_Roll_Normal_Sample3)[,'Sigma']
Sample3_Brent_Holdout_Mu <- as.data.frame(Brent_Roll_Normal_Sample3)[,'Mu']
Sample3_Brent_Holdout_Sigma <- as.data.frame(Brent_Roll_Normal_Sample3)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample3_WTI <- qpearson(signific_levels[j],WTI_Param_Pearson_Sample3)
  Z_Pearson_Sample3_Brent   <- qpearson(signific_levels[j],Brent_Param_Pearson_Sample3)
  
  Sample3_WTI_Holdout_VaR <- Sample3_WTI_Holdout_Mu + Sample3_WTI_Holdout_Sigma*Z_Pearson_Sample3_WTI
  Sample3_Brent_Holdout_VaR   <- Sample3_Brent_Holdout_Mu + Sample3_Brent_Holdout_Sigma*Z_Pearson_Sample3_Brent
  Sample3_WTI_Violations <- BacktestVaR(Returns_WTI[(WTI_samples[3]-250+1):(WTI_samples[3]-1)],Sample3_WTI_Holdout_VaR[-1],signific_levels[j])
  Sample3_Brent_Violations   <- BacktestVaR(Returns_Brent[(Brent_samples[3]-250+1):(Brent_samples[3]-1)],Sample3_Brent_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample3_WTI_Pearson[j,] <- c(Sample3_WTI_Violations$LRuc[1],Sample3_WTI_Violations$LRcc[1],Sample3_WTI_Violations$DQ$stat,Sample3_WTI_Violations$AE,Sample3_WTI_Violations$Loss$Loss)
  Holdout_Sample3_Brent_Pearson[j,]   <- c(Sample3_Brent_Violations$LRuc[1],Sample3_Brent_Violations$LRcc[1],Sample3_Brent_Violations$DQ$stat,Sample3_Brent_Violations$AE,Sample3_Brent_Violations$Loss$Loss)
  Holdout_Sample3_WTI_Pearson_pvalues[j,] <- c(Sample3_WTI_Violations$LRuc[2],Sample3_WTI_Violations$LRcc[2],Sample3_WTI_Violations$DQ$pvalue)
  Holdout_Sample3_Brent_Pearson_pvalues[j,] <- c(Sample3_Brent_Violations$LRuc[2],Sample3_Brent_Violations$LRcc[2],Sample3_Brent_Violations$DQ$pvalue)
}

## Store the volatility and mu forecast for Sample 4

Holdout_Sample4_WTI_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample4_Brent_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample4_WTI_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample4_Brent_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample4_WTI_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample4_Brent_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample4_WTI_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample4_Brent_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample4_WTI_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample4_Brent_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample4_WTI_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample4_Brent_Pearson_pvalues) <- c("1%","5%","10%")

Sample4_WTI_Holdout_Mu <- as.data.frame(WTI_Roll_Normal_Sample4)[,'Mu']
Sample4_WTI_Holdout_Sigma <- as.data.frame(WTI_Roll_Normal_Sample4)[,'Sigma']
Sample4_Brent_Holdout_Mu <- as.data.frame(Brent_Roll_Normal_Sample4)[,'Mu']
Sample4_Brent_Holdout_Sigma <- as.data.frame(Brent_Roll_Normal_Sample4)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample4_WTI <- qpearson(signific_levels[j],WTI_Param_Pearson_Sample4)
  Z_Pearson_Sample4_Brent   <- qpearson(signific_levels[j],Brent_Param_Pearson_Sample4)
  
  Sample4_WTI_Holdout_VaR <- Sample4_WTI_Holdout_Mu + Sample4_WTI_Holdout_Sigma*Z_Pearson_Sample4_WTI
  Sample4_Brent_Holdout_VaR   <- Sample4_Brent_Holdout_Mu + Sample4_Brent_Holdout_Sigma*Z_Pearson_Sample4_Brent
  Sample4_WTI_Violations <- BacktestVaR(Returns_WTI[(WTI_samples[4]-250+1):(WTI_samples[4]-1)],Sample4_WTI_Holdout_VaR[-1],signific_levels[j])
  Sample4_Brent_Violations   <- BacktestVaR(Returns_Brent[(Brent_samples[4]-250+1):(Brent_samples[4]-1)],Sample4_Brent_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample4_WTI_Pearson[j,] <- c(Sample4_WTI_Violations$LRuc[1],Sample4_WTI_Violations$LRcc[1],Sample4_WTI_Violations$DQ$stat,Sample4_WTI_Violations$AE,Sample4_WTI_Violations$Loss$Loss)
  Holdout_Sample4_Brent_Pearson[j,]   <- c(Sample4_Brent_Violations$LRuc[1],Sample4_Brent_Violations$LRcc[1],Sample4_Brent_Violations$DQ$stat,Sample4_Brent_Violations$AE,Sample4_Brent_Violations$Loss$Loss)
  Holdout_Sample4_WTI_Pearson_pvalues[j,] <- c(Sample4_WTI_Violations$LRuc[2],Sample4_WTI_Violations$LRcc[2],Sample4_WTI_Violations$DQ$pvalue)
  Holdout_Sample4_Brent_Pearson_pvalues[j,] <- c(Sample4_Brent_Violations$LRuc[2],Sample4_Brent_Violations$LRcc[2],Sample4_Brent_Violations$DQ$pvalue)
}

## Store the volatility and mu forecast for Sample 5

Holdout_Sample5_WTI_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample5_Brent_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample5_WTI_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample5_Brent_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample5_WTI_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample5_Brent_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample5_WTI_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample5_Brent_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample5_WTI_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample5_Brent_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample5_WTI_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample5_Brent_Pearson_pvalues) <- c("1%","5%","10%")

Sample5_WTI_Holdout_Mu <- as.data.frame(WTI_Roll_Normal_Sample5)[,'Mu']
Sample5_WTI_Holdout_Sigma <- as.data.frame(WTI_Roll_Normal_Sample5)[,'Sigma']
Sample5_Brent_Holdout_Mu <- as.data.frame(Brent_Roll_Normal_Sample5)[,'Mu']
Sample5_Brent_Holdout_Sigma <- as.data.frame(Brent_Roll_Normal_Sample5)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample5_WTI <- qpearson(signific_levels[j],WTI_Param_Pearson_Sample5)
  Z_Pearson_Sample5_Brent   <- qpearson(signific_levels[j],Brent_Param_Pearson_Sample5)
  
  Sample5_WTI_Holdout_VaR <- Sample5_WTI_Holdout_Mu + Sample5_WTI_Holdout_Sigma*Z_Pearson_Sample5_WTI
  Sample5_Brent_Holdout_VaR   <- Sample5_Brent_Holdout_Mu + Sample5_Brent_Holdout_Sigma*Z_Pearson_Sample5_Brent
  Sample5_WTI_Violations <- BacktestVaR(Returns_WTI[(length(Returns_WTI)-250+1):(length(Returns_WTI)-1)],Sample5_WTI_Holdout_VaR[-1],signific_levels[j])
  Sample5_Brent_Violations   <- BacktestVaR(Returns_Brent[(length(Returns_Brent)-250+1):(length(Returns_Brent)-1)],Sample5_Brent_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample5_WTI_Pearson[j,] <- c(Sample5_WTI_Violations$LRuc[1],Sample5_WTI_Violations$LRcc[1],Sample5_WTI_Violations$DQ$stat,Sample5_WTI_Violations$AE,Sample5_WTI_Violations$Loss$Loss)
  Holdout_Sample5_Brent_Pearson[j,]   <- c(Sample5_Brent_Violations$LRuc[1],Sample5_Brent_Violations$LRcc[1],Sample5_Brent_Violations$DQ$stat,Sample5_Brent_Violations$AE,Sample5_Brent_Violations$Loss$Loss)
  Holdout_Sample5_WTI_Pearson_pvalues[j,] <- c(Sample5_WTI_Violations$LRuc[2],Sample5_WTI_Violations$LRcc[2],Sample5_WTI_Violations$DQ$pvalue)
  Holdout_Sample5_Brent_Pearson_pvalues[j,] <- c(Sample5_Brent_Violations$LRuc[2],Sample5_Brent_Violations$LRcc[2],Sample5_Brent_Violations$DQ$pvalue)
}

### Plotting the VaR for subsamples

WTI_Sample1_Dates <- format(WTI$Date[(WTI_samples[1]-249):WTI_samples[1]],"%b-%y")
WTI_Sample2_Dates <- format(WTI$Date[(WTI_samples[2]-249):WTI_samples[2]],"%b-%y")
WTI_Sample3_Dates <- format(WTI$Date[(WTI_samples[3]-249):WTI_samples[3]],"%b-%y")
WTI_Sample4_Dates <- format(WTI$Date[(WTI_samples[4]-249):WTI_samples[4]],"%b-%y")
WTI_Sample5_Dates <- format(WTI$Date[(length(WTI_xts)-249):length(WTI_xts)],"%b-%y")

Brent_Sample1_Dates <- format(Brent$Date[(Brent_samples[1]-249):Brent_samples[1]],"%b-%y")
Brent_Sample2_Dates <- format(Brent$Date[(Brent_samples[2]-249):Brent_samples[2]],"%b-%y")
Brent_Sample3_Dates <- format(Brent$Date[(Brent_samples[3]-249):Brent_samples[3]],"%b-%y")
Brent_Sample4_Dates <- format(Brent$Date[(Brent_samples[4]-249):Brent_samples[4]],"%b-%y")
Brent_Sample5_Dates <- format(Brent$Date[(length(Brent_xts)-249):length(Brent_xts)],"%b-%y")

par(mfrow=c(2,3))
plot(Brent_Roll_Normal_Sample1@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample1",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Brent_Roll_Normal_Sample1@forecast$VaR[,3],type="l",col="yellow")
lines(Brent_Roll_SST_Sample1@forecast$VaR[,3],type="l",col="green")
lines(Brent_Roll_SGED_Sample1@forecast$VaR[,3],type="l",col="gray")
lines(Brent_Roll_JSU_Sample1@forecast$VaR[,3],type="l",col="brown")
lines(Sample1_Brent_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Brent_Sample1_Dates[c(1,50,100,150,200,250)],las=2)


plot(Brent_Roll_JSU_Sample2@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample2",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Brent_Roll_Normal_Sample2@forecast$VaR[,3],type="l",col="yellow")
lines(Brent_Roll_SST_Sample2@forecast$VaR[,3],type="l",col="green")
lines(Brent_Roll_SGED_Sample2@forecast$VaR[,3],type="l",col="gray")
lines(Brent_Roll_JSU_Sample2@forecast$VaR[,3],type="l",col="brown")
lines(Sample2_Brent_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Brent_Sample2_Dates[c(1,50,100,150,200,250)],las=2)

plot(Brent_Roll_JSU_Sample3@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample3",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Brent_Roll_Normal_Sample3@forecast$VaR[,3],type="l",col="yellow")
lines(Brent_Roll_SST_Sample3@forecast$VaR[,3],type="l",col="green")
lines(Brent_Roll_SGED_Sample3@forecast$VaR[,3],type="l",col="gray")
lines(Brent_Roll_JSU_Sample3@forecast$VaR[,3],type="l",col="brown")
lines(Sample3_Brent_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Brent_Sample3_Dates[c(1,50,100,150,200,250)],las=2)

plot(Brent_Roll_JSU_Sample4@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample4",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Brent_Roll_Normal_Sample4@forecast$VaR[,3],type="l",col="yellow")
lines(Brent_Roll_SST_Sample4@forecast$VaR[,3],type="l",col="green")
lines(Brent_Roll_SGED_Sample4@forecast$VaR[,3],type="l",col="gray")
lines(Brent_Roll_JSU_Sample4@forecast$VaR[,3],type="l",col="brown")
lines(Sample4_Brent_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Brent_Sample4_Dates[c(1,50,100,150,200,250)],las=2)

plot(Brent_Roll_JSU_Sample5@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample5",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Brent_Roll_Normal_Sample5@forecast$VaR[,3],type="l",col="yellow")
lines(Brent_Roll_SST_Sample5@forecast$VaR[,3],type="l",col="green")
lines(Brent_Roll_SGED_Sample5@forecast$VaR[,3],type="l",col="gray")
lines(Brent_Roll_JSU_Sample5@forecast$VaR[,3],type="l",col="brown")
lines(Sample5_Brent_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Brent_Sample5_Dates[c(1,50,100,150,200,250)],las=2)

### For WTI

par(mfrow=c(2,3))
plot(WTI_Roll_Normal_Sample1@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample1",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(WTI_Roll_Normal_Sample1@forecast$VaR[,3],type="l",col="yellow")
lines(WTI_Roll_SST_Sample1@forecast$VaR[,3],type="l",col="green")
lines(WTI_Roll_SGED_Sample1@forecast$VaR[,3],type="l",col="gray")
lines(WTI_Roll_JSU_Sample1@forecast$VaR[,3],type="l",col="brown")
lines(Sample1_WTI_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=WTI_Sample1_Dates[c(1,50,100,150,200,250)],las=2)


plot(WTI_Roll_JSU_Sample2@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample2",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(WTI_Roll_Normal_Sample2@forecast$VaR[,3],type="l",col="yellow")
lines(WTI_Roll_SST_Sample2@forecast$VaR[,3],type="l",col="green")
lines(WTI_Roll_SGED_Sample2@forecast$VaR[,3],type="l",col="gray")
lines(WTI_Roll_JSU_Sample2@forecast$VaR[,3],type="l",col="brown")
lines(Sample2_WTI_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=WTI_Sample2_Dates[c(1,50,100,150,200,250)],las=2)

plot(WTI_Roll_JSU_Sample3@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample3",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(WTI_Roll_Normal_Sample3@forecast$VaR[,3],type="l",col="yellow")
lines(WTI_Roll_SST_Sample3@forecast$VaR[,3],type="l",col="green")
lines(WTI_Roll_SGED_Sample3@forecast$VaR[,3],type="l",col="gray")
lines(WTI_Roll_JSU_Sample3@forecast$VaR[,3],type="l",col="brown")
lines(Sample3_WTI_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=WTI_Sample3_Dates[c(1,50,100,150,200,250)],las=2)

plot(WTI_Roll_JSU_Sample4@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample4",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(WTI_Roll_Normal_Sample4@forecast$VaR[,3],type="l",col="yellow")
lines(WTI_Roll_SST_Sample4@forecast$VaR[,3],type="l",col="green")
lines(WTI_Roll_SGED_Sample4@forecast$VaR[,3],type="l",col="gray")
lines(WTI_Roll_JSU_Sample4@forecast$VaR[,3],type="l",col="brown")
lines(Sample4_WTI_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=WTI_Sample4_Dates[c(1,50,100,150,200,250)],las=2)

plot(WTI_Roll_JSU_Sample5@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample5",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(WTI_Roll_Normal_Sample5@forecast$VaR[,3],type="l",col="yellow")
lines(WTI_Roll_SST_Sample5@forecast$VaR[,3],type="l",col="green")
lines(WTI_Roll_SGED_Sample5@forecast$VaR[,3],type="l",col="gray")
lines(WTI_Roll_JSU_Sample5@forecast$VaR[,3],type="l",col="brown")
lines(Sample5_WTI_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=WTI_Sample5_Dates[c(1,50,100,150,200,250)],las=2)


### Incorporating Structural Breaks in a different way using dummies###
# 

par(mfrow=c(1,2))
plot(WTI_Roll8@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample for WTI",ylab="Returns/ VaR at 1%",xlab="",xaxt="n")
lines(WTI_Roll8@forecast$VaR[,1],type="l",col="yellow")
lines(WTI_Roll5@forecast$VaR[,1],type="l",col="green")
lines(WTI_Roll6@forecast$VaR[,1],type="l",col="gray")
lines(WTI_Roll7@forecast$VaR[,1],type="l",col="brown")
lines(WTI_Holdout_VaR[-1],type="l",col="black")
legend(x=1,y=0.14,legend=c("Returns","Normal","SST","SGED","JSU","Pearson"),col=c("blue","yellow","green","gray","brown","black"),cex=0.8,lty=1,bty="n",y.intersp = 0.7,ncol=2,x.intersp = 0.49)
axis(1,at=c(0,200,400,600,800,950),labels=format(WTI_Date[c(4100,4300,4500,4700,4900,5100)],"%b-%y"),las=2)

### Keep the plot width as 800 pixels and save as .tiff extension

plot(Brent_Roll8@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample for Brent",ylab="Returns/ VaR at 1%",xlab="",xaxt="n")
lines(Brent_Roll8@forecast$VaR[,1],type="l",col="yellow")
lines(Brent_Roll5@forecast$VaR[,1],type="l",col="green")
lines(Brent_Roll6@forecast$VaR[,1],type="l",col="gray")
lines(Brent_Roll7@forecast$VaR[,1],type="l",col="brown")
lines(Brent_Holdout_VaR[-1],type="l",col="black")
legend(x=1,y=0.14,legend=c("Returns","Normal","SST","SGED","JSU","Pearson"),col=c("blue","yellow","green","gray","brown","black"),cex=0.8,lty=1,bty="n",y.intersp = 0.7,x.intersp = 0.49,ncol=2)
axis(1,at=c(0,200,400,600,800,950),labels=format(Brent_Date[c(4001,4200,4400,4600,4800,5000)],"%b-%y"),las=2)

## Storing Pearson subsample results ###

WTI_subsample_1pc <- rbind(Holdout_Sample1_WTI_Pearson[1,],Holdout_Sample2_WTI_Pearson[1,],Holdout_Sample3_WTI_Pearson[1,],Holdout_Sample4_WTI_Pearson[1,],Holdout_Sample5_WTI_Pearson[1,])
Brent_subsample_1pc   <- rbind(Holdout_Sample1_Brent_Pearson[1,],Holdout_Sample2_Brent_Pearson[1,],Holdout_Sample3_Brent_Pearson[1,],Holdout_Sample4_Brent_Pearson[1,],Holdout_Sample5_Brent_Pearson[1,])

WTI_subsample_5pc <- rbind(Holdout_Sample1_WTI_Pearson[2,],Holdout_Sample2_WTI_Pearson[2,],Holdout_Sample3_WTI_Pearson[2,],Holdout_Sample4_WTI_Pearson[2,],Holdout_Sample5_WTI_Pearson[2,])
Brent_subsample_5pc   <- rbind(Holdout_Sample1_Brent_Pearson[2,],Holdout_Sample2_Brent_Pearson[2,],Holdout_Sample3_Brent_Pearson[2,],Holdout_Sample4_Brent_Pearson[2,],Holdout_Sample5_Brent_Pearson[2,])

WTI_subsample_10pc <- rbind(Holdout_Sample1_WTI_Pearson[3,],Holdout_Sample2_WTI_Pearson[3,],Holdout_Sample3_WTI_Pearson[3,],Holdout_Sample4_WTI_Pearson[3,],Holdout_Sample5_WTI_Pearson[3,])
Brent_subsample_10pc   <- rbind(Holdout_Sample1_Brent_Pearson[3,],Holdout_Sample2_Brent_Pearson[3,],Holdout_Sample3_Brent_Pearson[3,],Holdout_Sample4_Brent_Pearson[3,],Holdout_Sample5_Brent_Pearson[3,])

WTI_subsample_pvalues_1pc <- rbind(Holdout_Sample1_WTI_Pearson_pvalues[1,],Holdout_Sample2_WTI_Pearson_pvalues[1,],Holdout_Sample3_WTI_Pearson_pvalues[1,],Holdout_Sample4_WTI_Pearson_pvalues[1,],Holdout_Sample5_WTI_Pearson_pvalues[1,])
Brent_subsample_pvalues_1pc   <- rbind(Holdout_Sample1_Brent_Pearson_pvalues[1,],Holdout_Sample2_Brent_Pearson_pvalues[1,],Holdout_Sample3_Brent_Pearson_pvalues[1,],Holdout_Sample4_Brent_Pearson_pvalues[1,],Holdout_Sample5_Brent_Pearson_pvalues[1,])

WTI_subsample_pvalues_5pc <- rbind(Holdout_Sample1_WTI_Pearson_pvalues[2,],Holdout_Sample2_WTI_Pearson_pvalues[2,],Holdout_Sample3_WTI_Pearson_pvalues[2,],Holdout_Sample4_WTI_Pearson_pvalues[2,],Holdout_Sample5_WTI_Pearson_pvalues[2,])
Brent_subsample_pvalues_5pc   <- rbind(Holdout_Sample1_Brent_Pearson_pvalues[2,],Holdout_Sample2_Brent_Pearson_pvalues[2,],Holdout_Sample3_Brent_Pearson_pvalues[2,],Holdout_Sample4_Brent_Pearson_pvalues[2,],Holdout_Sample5_Brent_Pearson_pvalues[2,])

WTI_subsample_pvalues_10pc <- rbind(Holdout_Sample1_WTI_Pearson_pvalues[3,],Holdout_Sample2_WTI_Pearson_pvalues[3,],Holdout_Sample3_WTI_Pearson_pvalues[3,],Holdout_Sample4_WTI_Pearson_pvalues[3,],Holdout_Sample5_WTI_Pearson_pvalues[3,])
Brent_subsample_pvalues_10pc   <- rbind(Holdout_Sample1_Brent_Pearson_pvalues[3,],Holdout_Sample2_Brent_Pearson_pvalues[3,],Holdout_Sample3_Brent_Pearson_pvalues[3,],Holdout_Sample4_Brent_Pearson_pvalues[3,],Holdout_Sample5_Brent_Pearson_pvalues[3,])

#### VaRDuration Test for subsamples

WTI_Dur_Normal_Sample1 <- VaRDurTest(0.10,WTI_Roll_Normal_Sample1@forecast$VaR[,4],WTI_Roll_Normal_Sample1@forecast$VaR[,3])
WTI_Dur_Normal_Sample2 <- VaRDurTest(0.10,WTI_Roll_Normal_Sample2@forecast$VaR[,4],WTI_Roll_Normal_Sample2@forecast$VaR[,3])
WTI_Dur_Normal_Sample3 <- VaRDurTest(0.10,WTI_Roll_Normal_Sample3@forecast$VaR[,4],WTI_Roll_Normal_Sample3@forecast$VaR[,3])
WTI_Dur_Normal_Sample4 <- VaRDurTest(0.10,WTI_Roll_Normal_Sample4@forecast$VaR[,4],WTI_Roll_Normal_Sample4@forecast$VaR[,3])
WTI_Dur_Normal_Sample5 <- VaRDurTest(0.10,WTI_Roll_Normal_Sample5@forecast$VaR[,4],WTI_Roll_Normal_Sample5@forecast$VaR[,3])

Brent_Dur_Normal_Sample1 <- VaRDurTest(0.10,Brent_Roll_Normal_Sample1@forecast$VaR[,4],Brent_Roll_Normal_Sample1@forecast$VaR[,3])
Brent_Dur_Normal_Sample2 <- VaRDurTest(0.10,Brent_Roll_Normal_Sample2@forecast$VaR[,4],Brent_Roll_Normal_Sample2@forecast$VaR[,3])
Brent_Dur_Normal_Sample3 <- VaRDurTest(0.10,Brent_Roll_Normal_Sample3@forecast$VaR[,4],Brent_Roll_Normal_Sample3@forecast$VaR[,3])
Brent_Dur_Normal_Sample4 <- VaRDurTest(0.10,Brent_Roll_Normal_Sample4@forecast$VaR[,4],Brent_Roll_Normal_Sample4@forecast$VaR[,3])
Brent_Dur_Normal_Sample5 <- VaRDurTest(0.10,Brent_Roll_Normal_Sample5@forecast$VaR[,4],Brent_Roll_Normal_Sample5@forecast$VaR[,3])

WTI_Dur_SST_Sample1 <- VaRDurTest(0.10,WTI_Roll_SST_Sample1@forecast$VaR[,4],WTI_Roll_SST_Sample1@forecast$VaR[,3])
WTI_Dur_SST_Sample2 <- VaRDurTest(0.10,WTI_Roll_SST_Sample2@forecast$VaR[,4],WTI_Roll_SST_Sample2@forecast$VaR[,3])
WTI_Dur_SST_Sample3 <- VaRDurTest(0.10,WTI_Roll_SST_Sample3@forecast$VaR[,4],WTI_Roll_SST_Sample3@forecast$VaR[,3])
WTI_Dur_SST_Sample4 <- VaRDurTest(0.10,WTI_Roll_SST_Sample4@forecast$VaR[,4],WTI_Roll_SST_Sample4@forecast$VaR[,3])
WTI_Dur_SST_Sample5 <- VaRDurTest(0.10,WTI_Roll_SST_Sample5@forecast$VaR[,4],WTI_Roll_SST_Sample5@forecast$VaR[,3])

Brent_Dur_SST_Sample1 <- VaRDurTest(0.10,Brent_Roll_SST_Sample1@forecast$VaR[,4],Brent_Roll_SST_Sample1@forecast$VaR[,3])
Brent_Dur_SST_Sample2 <- VaRDurTest(0.10,Brent_Roll_SST_Sample2@forecast$VaR[,4],Brent_Roll_SST_Sample2@forecast$VaR[,3])
Brent_Dur_SST_Sample3 <- VaRDurTest(0.10,Brent_Roll_SST_Sample3@forecast$VaR[,4],Brent_Roll_SST_Sample3@forecast$VaR[,3])
Brent_Dur_SST_Sample4 <- VaRDurTest(0.10,Brent_Roll_SST_Sample4@forecast$VaR[,4],Brent_Roll_SST_Sample4@forecast$VaR[,3])
Brent_Dur_SST_Sample5 <- VaRDurTest(0.10,Brent_Roll_SST_Sample5@forecast$VaR[,4],Brent_Roll_SST_Sample5@forecast$VaR[,3])

WTI_Dur_SGED_Sample1 <- VaRDurTest(0.10,WTI_Roll_SGED_Sample1@forecast$VaR[,4],WTI_Roll_SGED_Sample1@forecast$VaR[,3])
WTI_Dur_SGED_Sample2 <- VaRDurTest(0.10,WTI_Roll_SGED_Sample2@forecast$VaR[,4],WTI_Roll_SGED_Sample2@forecast$VaR[,3])
WTI_Dur_SGED_Sample3 <- VaRDurTest(0.10,WTI_Roll_SGED_Sample3@forecast$VaR[,4],WTI_Roll_SGED_Sample3@forecast$VaR[,3])
WTI_Dur_SGED_Sample4 <- VaRDurTest(0.10,WTI_Roll_SGED_Sample4@forecast$VaR[,4],WTI_Roll_SGED_Sample4@forecast$VaR[,3])
WTI_Dur_SGED_Sample5 <- VaRDurTest(0.10,WTI_Roll_SGED_Sample5@forecast$VaR[,4],WTI_Roll_SGED_Sample5@forecast$VaR[,3])

Brent_Dur_SGED_Sample1 <- VaRDurTest(0.10,Brent_Roll_SGED_Sample1@forecast$VaR[,4],Brent_Roll_SGED_Sample1@forecast$VaR[,3])
Brent_Dur_SGED_Sample2 <- VaRDurTest(0.10,Brent_Roll_SGED_Sample2@forecast$VaR[,4],Brent_Roll_SGED_Sample2@forecast$VaR[,3])
Brent_Dur_SGED_Sample3 <- VaRDurTest(0.10,Brent_Roll_SGED_Sample3@forecast$VaR[,4],Brent_Roll_SGED_Sample3@forecast$VaR[,3])
Brent_Dur_SGED_Sample4 <- VaRDurTest(0.10,Brent_Roll_SGED_Sample4@forecast$VaR[,4],Brent_Roll_SGED_Sample4@forecast$VaR[,3])
Brent_Dur_SGED_Sample5 <- VaRDurTest(0.10,Brent_Roll_SGED_Sample5@forecast$VaR[,4],Brent_Roll_SGED_Sample5@forecast$VaR[,3])

WTI_Dur_JSU_Sample1 <- VaRDurTest(0.10,WTI_Roll_JSU_Sample1@forecast$VaR[,4],WTI_Roll_JSU_Sample1@forecast$VaR[,3])
WTI_Dur_JSU_Sample2 <- VaRDurTest(0.10,WTI_Roll_JSU_Sample2@forecast$VaR[,4],WTI_Roll_JSU_Sample2@forecast$VaR[,3])
WTI_Dur_JSU_Sample3 <- VaRDurTest(0.10,WTI_Roll_JSU_Sample3@forecast$VaR[,4],WTI_Roll_JSU_Sample3@forecast$VaR[,3])
WTI_Dur_JSU_Sample4 <- VaRDurTest(0.10,WTI_Roll_JSU_Sample4@forecast$VaR[,4],WTI_Roll_JSU_Sample4@forecast$VaR[,3])
WTI_Dur_JSU_Sample5 <- VaRDurTest(0.10,WTI_Roll_JSU_Sample5@forecast$VaR[,4],WTI_Roll_JSU_Sample5@forecast$VaR[,3])

Brent_Dur_JSU_Sample1 <- VaRDurTest(0.10,Brent_Roll_JSU_Sample1@forecast$VaR[,4],Brent_Roll_JSU_Sample1@forecast$VaR[,3])
Brent_Dur_JSU_Sample2 <- VaRDurTest(0.10,Brent_Roll_JSU_Sample2@forecast$VaR[,4],Brent_Roll_JSU_Sample2@forecast$VaR[,3])
Brent_Dur_JSU_Sample3 <- VaRDurTest(0.10,Brent_Roll_JSU_Sample3@forecast$VaR[,4],Brent_Roll_JSU_Sample3@forecast$VaR[,3])
Brent_Dur_JSU_Sample4 <- VaRDurTest(0.10,Brent_Roll_JSU_Sample4@forecast$VaR[,4],Brent_Roll_JSU_Sample4@forecast$VaR[,3])
Brent_Dur_JSU_Sample5 <- VaRDurTest(0.10,Brent_Roll_JSU_Sample5@forecast$VaR[,4],Brent_Roll_JSU_Sample5@forecast$VaR[,3])

## Store the results

WTI_Dur_Subsample_stat <- matrix(c(WTI_Dur_Normal_Sample1$rLL,
                                 WTI_Dur_SST_Sample1$rLL,
                                 WTI_Dur_SGED_Sample1$rLL,
                                 WTI_Dur_JSU_Sample1$rLL,
                                 WTI_Dur_Normal_Sample2$rLL,
                                 WTI_Dur_SST_Sample2$rLL,
                                 WTI_Dur_SGED_Sample2$rLL,
                                 WTI_Dur_JSU_Sample2$rLL,
                                 WTI_Dur_Normal_Sample3$rLL,
                                 WTI_Dur_SST_Sample3$rLL,
                                 WTI_Dur_SGED_Sample3$rLL,
                                 WTI_Dur_JSU_Sample3$rLL,
                                 WTI_Dur_Normal_Sample4$rLL,
                                 WTI_Dur_SST_Sample4$rLL,
                                 WTI_Dur_SGED_Sample4$rLL,
                                 WTI_Dur_JSU_Sample4$rLL,
                                 WTI_Dur_Normal_Sample5$rLL,
                                 WTI_Dur_SST_Sample5$rLL,
                                 WTI_Dur_SGED_Sample5$rLL,
                                 WTI_Dur_JSU_Sample5$rLL),nrow=4)

Brent_Dur_Subsample_stat <- matrix(c(Brent_Dur_Normal_Sample1$rLL,
                                 Brent_Dur_SST_Sample1$rLL,
                                 Brent_Dur_SGED_Sample1$rLL,
                                 Brent_Dur_JSU_Sample1$rLL,
                                 Brent_Dur_Normal_Sample2$rLL,
                                 Brent_Dur_SST_Sample2$rLL,
                                 Brent_Dur_SGED_Sample2$rLL,
                                 Brent_Dur_JSU_Sample2$rLL,
                                 Brent_Dur_Normal_Sample3$rLL,
                                 Brent_Dur_SST_Sample3$rLL,
                                 Brent_Dur_SGED_Sample3$rLL,
                                 Brent_Dur_JSU_Sample3$rLL,
                                 Brent_Dur_Normal_Sample4$rLL,
                                 Brent_Dur_SST_Sample4$rLL,
                                 Brent_Dur_SGED_Sample4$rLL,
                                 Brent_Dur_JSU_Sample4$rLL,
                                 Brent_Dur_Normal_Sample5$rLL,
                                 Brent_Dur_SST_Sample5$rLL,
                                 Brent_Dur_SGED_Sample5$rLL,
                                 Brent_Dur_JSU_Sample5$rLL),nrow=4)

WTI_Dur_Subsample_pvalues <- matrix(c(WTI_Dur_Normal_Sample1$LRp,
                                 WTI_Dur_SST_Sample1$LRp,
                                 WTI_Dur_SGED_Sample1$LRp,
                                 WTI_Dur_JSU_Sample1$LRp,
                                 WTI_Dur_Normal_Sample2$LRp,
                                 WTI_Dur_SST_Sample2$LRp,
                                 WTI_Dur_SGED_Sample2$LRp,
                                 WTI_Dur_JSU_Sample2$LRp,
                                 WTI_Dur_Normal_Sample3$LRp,
                                 WTI_Dur_SST_Sample3$LRp,
                                 WTI_Dur_SGED_Sample3$LRp,
                                 WTI_Dur_JSU_Sample3$LRp,
                                 WTI_Dur_Normal_Sample4$LRp,
                                 WTI_Dur_SST_Sample4$LRp,
                                 WTI_Dur_SGED_Sample4$LRp,
                                 WTI_Dur_JSU_Sample4$LRp,
                                 WTI_Dur_Normal_Sample5$LRp,
                                 WTI_Dur_SST_Sample5$LRp,
                                 WTI_Dur_SGED_Sample5$LRp,
                                 WTI_Dur_JSU_Sample5$LRp),nrow=4)

Brent_Dur_Subsample_pvalues <- matrix(c(Brent_Dur_Normal_Sample1$LRp,
                               Brent_Dur_SST_Sample1$LRp,
                               Brent_Dur_SGED_Sample1$LRp,
                               Brent_Dur_JSU_Sample1$LRp,
                               Brent_Dur_Normal_Sample2$LRp,
                               Brent_Dur_SST_Sample2$LRp,
                               Brent_Dur_SGED_Sample2$LRp,
                               Brent_Dur_JSU_Sample2$LRp,
                               Brent_Dur_Normal_Sample3$LRp,
                               Brent_Dur_SST_Sample3$LRp,
                               Brent_Dur_SGED_Sample3$LRp,
                               Brent_Dur_JSU_Sample3$LRp,
                               Brent_Dur_Normal_Sample4$LRp,
                               Brent_Dur_SST_Sample4$LRp,
                               Brent_Dur_SGED_Sample4$LRp,
                               Brent_Dur_JSU_Sample4$LRp,
                               Brent_Dur_Normal_Sample5$LRp,
                               Brent_Dur_SST_Sample5$LRp,
                               Brent_Dur_SGED_Sample5$LRp,
                               Brent_Dur_JSU_Sample5$LRp),nrow=4)


### Expected Shortfall Calculations ####

## Assuming Normal Distribution

WTI_index_normal1 <- which(WTI_Roll8@forecast$VaR[,1]>WTI_Roll8@forecast$VaR[,4])
WTI_VaR_values_normal1 <- WTI_Roll8@forecast$VaR[WTI_index_normal1,1]
WTI_return_values_normal1 <- WTI_Roll8@forecast$VaR[WTI_index_normal1,4]
WTI_ES_normal1 <- (WTI_return_values_normal1/WTI_VaR_values_normal1)
WTI_ES_Measure1_normal1 <- mean(WTI_ES_normal1)
WTI_ES_Measure2_normal1 <- max(WTI_ES_normal1)

WTI_index_normal5 <- which(WTI_Roll8@forecast$VaR[,2]>WTI_Roll8@forecast$VaR[,4])
WTI_VaR_values_normal5 <- WTI_Roll8@forecast$VaR[WTI_index_normal5,2]
WTI_return_values_normal5 <- WTI_Roll8@forecast$VaR[WTI_index_normal5,4]
WTI_ES_normal5 <- (WTI_return_values_normal5/WTI_VaR_values_normal5)
WTI_ES_Measure1_normal5 <- mean(WTI_ES_normal5)
WTI_ES_Measure2_normal5 <- max(WTI_ES_normal5)

WTI_index_normal10 <- which(WTI_Roll8@forecast$VaR[,3]>WTI_Roll8@forecast$VaR[,4])
WTI_VaR_values_normal10 <- WTI_Roll8@forecast$VaR[WTI_index_normal10,3]
WTI_return_values_normal10 <- WTI_Roll8@forecast$VaR[WTI_index_normal10,4]
WTI_ES_normal10 <- (WTI_return_values_normal10/WTI_VaR_values_normal10)
WTI_ES_Measure1_normal10 <- mean(WTI_ES_normal10)
WTI_ES_Measure2_normal10 <- max(WTI_ES_normal10)

Brent_index_normal1 <- which(Brent_Roll8@forecast$VaR[,1]>Brent_Roll8@forecast$VaR[,4])
Brent_VaR_values_normal1 <- Brent_Roll8@forecast$VaR[Brent_index_normal1,1]
Brent_return_values_normal1 <- Brent_Roll8@forecast$VaR[Brent_index_normal1,4]
Brent_ES_normal1 <- (Brent_return_values_normal1/Brent_VaR_values_normal1)
Brent_ES_Measure1_normal1 <- mean(Brent_ES_normal1)
Brent_ES_Measure2_normal1 <- max(Brent_ES_normal1)

Brent_index_normal5 <- which(Brent_Roll8@forecast$VaR[,2]>Brent_Roll8@forecast$VaR[,4])
Brent_VaR_values_normal5 <- Brent_Roll8@forecast$VaR[Brent_index_normal5,2]
Brent_return_values_normal5 <- Brent_Roll8@forecast$VaR[Brent_index_normal5,4]
Brent_ES_normal5 <- (Brent_return_values_normal5/Brent_VaR_values_normal5)
Brent_ES_Measure1_normal5 <- mean(Brent_ES_normal5)
Brent_ES_Measure2_normal5 <- max(Brent_ES_normal5)

Brent_index_normal10 <- which(Brent_Roll8@forecast$VaR[,3]>Brent_Roll8@forecast$VaR[,4])
Brent_VaR_values_normal10 <- Brent_Roll8@forecast$VaR[Brent_index_normal10,3]
Brent_return_values_normal10 <- Brent_Roll8@forecast$VaR[Brent_index_normal10,4]
Brent_ES_normal10 <- (Brent_return_values_normal10/Brent_VaR_values_normal10)
Brent_ES_Measure1_normal10 <- mean(Brent_ES_normal10)
Brent_ES_Measure2_normal10 <- max(Brent_ES_normal10)

## Assuming SST Distribution

WTI_index_SST1 <- which(WTI_Roll5@forecast$VaR[,1]>WTI_Roll5@forecast$VaR[,4])
WTI_VaR_values_SST1 <- WTI_Roll5@forecast$VaR[WTI_index_SST1,1]
WTI_return_values_SST1 <- WTI_Roll5@forecast$VaR[WTI_index_SST1,4]
WTI_ES_SST1 <- (WTI_return_values_SST1/WTI_VaR_values_SST1)
WTI_ES_Measure1_SST1 <- mean(WTI_ES_SST1)
WTI_ES_Measure2_SST1 <- max(WTI_ES_SST1)

WTI_index_SST5 <- which(WTI_Roll5@forecast$VaR[,2]>WTI_Roll5@forecast$VaR[,4])
WTI_VaR_values_SST5 <- WTI_Roll5@forecast$VaR[WTI_index_SST5,2]
WTI_return_values_SST5 <- WTI_Roll5@forecast$VaR[WTI_index_SST5,4]
WTI_ES_SST5 <- (WTI_return_values_SST5/WTI_VaR_values_SST5)
WTI_ES_Measure1_SST5 <- mean(WTI_ES_SST5)
WTI_ES_Measure2_SST5 <- max(WTI_ES_SST5)

WTI_index_SST10 <- which(WTI_Roll5@forecast$VaR[,3]>WTI_Roll5@forecast$VaR[,4])
WTI_VaR_values_SST10 <- WTI_Roll5@forecast$VaR[WTI_index_SST10,3]
WTI_return_values_SST10 <- WTI_Roll5@forecast$VaR[WTI_index_SST10,4]
WTI_ES_SST10 <- (WTI_return_values_SST10/WTI_VaR_values_SST10)
WTI_ES_Measure1_SST10 <- mean(WTI_ES_SST10)
WTI_ES_Measure2_SST10 <- max(WTI_ES_SST10)

Brent_index_SST1 <- which(Brent_Roll5@forecast$VaR[,1]>Brent_Roll5@forecast$VaR[,4])
Brent_VaR_values_SST1 <- Brent_Roll5@forecast$VaR[Brent_index_SST1,1]
Brent_return_values_SST1 <- Brent_Roll5@forecast$VaR[Brent_index_SST1,4]
Brent_ES_SST1 <- (Brent_return_values_SST1/Brent_VaR_values_SST1)
Brent_ES_Measure1_SST1 <- mean(Brent_ES_SST1)
Brent_ES_Measure2_SST1 <- max(Brent_ES_SST1)

Brent_index_SST5 <- which(Brent_Roll5@forecast$VaR[,2]>Brent_Roll5@forecast$VaR[,4])
Brent_VaR_values_SST5 <- Brent_Roll5@forecast$VaR[Brent_index_SST5,2]
Brent_return_values_SST5 <- Brent_Roll5@forecast$VaR[Brent_index_SST5,4]
Brent_ES_SST5 <- (Brent_return_values_SST5/Brent_VaR_values_SST5)
Brent_ES_Measure1_SST5 <- mean(Brent_ES_SST5)
Brent_ES_Measure2_SST5 <- max(Brent_ES_SST5)

Brent_index_SST10 <- which(Brent_Roll5@forecast$VaR[,3]>Brent_Roll5@forecast$VaR[,4])
Brent_VaR_values_SST10 <- Brent_Roll5@forecast$VaR[Brent_index_SST10,3]
Brent_return_values_SST10 <- Brent_Roll5@forecast$VaR[Brent_index_SST10,4]
Brent_ES_SST10 <- (Brent_return_values_SST10/Brent_VaR_values_SST10)
Brent_ES_Measure1_SST10 <- mean(Brent_ES_SST10)
Brent_ES_Measure2_SST10 <- max(Brent_ES_SST10)


## Assuming SGED Distribution

WTI_index_SGED1 <- which(WTI_Roll6@forecast$VaR[,1]>WTI_Roll6@forecast$VaR[,4])
WTI_VaR_values_SGED1 <- WTI_Roll6@forecast$VaR[WTI_index_SGED1,1]
WTI_return_values_SGED1 <- WTI_Roll6@forecast$VaR[WTI_index_SGED1,4]
WTI_ES_SGED1 <- (WTI_return_values_SGED1/WTI_VaR_values_SGED1)
WTI_ES_Measure1_SGED1 <- mean(WTI_ES_SGED1)
WTI_ES_Measure2_SGED1 <- max(WTI_ES_SGED1)

WTI_index_SGED5 <- which(WTI_Roll6@forecast$VaR[,2]>WTI_Roll6@forecast$VaR[,4])
WTI_VaR_values_SGED5 <- WTI_Roll6@forecast$VaR[WTI_index_SGED5,2]
WTI_return_values_SGED5 <- WTI_Roll6@forecast$VaR[WTI_index_SGED5,4]
WTI_ES_SGED5 <- (WTI_return_values_SGED5/WTI_VaR_values_SGED5)
WTI_ES_Measure1_SGED5 <- mean(WTI_ES_SGED5)
WTI_ES_Measure2_SGED5 <- max(WTI_ES_SGED5)

WTI_index_SGED10 <- which(WTI_Roll6@forecast$VaR[,3]>WTI_Roll6@forecast$VaR[,4])
WTI_VaR_values_SGED10 <- WTI_Roll6@forecast$VaR[WTI_index_SGED10,3]
WTI_return_values_SGED10 <- WTI_Roll6@forecast$VaR[WTI_index_SGED10,4]
WTI_ES_SGED10 <- (WTI_return_values_SGED10/WTI_VaR_values_SGED10)
WTI_ES_Measure1_SGED10 <- mean(WTI_ES_SGED10)
WTI_ES_Measure2_SGED10 <- max(WTI_ES_SGED10)

Brent_index_SGED1 <- which(Brent_Roll6@forecast$VaR[,1]>Brent_Roll6@forecast$VaR[,4])
Brent_VaR_values_SGED1 <- Brent_Roll6@forecast$VaR[Brent_index_SGED1,1]
Brent_return_values_SGED1 <- Brent_Roll6@forecast$VaR[Brent_index_SGED1,4]
Brent_ES_SGED1 <- (Brent_return_values_SGED1/Brent_VaR_values_SGED1)
Brent_ES_Measure1_SGED1 <- mean(Brent_ES_SGED1)
Brent_ES_Measure2_SGED1 <- max(Brent_ES_SGED1)

Brent_index_SGED5 <- which(Brent_Roll6@forecast$VaR[,2]>Brent_Roll6@forecast$VaR[,4])
Brent_VaR_values_SGED5 <- Brent_Roll6@forecast$VaR[Brent_index_SGED5,2]
Brent_return_values_SGED5 <- Brent_Roll6@forecast$VaR[Brent_index_SGED5,4]
Brent_ES_SGED5 <- (Brent_return_values_SGED5/Brent_VaR_values_SGED5)
Brent_ES_Measure1_SGED5 <- mean(Brent_ES_SGED5)
Brent_ES_Measure2_SGED5 <- max(Brent_ES_SGED5)

Brent_index_SGED10 <- which(Brent_Roll6@forecast$VaR[,3]>Brent_Roll6@forecast$VaR[,4])
Brent_VaR_values_SGED10 <- Brent_Roll6@forecast$VaR[Brent_index_SGED10,3]
Brent_return_values_SGED10 <- Brent_Roll6@forecast$VaR[Brent_index_SGED10,4]
Brent_ES_SGED10 <- (Brent_return_values_SGED10/Brent_VaR_values_SGED10)
Brent_ES_Measure1_SGED10 <- mean(Brent_ES_SGED10)
Brent_ES_Measure2_SGED10 <- max(Brent_ES_SGED10)


## Assuming JSU Distribution

WTI_index_JSU1 <- which(WTI_Roll7@forecast$VaR[,1]>WTI_Roll7@forecast$VaR[,4])
WTI_VaR_values_JSU1 <- WTI_Roll7@forecast$VaR[WTI_index_JSU1,1]
WTI_return_values_JSU1 <- WTI_Roll7@forecast$VaR[WTI_index_JSU1,4]
WTI_ES_JSU1 <- (WTI_return_values_JSU1/WTI_VaR_values_JSU1)
WTI_ES_Measure1_JSU1 <- mean(WTI_ES_JSU1)
WTI_ES_Measure2_JSU1 <- max(WTI_ES_JSU1)

WTI_index_JSU5 <- which(WTI_Roll7@forecast$VaR[,2]>WTI_Roll7@forecast$VaR[,4])
WTI_VaR_values_JSU5 <- WTI_Roll7@forecast$VaR[WTI_index_JSU5,2]
WTI_return_values_JSU5 <- WTI_Roll7@forecast$VaR[WTI_index_JSU5,4]
WTI_ES_JSU5 <- (WTI_return_values_JSU5/WTI_VaR_values_JSU5)
WTI_ES_Measure1_JSU5 <- mean(WTI_ES_JSU5)
WTI_ES_Measure2_JSU5 <- max(WTI_ES_JSU5)

WTI_index_JSU10 <- which(WTI_Roll7@forecast$VaR[,3]>WTI_Roll7@forecast$VaR[,4])
WTI_VaR_values_JSU10 <- WTI_Roll7@forecast$VaR[WTI_index_JSU10,3]
WTI_return_values_JSU10 <- WTI_Roll7@forecast$VaR[WTI_index_JSU10,4]
WTI_ES_JSU10 <- (WTI_return_values_JSU10/WTI_VaR_values_JSU10)
WTI_ES_Measure1_JSU10 <- mean(WTI_ES_JSU10)
WTI_ES_Measure2_JSU10 <- max(WTI_ES_JSU10)

Brent_index_JSU1 <- which(Brent_Roll7@forecast$VaR[,1]>Brent_Roll7@forecast$VaR[,4])
Brent_VaR_values_JSU1 <- Brent_Roll7@forecast$VaR[Brent_index_JSU1,1]
Brent_return_values_JSU1 <- Brent_Roll7@forecast$VaR[Brent_index_JSU1,4]
Brent_ES_JSU1 <- (Brent_return_values_JSU1/Brent_VaR_values_JSU1)
Brent_ES_Measure1_JSU1 <- mean(Brent_ES_JSU1)
Brent_ES_Measure2_JSU1 <- max(Brent_ES_JSU1)

Brent_index_JSU5 <- which(Brent_Roll7@forecast$VaR[,2]>Brent_Roll7@forecast$VaR[,4])
Brent_VaR_values_JSU5 <- Brent_Roll7@forecast$VaR[Brent_index_JSU5,2]
Brent_return_values_JSU5 <- Brent_Roll7@forecast$VaR[Brent_index_JSU5,4]
Brent_ES_JSU5 <- (Brent_return_values_JSU5/Brent_VaR_values_JSU5)
Brent_ES_Measure1_JSU5 <- mean(Brent_ES_JSU5)
Brent_ES_Measure2_JSU5 <- max(Brent_ES_JSU5)

Brent_index_JSU10 <- which(Brent_Roll7@forecast$VaR[,3]>Brent_Roll7@forecast$VaR[,4])
Brent_VaR_values_JSU10 <- Brent_Roll7@forecast$VaR[Brent_index_JSU10,3]
Brent_return_values_JSU10 <- Brent_Roll7@forecast$VaR[Brent_index_JSU10,4]
Brent_ES_JSU10 <- (Brent_return_values_JSU10/Brent_VaR_values_JSU10)
Brent_ES_Measure1_JSU10 <- mean(Brent_ES_JSU10)
Brent_ES_Measure2_JSU10 <- max(Brent_ES_JSU10)

### For Pearson

WTI_index_Pearson1 <- which(WTI_Holdout_VaR[-1]>WTI_Roll7@forecast$VaR[1:999,4])
WTI_VaR_values_Pearson1 <- WTI_Holdout_VaR[WTI_index_Pearson1]
WTI_return_values_Pearson1 <- WTI_Roll7@forecast$VaR[WTI_index_Pearson1,4]
WTI_ES_Pearson1 <- (WTI_return_values_Pearson1/WTI_VaR_values_Pearson1)
WTI_ES_Measure1_Pearson1 <- mean(WTI_ES_Pearson1)
WTI_ES_Measure2_Pearson1 <- max(WTI_ES_Pearson1)

WTI_index_Pearson5 <- which(WTI_Holdout_VaR[-1]>WTI_Roll7@forecast$VaR[1:999,4])
WTI_VaR_values_Pearson5 <- WTI_Holdout_VaR[WTI_index_Pearson5]
WTI_return_values_Pearson5 <- WTI_Roll7@forecast$VaR[WTI_index_Pearson5,4]
WTI_ES_Pearson5 <- (WTI_return_values_Pearson5/WTI_VaR_values_Pearson5)
WTI_ES_Measure1_Pearson5 <- mean(WTI_ES_Pearson5)
WTI_ES_Measure2_Pearson5 <- max(WTI_ES_Pearson5)

WTI_index_Pearson10 <- which(WTI_Holdout_VaR[-1]>WTI_Roll7@forecast$VaR[1:999,4])
WTI_VaR_values_Pearson10 <- WTI_Holdout_VaR[WTI_index_Pearson10]
WTI_return_values_Pearson10 <- WTI_Roll7@forecast$VaR[WTI_index_Pearson10,4]
WTI_ES_Pearson10 <- (WTI_return_values_Pearson10/WTI_VaR_values_Pearson10)
WTI_ES_Measure1_Pearson10 <- mean(WTI_ES_Pearson10)
WTI_ES_Measure2_Pearson10 <- max(WTI_ES_Pearson10)

Brent_index_Pearson1 <- which(Brent_Holdout_VaR[-1]>Brent_Roll7@forecast$VaR[1:999,4])
Brent_VaR_values_Pearson1 <- Brent_Holdout_VaR[Brent_index_Pearson1]
Brent_return_values_Pearson1 <- Brent_Roll7@forecast$VaR[Brent_index_Pearson1,4]
Brent_ES_Pearson1 <- (Brent_return_values_Pearson1/Brent_VaR_values_Pearson1)
Brent_ES_Measure1_Pearson1 <- mean(Brent_ES_Pearson1)
Brent_ES_Measure2_Pearson1 <- max(Brent_ES_Pearson1)

Brent_index_Pearson5 <- which(Brent_Holdout_VaR[-1]>Brent_Roll7@forecast$VaR[1:999,4])
Brent_VaR_values_Pearson5 <- Brent_Holdout_VaR[Brent_index_Pearson5]
Brent_return_values_Pearson5 <- Brent_Roll7@forecast$VaR[Brent_index_Pearson5,4]
Brent_ES_Pearson5 <- (Brent_return_values_Pearson5/Brent_VaR_values_Pearson5)
Brent_ES_Measure1_Pearson5 <- mean(Brent_ES_Pearson5)
Brent_ES_Measure2_Pearson5 <- max(Brent_ES_Pearson5)

Brent_index_Pearson10 <- which(Brent_Holdout_VaR[-1]>Brent_Roll7@forecast$VaR[1:999,4])
Brent_VaR_values_Pearson10 <- Brent_Holdout_VaR[Brent_index_Pearson10]
Brent_return_values_Pearson10 <- Brent_Roll7@forecast$VaR[Brent_index_Pearson10,4]
Brent_ES_Pearson10 <- (Brent_return_values_Pearson10/Brent_VaR_values_Pearson10)
Brent_ES_Measure1_Pearson10 <- mean(Brent_ES_Pearson10)
Brent_ES_Measure2_Pearson10 <- max(Brent_ES_Pearson10)

### Store for Pearson

Pearson_Full_ES_Measure1 <- matrix(c(WTI_ES_Measure1_Pearson1,WTI_ES_Measure1_Pearson5,WTI_ES_Measure1_Pearson10,
                                     Brent_ES_Measure1_Pearson1,Brent_ES_Measure1_Pearson5,Brent_ES_Measure1_Pearson10),nrow=3)

Pearson_Full_ES_Measure2 <- matrix(c(WTI_ES_Measure2_Pearson1,WTI_ES_Measure2_Pearson5,WTI_ES_Measure2_Pearson10,
                                     Brent_ES_Measure2_Pearson1,Brent_ES_Measure2_Pearson5,Brent_ES_Measure2_Pearson10),nrow=3)

rownames(Pearson_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(Pearson_Full_ES_Measure2) <- c("1%","5%","10%")
colnames(Pearson_Full_ES_Measure1) <- c("WTI","Brent")
colnames(Pearson_Full_ES_Measure2) <- c("WTI","Brent")

#### Store the values

WTI_Full_ES_Measure1 <- matrix(c(WTI_ES_Measure1_normal1,WTI_ES_Measure1_SST1,WTI_ES_Measure1_SGED1,WTI_ES_Measure1_JSU1,
                                  WTI_ES_Measure1_normal5,WTI_ES_Measure1_SST5,WTI_ES_Measure1_SGED5,WTI_ES_Measure1_JSU5,
                                  WTI_ES_Measure1_normal10,WTI_ES_Measure1_SST10,WTI_ES_Measure1_SGED10,WTI_ES_Measure1_JSU10),nrow=3,byrow=TRUE)
  
  
Brent_Full_ES_Measure1 <- matrix(c(Brent_ES_Measure1_normal1,Brent_ES_Measure1_SST1,Brent_ES_Measure1_SGED1,Brent_ES_Measure1_JSU1,
                                   Brent_ES_Measure1_normal5,Brent_ES_Measure1_SST5,Brent_ES_Measure1_SGED5,Brent_ES_Measure1_JSU5,
                                   Brent_ES_Measure1_normal10,Brent_ES_Measure1_SST10,Brent_ES_Measure1_SGED10,Brent_ES_Measure1_JSU10),nrow=3,byrow=TRUE)

WTI_Full_ES_Measure2 <- matrix(c(WTI_ES_Measure2_normal1,WTI_ES_Measure2_SST1,WTI_ES_Measure2_SGED1,WTI_ES_Measure2_JSU1,
                                   WTI_ES_Measure2_normal5,WTI_ES_Measure2_SST5,WTI_ES_Measure2_SGED5,WTI_ES_Measure2_JSU5,
                                   WTI_ES_Measure2_normal10,WTI_ES_Measure2_SST10,WTI_ES_Measure2_SGED10,WTI_ES_Measure2_JSU10),nrow=3,byrow=TRUE)

Brent_Full_ES_Measure2 <- matrix(c(Brent_ES_Measure2_normal1,Brent_ES_Measure2_SST1,Brent_ES_Measure2_SGED1,Brent_ES_Measure2_JSU1,
                                   Brent_ES_Measure2_normal5,Brent_ES_Measure2_SST5,Brent_ES_Measure2_SGED5,Brent_ES_Measure2_JSU5,
                                   Brent_ES_Measure2_normal10,Brent_ES_Measure2_SST10,Brent_ES_Measure2_SGED10,Brent_ES_Measure2_JSU10),nrow=3,byrow=TRUE)

rownames(WTI_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(WTI_Full_ES_Measure2) <- c("1%","5%","10%")
rownames(Brent_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(Brent_Full_ES_Measure2) <- c("1%","5%","10%")

colnames(WTI_Full_ES_Measure1) <- c("Normal","SST","SGED","JSU")
colnames(WTI_Full_ES_Measure2) <- c("Normal","SST","SGED","JSU")
colnames(Brent_Full_ES_Measure1) <- c("Normal","SST","SGED","JSU")
colnames(Brent_Full_ES_Measure2) <- c("Normal","SST","SGED","JSU")

### MAPE

# WTI
check1 <- Holdout_Returns_WTI[-c(which(Holdout_Returns_WTI==0))]
check2 <- as.data.frame(x=WTI_Roll5)[-c(which(Holdout_Returns_WTI==0)),1]
check3 <- as.data.frame(x=WTI_Roll6)[-c(which(Holdout_Returns_WTI==0)),1]
check4 <- as.data.frame(x=WTI_Roll7)[-c(which(Holdout_Returns_WTI==0)),1]
check5 <- as.data.frame(x=WTI_Roll8)[-c(which(Holdout_Returns_WTI==0)),1]

MAPE_WTI <- c(mape(check1,check2),
mape(check1,check3),
mape(check1,check4),
mape(check1,check5))

## Brent
check1 <- Holdout_Returns_Brent[-c(which(Holdout_Returns_Brent==0))]
check2 <- as.data.frame(x=Brent_Roll5)[-c(which(Holdout_Returns_Brent==0)),1]
check3 <- as.data.frame(x=Brent_Roll6)[-c(which(Holdout_Returns_Brent==0)),1]
check4 <- as.data.frame(x=Brent_Roll7)[-c(which(Holdout_Returns_Brent==0)),1]
check5 <- as.data.frame(x=Brent_Roll8)[-c(which(Holdout_Returns_Brent==0)),1]

MAPE_Brent <- c(mape(check1,check2),
              mape(check1,check3),
              mape(check1,check4),
              mape(check1,check5))



#### End of all the code/lines ######  
  