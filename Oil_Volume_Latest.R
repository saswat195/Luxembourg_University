### Updating 2020-2021 data ###

## Energy Economics Paper

## Job is to compare the volatility forecasts for the oil markets under MDH ####
## using structural breaks and different distributions with volume intact

## Load the following libraries
library(forecast)
library(MCS)
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

Brent <- read.csv("~//Desktop//Latest_Research//Energy_Related//Brent_Volume.csv")
WTI <- read.csv("~//Desktop//Latest_Research//Energy_Related//WTI_Volume.csv")  
WTI_Date <- as.Date.character(WTI$Date,format = "%d/%m/%Y")
Brent_Date   <- as.Date.character(Brent$Date,format = "%d/%m/%Y")
Brent   <- data.frame(Brent_Date,Brent$Price,Brent$Volume)     
WTI <- data.frame(WTI_Date,WTI$Price,WTI$Volume)
colnames(Brent) <- c("Date","Price","Volume")
colnames(WTI) <- c("Date","Price","Volume")

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
WTI_sub5    <- WTI_xts[WTI_samples[4]:WTI_samples[5]]
WTI_sub6    <- WTI_xts[WTI_samples[5]:nrow(WTI_xts)]

Brent_sub1    <- Brent_xts[1:Brent_samples[1]]
Brent_sub2    <- Brent_xts[Brent_samples[1]:Brent_samples[2]]
Brent_sub3    <- Brent_xts[Brent_samples[2]:Brent_samples[3]]
Brent_sub4    <- Brent_xts[Brent_samples[3]:Brent_samples[4]]
Brent_sub5    <- Brent_xts[Brent_samples[4]:Brent_samples[5]]
Brent_sub6    <- Brent_xts[Brent_samples[5]:nrow(Brent_xts)]

###Dummy variables for structural breaks ####

WTI_Dummy1 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy2 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy3 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy4 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy5 <- array(data=0,dim=length(WTI_xts))
WTI_Dummy6 <- array(data=0,dim=length(WTI_xts))

WTI_Dummy1[1:WTI_samples[1]]  <- 1
WTI_Dummy2[(WTI_samples[1]+1):WTI_samples[2]] <- 1
WTI_Dummy3[(WTI_samples[2]+1):WTI_samples[3]] <- 1
WTI_Dummy4[(WTI_samples[3]+1):WTI_samples[4]] <- 1
WTI_Dummy5[(WTI_samples[4]+1):WTI_samples[5]] <- 1
WTI_Dummy6[(WTI_samples[5]+1):length(WTI_xts)] <- 1

Brent_Dummy1 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy2 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy3 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy4 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy5 <- array(data=0,dim=length(Brent_xts))
Brent_Dummy6 <- array(data=0,dim=length(Brent_xts))

Brent_Dummy1[1:Brent_samples[1]]  <- 1
Brent_Dummy2[(Brent_samples[1]+1):Brent_samples[2]] <- 1
Brent_Dummy3[(Brent_samples[2]+1):Brent_samples[3]] <- 1
Brent_Dummy4[(Brent_samples[3]+1):Brent_samples[4]] <- 1
Brent_Dummy5[(Brent_samples[4]+1):Brent_samples[5]] <- 1
Brent_Dummy6[(Brent_samples[5]+1):length(Brent_xts)] <- 1

WTI_Dummy <- cbind(WTI_Dummy1,WTI_Dummy2,WTI_Dummy3,WTI_Dummy4,WTI_Dummy5,log(as.integer(WTI$Volume[1:5537])))
Brent_Dummy <- cbind(Brent_Dummy1,Brent_Dummy2,Brent_Dummy3,Brent_Dummy4,Brent_Dummy5,log(Brent$Volume[1:5646]))

### Descriptive Statistics and Basic Tests ###

Descr_WTI <-  cbind(basicStats(WTI_sub1),basicStats(WTI_sub2),basicStats(WTI_sub3),basicStats(WTI_sub4),basicStats(WTI_sub5),basicStats(WTI_sub6))
Descr_Brent   <- cbind(basicStats(Brent_sub1),basicStats(Brent_sub2),basicStats(Brent_sub3),basicStats(Brent_sub4),basicStats(Brent_sub5),basicStats(Brent_sub6))
ADF_WTI <- cbind(adf.test(WTI_sub1),adf.test(WTI_sub2),adf.test(WTI_sub3),adf.test(WTI_sub4),adf.test(WTI_sub5),adf.test(WTI_sub6))
ADF_Brent   <- cbind(adf.test(Brent_sub1),adf.test(Brent_sub2),adf.test(Brent_sub3),adf.test(Brent_sub4),adf.test(Brent_sub5),adf.test(Brent_sub6))
pp_WTI <- cbind(pp.test(WTI_sub1),pp.test(WTI_sub2),pp.test(WTI_sub3),pp.test(WTI_sub4),pp.test(WTI_sub5),pp.test(WTI_sub6))
pp_Brent   <- cbind(pp.test(Brent_sub1),pp.test(Brent_sub2),pp.test(Brent_sub3),pp.test(Brent_sub4),pp.test(Brent_sub5),pp.test(Brent_sub6))
jb_WTI <- cbind(jarque.bera.test(WTI_sub1),jarque.bera.test(WTI_sub2),jarque.bera.test(WTI_sub3),jarque.bera.test(WTI_sub4),jarque.bera.test(WTI_sub5),jarque.bera.test(WTI_sub6))
jb_Brent   <- cbind(jarque.bera.test(Brent_sub1),jarque.bera.test(Brent_sub2),jarque.bera.test(Brent_sub3),jarque.bera.test(Brent_sub4),jarque.bera.test(Brent_sub5),jarque.bera.test(Brent_sub6))

LQ_WTI5 <- cbind(Box.test(WTI_sub1^2,lag=5,type="Lj"),Box.test(WTI_sub2^2,lag=5,type="Lj"),Box.test(WTI_sub3^2,lag=5,type="Lj"),Box.test(WTI_sub4^2,lag=5,type="Lj"),Box.test(WTI_sub5^2,lag=5,type="Lj"),Box.test(WTI_sub6^2,lag=5,type="Lj"))
LQ_Brent5   <- cbind(Box.test(Brent_sub1^2,lag=5,type="Lj"),Box.test(Brent_sub2^2,lag=5,type="Lj"),Box.test(Brent_sub3^2,lag=5,type="Lj"),Box.test(Brent_sub4^2,lag=5,type="Lj"),Box.test(Brent_sub5^2,lag=5,type="Lj"),Box.test(Brent_sub6^2,lag=5,type="Lj"))

LQ_WTI10 <- cbind(Box.test(WTI_sub1^2,lag=10,type="Lj"),Box.test(WTI_sub2^2,lag=10,type="Lj"),Box.test(WTI_sub3^2,lag=10,type="Lj"),Box.test(WTI_sub4^2,lag=10,type="Lj"),Box.test(WTI_sub5^2,lag=10,type="Lj"),Box.test(WTI_sub6^2,lag=10,type="Lj"))
LQ_Brent10   <- cbind(Box.test(Brent_sub1^2,lag=10,type="Lj"),Box.test(Brent_sub2^2,lag=10,type="Lj"),Box.test(Brent_sub3^2,lag=10,type="Lj"),Box.test(Brent_sub4^2,lag=10,type="Lj"),Box.test(Brent_sub5^2,lag=10,type="Lj"),Box.test(Brent_sub6^2,lag=10,type="Lj"))

#### Plot the Returns and Oil Prices with Breakpoints indicated

par(mfrow=c(1,2))
plot(WTI$Date,WTI$Price,type="l",col="blue",main="Structural Breaks for WTI Oil",ylab="Price/Barrel (USD)",xaxt="n",xlab="")
abline(v=WTI$Date[WTI_samples[1]],col="brown")
abline(v=WTI$Date[WTI_samples[2]],col="brown")
abline(v=WTI$Date[WTI_samples[3]],col="brown")
abline(v=WTI$Date[WTI_samples[4]],col="brown")
abline(v=WTI$Date[WTI_samples[5]],col="brown")
Date_formatted <- format(WTI$Date[c(1,WTI_samples,nrow(WTI))],"%b-%y")
axis(1,at=c(min(WTI$Date),WTI$Date[WTI_samples],max(WTI$Date)),labels=Date_formatted,las = 2)

plot(Brent$Date,Brent$Price,type="l",col="blue",main="Structural Breaks for Brent Oil",ylab="Price/Barrel (USD)",xlab="",xaxt="n")
abline(v=Brent$Date[Brent_samples[1]],col="brown")
abline(v=Brent$Date[Brent_samples[2]],col="brown")
abline(v=Brent$Date[Brent_samples[3]],col="brown")
abline(v=Brent$Date[Brent_samples[4]],col="brown")
abline(v=Brent$Date[Brent_samples[5]],col="brown")
Date_formatted <- format(Brent$Date[c(1,Brent_samples,nrow(Brent))],"%b-%y")
axis(1,at=c(min(Brent$Date),Brent$Date[Brent_samples],max(Brent$Date)),labels=Date_formatted,las = 2)

### Plotting the Density and Histogram for Returns

par(mfrow=c(1,2))
densityPlot(as.timeSeries(WTI_xts),hist = TRUE,title = FALSE,grid=FALSE)
title(main="Density Plot for WTI Oil",xlab="Returns",ylab="Density/Frequency",col.main="blue",cex.main=0.7)
legend(x=0.06,y=25,legend = c("Normal","Empirical"),bty="n",cex=0.7,text.col = c("gray","brown"))

densityPlot(as.timeSeries(WTI_xts),hist = TRUE,title = FALSE,grid=FALSE)
title(main="Density Plot for Brent Oil",xlab="Returns",ylab="Density/Frequency",col.main="blue",cex.main=0.7)
legend(x=0.06,y=25,legend = c("Normal","Empirical"),bty="n",cex=0.7,text.col = c("gray","brown"),col=c("gray","brown"))

#### Base Model without Volume ######

### Garch Models

WTI_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model1 <- ugarchfit(WTI_Spec1,WTI_xts,out.sample = 1000)
WTI_Model2 <- ugarchfit(WTI_Spec2,WTI_xts,out.sample = 1000)
WTI_Model3 <- ugarchfit(WTI_Spec3,WTI_xts,out.sample = 1000)
WTI_Model4 <- ugarchfit(WTI_Spec4,WTI_xts,out.sample = 1000)

Brent_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model1 <- ugarchfit(Brent_Spec1,Brent_xts,out.sample = 1000)
Brent_Model2 <- ugarchfit(Brent_Spec2,Brent_xts,out.sample = 1000)
Brent_Model3 <- ugarchfit(Brent_Spec3,Brent_xts,out.sample = 1000)
Brent_Model4 <- ugarchfit(Brent_Spec4,Brent_xts,out.sample = 1000)

#### egarch models

WTI_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model5 <- ugarchfit(WTI_Spec5,WTI_xts,out.sample = 1000)
WTI_Model6 <- ugarchfit(WTI_Spec6,WTI_xts,out.sample = 1000)
WTI_Model7 <- ugarchfit(WTI_Spec7,WTI_xts,out.sample = 1000)
WTI_Model8 <- ugarchfit(WTI_Spec8,WTI_xts,out.sample = 1000)

Brent_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model5 <- ugarchfit(Brent_Spec5,Brent_xts,out.sample = 1000)
Brent_Model6 <- ugarchfit(Brent_Spec6,Brent_xts,out.sample = 1000)
Brent_Model7 <- ugarchfit(Brent_Spec7,Brent_xts,out.sample = 1000)
Brent_Model8 <- ugarchfit(Brent_Spec8,Brent_xts,out.sample = 1000)

### gjrgarch models

WTI_Spec9   <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec10  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec11  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec12  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model9  <- ugarchfit(WTI_Spec9,WTI_xts,out.sample = 1000)
WTI_Model10 <- ugarchfit(WTI_Spec10,WTI_xts,out.sample = 1000)
WTI_Model11 <- ugarchfit(WTI_Spec11,WTI_xts,out.sample = 1000)
WTI_Model12 <- ugarchfit(WTI_Spec12,WTI_xts,out.sample = 1000)

Brent_Spec9  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec10 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec11 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec12 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:5]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

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

WTI_Roll5 <-  ugarchroll(WTI_Spec5,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll6 <-  ugarchroll(WTI_Spec6,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll7 <-  ugarchroll(WTI_Spec7,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll8 <-  ugarchroll(WTI_Spec8,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)

Brent_Roll5 <-  ugarchroll(Brent_Spec5,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll6 <-  ugarchroll(Brent_Spec6,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll7 <-  ugarchroll(Brent_Spec7,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll8 <-  ugarchroll(Brent_Spec8,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)

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

# ### JSU specification

WTI_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

WTI_Model_Sample_JSU1 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub1,out.sample = 250)
WTI_Model_Sample_JSU2 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub2,out.sample = 250)
WTI_Model_Sample_JSU3 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub3,out.sample = 250)
WTI_Model_Sample_JSU4 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub4,out.sample = 250)
WTI_Model_Sample_JSU5 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub5,out.sample = 250)
WTI_Model_Sample_JSU6 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

Brent_Model_Sample_JSU1 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub1,out.sample = 250)
Brent_Model_Sample_JSU2 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub2,out.sample = 250)
Brent_Model_Sample_JSU3 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub3,out.sample = 250)
Brent_Model_Sample_JSU4 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub4,out.sample = 250)
Brent_Model_Sample_JSU5 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub5,out.sample = 250)
Brent_Model_Sample_JSU6 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_JSU_Sample1 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample2 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample3 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample4 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample5 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample6 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_JSU_Sample1 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample2 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample3 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample4 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample5 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample6 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### SGED specification

WTI_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

WTI_Model_Sample_SGED1 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub1,out.sample = 250)
WTI_Model_Sample_SGED2 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub2,out.sample = 250)
WTI_Model_Sample_SGED3 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub3,out.sample = 250)
WTI_Model_Sample_SGED4 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub4,out.sample = 250)
WTI_Model_Sample_SGED5 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub5,out.sample = 250)
WTI_Model_Sample_SGED6 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

Brent_Model_Sample_SGED1 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub1,out.sample = 250)
Brent_Model_Sample_SGED2 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub2,out.sample = 250)
Brent_Model_Sample_SGED3 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub3,out.sample = 250)
Brent_Model_Sample_SGED4 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub4,out.sample = 250)
Brent_Model_Sample_SGED5 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub5,out.sample = 250)
Brent_Model_Sample_SGED6 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_SGED_Sample1 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample2 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample3 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample4 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample5 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample6 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_SGED_Sample1 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample2 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample3 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample4 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample5 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample6 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### SST Specification
# 
WTI_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
# 
WTI_Model_Sample_SST1 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub1,out.sample = 250)
WTI_Model_Sample_SST2 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub2,out.sample = 250)
WTI_Model_Sample_SST3 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub3,out.sample = 250)
WTI_Model_Sample_SST4 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub4,out.sample = 250)
WTI_Model_Sample_SST5 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub5,out.sample = 250)
WTI_Model_Sample_SST6 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")

Brent_Model_Sample_SST1 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub1,out.sample = 250)
Brent_Model_Sample_SST2 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub2,out.sample = 250)
Brent_Model_Sample_SST3 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub3,out.sample = 250)
Brent_Model_Sample_SST4 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub4,out.sample = 250)
Brent_Model_Sample_SST5 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub5,out.sample = 250)
Brent_Model_Sample_SST6 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_SST_Sample1 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample2 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample3 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample4 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample5 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample6 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_SST_Sample1 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample2 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample3 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample4 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample5 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample6 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### Normal Specification

WTI_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model_Sample_Normal1 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)
WTI_Model_Sample_Normal2 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)
WTI_Model_Sample_Normal3 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)
WTI_Model_Sample_Normal4 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)
WTI_Model_Sample_Normal5 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)
WTI_Model_Sample_Normal6 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model_Sample_Normal1 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)
Brent_Model_Sample_Normal2 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)
Brent_Model_Sample_Normal3 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)
Brent_Model_Sample_Normal4 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)
Brent_Model_Sample_Normal5 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)
Brent_Model_Sample_Normal6 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)

WTI_Roll_Normal_Sample1 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample2 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample3 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample4 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample5 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample6 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_Normal_Sample1 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample2 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample3 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample4 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample5 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample6 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### Sample Residuals Statistics

WTI_Res_Sample1 <- c(normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample2 <- c(normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample3 <- c(normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample4 <- c(normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample5 <- c(normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample6 <- c(normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="ks")@test$statistic)

Brent_Res_Sample1 <- c(normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample2 <- c(normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample3 <- c(normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample4 <- c(normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample5 <- c(normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample6 <- c(normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="ks")@test$statistic)

WTI_Res_Sample_Test <- rbind(WTI_Res_Sample1,WTI_Res_Sample2,WTI_Res_Sample3,WTI_Res_Sample4,WTI_Res_Sample5,WTI_Res_Sample6)
Brent_Res_Sample_Test   <- rbind(Brent_Res_Sample1,Brent_Res_Sample2,Brent_Res_Sample3,Brent_Res_Sample4,Brent_Res_Sample5,Brent_Res_Sample6)
rownames(WTI_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")
colnames(WTI_Res_Sample_Test) <- c("JB","SW","KS") 
rownames(Brent_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")
colnames(Brent_Res_Sample_Test) <- c("JB","SW","KS") 

### Pearson Specification

WTI_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample6 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$sigma)

WTI_Param_Pearson_Subsample <- rbind(WTI_Param_Pearson_Sample1,WTI_Param_Pearson_Sample2,WTI_Param_Pearson_Sample3,WTI_Param_Pearson_Sample4,WTI_Param_Pearson_Sample5,WTI_Param_Pearson_Sample6)
rownames(WTI_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

Brent_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample6 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$sigma)

Brent_Param_Pearson_Subsample <- rbind(Brent_Param_Pearson_Sample1,Brent_Param_Pearson_Sample2,Brent_Param_Pearson_Sample3,Brent_Param_Pearson_Sample4,Brent_Param_Pearson_Sample5,Brent_Param_Pearson_Sample6)
rownames(Brent_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

### johnson SU parameter


WTI_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample6 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$sigma,moment = "find")

WTI_Param_Johnson_Subsample <- rbind(WTI_Param_Johnson_Sample1,WTI_Param_Johnson_Sample2,WTI_Param_Johnson_Sample3,WTI_Param_Johnson_Sample4,WTI_Param_Johnson_Sample5,WTI_Param_Johnson_Sample6)
rownames(WTI_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

Brent_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample6 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$sigma,moment = "find")

Brent_Param_Johnson_Subsample <- rbind(Brent_Param_Johnson_Sample1,Brent_Param_Johnson_Sample2,Brent_Param_Johnson_Sample3,Brent_Param_Johnson_Sample4,Brent_Param_Johnson_Sample5,Brent_Param_Johnson_Sample6)
rownames(Brent_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

### Plotting the VaR for subsamples

WTI_Sample1_Dates <- format(WTI$Date[(WTI_samples[1]-249):WTI_samples[1]],"%b-%y")
WTI_Sample2_Dates <- format(WTI$Date[(WTI_samples[2]-249):WTI_samples[2]],"%b-%y")
WTI_Sample3_Dates <- format(WTI$Date[(WTI_samples[3]-249):WTI_samples[3]],"%b-%y")
WTI_Sample4_Dates <- format(WTI$Date[(WTI_samples[4]-249):WTI_samples[4]],"%b-%y")
WTI_Sample5_Dates <- format(WTI$Date[(WTI_samples[5]-249):WTI_samples[5]],"%b-%y")
WTI_Sample6_Dates <- format(WTI$Date[(length(WTI_xts)-249):length(WTI_xts)],"%b-%y")

Brent_Sample1_Dates <- format(Brent$Date[(Brent_samples[1]-249):Brent_samples[1]],"%b-%y")
Brent_Sample2_Dates <- format(Brent$Date[(Brent_samples[2]-249):Brent_samples[2]],"%b-%y")
Brent_Sample3_Dates <- format(Brent$Date[(Brent_samples[3]-249):Brent_samples[3]],"%b-%y")
Brent_Sample4_Dates <- format(Brent$Date[(Brent_samples[4]-249):Brent_samples[4]],"%b-%y")
Brent_Sample5_Dates <- format(Brent$Date[(Brent_samples[5]-249):Brent_samples[5]],"%b-%y")
Brent_Sample6_Dates <- format(Brent$Date[(length(Brent_xts)-249):length(Brent_xts)],"%b-%y")

### Forecast Performance Measures ###

### Loss functions ###

# In Sample:

Brent_Temp_Returns <- Returns_Brent[1:4646]
Brent_Realized_Vol <- array(data=0,dim=4500)

for(i in 1:4500)
{
  Brent_Realized_Vol[i] <- sum(Brent_Temp_Returns[i:(i+146)]^2)
  #+Brent_Temp_Returns[i+1]^2+Brent_Temp_Returns[i+2]^2+Brent_Temp_Returns[i+3]^2+
  #Brent_Temp_Returns[i+4]^2+Brent_Temp_Returns[i+5]^2+Brent_Temp_Returns[i+6]^2+Brent_Temp_Returns[i+7]^2+
  #Brent_Temp_Returns[i+8]^2+Brent_Temp_Returns[i+9]^2
}

Brent_Realized_Vol  <- sqrt(Brent_Realized_Vol)  

Brent_In_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(Brent_In_Sample) <- c("Normal","JSU","SGED","SST")
rownames(Brent_In_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
Brent_In_Sample[1,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1")))
Brent_In_Sample[2,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE2")))
Brent_In_Sample[3,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="AE1")))
Brent_In_Sample[4,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="AE2")))
Brent_In_Sample[5,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="R2LOG")))
Brent_In_Sample[6,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="QLIKE")))
Brent_In_Sample <- round(Brent_In_Sample,6)


### SPA test

SPA_Brent_base <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1"))
                               ,0.05,statistic = "Tmax")

Brent_Base <- Brent_Model5@fit$sigma[147:4646]

# FOR WTI

WTI_Temp_Returns <- Returns_WTI[1:4537]
WTI_Realized_Vol <- array(data=0,dim=4400)

for(i in 1:4400)
{
  WTI_Realized_Vol[i] <- sum(WTI_Temp_Returns[i:(i+137)]^2)
  #+WTI_Temp_Returns[i+1]^2+WTI_Temp_Returns[i+2]^2+WTI_Temp_Returns[i+3]^2+
  #WTI_Temp_Returns[i+4]^2+WTI_Temp_Returns[i+5]^2+WTI_Temp_Returns[i+6]^2+WTI_Temp_Returns[i+7]^2+
  #WTI_Temp_Returns[i+8]^2+WTI_Temp_Returns[i+9]^2
}

WTI_Realized_Vol  <- sqrt(WTI_Realized_Vol)  

WTI_In_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(WTI_In_Sample) <- c("Normal","JSU","SGED","SST")
rownames(WTI_In_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
WTI_In_Sample[1,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1")))
WTI_In_Sample[2,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE2")))
WTI_In_Sample[3,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="AE1")))
WTI_In_Sample[4,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="AE2")))
WTI_In_Sample[5,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="R2LOG")))
WTI_In_Sample[6,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="QLIKE")))
WTI_In_Sample <- round(WTI_In_Sample,6)


### SPA test

SPA_WTI_base <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1"))
                             ,0.05,statistic = "Tmax")

WTI_Base <- WTI_Model7@fit$sigma[138:4537]

# Out Sample

## For Brent
BrentOut_Temp_Returns <- Returns_Brent[(length(Returns_Brent)-999):length(Returns_Brent)]
Brent_Out_Realized_Vol <- array(data=0,dim=900)

for(i in 1:900)
{
  Brent_Out_Realized_Vol[i] <- sum(BrentOut_Temp_Returns[i:(i+99)]^2)
  #+BrentOut_Temp_Returns[i+1]^2+BrentOut_Temp_Returns[i+2]^2+BrentOut_Temp_Returns[i+3]^2+
  #BrentOut_Temp_Returns[i+4]^2+BrentOut_Temp_Returns[i+5]^2+BrentOut_Temp_Returns[i+6]^2+BrentOut_Temp_Returns[i+7]^2+
  #BrentOut_Temp_Returns[i+8]^2+BrentOut_Temp_Returns[i+9]^2
}

Brent_Out_Realized_Vol  <- sqrt(Brent_Out_Realized_Vol)  

Brent_Out_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(Brent_Out_Sample) <- c("Normal","JSU","SGED","SST")
rownames(Brent_Out_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
Brent_Out_Sample[1,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE1")))
Brent_Out_Sample[2,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE2")))
Brent_Out_Sample[3,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="AE1")))
Brent_Out_Sample[4,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="AE2")))
Brent_Out_Sample[5,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="R2LOG")))
Brent_Out_Sample[6,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="QLIKE")))
Brent_Out_Sample <- round(Brent_Out_Sample,6)


### SPA test

Brent_SPA_out <- MCSprocedure(cbind(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE1"))
                            ,0.05,statistic = "Tmax")

# FOR WTI

WTIOut_Temp_Returns <- Returns_WTI[(length(Returns_WTI)-999):length(Returns_WTI)]
WTI_Out_Realized_Vol <- array(data=0,dim=900)

for(i in 1:900)
{
  WTI_Out_Realized_Vol[i] <- sum(WTIOut_Temp_Returns[i:(i+99)]^2)
  #+WTIOut_Temp_Returns[i+1]^2+WTIOut_Temp_Returns[i+2]^2+WTIOut_Temp_Returns[i+3]^2+
  #WTIOut_Temp_Returns[i+4]^2+WTIOut_Temp_Returns[i+5]^2+WTIOut_Temp_Returns[i+6]^2+WTIOut_Temp_Returns[i+7]^2+
  #WTIOut_Temp_Returns[i+8]^2+WTIOut_Temp_Returns[i+9]^2
}

WTI_Out_Realized_Vol  <- sqrt(WTI_Out_Realized_Vol)  

WTI_Out_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(WTI_Out_Sample) <- c("Normal","JSU","SGED","SST")
rownames(WTI_Out_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
WTI_Out_Sample[1,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE1")))
WTI_Out_Sample[2,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE2")))
WTI_Out_Sample[3,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="AE1")))
WTI_Out_Sample[4,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="AE2")))
WTI_Out_Sample[5,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="R2LOG")))
WTI_Out_Sample[6,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="QLIKE")))
WTI_Out_Sample <- round(WTI_Out_Sample,6)


### SPA test

WTI_SPA_out <- MCSprocedure(cbind(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE1"),
                               LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE1"),
                               LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE1"),
                               LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE1"))
                         ,0.05,statistic = "Tmax")


#### Model with Lagged Volume ##### 

### Repeating the above codes ##

### Garch Models

WTI_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model1 <- ugarchfit(WTI_Spec1,WTI_xts,out.sample = 1000)
WTI_Model2 <- ugarchfit(WTI_Spec2,WTI_xts,out.sample = 1000)
WTI_Model3 <- ugarchfit(WTI_Spec3,WTI_xts,out.sample = 1000)
WTI_Model4 <- ugarchfit(WTI_Spec4,WTI_xts,out.sample = 1000)

Brent_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model1 <- ugarchfit(Brent_Spec1,Brent_xts,out.sample = 1000)
Brent_Model2 <- ugarchfit(Brent_Spec2,Brent_xts,out.sample = 1000)
Brent_Model3 <- ugarchfit(Brent_Spec3,Brent_xts,out.sample = 1000)
Brent_Model4 <- ugarchfit(Brent_Spec4,Brent_xts,out.sample = 1000)

#### egarch models

WTI_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model5 <- ugarchfit(WTI_Spec5,WTI_xts,out.sample = 1000)
WTI_Model6 <- ugarchfit(WTI_Spec6,WTI_xts,out.sample = 1000)
WTI_Model7 <- ugarchfit(WTI_Spec7,WTI_xts,out.sample = 1000)
WTI_Model8 <- ugarchfit(WTI_Spec8,WTI_xts,out.sample = 1000)

Brent_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model5 <- ugarchfit(Brent_Spec5,Brent_xts,out.sample = 1000)
Brent_Model6 <- ugarchfit(Brent_Spec6,Brent_xts,out.sample = 1000)
Brent_Model7 <- ugarchfit(Brent_Spec7,Brent_xts,out.sample = 1000)
Brent_Model8 <- ugarchfit(Brent_Spec8,Brent_xts,out.sample = 1000)

### gjrgarch models

WTI_Spec9   <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec10  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec11  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec12  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model9  <- ugarchfit(WTI_Spec9,WTI_xts,out.sample = 1000)
WTI_Model10 <- ugarchfit(WTI_Spec10,WTI_xts,out.sample = 1000)
WTI_Model11 <- ugarchfit(WTI_Spec11,WTI_xts,out.sample = 1000)
WTI_Model12 <- ugarchfit(WTI_Spec12,WTI_xts,out.sample = 1000)

Brent_Spec9  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec10 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec11 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec12 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

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

WTI_Roll5 <-  ugarchroll(WTI_Spec5,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll6 <-  ugarchroll(WTI_Spec6,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll7 <-  ugarchroll(WTI_Spec7,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll8 <-  ugarchroll(WTI_Spec8,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)

Brent_Roll5 <-  ugarchroll(Brent_Spec5,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll6 <-  ugarchroll(Brent_Spec6,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll7 <-  ugarchroll(Brent_Spec7,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll8 <-  ugarchroll(Brent_Spec8,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)

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

# ### JSU specification

WTI_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

WTI_Model_Sample_JSU1 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub1,out.sample = 250)
WTI_Model_Sample_JSU2 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub2,out.sample = 250)
WTI_Model_Sample_JSU3 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub3,out.sample = 250)
WTI_Model_Sample_JSU4 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub4,out.sample = 250)
WTI_Model_Sample_JSU5 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub5,out.sample = 250)
WTI_Model_Sample_JSU6 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

Brent_Model_Sample_JSU1 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub1,out.sample = 250)
Brent_Model_Sample_JSU2 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub2,out.sample = 250)
Brent_Model_Sample_JSU3 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub3,out.sample = 250)
Brent_Model_Sample_JSU4 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub4,out.sample = 250)
Brent_Model_Sample_JSU5 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub5,out.sample = 250)
Brent_Model_Sample_JSU6 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_JSU_Sample1 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample2 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample3 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample4 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample5 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample6 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_JSU_Sample1 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample2 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample3 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample4 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample5 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample6 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### SGED specification

WTI_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

WTI_Model_Sample_SGED1 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub1,out.sample = 250)
WTI_Model_Sample_SGED2 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub2,out.sample = 250)
WTI_Model_Sample_SGED3 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub3,out.sample = 250)
WTI_Model_Sample_SGED4 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub4,out.sample = 250)
WTI_Model_Sample_SGED5 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub5,out.sample = 250)
WTI_Model_Sample_SGED6 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

Brent_Model_Sample_SGED1 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub1,out.sample = 250)
Brent_Model_Sample_SGED2 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub2,out.sample = 250)
Brent_Model_Sample_SGED3 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub3,out.sample = 250)
Brent_Model_Sample_SGED4 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub4,out.sample = 250)
Brent_Model_Sample_SGED5 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub5,out.sample = 250)
Brent_Model_Sample_SGED6 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_SGED_Sample1 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample2 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample3 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample4 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample5 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample6 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_SGED_Sample1 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample2 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample3 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample4 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample5 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample6 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### SST Specification
# 
WTI_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
# 
WTI_Model_Sample_SST1 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub1,out.sample = 250)
WTI_Model_Sample_SST2 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub2,out.sample = 250)
WTI_Model_Sample_SST3 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub3,out.sample = 250)
WTI_Model_Sample_SST4 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub4,out.sample = 250)
WTI_Model_Sample_SST5 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub5,out.sample = 250)
WTI_Model_Sample_SST6 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")

Brent_Model_Sample_SST1 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub1,out.sample = 250)
Brent_Model_Sample_SST2 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub2,out.sample = 250)
Brent_Model_Sample_SST3 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub3,out.sample = 250)
Brent_Model_Sample_SST4 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub4,out.sample = 250)
Brent_Model_Sample_SST5 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub5,out.sample = 250)
Brent_Model_Sample_SST6 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_SST_Sample1 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample2 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample3 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample4 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample5 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample6 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_SST_Sample1 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample2 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample3 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample4 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample5 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample6 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### Normal Specification

WTI_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model_Sample_Normal1 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)
WTI_Model_Sample_Normal2 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)
WTI_Model_Sample_Normal3 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)
WTI_Model_Sample_Normal4 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)
WTI_Model_Sample_Normal5 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)
WTI_Model_Sample_Normal6 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model_Sample_Normal1 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)
Brent_Model_Sample_Normal2 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)
Brent_Model_Sample_Normal3 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)
Brent_Model_Sample_Normal4 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)
Brent_Model_Sample_Normal5 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)
Brent_Model_Sample_Normal6 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)

WTI_Roll_Normal_Sample1 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample2 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample3 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample4 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample5 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample6 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_Normal_Sample1 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample2 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample3 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample4 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample5 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample6 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### Sample Residuals Statistics

WTI_Res_Sample1 <- c(normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample2 <- c(normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample3 <- c(normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample4 <- c(normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample5 <- c(normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample6 <- c(normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="ks")@test$statistic)

Brent_Res_Sample1 <- c(normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample2 <- c(normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample3 <- c(normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample4 <- c(normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample5 <- c(normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample6 <- c(normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="ks")@test$statistic)

WTI_Res_Sample_Test <- rbind(WTI_Res_Sample1,WTI_Res_Sample2,WTI_Res_Sample3,WTI_Res_Sample4,WTI_Res_Sample5,WTI_Res_Sample6)
Brent_Res_Sample_Test   <- rbind(Brent_Res_Sample1,Brent_Res_Sample2,Brent_Res_Sample3,Brent_Res_Sample4,Brent_Res_Sample5,Brent_Res_Sample6)
rownames(WTI_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")
colnames(WTI_Res_Sample_Test) <- c("JB","SW","KS") 
rownames(Brent_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")
colnames(Brent_Res_Sample_Test) <- c("JB","SW","KS") 

### Pearson Specification

WTI_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample6 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$sigma)

WTI_Param_Pearson_Subsample <- rbind(WTI_Param_Pearson_Sample1,WTI_Param_Pearson_Sample2,WTI_Param_Pearson_Sample3,WTI_Param_Pearson_Sample4,WTI_Param_Pearson_Sample5,WTI_Param_Pearson_Sample6)
rownames(WTI_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

Brent_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample6 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$sigma)

Brent_Param_Pearson_Subsample <- rbind(Brent_Param_Pearson_Sample1,Brent_Param_Pearson_Sample2,Brent_Param_Pearson_Sample3,Brent_Param_Pearson_Sample4,Brent_Param_Pearson_Sample5,Brent_Param_Pearson_Sample6)
rownames(Brent_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

### johnson SU parameter


WTI_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample6 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$sigma,moment = "find")

WTI_Param_Johnson_Subsample <- rbind(WTI_Param_Johnson_Sample1,WTI_Param_Johnson_Sample2,WTI_Param_Johnson_Sample3,WTI_Param_Johnson_Sample4,WTI_Param_Johnson_Sample5,WTI_Param_Johnson_Sample6)
rownames(WTI_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

Brent_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample6 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$sigma,moment = "find")

Brent_Param_Johnson_Subsample <- rbind(Brent_Param_Johnson_Sample1,Brent_Param_Johnson_Sample2,Brent_Param_Johnson_Sample3,Brent_Param_Johnson_Sample4,Brent_Param_Johnson_Sample5,Brent_Param_Johnson_Sample6)
rownames(Brent_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

### Plotting the VaR for subsamples

WTI_Sample1_Dates <- format(WTI$Date[(WTI_samples[1]-249):WTI_samples[1]],"%b-%y")
WTI_Sample2_Dates <- format(WTI$Date[(WTI_samples[2]-249):WTI_samples[2]],"%b-%y")
WTI_Sample3_Dates <- format(WTI$Date[(WTI_samples[3]-249):WTI_samples[3]],"%b-%y")
WTI_Sample4_Dates <- format(WTI$Date[(WTI_samples[4]-249):WTI_samples[4]],"%b-%y")
WTI_Sample5_Dates <- format(WTI$Date[(WTI_samples[5]-249):WTI_samples[5]],"%b-%y")
WTI_Sample6_Dates <- format(WTI$Date[(length(WTI_xts)-249):length(WTI_xts)],"%b-%y")

Brent_Sample1_Dates <- format(Brent$Date[(Brent_samples[1]-249):Brent_samples[1]],"%b-%y")
Brent_Sample2_Dates <- format(Brent$Date[(Brent_samples[2]-249):Brent_samples[2]],"%b-%y")
Brent_Sample3_Dates <- format(Brent$Date[(Brent_samples[3]-249):Brent_samples[3]],"%b-%y")
Brent_Sample4_Dates <- format(Brent$Date[(Brent_samples[4]-249):Brent_samples[4]],"%b-%y")
Brent_Sample5_Dates <- format(Brent$Date[(Brent_samples[5]-249):Brent_samples[5]],"%b-%y")
Brent_Sample6_Dates <- format(Brent$Date[(length(Brent_xts)-249):length(Brent_xts)],"%b-%y")

### Forecast Performance Measures ###

### Loss functions ###

# In Sample:

Brent_Temp_Returns <- Returns_Brent[1:4646]
Brent_Realized_Vol <- array(data=0,dim=4500)

for(i in 1:4500)
{
  Brent_Realized_Vol[i] <- sum(Brent_Temp_Returns[i:(i+146)]^2)
  #+Brent_Temp_Returns[i+1]^2+Brent_Temp_Returns[i+2]^2+Brent_Temp_Returns[i+3]^2+
  #Brent_Temp_Returns[i+4]^2+Brent_Temp_Returns[i+5]^2+Brent_Temp_Returns[i+6]^2+Brent_Temp_Returns[i+7]^2+
  #Brent_Temp_Returns[i+8]^2+Brent_Temp_Returns[i+9]^2
}

Brent_Realized_Vol  <- sqrt(Brent_Realized_Vol)  

Brent_In_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(Brent_In_Sample) <- c("Normal","JSU","SGED","SST")
rownames(Brent_In_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
Brent_In_Sample[1,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1")))
Brent_In_Sample[2,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE2")))
Brent_In_Sample[3,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="AE1")))
Brent_In_Sample[4,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="AE2")))
Brent_In_Sample[5,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="R2LOG")))
Brent_In_Sample[6,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="QLIKE")))
Brent_In_Sample <- round(Brent_In_Sample,6)


### SPA test

SPA_Brent_base <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1"))
                               ,0.05,statistic = "Tmax")

Brent_Lag <- Brent_Model5@fit$sigma[147:4646]

# FOR WTI

WTI_Temp_Returns <- Returns_WTI[1:4537]
WTI_Realized_Vol <- array(data=0,dim=4400)

for(i in 1:4400)
{
  WTI_Realized_Vol[i] <- sum(WTI_Temp_Returns[i:(i+137)]^2)
  #+WTI_Temp_Returns[i+1]^2+WTI_Temp_Returns[i+2]^2+WTI_Temp_Returns[i+3]^2+
  #WTI_Temp_Returns[i+4]^2+WTI_Temp_Returns[i+5]^2+WTI_Temp_Returns[i+6]^2+WTI_Temp_Returns[i+7]^2+
  #WTI_Temp_Returns[i+8]^2+WTI_Temp_Returns[i+9]^2
}

WTI_Realized_Vol  <- sqrt(WTI_Realized_Vol)  

WTI_In_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(WTI_In_Sample) <- c("Normal","JSU","SGED","SST")
rownames(WTI_In_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
WTI_In_Sample[1,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1")))
WTI_In_Sample[2,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE2")))
WTI_In_Sample[3,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="AE1")))
WTI_In_Sample[4,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="AE2")))
WTI_In_Sample[5,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="R2LOG")))
WTI_In_Sample[6,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="QLIKE")))
WTI_In_Sample <- round(WTI_In_Sample,6)


### SPA test

SPA_WTI_base <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1"))
                             ,0.05,statistic = "Tmax")

WTI_Lag <- WTI_Model7@fit$sigma[138:4537]

# Out Sample

## For Brent
BrentOut_Temp_Returns <- Returns_Brent[(length(Returns_Brent)-999):length(Returns_Brent)]
Brent_Out_Realized_Vol <- array(data=0,dim=900)

for(i in 1:900)
{
  Brent_Out_Realized_Vol[i] <- sum(BrentOut_Temp_Returns[i:(i+99)]^2)
  #+BrentOut_Temp_Returns[i+1]^2+BrentOut_Temp_Returns[i+2]^2+BrentOut_Temp_Returns[i+3]^2+
  #BrentOut_Temp_Returns[i+4]^2+BrentOut_Temp_Returns[i+5]^2+BrentOut_Temp_Returns[i+6]^2+BrentOut_Temp_Returns[i+7]^2+
  #BrentOut_Temp_Returns[i+8]^2+BrentOut_Temp_Returns[i+9]^2
}

Brent_Out_Realized_Vol  <- sqrt(Brent_Out_Realized_Vol)  

Brent_Out_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(Brent_Out_Sample) <- c("Normal","JSU","SGED","SST")
rownames(Brent_Out_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
Brent_Out_Sample[1,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE1")))
Brent_Out_Sample[2,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE2")))
Brent_Out_Sample[3,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="AE1")))
Brent_Out_Sample[4,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="AE2")))
Brent_Out_Sample[5,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="R2LOG")))
Brent_Out_Sample[6,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="QLIKE")))
Brent_Out_Sample <- round(Brent_Out_Sample,6)


### SPA test

Brent_SPA_out <- MCSprocedure(cbind(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE1"),
                                    LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE1"),
                                    LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE1"),
                                    LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE1"))
                              ,0.05,statistic = "Tmax")

# FOR WTI

WTIOut_Temp_Returns <- Returns_WTI[(length(Returns_WTI)-999):length(Returns_WTI)]
WTI_Out_Realized_Vol <- array(data=0,dim=900)

for(i in 1:900)
{
  WTI_Out_Realized_Vol[i] <- sum(WTIOut_Temp_Returns[i:(i+99)]^2)
  #+WTIOut_Temp_Returns[i+1]^2+WTIOut_Temp_Returns[i+2]^2+WTIOut_Temp_Returns[i+3]^2+
  #WTIOut_Temp_Returns[i+4]^2+WTIOut_Temp_Returns[i+5]^2+WTIOut_Temp_Returns[i+6]^2+WTIOut_Temp_Returns[i+7]^2+
  #WTIOut_Temp_Returns[i+8]^2+WTIOut_Temp_Returns[i+9]^2
}

WTI_Out_Realized_Vol  <- sqrt(WTI_Out_Realized_Vol)  

WTI_Out_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(WTI_Out_Sample) <- c("Normal","JSU","SGED","SST")
rownames(WTI_Out_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
WTI_Out_Sample[1,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE1")))
WTI_Out_Sample[2,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE2")))
WTI_Out_Sample[3,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="AE1")))
WTI_Out_Sample[4,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="AE2")))
WTI_Out_Sample[5,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="R2LOG")))
WTI_Out_Sample[6,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="QLIKE")))
WTI_Out_Sample <- round(WTI_Out_Sample,6)


### SPA test

WTI_SPA_out <- MCSprocedure(cbind(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE1"))
                            ,0.05,statistic = "Tmax")


##### Model with Contemporaneous Volume #####


#### Structural Breaks for WTI and Brent #####

WTI_Dummy <- cbind(WTI_Dummy1,WTI_Dummy2,WTI_Dummy3,WTI_Dummy4,WTI_Dummy5,log(as.integer(WTI$Volume[2:5538])))
Brent_Dummy <- cbind(Brent_Dummy1,Brent_Dummy2,Brent_Dummy3,Brent_Dummy4,Brent_Dummy5,log(Brent$Volume[2:5647]))

### Repeating the code for this model now

### Garch Models

WTI_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model1 <- ugarchfit(WTI_Spec1,WTI_xts,out.sample = 1000)
WTI_Model2 <- ugarchfit(WTI_Spec2,WTI_xts,out.sample = 1000)
WTI_Model3 <- ugarchfit(WTI_Spec3,WTI_xts,out.sample = 1000)
WTI_Model4 <- ugarchfit(WTI_Spec4,WTI_xts,out.sample = 1000)

Brent_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model1 <- ugarchfit(Brent_Spec1,Brent_xts,out.sample = 1000)
Brent_Model2 <- ugarchfit(Brent_Spec2,Brent_xts,out.sample = 1000)
Brent_Model3 <- ugarchfit(Brent_Spec3,Brent_xts,out.sample = 1000)
Brent_Model4 <- ugarchfit(Brent_Spec4,Brent_xts,out.sample = 1000)

#### egarch models

WTI_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model5 <- ugarchfit(WTI_Spec5,WTI_xts,out.sample = 1000)
WTI_Model6 <- ugarchfit(WTI_Spec6,WTI_xts,out.sample = 1000)
WTI_Model7 <- ugarchfit(WTI_Spec7,WTI_xts,out.sample = 1000)
WTI_Model8 <- ugarchfit(WTI_Spec8,WTI_xts,out.sample = 1000)

Brent_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model5 <- ugarchfit(Brent_Spec5,Brent_xts,out.sample = 1000)
Brent_Model6 <- ugarchfit(Brent_Spec6,Brent_xts,out.sample = 1000)
Brent_Model7 <- ugarchfit(Brent_Spec7,Brent_xts,out.sample = 1000)
Brent_Model8 <- ugarchfit(Brent_Spec8,Brent_xts,out.sample = 1000)

### gjrgarch models

WTI_Spec9   <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
WTI_Spec10  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
WTI_Spec11  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
WTI_Spec12  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=WTI_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model9  <- ugarchfit(WTI_Spec9,WTI_xts,out.sample = 1000)
WTI_Model10 <- ugarchfit(WTI_Spec10,WTI_xts,out.sample = 1000)
WTI_Model11 <- ugarchfit(WTI_Spec11,WTI_xts,out.sample = 1000)
WTI_Model12 <- ugarchfit(WTI_Spec12,WTI_xts,out.sample = 1000)

Brent_Spec9  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Brent_Spec10 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Brent_Spec11 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Brent_Spec12 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Brent_Dummy[,1:6]),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

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

WTI_Roll5 <-  ugarchroll(WTI_Spec5,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll6 <-  ugarchroll(WTI_Spec6,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll7 <-  ugarchroll(WTI_Spec7,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
WTI_Roll8 <-  ugarchroll(WTI_Spec8,WTI_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)

Brent_Roll5 <-  ugarchroll(Brent_Spec5,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll6 <-  ugarchroll(Brent_Spec6,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll7 <-  ugarchroll(Brent_Spec7,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)
Brent_Roll8 <-  ugarchroll(Brent_Spec8,Brent_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",keep.coef = TRUE)

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

# ### JSU specification

WTI_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

WTI_Model_Sample_JSU1 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub1,out.sample = 250)
WTI_Model_Sample_JSU2 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub2,out.sample = 250)
WTI_Model_Sample_JSU3 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub3,out.sample = 250)
WTI_Model_Sample_JSU4 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub4,out.sample = 250)
WTI_Model_Sample_JSU5 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub5,out.sample = 250)
WTI_Model_Sample_JSU6 <- ugarchfit(WTI_Spec_Sample_JSU,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

Brent_Model_Sample_JSU1 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub1,out.sample = 250)
Brent_Model_Sample_JSU2 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub2,out.sample = 250)
Brent_Model_Sample_JSU3 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub3,out.sample = 250)
Brent_Model_Sample_JSU4 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub4,out.sample = 250)
Brent_Model_Sample_JSU5 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub5,out.sample = 250)
Brent_Model_Sample_JSU6 <- ugarchfit(Brent_Spec_Sample_JSU,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_JSU_Sample1 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample2 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample3 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample4 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample5 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_JSU_Sample6 <-  ugarchroll(WTI_Spec_Sample_JSU,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_JSU_Sample1 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample2 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample3 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample4 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample5 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_JSU_Sample6 <-  ugarchroll(Brent_Spec_Sample_JSU,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### SGED specification

WTI_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

WTI_Model_Sample_SGED1 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub1,out.sample = 250)
WTI_Model_Sample_SGED2 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub2,out.sample = 250)
WTI_Model_Sample_SGED3 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub3,out.sample = 250)
WTI_Model_Sample_SGED4 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub4,out.sample = 250)
WTI_Model_Sample_SGED5 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub5,out.sample = 250)
WTI_Model_Sample_SGED6 <- ugarchfit(WTI_Spec_Sample_SGED,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

Brent_Model_Sample_SGED1 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub1,out.sample = 250)
Brent_Model_Sample_SGED2 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub2,out.sample = 250)
Brent_Model_Sample_SGED3 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub3,out.sample = 250)
Brent_Model_Sample_SGED4 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub4,out.sample = 250)
Brent_Model_Sample_SGED5 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub5,out.sample = 250)
Brent_Model_Sample_SGED6 <- ugarchfit(Brent_Spec_Sample_SGED,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_SGED_Sample1 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample2 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample3 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample4 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample5 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SGED_Sample6 <-  ugarchroll(WTI_Spec_Sample_SGED,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_SGED_Sample1 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample2 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample3 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample4 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample5 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SGED_Sample6 <-  ugarchroll(Brent_Spec_Sample_SGED,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### SST Specification
# 
WTI_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
# 
WTI_Model_Sample_SST1 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub1,out.sample = 250)
WTI_Model_Sample_SST2 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub2,out.sample = 250)
WTI_Model_Sample_SST3 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub3,out.sample = 250)
WTI_Model_Sample_SST4 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub4,out.sample = 250)
WTI_Model_Sample_SST5 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub5,out.sample = 250)
WTI_Model_Sample_SST6 <- ugarchfit(WTI_Spec_Sample_SST,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")

Brent_Model_Sample_SST1 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub1,out.sample = 250)
Brent_Model_Sample_SST2 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub2,out.sample = 250)
Brent_Model_Sample_SST3 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub3,out.sample = 250)
Brent_Model_Sample_SST4 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub4,out.sample = 250)
Brent_Model_Sample_SST5 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub5,out.sample = 250)
Brent_Model_Sample_SST6 <- ugarchfit(Brent_Spec_Sample_SST,Brent_sub6,out.sample = 250)


## Rollover Analysis

WTI_Roll_SST_Sample1 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample2 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample3 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample4 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample5 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_SST_Sample6 <-  ugarchroll(WTI_Spec_Sample_SST,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_SST_Sample1 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample2 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample3 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample4 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample5 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_SST_Sample6 <-  ugarchroll(Brent_Spec_Sample_SST,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### Normal Specification

WTI_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

WTI_Model_Sample_Normal1 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)
WTI_Model_Sample_Normal2 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)
WTI_Model_Sample_Normal3 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)
WTI_Model_Sample_Normal4 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)
WTI_Model_Sample_Normal5 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)
WTI_Model_Sample_Normal6 <- ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)

Brent_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Brent_Model_Sample_Normal1 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)
Brent_Model_Sample_Normal2 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)
Brent_Model_Sample_Normal3 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)
Brent_Model_Sample_Normal4 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)
Brent_Model_Sample_Normal5 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)
Brent_Model_Sample_Normal6 <- ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)

WTI_Roll_Normal_Sample1 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample2 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample3 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample4 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample5 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
WTI_Roll_Normal_Sample6 <-  ugarchroll(WTI_Spec_Sample_Normal,WTI_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Brent_Roll_Normal_Sample1 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample2 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample3 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample4 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample5 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Brent_Roll_Normal_Sample6 <-  ugarchroll(Brent_Spec_Sample_Normal,Brent_sub6,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### Sample Residuals Statistics

WTI_Res_Sample1 <- c(normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal1@fit$residuals/WTI_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample2 <- c(normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal2@fit$residuals/WTI_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample3 <- c(normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal3@fit$residuals/WTI_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample4 <- c(normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal4@fit$residuals/WTI_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample5 <- c(normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal5@fit$residuals/WTI_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)
WTI_Res_Sample6 <- c(normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="jb")@test$statistic,normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="sw")@test$statistic,normalTest(WTI_Model_Sample_Normal6@fit$residuals/WTI_Model_Sample_Normal6@fit$sigma,method="ks")@test$statistic)

Brent_Res_Sample1 <- c(normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal1@fit$residuals/Brent_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample2 <- c(normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal2@fit$residuals/Brent_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample3 <- c(normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal3@fit$residuals/Brent_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample4 <- c(normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal4@fit$residuals/Brent_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample5 <- c(normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal5@fit$residuals/Brent_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)
Brent_Res_Sample6 <- c(normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="jb")@test$statistic,normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="sw")@test$statistic,normalTest(Brent_Model_Sample_Normal6@fit$residuals/Brent_Model_Sample_Normal6@fit$sigma,method="ks")@test$statistic)

WTI_Res_Sample_Test <- rbind(WTI_Res_Sample1,WTI_Res_Sample2,WTI_Res_Sample3,WTI_Res_Sample4,WTI_Res_Sample5,WTI_Res_Sample6)
Brent_Res_Sample_Test   <- rbind(Brent_Res_Sample1,Brent_Res_Sample2,Brent_Res_Sample3,Brent_Res_Sample4,Brent_Res_Sample5,Brent_Res_Sample6)
rownames(WTI_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")
colnames(WTI_Res_Sample_Test) <- c("JB","SW","KS") 
rownames(Brent_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")
colnames(Brent_Res_Sample_Test) <- c("JB","SW","KS") 

### Pearson Specification

WTI_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$sigma)
WTI_Param_Pearson_Sample6 <- pearsonFitML(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$sigma)

WTI_Param_Pearson_Subsample <- rbind(WTI_Param_Pearson_Sample1,WTI_Param_Pearson_Sample2,WTI_Param_Pearson_Sample3,WTI_Param_Pearson_Sample4,WTI_Param_Pearson_Sample5,WTI_Param_Pearson_Sample6)
rownames(WTI_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

Brent_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$sigma)
Brent_Param_Pearson_Sample6 <- pearsonFitML(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$sigma)

Brent_Param_Pearson_Subsample <- rbind(Brent_Param_Pearson_Sample1,Brent_Param_Pearson_Sample2,Brent_Param_Pearson_Sample3,Brent_Param_Pearson_Sample4,Brent_Param_Pearson_Sample5,Brent_Param_Pearson_Sample6)
rownames(Brent_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

### johnson SU parameter


WTI_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub1,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub2,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub3,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub4,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub5,out.sample = 250)@fit$sigma,moment = "find")
WTI_Param_Johnson_Sample6 <- JohnsonFit(ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$residuals/ugarchfit(WTI_Spec_Sample_Normal,WTI_sub6,out.sample = 250)@fit$sigma,moment = "find")

WTI_Param_Johnson_Subsample <- rbind(WTI_Param_Johnson_Sample1,WTI_Param_Johnson_Sample2,WTI_Param_Johnson_Sample3,WTI_Param_Johnson_Sample4,WTI_Param_Johnson_Sample5,WTI_Param_Johnson_Sample6)
rownames(WTI_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

Brent_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub1,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub2,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub3,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub4,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub5,out.sample = 250)@fit$sigma,moment = "find")
Brent_Param_Johnson_Sample6 <- JohnsonFit(ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$residuals/ugarchfit(Brent_Spec_Sample_Normal,Brent_sub6,out.sample = 250)@fit$sigma,moment = "find")

Brent_Param_Johnson_Subsample <- rbind(Brent_Param_Johnson_Sample1,Brent_Param_Johnson_Sample2,Brent_Param_Johnson_Sample3,Brent_Param_Johnson_Sample4,Brent_Param_Johnson_Sample5,Brent_Param_Johnson_Sample6)
rownames(Brent_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample6")

### Plotting the VaR for subsamples

WTI_Sample1_Dates <- format(WTI$Date[(WTI_samples[1]-249):WTI_samples[1]],"%b-%y")
WTI_Sample2_Dates <- format(WTI$Date[(WTI_samples[2]-249):WTI_samples[2]],"%b-%y")
WTI_Sample3_Dates <- format(WTI$Date[(WTI_samples[3]-249):WTI_samples[3]],"%b-%y")
WTI_Sample4_Dates <- format(WTI$Date[(WTI_samples[4]-249):WTI_samples[4]],"%b-%y")
WTI_Sample5_Dates <- format(WTI$Date[(WTI_samples[5]-249):WTI_samples[5]],"%b-%y")
WTI_Sample6_Dates <- format(WTI$Date[(length(WTI_xts)-249):length(WTI_xts)],"%b-%y")

Brent_Sample1_Dates <- format(Brent$Date[(Brent_samples[1]-249):Brent_samples[1]],"%b-%y")
Brent_Sample2_Dates <- format(Brent$Date[(Brent_samples[2]-249):Brent_samples[2]],"%b-%y")
Brent_Sample3_Dates <- format(Brent$Date[(Brent_samples[3]-249):Brent_samples[3]],"%b-%y")
Brent_Sample4_Dates <- format(Brent$Date[(Brent_samples[4]-249):Brent_samples[4]],"%b-%y")
Brent_Sample5_Dates <- format(Brent$Date[(Brent_samples[5]-249):Brent_samples[5]],"%b-%y")
Brent_Sample6_Dates <- format(Brent$Date[(length(Brent_xts)-249):length(Brent_xts)],"%b-%y")

### Forecast Performance Measures ###

### Loss functions ###

# In Sample:

Brent_Temp_Returns <- Returns_Brent[1:4646]
Brent_Realized_Vol <- array(data=0,dim=4500)

for(i in 1:4500)
{
  Brent_Realized_Vol[i] <- sum(Brent_Temp_Returns[i:(i+146)]^2)
  #+Brent_Temp_Returns[i+1]^2+Brent_Temp_Returns[i+2]^2+Brent_Temp_Returns[i+3]^2+
  #Brent_Temp_Returns[i+4]^2+Brent_Temp_Returns[i+5]^2+Brent_Temp_Returns[i+6]^2+Brent_Temp_Returns[i+7]^2+
  #Brent_Temp_Returns[i+8]^2+Brent_Temp_Returns[i+9]^2
}

Brent_Realized_Vol  <- sqrt(Brent_Realized_Vol)  

Brent_In_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(Brent_In_Sample) <- c("Normal","JSU","SGED","SST")
rownames(Brent_In_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
Brent_In_Sample[1,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1")))
Brent_In_Sample[2,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE2")))
Brent_In_Sample[3,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="AE1")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="AE1")))
Brent_In_Sample[4,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="AE2")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="AE2")))
Brent_In_Sample[5,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="R2LOG")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="R2LOG")))
Brent_In_Sample[6,] <- cbind(sum(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="QLIKE")),sum(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="QLIKE")))
Brent_In_Sample <- round(Brent_In_Sample,6)


### SPA test

SPA_Brent_base_SE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1"),
                                     LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1"))
                               ,0.10,statistic = "Tmax")

SPA_Brent_base_Q <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="QLIKE"),
                                     LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="QLIKE"),
                                     LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="QLIKE"),
                                     LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="QLIKE"))
                               ,0.10,statistic = "Tmax")

SPA_Brent_base_R2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="R2LOG"),
                                      LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="R2LOG"),
                                      LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="R2LOG"),
                                      LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="R2LOG"))
                                ,0.10,statistic = "Tmax")

SPA_Brent_base_AE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="AE2"),
                                       LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="AE2"),
                                       LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="AE2"),
                                       LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="AE2"))
                                 ,0.10,statistic = "Tmax")

SPA_Brent_base_AE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="AE1"),
                                       LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="AE1"),
                                       LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="AE1"),
                                       LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="AE1"))
                                 ,0.10,statistic = "Tmax")

SPA_Brent_base_SE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE2"),
                                       LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE2"),
                                       LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE2"),
                                       LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE2"))
                                 ,0.10,statistic = "Tmax")


## Pairwise SPA

Brent_SPA_comb1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1"))),alpha = 0.05,statistic = "Tmax")
Brent_SPA_comb2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1"))),alpha = 0.05,statistic = "Tmax")
Brent_SPA_comb3 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model8@fit$sigma[147:4646],which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1"))),alpha = 0.05,statistic = "Tmax")
Brent_SPA_comb4 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1"))),alpha = 0.05,statistic = "Tmax")
Brent_SPA_comb5 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model7@fit$sigma[147:4646],which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1"))),alpha = 0.05,statistic = "Tmax")
Brent_SPA_comb6 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Model6@fit$sigma[147:4646],which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Model5@fit$sigma[147:4646],which="SE1"))),alpha = 0.05,statistic = "Tmax")

#### Store the p-values

Brent_SPA_pvalue <- matrix(data=NA,nrow=4,ncol=4)
colnames(Brent_SPA_pvalue) <- c("Norm","JSU","SGED","SST")
rownames(Brent_SPA_pvalue) <- c("Norm","JSU","SGED","SST")

Brent_SPA_pvalue[2,1] <- Brent_SPA_comb1@Info$mcs_pvalue
Brent_SPA_pvalue[3,1] <- Brent_SPA_comb2@Info$mcs_pvalue
Brent_SPA_pvalue[4,1] <- Brent_SPA_comb3@Info$mcs_pvalue
Brent_SPA_pvalue[3,2] <- Brent_SPA_comb4@Info$mcs_pvalue
Brent_SPA_pvalue[4,2] <- Brent_SPA_comb5@Info$mcs_pvalue
Brent_SPA_pvalue[4,3] <- Brent_SPA_comb6@Info$mcs_pvalue

Brent_Cont <- Brent_Model5@fit$sigma[147:4646]

# FOR WTI

WTI_Temp_Returns <- Returns_WTI[1:4537]
WTI_Realized_Vol <- array(data=0,dim=4400)

for(i in 1:4400)
{
  WTI_Realized_Vol[i] <- sum(WTI_Temp_Returns[i:(i+137)]^2)
  #+WTI_Temp_Returns[i+1]^2+WTI_Temp_Returns[i+2]^2+WTI_Temp_Returns[i+3]^2+
  #WTI_Temp_Returns[i+4]^2+WTI_Temp_Returns[i+5]^2+WTI_Temp_Returns[i+6]^2+WTI_Temp_Returns[i+7]^2+
  #WTI_Temp_Returns[i+8]^2+WTI_Temp_Returns[i+9]^2
}

WTI_Realized_Vol  <- sqrt(WTI_Realized_Vol)  

WTI_In_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(WTI_In_Sample) <- c("Normal","JSU","SGED","SST")
rownames(WTI_In_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
WTI_In_Sample[1,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1")))
WTI_In_Sample[2,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE2")))
WTI_In_Sample[3,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="AE1")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="AE1")))
WTI_In_Sample[4,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="AE2")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="AE2")))
WTI_In_Sample[5,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="R2LOG")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="R2LOG")))
WTI_In_Sample[6,] <- cbind(sum(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="QLIKE")),sum(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="QLIKE")))
WTI_In_Sample <- round(WTI_In_Sample,6)


### SPA test

SPA_WTI_base_SE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1"))
                             ,0.10,statistic = "Tmax")

SPA_WTI_base_Q <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="QLIKE"),
                                   LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="QLIKE"),
                                   LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="QLIKE"),
                                   LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="QLIKE"))
                             ,0.10,statistic = "Tmax")

SPA_WTI_base_R2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="R2LOG"),
                                   LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="R2LOG"),
                                   LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="R2LOG"),
                                   LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="R2LOG"))
                             ,0.10,statistic = "Tmax")

SPA_WTI_base_AE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="AE2"),
                                   LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="AE2"),
                                   LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="AE2"),
                                   LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="AE2"))
                             ,0.10,statistic = "Tmax")

SPA_WTI_base_AE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="AE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="AE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="AE1"),
                                   LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="AE1"))
                             ,0.10,statistic = "Tmax")

SPA_WTI_base_SE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE2"),
                                   LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE2"),
                                   LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE2"),
                                   LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE2"))
                             ,0.10,statistic = "Tmax")


## Pairwise SPA ### 

WTI_SPA_comb1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1"))),alpha = 0.05,statistic = "Tmax")
WTI_SPA_comb2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1"))),alpha = 0.05,statistic = "Tmax")
WTI_SPA_comb3 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model8@fit$sigma[138:4537],which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1"))),alpha = 0.05,statistic = "Tmax")
WTI_SPA_comb4 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1"))),alpha = 0.05,statistic = "Tmax")
WTI_SPA_comb5 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model7@fit$sigma[138:4537],which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1"))),alpha = 0.05,statistic = "Tmax")
WTI_SPA_comb6 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Model6@fit$sigma[138:4537],which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Model5@fit$sigma[138:4537],which="SE1"))),alpha = 0.05,statistic = "Tmax")

#### Store the p-values

WTI_SPA_pvalue <- matrix(data=NA,nrow=4,ncol=4)
colnames(WTI_SPA_pvalue) <- c("Norm","JSU","SGED","SST")
rownames(WTI_SPA_pvalue) <- c("Norm","JSU","SGED","SST")

WTI_SPA_pvalue[2,1] <- WTI_SPA_comb1@Info$mcs_pvalue
WTI_SPA_pvalue[3,1] <- WTI_SPA_comb2@Info$mcs_pvalue
WTI_SPA_pvalue[4,1] <- WTI_SPA_comb3@Info$mcs_pvalue
WTI_SPA_pvalue[3,2] <- WTI_SPA_comb4@Info$mcs_pvalue
WTI_SPA_pvalue[4,2] <- WTI_SPA_comb5@Info$mcs_pvalue
WTI_SPA_pvalue[4,3] <- WTI_SPA_comb6@Info$mcs_pvalue

WTI_Cont <- WTI_Model7@fit$sigma[138:4537]

# Out Sample

## For Brent
BrentOut_Temp_Returns <- Returns_Brent[(length(Returns_Brent)-999):length(Returns_Brent)]
Brent_Out_Realized_Vol <- array(data=0,dim=900)

for(i in 1:900)
{
  Brent_Out_Realized_Vol[i] <- sum(BrentOut_Temp_Returns[i:(i+99)]^2)
  #+BrentOut_Temp_Returns[i+1]^2+BrentOut_Temp_Returns[i+2]^2+BrentOut_Temp_Returns[i+3]^2+
  #BrentOut_Temp_Returns[i+4]^2+BrentOut_Temp_Returns[i+5]^2+BrentOut_Temp_Returns[i+6]^2+BrentOut_Temp_Returns[i+7]^2+
  #BrentOut_Temp_Returns[i+8]^2+BrentOut_Temp_Returns[i+9]^2
}

Brent_Out_Realized_Vol  <- sqrt(Brent_Out_Realized_Vol)  

Brent_Out_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(Brent_Out_Sample) <- c("Normal","JSU","SGED","SST")
rownames(Brent_Out_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
Brent_Out_Sample[1,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE1")))
Brent_Out_Sample[2,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE2")))
Brent_Out_Sample[3,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="AE1")))
Brent_Out_Sample[4,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="AE2")))
Brent_Out_Sample[5,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="R2LOG")))
Brent_Out_Sample[6,] <- cbind(sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="QLIKE")))
Brent_Out_Sample <- round(Brent_Out_Sample,6)


### SPA test

Brent_SPA_out <- MCSprocedure(cbind(LossVol(Brent_Out_Realized_Vol,Brent_Roll8@forecast$density$Sigma[100:999],which="SE1"),
                                    LossVol(Brent_Out_Realized_Vol,Brent_Roll7@forecast$density$Sigma[100:999],which="SE1"),
                                    LossVol(Brent_Out_Realized_Vol,Brent_Roll6@forecast$density$Sigma[100:999],which="SE1"),
                                    LossVol(Brent_Out_Realized_Vol,Brent_Roll5@forecast$density$Sigma[100:999],which="SE1"))
                              ,0.05,statistic = "Tmax")

# FOR WTI

WTIOut_Temp_Returns <- Returns_WTI[(length(Returns_WTI)-999):length(Returns_WTI)]
WTI_Out_Realized_Vol <- array(data=0,dim=900)

for(i in 1:900)
{
  WTI_Out_Realized_Vol[i] <- sum(WTIOut_Temp_Returns[i:(i+99)]^2)
  #+WTIOut_Temp_Returns[i+1]^2+WTIOut_Temp_Returns[i+2]^2+WTIOut_Temp_Returns[i+3]^2+
  #WTIOut_Temp_Returns[i+4]^2+WTIOut_Temp_Returns[i+5]^2+WTIOut_Temp_Returns[i+6]^2+WTIOut_Temp_Returns[i+7]^2+
  #WTIOut_Temp_Returns[i+8]^2+WTIOut_Temp_Returns[i+9]^2
}

WTI_Out_Realized_Vol  <- sqrt(WTI_Out_Realized_Vol)  

WTI_Out_Sample <- matrix(data=0,nrow=6,ncol=4)
colnames(WTI_Out_Sample) <- c("Normal","JSU","SGED","SST")
rownames(WTI_Out_Sample) <- c("SE1","SE2","AE1","AE2","R2LOG","QLIKE")
WTI_Out_Sample[1,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE1")))
WTI_Out_Sample[2,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE2")))
WTI_Out_Sample[3,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="AE1")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="AE1")))
WTI_Out_Sample[4,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="AE2")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="AE2")))
WTI_Out_Sample[5,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="R2LOG")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="R2LOG")))
WTI_Out_Sample[6,] <- cbind(sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="QLIKE")),sum(LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="QLIKE")))
WTI_Out_Sample <- round(WTI_Out_Sample,6)


### SPA test

WTI_SPA_out <- MCSprocedure(cbind(LossVol(WTI_Out_Realized_Vol,WTI_Roll8@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(WTI_Out_Realized_Vol,WTI_Roll7@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(WTI_Out_Realized_Vol,WTI_Roll6@forecast$density$Sigma[100:999],which="SE1"),
                                  LossVol(WTI_Out_Realized_Vol,WTI_Roll5@forecast$density$Sigma[100:999],which="SE1"))
                            ,0.05,statistic = "Tmax")


#### SPA OVERALLL ####

## Pairwise SPA comparing models with/without volume### 

WTI_SPA_comb11 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Base,which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Cont,which="SE1"))),alpha = 0.05,statistic = "Tmax")
WTI_SPA_comb21 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Base,which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Lag,which="SE1"))),alpha = 0.05,statistic = "Tmax")
WTI_SPA_comb31 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Cont,which="SE1"),(LossVol(WTI_Realized_Vol,WTI_Lag,which="SE1"))),alpha = 0.05,statistic = "Tmax")

Brent_SPA_comb11 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Base,which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Cont,which="SE1"))),alpha = 0.05,statistic = "Tmax")
Brent_SPA_comb21 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Base,which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Lag,which="SE1"))),alpha = 0.05,statistic = "Tmax")
Brent_SPA_comb31 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Cont,which="SE1"),(LossVol(Brent_Realized_Vol,Brent_Lag,which="SE1"))),alpha = 0.05,statistic = "Tmax")


# ### CALCULATING MSE, MAE AND DAE for in-sample
# 
# In_Brent_MSE_SSTD <- sum((Brent_Roll5@forecast$density$Mu-Brent_Roll5@forecast$density$Realized)^2)/1000
# In_Brent_MSE_SGED <- sum((Brent_Roll6@forecast$density$Mu-Brent_Roll6@forecast$density$Realized)^2)/1000
# In_Brent_MSE_JSU  <- sum((Brent_Roll7@forecast$density$Mu-Brent_Roll7@forecast$density$Realized)^2)/1000
# In_Brent_MSE_Norm <- sum((Brent_Roll8@forecast$density$Mu-Brent_Roll8@forecast$density$Realized)^2)/1000
# 
# In_WTI_MSE_SSTD <- sum((WTI_Roll5@forecast$density$Mu-WTI_Roll5@forecast$density$Realized)^2)/1000
# In_WTI_MSE_SGED <- sum((WTI_Roll6@forecast$density$Mu-WTI_Roll6@forecast$density$Realized)^2)/1000
# In_WTI_MSE_JSU  <- sum((WTI_Roll7@forecast$density$Mu-WTI_Roll7@forecast$density$Realized)^2)/1000
# In_WTI_MSE_Norm <- sum((WTI_Roll8@forecast$density$Mu-WTI_Roll8@forecast$density$Realized)^2)/1000
# 
# In_Brent_MSE <- c(In_Brent_MSE_SSTD,In_Brent_MSE_SGED,In_Brent_MSE_JSU,In_Brent_MSE_Norm)
# In_WTI_MSE <- c(In_WTI_MSE_SSTD,In_WTI_MSE_SGED,In_WTI_MSE_JSU,In_WTI_MSE_Norm)
# In_MSE <- cbind(In_Brent_MSE,In_WTI_MSE)
# rownames(In_MSE) <- c("SSTD","SGED","JSU","Norm")
# 
# ## Calculating MAE
# 
# In_Brent_MAE_SSTD <- sum(abs(Brent_Roll5@forecast$density$Mu-Brent_Roll5@forecast$density$Realized))/1000
# In_Brent_MAE_SGED <- sum(abs(Brent_Roll6@forecast$density$Mu-Brent_Roll6@forecast$density$Realized))/1000
# In_Brent_MAE_JSU  <- sum(abs(Brent_Roll7@forecast$density$Mu-Brent_Roll7@forecast$density$Realized))/1000
# In_Brent_MAE_Norm <- sum(abs(Brent_Roll8@forecast$density$Mu-Brent_Roll8@forecast$density$Realized))/1000
# 
# In_WTI_MAE_SSTD <- sum(abs(WTI_Roll5@forecast$density$Mu-WTI_Roll5@forecast$density$Realized))/1000
# In_WTI_MAE_SGED <- sum(abs(WTI_Roll6@forecast$density$Mu-WTI_Roll6@forecast$density$Realized))/1000
# In_WTI_MAE_JSU  <- sum(abs(WTI_Roll7@forecast$density$Mu-WTI_Roll7@forecast$density$Realized))/1000
# In_WTI_MAE_Norm <- sum(abs(WTI_Roll8@forecast$density$Mu-WTI_Roll8@forecast$density$Realized))/1000
# 
# In_Brent_MAE <- c(In_Brent_MAE_SSTD,In_Brent_MAE_SGED,In_Brent_MAE_JSU,In_Brent_MAE_Norm)
# In_WTI_MAE <- c(In_WTI_MAE_SSTD,In_WTI_MAE_SGED,In_WTI_MAE_JSU,In_WTI_MAE_Norm)
# In_MAE <- cbind(In_Brent_MAE,In_WTI_MAE)
# rownames(In_MAE) <- c("SSTD","SGED","JSU","Norm")
# 
# 
# ### CALCULATING MSE, MAE AND DAE for out of sample
# 
# Brent_MSE_SSTD <- sum((Brent_Roll5@forecast$density$Mu-Brent_Roll5@forecast$density$Realized)^2)/1000
# Brent_MSE_SGED <- sum((Brent_Roll6@forecast$density$Mu-Brent_Roll6@forecast$density$Realized)^2)/1000
# Brent_MSE_JSU  <- sum((Brent_Roll7@forecast$density$Mu-Brent_Roll7@forecast$density$Realized)^2)/1000
# Brent_MSE_Norm <- sum((Brent_Roll8@forecast$density$Mu-Brent_Roll8@forecast$density$Realized)^2)/1000
# 
# WTI_MSE_SSTD <- sum((WTI_Roll5@forecast$density$Mu-WTI_Roll5@forecast$density$Realized)^2)/1000
# WTI_MSE_SGED <- sum((WTI_Roll6@forecast$density$Mu-WTI_Roll6@forecast$density$Realized)^2)/1000
# WTI_MSE_JSU  <- sum((WTI_Roll7@forecast$density$Mu-WTI_Roll7@forecast$density$Realized)^2)/1000
# WTI_MSE_Norm <- sum((WTI_Roll8@forecast$density$Mu-WTI_Roll8@forecast$density$Realized)^2)/1000
# 
# Brent_MSE <- c(Brent_MSE_SSTD,Brent_MSE_SGED,Brent_MSE_JSU,Brent_MSE_Norm)
# WTI_MSE <- c(WTI_MSE_SSTD,WTI_MSE_SGED,WTI_MSE_JSU,WTI_MSE_Norm)
# MSE <- cbind(Brent_MSE,WTI_MSE)
# rownames(MSE) <- c("SSTD","SGED","JSU","Norm")
# 
# ## Calculating MAE
# 
# Brent_MAE_SSTD <- sum(abs(Brent_Roll5@forecast$density$Mu-Brent_Roll5@forecast$density$Realized))/1000
# Brent_MAE_SGED <- sum(abs(Brent_Roll6@forecast$density$Mu-Brent_Roll6@forecast$density$Realized))/1000
# Brent_MAE_JSU  <- sum(abs(Brent_Roll7@forecast$density$Mu-Brent_Roll7@forecast$density$Realized))/1000
# Brent_MAE_Norm <- sum(abs(Brent_Roll8@forecast$density$Mu-Brent_Roll8@forecast$density$Realized))/1000
# 
# WTI_MAE_SSTD <- sum(abs(WTI_Roll5@forecast$density$Mu-WTI_Roll5@forecast$density$Realized))/1000
# WTI_MAE_SGED <- sum(abs(WTI_Roll6@forecast$density$Mu-WTI_Roll6@forecast$density$Realized))/1000
# WTI_MAE_JSU  <- sum(abs(WTI_Roll7@forecast$density$Mu-WTI_Roll7@forecast$density$Realized))/1000
# WTI_MAE_Norm <- sum(abs(WTI_Roll8@forecast$density$Mu-WTI_Roll8@forecast$density$Realized))/1000
# 
# Brent_MAE <- c(Brent_MAE_SSTD,Brent_MAE_SGED,Brent_MAE_JSU,Brent_MAE_Norm)
# WTI_MAE <- c(WTI_MAE_SSTD,WTI_MAE_SGED,WTI_MAE_JSU,WTI_MAE_Norm)
# MAE <- cbind(Brent_MAE,WTI_MAE)
# rownames(MAE) <- c("SSTD","SGED","JSU","Norm")
# 
# The same can be done by using report(,type="fpm")

report(Brent_Roll5,type="fpm")
report(Brent_Roll6,type="fpm")
report(Brent_Roll7,type="fpm")
report(Brent_Roll8,type="fpm")

report(WTI_Roll5,type="fpm")
report(WTI_Roll6,type="fpm")
report(WTI_Roll7,type="fpm")
report(WTI_Roll8,type="fpm")


#### Sub-sample Analysis of MCS ####

### FOR BRENT

## Sub-sample 1

size1 <- (length(Brent_Model_Sample_Normal1@fit$sigma)-100)
Brent_Temp_Returns_Sample1 <- Returns_Brent[1:(size1+100)]
Brent_Realized_Vol_Sample1 <- array(data=0,dim=size1)

for(i in 1:size1)
{
  Brent_Realized_Vol_Sample1[i] <- sum(Brent_Temp_Returns_Sample1[i:(i+99)]^2)
}

Brent_Realized_Vol_Sample1  <- sqrt(Brent_Realized_Vol_Sample1)  

### SPA test

Brent_SPA_Sample1_SE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="SE1"),
                                    LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="SE1"),
                                    LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="SE1"),
                                    LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="SE1"))
                              ,0.10,statistic = "Tmax")

Brent_SPA_Sample1_SE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample1_AE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample1_AE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample1_R2LOG <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="R2LOG"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="R2LOG"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="R2LOG"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="R2LOG"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample1_QLIKE <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="QLIKE"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="QLIKE"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="QLIKE"),
                                            LossVol(Brent_Realized_Vol_Sample1,Brent_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="QLIKE"))
                                      ,0.10,statistic = "Tmax")

## Sub-sample2

size2 <- (length(Brent_Model_Sample_Normal2@fit$sigma)-100)
Brent_Temp_Returns_Sample2 <- Returns_Brent[(size1+350+1):(size1+350+size2+100)]
Brent_Realized_Vol_Sample2 <- array(data=0,dim=size2)

for(i in 1:size2)
{
  Brent_Realized_Vol_Sample2[i] <- sum(Brent_Temp_Returns_Sample2[i:(i+99)]^2)
}

Brent_Realized_Vol_Sample2  <- sqrt(Brent_Realized_Vol_Sample2)  

### SPA test

Brent_SPA_Sample2_SE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_Normal2@fit$sigma[100:(size2+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


Brent_SPA_Sample2_SE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample2_AE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample2_AE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample2_R2LOG <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

Brent_SPA_Sample2_QLIKE <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample2,Brent_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")


## Sub-sample3

size3 <- (length(Brent_Model_Sample_Normal3@fit$sigma)-100)
Brent_Temp_Returns_Sample3 <- Returns_Brent[(size2+350+1):(size2+350+size3+100)]
Brent_Realized_Vol_Sample3 <- array(data=0,dim=size3)

for(i in 1:size3)
{
  Brent_Realized_Vol_Sample3[i] <- sum(Brent_Temp_Returns_Sample3[i:(i+99)]^2)
}

Brent_Realized_Vol_Sample3  <- sqrt(Brent_Realized_Vol_Sample3)  

### SPA test

Brent_SPA_Sample3_SE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_Normal3@fit$sigma[100:(size3+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


Brent_SPA_Sample3_SE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample3_AE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample3_AE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample3_R2LOG <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

Brent_SPA_Sample3_QLIKE <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample3,Brent_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")

## Sub-sample4

size4 <- (length(Brent_Model_Sample_Normal4@fit$sigma)-100)
Brent_Temp_Returns_Sample4 <- Returns_Brent[(size3+350+1):(size3+350+size4+100)]
Brent_Realized_Vol_Sample4 <- array(data=0,dim=size4)

for(i in 1:size4)
{
  Brent_Realized_Vol_Sample4[i] <- sum(Brent_Temp_Returns_Sample4[i:(i+99)]^2)
}

Brent_Realized_Vol_Sample4  <- sqrt(Brent_Realized_Vol_Sample4)  

### SPA test

Brent_SPA_Sample4_SE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


Brent_SPA_Sample4_SE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample4_AE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample4_AE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample4_R2LOG <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

Brent_SPA_Sample4_QLIKE <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample4,Brent_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")

## Sub-Sample5

size5 <- (length(Brent_Model_Sample_Normal5@fit$sigma)-100)
Brent_Temp_Returns_Sample5 <- Returns_Brent[(size1+350+1):(size1+350+size5+100)]
Brent_Realized_Vol_Sample5 <- array(data=0,dim=size5)

for(i in 1:size5)
{
  Brent_Realized_Vol_Sample5[i] <- sum(Brent_Temp_Returns_Sample5[i:(i+99)]^2)
}

Brent_Realized_Vol_Sample5  <- sqrt(Brent_Realized_Vol_Sample5)  

### SPA test

Brent_SPA_Sample5_SE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


Brent_SPA_Sample5_SE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample5_AE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample5_AE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample5_R2LOG <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

Brent_SPA_Sample5_QLIKE <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample5,Brent_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")


## Sub-Sample6

size6 <- (length(Brent_Model_Sample_Normal6@fit$sigma)-100)
Brent_Temp_Returns_Sample6 <- Returns_Brent[(size5+350+1):(size5+350+size6+100)]
Brent_Realized_Vol_Sample6 <- array(data=0,dim=size6)

for(i in 1:size6)
{
  Brent_Realized_Vol_Sample6[i] <- sum(Brent_Temp_Returns_Sample6[i:(i+99)]^2)
}

Brent_Realized_Vol_Sample6  <- sqrt(Brent_Realized_Vol_Sample6)  

### SPA test

Brent_SPA_Sample6_SE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="SE1"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


Brent_SPA_Sample6_SE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="SE2"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample6_AE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="AE1"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample6_AE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="AE2"),
                                            LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

Brent_SPA_Sample6_R2LOG <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="R2LOG"),
                                              LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

Brent_SPA_Sample6_QLIKE <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="QLIKE"),
                                              LossVol(Brent_Realized_Vol_Sample6,Brent_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")



### FOR WTI

## Sub-sample 1

size1 <- (length(WTI_Model_Sample_Normal1@fit$sigma)-100)
WTI_Temp_Returns_Sample1 <- Returns_WTI[1:(size1+100)]
WTI_Realized_Vol_Sample1 <- array(data=0,dim=size1)

for(i in 1:size1)
{
  WTI_Realized_Vol_Sample1[i] <- sum(WTI_Temp_Returns_Sample1[i:(i+99)]^2)
}

WTI_Realized_Vol_Sample1  <- sqrt(WTI_Realized_Vol_Sample1)  

### SPA test

WTI_SPA_Sample1_SE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample1_SE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample1_AE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample1_AE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample1_R2LOG <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

WTI_SPA_Sample1_QLIKE <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_Normal1@fit$sigma[100:(size1+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_JSU1@fit$sigma[100:(size1+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SGED1@fit$sigma[100:(size1+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample1,WTI_Model_Sample_SST1@fit$sigma[100:(size1+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")

## Sub-sample2

size2 <- (length(WTI_Model_Sample_Normal2@fit$sigma)-100)
WTI_Temp_Returns_Sample2 <- Returns_WTI[(size1+350+1):(size1+350+size2+100)]
WTI_Realized_Vol_Sample2 <- array(data=0,dim=size2)

for(i in 1:size2)
{
  WTI_Realized_Vol_Sample2[i] <- sum(WTI_Temp_Returns_Sample2[i:(i+99)]^2)
}

WTI_Realized_Vol_Sample2  <- sqrt(WTI_Realized_Vol_Sample2)  

### SPA test

WTI_SPA_Sample2_SE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_Normal2@fit$sigma[100:(size2+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


WTI_SPA_Sample2_SE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample2_AE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample2_AE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample2_R2LOG <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

WTI_SPA_Sample2_QLIKE <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_Normal1@fit$sigma[100:(size2+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_JSU2@fit$sigma[100:(size2+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SGED2@fit$sigma[100:(size2+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample2,WTI_Model_Sample_SST2@fit$sigma[100:(size2+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")


## Sub-sample3

size3 <- (length(WTI_Model_Sample_Normal3@fit$sigma)-100)
WTI_Temp_Returns_Sample3 <- Returns_WTI[(size2+350+1):(size2+350+size3+100)]
WTI_Realized_Vol_Sample3 <- array(data=0,dim=size3)

for(i in 1:size3)
{
  WTI_Realized_Vol_Sample3[i] <- sum(WTI_Temp_Returns_Sample3[i:(i+99)]^2)
}

WTI_Realized_Vol_Sample3  <- sqrt(WTI_Realized_Vol_Sample3)  

### SPA test

WTI_SPA_Sample3_SE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_Normal3@fit$sigma[100:(size3+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


WTI_SPA_Sample3_SE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample3_AE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample3_AE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample3_R2LOG <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

WTI_SPA_Sample3_QLIKE <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_Normal1@fit$sigma[100:(size3+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_JSU3@fit$sigma[100:(size3+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SGED3@fit$sigma[100:(size3+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample3,WTI_Model_Sample_SST3@fit$sigma[100:(size3+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")

## Sub-sample4

size4 <- (length(WTI_Model_Sample_Normal4@fit$sigma)-100)
WTI_Temp_Returns_Sample4 <- Returns_WTI[(size3+350+1):(size3+350+size4+100)]
WTI_Realized_Vol_Sample4 <- array(data=0,dim=size4)

for(i in 1:size4)
{
  WTI_Realized_Vol_Sample4[i] <- sum(WTI_Temp_Returns_Sample4[i:(i+99)]^2)
}

WTI_Realized_Vol_Sample4  <- sqrt(WTI_Realized_Vol_Sample4)  

### SPA test

WTI_SPA_Sample4_SE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


WTI_SPA_Sample4_SE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample4_AE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample4_AE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample4_R2LOG <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

WTI_SPA_Sample4_QLIKE <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_Normal4@fit$sigma[100:(size4+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_JSU4@fit$sigma[100:(size4+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SGED4@fit$sigma[100:(size4+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample4,WTI_Model_Sample_SST4@fit$sigma[100:(size4+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")

## Sub-Sample5

size5 <- (length(WTI_Model_Sample_Normal5@fit$sigma)-100)
WTI_Temp_Returns_Sample5 <- Returns_WTI[(size1+350+1):(size1+350+size5+100)]
WTI_Realized_Vol_Sample5 <- array(data=0,dim=size5)

for(i in 1:size5)
{
  WTI_Realized_Vol_Sample5[i] <- sum(WTI_Temp_Returns_Sample5[i:(i+99)]^2)
}

WTI_Realized_Vol_Sample5  <- sqrt(WTI_Realized_Vol_Sample5)  

### SPA test

WTI_SPA_Sample5_SE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


WTI_SPA_Sample5_SE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample5_AE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample5_AE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample5_R2LOG <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

WTI_SPA_Sample5_QLIKE <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_Normal5@fit$sigma[100:(size5+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_JSU5@fit$sigma[100:(size5+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SGED5@fit$sigma[100:(size5+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample5,WTI_Model_Sample_SST5@fit$sigma[100:(size5+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")


## Sub-Sample6

size6 <- (length(WTI_Model_Sample_Normal6@fit$sigma)-100)
WTI_Temp_Returns_Sample6 <- Returns_WTI[(size5+350+1):(size5+350+size6+100)]
WTI_Realized_Vol_Sample6 <- array(data=0,dim=size6)

for(i in 1:size6)
{
  WTI_Realized_Vol_Sample6[i] <- sum(WTI_Temp_Returns_Sample6[i:(i+99)]^2)
}

WTI_Realized_Vol_Sample6  <- sqrt(WTI_Realized_Vol_Sample6)  

### SPA test

WTI_SPA_Sample6_SE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="SE1"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="SE1"))
                                      ,0.10,statistic = "Tmax")


WTI_SPA_Sample6_SE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="SE2"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="SE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample6_AE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="AE1"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="AE1"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample6_AE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="AE2"),
                                            LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="AE2"))
                                      ,0.10,statistic = "Tmax")

WTI_SPA_Sample6_R2LOG <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="R2LOG"),
                                              LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="R2LOG"))
                                        ,0.10,statistic = "Tmax")

WTI_SPA_Sample6_QLIKE <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_Normal6@fit$sigma[100:(size6+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_JSU6@fit$sigma[100:(size6+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SGED6@fit$sigma[100:(size6+100-1)],which="QLIKE"),
                                              LossVol(WTI_Realized_Vol_Sample6,WTI_Model_Sample_SST6@fit$sigma[100:(size6+100-1)],which="QLIKE"))
                                        ,0.10,statistic = "Tmax")


### Overall MCS for models ###

### For Brent

Brent_model_SE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Base,which="SE1"),
                                    LossVol(Brent_Realized_Vol,Brent_Cont,which="SE1"),
                                    LossVol(Brent_Realized_Vol,Brent_Lag,which="SE1")))

Brent_model_SE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Base,which="SE2"),
                                    LossVol(Brent_Realized_Vol,Brent_Cont,which="SE2"),
                                    LossVol(Brent_Realized_Vol,Brent_Lag,which="SE2")))

Brent_model_AE1 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Base,which="AE1"),
                                    LossVol(Brent_Realized_Vol,Brent_Cont,which="AE1"),
                                    LossVol(Brent_Realized_Vol,Brent_Lag,which="AE1")))

Brent_model_AE2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Base,which="AE2"),
                                    LossVol(Brent_Realized_Vol,Brent_Cont,which="AE2"),
                                    LossVol(Brent_Realized_Vol,Brent_Lag,which="AE2")))

Brent_model_R2 <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Base,which="R2LOG"),
                                   LossVol(Brent_Realized_Vol,Brent_Cont,which="R2LOG"),
                                   LossVol(Brent_Realized_Vol,Brent_Lag,which="R2LOG")))

Brent_model_Q <- MCSprocedure(cbind(LossVol(Brent_Realized_Vol,Brent_Base,which="QLIKE"),
                                  LossVol(Brent_Realized_Vol,Brent_Cont,which="QLIKE"),
                                  LossVol(Brent_Realized_Vol,Brent_Lag,which="QLIKE")))

## For WTI 
WTI_model_SE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Base,which="SE1"),
                   LossVol(WTI_Realized_Vol,WTI_Cont,which="SE1"),
                   LossVol(WTI_Realized_Vol,WTI_Lag,which="SE1")))

WTI_model_SE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Base,which="SE2"),
                   LossVol(WTI_Realized_Vol,WTI_Cont,which="SE2"),
                   LossVol(WTI_Realized_Vol,WTI_Lag,which="SE2")))

WTI_model_AE1 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Base,which="AE1"),
                   LossVol(WTI_Realized_Vol,WTI_Cont,which="AE1"),
                   LossVol(WTI_Realized_Vol,WTI_Lag,which="AE1")))

WTI_model_AE2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Base,which="AE2"),
                   LossVol(WTI_Realized_Vol,WTI_Cont,which="AE2"),
                   LossVol(WTI_Realized_Vol,WTI_Lag,which="AE2")))

WTI_model_R2 <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Base,which="R2LOG"),
                   LossVol(WTI_Realized_Vol,WTI_Cont,which="R2LOG"),
                   LossVol(WTI_Realized_Vol,WTI_Lag,which="R2LOG")))

WTI_model_Q <- MCSprocedure(cbind(LossVol(WTI_Realized_Vol,WTI_Base,which="QLIKE"),
                   LossVol(WTI_Realized_Vol,WTI_Cont,which="QLIKE"),
                   LossVol(WTI_Realized_Vol,WTI_Lag,which="QLIKE")))

#### End of all the code/lines ######  


