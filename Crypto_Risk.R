#### IRFA VaR for Crypto Markets #####

### Title: Estimating risk in crypto-currencies with structural breaks: The role of distributions ###

## Job is to find the 1-day ahead VaR for the crypto-currrency markets ####
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

### We will use both Bitcoin and Litecoin Oil data series ###

Bitcoin <- read.csv("~//Desktop//Woxsen//Bitcoin.csv")
Litecoin <- read.csv("~//Desktop//Woxsen//Litecoin.csv")  
Litecoin_Date <- as.Date.character(Litecoin$Date,format = "%m/%d/%Y")
Bitcoin_Date   <- as.Date.character(Bitcoin$Date,format = "%m/%d/%Y")
Bitcoin   <- data.frame(Bitcoin_Date,Bitcoin$Bitcoin)     
Litecoin <- data.frame(Litecoin_Date,Litecoin$Litecoin)
colnames(Bitcoin) <- c("Date","Price")
colnames(Litecoin) <- c("Date","Price")

#### Structural Breaks for Litecoin and Bitcoin #####

Litecoin_breaks <- breakpoints(Litecoin$Price~1)
Bitcoin_breaks   <- breakpoints(Bitcoin$Price~1)
Litecoin_samples <- Litecoin_breaks$breakpoints
Bitcoin_samples <- Bitcoin_breaks$breakpoints

Returns_Bitcoin   <- log(Bitcoin$Price[-1]/Bitcoin$Price[1:(nrow(Bitcoin)-1)])
Returns_Litecoin <- log(Litecoin$Price[-1]/Litecoin$Price[1:(nrow(Litecoin)-1)])
Des_Returns_Bitcoin <- basicStats(Returns_Bitcoin)
Des_Returns_Litecoin <- basicStats(Returns_Litecoin)

Holdout_Returns_Bitcoin   <- Returns_Bitcoin[(length(Returns_Bitcoin)-499):length(Returns_Bitcoin)] 
Holdout_Returns_Litecoin <- Returns_Litecoin[(length(Returns_Litecoin)-499):length(Returns_Litecoin)]

#### Convert into xts objects

Litecoin_xts <- xts(Returns_Litecoin,order.by = Litecoin_Date[-1])
Bitcoin_xts   <- xts(Returns_Bitcoin,order.by = Bitcoin_Date[-1])
names(Litecoin_xts)    <- "Returns"
names(Bitcoin_xts)   <- "Returns"

### Subsamples XTS  ###

Litecoin_sub1    <- Litecoin_xts[1:Litecoin_samples[1]]
Litecoin_sub2    <- Litecoin_xts[Litecoin_samples[1]:Litecoin_samples[2]]
Litecoin_sub3    <- Litecoin_xts[Litecoin_samples[2]:Litecoin_samples[3]]
Litecoin_sub4    <- Litecoin_xts[Litecoin_samples[3]:nrow(Litecoin_xts)]

Bitcoin_sub1    <- Bitcoin_xts[1:Bitcoin_samples[1]]
Bitcoin_sub2    <- Bitcoin_xts[Bitcoin_samples[1]:Bitcoin_samples[2]]
Bitcoin_sub3    <- Bitcoin_xts[Bitcoin_samples[2]:nrow(Bitcoin_xts)]

###Dummy variables for structural breaks ####

Litecoin_Dummy1 <- array(data=0,dim=length(Litecoin_xts))
Litecoin_Dummy2 <- array(data=0,dim=length(Litecoin_xts))
Litecoin_Dummy3 <- array(data=0,dim=length(Litecoin_xts))
Litecoin_Dummy4 <- array(data=0,dim=length(Litecoin_xts))

Litecoin_Dummy1[1:Litecoin_samples[1]]  <- 1
Litecoin_Dummy2[(Litecoin_samples[1]+1):Litecoin_samples[2]] <- 1
Litecoin_Dummy3[(Litecoin_samples[2]+1):Litecoin_samples[3]] <- 1
Litecoin_Dummy4[(Litecoin_samples[3]+1):length(Litecoin_xts)] <- 1

Bitcoin_Dummy1 <- array(data=0,dim=length(Bitcoin_xts))
Bitcoin_Dummy2 <- array(data=0,dim=length(Bitcoin_xts))
Bitcoin_Dummy3 <- array(data=0,dim=length(Bitcoin_xts))

Bitcoin_Dummy1[1:Bitcoin_samples[1]]  <- 1
Bitcoin_Dummy2[(Bitcoin_samples[1]+1):Bitcoin_samples[2]] <- 1
Bitcoin_Dummy3[(Bitcoin_samples[2]+1):length(Bitcoin_xts)] <- 1

Litecoin_Dummy <- cbind(Litecoin_Dummy1,Litecoin_Dummy2,Litecoin_Dummy3)
Bitcoin_Dummy <- cbind(Bitcoin_Dummy1,Bitcoin_Dummy2)

### Descriptive Statistics and Basic Tests ###

Descr_Litecoin <-  cbind(basicStats(Litecoin_sub1),basicStats(Litecoin_sub2),basicStats(Litecoin_sub3),basicStats(Litecoin_sub4))
Descr_Bitcoin   <- cbind(basicStats(Bitcoin_sub1),basicStats(Bitcoin_sub2),basicStats(Bitcoin_sub3))
ADF_Litecoin <- cbind(adf.test(Litecoin_sub1),adf.test(Litecoin_sub2),adf.test(Litecoin_sub3),adf.test(Litecoin_sub4))
ADF_Bitcoin   <- cbind(adf.test(Bitcoin_sub1),adf.test(Bitcoin_sub2),adf.test(Bitcoin_sub3))
pp_Litecoin <- cbind(pp.test(Litecoin_sub1),pp.test(Litecoin_sub2),pp.test(Litecoin_sub3),pp.test(Litecoin_sub4))
pp_Bitcoin   <- cbind(pp.test(Bitcoin_sub1),pp.test(Bitcoin_sub2),pp.test(Bitcoin_sub3))
jb_Litecoin <- cbind(jarque.bera.test(Litecoin_sub1),jarque.bera.test(Litecoin_sub2),jarque.bera.test(Litecoin_sub3),jarque.bera.test(Litecoin_sub4))
jb_Bitcoin   <- cbind(jarque.bera.test(Bitcoin_sub1),jarque.bera.test(Bitcoin_sub2),jarque.bera.test(Bitcoin_sub3))

LQ_Litecoin5 <- cbind(Box.test(Litecoin_sub1^2,lag=5,type="Lj"),Box.test(Litecoin_sub2^2,lag=5,type="Lj"),Box.test(Litecoin_sub3^2,lag=5,type="Lj"),Box.test(Litecoin_sub4^2,lag=5,type="Lj"))
LQ_Bitcoin5   <- cbind(Box.test(Bitcoin_sub1^2,lag=5,type="Lj"),Box.test(Bitcoin_sub2^2,lag=5,type="Lj"),Box.test(Bitcoin_sub3^2,lag=5,type="Lj"))

LQ_Litecoin10 <- cbind(Box.test(Litecoin_sub1^2,lag=10,type="Lj"),Box.test(Litecoin_sub2^2,lag=10,type="Lj"),Box.test(Litecoin_sub3^2,lag=10,type="Lj"),Box.test(Litecoin_sub4^2,lag=10,type="Lj"))
LQ_Bitcoin10   <- cbind(Box.test(Bitcoin_sub1^2,lag=10,type="Lj"),Box.test(Bitcoin_sub2^2,lag=10,type="Lj"),Box.test(Bitcoin_sub3^2,lag=10,type="Lj"))

#### Plot the Returns and Oil Prices with Breakpoints indicated

par(mfrow=c(1,2))
plot(Litecoin,type="l",col="blue",main="Structural Breaks for Litecoin",ylab="Price (USD)",xaxt="n",xlab="")
abline(v=Litecoin$Date[Litecoin_samples[1]],col="brown")
abline(v=Litecoin$Date[Litecoin_samples[2]],col="brown")
abline(v=Litecoin$Date[Litecoin_samples[3]],col="brown")
Date_formatted <- format(Litecoin$Date[c(1,Litecoin_samples,nrow(Litecoin))],"%b-%y")
axis(1,at=c(min(Litecoin$Date),Litecoin$Date[Litecoin_samples],max(Litecoin$Date)),labels=Date_formatted,las = 2)

plot(Bitcoin,type="l",col="blue",main="Structural Breaks for Bitcoin",ylab="Price (USD)",xlab="",xaxt="n")
abline(v=Bitcoin$Date[Bitcoin_samples[1]],col="brown")
abline(v=Bitcoin$Date[Bitcoin_samples[2]],col="brown")
Date_formatted <- format(Bitcoin$Date[c(1,Bitcoin_samples,nrow(Bitcoin))],"%b-%y")
axis(1,at=c(min(Bitcoin$Date),Bitcoin$Date[Bitcoin_samples],max(Bitcoin$Date)),labels=Date_formatted,las = 2)

### Plotting the Density and Histogram for Returns

par(mfrow=c(1,2))
densityPlot(as.timeSeries(Litecoin_xts),hist = TRUE,title = FALSE,grid=FALSE)
title(main="Density Plot for Litecoin",xlab="Returns",ylab="Density/Frequency",col.main="blue",cex.main=0.9)
legend(x=0.4,y=15,legend = c("Normal","Empirical"),bty="n",cex=0.7,text.col = c("gray","brown"))

densityPlot(as.timeSeries(Bitcoin_xts),hist = TRUE,title = FALSE,grid=FALSE)
title(main="Density Plot for Bitcoin",xlab="Returns",ylab="Density/Frequency",col.main="blue",cex.main=0.9)
legend(x=0.04,y=15,legend = c("Normal","Empirical"),bty="n",cex=0.7,text.col = c("gray","brown"),col=c("gray","brown"))

### Garch Models

Litecoin_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Litecoin_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Litecoin_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Litecoin_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Litecoin_Model1 <- ugarchfit(Litecoin_Spec1,Litecoin_xts,out.sample = 500)
Litecoin_Model2 <- ugarchfit(Litecoin_Spec2,Litecoin_xts,out.sample = 500)
Litecoin_Model3 <- ugarchfit(Litecoin_Spec3,Litecoin_xts,out.sample = 500)
Litecoin_Model4 <- ugarchfit(Litecoin_Spec4,Litecoin_xts,out.sample = 500)

Bitcoin_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Bitcoin_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Bitcoin_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Bitcoin_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Bitcoin_Model1 <- ugarchfit(Bitcoin_Spec1,Bitcoin_xts,out.sample = 500)
Bitcoin_Model2 <- ugarchfit(Bitcoin_Spec2,Bitcoin_xts,out.sample = 500)
Bitcoin_Model3 <- ugarchfit(Bitcoin_Spec3,Bitcoin_xts,out.sample = 500)
Bitcoin_Model4 <- ugarchfit(Bitcoin_Spec4,Bitcoin_xts,out.sample = 500)

#### egarch models

Litecoin_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Litecoin_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Litecoin_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Litecoin_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Litecoin_Model5 <- ugarchfit(Litecoin_Spec5,Litecoin_xts,out.sample = 500)
Litecoin_Model6 <- ugarchfit(Litecoin_Spec6,Litecoin_xts,out.sample = 500)
Litecoin_Model7 <- ugarchfit(Litecoin_Spec7,Litecoin_xts,out.sample = 500)
Litecoin_Model8 <- ugarchfit(Litecoin_Spec8,Litecoin_xts,out.sample = 500)

Bitcoin_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Bitcoin_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Bitcoin_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Bitcoin_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Bitcoin_Model5 <- ugarchfit(Bitcoin_Spec5,Bitcoin_xts,out.sample = 500)
Bitcoin_Model6 <- ugarchfit(Bitcoin_Spec6,Bitcoin_xts,out.sample = 500)
Bitcoin_Model7 <- ugarchfit(Bitcoin_Spec7,Bitcoin_xts,out.sample = 500)
Bitcoin_Model8 <- ugarchfit(Bitcoin_Spec8,Bitcoin_xts,out.sample = 500)

### gjrgarch models

Litecoin_Spec9   <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Litecoin_Spec10  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Litecoin_Spec11  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Litecoin_Spec12  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Litecoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Litecoin_Model9  <- ugarchfit(Litecoin_Spec9,Litecoin_xts,out.sample = 500)
Litecoin_Model10 <- ugarchfit(Litecoin_Spec10,Litecoin_xts,out.sample = 500)
Litecoin_Model11 <- ugarchfit(Litecoin_Spec11,Litecoin_xts,out.sample = 500)
Litecoin_Model12 <- ugarchfit(Litecoin_Spec12,Litecoin_xts,out.sample = 500)

Bitcoin_Spec9  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Bitcoin_Spec10 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Bitcoin_Spec11 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Bitcoin_Spec12 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Bitcoin_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Bitcoin_Model9  <- ugarchfit(Bitcoin_Spec9,Bitcoin_xts,out.sample = 500)
Bitcoin_Model10 <- ugarchfit(Bitcoin_Spec10,Bitcoin_xts,out.sample = 500)
Bitcoin_Model11 <- ugarchfit(Bitcoin_Spec11,Bitcoin_xts,out.sample = 500)
Bitcoin_Model12 <- ugarchfit(Bitcoin_Spec12,Bitcoin_xts,out.sample = 500)

### Store the loglikelihood values, AIC and BIC ####

Llhood  <- matrix(data=0,nrow=12,ncol=2)
colnames(Llhood) <- c("Bitcoin","Litecoin")
rownames(Llhood) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")

Llhood[1,1]  <-  Bitcoin_Model1@fit$LLH
Llhood[2,1]  <-  Bitcoin_Model2@fit$LLH
Llhood[3,1]  <-  Bitcoin_Model3@fit$LLH
Llhood[4,1]  <-  Bitcoin_Model4@fit$LLH
Llhood[5,1]  <-  Bitcoin_Model5@fit$LLH
Llhood[6,1]  <-  Bitcoin_Model6@fit$LLH
Llhood[7,1]  <-  Bitcoin_Model7@fit$LLH
Llhood[8,1]  <-  Bitcoin_Model8@fit$LLH
Llhood[9,1]  <-  Bitcoin_Model9@fit$LLH
Llhood[10,1] <- Bitcoin_Model10@fit$LLH
Llhood[11,1] <- Bitcoin_Model11@fit$LLH
Llhood[12,1] <- Bitcoin_Model12@fit$LLH

Llhood[1,2]  <-  Litecoin_Model1@fit$LLH
Llhood[2,2]  <-  Litecoin_Model2@fit$LLH
Llhood[3,2]  <-  Litecoin_Model3@fit$LLH
Llhood[4,2]  <-  Litecoin_Model4@fit$LLH
Llhood[5,2]  <-  Litecoin_Model5@fit$LLH
Llhood[6,2]  <-  Litecoin_Model6@fit$LLH
Llhood[7,2]  <-  Litecoin_Model7@fit$LLH
Llhood[8,2]  <-  Litecoin_Model8@fit$LLH
Llhood[9,2]  <-  Litecoin_Model9@fit$LLH
Llhood[10,2] <- Litecoin_Model10@fit$LLH
Llhood[11,2] <- Litecoin_Model11@fit$LLH
Llhood[12,2] <- Litecoin_Model12@fit$LLH

BIC  <- matrix(data=0,nrow=12,ncol=2)
colnames(BIC) <- c("Bitcoin","Litecoin")
rownames(BIC) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")

BIC[1,1]  <-  infocriteria(Bitcoin_Model1)[2]
BIC[2,1]  <-  infocriteria(Bitcoin_Model2)[2]
BIC[3,1]  <-  infocriteria(Bitcoin_Model3)[2]
BIC[4,1]  <-  infocriteria(Bitcoin_Model4)[2]
BIC[5,1]  <-  infocriteria(Bitcoin_Model5)[2]
BIC[6,1]  <-  infocriteria(Bitcoin_Model6)[2]
BIC[7,1]  <-  infocriteria(Bitcoin_Model7)[2]
BIC[8,1]  <-  infocriteria(Bitcoin_Model8)[2]
BIC[9,1]  <-  infocriteria(Bitcoin_Model9)[2]
BIC[10,1] <- infocriteria(Bitcoin_Model10)[2]
BIC[11,1] <- infocriteria(Bitcoin_Model11)[2]
BIC[12,1] <- infocriteria(Bitcoin_Model12)[2]

BIC[1,2]  <-  infocriteria(Litecoin_Model1)[2]
BIC[2,2]  <-  infocriteria(Litecoin_Model2)[2]
BIC[3,2]  <-  infocriteria(Litecoin_Model3)[2]
BIC[4,2]  <-  infocriteria(Litecoin_Model4)[2]
BIC[5,2]  <-  infocriteria(Litecoin_Model5)[2]
BIC[6,2]  <-  infocriteria(Litecoin_Model6)[2]
BIC[7,2]  <-  infocriteria(Litecoin_Model7)[2]
BIC[8,2]  <-  infocriteria(Litecoin_Model8)[2]
BIC[9,2]  <-  infocriteria(Litecoin_Model9)[2]
BIC[10,2] <- infocriteria(Litecoin_Model10)[2]
BIC[11,2] <- infocriteria(Litecoin_Model11)[2]
BIC[12,2] <- infocriteria(Litecoin_Model12)[2]

AIC  <- matrix(data=0,nrow=12,ncol=2)
colnames(AIC) <- c("Bitcoin","Litecoin")
rownames(AIC) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")

AIC[1,1]  <-  infocriteria(Bitcoin_Model1)[1]
AIC[2,1]  <-  infocriteria(Bitcoin_Model2)[1]
AIC[3,1]  <-  infocriteria(Bitcoin_Model3)[1]
AIC[4,1]  <-  infocriteria(Bitcoin_Model4)[1]
AIC[5,1]  <-  infocriteria(Bitcoin_Model5)[1]
AIC[6,1]  <-  infocriteria(Bitcoin_Model6)[1]
AIC[7,1]  <-  infocriteria(Bitcoin_Model7)[1]
AIC[8,1]  <-  infocriteria(Bitcoin_Model8)[1]
AIC[9,1]  <-  infocriteria(Bitcoin_Model9)[1]
AIC[10,1] <- infocriteria(Bitcoin_Model10)[1]
AIC[11,1] <- infocriteria(Bitcoin_Model11)[1]
AIC[12,1] <- infocriteria(Bitcoin_Model12)[1]

AIC[1,2]  <-  infocriteria(Litecoin_Model1)[1]
AIC[2,2]  <-  infocriteria(Litecoin_Model2)[1]
AIC[3,2]  <-  infocriteria(Litecoin_Model3)[1]
AIC[4,2]  <-  infocriteria(Litecoin_Model4)[1]
AIC[5,2]  <-  infocriteria(Litecoin_Model5)[1]
AIC[6,2]  <-  infocriteria(Litecoin_Model6)[1]
AIC[7,2]  <-  infocriteria(Litecoin_Model7)[1]
AIC[8,2]  <-  infocriteria(Litecoin_Model8)[1]
AIC[9,2]  <-  infocriteria(Litecoin_Model9)[1]
AIC[10,2] <- infocriteria(Litecoin_Model10)[1]
AIC[11,2] <- infocriteria(Litecoin_Model11)[1]
AIC[12,2] <- infocriteria(Litecoin_Model12)[1]

## Conclusion: eGarch with SST and JSU are the best models #####
## Model5 for SST and Model7 for JSU, Model8 for Normal required to calculate Pearson 

### Store the coefficients of the models

Litecoin_Coef5   <- round(coef(Litecoin_Model5),5)
Litecoin_Coef6   <- round(coef(Litecoin_Model6),5)
Litecoin_Coef7   <- round(coef(Litecoin_Model7),5)
Litecoin_Coef8   <- round(coef(Litecoin_Model8),5)
Litecoin_Coef_Matrix   <- list(Litecoin_Coef5,Litecoin_Coef6,Litecoin_Coef7,Litecoin_Coef8)
names(Litecoin_Coef_Matrix) <- c("SST","SGED","JSU","Normal")

Litecoin_Robust_t5   <- round(Litecoin_Model5@fit$robust.tval,2)
Litecoin_Robust_t6   <- round(Litecoin_Model6@fit$robust.tval,2)
Litecoin_Robust_t7   <- round(Litecoin_Model7@fit$robust.tval,2)
Litecoin_Robust_t8   <- round(Litecoin_Model8@fit$robust.tval,2)

Litecoin_Robust_t   <- list(Litecoin_Robust_t5,Litecoin_Robust_t6,Litecoin_Robust_t7,Litecoin_Robust_t8)
names(Litecoin_Robust_t) <- c("SST","SGED","JSU","Normal")

Bitcoin_Coef5   <- round(coef(Bitcoin_Model5),5)
Bitcoin_Coef6   <- round(coef(Bitcoin_Model6),5)
Bitcoin_Coef7   <- round(coef(Bitcoin_Model7),5)
Bitcoin_Coef8   <- round(coef(Bitcoin_Model8),5)
Bitcoin_Coef_Matrix   <- list(Bitcoin_Coef5,Bitcoin_Coef6,Bitcoin_Coef7,Bitcoin_Coef8)
names(Bitcoin_Coef_Matrix) <- c("SST","SGED","JSU","Normal")

Bitcoin_Robust_t5   <- round(Bitcoin_Model5@fit$robust.tval,2)
Bitcoin_Robust_t6   <- round(Bitcoin_Model6@fit$robust.tval,2)
Bitcoin_Robust_t7   <- round(Bitcoin_Model7@fit$robust.tval,2)
Bitcoin_Robust_t8   <- round(Bitcoin_Model8@fit$robust.tval,2)

Bitcoin_Robust_t   <- list(Bitcoin_Robust_t5,Bitcoin_Robust_t6,Bitcoin_Robust_t7,Bitcoin_Robust_t8)
names(Bitcoin_Robust_t) <- c("SST","SGED","JSU","Normal")

### Ljung-box test statistic

Litecoin5_LQ1  <- Box.test(Litecoin_Model5@fit$residuals/Litecoin_Model5@fit$sigma,lag=1,type="Lj")
Litecoin5_LQ5  <- Box.test(Litecoin_Model5@fit$residuals/Litecoin_Model5@fit$sigma,lag=5,type="Lj")
Litecoin5_LQ10 <- Box.test(Litecoin_Model5@fit$residuals/Litecoin_Model5@fit$sigma,lag=10,type="Lj")

Litecoin6_LQ1  <- Box.test(Litecoin_Model6@fit$residuals/Litecoin_Model6@fit$sigma,lag=1,type="Lj")
Litecoin6_LQ5  <- Box.test(Litecoin_Model6@fit$residuals/Litecoin_Model6@fit$sigma,lag=5,type="Lj")
Litecoin6_LQ10 <- Box.test(Litecoin_Model6@fit$residuals/Litecoin_Model6@fit$sigma,lag=10,type="Lj")

Litecoin7_LQ1  <- Box.test(Litecoin_Model7@fit$residuals/Litecoin_Model7@fit$sigma,lag=1,type="Lj")
Litecoin7_LQ5  <- Box.test(Litecoin_Model7@fit$residuals/Litecoin_Model7@fit$sigma,lag=5,type="Lj")
Litecoin7_LQ10 <- Box.test(Litecoin_Model7@fit$residuals/Litecoin_Model7@fit$sigma,lag=10,type="Lj")

Litecoin8_LQ1  <- Box.test(Litecoin_Model8@fit$residuals/Litecoin_Model8@fit$sigma,lag=1,type="Lj")
Litecoin8_LQ5  <- Box.test(Litecoin_Model8@fit$residuals/Litecoin_Model8@fit$sigma,lag=5,type="Lj")
Litecoin8_LQ10 <- Box.test(Litecoin_Model8@fit$residuals/Litecoin_Model8@fit$sigma,lag=10,type="Lj")

Litecoin5_LQ <- c(Litecoin5_LQ1$statistic,Litecoin5_LQ5$statistic,Litecoin5_LQ10$statistic)
Litecoin6_LQ <- c(Litecoin6_LQ1$statistic,Litecoin6_LQ5$statistic,Litecoin6_LQ10$statistic)
Litecoin7_LQ <- c(Litecoin7_LQ1$statistic,Litecoin7_LQ5$statistic,Litecoin7_LQ10$statistic)
Litecoin8_LQ <- c(Litecoin8_LQ1$statistic,Litecoin8_LQ5$statistic,Litecoin8_LQ10$statistic)

Bitcoin5_LQ1  <- Box.test(Bitcoin_Model5@fit$residuals/Bitcoin_Model5@fit$sigma,lag=1,type="Lj")
Bitcoin5_LQ5  <- Box.test(Bitcoin_Model5@fit$residuals/Bitcoin_Model5@fit$sigma,lag=5,type="Lj")
Bitcoin5_LQ10 <- Box.test(Bitcoin_Model5@fit$residuals/Bitcoin_Model5@fit$sigma,lag=10,type="Lj")

Bitcoin6_LQ1  <- Box.test(Bitcoin_Model6@fit$residuals/Bitcoin_Model6@fit$sigma,lag=1,type="Lj")
Bitcoin6_LQ5  <- Box.test(Bitcoin_Model6@fit$residuals/Bitcoin_Model6@fit$sigma,lag=5,type="Lj")
Bitcoin6_LQ10 <- Box.test(Bitcoin_Model6@fit$residuals/Bitcoin_Model6@fit$sigma,lag=10,type="Lj")

Bitcoin7_LQ1  <- Box.test(Bitcoin_Model7@fit$residuals/Bitcoin_Model7@fit$sigma,lag=1,type="Lj")
Bitcoin7_LQ5  <- Box.test(Bitcoin_Model7@fit$residuals/Bitcoin_Model7@fit$sigma,lag=5,type="Lj")
Bitcoin7_LQ10 <- Box.test(Bitcoin_Model7@fit$residuals/Bitcoin_Model7@fit$sigma,lag=10,type="Lj")

Bitcoin8_LQ1  <- Box.test(Bitcoin_Model8@fit$residuals/Bitcoin_Model8@fit$sigma,lag=1,type="Lj")
Bitcoin8_LQ5  <- Box.test(Bitcoin_Model8@fit$residuals/Bitcoin_Model8@fit$sigma,lag=5,type="Lj")
Bitcoin8_LQ10 <- Box.test(Bitcoin_Model8@fit$residuals/Bitcoin_Model8@fit$sigma,lag=10,type="Lj")

Bitcoin5_LQ <- c(Bitcoin5_LQ1$statistic,Bitcoin5_LQ5$statistic,Bitcoin5_LQ10$statistic)
Bitcoin6_LQ <- c(Bitcoin6_LQ1$statistic,Bitcoin6_LQ5$statistic,Bitcoin6_LQ10$statistic)
Bitcoin7_LQ <- c(Bitcoin7_LQ1$statistic,Bitcoin7_LQ5$statistic,Bitcoin7_LQ10$statistic)
Bitcoin8_LQ <- c(Bitcoin8_LQ1$statistic,Bitcoin8_LQ5$statistic,Bitcoin8_LQ10$statistic)

Litecoin_res_stat <-   cbind(Litecoin5_LQ,Litecoin6_LQ,Litecoin7_LQ,Litecoin8_LQ)
rownames(Litecoin_res_stat) <- c("LQ1","LQ5","LQ10")
colnames(Litecoin_res_stat) <- c("SST","SGED","JSU","Normal")

Bitcoin_res_stat <-   cbind(Bitcoin5_LQ,Bitcoin6_LQ,Bitcoin7_LQ,Bitcoin8_LQ)
rownames(Bitcoin_res_stat) <- c("LQ1","LQ5","LQ10")
colnames(Bitcoin_res_stat) <- c("SST","SGED","JSU","Normal")

#### 1-day ahead Rolling Forecasts ####
## Full sample for roll i.e. last 500 points for forecast and rest for model

Litecoin_Roll5 <-  ugarchroll(Litecoin_Spec5,Litecoin_xts,forecast.length = 500,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll6 <-  ugarchroll(Litecoin_Spec6,Litecoin_xts,forecast.length = 500,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll7 <-  ugarchroll(Litecoin_Spec7,Litecoin_xts,forecast.length = 500,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll8 <-  ugarchroll(Litecoin_Spec8,Litecoin_xts,forecast.length = 500,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Bitcoin_Roll5 <-  ugarchroll(Bitcoin_Spec5,Bitcoin_xts,forecast.length = 500,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll6 <-  ugarchroll(Bitcoin_Spec6,Bitcoin_xts,forecast.length = 500,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll7 <-  ugarchroll(Bitcoin_Spec7,Bitcoin_xts,forecast.length = 500,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll8 <-  ugarchroll(Bitcoin_Spec8,Bitcoin_xts,forecast.length = 500,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


#### Store the VaR for the out sample using rolling window forecasts

Litecoin_VaR_SST   <- Litecoin_Roll5@forecast$VaR
Litecoin_VaR_SGED  <- Litecoin_Roll6@forecast$VaR
Litecoin_VaR_JSU   <- Litecoin_Roll7@forecast$VaR
Litecoin_VaR_Norm  <- Litecoin_Roll8@forecast$VaR

Bitcoin_VaR_SST   <- Bitcoin_Roll5@forecast$VaR
Bitcoin_VaR_SGED  <- Bitcoin_Roll6@forecast$VaR
Bitcoin_VaR_JSU   <- Bitcoin_Roll7@forecast$VaR
Bitcoin_VaR_Norm  <- Bitcoin_Roll8@forecast$VaR

## Significance Levels

signific_levels <- c(0.01,0.05,0.10)

### Backtest the results

Litecoin_Test_SST_1  <- BacktestVaR(Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,1],signific_levels[1])
Litecoin_Test_SST_5  <- BacktestVaR(Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,2],signific_levels[2])
Litecoin_Test_SST_10 <- BacktestVaR(Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,3],signific_levels[3])

Litecoin_Test_SGED_1  <- BacktestVaR(Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,1],signific_levels[1])
Litecoin_Test_SGED_5  <- BacktestVaR(Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,2],signific_levels[2])
Litecoin_Test_SGED_10 <- BacktestVaR(Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,3],signific_levels[3])

Litecoin_Test_JSU_1  <- BacktestVaR(Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,1],signific_levels[1])
Litecoin_Test_JSU_5  <- BacktestVaR(Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,2],signific_levels[2])
Litecoin_Test_JSU_10 <- BacktestVaR(Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,3],signific_levels[3])

Litecoin_Test_Norm_1  <- BacktestVaR(Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,1],signific_levels[1])
Litecoin_Test_Norm_5  <- BacktestVaR(Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,2],signific_levels[2])
Litecoin_Test_Norm_10 <- BacktestVaR(Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,3],signific_levels[3])

Bitcoin_Test_SST_1  <- BacktestVaR(Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,1],signific_levels[1])
Bitcoin_Test_SST_5  <- BacktestVaR(Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,2],signific_levels[2])
Bitcoin_Test_SST_10 <- BacktestVaR(Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,3],signific_levels[3])

Bitcoin_Test_SGED_1  <- BacktestVaR(Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,1],signific_levels[1])
Bitcoin_Test_SGED_5  <- BacktestVaR(Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,2],signific_levels[2])
Bitcoin_Test_SGED_10 <- BacktestVaR(Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,3],signific_levels[3])

Bitcoin_Test_JSU_1  <- BacktestVaR(Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,1],signific_levels[1])
Bitcoin_Test_JSU_5  <- BacktestVaR(Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,2],signific_levels[2])
Bitcoin_Test_JSU_10 <- BacktestVaR(Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,3],signific_levels[3])

Bitcoin_Test_Norm_1  <- BacktestVaR(Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,1],signific_levels[1])
Bitcoin_Test_Norm_5  <- BacktestVaR(Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,2],signific_levels[2])
Bitcoin_Test_Norm_10 <- BacktestVaR(Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,3],signific_levels[3])

#### Saving the test stats and pvalues for 1%, 5% and 10% VaR

## 1%
Litecoin_Backtest_teststat1 <- matrix(c(Litecoin_Test_SST_1$LRuc[1],Litecoin_Test_SST_1$LRcc[1],Litecoin_Test_SST_1$DQ$stat,
                                   Litecoin_Test_SGED_1$LRuc[1],Litecoin_Test_SGED_1$LRcc[1],Litecoin_Test_SGED_1$DQ$stat,
                                   Litecoin_Test_JSU_1$LRuc[1],Litecoin_Test_JSU_1$LRcc[1],Litecoin_Test_JSU_1$DQ$stat,
                                   Litecoin_Test_Norm_1$LRuc[1],Litecoin_Test_Norm_1$LRcc[1],Litecoin_Test_Norm_1$DQ$stat),nrow=4,byrow = TRUE)
rownames(Litecoin_Backtest_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_Backtest_teststat1)  <- c("Kupiec","Christoffersen","DQ")

Bitcoin_Backtest_teststat1  <- matrix(c(Bitcoin_Test_SST_1$LRuc[1],Bitcoin_Test_SST_1$LRcc[1],Bitcoin_Test_SST_1$DQ$stat,
                                      Bitcoin_Test_SGED_1$LRuc[1],Bitcoin_Test_SGED_1$LRcc[1],Bitcoin_Test_SGED_1$DQ$stat,
                                      Bitcoin_Test_JSU_1$LRuc[1],Bitcoin_Test_JSU_1$LRcc[1],Bitcoin_Test_JSU_1$DQ$stat,
                                      Bitcoin_Test_Norm_1$LRuc[1],Bitcoin_Test_Norm_1$LRcc[1],Bitcoin_Test_Norm_1$DQ$stat),nrow=4,byrow = TRUE)
rownames(Bitcoin_Backtest_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_Backtest_teststat1)  <- c("Kupiec","Christoffersen","DQ")

Litecoin_Backtest_pvalues1 <- matrix(c(Litecoin_Test_SST_1$LRuc[2],Litecoin_Test_SST_1$LRcc[2],Litecoin_Test_SST_1$DQ$pvalue,
                                  Litecoin_Test_SGED_1$LRuc[2],Litecoin_Test_SGED_1$LRcc[2],Litecoin_Test_SGED_1$DQ$pvalue,
                                  Litecoin_Test_JSU_1$LRuc[2],Litecoin_Test_JSU_1$LRcc[2],Litecoin_Test_JSU_1$DQ$pvalue,
                                  Litecoin_Test_Norm_1$LRuc[2],Litecoin_Test_Norm_1$LRcc[2],Litecoin_Test_Norm_1$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Litecoin_Backtest_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_Backtest_pvalues1)  <- c("Kupiec","Christoffersen","DQ")

Bitcoin_Backtest_pvalues1  <- matrix(c(Bitcoin_Test_SST_1$LRuc[2],Bitcoin_Test_SST_1$LRcc[2],Bitcoin_Test_SST_1$DQ$pvalue,
                                     Bitcoin_Test_SGED_1$LRuc[2],Bitcoin_Test_SGED_1$LRcc[2],Bitcoin_Test_SGED_1$DQ$pvalue,
                                     Bitcoin_Test_JSU_1$LRuc[2],Bitcoin_Test_JSU_1$LRcc[2],Bitcoin_Test_JSU_1$DQ$pvalue,
                                     Bitcoin_Test_Norm_1$LRuc[2],Bitcoin_Test_Norm_1$LRcc[2],Bitcoin_Test_Norm_1$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Bitcoin_Backtest_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_Backtest_pvalues1)  <- c("Kupiec","Christoffersen","DQ")

### 5%
Litecoin_Backtest_teststat5 <- matrix(c(Litecoin_Test_SST_5$LRuc[1],Litecoin_Test_SST_5$LRcc[1],Litecoin_Test_SST_5$DQ$stat,
                                   Litecoin_Test_SGED_5$LRuc[1],Litecoin_Test_SGED_5$LRcc[1],Litecoin_Test_SGED_5$DQ$stat,
                                   Litecoin_Test_JSU_5$LRuc[1],Litecoin_Test_JSU_5$LRcc[1],Litecoin_Test_JSU_5$DQ$stat,
                                   Litecoin_Test_Norm_5$LRuc[1],Litecoin_Test_Norm_5$LRcc[1],Litecoin_Test_Norm_5$DQ$stat),nrow=4,byrow = TRUE)
rownames(Litecoin_Backtest_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_Backtest_teststat5)  <- c("Kupiec","Christoffersen","DQ")

Bitcoin_Backtest_teststat5  <- matrix(c(Bitcoin_Test_SST_5$LRuc[1],Bitcoin_Test_SST_5$LRcc[1],Bitcoin_Test_SST_5$DQ$stat,
                                      Bitcoin_Test_SGED_5$LRuc[1],Bitcoin_Test_SGED_5$LRcc[1],Bitcoin_Test_SGED_5$DQ$stat,
                                      Bitcoin_Test_JSU_5$LRuc[1],Bitcoin_Test_JSU_5$LRcc[1],Bitcoin_Test_JSU_5$DQ$stat,
                                      Bitcoin_Test_Norm_5$LRuc[1],Bitcoin_Test_Norm_5$LRcc[1],Bitcoin_Test_Norm_5$DQ$stat),nrow=4,byrow = TRUE)
rownames(Bitcoin_Backtest_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_Backtest_teststat5)  <- c("Kupiec","Christoffersen","DQ")

Litecoin_Backtest_pvalues5 <- matrix(c(Litecoin_Test_SST_5$LRuc[2],Litecoin_Test_SST_5$LRcc[2],Litecoin_Test_SST_5$DQ$pvalue,
                                  Litecoin_Test_SGED_5$LRuc[2],Litecoin_Test_SGED_5$LRcc[2],Litecoin_Test_SGED_5$DQ$pvalue,
                                  Litecoin_Test_JSU_5$LRuc[2],Litecoin_Test_JSU_5$LRcc[2],Litecoin_Test_JSU_5$DQ$pvalue,
                                  Litecoin_Test_Norm_5$LRuc[2],Litecoin_Test_Norm_5$LRcc[2],Litecoin_Test_Norm_5$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Litecoin_Backtest_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_Backtest_pvalues5)  <- c("Kupiec","Christoffersen","DQ")

Bitcoin_Backtest_pvalues5  <- matrix(c(Bitcoin_Test_SST_5$LRuc[2],Bitcoin_Test_SST_5$LRcc[2],Bitcoin_Test_SST_5$DQ$pvalue,
                                     Bitcoin_Test_SGED_5$LRuc[2],Bitcoin_Test_SGED_5$LRcc[2],Bitcoin_Test_SGED_5$DQ$pvalue,
                                     Bitcoin_Test_JSU_5$LRuc[2],Bitcoin_Test_JSU_5$LRcc[2],Bitcoin_Test_JSU_5$DQ$pvalue,
                                     Bitcoin_Test_Norm_5$LRuc[2],Bitcoin_Test_Norm_5$LRcc[2],Bitcoin_Test_Norm_5$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Bitcoin_Backtest_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_Backtest_pvalues5)  <- c("Kupiec","Christoffersen","DQ")

### 10%

Litecoin_Backtest_teststat10 <- matrix(c(Litecoin_Test_SST_10$LRuc[1],Litecoin_Test_SST_10$LRcc[1],Litecoin_Test_SST_10$DQ$stat,
                                    Litecoin_Test_SGED_10$LRuc[1],Litecoin_Test_SGED_10$LRcc[1],Litecoin_Test_SGED_10$DQ$stat,
                                    Litecoin_Test_JSU_10$LRuc[1],Litecoin_Test_JSU_10$LRcc[1],Litecoin_Test_JSU_10$DQ$stat,
                                    Litecoin_Test_Norm_10$LRuc[1],Litecoin_Test_Norm_10$LRcc[1],Litecoin_Test_Norm_10$DQ$stat),nrow=4,byrow = TRUE)
rownames(Litecoin_Backtest_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_Backtest_teststat10)  <- c("Kupiec","Christoffersen","DQ")

Bitcoin_Backtest_teststat10  <- matrix(c(Bitcoin_Test_SST_10$LRuc[1],Bitcoin_Test_SST_10$LRcc[1],Bitcoin_Test_SST_10$DQ$stat,
                                       Bitcoin_Test_SGED_10$LRuc[1],Bitcoin_Test_SGED_10$LRcc[1],Bitcoin_Test_SGED_10$DQ$stat,
                                       Bitcoin_Test_JSU_10$LRuc[1],Bitcoin_Test_JSU_10$LRcc[1],Bitcoin_Test_JSU_10$DQ$stat,
                                       Bitcoin_Test_Norm_10$LRuc[1],Bitcoin_Test_Norm_10$LRcc[1],Bitcoin_Test_Norm_10$DQ$stat),nrow=4,byrow = TRUE)
rownames(Bitcoin_Backtest_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_Backtest_teststat10)  <- c("Kupiec","Christoffersen","DQ")

Litecoin_Backtest_pvalues10 <- matrix(c(Litecoin_Test_SST_10$LRuc[2],Litecoin_Test_SST_10$LRcc[2],Litecoin_Test_SST_10$DQ$pvalue,
                                   Litecoin_Test_SGED_10$LRuc[2],Litecoin_Test_SGED_10$LRcc[2],Litecoin_Test_SGED_10$DQ$pvalue,
                                   Litecoin_Test_JSU_10$LRuc[2],Litecoin_Test_JSU_10$LRcc[2],Litecoin_Test_JSU_10$DQ$pvalue,
                                   Litecoin_Test_Norm_10$LRuc[2],Litecoin_Test_Norm_10$LRcc[2],Litecoin_Test_Norm_10$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Litecoin_Backtest_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_Backtest_pvalues10)  <- c("Kupiec","Christoffersen","DQ")

Bitcoin_Backtest_pvalues10  <- matrix(c(Bitcoin_Test_SST_10$LRuc[2],Bitcoin_Test_SST_10$LRcc[2],Bitcoin_Test_SST_10$DQ$pvalue,
                                      Bitcoin_Test_SGED_10$LRuc[2],Bitcoin_Test_SGED_10$LRcc[2],Bitcoin_Test_SGED_10$DQ$pvalue,
                                      Bitcoin_Test_JSU_10$LRuc[2],Bitcoin_Test_JSU_10$LRcc[2],Bitcoin_Test_JSU_10$DQ$pvalue,
                                      Bitcoin_Test_Norm_10$LRuc[2],Bitcoin_Test_Norm_10$LRcc[2],Bitcoin_Test_Norm_10$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Bitcoin_Backtest_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_Backtest_pvalues10)  <- c("Kupiec","Christoffersen","DQ")

Litecoin_Backtest_pvalues1  <- round(Litecoin_Backtest_pvalues1,4)
Litecoin_Backtest_pvalues5  <- round(Litecoin_Backtest_pvalues5,4)
Litecoin_Backtest_pvalues10 <- round(Litecoin_Backtest_pvalues10,4)

Bitcoin_Backtest_pvalues1  <- round(Bitcoin_Backtest_pvalues1,4)
Bitcoin_Backtest_pvalues5  <- round(Bitcoin_Backtest_pvalues5,4)
Bitcoin_Backtest_pvalues10 <- round(Bitcoin_Backtest_pvalues10,4)

Litecoin_Backtest_teststat1  <- round(Litecoin_Backtest_teststat1,4)
Litecoin_Backtest_teststat5  <- round(Litecoin_Backtest_teststat5,4)
Litecoin_Backtest_teststat10 <- round(Litecoin_Backtest_teststat10,4)

Bitcoin_Backtest_teststat1  <- round(Bitcoin_Backtest_teststat1,4)
Bitcoin_Backtest_teststat5  <- round(Bitcoin_Backtest_teststat5,4)
Bitcoin_Backtest_teststat10 <- round(Bitcoin_Backtest_teststat10,4)

### Using VaRTest and VaRDuration Test individually ###

## VaR Test in R

Litecoin_Norm_VaR_Test1  <- VaRTest(0.01,Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,1])
Litecoin_Norm_VaR_Test5  <- VaRTest(0.05,Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,2])
Litecoin_Norm_VaR_Test10 <- VaRTest(0.1,Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,3])

Litecoin_JSU_VaR_Test1  <- VaRTest(0.01,Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,1])
Litecoin_JSU_VaR_Test5  <- VaRTest(0.05,Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,2])
Litecoin_JSU_VaR_Test10 <- VaRTest(0.1,Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,3])

Litecoin_SGED_VaR_Test1  <- VaRTest(0.01,Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,1])
Litecoin_SGED_VaR_Test5  <- VaRTest(0.05,Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,2])
Litecoin_SGED_VaR_Test10 <- VaRTest(0.1,Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,3])

Litecoin_SST_VaR_Test1  <- VaRTest(0.01,Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,1])
Litecoin_SST_VaR_Test5  <- VaRTest(0.05,Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,2])
Litecoin_SST_VaR_Test10 <- VaRTest(0.1,Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,3])

Bitcoin_Norm_VaR_Test1  <- VaRTest(0.01,Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,1])
Bitcoin_Norm_VaR_Test5  <- VaRTest(0.05,Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,2])
Bitcoin_Norm_VaR_Test10 <- VaRTest(0.1,Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,3])

Bitcoin_JSU_VaR_Test1  <- VaRTest(0.01,Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,1])
Bitcoin_JSU_VaR_Test5  <- VaRTest(0.05,Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,2])
Bitcoin_JSU_VaR_Test10 <- VaRTest(0.1,Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,3])

Bitcoin_SGED_VaR_Test1  <- VaRTest(0.01,Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,1])
Bitcoin_SGED_VaR_Test5  <- VaRTest(0.05,Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,2])
Bitcoin_SGED_VaR_Test10 <- VaRTest(0.1,Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,3])

Bitcoin_SST_VaR_Test1  <- VaRTest(0.01,Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,1])
Bitcoin_SST_VaR_Test5  <- VaRTest(0.05,Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,2])
Bitcoin_SST_VaR_Test10 <- VaRTest(0.1,Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,3])

## Store the pvalues at a place

Litecoin_VaR_Test_pvalues1 <- matrix(c(Litecoin_SST_VaR_Test1$uc.LRp,Litecoin_SST_VaR_Test1$cc.LRp,
                                  Litecoin_SGED_VaR_Test1$uc.LRp,Litecoin_SGED_VaR_Test1$cc.LRp,
                                  Litecoin_JSU_VaR_Test1$uc.LRp,Litecoin_JSU_VaR_Test1$cc.LRp,
                                  Litecoin_Norm_VaR_Test1$uc.LRp,Litecoin_Norm_VaR_Test1$cc.LRp),ncol=2,byrow = TRUE)

Litecoin_VaR_Test_pvalues5 <- matrix(c(Litecoin_SST_VaR_Test5$uc.LRp,Litecoin_SST_VaR_Test5$cc.LRp,
                                  Litecoin_SGED_VaR_Test5$uc.LRp,Litecoin_SGED_VaR_Test5$cc.LRp,
                                  Litecoin_JSU_VaR_Test5$uc.LRp,Litecoin_JSU_VaR_Test5$cc.LRp,
                                  Litecoin_Norm_VaR_Test5$uc.LRp,Litecoin_Norm_VaR_Test5$cc.LRp),ncol=2,byrow = TRUE)

Litecoin_VaR_Test_pvalues10 <- matrix(c(Litecoin_SST_VaR_Test10$uc.LRp,Litecoin_SST_VaR_Test10$cc.LRp,
                                   Litecoin_SGED_VaR_Test10$uc.LRp,Litecoin_SGED_VaR_Test10$cc.LRp,
                                   Litecoin_JSU_VaR_Test10$uc.LRp,Litecoin_JSU_VaR_Test10$cc.LRp,
                                   Litecoin_Norm_VaR_Test10$uc.LRp,Litecoin_Norm_VaR_Test10$cc.LRp),ncol=2,byrow = TRUE)

Bitcoin_VaR_Test_pvalues1 <- matrix(c(Bitcoin_SST_VaR_Test1$uc.LRp,Bitcoin_SST_VaR_Test1$cc.LRp,
                                    Bitcoin_SGED_VaR_Test1$uc.LRp,Bitcoin_SGED_VaR_Test1$cc.LRp,
                                    Bitcoin_JSU_VaR_Test1$uc.LRp,Bitcoin_JSU_VaR_Test1$cc.LRp,
                                    Bitcoin_Norm_VaR_Test1$uc.LRp,Bitcoin_Norm_VaR_Test1$cc.LRp),ncol=2,byrow = TRUE)

Bitcoin_VaR_Test_pvalues5 <- matrix(c(Bitcoin_SST_VaR_Test5$uc.LRp,Bitcoin_SST_VaR_Test5$cc.LRp,
                                    Bitcoin_SGED_VaR_Test5$uc.LRp,Bitcoin_SGED_VaR_Test5$cc.LRp,
                                    Bitcoin_JSU_VaR_Test5$uc.LRp,Bitcoin_JSU_VaR_Test5$cc.LRp,
                                    Bitcoin_Norm_VaR_Test5$uc.LRp,Bitcoin_Norm_VaR_Test5$cc.LRp),ncol=2,byrow = TRUE)

Bitcoin_VaR_Test_pvalues10 <- matrix(c(Bitcoin_SST_VaR_Test10$uc.LRp,Bitcoin_SST_VaR_Test10$cc.LRp,
                                     Bitcoin_SGED_VaR_Test10$uc.LRp,Bitcoin_SGED_VaR_Test10$cc.LRp,
                                     Bitcoin_JSU_VaR_Test10$uc.LRp,Bitcoin_JSU_VaR_Test10$cc.LRp,
                                     Bitcoin_Norm_VaR_Test10$uc.LRp,Bitcoin_Norm_VaR_Test10$cc.LRp),ncol=2,byrow = TRUE)

## Storing the test stats
Litecoin_VaR_Test_teststat1 <- matrix(c(Litecoin_SST_VaR_Test1$uc.LRstat,Litecoin_SST_VaR_Test1$cc.LRstat,
                                   Litecoin_SGED_VaR_Test1$uc.LRstat,Litecoin_SGED_VaR_Test1$cc.LRstat,
                                   Litecoin_JSU_VaR_Test1$uc.LRstat,Litecoin_JSU_VaR_Test1$cc.LRstat,
                                   Litecoin_Norm_VaR_Test1$uc.LRstat,Litecoin_Norm_VaR_Test1$cc.LRstat),ncol=2,byrow = TRUE)

Litecoin_VaR_Test_teststat5 <- matrix(c(Litecoin_SST_VaR_Test5$uc.LRstat,Litecoin_SST_VaR_Test5$cc.LRstat,
                                   Litecoin_SGED_VaR_Test5$uc.LRstat,Litecoin_SGED_VaR_Test5$cc.LRstat,
                                   Litecoin_JSU_VaR_Test5$uc.LRstat,Litecoin_JSU_VaR_Test5$cc.LRstat,
                                   Litecoin_Norm_VaR_Test5$uc.LRstat,Litecoin_Norm_VaR_Test5$cc.LRstat),ncol=2,byrow = TRUE)

Litecoin_VaR_Test_teststat10 <- matrix(c(Litecoin_SST_VaR_Test10$uc.LRstat,Litecoin_SST_VaR_Test10$cc.LRstat,
                                    Litecoin_SGED_VaR_Test10$uc.LRstat,Litecoin_SGED_VaR_Test10$cc.LRstat,
                                    Litecoin_JSU_VaR_Test10$uc.LRstat,Litecoin_JSU_VaR_Test10$cc.LRstat,
                                    Litecoin_Norm_VaR_Test10$uc.LRstat,Litecoin_Norm_VaR_Test10$cc.LRstat),ncol=2,byrow = TRUE)

Bitcoin_VaR_Test_teststat1 <- matrix(c(Bitcoin_SST_VaR_Test1$uc.LRstat,Bitcoin_SST_VaR_Test1$cc.LRstat,
                                     Bitcoin_SGED_VaR_Test1$uc.LRstat,Bitcoin_SGED_VaR_Test1$cc.LRstat,
                                     Bitcoin_JSU_VaR_Test1$uc.LRstat,Bitcoin_JSU_VaR_Test1$cc.LRstat,
                                     Bitcoin_Norm_VaR_Test1$uc.LRstat,Bitcoin_Norm_VaR_Test1$cc.LRstat),ncol=2,byrow = TRUE)

Bitcoin_VaR_Test_teststat5 <- matrix(c(Bitcoin_SST_VaR_Test5$uc.LRstat,Bitcoin_SST_VaR_Test5$cc.LRstat,
                                     Bitcoin_SGED_VaR_Test5$uc.LRstat,Bitcoin_SGED_VaR_Test5$cc.LRstat,
                                     Bitcoin_JSU_VaR_Test5$uc.LRstat,Bitcoin_JSU_VaR_Test5$cc.LRstat,
                                     Bitcoin_Norm_VaR_Test5$uc.LRstat,Bitcoin_Norm_VaR_Test5$cc.LRstat),ncol=2,byrow = TRUE)

Bitcoin_VaR_Test_teststat10 <- matrix(c(Bitcoin_SST_VaR_Test10$uc.LRstat,Bitcoin_SST_VaR_Test10$cc.LRstat,
                                      Bitcoin_SGED_VaR_Test10$uc.LRstat,Bitcoin_SGED_VaR_Test10$cc.LRstat,
                                      Bitcoin_JSU_VaR_Test10$uc.LRstat,Bitcoin_JSU_VaR_Test10$cc.LRstat,
                                      Bitcoin_Norm_VaR_Test10$uc.LRstat,Bitcoin_Norm_VaR_Test10$cc.LRstat),ncol=2,byrow = TRUE)
## Round them up to 4 decimals

Litecoin_VaR_Test_pvalues1  <- round(Litecoin_VaR_Test_pvalues1,4)
Litecoin_VaR_Test_pvalues5  <- round(Litecoin_VaR_Test_pvalues5,4)
Litecoin_VaR_Test_pvalues10 <- round(Litecoin_VaR_Test_pvalues10,4)

Bitcoin_VaR_Test_pvalues1  <- round(Bitcoin_VaR_Test_pvalues1,4)
Bitcoin_VaR_Test_pvalues5  <- round(Bitcoin_VaR_Test_pvalues5,4)
Bitcoin_VaR_Test_pvalues10 <- round(Bitcoin_VaR_Test_pvalues10,4)

Litecoin_VaR_Test_teststat1  <- round(Litecoin_VaR_Test_teststat1,4)
Litecoin_VaR_Test_teststat5  <- round(Litecoin_VaR_Test_teststat5,4)
Litecoin_VaR_Test_teststat10 <- round(Litecoin_VaR_Test_teststat10,4)

Bitcoin_VaR_Test_teststat1  <- round(Bitcoin_VaR_Test_teststat1,4)
Bitcoin_VaR_Test_teststat5  <- round(Bitcoin_VaR_Test_teststat5,4)
Bitcoin_VaR_Test_teststat10 <- round(Bitcoin_VaR_Test_teststat10,4)

rownames(Litecoin_VaR_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Test_pvalues1)  <- c("Kupiec","Christoffersen")
rownames(Litecoin_VaR_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Test_pvalues5)  <- c("Kupiec","Christoffersen")
rownames(Litecoin_VaR_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Test_pvalues10)  <- c("Kupiec","Christoffersen")

rownames(Bitcoin_VaR_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Test_pvalues1)  <- c("Kupiec","Christoffersen")
rownames(Bitcoin_VaR_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Test_pvalues5)  <- c("Kupiec","Christoffersen")
rownames(Bitcoin_VaR_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Test_pvalues10)  <- c("Kupiec","Christoffersen")

rownames(Litecoin_VaR_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Test_teststat1)  <- c("Kupiec","Christoffersen")
rownames(Litecoin_VaR_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Test_teststat5)  <- c("Kupiec","Christoffersen")
rownames(Litecoin_VaR_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Test_teststat10)  <- c("Kupiec","Christoffersen")

rownames(Bitcoin_VaR_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Test_teststat1)  <- c("Kupiec","Christoffersen")
rownames(Bitcoin_VaR_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Test_teststat5)  <- c("Kupiec","Christoffersen")
rownames(Bitcoin_VaR_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Test_teststat10)  <- c("Kupiec","Christoffersen")


## VaR Duration test

Litecoin_Norm_VaR_Dur_Test1  <- VaRDurTest(0.01,Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,1])
Litecoin_Norm_VaR_Dur_Test5  <- VaRDurTest(0.05,Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,2])
Litecoin_Norm_VaR_Dur_Test10 <- VaRDurTest(0.1,Litecoin_VaR_Norm[,4],Litecoin_VaR_Norm[,3])

Litecoin_JSU_VaR_Dur_Test1  <- VaRDurTest(0.01,Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,1])
Litecoin_JSU_VaR_Dur_Test5  <- VaRDurTest(0.05,Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,2])
Litecoin_JSU_VaR_Dur_Test10 <- VaRDurTest(0.1,Litecoin_VaR_JSU[,4],Litecoin_VaR_JSU[,3])

Litecoin_SGED_VaR_Dur_Test1  <- VaRDurTest(0.01,Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,1])
Litecoin_SGED_VaR_Dur_Test5  <- VaRDurTest(0.05,Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,2])
Litecoin_SGED_VaR_Dur_Test10 <- VaRDurTest(0.1,Litecoin_VaR_SGED[,4],Litecoin_VaR_SGED[,3])

Litecoin_SST_VaR_Dur_Test1  <- VaRDurTest(0.01,Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,1])
Litecoin_SST_VaR_Dur_Test5  <- VaRDurTest(0.05,Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,2])
Litecoin_SST_VaR_Dur_Test10 <- VaRDurTest(0.1,Litecoin_VaR_SST[,4],Litecoin_VaR_SST[,3])

Bitcoin_Norm_VaR_Dur_Test1  <- VaRDurTest(0.01,Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,1])
Bitcoin_Norm_VaR_Dur_Test5  <- VaRDurTest(0.05,Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,2])
Bitcoin_Norm_VaR_Dur_Test10 <- VaRDurTest(0.1,Bitcoin_VaR_Norm[,4],Bitcoin_VaR_Norm[,3])

Bitcoin_JSU_VaR_Dur_Test1  <- VaRDurTest(0.01,Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,1])
Bitcoin_JSU_VaR_Dur_Test5  <- VaRDurTest(0.05,Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,2])
Bitcoin_JSU_VaR_Dur_Test10 <- VaRDurTest(0.1,Bitcoin_VaR_JSU[,4],Bitcoin_VaR_JSU[,3])

Bitcoin_SGED_VaR_Dur_Test1  <- VaRDurTest(0.01,Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,1])
Bitcoin_SGED_VaR_Dur_Test5  <- VaRDurTest(0.05,Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,2])
Bitcoin_SGED_VaR_Dur_Test10 <- VaRDurTest(0.1,Bitcoin_VaR_SGED[,4],Bitcoin_VaR_SGED[,3])

Bitcoin_SST_VaR_Dur_Test1  <- VaRDurTest(0.01,Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,1])
Bitcoin_SST_VaR_Dur_Test5  <- VaRDurTest(0.05,Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,2])
Bitcoin_SST_VaR_Dur_Test10 <- VaRDurTest(0.1,Bitcoin_VaR_SST[,4],Bitcoin_VaR_SST[,3])

### Store the test results and p=values
Litecoin_VaR_Dur_Test_pvalues1 <- matrix(c(Litecoin_SST_VaR_Dur_Test1$LRp,Litecoin_SGED_VaR_Dur_Test1$LRp,
                                      Litecoin_JSU_VaR_Dur_Test1$LRp,Litecoin_Norm_VaR_Dur_Test1$LRp),nrow=4,byrow = TRUE)
Litecoin_VaR_Dur_Test_pvalues5 <- matrix(c(Litecoin_SST_VaR_Dur_Test5$LRp,Litecoin_SGED_VaR_Dur_Test5$LRp,
                                      Litecoin_JSU_VaR_Dur_Test5$LRp,Litecoin_Norm_VaR_Dur_Test5$LRp),nrow=4,byrow = TRUE)
Litecoin_VaR_Dur_Test_pvalues10 <- matrix(c(Litecoin_SST_VaR_Dur_Test10$LRp,Litecoin_SGED_VaR_Dur_Test10$LRp,
                                       Litecoin_JSU_VaR_Dur_Test10$LRp,Litecoin_Norm_VaR_Dur_Test10$LRp),nrow=4,byrow = TRUE)

Bitcoin_VaR_Dur_Test_pvalues1 <- matrix(c(Bitcoin_SST_VaR_Dur_Test1$LRp,Bitcoin_SGED_VaR_Dur_Test1$LRp,
                                        Bitcoin_JSU_VaR_Dur_Test1$LRp,Bitcoin_Norm_VaR_Dur_Test1$LRp),nrow=4,byrow = TRUE)
Bitcoin_VaR_Dur_Test_pvalues5 <- matrix(c(Bitcoin_SST_VaR_Dur_Test5$LRp,Bitcoin_SGED_VaR_Dur_Test5$LRp,
                                        Bitcoin_JSU_VaR_Dur_Test5$LRp,Bitcoin_Norm_VaR_Dur_Test5$LRp),nrow=4,byrow = TRUE)
Bitcoin_VaR_Dur_Test_pvalues10 <- matrix(c(Bitcoin_SST_VaR_Dur_Test10$LRp,Bitcoin_SGED_VaR_Dur_Test10$LRp,
                                         Bitcoin_JSU_VaR_Dur_Test10$LRp,Bitcoin_Norm_VaR_Dur_Test10$LRp),nrow=4,byrow = TRUE)

Litecoin_VaR_Dur_Test_teststat1 <- matrix(c(Litecoin_SST_VaR_Dur_Test1$uLL,Litecoin_SST_VaR_Dur_Test1$rLL,
                                       Litecoin_SGED_VaR_Dur_Test1$uLL,Litecoin_SGED_VaR_Dur_Test1$rLL,
                                       Litecoin_JSU_VaR_Dur_Test1$uLL,Litecoin_JSU_VaR_Dur_Test1$rLL,
                                       Litecoin_Norm_VaR_Dur_Test1$uLL,Litecoin_Norm_VaR_Dur_Test1$rLL),nrow=4,byrow = TRUE)
Litecoin_VaR_Dur_Test_teststat5 <- matrix(c(Litecoin_SST_VaR_Dur_Test5$uLL,Litecoin_SST_VaR_Dur_Test5$rLL,
                                       Litecoin_SGED_VaR_Dur_Test5$uLL,Litecoin_SGED_VaR_Dur_Test5$rLL,
                                       Litecoin_JSU_VaR_Dur_Test5$uLL,Litecoin_JSU_VaR_Dur_Test5$rLL,
                                       Litecoin_Norm_VaR_Dur_Test5$uLL,Litecoin_Norm_VaR_Dur_Test5$rLL),nrow=4,byrow = TRUE)
Litecoin_VaR_Dur_Test_teststat10 <- matrix(c(Litecoin_SST_VaR_Dur_Test10$uLL,Litecoin_SST_VaR_Dur_Test10$rLL,
                                        Litecoin_SGED_VaR_Dur_Test10$uLL,Litecoin_SGED_VaR_Dur_Test10$rLL,
                                        Litecoin_JSU_VaR_Dur_Test10$uLL,Litecoin_JSU_VaR_Dur_Test10$rLL,
                                        Litecoin_Norm_VaR_Dur_Test10$uLL,Litecoin_Norm_VaR_Dur_Test10$rLL),nrow=4,byrow = TRUE)

Bitcoin_VaR_Dur_Test_teststat1 <- matrix(c(Bitcoin_SST_VaR_Dur_Test1$uLL,Bitcoin_SST_VaR_Dur_Test1$rLL,
                                         Bitcoin_SGED_VaR_Dur_Test1$uLL,Bitcoin_SGED_VaR_Dur_Test1$rLL,
                                         Bitcoin_JSU_VaR_Dur_Test1$uLL,Bitcoin_JSU_VaR_Dur_Test1$rLL,
                                         Bitcoin_Norm_VaR_Dur_Test1$uLL,Bitcoin_Norm_VaR_Dur_Test1$rLL),nrow=4,byrow = TRUE)
Bitcoin_VaR_Dur_Test_teststat5 <- matrix(c(Bitcoin_SST_VaR_Dur_Test5$uLL,Bitcoin_SST_VaR_Dur_Test5$rLL,
                                         Bitcoin_SGED_VaR_Dur_Test5$uLL,Bitcoin_SGED_VaR_Dur_Test5$rLL,
                                         Bitcoin_JSU_VaR_Dur_Test5$uLL,Bitcoin_JSU_VaR_Dur_Test5$rLL,
                                         Bitcoin_Norm_VaR_Dur_Test5$uLL,Bitcoin_Norm_VaR_Dur_Test5$rLL),nrow=4,byrow = TRUE)
Bitcoin_VaR_Dur_Test_teststat10 <- matrix(c(Bitcoin_SST_VaR_Dur_Test10$uLL,Bitcoin_SST_VaR_Dur_Test10$rLL,
                                          Bitcoin_SGED_VaR_Dur_Test10$uLL,Bitcoin_SGED_VaR_Dur_Test10$rLL,
                                          Bitcoin_JSU_VaR_Dur_Test10$uLL,Bitcoin_JSU_VaR_Dur_Test10$rLL,
                                          Bitcoin_Norm_VaR_Dur_Test10$uLL,Bitcoin_Norm_VaR_Dur_Test10$rLL),nrow=4,byrow = TRUE)

## Round them up to 4 decimals

Litecoin_VaR_Dur_Test_pvalues1  <- round(Litecoin_VaR_Dur_Test_pvalues1,4)
Litecoin_VaR_Dur_Test_pvalues5  <- round(Litecoin_VaR_Dur_Test_pvalues5,4)
Litecoin_VaR_Dur_Test_pvalues10 <- round(Litecoin_VaR_Dur_Test_pvalues10,4)

Bitcoin_VaR_Dur_Test_pvalues1  <- round(Bitcoin_VaR_Dur_Test_pvalues1,4)
Bitcoin_VaR_Dur_Test_pvalues5  <- round(Bitcoin_VaR_Dur_Test_pvalues5,4)
Bitcoin_VaR_Dur_Test_pvalues10 <- round(Bitcoin_VaR_Dur_Test_pvalues10,4)

Litecoin_VaR_Dur_Test_teststat1  <- round(Litecoin_VaR_Dur_Test_teststat1,4)
Litecoin_VaR_Dur_Test_teststat5  <- round(Litecoin_VaR_Dur_Test_teststat5,4)
Litecoin_VaR_Dur_Test_teststat10 <- round(Litecoin_VaR_Dur_Test_teststat10,4)

Bitcoin_VaR_Dur_Test_teststat1  <- round(Bitcoin_VaR_Dur_Test_teststat1,4)
Bitcoin_VaR_Dur_Test_teststat5  <- round(Bitcoin_VaR_Dur_Test_teststat5,4)
Bitcoin_VaR_Dur_Test_teststat10 <- round(Bitcoin_VaR_Dur_Test_teststat10,4)

rownames(Litecoin_VaR_Dur_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
rownames(Litecoin_VaR_Dur_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
rownames(Litecoin_VaR_Dur_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")

rownames(Bitcoin_VaR_Dur_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
rownames(Bitcoin_VaR_Dur_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
rownames(Bitcoin_VaR_Dur_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")

rownames(Litecoin_VaR_Dur_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Dur_Test_teststat1)  <- c("Unrestricted","Restricted")
rownames(Litecoin_VaR_Dur_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Dur_Test_teststat5)  <- c("Unrestricted","Restricted")
rownames(Litecoin_VaR_Dur_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Litecoin_VaR_Dur_Test_teststat10)  <- c("Unrestricted","Restricted")

rownames(Bitcoin_VaR_Dur_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Dur_Test_teststat1)  <- c("Unrestricted","Restricted")
rownames(Bitcoin_VaR_Dur_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Dur_Test_teststat5)  <- c("Unrestricted","Restricted")
rownames(Bitcoin_VaR_Dur_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Bitcoin_VaR_Dur_Test_teststat10)  <- c("Unrestricted","Restricted")


## QL ratio and FZL ratio

QL_SGED_1  <- round(Bitcoin_Test_SGED_1$Loss$Loss/Bitcoin_Test_Norm_1$Loss$Loss,4)
QL_SGED_5  <- round(Bitcoin_Test_SGED_5$Loss$Loss/Bitcoin_Test_Norm_5$Loss$Loss,4)
QL_SGED_10 <- round(Bitcoin_Test_SGED_10$Loss$Loss/Bitcoin_Test_Norm_10$Loss$Loss,4)

QL_SST_1  <- round(Bitcoin_Test_SST_1$Loss$Loss/Bitcoin_Test_Norm_1$Loss$Loss,4)
QL_SST_5  <- round(Bitcoin_Test_SST_5$Loss$Loss/Bitcoin_Test_Norm_5$Loss$Loss,4)
QL_SST_10 <- round(Bitcoin_Test_SST_10$Loss$Loss/Bitcoin_Test_Norm_10$Loss$Loss,4)

QL_JSU_1  <- round(Bitcoin_Test_JSU_1$Loss$Loss/Bitcoin_Test_Norm_1$Loss$Loss,4)
QL_JSU_5  <- round(Bitcoin_Test_JSU_5$Loss$Loss/Bitcoin_Test_Norm_5$Loss$Loss,4)
QL_JSU_10 <- round(Bitcoin_Test_JSU_10$Loss$Loss/Bitcoin_Test_Norm_10$Loss$Loss,4)

QL <- matrix(c(QL_SST_1,QL_SST_5,QL_SST_10,QL_SGED_1,QL_SGED_5,QL_SGED_10,QL_JSU_1,QL_JSU_5,QL_JSU_10),nrow=3,byrow=TRUE)
rownames(QL) <- c("SST","SGED","JSU")
colnames(QL) <- c("1%","5%","10%")

AE_SST_1  <- round(Bitcoin_Test_SST_1$AE,4)
AE_SST_5  <- round(Bitcoin_Test_SST_5$AE,4)
AE_SST_10 <- round(Bitcoin_Test_SST_10$AE,4)

AE_SGED_1  <- round(Bitcoin_Test_SGED_1$AE,4)
AE_SGED_5  <- round(Bitcoin_Test_SGED_5$AE,4)
AE_SGED_10 <- round(Bitcoin_Test_SGED_10$AE,4)

AE_JSU_1  <- round(Bitcoin_Test_JSU_1$AE,4)
AE_JSU_5  <- round(Bitcoin_Test_JSU_5$AE,4)
AE_JSU_10 <- round(Bitcoin_Test_JSU_10$AE,4)

AE_Norm_1  <- round(Bitcoin_Test_Norm_1$AE,4)
AE_Norm_5  <- round(Bitcoin_Test_Norm_5$AE,4)
AE_Norm_10 <- round(Bitcoin_Test_Norm_10$AE,4)

AE <- matrix(c(AE_SST_1,AE_SST_5,AE_SST_10,AE_SGED_1,AE_SGED_5,AE_SGED_10,AE_JSU_1,AE_JSU_5,AE_JSU_10,AE_Norm_1,AE_Norm_5,AE_Norm_10),nrow=4,byrow=TRUE)
rownames(AE) <- c("SST","SGED","JSU","Norm")
colnames(AE) <- c("1%","5%","10%")

### Check the summary stats of the standardized residuals

Litecoin_Res8 <- Litecoin_Model8@fit$residuals/Litecoin_Model8@fit$sigma
Bitcoin_Res8   <- Bitcoin_Model8@fit$residuals/Bitcoin_Model8@fit$sigma

Litecoin_Descr8 <- basicStats(Litecoin_Res8)
Bitcoin_Descr8 <- basicStats(Bitcoin_Res8)

### Fitting distributions on the standardized residuals 

Litecoin_Dist8 <- pearsonFitML(Litecoin_Res8)
Bitcoin_Dist8   <- pearsonFitML(Bitcoin_Res8)
Litecoin_JSU_param <- JohnsonFit(Litecoin_Res8,moment ="find")
Bitcoin_JSU_param <- JohnsonFit(Bitcoin_Res8,moment ="find")

### Store the parameters for different distributions  

Litecoin_P4_Param  <- Litecoin_Dist8[2:5]
Bitcoin_P4_Param    <- Bitcoin_Dist8[2:5]

### Checking the fit using KS tests

Litecoin_test1 <- ks.test(Litecoin_Res8,rpearson(length(Litecoin_Res8),params = Litecoin_Dist8))
Bitcoin_test1 <- ks.test(Bitcoin_Res8,rpearson(length(Bitcoin_Res8),params = Bitcoin_Dist8))


### Calculate the quantiles

Litecoin_Exceed_Pearson  <- array(data=0,dim=length(signific_levels))
Bitcoin_Exceed_Pearson    <- array(data=0,dim=length(signific_levels))

for(j in 1:3)
{
  Z_Pearson_Litecoin <- qpearson(signific_levels[j],Litecoin_Dist8)
  Z_Pearson_Bitcoin   <- qpearson(signific_levels[j],Bitcoin_Dist8)
  
  mu_Litecoin   <- array(data=0,dim=length(Returns_Litecoin)-500)
  VaR_Litecoin   <- array(data=0,dim=length(Returns_Litecoin)-501)
  mu_Bitcoin   <- array(data=0,dim=length(Returns_Bitcoin)-500)
  VaR_Bitcoin   <- array(data=0,dim=length(Returns_Bitcoin)-501)
  
  for(i in 2:(length(Returns_Litecoin)-500))
  {
    mu_Litecoin[i]    <- Litecoin_Model8@fit$coef[1]+Litecoin_Model8@fit$coef[2]*Returns_Litecoin[i-1]
    VaR_Litecoin[i-1] <- mu_Litecoin[i]+Litecoin_Model8@fit$sigma[i]*Z_Pearson_Litecoin
  }
  
  for(i in 2:(length(Returns_Bitcoin)-500))
  {
    mu_Bitcoin[i] <- Bitcoin_Model8@fit$coef[1]+Bitcoin_Model8@fit$coef[2]*Returns_Bitcoin[i-1]
    VaR_Bitcoin[i-1] <- mu_Bitcoin[i]+Bitcoin_Model8@fit$sigma[i]*Z_Pearson_Bitcoin
  }
  
  Litecoin_Violations <- VaRTest(signific_levels[j],Returns_Litecoin[1:(length(VaR_Litecoin))],VaR_Litecoin)
  Bitcoin_Violations   <- VaRTest(signific_levels[j],Returns_Bitcoin[1:(length(VaR_Bitcoin))],VaR_Bitcoin)
  Litecoin_Exceed_Pearson[j] <- Litecoin_Violations$actual.exceed
  Bitcoin_Exceed_Pearson[j] <- Bitcoin_Violations$actual.exceed
  
}


#### Out of sample VaR forecast  #####

#### Using Hold out sample

Holdout_Litecoin_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Bitcoin_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Litecoin_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Bitcoin_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Litecoin_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Bitcoin_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Litecoin_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Bitcoin_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Litecoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Bitcoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Litecoin_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Bitcoin_Pearson_pvalues) <- c("1%","5%","10%")

## Store the volatility and mu forecast
#signific_levels <- c(0.1,0.05,0.01)

Litecoin_Holdout_Mu <- as.data.frame(Litecoin_Roll8)[,'Mu']
Litecoin_Holdout_Sigma <- as.data.frame(Litecoin_Roll8)[,'Sigma']
Bitcoin_Holdout_Mu <- as.data.frame(Bitcoin_Roll8)[,'Mu']
Bitcoin_Holdout_Sigma <- as.data.frame(Bitcoin_Roll8)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Litecoin <- qpearson(signific_levels[j],Litecoin_Dist8)
  Z_Pearson_Bitcoin   <- qpearson(signific_levels[j],Bitcoin_Dist8)
  
  Litecoin_Holdout_VaR <- Litecoin_Holdout_Mu + Litecoin_Holdout_Sigma*Z_Pearson_Litecoin
  Bitcoin_Holdout_VaR   <- Bitcoin_Holdout_Mu + Bitcoin_Holdout_Sigma*Z_Pearson_Bitcoin
  Litecoin_Violations <- BacktestVaR(Holdout_Returns_Litecoin[1:(length(Holdout_Returns_Litecoin)-1)],Litecoin_Holdout_VaR[-1],signific_levels[j])
  Bitcoin_Violations   <- BacktestVaR(Holdout_Returns_Bitcoin[1:(length(Holdout_Returns_Bitcoin)-1)],Bitcoin_Holdout_VaR[-1],signific_levels[j])
  Holdout_Litecoin_Pearson[j,] <- c(Litecoin_Violations$LRuc[1],Litecoin_Violations$LRcc[1],Litecoin_Violations$DQ$stat,Litecoin_Violations$AE,Litecoin_Violations$Loss$Loss)
  Holdout_Bitcoin_Pearson[j,]   <- c(Bitcoin_Violations$LRuc[1],Bitcoin_Violations$LRcc[1],Bitcoin_Violations$DQ$stat,Bitcoin_Violations$AE,Bitcoin_Violations$Loss$Loss)
  Holdout_Litecoin_Pearson_pvalues[j,] <- c(Litecoin_Violations$LRuc[2],Litecoin_Violations$LRcc[2],Litecoin_Violations$DQ$pvalue)
  Holdout_Bitcoin_Pearson_pvalues[j,] <- c(Bitcoin_Violations$LRuc[2],Bitcoin_Violations$LRcc[2],Bitcoin_Violations$DQ$pvalue)
}


#### Subsample Analysis ####

# ### JSU specification

Litecoin_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

Litecoin_Model_Sample_JSU1 <- ugarchfit(Litecoin_Spec_Sample_JSU,Litecoin_sub1,out.sample = 250)
Litecoin_Model_Sample_JSU2 <- ugarchfit(Litecoin_Spec_Sample_JSU,Litecoin_sub2,out.sample = 100)
Litecoin_Model_Sample_JSU3 <- ugarchfit(Litecoin_Spec_Sample_JSU,Litecoin_sub3,out.sample = 100)
Litecoin_Model_Sample_JSU4 <- ugarchfit(Litecoin_Spec_Sample_JSU,Litecoin_sub4,out.sample = 100)

Bitcoin_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

Bitcoin_Model_Sample_JSU1 <- ugarchfit(Bitcoin_Spec_Sample_JSU,Bitcoin_sub1,out.sample = 250)
Bitcoin_Model_Sample_JSU2 <- ugarchfit(Bitcoin_Spec_Sample_JSU,Bitcoin_sub2,out.sample = 250)
Bitcoin_Model_Sample_JSU3 <- ugarchfit(Bitcoin_Spec_Sample_JSU,Bitcoin_sub3,out.sample = 100)


## Rollover Analysis

Litecoin_Roll_JSU_Sample1 <-  ugarchroll(Litecoin_Spec_Sample_JSU,Litecoin_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_JSU_Sample2 <-  ugarchroll(Litecoin_Spec_Sample_JSU,Litecoin_sub2,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_JSU_Sample3 <-  ugarchroll(Litecoin_Spec_Sample_JSU,Litecoin_sub3,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_JSU_Sample4 <-  ugarchroll(Litecoin_Spec_Sample_JSU,Litecoin_sub4,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Bitcoin_Roll_JSU_Sample1 <-  ugarchroll(Bitcoin_Spec_Sample_JSU,Bitcoin_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll_JSU_Sample2 <-  ugarchroll(Bitcoin_Spec_Sample_JSU,Bitcoin_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll_JSU_Sample3 <-  ugarchroll(Bitcoin_Spec_Sample_JSU,Bitcoin_sub3,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### SGED specification

Litecoin_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

Litecoin_Model_Sample_SGED1 <- ugarchfit(Litecoin_Spec_Sample_SGED,Litecoin_sub1,out.sample = 250)
Litecoin_Model_Sample_SGED2 <- ugarchfit(Litecoin_Spec_Sample_SGED,Litecoin_sub2,out.sample = 100)
Litecoin_Model_Sample_SGED3 <- ugarchfit(Litecoin_Spec_Sample_SGED,Litecoin_sub3,out.sample = 100)
Litecoin_Model_Sample_SGED4 <- ugarchfit(Litecoin_Spec_Sample_SGED,Litecoin_sub4,out.sample = 100)

Bitcoin_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

Bitcoin_Model_Sample_SGED1 <- ugarchfit(Bitcoin_Spec_Sample_SGED,Bitcoin_sub1,out.sample = 250)
Bitcoin_Model_Sample_SGED2 <- ugarchfit(Bitcoin_Spec_Sample_SGED,Bitcoin_sub2,out.sample = 250)
Bitcoin_Model_Sample_SGED3 <- ugarchfit(Bitcoin_Spec_Sample_SGED,Bitcoin_sub3,out.sample = 100)


## Rollover Analysis

Litecoin_Roll_SGED_Sample1 <-  ugarchroll(Litecoin_Spec_Sample_SGED,Litecoin_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_SGED_Sample2 <-  ugarchroll(Litecoin_Spec_Sample_SGED,Litecoin_sub2,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_SGED_Sample3 <-  ugarchroll(Litecoin_Spec_Sample_SGED,Litecoin_sub3,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_SGED_Sample4 <-  ugarchroll(Litecoin_Spec_Sample_SGED,Litecoin_sub4,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Bitcoin_Roll_SGED_Sample1 <-  ugarchroll(Bitcoin_Spec_Sample_SGED,Bitcoin_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll_SGED_Sample2 <-  ugarchroll(Bitcoin_Spec_Sample_SGED,Bitcoin_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll_SGED_Sample3 <-  ugarchroll(Bitcoin_Spec_Sample_SGED,Bitcoin_sub3,forecast.length = 100, refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### SST Specification
# 
Litecoin_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
# 
Litecoin_Model_Sample_SST1 <- ugarchfit(Litecoin_Spec_Sample_SST,Litecoin_sub1,out.sample = 250)
Litecoin_Model_Sample_SST2 <- ugarchfit(Litecoin_Spec_Sample_SST,Litecoin_sub2,out.sample = 100)
Litecoin_Model_Sample_SST3 <- ugarchfit(Litecoin_Spec_Sample_SST,Litecoin_sub3,out.sample = 100)
Litecoin_Model_Sample_SST4 <- ugarchfit(Litecoin_Spec_Sample_SST,Litecoin_sub4,out.sample = 100)

Bitcoin_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")

Bitcoin_Model_Sample_SST1 <- ugarchfit(Bitcoin_Spec_Sample_SST,Bitcoin_sub1,out.sample = 250)
Bitcoin_Model_Sample_SST2 <- ugarchfit(Bitcoin_Spec_Sample_SST,Bitcoin_sub2,out.sample = 250)
Bitcoin_Model_Sample_SST3 <- ugarchfit(Bitcoin_Spec_Sample_SST,Bitcoin_sub3,out.sample = 100)


## Rollover Analysis

Litecoin_Roll_SST_Sample1 <-  ugarchroll(Litecoin_Spec_Sample_SST,Litecoin_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_SST_Sample2 <-  ugarchroll(Litecoin_Spec_Sample_SST,Litecoin_sub2,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_SST_Sample3 <-  ugarchroll(Litecoin_Spec_Sample_SST,Litecoin_sub3,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_SST_Sample4 <-  ugarchroll(Litecoin_Spec_Sample_SST,Litecoin_sub4,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Bitcoin_Roll_SST_Sample1 <-  ugarchroll(Bitcoin_Spec_Sample_SST,Bitcoin_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll_SST_Sample2 <-  ugarchroll(Bitcoin_Spec_Sample_SST,Bitcoin_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll_SST_Sample3 <-  ugarchroll(Bitcoin_Spec_Sample_SST,Bitcoin_sub3,forecast.length = 100, refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### Normal Specification

Litecoin_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Litecoin_Model_Sample_Normal1 <- ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub1,out.sample = 250)
Litecoin_Model_Sample_Normal2 <- ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub2,out.sample = 100)
Litecoin_Model_Sample_Normal3 <- ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub3,out.sample = 100)
Litecoin_Model_Sample_Normal4 <- ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub4,out.sample = 100)

Bitcoin_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Bitcoin_Model_Sample_Normal1 <- ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub1,out.sample = 250)
Bitcoin_Model_Sample_Normal2 <- ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub2,out.sample = 250)
Bitcoin_Model_Sample_Normal3 <- ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub3,out.sample = 100)

Litecoin_Roll_Normal_Sample1 <-  ugarchroll(Litecoin_Spec_Sample_Normal,Litecoin_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_Normal_Sample2 <-  ugarchroll(Litecoin_Spec_Sample_Normal,Litecoin_sub2,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_Normal_Sample3 <-  ugarchroll(Litecoin_Spec_Sample_Normal,Litecoin_sub3,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Litecoin_Roll_Normal_Sample4 <-  ugarchroll(Litecoin_Spec_Sample_Normal,Litecoin_sub4,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Bitcoin_Roll_Normal_Sample1 <-  ugarchroll(Bitcoin_Spec_Sample_Normal,Bitcoin_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll_Normal_Sample2 <-  ugarchroll(Bitcoin_Spec_Sample_Normal,Bitcoin_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Bitcoin_Roll_Normal_Sample3 <-  ugarchroll(Bitcoin_Spec_Sample_Normal,Bitcoin_sub3,forecast.length = 100,refit.every = 10,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### Sample Residuals Statistics

Litecoin_Res_Sample1 <- c(normalTest(Litecoin_Model_Sample_Normal1@fit$residuals/Litecoin_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(Litecoin_Model_Sample_Normal1@fit$residuals/Litecoin_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(Litecoin_Model_Sample_Normal1@fit$residuals/Litecoin_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
Litecoin_Res_Sample2 <- c(normalTest(Litecoin_Model_Sample_Normal2@fit$residuals/Litecoin_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(Litecoin_Model_Sample_Normal2@fit$residuals/Litecoin_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(Litecoin_Model_Sample_Normal2@fit$residuals/Litecoin_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
Litecoin_Res_Sample3 <- c(normalTest(Litecoin_Model_Sample_Normal3@fit$residuals/Litecoin_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(Litecoin_Model_Sample_Normal3@fit$residuals/Litecoin_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(Litecoin_Model_Sample_Normal3@fit$residuals/Litecoin_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
Litecoin_Res_Sample4 <- c(normalTest(Litecoin_Model_Sample_Normal4@fit$residuals/Litecoin_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(Litecoin_Model_Sample_Normal4@fit$residuals/Litecoin_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(Litecoin_Model_Sample_Normal4@fit$residuals/Litecoin_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)

Bitcoin_Res_Sample1 <- c(normalTest(Bitcoin_Model_Sample_Normal1@fit$residuals/Bitcoin_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(Bitcoin_Model_Sample_Normal1@fit$residuals/Bitcoin_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(Bitcoin_Model_Sample_Normal1@fit$residuals/Bitcoin_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
Bitcoin_Res_Sample2 <- c(normalTest(Bitcoin_Model_Sample_Normal2@fit$residuals/Bitcoin_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(Bitcoin_Model_Sample_Normal2@fit$residuals/Bitcoin_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(Bitcoin_Model_Sample_Normal2@fit$residuals/Bitcoin_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
Bitcoin_Res_Sample3 <- c(normalTest(Bitcoin_Model_Sample_Normal3@fit$residuals/Bitcoin_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(Bitcoin_Model_Sample_Normal3@fit$residuals/Bitcoin_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(Bitcoin_Model_Sample_Normal3@fit$residuals/Bitcoin_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)

Litecoin_Res_Sample_Test <- rbind(Litecoin_Res_Sample1,Litecoin_Res_Sample2,Litecoin_Res_Sample3,Litecoin_Res_Sample4)
Bitcoin_Res_Sample_Test   <- rbind(Bitcoin_Res_Sample1,Bitcoin_Res_Sample2,Bitcoin_Res_Sample3)
rownames(Litecoin_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4")
colnames(Litecoin_Res_Sample_Test) <- c("JB","SW","KS") 
rownames(Bitcoin_Res_Sample_Test) <- c("Sample1","Sample2","Sample3")
colnames(Bitcoin_Res_Sample_Test) <- c("JB","SW","KS") 

### Pearson Specification

Litecoin_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub1,out.sample = 250)@fit$residuals/ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub1,out.sample = 250)@fit$sigma)
Litecoin_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub2,out.sample = 100)@fit$residuals/ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub2,out.sample = 100)@fit$sigma)
Litecoin_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub3,out.sample = 100)@fit$residuals/ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub3,out.sample = 100)@fit$sigma)
Litecoin_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub4,out.sample = 100)@fit$residuals/ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub4,out.sample = 100)@fit$sigma)

Litecoin_Param_Pearson_Subsample <- rbind(Litecoin_Param_Pearson_Sample1,Litecoin_Param_Pearson_Sample2,Litecoin_Param_Pearson_Sample3,Litecoin_Param_Pearson_Sample4)
rownames(Litecoin_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4")

Bitcoin_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub1,out.sample = 250)@fit$residuals/ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub1,out.sample = 250)@fit$sigma)
Bitcoin_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub2,out.sample = 250)@fit$residuals/ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub2,out.sample = 250)@fit$sigma)
Bitcoin_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub3,out.sample = 100)@fit$residuals/ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub3,out.sample = 100)@fit$sigma)

Bitcoin_Param_Pearson_Subsample <- rbind(Bitcoin_Param_Pearson_Sample1,Bitcoin_Param_Pearson_Sample2,Bitcoin_Param_Pearson_Sample3)
rownames(Bitcoin_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3")

### johnson SU parameter


Litecoin_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub1,out.sample = 250)@fit$residuals/ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub1,out.sample = 250)@fit$sigma,moment = "find")
Litecoin_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub2,out.sample = 100)@fit$residuals/ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub2,out.sample = 100)@fit$sigma,moment = "find")
Litecoin_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub3,out.sample = 100)@fit$residuals/ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub3,out.sample = 100)@fit$sigma,moment = "find")
Litecoin_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub4,out.sample = 100)@fit$residuals/ugarchfit(Litecoin_Spec_Sample_Normal,Litecoin_sub4,out.sample = 100)@fit$sigma,moment = "find")

Litecoin_Param_Johnson_Subsample <- rbind(Litecoin_Param_Johnson_Sample1,Litecoin_Param_Johnson_Sample2,Litecoin_Param_Johnson_Sample3,Litecoin_Param_Johnson_Sample4)
rownames(Litecoin_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4")

Bitcoin_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub1,out.sample = 250)@fit$residuals/ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub1,out.sample = 250)@fit$sigma,moment = "find")
Bitcoin_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub2,out.sample = 250)@fit$residuals/ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub2,out.sample = 250)@fit$sigma,moment = "find")
Bitcoin_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub3,out.sample = 100)@fit$residuals/ugarchfit(Bitcoin_Spec_Sample_Normal,Bitcoin_sub3,out.sample = 100)@fit$sigma,moment = "find")

Bitcoin_Param_Johnson_Subsample <- rbind(Bitcoin_Param_Johnson_Sample1,Bitcoin_Param_Johnson_Sample2,Bitcoin_Param_Johnson_Sample3)
rownames(Bitcoin_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3")


#### Backtesting and Storing the results

Litecoin_Test_SST1_Sample1 <- BacktestVaR(Litecoin_Roll_SST_Sample1@forecast$VaR[,4],Litecoin_Roll_SST_Sample1@forecast$VaR[,1],0.01)
Litecoin_Test_SST5_Sample1 <- BacktestVaR(Litecoin_Roll_SST_Sample1@forecast$VaR[,4],Litecoin_Roll_SST_Sample1@forecast$VaR[,2],0.05)
Litecoin_Test_SST10_Sample1 <- BacktestVaR(Litecoin_Roll_SST_Sample1@forecast$VaR[,4],Litecoin_Roll_SST_Sample1@forecast$VaR[,3],0.1)

Bitcoin_Test_SST1_Sample1 <- BacktestVaR(Bitcoin_Roll_SST_Sample1@forecast$VaR[,4],Bitcoin_Roll_SST_Sample1@forecast$VaR[,1],0.01)
Bitcoin_Test_SST5_Sample1 <- BacktestVaR(Bitcoin_Roll_SST_Sample1@forecast$VaR[,4],Bitcoin_Roll_SST_Sample1@forecast$VaR[,2],0.05)
Bitcoin_Test_SST10_Sample1 <- BacktestVaR(Bitcoin_Roll_SST_Sample1@forecast$VaR[,4],Bitcoin_Roll_SST_Sample1@forecast$VaR[,3],0.1)

Litecoin_Test_SGED1_Sample1 <- BacktestVaR(Litecoin_Roll_SGED_Sample1@forecast$VaR[,4],Litecoin_Roll_SGED_Sample1@forecast$VaR[,1],0.01)
Litecoin_Test_SGED5_Sample1 <- BacktestVaR(Litecoin_Roll_SGED_Sample1@forecast$VaR[,4],Litecoin_Roll_SGED_Sample1@forecast$VaR[,2],0.05)
Litecoin_Test_SGED10_Sample1 <- BacktestVaR(Litecoin_Roll_SGED_Sample1@forecast$VaR[,4],Litecoin_Roll_SGED_Sample1@forecast$VaR[,3],0.1)

Bitcoin_Test_SGED1_Sample1 <- BacktestVaR(Bitcoin_Roll_SGED_Sample1@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample1@forecast$VaR[,1],0.01)
Bitcoin_Test_SGED5_Sample1 <- BacktestVaR(Bitcoin_Roll_SGED_Sample1@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample1@forecast$VaR[,2],0.05)
Bitcoin_Test_SGED10_Sample1 <- BacktestVaR(Bitcoin_Roll_SGED_Sample1@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample1@forecast$VaR[,3],0.1)

Litecoin_Test_JSU1_Sample1 <- BacktestVaR(Litecoin_Roll_JSU_Sample1@forecast$VaR[,4],Litecoin_Roll_JSU_Sample1@forecast$VaR[,1],0.01)
Litecoin_Test_JSU5_Sample1 <- BacktestVaR(Litecoin_Roll_JSU_Sample1@forecast$VaR[,4],Litecoin_Roll_JSU_Sample1@forecast$VaR[,2],0.05)
Litecoin_Test_JSU10_Sample1 <- BacktestVaR(Litecoin_Roll_JSU_Sample1@forecast$VaR[,4],Litecoin_Roll_JSU_Sample1@forecast$VaR[,3],0.1)

Bitcoin_Test_JSU1_Sample1 <- BacktestVaR(Bitcoin_Roll_JSU_Sample1@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample1@forecast$VaR[,1],0.01)
Bitcoin_Test_JSU5_Sample1 <- BacktestVaR(Bitcoin_Roll_JSU_Sample1@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample1@forecast$VaR[,2],0.05)
Bitcoin_Test_JSU10_Sample1 <- BacktestVaR(Bitcoin_Roll_JSU_Sample1@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample1@forecast$VaR[,3],0.1)

Litecoin_Test_Normal1_Sample1 <- BacktestVaR(Litecoin_Roll_Normal_Sample1@forecast$VaR[,4],Litecoin_Roll_Normal_Sample1@forecast$VaR[,1],0.01)
Litecoin_Test_Normal5_Sample1 <- BacktestVaR(Litecoin_Roll_Normal_Sample1@forecast$VaR[,4],Litecoin_Roll_Normal_Sample1@forecast$VaR[,2],0.05)
Litecoin_Test_Normal10_Sample1 <- BacktestVaR(Litecoin_Roll_Normal_Sample1@forecast$VaR[,4],Litecoin_Roll_Normal_Sample1@forecast$VaR[,3],0.1)

Bitcoin_Test_Normal1_Sample1 <- BacktestVaR(Bitcoin_Roll_Normal_Sample1@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample1@forecast$VaR[,1],0.01)
Bitcoin_Test_Normal5_Sample1 <- BacktestVaR(Bitcoin_Roll_Normal_Sample1@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample1@forecast$VaR[,2],0.05)
Bitcoin_Test_Normal10_Sample1 <- BacktestVaR(Bitcoin_Roll_Normal_Sample1@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample1@forecast$VaR[,3],0.1)

##### Sample 2


Litecoin_Test_SST1_Sample2 <- BacktestVaR(Litecoin_Roll_SST_Sample2@forecast$VaR[,4],Litecoin_Roll_SST_Sample2@forecast$VaR[,1],0.01)
Litecoin_Test_SST5_Sample2 <- BacktestVaR(Litecoin_Roll_SST_Sample2@forecast$VaR[,4],Litecoin_Roll_SST_Sample2@forecast$VaR[,2],0.05)
Litecoin_Test_SST10_Sample2 <- BacktestVaR(Litecoin_Roll_SST_Sample2@forecast$VaR[,4],Litecoin_Roll_SST_Sample2@forecast$VaR[,3],0.1)

Bitcoin_Test_SST1_Sample2 <- BacktestVaR(Bitcoin_Roll_SST_Sample2@forecast$VaR[,4],Bitcoin_Roll_SST_Sample2@forecast$VaR[,1],0.01)
Bitcoin_Test_SST5_Sample2 <- BacktestVaR(Bitcoin_Roll_SST_Sample2@forecast$VaR[,4],Bitcoin_Roll_SST_Sample2@forecast$VaR[,2],0.05)
Bitcoin_Test_SST10_Sample2 <- BacktestVaR(Bitcoin_Roll_SST_Sample2@forecast$VaR[,4],Bitcoin_Roll_SST_Sample2@forecast$VaR[,3],0.1)

Litecoin_Test_SGED1_Sample2 <- BacktestVaR(Litecoin_Roll_SGED_Sample2@forecast$VaR[,4],Litecoin_Roll_SGED_Sample2@forecast$VaR[,1],0.01)
Litecoin_Test_SGED5_Sample2 <- BacktestVaR(Litecoin_Roll_SGED_Sample2@forecast$VaR[,4],Litecoin_Roll_SGED_Sample2@forecast$VaR[,2],0.05)
Litecoin_Test_SGED10_Sample2 <- BacktestVaR(Litecoin_Roll_SGED_Sample2@forecast$VaR[,4],Litecoin_Roll_SGED_Sample2@forecast$VaR[,3],0.1)

Bitcoin_Test_SGED1_Sample2 <- BacktestVaR(Bitcoin_Roll_SGED_Sample2@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample2@forecast$VaR[,1],0.01)
Bitcoin_Test_SGED5_Sample2 <- BacktestVaR(Bitcoin_Roll_SGED_Sample2@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample2@forecast$VaR[,2],0.05)
Bitcoin_Test_SGED10_Sample2 <- BacktestVaR(Bitcoin_Roll_SGED_Sample2@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample2@forecast$VaR[,3],0.1)

Litecoin_Test_JSU1_Sample2 <- BacktestVaR(Litecoin_Roll_JSU_Sample2@forecast$VaR[,4],Litecoin_Roll_JSU_Sample2@forecast$VaR[,1],0.01)
Litecoin_Test_JSU5_Sample2 <- BacktestVaR(Litecoin_Roll_JSU_Sample2@forecast$VaR[,4],Litecoin_Roll_JSU_Sample2@forecast$VaR[,2],0.05)
Litecoin_Test_JSU10_Sample2 <- BacktestVaR(Litecoin_Roll_JSU_Sample2@forecast$VaR[,4],Litecoin_Roll_JSU_Sample2@forecast$VaR[,3],0.1)

Bitcoin_Test_JSU1_Sample2 <- BacktestVaR(Bitcoin_Roll_JSU_Sample2@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample2@forecast$VaR[,1],0.01)
Bitcoin_Test_JSU5_Sample2 <- BacktestVaR(Bitcoin_Roll_JSU_Sample2@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample2@forecast$VaR[,2],0.05)
Bitcoin_Test_JSU10_Sample2 <- BacktestVaR(Bitcoin_Roll_JSU_Sample2@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample2@forecast$VaR[,3],0.1)

Litecoin_Test_Normal1_Sample2 <- BacktestVaR(Litecoin_Roll_Normal_Sample2@forecast$VaR[,4],Litecoin_Roll_Normal_Sample2@forecast$VaR[,1],0.01)
Litecoin_Test_Normal5_Sample2 <- BacktestVaR(Litecoin_Roll_Normal_Sample2@forecast$VaR[,4],Litecoin_Roll_Normal_Sample2@forecast$VaR[,2],0.05)
Litecoin_Test_Normal10_Sample2 <- BacktestVaR(Litecoin_Roll_Normal_Sample2@forecast$VaR[,4],Litecoin_Roll_Normal_Sample2@forecast$VaR[,3],0.1)

Bitcoin_Test_Normal1_Sample2 <- BacktestVaR(Bitcoin_Roll_Normal_Sample2@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample2@forecast$VaR[,1],0.01)
Bitcoin_Test_Normal5_Sample2 <- BacktestVaR(Bitcoin_Roll_Normal_Sample2@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample2@forecast$VaR[,2],0.05)
Bitcoin_Test_Normal10_Sample2 <- BacktestVaR(Bitcoin_Roll_Normal_Sample2@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample2@forecast$VaR[,3],0.1)

### Sample3


Litecoin_Test_SST1_Sample3 <- BacktestVaR(Litecoin_Roll_SST_Sample3@forecast$VaR[,4],Litecoin_Roll_SST_Sample3@forecast$VaR[,1],0.01)
Litecoin_Test_SST5_Sample3 <- BacktestVaR(Litecoin_Roll_SST_Sample3@forecast$VaR[,4],Litecoin_Roll_SST_Sample3@forecast$VaR[,2],0.05)
Litecoin_Test_SST10_Sample3 <- BacktestVaR(Litecoin_Roll_SST_Sample3@forecast$VaR[,4],Litecoin_Roll_SST_Sample3@forecast$VaR[,3],0.1)

Bitcoin_Test_SST1_Sample3 <- BacktestVaR(Bitcoin_Roll_SST_Sample3@forecast$VaR[,4],Bitcoin_Roll_SST_Sample3@forecast$VaR[,1],0.01)
Bitcoin_Test_SST5_Sample3 <- BacktestVaR(Bitcoin_Roll_SST_Sample3@forecast$VaR[,4],Bitcoin_Roll_SST_Sample3@forecast$VaR[,2],0.05)
Bitcoin_Test_SST10_Sample3 <- BacktestVaR(Bitcoin_Roll_SST_Sample3@forecast$VaR[,4],Bitcoin_Roll_SST_Sample3@forecast$VaR[,3],0.1)

Litecoin_Test_SGED1_Sample3 <- BacktestVaR(Litecoin_Roll_SGED_Sample3@forecast$VaR[,4],Litecoin_Roll_SGED_Sample3@forecast$VaR[,1],0.01)
Litecoin_Test_SGED5_Sample3 <- BacktestVaR(Litecoin_Roll_SGED_Sample3@forecast$VaR[,4],Litecoin_Roll_SGED_Sample3@forecast$VaR[,2],0.05)
Litecoin_Test_SGED10_Sample3 <- BacktestVaR(Litecoin_Roll_SGED_Sample3@forecast$VaR[,4],Litecoin_Roll_SGED_Sample3@forecast$VaR[,3],0.1)

Bitcoin_Test_SGED1_Sample3 <- BacktestVaR(Bitcoin_Roll_SGED_Sample3@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample3@forecast$VaR[,1],0.01)
Bitcoin_Test_SGED5_Sample3 <- BacktestVaR(Bitcoin_Roll_SGED_Sample3@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample3@forecast$VaR[,2],0.05)
Bitcoin_Test_SGED10_Sample3 <- BacktestVaR(Bitcoin_Roll_SGED_Sample3@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample3@forecast$VaR[,3],0.1)

Litecoin_Test_JSU1_Sample3 <- BacktestVaR(Litecoin_Roll_JSU_Sample3@forecast$VaR[,4],Litecoin_Roll_JSU_Sample3@forecast$VaR[,1],0.01)
Litecoin_Test_JSU5_Sample3 <- BacktestVaR(Litecoin_Roll_JSU_Sample3@forecast$VaR[,4],Litecoin_Roll_JSU_Sample3@forecast$VaR[,2],0.05)
Litecoin_Test_JSU10_Sample3 <- BacktestVaR(Litecoin_Roll_JSU_Sample3@forecast$VaR[,4],Litecoin_Roll_JSU_Sample3@forecast$VaR[,3],0.1)

Bitcoin_Test_JSU1_Sample3 <- BacktestVaR(Bitcoin_Roll_JSU_Sample3@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample3@forecast$VaR[,1],0.01)
Bitcoin_Test_JSU5_Sample3 <- BacktestVaR(Bitcoin_Roll_JSU_Sample3@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample3@forecast$VaR[,2],0.05)
Bitcoin_Test_JSU10_Sample3 <- BacktestVaR(Bitcoin_Roll_JSU_Sample3@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample3@forecast$VaR[,3],0.1)

Litecoin_Test_Normal1_Sample3 <- BacktestVaR(Litecoin_Roll_Normal_Sample3@forecast$VaR[,4],Litecoin_Roll_Normal_Sample3@forecast$VaR[,1],0.01)
Litecoin_Test_Normal5_Sample3 <- BacktestVaR(Litecoin_Roll_Normal_Sample3@forecast$VaR[,4],Litecoin_Roll_Normal_Sample3@forecast$VaR[,2],0.05)
Litecoin_Test_Normal10_Sample3 <- BacktestVaR(Litecoin_Roll_Normal_Sample3@forecast$VaR[,4],Litecoin_Roll_Normal_Sample3@forecast$VaR[,3],0.1)

Bitcoin_Test_Normal1_Sample3 <- BacktestVaR(Bitcoin_Roll_Normal_Sample3@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample3@forecast$VaR[,1],0.01)
Bitcoin_Test_Normal5_Sample3 <- BacktestVaR(Bitcoin_Roll_Normal_Sample3@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample3@forecast$VaR[,2],0.05)
Bitcoin_Test_Normal10_Sample3 <- BacktestVaR(Bitcoin_Roll_Normal_Sample3@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample3@forecast$VaR[,3],0.1)

### Sample 4


Litecoin_Test_SST1_Sample4 <- BacktestVaR(Litecoin_Roll_SST_Sample4@forecast$VaR[,4],Litecoin_Roll_SST_Sample4@forecast$VaR[,1],0.01)
Litecoin_Test_SST5_Sample4 <- BacktestVaR(Litecoin_Roll_SST_Sample4@forecast$VaR[,4],Litecoin_Roll_SST_Sample4@forecast$VaR[,2],0.05)
Litecoin_Test_SST10_Sample4 <- BacktestVaR(Litecoin_Roll_SST_Sample4@forecast$VaR[,4],Litecoin_Roll_SST_Sample4@forecast$VaR[,3],0.1)

Litecoin_Test_SGED1_Sample4 <- BacktestVaR(Litecoin_Roll_SGED_Sample4@forecast$VaR[,4],Litecoin_Roll_SGED_Sample4@forecast$VaR[,1],0.01)
Litecoin_Test_SGED5_Sample4 <- BacktestVaR(Litecoin_Roll_SGED_Sample4@forecast$VaR[,4],Litecoin_Roll_SGED_Sample4@forecast$VaR[,2],0.05)
Litecoin_Test_SGED10_Sample4 <- BacktestVaR(Litecoin_Roll_SGED_Sample4@forecast$VaR[,4],Litecoin_Roll_SGED_Sample4@forecast$VaR[,3],0.1)

Litecoin_Test_JSU1_Sample4 <- BacktestVaR(Litecoin_Roll_JSU_Sample4@forecast$VaR[,4],Litecoin_Roll_JSU_Sample4@forecast$VaR[,1],0.01)
Litecoin_Test_JSU5_Sample4 <- BacktestVaR(Litecoin_Roll_JSU_Sample4@forecast$VaR[,4],Litecoin_Roll_JSU_Sample4@forecast$VaR[,2],0.05)
Litecoin_Test_JSU10_Sample4 <- BacktestVaR(Litecoin_Roll_JSU_Sample4@forecast$VaR[,4],Litecoin_Roll_JSU_Sample4@forecast$VaR[,3],0.1)

Litecoin_Test_Normal1_Sample4 <- BacktestVaR(Litecoin_Roll_Normal_Sample4@forecast$VaR[,4],Litecoin_Roll_Normal_Sample4@forecast$VaR[,1],0.01)
Litecoin_Test_Normal5_Sample4 <- BacktestVaR(Litecoin_Roll_Normal_Sample4@forecast$VaR[,4],Litecoin_Roll_Normal_Sample4@forecast$VaR[,2],0.05)
Litecoin_Test_Normal10_Sample4 <- BacktestVaR(Litecoin_Roll_Normal_Sample4@forecast$VaR[,4],Litecoin_Roll_Normal_Sample4@forecast$VaR[,3],0.1)


Litecoin_Sample1 <- matrix(c(Litecoin_Test_Normal1_Sample1$LRuc[1],Litecoin_Test_SST1_Sample1$LRuc[1],Litecoin_Test_SGED1_Sample1$LRuc[1],Litecoin_Test_JSU1_Sample1$LRuc[1],
                        Litecoin_Test_Normal1_Sample1$LRcc[1],Litecoin_Test_SST1_Sample1$LRcc[1],Litecoin_Test_SGED1_Sample1$LRcc[1],Litecoin_Test_JSU1_Sample1$LRcc[1],
                        Litecoin_Test_Normal1_Sample1$DQ$stat,Litecoin_Test_SST1_Sample1$DQ$stat,Litecoin_Test_SGED1_Sample1$DQ$stat,Litecoin_Test_JSU1_Sample1$DQ$stat,
                        Litecoin_Test_Normal1_Sample1$AE,Litecoin_Test_SST1_Sample1$AE,Litecoin_Test_SGED1_Sample1$AE,Litecoin_Test_JSU1_Sample1$AE,
                        (Litecoin_Test_Normal1_Sample1$Loss$Loss/Litecoin_Test_Normal1_Sample1$Loss$Loss),(Litecoin_Test_SST1_Sample1$Loss$Loss/Litecoin_Test_Normal1_Sample1$Loss$Loss),
                        (Litecoin_Test_SGED1_Sample1$Loss$Loss/Litecoin_Test_Normal1_Sample1$Loss$Loss),(Litecoin_Test_JSU1_Sample1$Loss$Loss/Litecoin_Test_Normal1_Sample1$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Litecoin_Sample1) <- c("Normal","SST","SGED","JSU")    
rownames(Litecoin_Sample1) <- c("UC","CC","DQ","AE","QL")    

Litecoin_Sample2 <- matrix(c(Litecoin_Test_Normal1_Sample2$LRuc[1],Litecoin_Test_SST1_Sample2$LRuc[1],Litecoin_Test_SGED1_Sample2$LRuc[1],Litecoin_Test_JSU1_Sample2$LRuc[1],
                        Litecoin_Test_Normal1_Sample2$LRcc[1],Litecoin_Test_SST1_Sample2$LRcc[1],Litecoin_Test_SGED1_Sample2$LRcc[1],Litecoin_Test_JSU1_Sample2$LRcc[1],
                        Litecoin_Test_Normal1_Sample2$DQ$stat,Litecoin_Test_SST1_Sample2$DQ$stat,Litecoin_Test_SGED1_Sample2$DQ$stat,Litecoin_Test_JSU1_Sample2$DQ$stat,
                        Litecoin_Test_Normal1_Sample2$AE,Litecoin_Test_SST1_Sample2$AE,Litecoin_Test_SGED1_Sample2$AE,Litecoin_Test_JSU1_Sample2$AE,
                        (Litecoin_Test_Normal1_Sample2$Loss$Loss/Litecoin_Test_Normal1_Sample2$Loss$Loss),(Litecoin_Test_SST1_Sample2$Loss$Loss/Litecoin_Test_Normal1_Sample2$Loss$Loss),
                        (Litecoin_Test_SGED1_Sample2$Loss$Loss/Litecoin_Test_Normal1_Sample2$Loss$Loss),(Litecoin_Test_JSU1_Sample2$Loss$Loss/Litecoin_Test_Normal1_Sample2$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Litecoin_Sample2) <- c("Normal","SST","SGED","JSU")    
rownames(Litecoin_Sample2) <- c("UC","CC","DQ","AE","QL")    

Litecoin_Sample3 <- matrix(c(Litecoin_Test_Normal1_Sample3$LRuc[1],Litecoin_Test_SST1_Sample3$LRuc[1],Litecoin_Test_SGED1_Sample3$LRuc[1],Litecoin_Test_JSU1_Sample3$LRuc[1],
                        Litecoin_Test_Normal1_Sample3$LRcc[1],Litecoin_Test_SST1_Sample3$LRcc[1],Litecoin_Test_SGED1_Sample3$LRcc[1],Litecoin_Test_JSU1_Sample3$LRcc[1],
                        Litecoin_Test_Normal1_Sample3$DQ$stat,Litecoin_Test_SST1_Sample3$DQ$stat,Litecoin_Test_SGED1_Sample3$DQ$stat,Litecoin_Test_JSU1_Sample3$DQ$stat,
                        Litecoin_Test_Normal1_Sample3$AE,Litecoin_Test_SST1_Sample3$AE,Litecoin_Test_SGED1_Sample3$AE,Litecoin_Test_JSU1_Sample3$AE,
                        (Litecoin_Test_Normal1_Sample3$Loss$Loss/Litecoin_Test_Normal1_Sample3$Loss$Loss),(Litecoin_Test_SST1_Sample3$Loss$Loss/Litecoin_Test_Normal1_Sample3$Loss$Loss),
                        (Litecoin_Test_SGED1_Sample3$Loss$Loss/Litecoin_Test_Normal1_Sample3$Loss$Loss),(Litecoin_Test_JSU1_Sample3$Loss$Loss/Litecoin_Test_Normal1_Sample3$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Litecoin_Sample3) <- c("Normal","SST","SGED","JSU")    
rownames(Litecoin_Sample3) <- c("UC","CC","DQ","AE","QL")    

Litecoin_Sample4 <- matrix(c(Litecoin_Test_Normal1_Sample4$LRuc[1],Litecoin_Test_SST1_Sample4$LRuc[1],Litecoin_Test_SGED1_Sample4$LRuc[1],Litecoin_Test_JSU1_Sample4$LRuc[1],
                        Litecoin_Test_Normal1_Sample4$LRcc[1],Litecoin_Test_SST1_Sample4$LRcc[1],Litecoin_Test_SGED1_Sample4$LRcc[1],Litecoin_Test_JSU1_Sample4$LRcc[1],
                        Litecoin_Test_Normal1_Sample4$DQ$stat,Litecoin_Test_SST1_Sample4$DQ$stat,Litecoin_Test_SGED1_Sample4$DQ$stat,Litecoin_Test_JSU1_Sample4$DQ$stat,
                        Litecoin_Test_Normal1_Sample4$AE,Litecoin_Test_SST1_Sample4$AE,Litecoin_Test_SGED1_Sample4$AE,Litecoin_Test_JSU1_Sample4$AE,
                        (Litecoin_Test_Normal1_Sample4$Loss$Loss/Litecoin_Test_Normal1_Sample4$Loss$Loss),(Litecoin_Test_SST1_Sample4$Loss$Loss/Litecoin_Test_Normal1_Sample4$Loss$Loss),
                        (Litecoin_Test_SGED1_Sample4$Loss$Loss/Litecoin_Test_Normal1_Sample4$Loss$Loss),(Litecoin_Test_JSU1_Sample4$Loss$Loss/Litecoin_Test_Normal1_Sample4$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Litecoin_Sample4) <- c("Normal","SST","SGED","JSU")    
rownames(Litecoin_Sample4) <- c("UC","CC","DQ","AE","QL")    


Bitcoin_Sample1 <- matrix(c(Bitcoin_Test_Normal1_Sample1$LRuc[1],Bitcoin_Test_SST1_Sample1$LRuc[1],Bitcoin_Test_SGED1_Sample1$LRuc[1],Bitcoin_Test_JSU1_Sample1$LRuc[1],
                          Bitcoin_Test_Normal1_Sample1$LRcc[1],Bitcoin_Test_SST1_Sample1$LRcc[1],Bitcoin_Test_SGED1_Sample1$LRcc[1],Bitcoin_Test_JSU1_Sample1$LRcc[1],
                          Bitcoin_Test_Normal1_Sample1$DQ$stat,Bitcoin_Test_SST1_Sample1$DQ$stat,Bitcoin_Test_SGED1_Sample1$DQ$stat,Bitcoin_Test_JSU1_Sample1$DQ$stat,
                          Bitcoin_Test_Normal1_Sample1$AE,Bitcoin_Test_SST1_Sample1$AE,Bitcoin_Test_SGED1_Sample1$AE,Bitcoin_Test_JSU1_Sample1$AE,
                          (Bitcoin_Test_Normal1_Sample1$Loss$Loss/Bitcoin_Test_Normal1_Sample1$Loss$Loss),(Bitcoin_Test_SST1_Sample1$Loss$Loss/Bitcoin_Test_Normal1_Sample1$Loss$Loss),
                          (Bitcoin_Test_SGED1_Sample1$Loss$Loss/Bitcoin_Test_Normal1_Sample1$Loss$Loss),(Bitcoin_Test_JSU1_Sample1$Loss$Loss/Bitcoin_Test_Normal1_Sample1$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Bitcoin_Sample1) <- c("Normal","SST","SGED","JSU")    
rownames(Bitcoin_Sample1) <- c("UC","CC","DQ","AE","QL")    

Bitcoin_Sample2 <- matrix(c(Bitcoin_Test_Normal1_Sample2$LRuc[1],Bitcoin_Test_SST1_Sample2$LRuc[1],Bitcoin_Test_SGED1_Sample2$LRuc[1],Bitcoin_Test_JSU1_Sample2$LRuc[1],
                          Bitcoin_Test_Normal1_Sample2$LRcc[1],Bitcoin_Test_SST1_Sample2$LRcc[1],Bitcoin_Test_SGED1_Sample2$LRcc[1],Bitcoin_Test_JSU1_Sample2$LRcc[1],
                          Bitcoin_Test_Normal1_Sample2$DQ$stat,Bitcoin_Test_SST1_Sample2$DQ$stat,Bitcoin_Test_SGED1_Sample2$DQ$stat,Bitcoin_Test_JSU1_Sample2$DQ$stat,
                          Bitcoin_Test_Normal1_Sample2$AE,Bitcoin_Test_SST1_Sample2$AE,Bitcoin_Test_SGED1_Sample2$AE,Bitcoin_Test_JSU1_Sample2$AE,
                          (Bitcoin_Test_Normal1_Sample2$Loss$Loss/Bitcoin_Test_Normal1_Sample2$Loss$Loss),(Bitcoin_Test_SST1_Sample2$Loss$Loss/Bitcoin_Test_Normal1_Sample2$Loss$Loss),
                          (Bitcoin_Test_SGED1_Sample2$Loss$Loss/Bitcoin_Test_Normal1_Sample2$Loss$Loss),(Bitcoin_Test_JSU1_Sample2$Loss$Loss/Bitcoin_Test_Normal1_Sample2$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Bitcoin_Sample2) <- c("Normal","SST","SGED","JSU")    
rownames(Bitcoin_Sample2) <- c("UC","CC","DQ","AE","QL")    

Bitcoin_Sample3 <- matrix(c(Bitcoin_Test_Normal1_Sample3$LRuc[1],Bitcoin_Test_SST1_Sample3$LRuc[1],Bitcoin_Test_SGED1_Sample3$LRuc[1],Bitcoin_Test_JSU1_Sample3$LRuc[1],
                          Bitcoin_Test_Normal1_Sample3$LRcc[1],Bitcoin_Test_SST1_Sample3$LRcc[1],Bitcoin_Test_SGED1_Sample3$LRcc[1],Bitcoin_Test_JSU1_Sample3$LRcc[1],
                          Bitcoin_Test_Normal1_Sample3$DQ$stat,Bitcoin_Test_SST1_Sample3$DQ$stat,Bitcoin_Test_SGED1_Sample3$DQ$stat,Bitcoin_Test_JSU1_Sample3$DQ$stat,
                          Bitcoin_Test_Normal1_Sample3$AE,Bitcoin_Test_SST1_Sample3$AE,Bitcoin_Test_SGED1_Sample3$AE,Bitcoin_Test_JSU1_Sample3$AE,
                          (Bitcoin_Test_Normal1_Sample3$Loss$Loss/Bitcoin_Test_Normal1_Sample3$Loss$Loss),(Bitcoin_Test_SST1_Sample3$Loss$Loss/Bitcoin_Test_Normal1_Sample3$Loss$Loss),
                          (Bitcoin_Test_SGED1_Sample3$Loss$Loss/Bitcoin_Test_Normal1_Sample3$Loss$Loss),(Bitcoin_Test_JSU1_Sample3$Loss$Loss/Bitcoin_Test_Normal1_Sample3$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Bitcoin_Sample3) <- c("Normal","SST","SGED","JSU")    
rownames(Bitcoin_Sample3) <- c("UC","CC","DQ","AE","QL")    

### Storing the p-values for subsamples 

# 1 percent
Litecoin_pvalues_1pc_Sample1 <- matrix(c(Litecoin_Test_Normal1_Sample1$LRuc[2],Litecoin_Test_SST1_Sample1$LRuc[2],Litecoin_Test_SGED1_Sample1$LRuc[2],Litecoin_Test_JSU1_Sample1$LRuc[2],
                                    Litecoin_Test_Normal1_Sample1$LRcc[2],Litecoin_Test_SST1_Sample1$LRcc[2],Litecoin_Test_SGED1_Sample1$LRcc[2],Litecoin_Test_JSU1_Sample1$LRcc[2],
                                    Litecoin_Test_Normal1_Sample1$DQ$pvalue,Litecoin_Test_SST1_Sample1$DQ$pvalue,Litecoin_Test_SGED1_Sample1$DQ$pvalue,Litecoin_Test_JSU1_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_1pc_Sample2 <- matrix(c(Litecoin_Test_Normal1_Sample2$LRuc[2],Litecoin_Test_SST1_Sample2$LRuc[2],Litecoin_Test_SGED1_Sample2$LRuc[2],Litecoin_Test_JSU1_Sample2$LRuc[2],
                                    Litecoin_Test_Normal1_Sample2$LRcc[2],Litecoin_Test_SST1_Sample2$LRcc[2],Litecoin_Test_SGED1_Sample2$LRcc[2],Litecoin_Test_JSU1_Sample2$LRcc[2],
                                    Litecoin_Test_Normal1_Sample2$DQ$pvalue,Litecoin_Test_SST1_Sample2$DQ$pvalue,Litecoin_Test_SGED1_Sample2$DQ$pvalue,Litecoin_Test_JSU1_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_1pc_Sample3 <- matrix(c(Litecoin_Test_Normal1_Sample3$LRuc[2],Litecoin_Test_SST1_Sample3$LRuc[2],Litecoin_Test_SGED1_Sample3$LRuc[2],Litecoin_Test_JSU1_Sample3$LRuc[2],
                                    Litecoin_Test_Normal1_Sample3$LRcc[2],Litecoin_Test_SST1_Sample3$LRcc[2],Litecoin_Test_SGED1_Sample3$LRcc[2],Litecoin_Test_JSU1_Sample3$LRcc[2],
                                    Litecoin_Test_Normal1_Sample3$DQ$pvalue,Litecoin_Test_SST1_Sample3$DQ$pvalue,Litecoin_Test_SGED1_Sample3$DQ$pvalue,Litecoin_Test_JSU1_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_1pc_Sample4 <- matrix(c(Litecoin_Test_Normal1_Sample4$LRuc[2],Litecoin_Test_SST1_Sample4$LRuc[2],Litecoin_Test_SGED1_Sample4$LRuc[2],Litecoin_Test_JSU1_Sample4$LRuc[2],
                                    Litecoin_Test_Normal1_Sample4$LRcc[2],Litecoin_Test_SST1_Sample4$LRcc[2],Litecoin_Test_SGED1_Sample4$LRcc[2],Litecoin_Test_JSU1_Sample4$LRcc[2],
                                    Litecoin_Test_Normal1_Sample4$DQ$pvalue,Litecoin_Test_SST1_Sample4$DQ$pvalue,Litecoin_Test_SGED1_Sample4$DQ$pvalue,Litecoin_Test_JSU1_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Bitcoin_pvalues_1pc_Sample1 <- matrix(c(Bitcoin_Test_Normal1_Sample1$LRuc[2],Bitcoin_Test_SST1_Sample1$LRuc[2],Bitcoin_Test_SGED1_Sample1$LRuc[2],Bitcoin_Test_JSU1_Sample1$LRuc[2],
                                      Bitcoin_Test_Normal1_Sample1$LRcc[2],Bitcoin_Test_SST1_Sample1$LRcc[2],Bitcoin_Test_SGED1_Sample1$LRcc[2],Bitcoin_Test_JSU1_Sample1$LRcc[2],
                                      Bitcoin_Test_Normal1_Sample1$DQ$pvalue,Bitcoin_Test_SST1_Sample1$DQ$pvalue,Bitcoin_Test_SGED1_Sample1$DQ$pvalue,Bitcoin_Test_JSU1_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Bitcoin_pvalues_1pc_Sample2 <- matrix(c(Bitcoin_Test_Normal1_Sample2$LRuc[2],Bitcoin_Test_SST1_Sample2$LRuc[2],Bitcoin_Test_SGED1_Sample2$LRuc[2],Bitcoin_Test_JSU1_Sample2$LRuc[2],
                                      Bitcoin_Test_Normal1_Sample2$LRcc[2],Bitcoin_Test_SST1_Sample2$LRcc[2],Bitcoin_Test_SGED1_Sample2$LRcc[2],Bitcoin_Test_JSU1_Sample2$LRcc[2],
                                      Bitcoin_Test_Normal1_Sample2$DQ$pvalue,Bitcoin_Test_SST1_Sample2$DQ$pvalue,Bitcoin_Test_SGED1_Sample2$DQ$pvalue,Bitcoin_Test_JSU1_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Bitcoin_pvalues_1pc_Sample3 <- matrix(c(Bitcoin_Test_Normal1_Sample3$LRuc[2],Bitcoin_Test_SST1_Sample3$LRuc[2],Bitcoin_Test_SGED1_Sample3$LRuc[2],Bitcoin_Test_JSU1_Sample3$LRuc[2],
                                      Bitcoin_Test_Normal1_Sample3$LRcc[2],Bitcoin_Test_SST1_Sample3$LRcc[2],Bitcoin_Test_SGED1_Sample3$LRcc[2],Bitcoin_Test_JSU1_Sample3$LRcc[2],
                                      Bitcoin_Test_Normal1_Sample3$DQ$pvalue,Bitcoin_Test_SST1_Sample3$DQ$pvalue,Bitcoin_Test_SGED1_Sample3$DQ$pvalue,Bitcoin_Test_JSU1_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)


# 5 percent

Litecoin_pvalues_5pc_Sample1 <- matrix(c(Litecoin_Test_Normal5_Sample1$LRuc[2],Litecoin_Test_SST5_Sample1$LRuc[2],Litecoin_Test_SGED5_Sample1$LRuc[2],Litecoin_Test_JSU5_Sample1$LRuc[2],
                                    Litecoin_Test_Normal5_Sample1$LRcc[2],Litecoin_Test_SST5_Sample1$LRcc[2],Litecoin_Test_SGED5_Sample1$LRcc[2],Litecoin_Test_JSU5_Sample1$LRcc[2],
                                    Litecoin_Test_Normal5_Sample1$DQ$pvalue,Litecoin_Test_SST5_Sample1$DQ$pvalue,Litecoin_Test_SGED5_Sample1$DQ$pvalue,Litecoin_Test_JSU5_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_5pc_Sample2 <- matrix(c(Litecoin_Test_Normal5_Sample2$LRuc[2],Litecoin_Test_SST5_Sample2$LRuc[2],Litecoin_Test_SGED5_Sample2$LRuc[2],Litecoin_Test_JSU5_Sample2$LRuc[2],
                                    Litecoin_Test_Normal5_Sample2$LRcc[2],Litecoin_Test_SST5_Sample2$LRcc[2],Litecoin_Test_SGED5_Sample2$LRcc[2],Litecoin_Test_JSU5_Sample2$LRcc[2],
                                    Litecoin_Test_Normal5_Sample2$DQ$pvalue,Litecoin_Test_SST5_Sample2$DQ$pvalue,Litecoin_Test_SGED5_Sample2$DQ$pvalue,Litecoin_Test_JSU5_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_5pc_Sample3 <- matrix(c(Litecoin_Test_Normal5_Sample3$LRuc[2],Litecoin_Test_SST5_Sample3$LRuc[2],Litecoin_Test_SGED5_Sample3$LRuc[2],Litecoin_Test_JSU5_Sample3$LRuc[2],
                                    Litecoin_Test_Normal5_Sample3$LRcc[2],Litecoin_Test_SST5_Sample3$LRcc[2],Litecoin_Test_SGED5_Sample3$LRcc[2],Litecoin_Test_JSU5_Sample3$LRcc[2],
                                    Litecoin_Test_Normal5_Sample3$DQ$pvalue,Litecoin_Test_SST5_Sample3$DQ$pvalue,Litecoin_Test_SGED5_Sample3$DQ$pvalue,Litecoin_Test_JSU5_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_5pc_Sample4 <- matrix(c(Litecoin_Test_Normal5_Sample4$LRuc[2],Litecoin_Test_SST5_Sample4$LRuc[2],Litecoin_Test_SGED5_Sample4$LRuc[2],Litecoin_Test_JSU5_Sample4$LRuc[2],
                                    Litecoin_Test_Normal5_Sample4$LRcc[2],Litecoin_Test_SST5_Sample4$LRcc[2],Litecoin_Test_SGED5_Sample4$LRcc[2],Litecoin_Test_JSU5_Sample4$LRcc[2],
                                    Litecoin_Test_Normal5_Sample4$DQ$pvalue,Litecoin_Test_SST5_Sample4$DQ$pvalue,Litecoin_Test_SGED5_Sample4$DQ$pvalue,Litecoin_Test_JSU5_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)


Bitcoin_pvalues_5pc_Sample1 <- matrix(c(Bitcoin_Test_Normal5_Sample1$LRuc[2],Bitcoin_Test_SST5_Sample1$LRuc[2],Bitcoin_Test_SGED5_Sample1$LRuc[2],Bitcoin_Test_JSU5_Sample1$LRuc[2],
                                      Bitcoin_Test_Normal5_Sample1$LRcc[2],Bitcoin_Test_SST5_Sample1$LRcc[2],Bitcoin_Test_SGED5_Sample1$LRcc[2],Bitcoin_Test_JSU5_Sample1$LRcc[2],
                                      Bitcoin_Test_Normal5_Sample1$DQ$pvalue,Bitcoin_Test_SST5_Sample1$DQ$pvalue,Bitcoin_Test_SGED5_Sample1$DQ$pvalue,Bitcoin_Test_JSU5_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Bitcoin_pvalues_5pc_Sample2 <- matrix(c(Bitcoin_Test_Normal5_Sample2$LRuc[2],Bitcoin_Test_SST5_Sample2$LRuc[2],Bitcoin_Test_SGED5_Sample2$LRuc[2],Bitcoin_Test_JSU5_Sample2$LRuc[2],
                                      Bitcoin_Test_Normal5_Sample2$LRcc[2],Bitcoin_Test_SST5_Sample2$LRcc[2],Bitcoin_Test_SGED5_Sample2$LRcc[2],Bitcoin_Test_JSU5_Sample2$LRcc[2],
                                      Bitcoin_Test_Normal5_Sample2$DQ$pvalue,Bitcoin_Test_SST5_Sample2$DQ$pvalue,Bitcoin_Test_SGED5_Sample2$DQ$pvalue,Bitcoin_Test_JSU5_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Bitcoin_pvalues_5pc_Sample3 <- matrix(c(Bitcoin_Test_Normal5_Sample3$LRuc[2],Bitcoin_Test_SST5_Sample3$LRuc[2],Bitcoin_Test_SGED5_Sample3$LRuc[2],Bitcoin_Test_JSU5_Sample3$LRuc[2],
                                      Bitcoin_Test_Normal5_Sample3$LRcc[2],Bitcoin_Test_SST5_Sample3$LRcc[2],Bitcoin_Test_SGED5_Sample3$LRcc[2],Bitcoin_Test_JSU5_Sample3$LRcc[2],
                                      Bitcoin_Test_Normal5_Sample3$DQ$pvalue,Bitcoin_Test_SST5_Sample3$DQ$pvalue,Bitcoin_Test_SGED5_Sample3$DQ$pvalue,Bitcoin_Test_JSU5_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

# 10 percent

Litecoin_pvalues_10pc_Sample1 <- matrix(c(Litecoin_Test_Normal10_Sample1$LRuc[2],Litecoin_Test_SST10_Sample1$LRuc[2],Litecoin_Test_SGED10_Sample1$LRuc[2],Litecoin_Test_JSU10_Sample1$LRuc[2],
                                     Litecoin_Test_Normal10_Sample1$LRcc[2],Litecoin_Test_SST10_Sample1$LRcc[2],Litecoin_Test_SGED10_Sample1$LRcc[2],Litecoin_Test_JSU10_Sample1$LRcc[2],
                                     Litecoin_Test_Normal10_Sample1$DQ$pvalue,Litecoin_Test_SST10_Sample1$DQ$pvalue,Litecoin_Test_SGED10_Sample1$DQ$pvalue,Litecoin_Test_JSU10_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_10pc_Sample2 <- matrix(c(Litecoin_Test_Normal10_Sample2$LRuc[2],Litecoin_Test_SST10_Sample2$LRuc[2],Litecoin_Test_SGED10_Sample2$LRuc[2],Litecoin_Test_JSU10_Sample2$LRuc[2],
                                     Litecoin_Test_Normal10_Sample2$LRcc[2],Litecoin_Test_SST10_Sample2$LRcc[2],Litecoin_Test_SGED10_Sample2$LRcc[2],Litecoin_Test_JSU10_Sample2$LRcc[2],
                                     Litecoin_Test_Normal10_Sample2$DQ$pvalue,Litecoin_Test_SST10_Sample2$DQ$pvalue,Litecoin_Test_SGED10_Sample2$DQ$pvalue,Litecoin_Test_JSU10_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_10pc_Sample3 <- matrix(c(Litecoin_Test_Normal10_Sample3$LRuc[2],Litecoin_Test_SST10_Sample3$LRuc[2],Litecoin_Test_SGED10_Sample3$LRuc[2],Litecoin_Test_JSU10_Sample3$LRuc[2],
                                     Litecoin_Test_Normal10_Sample3$LRcc[2],Litecoin_Test_SST10_Sample3$LRcc[2],Litecoin_Test_SGED10_Sample3$LRcc[2],Litecoin_Test_JSU10_Sample3$LRcc[2],
                                     Litecoin_Test_Normal10_Sample3$DQ$pvalue,Litecoin_Test_SST10_Sample3$DQ$pvalue,Litecoin_Test_SGED10_Sample3$DQ$pvalue,Litecoin_Test_JSU10_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Litecoin_pvalues_10pc_Sample4 <- matrix(c(Litecoin_Test_Normal10_Sample4$LRuc[2],Litecoin_Test_SST10_Sample4$LRuc[2],Litecoin_Test_SGED10_Sample4$LRuc[2],Litecoin_Test_JSU10_Sample4$LRuc[2],
                                     Litecoin_Test_Normal10_Sample4$LRcc[2],Litecoin_Test_SST10_Sample4$LRcc[2],Litecoin_Test_SGED10_Sample4$LRcc[2],Litecoin_Test_JSU10_Sample4$LRcc[2],
                                     Litecoin_Test_Normal10_Sample4$DQ$pvalue,Litecoin_Test_SST10_Sample4$DQ$pvalue,Litecoin_Test_SGED10_Sample4$DQ$pvalue,Litecoin_Test_JSU10_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)


Bitcoin_pvalues_10pc_Sample1 <- matrix(c(Bitcoin_Test_Normal10_Sample1$LRuc[2],Bitcoin_Test_SST10_Sample1$LRuc[2],Bitcoin_Test_SGED10_Sample1$LRuc[2],Bitcoin_Test_JSU10_Sample1$LRuc[2],
                                       Bitcoin_Test_Normal10_Sample1$LRcc[2],Bitcoin_Test_SST10_Sample1$LRcc[2],Bitcoin_Test_SGED10_Sample1$LRcc[2],Bitcoin_Test_JSU10_Sample1$LRcc[2],
                                       Bitcoin_Test_Normal10_Sample1$DQ$pvalue,Bitcoin_Test_SST10_Sample1$DQ$pvalue,Bitcoin_Test_SGED10_Sample1$DQ$pvalue,Bitcoin_Test_JSU10_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Bitcoin_pvalues_10pc_Sample2 <- matrix(c(Bitcoin_Test_Normal10_Sample2$LRuc[2],Bitcoin_Test_SST10_Sample2$LRuc[2],Bitcoin_Test_SGED10_Sample2$LRuc[2],Bitcoin_Test_JSU10_Sample2$LRuc[2],
                                       Bitcoin_Test_Normal10_Sample2$LRcc[2],Bitcoin_Test_SST10_Sample2$LRcc[2],Bitcoin_Test_SGED10_Sample2$LRcc[2],Bitcoin_Test_JSU10_Sample2$LRcc[2],
                                       Bitcoin_Test_Normal10_Sample2$DQ$pvalue,Bitcoin_Test_SST10_Sample2$DQ$pvalue,Bitcoin_Test_SGED10_Sample2$DQ$pvalue,Bitcoin_Test_JSU10_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Bitcoin_pvalues_10pc_Sample3 <- matrix(c(Bitcoin_Test_Normal10_Sample3$LRuc[2],Bitcoin_Test_SST10_Sample3$LRuc[2],Bitcoin_Test_SGED10_Sample3$LRuc[2],Bitcoin_Test_JSU10_Sample3$LRuc[2],
                                       Bitcoin_Test_Normal10_Sample3$LRcc[2],Bitcoin_Test_SST10_Sample3$LRcc[2],Bitcoin_Test_SGED10_Sample3$LRcc[2],Bitcoin_Test_JSU10_Sample3$LRcc[2],
                                       Bitcoin_Test_Normal10_Sample3$DQ$pvalue,Bitcoin_Test_SST10_Sample3$DQ$pvalue,Bitcoin_Test_SGED10_Sample3$DQ$pvalue,Bitcoin_Test_JSU10_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)


Litecoin_pvalues_1pc   <- rbind(Litecoin_pvalues_1pc_Sample1,Litecoin_pvalues_1pc_Sample2,Litecoin_pvalues_1pc_Sample3,Litecoin_pvalues_1pc_Sample4)
Litecoin_pvalues_5pc   <- rbind(Litecoin_pvalues_5pc_Sample1,Litecoin_pvalues_5pc_Sample2,Litecoin_pvalues_5pc_Sample3,Litecoin_pvalues_5pc_Sample4)
Litecoin_pvalues_10pc  <- rbind(Litecoin_pvalues_10pc_Sample1,Litecoin_pvalues_10pc_Sample2,Litecoin_pvalues_10pc_Sample3,Litecoin_pvalues_10pc_Sample4)

Bitcoin_pvalues_1pc   <- rbind(Bitcoin_pvalues_1pc_Sample1,Bitcoin_pvalues_1pc_Sample2,Bitcoin_pvalues_1pc_Sample3)
Bitcoin_pvalues_5pc   <- rbind(Bitcoin_pvalues_5pc_Sample1,Bitcoin_pvalues_5pc_Sample2,Bitcoin_pvalues_5pc_Sample3)
Bitcoin_pvalues_10pc  <- rbind(Bitcoin_pvalues_10pc_Sample1,Bitcoin_pvalues_10pc_Sample2,Bitcoin_pvalues_10pc_Sample3)

rownames(Litecoin_pvalues_1pc) <- c(rep(c("UC","CC","DQ"),4))
colnames(Litecoin_pvalues_1pc) <- c("Normal","SST","SGED","JSU")

rownames(Litecoin_pvalues_5pc) <- c(rep(c("UC","CC","DQ"),4))
colnames(Litecoin_pvalues_5pc) <- c("Normal","SST","SGED","JSU")

rownames(Litecoin_pvalues_10pc) <- c(rep(c("UC","CC","DQ"),4))
colnames(Litecoin_pvalues_10pc) <- c("Normal","SST","SGED","JSU")

rownames(Bitcoin_pvalues_1pc) <- c(rep(c("UC","CC","DQ"),3))
colnames(Bitcoin_pvalues_1pc) <- c("Normal","SST","SGED","JSU")

rownames(Bitcoin_pvalues_5pc) <- c(rep(c("UC","CC","DQ"),3))
colnames(Bitcoin_pvalues_5pc) <- c("Normal","SST","SGED","JSU")

rownames(Bitcoin_pvalues_10pc) <- c(rep(c("UC","CC","DQ"),3))
colnames(Bitcoin_pvalues_10pc) <- c("Normal","SST","SGED","JSU")

#### Repeating for Pearson Distribution  #####

#signific_levels <- c(0.1,0.05,0.01)
## Store the volatility and mu forecast for Sample 1

Holdout_Sample1_Litecoin_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample1_Bitcoin_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample1_Litecoin_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample1_Bitcoin_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample1_Litecoin_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample1_Bitcoin_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample1_Litecoin_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample1_Bitcoin_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample1_Litecoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample1_Bitcoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample1_Litecoin_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample1_Bitcoin_Pearson_pvalues) <- c("1%","5%","10%")

Sample1_Litecoin_Holdout_Mu <- as.data.frame(Litecoin_Roll_Normal_Sample1)[,'Mu']
Sample1_Litecoin_Holdout_Sigma <- as.data.frame(Litecoin_Roll_Normal_Sample1)[,'Sigma']
Sample1_Bitcoin_Holdout_Mu <- as.data.frame(Bitcoin_Roll_Normal_Sample1)[,'Mu']
Sample1_Bitcoin_Holdout_Sigma <- as.data.frame(Bitcoin_Roll_Normal_Sample1)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample1_Litecoin <- qpearson(signific_levels[j],Litecoin_Param_Pearson_Sample1)
  Z_Pearson_Sample1_Bitcoin   <- qpearson(signific_levels[j],Bitcoin_Param_Pearson_Sample1)
  
  Sample1_Litecoin_Holdout_VaR <- Sample1_Litecoin_Holdout_Mu + Sample1_Litecoin_Holdout_Sigma*Z_Pearson_Sample1_Litecoin
  Sample1_Bitcoin_Holdout_VaR   <- Sample1_Bitcoin_Holdout_Mu + Sample1_Bitcoin_Holdout_Sigma*Z_Pearson_Sample1_Bitcoin
  Sample1_Litecoin_Violations <- BacktestVaR(Returns_Litecoin[(Litecoin_samples[1]-250+1):(Litecoin_samples[1]-1)],Sample1_Litecoin_Holdout_VaR[-1],signific_levels[j])
  Sample1_Bitcoin_Violations   <- BacktestVaR(Returns_Bitcoin[(Bitcoin_samples[1]-250+1):(Bitcoin_samples[1]-1)],Sample1_Bitcoin_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample1_Litecoin_Pearson[j,] <- c(Sample1_Litecoin_Violations$LRuc[1],Sample1_Litecoin_Violations$LRcc[1],Sample1_Litecoin_Violations$DQ$stat,Sample1_Litecoin_Violations$AE,Sample1_Litecoin_Violations$Loss$Loss)
  Holdout_Sample1_Bitcoin_Pearson[j,]   <- c(Sample1_Bitcoin_Violations$LRuc[1],Sample1_Bitcoin_Violations$LRcc[1],Sample1_Bitcoin_Violations$DQ$stat,Sample1_Bitcoin_Violations$AE,Sample1_Bitcoin_Violations$Loss$Loss)
  Holdout_Sample1_Litecoin_Pearson_pvalues[j,] <- c(Sample1_Litecoin_Violations$LRuc[2],Sample1_Litecoin_Violations$LRcc[2],Sample1_Litecoin_Violations$DQ$pvalue)
  Holdout_Sample1_Bitcoin_Pearson_pvalues[j,] <- c(Sample1_Bitcoin_Violations$LRuc[2],Sample1_Bitcoin_Violations$LRcc[2],Sample1_Bitcoin_Violations$DQ$pvalue)
}


## Store the volatility and mu forecast for Sample 2

Holdout_Sample2_Litecoin_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample2_Bitcoin_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample2_Litecoin_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample2_Bitcoin_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample2_Litecoin_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample2_Bitcoin_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample2_Litecoin_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample2_Bitcoin_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample2_Litecoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample2_Bitcoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample2_Litecoin_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample2_Bitcoin_Pearson_pvalues) <- c("1%","5%","10%")

Sample2_Litecoin_Holdout_Mu <- as.data.frame(Litecoin_Roll_Normal_Sample2)[,'Mu']
Sample2_Litecoin_Holdout_Sigma <- as.data.frame(Litecoin_Roll_Normal_Sample2)[,'Sigma']
Sample2_Bitcoin_Holdout_Mu <- as.data.frame(Bitcoin_Roll_Normal_Sample2)[,'Mu']
Sample2_Bitcoin_Holdout_Sigma <- as.data.frame(Bitcoin_Roll_Normal_Sample2)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample2_Litecoin <- qpearson(signific_levels[j],Litecoin_Param_Pearson_Sample2)
  Z_Pearson_Sample2_Bitcoin   <- qpearson(signific_levels[j],Bitcoin_Param_Pearson_Sample2)
  
  Sample2_Litecoin_Holdout_VaR <- Sample2_Litecoin_Holdout_Mu + Sample2_Litecoin_Holdout_Sigma*Z_Pearson_Sample2_Litecoin
  Sample2_Bitcoin_Holdout_VaR   <- Sample2_Bitcoin_Holdout_Mu + Sample2_Bitcoin_Holdout_Sigma*Z_Pearson_Sample2_Bitcoin
  Sample2_Litecoin_Violations <- BacktestVaR(Returns_Litecoin[(Litecoin_samples[2]-100+1):(Litecoin_samples[2]-1)],Sample2_Litecoin_Holdout_VaR[-1],signific_levels[j])
  Sample2_Bitcoin_Violations   <- BacktestVaR(Returns_Bitcoin[(Bitcoin_samples[2]-250+1):(Bitcoin_samples[2]-1)],Sample2_Bitcoin_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample2_Litecoin_Pearson[j,] <- c(Sample2_Litecoin_Violations$LRuc[1],Sample2_Litecoin_Violations$LRcc[1],Sample2_Litecoin_Violations$DQ$stat,Sample2_Litecoin_Violations$AE,Sample2_Litecoin_Violations$Loss$Loss)
  Holdout_Sample2_Bitcoin_Pearson[j,]   <- c(Sample2_Bitcoin_Violations$LRuc[1],Sample2_Bitcoin_Violations$LRcc[1],Sample2_Bitcoin_Violations$DQ$stat,Sample2_Bitcoin_Violations$AE,Sample2_Bitcoin_Violations$Loss$Loss)
  Holdout_Sample2_Litecoin_Pearson_pvalues[j,] <- c(Sample2_Litecoin_Violations$LRuc[2],Sample2_Litecoin_Violations$LRcc[2],Sample2_Litecoin_Violations$DQ$pvalue)
  Holdout_Sample2_Bitcoin_Pearson_pvalues[j,] <- c(Sample2_Bitcoin_Violations$LRuc[2],Sample2_Bitcoin_Violations$LRcc[2],Sample2_Bitcoin_Violations$DQ$pvalue)
}

## Store the volatility and mu forecast for Sample 3

Holdout_Sample3_Litecoin_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample3_Bitcoin_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample3_Litecoin_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample3_Bitcoin_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample3_Litecoin_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample3_Bitcoin_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample3_Litecoin_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample3_Bitcoin_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample3_Litecoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample3_Bitcoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample3_Litecoin_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample3_Bitcoin_Pearson_pvalues) <- c("1%","5%","10%")

Sample3_Litecoin_Holdout_Mu <- as.data.frame(Litecoin_Roll_Normal_Sample3)[,'Mu']
Sample3_Litecoin_Holdout_Sigma <- as.data.frame(Litecoin_Roll_Normal_Sample3)[,'Sigma']
Sample3_Bitcoin_Holdout_Mu <- as.data.frame(Bitcoin_Roll_Normal_Sample3)[,'Mu']
Sample3_Bitcoin_Holdout_Sigma <- as.data.frame(Bitcoin_Roll_Normal_Sample3)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample3_Litecoin <- qpearson(signific_levels[j],Litecoin_Param_Pearson_Sample3)
  Z_Pearson_Sample3_Bitcoin   <- qpearson(signific_levels[j],Bitcoin_Param_Pearson_Sample3)
  
  Sample3_Litecoin_Holdout_VaR <- Sample3_Litecoin_Holdout_Mu + Sample3_Litecoin_Holdout_Sigma*Z_Pearson_Sample3_Litecoin
  Sample3_Bitcoin_Holdout_VaR   <- Sample3_Bitcoin_Holdout_Mu + Sample3_Bitcoin_Holdout_Sigma*Z_Pearson_Sample3_Bitcoin
  Sample3_Litecoin_Violations <- BacktestVaR(Returns_Litecoin[(Litecoin_samples[3]-100+1):(Litecoin_samples[3]-1)],Sample3_Litecoin_Holdout_VaR[-1],signific_levels[j])
  Sample3_Bitcoin_Violations   <- BacktestVaR(Returns_Bitcoin[(length(Returns_Bitcoin)-100+1):(length(Returns_Bitcoin)-1)],Sample3_Bitcoin_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample3_Litecoin_Pearson[j,] <- c(Sample3_Litecoin_Violations$LRuc[1],Sample3_Litecoin_Violations$LRcc[1],Sample3_Litecoin_Violations$DQ$stat,Sample3_Litecoin_Violations$AE,Sample3_Litecoin_Violations$Loss$Loss)
  Holdout_Sample3_Bitcoin_Pearson[j,]   <- c(Sample3_Bitcoin_Violations$LRuc[1],Sample3_Bitcoin_Violations$LRcc[1],Sample3_Bitcoin_Violations$DQ$stat,Sample3_Bitcoin_Violations$AE,Sample3_Bitcoin_Violations$Loss$Loss)
  Holdout_Sample3_Litecoin_Pearson_pvalues[j,] <- c(Sample3_Litecoin_Violations$LRuc[2],Sample3_Litecoin_Violations$LRcc[2],Sample3_Litecoin_Violations$DQ$pvalue)
  Holdout_Sample3_Bitcoin_Pearson_pvalues[j,] <- c(Sample3_Bitcoin_Violations$LRuc[2],Sample3_Bitcoin_Violations$LRcc[2],Sample3_Bitcoin_Violations$DQ$pvalue)
}

## Store the volatility and mu forecast for Sample 4

Holdout_Sample4_Litecoin_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample4_Litecoin_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample4_Litecoin_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample4_Litecoin_Pearson_pvalues) <- c("UC","CC","DQ")
rownames(Holdout_Sample4_Litecoin_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample4_Litecoin_Pearson_pvalues) <- c("1%","5%","10%")

Sample4_Litecoin_Holdout_Mu <- as.data.frame(Litecoin_Roll_Normal_Sample4)[,'Mu']
Sample4_Litecoin_Holdout_Sigma <- as.data.frame(Litecoin_Roll_Normal_Sample4)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample4_Litecoin <- qpearson(signific_levels[j],Litecoin_Param_Pearson_Sample4)
  
  Sample4_Litecoin_Holdout_VaR <- Sample4_Litecoin_Holdout_Mu + Sample4_Litecoin_Holdout_Sigma*Z_Pearson_Sample4_Litecoin
  Sample4_Litecoin_Violations  <- BacktestVaR(Returns_Litecoin[(length(Returns_Litecoin)-100+1):(length(Returns_Litecoin)-1)],Sample4_Litecoin_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample4_Litecoin_Pearson[j,] <- c(Sample4_Litecoin_Violations$LRuc[1],Sample4_Litecoin_Violations$LRcc[1],Sample4_Litecoin_Violations$DQ$stat,Sample4_Litecoin_Violations$AE,Sample4_Litecoin_Violations$Loss$Loss)
  Holdout_Sample4_Litecoin_Pearson_pvalues[j,] <- c(Sample4_Litecoin_Violations$LRuc[2],Sample4_Litecoin_Violations$LRcc[2],Sample4_Litecoin_Violations$DQ$pvalue)
  }


### Plotting the VaR for subsamples

Litecoin_Sample1_Dates <- format(Litecoin$Date[(Litecoin_samples[1]-249):Litecoin_samples[1]],"%b-%y")
Litecoin_Sample2_Dates <- format(Litecoin$Date[(Litecoin_samples[2]-49):Litecoin_samples[2]],"%d-%b-%y")
Litecoin_Sample3_Dates <- format(Litecoin$Date[(Litecoin_samples[3]-49):Litecoin_samples[3]],"%d-%b-%y")
Litecoin_Sample4_Dates <- format(Litecoin$Date[(length(Litecoin_xts)-49):length(Litecoin_xts)],"%d-%b-%y")

Bitcoin_Sample1_Dates <- format(Bitcoin$Date[(Bitcoin_samples[1]-249):Bitcoin_samples[1]],"%b-%y")
Bitcoin_Sample2_Dates <- format(Bitcoin$Date[(Bitcoin_samples[2]-249):Bitcoin_samples[2]],"%b-%y")
Bitcoin_Sample3_Dates <- format(Bitcoin$Date[(length(Bitcoin_xts)-49):length(Bitcoin_xts)],"%d-%b-%y")

par(mfrow=c(2,3))
plot(Bitcoin_Roll_Normal_Sample1@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample1",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Bitcoin_Roll_Normal_Sample1@forecast$VaR[,3],type="l",col="yellow")
lines(Bitcoin_Roll_SST_Sample1@forecast$VaR[,3],type="l",col="green")
lines(Bitcoin_Roll_SGED_Sample1@forecast$VaR[,3],type="l",col="gray")
lines(Bitcoin_Roll_JSU_Sample1@forecast$VaR[,3],type="l",col="brown")
lines(Sample1_Bitcoin_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Bitcoin_Sample1_Dates[c(1,50,100,150,200,250)],las=2)


plot(Bitcoin_Roll_JSU_Sample2@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample2",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Bitcoin_Roll_Normal_Sample2@forecast$VaR[,3],type="l",col="yellow")
lines(Bitcoin_Roll_SST_Sample2@forecast$VaR[,3],type="l",col="green")
lines(Bitcoin_Roll_SGED_Sample2@forecast$VaR[,3],type="l",col="gray")
lines(Bitcoin_Roll_JSU_Sample2@forecast$VaR[,3],type="l",col="brown")
lines(Sample2_Bitcoin_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Bitcoin_Sample2_Dates[c(1,50,100,150,200,250)],las=2)

plot(Bitcoin_Roll_JSU_Sample3@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample3",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Bitcoin_Roll_Normal_Sample3@forecast$VaR[,3],type="l",col="yellow")
lines(Bitcoin_Roll_SST_Sample3@forecast$VaR[,3],type="l",col="green")
lines(Bitcoin_Roll_SGED_Sample3@forecast$VaR[,3],type="l",col="gray")
lines(Bitcoin_Roll_JSU_Sample3@forecast$VaR[,3],type="l",col="brown")
lines(Sample3_Bitcoin_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,10,20,30,40,50),labels=Bitcoin_Sample3_Dates[c(1,10,20,30,40,50)],las=2)



### For Litecoin

par(mfrow=c(2,3))
plot(Litecoin_Roll_Normal_Sample1@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample1",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Litecoin_Roll_Normal_Sample1@forecast$VaR[,3],type="l",col="yellow")
lines(Litecoin_Roll_SST_Sample1@forecast$VaR[,3],type="l",col="green")
lines(Litecoin_Roll_SGED_Sample1@forecast$VaR[,3],type="l",col="gray")
lines(Litecoin_Roll_JSU_Sample1@forecast$VaR[,3],type="l",col="brown")
lines(Sample1_Litecoin_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Litecoin_Sample1_Dates[c(1,50,100,150,200,250)],las=2)


plot(Litecoin_Roll_JSU_Sample2@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample2",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Litecoin_Roll_Normal_Sample2@forecast$VaR[,3],type="l",col="yellow")
lines(Litecoin_Roll_SST_Sample2@forecast$VaR[,3],type="l",col="green")
lines(Litecoin_Roll_SGED_Sample2@forecast$VaR[,3],type="l",col="gray")
lines(Litecoin_Roll_JSU_Sample2@forecast$VaR[,3],type="l",col="brown")
lines(Sample2_Litecoin_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,10,20,30,40,50),labels=Litecoin_Sample2_Dates[c(1,10,20,30,40,50)],las=2)

plot(Litecoin_Roll_JSU_Sample3@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample3",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Litecoin_Roll_Normal_Sample3@forecast$VaR[,3],type="l",col="yellow")
lines(Litecoin_Roll_SST_Sample3@forecast$VaR[,3],type="l",col="green")
lines(Litecoin_Roll_SGED_Sample3@forecast$VaR[,3],type="l",col="gray")
lines(Litecoin_Roll_JSU_Sample3@forecast$VaR[,3],type="l",col="brown")
lines(Sample3_Litecoin_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,10,20,30,40,50),labels=Litecoin_Sample3_Dates[c(1,10,20,30,40,50)],las=2)

plot(Litecoin_Roll_JSU_Sample4@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample4",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Litecoin_Roll_Normal_Sample4@forecast$VaR[,3],type="l",col="yellow")
lines(Litecoin_Roll_SST_Sample4@forecast$VaR[,3],type="l",col="green")
lines(Litecoin_Roll_SGED_Sample4@forecast$VaR[,3],type="l",col="gray")
lines(Litecoin_Roll_JSU_Sample4@forecast$VaR[,3],type="l",col="brown")
lines(Sample4_Litecoin_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,10,20,30,40,50),labels=Litecoin_Sample4_Dates[c(1,10,20,30,40,50)],las=2)


### Incorporating Structural Breaks in a different way using dummies###
# 

par(mfrow=c(1,2))
plot(Litecoin_Roll8@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample for Litecoin",ylab="Returns/ VaR at 1%",xlab="",xaxt="n")
lines(Litecoin_Roll8@forecast$VaR[,1],type="l",col="yellow")
lines(Litecoin_Roll5@forecast$VaR[,1],type="l",col="green")
lines(Litecoin_Roll6@forecast$VaR[,1],type="l",col="gray")
lines(Litecoin_Roll7@forecast$VaR[,1],type="l",col="brown")
legend(x=1,y=0.34,legend=c("Returns","Normal","SST","SGED","JSU"),col=c("blue","yellow","green","gray","brown"),cex=0.8,lty=1,bty="n",y.intersp = 0.5,ncol=3,x.intersp = 0.49)
axis(1,at=c(0,100,200,300,400,450),labels=format(Litecoin_Date[c(2500,2600,2700,2800,2900,2950)],"%b-%y"),las=2)

### Keep the plot width as 800 pixels and save as .tiff extension

plot(Bitcoin_Roll8@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample for Bitcoin",ylab="Returns/ VaR at 1%",xlab="",xaxt="n")
lines(Bitcoin_Roll8@forecast$VaR[,1],type="l",col="yellow")
lines(Bitcoin_Roll5@forecast$VaR[,1],type="l",col="green")
lines(Bitcoin_Roll6@forecast$VaR[,1],type="l",col="gray")
lines(Bitcoin_Roll7@forecast$VaR[,1],type="l",col="brown")
lines(Bitcoin_Holdout_VaR[-1],type="l",col="black")
legend(x=1,y=0.14,legend=c("Returns","Normal","SST","SGED","JSU","Pearson"),col=c("blue","yellow","green","gray","brown","black"),cex=0.8,lty=1,bty="n",y.intersp = 0.7,x.intersp = 0.49,ncol=2)
axis(1,at=c(0,200,400,600,800,950),labels=format(Bitcoin_Date[c(4001,4200,4400,4600,4800,5000)],"%b-%y"),las=2)

#### Checked on Monday, 12th July, 2021 ###
### Storing Pearson subsample results ###

Litecoin_subsample_1pc <- rbind(Holdout_Sample1_Litecoin_Pearson[1,],Holdout_Sample2_Litecoin_Pearson[1,],Holdout_Sample3_Litecoin_Pearson[1,],Holdout_Sample4_Litecoin_Pearson[1,],Holdout_Sample5_Litecoin_Pearson[1,])
Bitcoin_subsample_1pc   <- rbind(Holdout_Sample1_Bitcoin_Pearson[1,],Holdout_Sample2_Bitcoin_Pearson[1,],Holdout_Sample3_Bitcoin_Pearson[1,],Holdout_Sample4_Bitcoin_Pearson[1,],Holdout_Sample5_Bitcoin_Pearson[1,])

Litecoin_subsample_5pc <- rbind(Holdout_Sample1_Litecoin_Pearson[2,],Holdout_Sample2_Litecoin_Pearson[2,],Holdout_Sample3_Litecoin_Pearson[2,],Holdout_Sample4_Litecoin_Pearson[2,],Holdout_Sample5_Litecoin_Pearson[2,])
Bitcoin_subsample_5pc   <- rbind(Holdout_Sample1_Bitcoin_Pearson[2,],Holdout_Sample2_Bitcoin_Pearson[2,],Holdout_Sample3_Bitcoin_Pearson[2,],Holdout_Sample4_Bitcoin_Pearson[2,],Holdout_Sample5_Bitcoin_Pearson[2,])

Litecoin_subsample_10pc <- rbind(Holdout_Sample1_Litecoin_Pearson[3,],Holdout_Sample2_Litecoin_Pearson[3,],Holdout_Sample3_Litecoin_Pearson[3,],Holdout_Sample4_Litecoin_Pearson[3,],Holdout_Sample5_Litecoin_Pearson[3,])
Bitcoin_subsample_10pc   <- rbind(Holdout_Sample1_Bitcoin_Pearson[3,],Holdout_Sample2_Bitcoin_Pearson[3,],Holdout_Sample3_Bitcoin_Pearson[3,],Holdout_Sample4_Bitcoin_Pearson[3,],Holdout_Sample5_Bitcoin_Pearson[3,])

Litecoin_subsample_pvalues_1pc <- rbind(Holdout_Sample1_Litecoin_Pearson_pvalues[1,],Holdout_Sample2_Litecoin_Pearson_pvalues[1,],Holdout_Sample3_Litecoin_Pearson_pvalues[1,],Holdout_Sample4_Litecoin_Pearson_pvalues[1,],Holdout_Sample5_Litecoin_Pearson_pvalues[1,])
Bitcoin_subsample_pvalues_1pc   <- rbind(Holdout_Sample1_Bitcoin_Pearson_pvalues[1,],Holdout_Sample2_Bitcoin_Pearson_pvalues[1,],Holdout_Sample3_Bitcoin_Pearson_pvalues[1,],Holdout_Sample4_Bitcoin_Pearson_pvalues[1,],Holdout_Sample5_Bitcoin_Pearson_pvalues[1,])

Litecoin_subsample_pvalues_5pc <- rbind(Holdout_Sample1_Litecoin_Pearson_pvalues[2,],Holdout_Sample2_Litecoin_Pearson_pvalues[2,],Holdout_Sample3_Litecoin_Pearson_pvalues[2,],Holdout_Sample4_Litecoin_Pearson_pvalues[2,],Holdout_Sample5_Litecoin_Pearson_pvalues[2,])
Bitcoin_subsample_pvalues_5pc   <- rbind(Holdout_Sample1_Bitcoin_Pearson_pvalues[2,],Holdout_Sample2_Bitcoin_Pearson_pvalues[2,],Holdout_Sample3_Bitcoin_Pearson_pvalues[2,],Holdout_Sample4_Bitcoin_Pearson_pvalues[2,],Holdout_Sample5_Bitcoin_Pearson_pvalues[2,])

Litecoin_subsample_pvalues_10pc <- rbind(Holdout_Sample1_Litecoin_Pearson_pvalues[3,],Holdout_Sample2_Litecoin_Pearson_pvalues[3,],Holdout_Sample3_Litecoin_Pearson_pvalues[3,],Holdout_Sample4_Litecoin_Pearson_pvalues[3,],Holdout_Sample5_Litecoin_Pearson_pvalues[3,])
Bitcoin_subsample_pvalues_10pc   <- rbind(Holdout_Sample1_Bitcoin_Pearson_pvalues[3,],Holdout_Sample2_Bitcoin_Pearson_pvalues[3,],Holdout_Sample3_Bitcoin_Pearson_pvalues[3,],Holdout_Sample4_Bitcoin_Pearson_pvalues[3,],Holdout_Sample5_Bitcoin_Pearson_pvalues[3,])

#### VaRDuration Test for subsamples

Litecoin_Dur_Normal_Sample1 <- VaRDurTest(0.10,Litecoin_Roll_Normal_Sample1@forecast$VaR[,4],Litecoin_Roll_Normal_Sample1@forecast$VaR[,3])
Litecoin_Dur_Normal_Sample2 <- VaRDurTest(0.10,Litecoin_Roll_Normal_Sample2@forecast$VaR[,4],Litecoin_Roll_Normal_Sample2@forecast$VaR[,3])
Litecoin_Dur_Normal_Sample3 <- VaRDurTest(0.10,Litecoin_Roll_Normal_Sample3@forecast$VaR[,4],Litecoin_Roll_Normal_Sample3@forecast$VaR[,3])
Litecoin_Dur_Normal_Sample4 <- VaRDurTest(0.10,Litecoin_Roll_Normal_Sample4@forecast$VaR[,4],Litecoin_Roll_Normal_Sample4@forecast$VaR[,3])
Litecoin_Dur_Normal_Sample5 <- VaRDurTest(0.10,Litecoin_Roll_Normal_Sample5@forecast$VaR[,4],Litecoin_Roll_Normal_Sample5@forecast$VaR[,3])

Bitcoin_Dur_Normal_Sample1 <- VaRDurTest(0.10,Bitcoin_Roll_Normal_Sample1@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample1@forecast$VaR[,3])
Bitcoin_Dur_Normal_Sample2 <- VaRDurTest(0.10,Bitcoin_Roll_Normal_Sample2@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample2@forecast$VaR[,3])
Bitcoin_Dur_Normal_Sample3 <- VaRDurTest(0.10,Bitcoin_Roll_Normal_Sample3@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample3@forecast$VaR[,3])
Bitcoin_Dur_Normal_Sample4 <- VaRDurTest(0.10,Bitcoin_Roll_Normal_Sample4@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample4@forecast$VaR[,3])
Bitcoin_Dur_Normal_Sample5 <- VaRDurTest(0.10,Bitcoin_Roll_Normal_Sample5@forecast$VaR[,4],Bitcoin_Roll_Normal_Sample5@forecast$VaR[,3])

Litecoin_Dur_SST_Sample1 <- VaRDurTest(0.10,Litecoin_Roll_SST_Sample1@forecast$VaR[,4],Litecoin_Roll_SST_Sample1@forecast$VaR[,3])
Litecoin_Dur_SST_Sample2 <- VaRDurTest(0.10,Litecoin_Roll_SST_Sample2@forecast$VaR[,4],Litecoin_Roll_SST_Sample2@forecast$VaR[,3])
Litecoin_Dur_SST_Sample3 <- VaRDurTest(0.10,Litecoin_Roll_SST_Sample3@forecast$VaR[,4],Litecoin_Roll_SST_Sample3@forecast$VaR[,3])
Litecoin_Dur_SST_Sample4 <- VaRDurTest(0.10,Litecoin_Roll_SST_Sample4@forecast$VaR[,4],Litecoin_Roll_SST_Sample4@forecast$VaR[,3])
Litecoin_Dur_SST_Sample5 <- VaRDurTest(0.10,Litecoin_Roll_SST_Sample5@forecast$VaR[,4],Litecoin_Roll_SST_Sample5@forecast$VaR[,3])

Bitcoin_Dur_SST_Sample1 <- VaRDurTest(0.10,Bitcoin_Roll_SST_Sample1@forecast$VaR[,4],Bitcoin_Roll_SST_Sample1@forecast$VaR[,3])
Bitcoin_Dur_SST_Sample2 <- VaRDurTest(0.10,Bitcoin_Roll_SST_Sample2@forecast$VaR[,4],Bitcoin_Roll_SST_Sample2@forecast$VaR[,3])
Bitcoin_Dur_SST_Sample3 <- VaRDurTest(0.10,Bitcoin_Roll_SST_Sample3@forecast$VaR[,4],Bitcoin_Roll_SST_Sample3@forecast$VaR[,3])
Bitcoin_Dur_SST_Sample4 <- VaRDurTest(0.10,Bitcoin_Roll_SST_Sample4@forecast$VaR[,4],Bitcoin_Roll_SST_Sample4@forecast$VaR[,3])
Bitcoin_Dur_SST_Sample5 <- VaRDurTest(0.10,Bitcoin_Roll_SST_Sample5@forecast$VaR[,4],Bitcoin_Roll_SST_Sample5@forecast$VaR[,3])

Litecoin_Dur_SGED_Sample1 <- VaRDurTest(0.10,Litecoin_Roll_SGED_Sample1@forecast$VaR[,4],Litecoin_Roll_SGED_Sample1@forecast$VaR[,3])
Litecoin_Dur_SGED_Sample2 <- VaRDurTest(0.10,Litecoin_Roll_SGED_Sample2@forecast$VaR[,4],Litecoin_Roll_SGED_Sample2@forecast$VaR[,3])
Litecoin_Dur_SGED_Sample3 <- VaRDurTest(0.10,Litecoin_Roll_SGED_Sample3@forecast$VaR[,4],Litecoin_Roll_SGED_Sample3@forecast$VaR[,3])
Litecoin_Dur_SGED_Sample4 <- VaRDurTest(0.10,Litecoin_Roll_SGED_Sample4@forecast$VaR[,4],Litecoin_Roll_SGED_Sample4@forecast$VaR[,3])
Litecoin_Dur_SGED_Sample5 <- VaRDurTest(0.10,Litecoin_Roll_SGED_Sample5@forecast$VaR[,4],Litecoin_Roll_SGED_Sample5@forecast$VaR[,3])

Bitcoin_Dur_SGED_Sample1 <- VaRDurTest(0.10,Bitcoin_Roll_SGED_Sample1@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample1@forecast$VaR[,3])
Bitcoin_Dur_SGED_Sample2 <- VaRDurTest(0.10,Bitcoin_Roll_SGED_Sample2@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample2@forecast$VaR[,3])
Bitcoin_Dur_SGED_Sample3 <- VaRDurTest(0.10,Bitcoin_Roll_SGED_Sample3@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample3@forecast$VaR[,3])
Bitcoin_Dur_SGED_Sample4 <- VaRDurTest(0.10,Bitcoin_Roll_SGED_Sample4@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample4@forecast$VaR[,3])
Bitcoin_Dur_SGED_Sample5 <- VaRDurTest(0.10,Bitcoin_Roll_SGED_Sample5@forecast$VaR[,4],Bitcoin_Roll_SGED_Sample5@forecast$VaR[,3])

Litecoin_Dur_JSU_Sample1 <- VaRDurTest(0.10,Litecoin_Roll_JSU_Sample1@forecast$VaR[,4],Litecoin_Roll_JSU_Sample1@forecast$VaR[,3])
Litecoin_Dur_JSU_Sample2 <- VaRDurTest(0.10,Litecoin_Roll_JSU_Sample2@forecast$VaR[,4],Litecoin_Roll_JSU_Sample2@forecast$VaR[,3])
Litecoin_Dur_JSU_Sample3 <- VaRDurTest(0.10,Litecoin_Roll_JSU_Sample3@forecast$VaR[,4],Litecoin_Roll_JSU_Sample3@forecast$VaR[,3])
Litecoin_Dur_JSU_Sample4 <- VaRDurTest(0.10,Litecoin_Roll_JSU_Sample4@forecast$VaR[,4],Litecoin_Roll_JSU_Sample4@forecast$VaR[,3])
Litecoin_Dur_JSU_Sample5 <- VaRDurTest(0.10,Litecoin_Roll_JSU_Sample5@forecast$VaR[,4],Litecoin_Roll_JSU_Sample5@forecast$VaR[,3])

Bitcoin_Dur_JSU_Sample1 <- VaRDurTest(0.10,Bitcoin_Roll_JSU_Sample1@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample1@forecast$VaR[,3])
Bitcoin_Dur_JSU_Sample2 <- VaRDurTest(0.10,Bitcoin_Roll_JSU_Sample2@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample2@forecast$VaR[,3])
Bitcoin_Dur_JSU_Sample3 <- VaRDurTest(0.10,Bitcoin_Roll_JSU_Sample3@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample3@forecast$VaR[,3])
Bitcoin_Dur_JSU_Sample4 <- VaRDurTest(0.10,Bitcoin_Roll_JSU_Sample4@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample4@forecast$VaR[,3])
Bitcoin_Dur_JSU_Sample5 <- VaRDurTest(0.10,Bitcoin_Roll_JSU_Sample5@forecast$VaR[,4],Bitcoin_Roll_JSU_Sample5@forecast$VaR[,3])

## Store the results

Litecoin_Dur_Subsample_stat <- matrix(c(Litecoin_Dur_Normal_Sample1$rLL,
                                   Litecoin_Dur_SST_Sample1$rLL,
                                   Litecoin_Dur_SGED_Sample1$rLL,
                                   Litecoin_Dur_JSU_Sample1$rLL,
                                   Litecoin_Dur_Normal_Sample2$rLL,
                                   Litecoin_Dur_SST_Sample2$rLL,
                                   Litecoin_Dur_SGED_Sample2$rLL,
                                   Litecoin_Dur_JSU_Sample2$rLL,
                                   Litecoin_Dur_Normal_Sample3$rLL,
                                   Litecoin_Dur_SST_Sample3$rLL,
                                   Litecoin_Dur_SGED_Sample3$rLL,
                                   Litecoin_Dur_JSU_Sample3$rLL,
                                   Litecoin_Dur_Normal_Sample4$rLL,
                                   Litecoin_Dur_SST_Sample4$rLL,
                                   Litecoin_Dur_SGED_Sample4$rLL,
                                   Litecoin_Dur_JSU_Sample4$rLL,
                                   Litecoin_Dur_Normal_Sample5$rLL,
                                   Litecoin_Dur_SST_Sample5$rLL,
                                   Litecoin_Dur_SGED_Sample5$rLL,
                                   Litecoin_Dur_JSU_Sample5$rLL),nrow=4)

Bitcoin_Dur_Subsample_stat <- matrix(c(Bitcoin_Dur_Normal_Sample1$rLL,
                                     Bitcoin_Dur_SST_Sample1$rLL,
                                     Bitcoin_Dur_SGED_Sample1$rLL,
                                     Bitcoin_Dur_JSU_Sample1$rLL,
                                     Bitcoin_Dur_Normal_Sample2$rLL,
                                     Bitcoin_Dur_SST_Sample2$rLL,
                                     Bitcoin_Dur_SGED_Sample2$rLL,
                                     Bitcoin_Dur_JSU_Sample2$rLL,
                                     Bitcoin_Dur_Normal_Sample3$rLL,
                                     Bitcoin_Dur_SST_Sample3$rLL,
                                     Bitcoin_Dur_SGED_Sample3$rLL,
                                     Bitcoin_Dur_JSU_Sample3$rLL,
                                     Bitcoin_Dur_Normal_Sample4$rLL,
                                     Bitcoin_Dur_SST_Sample4$rLL,
                                     Bitcoin_Dur_SGED_Sample4$rLL,
                                     Bitcoin_Dur_JSU_Sample4$rLL,
                                     Bitcoin_Dur_Normal_Sample5$rLL,
                                     Bitcoin_Dur_SST_Sample5$rLL,
                                     Bitcoin_Dur_SGED_Sample5$rLL,
                                     Bitcoin_Dur_JSU_Sample5$rLL),nrow=4)

Litecoin_Dur_Subsample_pvalues <- matrix(c(Litecoin_Dur_Normal_Sample1$LRp,
                                      Litecoin_Dur_SST_Sample1$LRp,
                                      Litecoin_Dur_SGED_Sample1$LRp,
                                      Litecoin_Dur_JSU_Sample1$LRp,
                                      Litecoin_Dur_Normal_Sample2$LRp,
                                      Litecoin_Dur_SST_Sample2$LRp,
                                      Litecoin_Dur_SGED_Sample2$LRp,
                                      Litecoin_Dur_JSU_Sample2$LRp,
                                      Litecoin_Dur_Normal_Sample3$LRp,
                                      Litecoin_Dur_SST_Sample3$LRp,
                                      Litecoin_Dur_SGED_Sample3$LRp,
                                      Litecoin_Dur_JSU_Sample3$LRp,
                                      Litecoin_Dur_Normal_Sample4$LRp,
                                      Litecoin_Dur_SST_Sample4$LRp,
                                      Litecoin_Dur_SGED_Sample4$LRp,
                                      Litecoin_Dur_JSU_Sample4$LRp,
                                      Litecoin_Dur_Normal_Sample5$LRp,
                                      Litecoin_Dur_SST_Sample5$LRp,
                                      Litecoin_Dur_SGED_Sample5$LRp,
                                      Litecoin_Dur_JSU_Sample5$LRp),nrow=4)

Bitcoin_Dur_Subsample_pvalues <- matrix(c(Bitcoin_Dur_Normal_Sample1$LRp,
                                        Bitcoin_Dur_SST_Sample1$LRp,
                                        Bitcoin_Dur_SGED_Sample1$LRp,
                                        Bitcoin_Dur_JSU_Sample1$LRp,
                                        Bitcoin_Dur_Normal_Sample2$LRp,
                                        Bitcoin_Dur_SST_Sample2$LRp,
                                        Bitcoin_Dur_SGED_Sample2$LRp,
                                        Bitcoin_Dur_JSU_Sample2$LRp,
                                        Bitcoin_Dur_Normal_Sample3$LRp,
                                        Bitcoin_Dur_SST_Sample3$LRp,
                                        Bitcoin_Dur_SGED_Sample3$LRp,
                                        Bitcoin_Dur_JSU_Sample3$LRp,
                                        Bitcoin_Dur_Normal_Sample4$LRp,
                                        Bitcoin_Dur_SST_Sample4$LRp,
                                        Bitcoin_Dur_SGED_Sample4$LRp,
                                        Bitcoin_Dur_JSU_Sample4$LRp,
                                        Bitcoin_Dur_Normal_Sample5$LRp,
                                        Bitcoin_Dur_SST_Sample5$LRp,
                                        Bitcoin_Dur_SGED_Sample5$LRp,
                                        Bitcoin_Dur_JSU_Sample5$LRp),nrow=4)


### Expected Shortfall Calculations ####

## Assuming Normal Distribution

Litecoin_index_normal1 <- which(Litecoin_Roll8@forecast$VaR[,1]>Litecoin_Roll8@forecast$VaR[,4])
Litecoin_VaR_values_normal1 <- Litecoin_Roll8@forecast$VaR[Litecoin_index_normal1,1]
Litecoin_return_values_normal1 <- Litecoin_Roll8@forecast$VaR[Litecoin_index_normal1,4]
Litecoin_ES_normal1 <- (Litecoin_return_values_normal1/Litecoin_VaR_values_normal1)
Litecoin_ES_Measure1_normal1 <- mean(Litecoin_ES_normal1)
Litecoin_ES_Measure2_normal1 <- max(Litecoin_ES_normal1)

Litecoin_index_normal5 <- which(Litecoin_Roll8@forecast$VaR[,2]>Litecoin_Roll8@forecast$VaR[,4])
Litecoin_VaR_values_normal5 <- Litecoin_Roll8@forecast$VaR[Litecoin_index_normal5,2]
Litecoin_return_values_normal5 <- Litecoin_Roll8@forecast$VaR[Litecoin_index_normal5,4]
Litecoin_ES_normal5 <- (Litecoin_return_values_normal5/Litecoin_VaR_values_normal5)
Litecoin_ES_Measure1_normal5 <- mean(Litecoin_ES_normal5)
Litecoin_ES_Measure2_normal5 <- max(Litecoin_ES_normal5)

Litecoin_index_normal10 <- which(Litecoin_Roll8@forecast$VaR[,3]>Litecoin_Roll8@forecast$VaR[,4])
Litecoin_VaR_values_normal10 <- Litecoin_Roll8@forecast$VaR[Litecoin_index_normal10,3]
Litecoin_return_values_normal10 <- Litecoin_Roll8@forecast$VaR[Litecoin_index_normal10,4]
Litecoin_ES_normal10 <- (Litecoin_return_values_normal10/Litecoin_VaR_values_normal10)
Litecoin_ES_Measure1_normal10 <- mean(Litecoin_ES_normal10)
Litecoin_ES_Measure2_normal10 <- max(Litecoin_ES_normal10)

Bitcoin_index_normal1 <- which(Bitcoin_Roll8@forecast$VaR[,1]>Bitcoin_Roll8@forecast$VaR[,4])
Bitcoin_VaR_values_normal1 <- Bitcoin_Roll8@forecast$VaR[Bitcoin_index_normal1,1]
Bitcoin_return_values_normal1 <- Bitcoin_Roll8@forecast$VaR[Bitcoin_index_normal1,4]
Bitcoin_ES_normal1 <- (Bitcoin_return_values_normal1/Bitcoin_VaR_values_normal1)
Bitcoin_ES_Measure1_normal1 <- mean(Bitcoin_ES_normal1)
Bitcoin_ES_Measure2_normal1 <- max(Bitcoin_ES_normal1)

Bitcoin_index_normal5 <- which(Bitcoin_Roll8@forecast$VaR[,2]>Bitcoin_Roll8@forecast$VaR[,4])
Bitcoin_VaR_values_normal5 <- Bitcoin_Roll8@forecast$VaR[Bitcoin_index_normal5,2]
Bitcoin_return_values_normal5 <- Bitcoin_Roll8@forecast$VaR[Bitcoin_index_normal5,4]
Bitcoin_ES_normal5 <- (Bitcoin_return_values_normal5/Bitcoin_VaR_values_normal5)
Bitcoin_ES_Measure1_normal5 <- mean(Bitcoin_ES_normal5)
Bitcoin_ES_Measure2_normal5 <- max(Bitcoin_ES_normal5)

Bitcoin_index_normal10 <- which(Bitcoin_Roll8@forecast$VaR[,3]>Bitcoin_Roll8@forecast$VaR[,4])
Bitcoin_VaR_values_normal10 <- Bitcoin_Roll8@forecast$VaR[Bitcoin_index_normal10,3]
Bitcoin_return_values_normal10 <- Bitcoin_Roll8@forecast$VaR[Bitcoin_index_normal10,4]
Bitcoin_ES_normal10 <- (Bitcoin_return_values_normal10/Bitcoin_VaR_values_normal10)
Bitcoin_ES_Measure1_normal10 <- mean(Bitcoin_ES_normal10)
Bitcoin_ES_Measure2_normal10 <- max(Bitcoin_ES_normal10)

## Assuming SST Distribution

Litecoin_index_SST1 <- which(Litecoin_Roll5@forecast$VaR[,1]>Litecoin_Roll5@forecast$VaR[,4])
Litecoin_VaR_values_SST1 <- Litecoin_Roll5@forecast$VaR[Litecoin_index_SST1,1]
Litecoin_return_values_SST1 <- Litecoin_Roll5@forecast$VaR[Litecoin_index_SST1,4]
Litecoin_ES_SST1 <- (Litecoin_return_values_SST1/Litecoin_VaR_values_SST1)
Litecoin_ES_Measure1_SST1 <- mean(Litecoin_ES_SST1)
Litecoin_ES_Measure2_SST1 <- max(Litecoin_ES_SST1)

Litecoin_index_SST5 <- which(Litecoin_Roll5@forecast$VaR[,2]>Litecoin_Roll5@forecast$VaR[,4])
Litecoin_VaR_values_SST5 <- Litecoin_Roll5@forecast$VaR[Litecoin_index_SST5,2]
Litecoin_return_values_SST5 <- Litecoin_Roll5@forecast$VaR[Litecoin_index_SST5,4]
Litecoin_ES_SST5 <- (Litecoin_return_values_SST5/Litecoin_VaR_values_SST5)
Litecoin_ES_Measure1_SST5 <- mean(Litecoin_ES_SST5)
Litecoin_ES_Measure2_SST5 <- max(Litecoin_ES_SST5)

Litecoin_index_SST10 <- which(Litecoin_Roll5@forecast$VaR[,3]>Litecoin_Roll5@forecast$VaR[,4])
Litecoin_VaR_values_SST10 <- Litecoin_Roll5@forecast$VaR[Litecoin_index_SST10,3]
Litecoin_return_values_SST10 <- Litecoin_Roll5@forecast$VaR[Litecoin_index_SST10,4]
Litecoin_ES_SST10 <- (Litecoin_return_values_SST10/Litecoin_VaR_values_SST10)
Litecoin_ES_Measure1_SST10 <- mean(Litecoin_ES_SST10)
Litecoin_ES_Measure2_SST10 <- max(Litecoin_ES_SST10)

Bitcoin_index_SST1 <- which(Bitcoin_Roll5@forecast$VaR[,1]>Bitcoin_Roll5@forecast$VaR[,4])
Bitcoin_VaR_values_SST1 <- Bitcoin_Roll5@forecast$VaR[Bitcoin_index_SST1,1]
Bitcoin_return_values_SST1 <- Bitcoin_Roll5@forecast$VaR[Bitcoin_index_SST1,4]
Bitcoin_ES_SST1 <- (Bitcoin_return_values_SST1/Bitcoin_VaR_values_SST1)
Bitcoin_ES_Measure1_SST1 <- mean(Bitcoin_ES_SST1)
Bitcoin_ES_Measure2_SST1 <- max(Bitcoin_ES_SST1)

Bitcoin_index_SST5 <- which(Bitcoin_Roll5@forecast$VaR[,2]>Bitcoin_Roll5@forecast$VaR[,4])
Bitcoin_VaR_values_SST5 <- Bitcoin_Roll5@forecast$VaR[Bitcoin_index_SST5,2]
Bitcoin_return_values_SST5 <- Bitcoin_Roll5@forecast$VaR[Bitcoin_index_SST5,4]
Bitcoin_ES_SST5 <- (Bitcoin_return_values_SST5/Bitcoin_VaR_values_SST5)
Bitcoin_ES_Measure1_SST5 <- mean(Bitcoin_ES_SST5)
Bitcoin_ES_Measure2_SST5 <- max(Bitcoin_ES_SST5)

Bitcoin_index_SST10 <- which(Bitcoin_Roll5@forecast$VaR[,3]>Bitcoin_Roll5@forecast$VaR[,4])
Bitcoin_VaR_values_SST10 <- Bitcoin_Roll5@forecast$VaR[Bitcoin_index_SST10,3]
Bitcoin_return_values_SST10 <- Bitcoin_Roll5@forecast$VaR[Bitcoin_index_SST10,4]
Bitcoin_ES_SST10 <- (Bitcoin_return_values_SST10/Bitcoin_VaR_values_SST10)
Bitcoin_ES_Measure1_SST10 <- mean(Bitcoin_ES_SST10)
Bitcoin_ES_Measure2_SST10 <- max(Bitcoin_ES_SST10)


## Assuming SGED Distribution

Litecoin_index_SGED1 <- which(Litecoin_Roll6@forecast$VaR[,1]>Litecoin_Roll6@forecast$VaR[,4])
Litecoin_VaR_values_SGED1 <- Litecoin_Roll6@forecast$VaR[Litecoin_index_SGED1,1]
Litecoin_return_values_SGED1 <- Litecoin_Roll6@forecast$VaR[Litecoin_index_SGED1,4]
Litecoin_ES_SGED1 <- (Litecoin_return_values_SGED1/Litecoin_VaR_values_SGED1)
Litecoin_ES_Measure1_SGED1 <- mean(Litecoin_ES_SGED1)
Litecoin_ES_Measure2_SGED1 <- max(Litecoin_ES_SGED1)

Litecoin_index_SGED5 <- which(Litecoin_Roll6@forecast$VaR[,2]>Litecoin_Roll6@forecast$VaR[,4])
Litecoin_VaR_values_SGED5 <- Litecoin_Roll6@forecast$VaR[Litecoin_index_SGED5,2]
Litecoin_return_values_SGED5 <- Litecoin_Roll6@forecast$VaR[Litecoin_index_SGED5,4]
Litecoin_ES_SGED5 <- (Litecoin_return_values_SGED5/Litecoin_VaR_values_SGED5)
Litecoin_ES_Measure1_SGED5 <- mean(Litecoin_ES_SGED5)
Litecoin_ES_Measure2_SGED5 <- max(Litecoin_ES_SGED5)

Litecoin_index_SGED10 <- which(Litecoin_Roll6@forecast$VaR[,3]>Litecoin_Roll6@forecast$VaR[,4])
Litecoin_VaR_values_SGED10 <- Litecoin_Roll6@forecast$VaR[Litecoin_index_SGED10,3]
Litecoin_return_values_SGED10 <- Litecoin_Roll6@forecast$VaR[Litecoin_index_SGED10,4]
Litecoin_ES_SGED10 <- (Litecoin_return_values_SGED10/Litecoin_VaR_values_SGED10)
Litecoin_ES_Measure1_SGED10 <- mean(Litecoin_ES_SGED10)
Litecoin_ES_Measure2_SGED10 <- max(Litecoin_ES_SGED10)

Bitcoin_index_SGED1 <- which(Bitcoin_Roll6@forecast$VaR[,1]>Bitcoin_Roll6@forecast$VaR[,4])
Bitcoin_VaR_values_SGED1 <- Bitcoin_Roll6@forecast$VaR[Bitcoin_index_SGED1,1]
Bitcoin_return_values_SGED1 <- Bitcoin_Roll6@forecast$VaR[Bitcoin_index_SGED1,4]
Bitcoin_ES_SGED1 <- (Bitcoin_return_values_SGED1/Bitcoin_VaR_values_SGED1)
Bitcoin_ES_Measure1_SGED1 <- mean(Bitcoin_ES_SGED1)
Bitcoin_ES_Measure2_SGED1 <- max(Bitcoin_ES_SGED1)

Bitcoin_index_SGED5 <- which(Bitcoin_Roll6@forecast$VaR[,2]>Bitcoin_Roll6@forecast$VaR[,4])
Bitcoin_VaR_values_SGED5 <- Bitcoin_Roll6@forecast$VaR[Bitcoin_index_SGED5,2]
Bitcoin_return_values_SGED5 <- Bitcoin_Roll6@forecast$VaR[Bitcoin_index_SGED5,4]
Bitcoin_ES_SGED5 <- (Bitcoin_return_values_SGED5/Bitcoin_VaR_values_SGED5)
Bitcoin_ES_Measure1_SGED5 <- mean(Bitcoin_ES_SGED5)
Bitcoin_ES_Measure2_SGED5 <- max(Bitcoin_ES_SGED5)

Bitcoin_index_SGED10 <- which(Bitcoin_Roll6@forecast$VaR[,3]>Bitcoin_Roll6@forecast$VaR[,4])
Bitcoin_VaR_values_SGED10 <- Bitcoin_Roll6@forecast$VaR[Bitcoin_index_SGED10,3]
Bitcoin_return_values_SGED10 <- Bitcoin_Roll6@forecast$VaR[Bitcoin_index_SGED10,4]
Bitcoin_ES_SGED10 <- (Bitcoin_return_values_SGED10/Bitcoin_VaR_values_SGED10)
Bitcoin_ES_Measure1_SGED10 <- mean(Bitcoin_ES_SGED10)
Bitcoin_ES_Measure2_SGED10 <- max(Bitcoin_ES_SGED10)


## Assuming JSU Distribution

Litecoin_index_JSU1 <- which(Litecoin_Roll7@forecast$VaR[,1]>Litecoin_Roll7@forecast$VaR[,4])
Litecoin_VaR_values_JSU1 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_JSU1,1]
Litecoin_return_values_JSU1 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_JSU1,4]
Litecoin_ES_JSU1 <- (Litecoin_return_values_JSU1/Litecoin_VaR_values_JSU1)
Litecoin_ES_Measure1_JSU1 <- mean(Litecoin_ES_JSU1)
Litecoin_ES_Measure2_JSU1 <- max(Litecoin_ES_JSU1)

Litecoin_index_JSU5 <- which(Litecoin_Roll7@forecast$VaR[,2]>Litecoin_Roll7@forecast$VaR[,4])
Litecoin_VaR_values_JSU5 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_JSU5,2]
Litecoin_return_values_JSU5 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_JSU5,4]
Litecoin_ES_JSU5 <- (Litecoin_return_values_JSU5/Litecoin_VaR_values_JSU5)
Litecoin_ES_Measure1_JSU5 <- mean(Litecoin_ES_JSU5)
Litecoin_ES_Measure2_JSU5 <- max(Litecoin_ES_JSU5)

Litecoin_index_JSU10 <- which(Litecoin_Roll7@forecast$VaR[,3]>Litecoin_Roll7@forecast$VaR[,4])
Litecoin_VaR_values_JSU10 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_JSU10,3]
Litecoin_return_values_JSU10 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_JSU10,4]
Litecoin_ES_JSU10 <- (Litecoin_return_values_JSU10/Litecoin_VaR_values_JSU10)
Litecoin_ES_Measure1_JSU10 <- mean(Litecoin_ES_JSU10)
Litecoin_ES_Measure2_JSU10 <- max(Litecoin_ES_JSU10)

Bitcoin_index_JSU1 <- which(Bitcoin_Roll7@forecast$VaR[,1]>Bitcoin_Roll7@forecast$VaR[,4])
Bitcoin_VaR_values_JSU1 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_JSU1,1]
Bitcoin_return_values_JSU1 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_JSU1,4]
Bitcoin_ES_JSU1 <- (Bitcoin_return_values_JSU1/Bitcoin_VaR_values_JSU1)
Bitcoin_ES_Measure1_JSU1 <- mean(Bitcoin_ES_JSU1)
Bitcoin_ES_Measure2_JSU1 <- max(Bitcoin_ES_JSU1)

Bitcoin_index_JSU5 <- which(Bitcoin_Roll7@forecast$VaR[,2]>Bitcoin_Roll7@forecast$VaR[,4])
Bitcoin_VaR_values_JSU5 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_JSU5,2]
Bitcoin_return_values_JSU5 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_JSU5,4]
Bitcoin_ES_JSU5 <- (Bitcoin_return_values_JSU5/Bitcoin_VaR_values_JSU5)
Bitcoin_ES_Measure1_JSU5 <- mean(Bitcoin_ES_JSU5)
Bitcoin_ES_Measure2_JSU5 <- max(Bitcoin_ES_JSU5)

Bitcoin_index_JSU10 <- which(Bitcoin_Roll7@forecast$VaR[,3]>Bitcoin_Roll7@forecast$VaR[,4])
Bitcoin_VaR_values_JSU10 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_JSU10,3]
Bitcoin_return_values_JSU10 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_JSU10,4]
Bitcoin_ES_JSU10 <- (Bitcoin_return_values_JSU10/Bitcoin_VaR_values_JSU10)
Bitcoin_ES_Measure1_JSU10 <- mean(Bitcoin_ES_JSU10)
Bitcoin_ES_Measure2_JSU10 <- max(Bitcoin_ES_JSU10)

### For Pearson

Litecoin_index_Pearson1 <- which(Litecoin_Holdout_VaR[-1]>Litecoin_Roll7@forecast$VaR[1:999,4])
Litecoin_VaR_values_Pearson1 <- Litecoin_Holdout_VaR[Litecoin_index_Pearson1]
Litecoin_return_values_Pearson1 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_Pearson1,4]
Litecoin_ES_Pearson1 <- (Litecoin_return_values_Pearson1/Litecoin_VaR_values_Pearson1)
Litecoin_ES_Measure1_Pearson1 <- mean(Litecoin_ES_Pearson1)
Litecoin_ES_Measure2_Pearson1 <- max(Litecoin_ES_Pearson1)

Litecoin_index_Pearson5 <- which(Litecoin_Holdout_VaR[-1]>Litecoin_Roll7@forecast$VaR[1:999,4])
Litecoin_VaR_values_Pearson5 <- Litecoin_Holdout_VaR[Litecoin_index_Pearson5]
Litecoin_return_values_Pearson5 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_Pearson5,4]
Litecoin_ES_Pearson5 <- (Litecoin_return_values_Pearson5/Litecoin_VaR_values_Pearson5)
Litecoin_ES_Measure1_Pearson5 <- mean(Litecoin_ES_Pearson5)
Litecoin_ES_Measure2_Pearson5 <- max(Litecoin_ES_Pearson5)

Litecoin_index_Pearson10 <- which(Litecoin_Holdout_VaR[-1]>Litecoin_Roll7@forecast$VaR[1:999,4])
Litecoin_VaR_values_Pearson10 <- Litecoin_Holdout_VaR[Litecoin_index_Pearson10]
Litecoin_return_values_Pearson10 <- Litecoin_Roll7@forecast$VaR[Litecoin_index_Pearson10,4]
Litecoin_ES_Pearson10 <- (Litecoin_return_values_Pearson10/Litecoin_VaR_values_Pearson10)
Litecoin_ES_Measure1_Pearson10 <- mean(Litecoin_ES_Pearson10)
Litecoin_ES_Measure2_Pearson10 <- max(Litecoin_ES_Pearson10)

Bitcoin_index_Pearson1 <- which(Bitcoin_Holdout_VaR[-1]>Bitcoin_Roll7@forecast$VaR[1:999,4])
Bitcoin_VaR_values_Pearson1 <- Bitcoin_Holdout_VaR[Bitcoin_index_Pearson1]
Bitcoin_return_values_Pearson1 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_Pearson1,4]
Bitcoin_ES_Pearson1 <- (Bitcoin_return_values_Pearson1/Bitcoin_VaR_values_Pearson1)
Bitcoin_ES_Measure1_Pearson1 <- mean(Bitcoin_ES_Pearson1)
Bitcoin_ES_Measure2_Pearson1 <- max(Bitcoin_ES_Pearson1)

Bitcoin_index_Pearson5 <- which(Bitcoin_Holdout_VaR[-1]>Bitcoin_Roll7@forecast$VaR[1:999,4])
Bitcoin_VaR_values_Pearson5 <- Bitcoin_Holdout_VaR[Bitcoin_index_Pearson5]
Bitcoin_return_values_Pearson5 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_Pearson5,4]
Bitcoin_ES_Pearson5 <- (Bitcoin_return_values_Pearson5/Bitcoin_VaR_values_Pearson5)
Bitcoin_ES_Measure1_Pearson5 <- mean(Bitcoin_ES_Pearson5)
Bitcoin_ES_Measure2_Pearson5 <- max(Bitcoin_ES_Pearson5)

Bitcoin_index_Pearson10 <- which(Bitcoin_Holdout_VaR[-1]>Bitcoin_Roll7@forecast$VaR[1:999,4])
Bitcoin_VaR_values_Pearson10 <- Bitcoin_Holdout_VaR[Bitcoin_index_Pearson10]
Bitcoin_return_values_Pearson10 <- Bitcoin_Roll7@forecast$VaR[Bitcoin_index_Pearson10,4]
Bitcoin_ES_Pearson10 <- (Bitcoin_return_values_Pearson10/Bitcoin_VaR_values_Pearson10)
Bitcoin_ES_Measure1_Pearson10 <- mean(Bitcoin_ES_Pearson10)
Bitcoin_ES_Measure2_Pearson10 <- max(Bitcoin_ES_Pearson10)

### Store for Pearson

Pearson_Full_ES_Measure1 <- matrix(c(Litecoin_ES_Measure1_Pearson1,Litecoin_ES_Measure1_Pearson5,Litecoin_ES_Measure1_Pearson10,
                                     Bitcoin_ES_Measure1_Pearson1,Bitcoin_ES_Measure1_Pearson5,Bitcoin_ES_Measure1_Pearson10),nrow=3)

Pearson_Full_ES_Measure2 <- matrix(c(Litecoin_ES_Measure2_Pearson1,Litecoin_ES_Measure2_Pearson5,Litecoin_ES_Measure2_Pearson10,
                                     Bitcoin_ES_Measure2_Pearson1,Bitcoin_ES_Measure2_Pearson5,Bitcoin_ES_Measure2_Pearson10),nrow=3)

rownames(Pearson_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(Pearson_Full_ES_Measure2) <- c("1%","5%","10%")
colnames(Pearson_Full_ES_Measure1) <- c("Litecoin","Bitcoin")
colnames(Pearson_Full_ES_Measure2) <- c("Litecoin","Bitcoin")

#### Store the values

Litecoin_Full_ES_Measure1 <- matrix(c(Litecoin_ES_Measure1_normal1,Litecoin_ES_Measure1_SST1,Litecoin_ES_Measure1_SGED1,Litecoin_ES_Measure1_JSU1,
                                 Litecoin_ES_Measure1_normal5,Litecoin_ES_Measure1_SST5,Litecoin_ES_Measure1_SGED5,Litecoin_ES_Measure1_JSU5,
                                 Litecoin_ES_Measure1_normal10,Litecoin_ES_Measure1_SST10,Litecoin_ES_Measure1_SGED10,Litecoin_ES_Measure1_JSU10),nrow=3,byrow=TRUE)


Bitcoin_Full_ES_Measure1 <- matrix(c(Bitcoin_ES_Measure1_normal1,Bitcoin_ES_Measure1_SST1,Bitcoin_ES_Measure1_SGED1,Bitcoin_ES_Measure1_JSU1,
                                   Bitcoin_ES_Measure1_normal5,Bitcoin_ES_Measure1_SST5,Bitcoin_ES_Measure1_SGED5,Bitcoin_ES_Measure1_JSU5,
                                   Bitcoin_ES_Measure1_normal10,Bitcoin_ES_Measure1_SST10,Bitcoin_ES_Measure1_SGED10,Bitcoin_ES_Measure1_JSU10),nrow=3,byrow=TRUE)

Litecoin_Full_ES_Measure2 <- matrix(c(Litecoin_ES_Measure2_normal1,Litecoin_ES_Measure2_SST1,Litecoin_ES_Measure2_SGED1,Litecoin_ES_Measure2_JSU1,
                                 Litecoin_ES_Measure2_normal5,Litecoin_ES_Measure2_SST5,Litecoin_ES_Measure2_SGED5,Litecoin_ES_Measure2_JSU5,
                                 Litecoin_ES_Measure2_normal10,Litecoin_ES_Measure2_SST10,Litecoin_ES_Measure2_SGED10,Litecoin_ES_Measure2_JSU10),nrow=3,byrow=TRUE)

Bitcoin_Full_ES_Measure2 <- matrix(c(Bitcoin_ES_Measure2_normal1,Bitcoin_ES_Measure2_SST1,Bitcoin_ES_Measure2_SGED1,Bitcoin_ES_Measure2_JSU1,
                                   Bitcoin_ES_Measure2_normal5,Bitcoin_ES_Measure2_SST5,Bitcoin_ES_Measure2_SGED5,Bitcoin_ES_Measure2_JSU5,
                                   Bitcoin_ES_Measure2_normal10,Bitcoin_ES_Measure2_SST10,Bitcoin_ES_Measure2_SGED10,Bitcoin_ES_Measure2_JSU10),nrow=3,byrow=TRUE)

rownames(Litecoin_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(Litecoin_Full_ES_Measure2) <- c("1%","5%","10%")
rownames(Bitcoin_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(Bitcoin_Full_ES_Measure2) <- c("1%","5%","10%")

colnames(Litecoin_Full_ES_Measure1) <- c("Normal","SST","SGED","JSU")
colnames(Litecoin_Full_ES_Measure2) <- c("Normal","SST","SGED","JSU")
colnames(Bitcoin_Full_ES_Measure1) <- c("Normal","SST","SGED","JSU")
colnames(Bitcoin_Full_ES_Measure2) <- c("Normal","SST","SGED","JSU")

### MAPE

# Litecoin
check1 <- Holdout_Returns_Litecoin[-c(which(Holdout_Returns_Litecoin==0))]
check2 <- as.data.frame(x=Litecoin_Roll5)[-c(which(Holdout_Returns_Litecoin==0)),1]
check3 <- as.data.frame(x=Litecoin_Roll6)[-c(which(Holdout_Returns_Litecoin==0)),1]
check4 <- as.data.frame(x=Litecoin_Roll7)[-c(which(Holdout_Returns_Litecoin==0)),1]
check5 <- as.data.frame(x=Litecoin_Roll8)[-c(which(Holdout_Returns_Litecoin==0)),1]

MAPE_Litecoin <- c(mape(check1,check2),
              mape(check1,check3),
              mape(check1,check4),
              mape(check1,check5))

## Bitcoin
check1 <- Holdout_Returns_Bitcoin[-c(which(Holdout_Returns_Bitcoin==0))]
check2 <- as.data.frame(x=Bitcoin_Roll5)[-c(which(Holdout_Returns_Bitcoin==0)),1]
check3 <- as.data.frame(x=Bitcoin_Roll6)[-c(which(Holdout_Returns_Bitcoin==0)),1]
check4 <- as.data.frame(x=Bitcoin_Roll7)[-c(which(Holdout_Returns_Bitcoin==0)),1]
check5 <- as.data.frame(x=Bitcoin_Roll8)[-c(which(Holdout_Returns_Bitcoin==0)),1]

MAPE_Bitcoin <- c(mape(check1,check2),
                mape(check1,check3),
                mape(check1,check4),
                mape(check1,check5))



#### End of all the code/lines ######  
