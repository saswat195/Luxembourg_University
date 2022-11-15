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

### We will use both Gas and Heat Oil data series ###

Gas <- read.csv("~//Desktop//Research//Energy_Related//Natural_Gas.csv")
Heat <- read.csv("~//Desktop//Research//Energy_Related//Heating_Oil.csv")  
Heat_Date <- as.Date.character(Heat$Date,format = "%m/%d/%Y")
Gas_Date   <- as.Date.character(Gas$Date,format = "%m/%d/%Y")
Gas   <- data.frame(Gas_Date,Gas$Price)     
Heat <- data.frame(Heat_Date,Heat$Price)
colnames(Gas) <- c("Date","Price")
colnames(Heat) <- c("Date","Price")

#### Structural Breaks for Heat and Gas #####

Heat_breaks <- breakpoints(Heat$Price~1)
Gas_breaks   <- breakpoints(Gas$Price~1)
Heat_samples <- Heat_breaks$breakpoints
Gas_samples <- Gas_breaks$breakpoints

Returns_Gas   <- log(Gas$Price[-1]/Gas$Price[1:(nrow(Gas)-1)])
Returns_Heat <- log(Heat$Price[-1]/Heat$Price[1:(nrow(Heat)-1)])
Des_Returns_Gas <- basicStats(Returns_Gas)
Des_Returns_Heat <- basicStats(Returns_Heat)

Holdout_Returns_Gas   <- Returns_Gas[(length(Returns_Gas)-999):length(Returns_Gas)] 
Holdout_Returns_Heat <- Returns_Heat[(length(Returns_Heat)-999):length(Returns_Heat)]

#### Convert into xts objects

Heat_xts <- xts(Returns_Heat,order.by = Heat_Date[-1])
Gas_xts   <- xts(Returns_Gas,order.by = Gas_Date[-1])
names(Heat_xts)    <- "Returns"
names(Gas_xts)   <- "Returns"

### Subsamples XTS  ###

Heat_sub1    <- Heat_xts[1:Heat_samples[1]]
Heat_sub2    <- Heat_xts[Heat_samples[1]:Heat_samples[2]]
Heat_sub3    <- Heat_xts[Heat_samples[2]:Heat_samples[3]]
Heat_sub4    <- Heat_xts[Heat_samples[3]:Heat_samples[4]]
Heat_sub5    <- Heat_xts[Heat_samples[4]:nrow(Heat_xts)]

Gas_sub1    <- Gas_xts[1:Gas_samples[1]]
Gas_sub2    <- Gas_xts[Gas_samples[1]:Gas_samples[2]]
Gas_sub3    <- Gas_xts[Gas_samples[2]:Gas_samples[3]]
Gas_sub4    <- Gas_xts[Gas_samples[3]:Gas_samples[4]]
Gas_sub5    <- Gas_xts[Gas_samples[4]:nrow(Gas_xts)]

###Dummy variables for structural breaks ####

Heat_Dummy1 <- array(data=0,dim=length(Heat_xts))
Heat_Dummy2 <- array(data=0,dim=length(Heat_xts))
Heat_Dummy3 <- array(data=0,dim=length(Heat_xts))
Heat_Dummy4 <- array(data=0,dim=length(Heat_xts))
Heat_Dummy5 <- array(data=0,dim=length(Heat_xts))

Heat_Dummy1[1:Heat_samples[1]]  <- 1
Heat_Dummy2[(Heat_samples[1]+1):Heat_samples[2]] <- 1
Heat_Dummy3[(Heat_samples[2]+1):Heat_samples[3]] <- 1
Heat_Dummy4[(Heat_samples[3]+1):Heat_samples[4]] <- 1
Heat_Dummy5[(Heat_samples[4]+1):length(Heat_xts)] <- 1

Gas_Dummy1 <- array(data=0,dim=length(Gas_xts))
Gas_Dummy2 <- array(data=0,dim=length(Gas_xts))
Gas_Dummy3 <- array(data=0,dim=length(Gas_xts))
Gas_Dummy4 <- array(data=0,dim=length(Gas_xts))
Gas_Dummy5 <- array(data=0,dim=length(Gas_xts))

Gas_Dummy1[1:Gas_samples[1]]  <- 1
Gas_Dummy2[(Gas_samples[1]+1):Gas_samples[2]] <- 1
Gas_Dummy3[(Gas_samples[2]+1):Gas_samples[3]] <- 1
Gas_Dummy4[(Gas_samples[3]+1):Gas_samples[4]] <- 1
Gas_Dummy5[(Gas_samples[4]+1):length(Gas_xts)] <- 1

Heat_Dummy <- cbind(Heat_Dummy1,Heat_Dummy2,Heat_Dummy3,Heat_Dummy4)
Gas_Dummy <- cbind(Gas_Dummy1,Gas_Dummy2,Gas_Dummy3,Gas_Dummy4)

### Descriptive Statistics and Basic Tests ###

Descr_Heat <-  cbind(basicStats(Heat_sub1),basicStats(Heat_sub2),basicStats(Heat_sub3),basicStats(Heat_sub4),basicStats(Heat_sub5))
Descr_Gas   <- cbind(basicStats(Gas_sub1),basicStats(Gas_sub2),basicStats(Gas_sub3),basicStats(Gas_sub4),basicStats(Gas_sub5))
ADF_Heat <- cbind(adf.test(Heat_sub1),adf.test(Heat_sub2),adf.test(Heat_sub3),adf.test(Heat_sub4),adf.test(Heat_sub5))
ADF_Gas   <- cbind(adf.test(Gas_sub1),adf.test(Gas_sub2),adf.test(Gas_sub3),adf.test(Gas_sub4),adf.test(Gas_sub5))
pp_Heat <- cbind(pp.test(Heat_sub1),pp.test(Heat_sub2),pp.test(Heat_sub3),pp.test(Heat_sub4),pp.test(Heat_sub5))
pp_Gas   <- cbind(pp.test(Gas_sub1),pp.test(Gas_sub2),pp.test(Gas_sub3),pp.test(Gas_sub4),pp.test(Gas_sub5))
jb_Heat <- cbind(jarque.bera.test(Heat_sub1),jarque.bera.test(Heat_sub2),jarque.bera.test(Heat_sub3),jarque.bera.test(Heat_sub4),jarque.bera.test(Heat_sub5))
jb_Gas   <- cbind(jarque.bera.test(Gas_sub1),jarque.bera.test(Gas_sub2),jarque.bera.test(Gas_sub3),jarque.bera.test(Gas_sub4),jarque.bera.test(Gas_sub5))

LQ_Heat5 <- cbind(Box.test(Heat_sub1^2,lag=5,type="Lj"),Box.test(Heat_sub2^2,lag=5,type="Lj"),Box.test(Heat_sub3^2,lag=5,type="Lj"),Box.test(Heat_sub4^2,lag=5,type="Lj"),Box.test(Heat_sub5^2,lag=5,type="Lj"))
LQ_Gas5   <- cbind(Box.test(Gas_sub1^2,lag=5,type="Lj"),Box.test(Gas_sub2^2,lag=5,type="Lj"),Box.test(Gas_sub3^2,lag=5,type="Lj"),Box.test(Gas_sub4^2,lag=5,type="Lj"),Box.test(Gas_sub5^2,lag=5,type="Lj"))

LQ_Heat10 <- cbind(Box.test(Heat_sub1^2,lag=10,type="Lj"),Box.test(Heat_sub2^2,lag=10,type="Lj"),Box.test(Heat_sub3^2,lag=10,type="Lj"),Box.test(Heat_sub4^2,lag=10,type="Lj"),Box.test(Heat_sub5^2,lag=10,type="Lj"))
LQ_Gas10   <- cbind(Box.test(Gas_sub1^2,lag=10,type="Lj"),Box.test(Gas_sub2^2,lag=10,type="Lj"),Box.test(Gas_sub3^2,lag=10,type="Lj"),Box.test(Gas_sub4^2,lag=10,type="Lj"),Box.test(Gas_sub5^2,lag=10,type="Lj"))

#### Plot the Returns and Oil Prices with Breakpoints indicated

par(mfrow=c(1,2))
plot(Heat,type="l",col="blue",main="Structural Breaks for Heat Oil",ylab="Price/Barrel (USD)",xaxt="n",xlab="")
abline(v=Heat$Date[Heat_samples[1]],col="brown")
abline(v=Heat$Date[Heat_samples[2]],col="brown")
abline(v=Heat$Date[Heat_samples[3]],col="brown")
abline(v=Heat$Date[Heat_samples[4]],col="brown")
Date_formatted <- format(Heat$Date[c(1,Heat_samples,nrow(Heat))],"%b-%y")
axis(1,at=c(min(Heat$Date),Heat$Date[Heat_samples],max(Heat$Date)),labels=Date_formatted,las = 2)

plot(Gas,type="l",col="blue",main="Structural Breaks for Gas Oil",ylab="Price/Barrel (USD)",xlab="",xaxt="n")
abline(v=Gas$Date[Gas_samples[1]],col="brown")
abline(v=Gas$Date[Gas_samples[2]],col="brown")
abline(v=Gas$Date[Gas_samples[3]],col="brown")
abline(v=Gas$Date[Gas_samples[4]],col="brown")
Date_formatted <- format(Gas$Date[c(1,Gas_samples,nrow(Gas))],"%b-%y")
axis(1,at=c(min(Gas$Date),Gas$Date[Gas_samples],max(Gas$Date)),labels=Date_formatted,las = 2)

### Plotting the Density and Histogram for Returns

par(mfrow=c(1,2))
densityPlot(as.timeSeries(Heat_xts),hist = TRUE,title = FALSE,grid=FALSE)
title(main="Density Plot for Heating Oil",xlab="Returns",ylab="Density/Frequency",col.main="blue",cex.main=0.77)
legend(x=-0.4,y=25,legend = c("Normal","Empirical"),bty="n",cex=0.7,text.col = c("gray","brown"))

densityPlot(as.timeSeries(Heat_xts),hist = TRUE,title = FALSE,grid=FALSE)
title(main="Density Plot for Natural Gas ",xlab="Returns",ylab="Density/Frequency",col.main="blue",cex.main=0.77)
legend(x=-0.4,y=25,legend = c("Normal","Empirical"),bty="n",cex=0.7,text.col = c("gray","brown"),col=c("gray","brown"))

### Garch Models

Heat_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Heat_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Heat_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Heat_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Heat_Model1 <- ugarchfit(Heat_Spec1,Heat_xts,out.sample = 1000)
Heat_Model2 <- ugarchfit(Heat_Spec2,Heat_xts,out.sample = 1000)
Heat_Model3 <- ugarchfit(Heat_Spec3,Heat_xts,out.sample = 1000)
Heat_Model4 <- ugarchfit(Heat_Spec4,Heat_xts,out.sample = 1000)

Gas_Spec1 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Gas_Spec2 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Gas_Spec3 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Gas_Spec4 <- ugarchspec(variance.model=list(model="sGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Gas_Model1 <- ugarchfit(Gas_Spec1,Gas_xts,out.sample = 1000)
Gas_Model2 <- ugarchfit(Gas_Spec2,Gas_xts,out.sample = 1000)
Gas_Model3 <- ugarchfit(Gas_Spec3,Gas_xts,out.sample = 1000)
Gas_Model4 <- ugarchfit(Gas_Spec4,Gas_xts,out.sample = 1000)

#### egarch models

Heat_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Heat_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Heat_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Heat_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Heat_Model5 <- ugarchfit(Heat_Spec5,Heat_xts,out.sample = 1000)
Heat_Model6 <- ugarchfit(Heat_Spec6,Heat_xts,out.sample = 1000)
Heat_Model7 <- ugarchfit(Heat_Spec7,Heat_xts,out.sample = 1000)
Heat_Model8 <- ugarchfit(Heat_Spec8,Heat_xts,out.sample = 1000)

Gas_Spec5 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Gas_Spec6 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Gas_Spec7 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Gas_Spec8 <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Gas_Model5 <- ugarchfit(Gas_Spec5,Gas_xts,out.sample = 1000)
Gas_Model6 <- ugarchfit(Gas_Spec6,Gas_xts,out.sample = 1000)
Gas_Model7 <- ugarchfit(Gas_Spec7,Gas_xts,out.sample = 1000)
Gas_Model8 <- ugarchfit(Gas_Spec8,Gas_xts,out.sample = 1000)

### gjrgarch models

Heat_Spec9   <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Heat_Spec10  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Heat_Spec11  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Heat_Spec12  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Heat_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Heat_Model9  <- ugarchfit(Heat_Spec9,Heat_xts,out.sample = 1000)
Heat_Model10 <- ugarchfit(Heat_Spec10,Heat_xts,out.sample = 1000)
Heat_Model11 <- ugarchfit(Heat_Spec11,Heat_xts,out.sample = 1000)
Heat_Model12 <- ugarchfit(Heat_Spec12,Heat_xts,out.sample = 1000)

Gas_Spec9  <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
Gas_Spec10 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")
Gas_Spec11 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")
Gas_Spec12 <- ugarchspec(variance.model=list(model="gjrGARCH",garchOrder=c(1,1),external.regressors=Gas_Dummy),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Gas_Model9  <- ugarchfit(Gas_Spec9,Gas_xts,out.sample = 1000)
Gas_Model10 <- ugarchfit(Gas_Spec10,Gas_xts,out.sample = 1000)
Gas_Model11 <- ugarchfit(Gas_Spec11,Gas_xts,out.sample = 1000)
Gas_Model12 <- ugarchfit(Gas_Spec12,Gas_xts,out.sample = 1000)

### Store the loglikelihood values, AIC and BIC ####

Llhood  <- matrix(data=0,nrow=12,ncol=2)
colnames(Llhood) <- c("Gas","Heat")
rownames(Llhood) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")
  
Llhood[1,1]  <-  Gas_Model1@fit$LLH
Llhood[2,1]  <-  Gas_Model2@fit$LLH
Llhood[3,1]  <-  Gas_Model3@fit$LLH
Llhood[4,1]  <-  Gas_Model4@fit$LLH
Llhood[5,1]  <-  Gas_Model5@fit$LLH
Llhood[6,1]  <-  Gas_Model6@fit$LLH
Llhood[7,1]  <-  Gas_Model7@fit$LLH
Llhood[8,1]  <-  Gas_Model8@fit$LLH
Llhood[9,1]  <-  Gas_Model9@fit$LLH
Llhood[10,1] <- Gas_Model10@fit$LLH
Llhood[11,1] <- Gas_Model11@fit$LLH
Llhood[12,1] <- Gas_Model12@fit$LLH

Llhood[1,2]  <-  Heat_Model1@fit$LLH
Llhood[2,2]  <-  Heat_Model2@fit$LLH
Llhood[3,2]  <-  Heat_Model3@fit$LLH
Llhood[4,2]  <-  Heat_Model4@fit$LLH
Llhood[5,2]  <-  Heat_Model5@fit$LLH
Llhood[6,2]  <-  Heat_Model6@fit$LLH
Llhood[7,2]  <-  Heat_Model7@fit$LLH
Llhood[8,2]  <-  Heat_Model8@fit$LLH
Llhood[9,2]  <-  Heat_Model9@fit$LLH
Llhood[10,2] <- Heat_Model10@fit$LLH
Llhood[11,2] <- Heat_Model11@fit$LLH
Llhood[12,2] <- Heat_Model12@fit$LLH

BIC  <- matrix(data=0,nrow=12,ncol=2)
colnames(BIC) <- c("Gas","Heat")
rownames(BIC) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")

BIC[1,1]  <-  infocriteria(Gas_Model1)[2]
BIC[2,1]  <-  infocriteria(Gas_Model2)[2]
BIC[3,1]  <-  infocriteria(Gas_Model3)[2]
BIC[4,1]  <-  infocriteria(Gas_Model4)[2]
BIC[5,1]  <-  infocriteria(Gas_Model5)[2]
BIC[6,1]  <-  infocriteria(Gas_Model6)[2]
BIC[7,1]  <-  infocriteria(Gas_Model7)[2]
BIC[8,1]  <-  infocriteria(Gas_Model8)[2]
BIC[9,1]  <-  infocriteria(Gas_Model9)[2]
BIC[10,1] <- infocriteria(Gas_Model10)[2]
BIC[11,1] <- infocriteria(Gas_Model11)[2]
BIC[12,1] <- infocriteria(Gas_Model12)[2]

BIC[1,2]  <-  infocriteria(Heat_Model1)[2]
BIC[2,2]  <-  infocriteria(Heat_Model2)[2]
BIC[3,2]  <-  infocriteria(Heat_Model3)[2]
BIC[4,2]  <-  infocriteria(Heat_Model4)[2]
BIC[5,2]  <-  infocriteria(Heat_Model5)[2]
BIC[6,2]  <-  infocriteria(Heat_Model6)[2]
BIC[7,2]  <-  infocriteria(Heat_Model7)[2]
BIC[8,2]  <-  infocriteria(Heat_Model8)[2]
BIC[9,2]  <-  infocriteria(Heat_Model9)[2]
BIC[10,2] <- infocriteria(Heat_Model10)[2]
BIC[11,2] <- infocriteria(Heat_Model11)[2]
BIC[12,2] <- infocriteria(Heat_Model12)[2]

AIC  <- matrix(data=0,nrow=12,ncol=2)
colnames(AIC) <- c("Gas","Heat")
rownames(AIC) <- c("GarchSST","GarchSGED","GarchJSU","GarchNorm","eGarchSST","eGarchSGED","eGarchJSU","eGarchNorm","gjrGarchSST","gjrGarchSGED","gjrGarchJSU","gjrGarchNorm")

AIC[1,1]  <-  infocriteria(Gas_Model1)[1]
AIC[2,1]  <-  infocriteria(Gas_Model2)[1]
AIC[3,1]  <-  infocriteria(Gas_Model3)[1]
AIC[4,1]  <-  infocriteria(Gas_Model4)[1]
AIC[5,1]  <-  infocriteria(Gas_Model5)[1]
AIC[6,1]  <-  infocriteria(Gas_Model6)[1]
AIC[7,1]  <-  infocriteria(Gas_Model7)[1]
AIC[8,1]  <-  infocriteria(Gas_Model8)[1]
AIC[9,1]  <-  infocriteria(Gas_Model9)[1]
AIC[10,1] <- infocriteria(Gas_Model10)[1]
AIC[11,1] <- infocriteria(Gas_Model11)[1]
AIC[12,1] <- infocriteria(Gas_Model12)[1]

AIC[1,2]  <-  infocriteria(Heat_Model1)[1]
AIC[2,2]  <-  infocriteria(Heat_Model2)[1]
AIC[3,2]  <-  infocriteria(Heat_Model3)[1]
AIC[4,2]  <-  infocriteria(Heat_Model4)[1]
AIC[5,2]  <-  infocriteria(Heat_Model5)[1]
AIC[6,2]  <-  infocriteria(Heat_Model6)[1]
AIC[7,2]  <-  infocriteria(Heat_Model7)[1]
AIC[8,2]  <-  infocriteria(Heat_Model8)[1]
AIC[9,2]  <-  infocriteria(Heat_Model9)[1]
AIC[10,2] <- infocriteria(Heat_Model10)[1]
AIC[11,2] <- infocriteria(Heat_Model11)[1]
AIC[12,2] <- infocriteria(Heat_Model12)[1]

## Conclusion: eGarch with SST and JSU are the best models #####
## Model5 for SST and Model7 for JSU, Model8 for Normal required to calculate Pearson 

### Store the coefficients of the models

Heat_Coef5   <- round(coef(Heat_Model5),5)
Heat_Coef6   <- round(coef(Heat_Model6),5)
Heat_Coef7   <- round(coef(Heat_Model7),5)
Heat_Coef8   <- round(coef(Heat_Model8),5)
Heat_Coef_Matrix   <- list(Heat_Coef5,Heat_Coef6,Heat_Coef7,Heat_Coef8)
names(Heat_Coef_Matrix) <- c("SST","SGED","JSU","Normal")

Heat_Robust_t5   <- round(Heat_Model5@fit$robust.tval,2)
Heat_Robust_t6   <- round(Heat_Model6@fit$robust.tval,2)
Heat_Robust_t7   <- round(Heat_Model7@fit$robust.tval,2)
Heat_Robust_t8   <- round(Heat_Model8@fit$robust.tval,2)

Heat_Robust_t   <- list(Heat_Robust_t5,Heat_Robust_t6,Heat_Robust_t7,Heat_Robust_t8)
names(Heat_Robust_t) <- c("SST","SGED","JSU","Normal")

Gas_Coef5   <- round(coef(Gas_Model5),5)
Gas_Coef6   <- round(coef(Gas_Model6),5)
Gas_Coef7   <- round(coef(Gas_Model7),5)
Gas_Coef8   <- round(coef(Gas_Model8),5)
Gas_Coef_Matrix   <- list(Gas_Coef5,Gas_Coef6,Gas_Coef7,Gas_Coef8)
names(Gas_Coef_Matrix) <- c("SST","SGED","JSU","Normal")

Gas_Robust_t5   <- round(Gas_Model5@fit$robust.tval,2)
Gas_Robust_t6   <- round(Gas_Model6@fit$robust.tval,2)
Gas_Robust_t7   <- round(Gas_Model7@fit$robust.tval,2)
Gas_Robust_t8   <- round(Gas_Model8@fit$robust.tval,2)

Gas_Robust_t   <- list(Gas_Robust_t5,Gas_Robust_t6,Gas_Robust_t7,Gas_Robust_t8)
names(Gas_Robust_t) <- c("SST","SGED","JSU","Normal")

### Ljung-box test statistic

Heat5_LQ1  <- Box.test(Heat_Model5@fit$residuals/Heat_Model5@fit$sigma,lag=1,type="Lj")
Heat5_LQ5  <- Box.test(Heat_Model5@fit$residuals/Heat_Model5@fit$sigma,lag=5,type="Lj")
Heat5_LQ10 <- Box.test(Heat_Model5@fit$residuals/Heat_Model5@fit$sigma,lag=10,type="Lj")

Heat6_LQ1  <- Box.test(Heat_Model6@fit$residuals/Heat_Model6@fit$sigma,lag=1,type="Lj")
Heat6_LQ5  <- Box.test(Heat_Model6@fit$residuals/Heat_Model6@fit$sigma,lag=5,type="Lj")
Heat6_LQ10 <- Box.test(Heat_Model6@fit$residuals/Heat_Model6@fit$sigma,lag=10,type="Lj")

Heat7_LQ1  <- Box.test(Heat_Model7@fit$residuals/Heat_Model7@fit$sigma,lag=1,type="Lj")
Heat7_LQ5  <- Box.test(Heat_Model7@fit$residuals/Heat_Model7@fit$sigma,lag=5,type="Lj")
Heat7_LQ10 <- Box.test(Heat_Model7@fit$residuals/Heat_Model7@fit$sigma,lag=10,type="Lj")

Heat8_LQ1  <- Box.test(Heat_Model8@fit$residuals/Heat_Model8@fit$sigma,lag=1,type="Lj")
Heat8_LQ5  <- Box.test(Heat_Model8@fit$residuals/Heat_Model8@fit$sigma,lag=5,type="Lj")
Heat8_LQ10 <- Box.test(Heat_Model8@fit$residuals/Heat_Model8@fit$sigma,lag=10,type="Lj")

Heat5_LQ <- c(Heat5_LQ1$statistic,Heat5_LQ5$statistic,Heat5_LQ10$statistic)
Heat6_LQ <- c(Heat6_LQ1$statistic,Heat6_LQ5$statistic,Heat6_LQ10$statistic)
Heat7_LQ <- c(Heat7_LQ1$statistic,Heat7_LQ5$statistic,Heat7_LQ10$statistic)
Heat8_LQ <- c(Heat8_LQ1$statistic,Heat8_LQ5$statistic,Heat8_LQ10$statistic)

Gas5_LQ1  <- Box.test(Gas_Model5@fit$residuals/Gas_Model5@fit$sigma,lag=1,type="Lj")
Gas5_LQ5  <- Box.test(Gas_Model5@fit$residuals/Gas_Model5@fit$sigma,lag=5,type="Lj")
Gas5_LQ10 <- Box.test(Gas_Model5@fit$residuals/Gas_Model5@fit$sigma,lag=10,type="Lj")

Gas6_LQ1  <- Box.test(Gas_Model6@fit$residuals/Gas_Model6@fit$sigma,lag=1,type="Lj")
Gas6_LQ5  <- Box.test(Gas_Model6@fit$residuals/Gas_Model6@fit$sigma,lag=5,type="Lj")
Gas6_LQ10 <- Box.test(Gas_Model6@fit$residuals/Gas_Model6@fit$sigma,lag=10,type="Lj")

Gas7_LQ1  <- Box.test(Gas_Model7@fit$residuals/Gas_Model7@fit$sigma,lag=1,type="Lj")
Gas7_LQ5  <- Box.test(Gas_Model7@fit$residuals/Gas_Model7@fit$sigma,lag=5,type="Lj")
Gas7_LQ10 <- Box.test(Gas_Model7@fit$residuals/Gas_Model7@fit$sigma,lag=10,type="Lj")

Gas8_LQ1  <- Box.test(Gas_Model8@fit$residuals/Gas_Model8@fit$sigma,lag=1,type="Lj")
Gas8_LQ5  <- Box.test(Gas_Model8@fit$residuals/Gas_Model8@fit$sigma,lag=5,type="Lj")
Gas8_LQ10 <- Box.test(Gas_Model8@fit$residuals/Gas_Model8@fit$sigma,lag=10,type="Lj")

Gas5_LQ <- c(Gas5_LQ1$statistic,Gas5_LQ5$statistic,Gas5_LQ10$statistic)
Gas6_LQ <- c(Gas6_LQ1$statistic,Gas6_LQ5$statistic,Gas6_LQ10$statistic)
Gas7_LQ <- c(Gas7_LQ1$statistic,Gas7_LQ5$statistic,Gas7_LQ10$statistic)
Gas8_LQ <- c(Gas8_LQ1$statistic,Gas8_LQ5$statistic,Gas8_LQ10$statistic)

Heat_res_stat <-   cbind(Heat5_LQ,Heat6_LQ,Heat7_LQ,Heat8_LQ)
rownames(Heat_res_stat) <- c("LQ1","LQ5","LQ10")
colnames(Heat_res_stat) <- c("SST","SGED","JSU","Normal")

Gas_res_stat <-   cbind(Gas5_LQ,Gas6_LQ,Gas7_LQ,Gas8_LQ)
rownames(Gas_res_stat) <- c("LQ1","LQ5","LQ10")
colnames(Gas_res_stat) <- c("SST","SGED","JSU","Normal")

#### 1-day ahead Rolling Forecasts ####
## Full sample for roll i.e. last 1000 points for forecast and rest for model

Heat_Roll5 <-  ugarchroll(Heat_Spec5,Heat_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll6 <-  ugarchroll(Heat_Spec6,Heat_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll7 <-  ugarchroll(Heat_Spec7,Heat_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll8 <-  ugarchroll(Heat_Spec8,Heat_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Gas_Roll5 <-  ugarchroll(Gas_Spec5,Gas_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll6 <-  ugarchroll(Gas_Spec6,Gas_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll7 <-  ugarchroll(Gas_Spec7,Gas_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll8 <-  ugarchroll(Gas_Spec8,Gas_xts,forecast.length = 1000,refit.every = 50,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


#### Store the VaR for the out sample using rolling window forecasts

Heat_VaR_SST   <- Heat_Roll5@forecast$VaR
Heat_VaR_SGED  <- Heat_Roll6@forecast$VaR
Heat_VaR_JSU   <- Heat_Roll7@forecast$VaR
Heat_VaR_Norm  <- Heat_Roll8@forecast$VaR

Gas_VaR_SST   <- Gas_Roll5@forecast$VaR
Gas_VaR_SGED  <- Gas_Roll6@forecast$VaR
Gas_VaR_JSU   <- Gas_Roll7@forecast$VaR
Gas_VaR_Norm  <- Gas_Roll8@forecast$VaR

## Significance Levels

signific_levels <- c(0.01,0.05,0.10)

### Backtest the results

Heat_Test_SST_1  <- BacktestVaR(Heat_VaR_SST[,4],Heat_VaR_SST[,1],signific_levels[1])
Heat_Test_SST_5  <- BacktestVaR(Heat_VaR_SST[,4],Heat_VaR_SST[,2],signific_levels[2])
Heat_Test_SST_10 <- BacktestVaR(Heat_VaR_SST[,4],Heat_VaR_SST[,3],signific_levels[3])

Heat_Test_SGED_1  <- BacktestVaR(Heat_VaR_SGED[,4],Heat_VaR_SGED[,1],signific_levels[1])
Heat_Test_SGED_5  <- BacktestVaR(Heat_VaR_SGED[,4],Heat_VaR_SGED[,2],signific_levels[2])
Heat_Test_SGED_10 <- BacktestVaR(Heat_VaR_SGED[,4],Heat_VaR_SGED[,3],signific_levels[3])

Heat_Test_JSU_1  <- BacktestVaR(Heat_VaR_JSU[,4],Heat_VaR_JSU[,1],signific_levels[1])
Heat_Test_JSU_5  <- BacktestVaR(Heat_VaR_JSU[,4],Heat_VaR_JSU[,2],signific_levels[2])
Heat_Test_JSU_10 <- BacktestVaR(Heat_VaR_JSU[,4],Heat_VaR_JSU[,3],signific_levels[3])

Heat_Test_Norm_1  <- BacktestVaR(Heat_VaR_Norm[,4],Heat_VaR_Norm[,1],signific_levels[1])
Heat_Test_Norm_5  <- BacktestVaR(Heat_VaR_Norm[,4],Heat_VaR_Norm[,2],signific_levels[2])
Heat_Test_Norm_10 <- BacktestVaR(Heat_VaR_Norm[,4],Heat_VaR_Norm[,3],signific_levels[3])

Gas_Test_SST_1  <- BacktestVaR(Gas_VaR_SST[,4],Gas_VaR_SST[,1],signific_levels[1])
Gas_Test_SST_5  <- BacktestVaR(Gas_VaR_SST[,4],Gas_VaR_SST[,2],signific_levels[2])
Gas_Test_SST_10 <- BacktestVaR(Gas_VaR_SST[,4],Gas_VaR_SST[,3],signific_levels[3])

Gas_Test_SGED_1  <- BacktestVaR(Gas_VaR_SGED[,4],Gas_VaR_SGED[,1],signific_levels[1])
Gas_Test_SGED_5  <- BacktestVaR(Gas_VaR_SGED[,4],Gas_VaR_SGED[,2],signific_levels[2])
Gas_Test_SGED_10 <- BacktestVaR(Gas_VaR_SGED[,4],Gas_VaR_SGED[,3],signific_levels[3])

Gas_Test_JSU_1  <- BacktestVaR(Gas_VaR_JSU[,4],Gas_VaR_JSU[,1],signific_levels[1])
Gas_Test_JSU_5  <- BacktestVaR(Gas_VaR_JSU[,4],Gas_VaR_JSU[,2],signific_levels[2])
Gas_Test_JSU_10 <- BacktestVaR(Gas_VaR_JSU[,4],Gas_VaR_JSU[,3],signific_levels[3])

Gas_Test_Norm_1  <- BacktestVaR(Gas_VaR_Norm[,4],Gas_VaR_Norm[,1],signific_levels[1])
Gas_Test_Norm_5  <- BacktestVaR(Gas_VaR_Norm[,4],Gas_VaR_Norm[,2],signific_levels[2])
Gas_Test_Norm_10 <- BacktestVaR(Gas_VaR_Norm[,4],Gas_VaR_Norm[,3],signific_levels[3])

#### Saving the test stats and pvalues for 1%, 5% and 10% VaR

## 1%
Heat_Backtest_teststat1 <- matrix(c(Heat_Test_SST_1$LRuc[1],Heat_Test_SST_1$LRcc[1],Heat_Test_SST_1$DQ$stat,
                            Heat_Test_SGED_1$LRuc[1],Heat_Test_SGED_1$LRcc[1],Heat_Test_SGED_1$DQ$stat,
                            Heat_Test_JSU_1$LRuc[1],Heat_Test_JSU_1$LRcc[1],Heat_Test_JSU_1$DQ$stat,
                            Heat_Test_Norm_1$LRuc[1],Heat_Test_Norm_1$LRcc[1],Heat_Test_Norm_1$DQ$stat),nrow=4,byrow = TRUE)
rownames(Heat_Backtest_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_Backtest_teststat1)  <- c("Kupiec","Christoffersen","DQ")

Gas_Backtest_teststat1  <- matrix(c(Gas_Test_SST_1$LRuc[1],Gas_Test_SST_1$LRcc[1],Gas_Test_SST_1$DQ$stat,
                            Gas_Test_SGED_1$LRuc[1],Gas_Test_SGED_1$LRcc[1],Gas_Test_SGED_1$DQ$stat,
                            Gas_Test_JSU_1$LRuc[1],Gas_Test_JSU_1$LRcc[1],Gas_Test_JSU_1$DQ$stat,
                            Gas_Test_Norm_1$LRuc[1],Gas_Test_Norm_1$LRcc[1],Gas_Test_Norm_1$DQ$stat),nrow=4,byrow = TRUE)
rownames(Gas_Backtest_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_Backtest_teststat1)  <- c("Kupiec","Christoffersen","DQ")

Heat_Backtest_pvalues1 <- matrix(c(Heat_Test_SST_1$LRuc[2],Heat_Test_SST_1$LRcc[2],Heat_Test_SST_1$DQ$pvalue,
                             Heat_Test_SGED_1$LRuc[2],Heat_Test_SGED_1$LRcc[2],Heat_Test_SGED_1$DQ$pvalue,
                             Heat_Test_JSU_1$LRuc[2],Heat_Test_JSU_1$LRcc[2],Heat_Test_JSU_1$DQ$pvalue,
                             Heat_Test_Norm_1$LRuc[2],Heat_Test_Norm_1$LRcc[2],Heat_Test_Norm_1$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Heat_Backtest_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_Backtest_pvalues1)  <- c("Kupiec","Christoffersen","DQ")

Gas_Backtest_pvalues1  <- matrix(c(Gas_Test_SST_1$LRuc[2],Gas_Test_SST_1$LRcc[2],Gas_Test_SST_1$DQ$pvalue,
                            Gas_Test_SGED_1$LRuc[2],Gas_Test_SGED_1$LRcc[2],Gas_Test_SGED_1$DQ$pvalue,
                            Gas_Test_JSU_1$LRuc[2],Gas_Test_JSU_1$LRcc[2],Gas_Test_JSU_1$DQ$pvalue,
                            Gas_Test_Norm_1$LRuc[2],Gas_Test_Norm_1$LRcc[2],Gas_Test_Norm_1$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Gas_Backtest_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_Backtest_pvalues1)  <- c("Kupiec","Christoffersen","DQ")

### 5%
Heat_Backtest_teststat5 <- matrix(c(Heat_Test_SST_5$LRuc[1],Heat_Test_SST_5$LRcc[1],Heat_Test_SST_5$DQ$stat,
                                     Heat_Test_SGED_5$LRuc[1],Heat_Test_SGED_5$LRcc[1],Heat_Test_SGED_5$DQ$stat,
                                     Heat_Test_JSU_5$LRuc[1],Heat_Test_JSU_5$LRcc[1],Heat_Test_JSU_5$DQ$stat,
                                     Heat_Test_Norm_5$LRuc[1],Heat_Test_Norm_5$LRcc[1],Heat_Test_Norm_5$DQ$stat),nrow=4,byrow = TRUE)
rownames(Heat_Backtest_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_Backtest_teststat5)  <- c("Kupiec","Christoffersen","DQ")

Gas_Backtest_teststat5  <- matrix(c(Gas_Test_SST_5$LRuc[1],Gas_Test_SST_5$LRcc[1],Gas_Test_SST_5$DQ$stat,
                                    Gas_Test_SGED_5$LRuc[1],Gas_Test_SGED_5$LRcc[1],Gas_Test_SGED_5$DQ$stat,
                                    Gas_Test_JSU_5$LRuc[1],Gas_Test_JSU_5$LRcc[1],Gas_Test_JSU_5$DQ$stat,
                                    Gas_Test_Norm_5$LRuc[1],Gas_Test_Norm_5$LRcc[1],Gas_Test_Norm_5$DQ$stat),nrow=4,byrow = TRUE)
rownames(Gas_Backtest_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_Backtest_teststat5)  <- c("Kupiec","Christoffersen","DQ")

Heat_Backtest_pvalues5 <- matrix(c(Heat_Test_SST_5$LRuc[2],Heat_Test_SST_5$LRcc[2],Heat_Test_SST_5$DQ$pvalue,
                                    Heat_Test_SGED_5$LRuc[2],Heat_Test_SGED_5$LRcc[2],Heat_Test_SGED_5$DQ$pvalue,
                                    Heat_Test_JSU_5$LRuc[2],Heat_Test_JSU_5$LRcc[2],Heat_Test_JSU_5$DQ$pvalue,
                                    Heat_Test_Norm_5$LRuc[2],Heat_Test_Norm_5$LRcc[2],Heat_Test_Norm_5$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Heat_Backtest_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_Backtest_pvalues5)  <- c("Kupiec","Christoffersen","DQ")

Gas_Backtest_pvalues5  <- matrix(c(Gas_Test_SST_5$LRuc[2],Gas_Test_SST_5$LRcc[2],Gas_Test_SST_5$DQ$pvalue,
                                   Gas_Test_SGED_5$LRuc[2],Gas_Test_SGED_5$LRcc[2],Gas_Test_SGED_5$DQ$pvalue,
                                   Gas_Test_JSU_5$LRuc[2],Gas_Test_JSU_5$LRcc[2],Gas_Test_JSU_5$DQ$pvalue,
                                   Gas_Test_Norm_5$LRuc[2],Gas_Test_Norm_5$LRcc[2],Gas_Test_Norm_5$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Gas_Backtest_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_Backtest_pvalues5)  <- c("Kupiec","Christoffersen","DQ")

### 10%

Heat_Backtest_teststat10 <- matrix(c(Heat_Test_SST_10$LRuc[1],Heat_Test_SST_10$LRcc[1],Heat_Test_SST_10$DQ$stat,
                                     Heat_Test_SGED_10$LRuc[1],Heat_Test_SGED_10$LRcc[1],Heat_Test_SGED_10$DQ$stat,
                                     Heat_Test_JSU_10$LRuc[1],Heat_Test_JSU_10$LRcc[1],Heat_Test_JSU_10$DQ$stat,
                                     Heat_Test_Norm_10$LRuc[1],Heat_Test_Norm_10$LRcc[1],Heat_Test_Norm_10$DQ$stat),nrow=4,byrow = TRUE)
rownames(Heat_Backtest_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_Backtest_teststat10)  <- c("Kupiec","Christoffersen","DQ")

Gas_Backtest_teststat10  <- matrix(c(Gas_Test_SST_10$LRuc[1],Gas_Test_SST_10$LRcc[1],Gas_Test_SST_10$DQ$stat,
                                    Gas_Test_SGED_10$LRuc[1],Gas_Test_SGED_10$LRcc[1],Gas_Test_SGED_10$DQ$stat,
                                    Gas_Test_JSU_10$LRuc[1],Gas_Test_JSU_10$LRcc[1],Gas_Test_JSU_10$DQ$stat,
                                    Gas_Test_Norm_10$LRuc[1],Gas_Test_Norm_10$LRcc[1],Gas_Test_Norm_10$DQ$stat),nrow=4,byrow = TRUE)
rownames(Gas_Backtest_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_Backtest_teststat10)  <- c("Kupiec","Christoffersen","DQ")

Heat_Backtest_pvalues10 <- matrix(c(Heat_Test_SST_10$LRuc[2],Heat_Test_SST_10$LRcc[2],Heat_Test_SST_10$DQ$pvalue,
                                    Heat_Test_SGED_10$LRuc[2],Heat_Test_SGED_10$LRcc[2],Heat_Test_SGED_10$DQ$pvalue,
                                    Heat_Test_JSU_10$LRuc[2],Heat_Test_JSU_10$LRcc[2],Heat_Test_JSU_10$DQ$pvalue,
                                    Heat_Test_Norm_10$LRuc[2],Heat_Test_Norm_10$LRcc[2],Heat_Test_Norm_10$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Heat_Backtest_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_Backtest_pvalues10)  <- c("Kupiec","Christoffersen","DQ")

Gas_Backtest_pvalues10  <- matrix(c(Gas_Test_SST_10$LRuc[2],Gas_Test_SST_10$LRcc[2],Gas_Test_SST_10$DQ$pvalue,
                                   Gas_Test_SGED_10$LRuc[2],Gas_Test_SGED_10$LRcc[2],Gas_Test_SGED_10$DQ$pvalue,
                                   Gas_Test_JSU_10$LRuc[2],Gas_Test_JSU_10$LRcc[2],Gas_Test_JSU_10$DQ$pvalue,
                                   Gas_Test_Norm_10$LRuc[2],Gas_Test_Norm_10$LRcc[2],Gas_Test_Norm_10$DQ$pvalue),nrow=4,byrow = TRUE)
rownames(Gas_Backtest_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_Backtest_pvalues10)  <- c("Kupiec","Christoffersen","DQ")

Heat_Backtest_pvalues1  <- round(Heat_Backtest_pvalues1,4)
Heat_Backtest_pvalues5  <- round(Heat_Backtest_pvalues5,4)
Heat_Backtest_pvalues10 <- round(Heat_Backtest_pvalues10,4)

Gas_Backtest_pvalues1  <- round(Gas_Backtest_pvalues1,4)
Gas_Backtest_pvalues5  <- round(Gas_Backtest_pvalues5,4)
Gas_Backtest_pvalues10 <- round(Gas_Backtest_pvalues10,4)

Heat_Backtest_teststat1  <- round(Heat_Backtest_teststat1,4)
Heat_Backtest_teststat5  <- round(Heat_Backtest_teststat5,4)
Heat_Backtest_teststat10 <- round(Heat_Backtest_teststat10,4)

Gas_Backtest_teststat1  <- round(Gas_Backtest_teststat1,4)
Gas_Backtest_teststat5  <- round(Gas_Backtest_teststat5,4)
Gas_Backtest_teststat10 <- round(Gas_Backtest_teststat10,4)

### Using VaRTest and VaRDuration Test individually ###

## VaR Test in R

Heat_Norm_VaR_Test1  <- VaRTest(0.01,Heat_VaR_Norm[,4],Heat_VaR_Norm[,1])
Heat_Norm_VaR_Test5  <- VaRTest(0.05,Heat_VaR_Norm[,4],Heat_VaR_Norm[,2])
Heat_Norm_VaR_Test10 <- VaRTest(0.1,Heat_VaR_Norm[,4],Heat_VaR_Norm[,3])

Heat_JSU_VaR_Test1  <- VaRTest(0.01,Heat_VaR_JSU[,4],Heat_VaR_JSU[,1])
Heat_JSU_VaR_Test5  <- VaRTest(0.05,Heat_VaR_JSU[,4],Heat_VaR_JSU[,2])
Heat_JSU_VaR_Test10 <- VaRTest(0.1,Heat_VaR_JSU[,4],Heat_VaR_JSU[,3])

Heat_SGED_VaR_Test1  <- VaRTest(0.01,Heat_VaR_SGED[,4],Heat_VaR_SGED[,1])
Heat_SGED_VaR_Test5  <- VaRTest(0.05,Heat_VaR_SGED[,4],Heat_VaR_SGED[,2])
Heat_SGED_VaR_Test10 <- VaRTest(0.1,Heat_VaR_SGED[,4],Heat_VaR_SGED[,3])

Heat_SST_VaR_Test1  <- VaRTest(0.01,Heat_VaR_SST[,4],Heat_VaR_SST[,1])
Heat_SST_VaR_Test5  <- VaRTest(0.05,Heat_VaR_SST[,4],Heat_VaR_SST[,2])
Heat_SST_VaR_Test10 <- VaRTest(0.1,Heat_VaR_SST[,4],Heat_VaR_SST[,3])

Gas_Norm_VaR_Test1  <- VaRTest(0.01,Gas_VaR_Norm[,4],Gas_VaR_Norm[,1])
Gas_Norm_VaR_Test5  <- VaRTest(0.05,Gas_VaR_Norm[,4],Gas_VaR_Norm[,2])
Gas_Norm_VaR_Test10 <- VaRTest(0.1,Gas_VaR_Norm[,4],Gas_VaR_Norm[,3])

Gas_JSU_VaR_Test1  <- VaRTest(0.01,Gas_VaR_JSU[,4],Gas_VaR_JSU[,1])
Gas_JSU_VaR_Test5  <- VaRTest(0.05,Gas_VaR_JSU[,4],Gas_VaR_JSU[,2])
Gas_JSU_VaR_Test10 <- VaRTest(0.1,Gas_VaR_JSU[,4],Gas_VaR_JSU[,3])

Gas_SGED_VaR_Test1  <- VaRTest(0.01,Gas_VaR_SGED[,4],Gas_VaR_SGED[,1])
Gas_SGED_VaR_Test5  <- VaRTest(0.05,Gas_VaR_SGED[,4],Gas_VaR_SGED[,2])
Gas_SGED_VaR_Test10 <- VaRTest(0.1,Gas_VaR_SGED[,4],Gas_VaR_SGED[,3])

Gas_SST_VaR_Test1  <- VaRTest(0.01,Gas_VaR_SST[,4],Gas_VaR_SST[,1])
Gas_SST_VaR_Test5  <- VaRTest(0.05,Gas_VaR_SST[,4],Gas_VaR_SST[,2])
Gas_SST_VaR_Test10 <- VaRTest(0.1,Gas_VaR_SST[,4],Gas_VaR_SST[,3])

## Store the pvalues at a place

Heat_VaR_Test_pvalues1 <- matrix(c(Heat_SST_VaR_Test1$uc.LRp,Heat_SST_VaR_Test1$cc.LRp,
                                    Heat_SGED_VaR_Test1$uc.LRp,Heat_SGED_VaR_Test1$cc.LRp,
                                    Heat_JSU_VaR_Test1$uc.LRp,Heat_JSU_VaR_Test1$cc.LRp,
                                    Heat_Norm_VaR_Test1$uc.LRp,Heat_Norm_VaR_Test1$cc.LRp),ncol=2,byrow = TRUE)

Heat_VaR_Test_pvalues5 <- matrix(c(Heat_SST_VaR_Test5$uc.LRp,Heat_SST_VaR_Test5$cc.LRp,
                                    Heat_SGED_VaR_Test5$uc.LRp,Heat_SGED_VaR_Test5$cc.LRp,
                                    Heat_JSU_VaR_Test5$uc.LRp,Heat_JSU_VaR_Test5$cc.LRp,
                                    Heat_Norm_VaR_Test5$uc.LRp,Heat_Norm_VaR_Test5$cc.LRp),ncol=2,byrow = TRUE)

Heat_VaR_Test_pvalues10 <- matrix(c(Heat_SST_VaR_Test10$uc.LRp,Heat_SST_VaR_Test10$cc.LRp,
                                    Heat_SGED_VaR_Test10$uc.LRp,Heat_SGED_VaR_Test10$cc.LRp,
                                    Heat_JSU_VaR_Test10$uc.LRp,Heat_JSU_VaR_Test10$cc.LRp,
                                    Heat_Norm_VaR_Test10$uc.LRp,Heat_Norm_VaR_Test10$cc.LRp),ncol=2,byrow = TRUE)

Gas_VaR_Test_pvalues1 <- matrix(c(Gas_SST_VaR_Test1$uc.LRp,Gas_SST_VaR_Test1$cc.LRp,
                                    Gas_SGED_VaR_Test1$uc.LRp,Gas_SGED_VaR_Test1$cc.LRp,
                                    Gas_JSU_VaR_Test1$uc.LRp,Gas_JSU_VaR_Test1$cc.LRp,
                                    Gas_Norm_VaR_Test1$uc.LRp,Gas_Norm_VaR_Test1$cc.LRp),ncol=2,byrow = TRUE)

Gas_VaR_Test_pvalues5 <- matrix(c(Gas_SST_VaR_Test5$uc.LRp,Gas_SST_VaR_Test5$cc.LRp,
                                    Gas_SGED_VaR_Test5$uc.LRp,Gas_SGED_VaR_Test5$cc.LRp,
                                    Gas_JSU_VaR_Test5$uc.LRp,Gas_JSU_VaR_Test5$cc.LRp,
                                    Gas_Norm_VaR_Test5$uc.LRp,Gas_Norm_VaR_Test5$cc.LRp),ncol=2,byrow = TRUE)

Gas_VaR_Test_pvalues10 <- matrix(c(Gas_SST_VaR_Test10$uc.LRp,Gas_SST_VaR_Test10$cc.LRp,
                                     Gas_SGED_VaR_Test10$uc.LRp,Gas_SGED_VaR_Test10$cc.LRp,
                                     Gas_JSU_VaR_Test10$uc.LRp,Gas_JSU_VaR_Test10$cc.LRp,
                                     Gas_Norm_VaR_Test10$uc.LRp,Gas_Norm_VaR_Test10$cc.LRp),ncol=2,byrow = TRUE)

## Storing the test stats
Heat_VaR_Test_teststat1 <- matrix(c(Heat_SST_VaR_Test1$uc.LRstat,Heat_SST_VaR_Test1$cc.LRstat,
                                    Heat_SGED_VaR_Test1$uc.LRstat,Heat_SGED_VaR_Test1$cc.LRstat,
                                    Heat_JSU_VaR_Test1$uc.LRstat,Heat_JSU_VaR_Test1$cc.LRstat,
                                    Heat_Norm_VaR_Test1$uc.LRstat,Heat_Norm_VaR_Test1$cc.LRstat),ncol=2,byrow = TRUE)

Heat_VaR_Test_teststat5 <- matrix(c(Heat_SST_VaR_Test5$uc.LRstat,Heat_SST_VaR_Test5$cc.LRstat,
                                    Heat_SGED_VaR_Test5$uc.LRstat,Heat_SGED_VaR_Test5$cc.LRstat,
                                    Heat_JSU_VaR_Test5$uc.LRstat,Heat_JSU_VaR_Test5$cc.LRstat,
                                    Heat_Norm_VaR_Test5$uc.LRstat,Heat_Norm_VaR_Test5$cc.LRstat),ncol=2,byrow = TRUE)

Heat_VaR_Test_teststat10 <- matrix(c(Heat_SST_VaR_Test10$uc.LRstat,Heat_SST_VaR_Test10$cc.LRstat,
                                     Heat_SGED_VaR_Test10$uc.LRstat,Heat_SGED_VaR_Test10$cc.LRstat,
                                     Heat_JSU_VaR_Test10$uc.LRstat,Heat_JSU_VaR_Test10$cc.LRstat,
                                     Heat_Norm_VaR_Test10$uc.LRstat,Heat_Norm_VaR_Test10$cc.LRstat),ncol=2,byrow = TRUE)

Gas_VaR_Test_teststat1 <- matrix(c(Gas_SST_VaR_Test1$uc.LRstat,Gas_SST_VaR_Test1$cc.LRstat,
                                  Gas_SGED_VaR_Test1$uc.LRstat,Gas_SGED_VaR_Test1$cc.LRstat,
                                  Gas_JSU_VaR_Test1$uc.LRstat,Gas_JSU_VaR_Test1$cc.LRstat,
                                  Gas_Norm_VaR_Test1$uc.LRstat,Gas_Norm_VaR_Test1$cc.LRstat),ncol=2,byrow = TRUE)

Gas_VaR_Test_teststat5 <- matrix(c(Gas_SST_VaR_Test5$uc.LRstat,Gas_SST_VaR_Test5$cc.LRstat,
                                  Gas_SGED_VaR_Test5$uc.LRstat,Gas_SGED_VaR_Test5$cc.LRstat,
                                  Gas_JSU_VaR_Test5$uc.LRstat,Gas_JSU_VaR_Test5$cc.LRstat,
                                  Gas_Norm_VaR_Test5$uc.LRstat,Gas_Norm_VaR_Test5$cc.LRstat),ncol=2,byrow = TRUE)

Gas_VaR_Test_teststat10 <- matrix(c(Gas_SST_VaR_Test10$uc.LRstat,Gas_SST_VaR_Test10$cc.LRstat,
                                   Gas_SGED_VaR_Test10$uc.LRstat,Gas_SGED_VaR_Test10$cc.LRstat,
                                   Gas_JSU_VaR_Test10$uc.LRstat,Gas_JSU_VaR_Test10$cc.LRstat,
                                   Gas_Norm_VaR_Test10$uc.LRstat,Gas_Norm_VaR_Test10$cc.LRstat),ncol=2,byrow = TRUE)
## Round them up to 4 decimals

Heat_VaR_Test_pvalues1  <- round(Heat_VaR_Test_pvalues1,4)
Heat_VaR_Test_pvalues5  <- round(Heat_VaR_Test_pvalues5,4)
Heat_VaR_Test_pvalues10 <- round(Heat_VaR_Test_pvalues10,4)

Gas_VaR_Test_pvalues1  <- round(Gas_VaR_Test_pvalues1,4)
Gas_VaR_Test_pvalues5  <- round(Gas_VaR_Test_pvalues5,4)
Gas_VaR_Test_pvalues10 <- round(Gas_VaR_Test_pvalues10,4)

Heat_VaR_Test_teststat1  <- round(Heat_VaR_Test_teststat1,4)
Heat_VaR_Test_teststat5  <- round(Heat_VaR_Test_teststat5,4)
Heat_VaR_Test_teststat10 <- round(Heat_VaR_Test_teststat10,4)

Gas_VaR_Test_teststat1  <- round(Gas_VaR_Test_teststat1,4)
Gas_VaR_Test_teststat5  <- round(Gas_VaR_Test_teststat5,4)
Gas_VaR_Test_teststat10 <- round(Gas_VaR_Test_teststat10,4)

rownames(Heat_VaR_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Test_pvalues1)  <- c("Kupiec","Christoffersen")
rownames(Heat_VaR_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Test_pvalues5)  <- c("Kupiec","Christoffersen")
rownames(Heat_VaR_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Test_pvalues10)  <- c("Kupiec","Christoffersen")

rownames(Gas_VaR_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Test_pvalues1)  <- c("Kupiec","Christoffersen")
rownames(Gas_VaR_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Test_pvalues5)  <- c("Kupiec","Christoffersen")
rownames(Gas_VaR_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Test_pvalues10)  <- c("Kupiec","Christoffersen")

rownames(Heat_VaR_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Test_teststat1)  <- c("Kupiec","Christoffersen")
rownames(Heat_VaR_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Test_teststat5)  <- c("Kupiec","Christoffersen")
rownames(Heat_VaR_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Test_teststat10)  <- c("Kupiec","Christoffersen")

rownames(Gas_VaR_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Test_teststat1)  <- c("Kupiec","Christoffersen")
rownames(Gas_VaR_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Test_teststat5)  <- c("Kupiec","Christoffersen")
rownames(Gas_VaR_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Test_teststat10)  <- c("Kupiec","Christoffersen")


## VaR Duration test

Heat_Norm_VaR_Dur_Test1  <- VaRDurTest(0.01,Heat_VaR_Norm[,4],Heat_VaR_Norm[,1])
Heat_Norm_VaR_Dur_Test5  <- VaRDurTest(0.05,Heat_VaR_Norm[,4],Heat_VaR_Norm[,2])
Heat_Norm_VaR_Dur_Test10 <- VaRDurTest(0.1,Heat_VaR_Norm[,4],Heat_VaR_Norm[,3])

Heat_JSU_VaR_Dur_Test1  <- VaRDurTest(0.01,Heat_VaR_JSU[,4],Heat_VaR_JSU[,1])
Heat_JSU_VaR_Dur_Test5  <- VaRDurTest(0.05,Heat_VaR_JSU[,4],Heat_VaR_JSU[,2])
Heat_JSU_VaR_Dur_Test10 <- VaRDurTest(0.1,Heat_VaR_JSU[,4],Heat_VaR_JSU[,3])

Heat_SGED_VaR_Dur_Test1  <- VaRDurTest(0.01,Heat_VaR_SGED[,4],Heat_VaR_SGED[,1])
Heat_SGED_VaR_Dur_Test5  <- VaRDurTest(0.05,Heat_VaR_SGED[,4],Heat_VaR_SGED[,2])
Heat_SGED_VaR_Dur_Test10 <- VaRDurTest(0.1,Heat_VaR_SGED[,4],Heat_VaR_SGED[,3])

Heat_SST_VaR_Dur_Test1  <- VaRDurTest(0.01,Heat_VaR_SST[,4],Heat_VaR_SST[,1])
Heat_SST_VaR_Dur_Test5  <- VaRDurTest(0.05,Heat_VaR_SST[,4],Heat_VaR_SST[,2])
Heat_SST_VaR_Dur_Test10 <- VaRDurTest(0.1,Heat_VaR_SST[,4],Heat_VaR_SST[,3])

Gas_Norm_VaR_Dur_Test1  <- VaRDurTest(0.01,Gas_VaR_Norm[,4],Gas_VaR_Norm[,1])
Gas_Norm_VaR_Dur_Test5  <- VaRDurTest(0.05,Gas_VaR_Norm[,4],Gas_VaR_Norm[,2])
Gas_Norm_VaR_Dur_Test10 <- VaRDurTest(0.1,Gas_VaR_Norm[,4],Gas_VaR_Norm[,3])

Gas_JSU_VaR_Dur_Test1  <- VaRDurTest(0.01,Gas_VaR_JSU[,4],Gas_VaR_JSU[,1])
Gas_JSU_VaR_Dur_Test5  <- VaRDurTest(0.05,Gas_VaR_JSU[,4],Gas_VaR_JSU[,2])
Gas_JSU_VaR_Dur_Test10 <- VaRDurTest(0.1,Gas_VaR_JSU[,4],Gas_VaR_JSU[,3])

Gas_SGED_VaR_Dur_Test1  <- VaRDurTest(0.01,Gas_VaR_SGED[,4],Gas_VaR_SGED[,1])
Gas_SGED_VaR_Dur_Test5  <- VaRDurTest(0.05,Gas_VaR_SGED[,4],Gas_VaR_SGED[,2])
Gas_SGED_VaR_Dur_Test10 <- VaRDurTest(0.1,Gas_VaR_SGED[,4],Gas_VaR_SGED[,3])

Gas_SST_VaR_Dur_Test1  <- VaRDurTest(0.01,Gas_VaR_SST[,4],Gas_VaR_SST[,1])
Gas_SST_VaR_Dur_Test5  <- VaRDurTest(0.05,Gas_VaR_SST[,4],Gas_VaR_SST[,2])
Gas_SST_VaR_Dur_Test10 <- VaRDurTest(0.1,Gas_VaR_SST[,4],Gas_VaR_SST[,3])

### Store the test results and p=values
Heat_VaR_Dur_Test_pvalues1 <- matrix(c(Heat_SST_VaR_Dur_Test1$LRp,Heat_SGED_VaR_Dur_Test1$LRp,
                                        Heat_JSU_VaR_Dur_Test1$LRp,Heat_Norm_VaR_Dur_Test1$LRp),nrow=4,byrow = TRUE)
Heat_VaR_Dur_Test_pvalues5 <- matrix(c(Heat_SST_VaR_Dur_Test5$LRp,Heat_SGED_VaR_Dur_Test5$LRp,
                                        Heat_JSU_VaR_Dur_Test5$LRp,Heat_Norm_VaR_Dur_Test5$LRp),nrow=4,byrow = TRUE)
Heat_VaR_Dur_Test_pvalues10 <- matrix(c(Heat_SST_VaR_Dur_Test10$LRp,Heat_SGED_VaR_Dur_Test10$LRp,
                                        Heat_JSU_VaR_Dur_Test10$LRp,Heat_Norm_VaR_Dur_Test10$LRp),nrow=4,byrow = TRUE)

Gas_VaR_Dur_Test_pvalues1 <- matrix(c(Gas_SST_VaR_Dur_Test1$LRp,Gas_SGED_VaR_Dur_Test1$LRp,
                                        Gas_JSU_VaR_Dur_Test1$LRp,Gas_Norm_VaR_Dur_Test1$LRp),nrow=4,byrow = TRUE)
Gas_VaR_Dur_Test_pvalues5 <- matrix(c(Gas_SST_VaR_Dur_Test5$LRp,Gas_SGED_VaR_Dur_Test5$LRp,
                                        Gas_JSU_VaR_Dur_Test5$LRp,Gas_Norm_VaR_Dur_Test5$LRp),nrow=4,byrow = TRUE)
Gas_VaR_Dur_Test_pvalues10 <- matrix(c(Gas_SST_VaR_Dur_Test10$LRp,Gas_SGED_VaR_Dur_Test10$LRp,
                                         Gas_JSU_VaR_Dur_Test10$LRp,Gas_Norm_VaR_Dur_Test10$LRp),nrow=4,byrow = TRUE)

Heat_VaR_Dur_Test_teststat1 <- matrix(c(Heat_SST_VaR_Dur_Test1$uLL,Heat_SST_VaR_Dur_Test1$rLL,
                                         Heat_SGED_VaR_Dur_Test1$uLL,Heat_SGED_VaR_Dur_Test1$rLL,
                                        Heat_JSU_VaR_Dur_Test1$uLL,Heat_JSU_VaR_Dur_Test1$rLL,
                                        Heat_Norm_VaR_Dur_Test1$uLL,Heat_Norm_VaR_Dur_Test1$rLL),nrow=4,byrow = TRUE)
Heat_VaR_Dur_Test_teststat5 <- matrix(c(Heat_SST_VaR_Dur_Test5$uLL,Heat_SST_VaR_Dur_Test5$rLL,
                                         Heat_SGED_VaR_Dur_Test5$uLL,Heat_SGED_VaR_Dur_Test5$rLL,
                                         Heat_JSU_VaR_Dur_Test5$uLL,Heat_JSU_VaR_Dur_Test5$rLL,
                                         Heat_Norm_VaR_Dur_Test5$uLL,Heat_Norm_VaR_Dur_Test5$rLL),nrow=4,byrow = TRUE)
Heat_VaR_Dur_Test_teststat10 <- matrix(c(Heat_SST_VaR_Dur_Test10$uLL,Heat_SST_VaR_Dur_Test10$rLL,
                                         Heat_SGED_VaR_Dur_Test10$uLL,Heat_SGED_VaR_Dur_Test10$rLL,
                                         Heat_JSU_VaR_Dur_Test10$uLL,Heat_JSU_VaR_Dur_Test10$rLL,
                                         Heat_Norm_VaR_Dur_Test10$uLL,Heat_Norm_VaR_Dur_Test10$rLL),nrow=4,byrow = TRUE)

Gas_VaR_Dur_Test_teststat1 <- matrix(c(Gas_SST_VaR_Dur_Test1$uLL,Gas_SST_VaR_Dur_Test1$rLL,
                                         Gas_SGED_VaR_Dur_Test1$uLL,Gas_SGED_VaR_Dur_Test1$rLL,
                                         Gas_JSU_VaR_Dur_Test1$uLL,Gas_JSU_VaR_Dur_Test1$rLL,
                                         Gas_Norm_VaR_Dur_Test1$uLL,Gas_Norm_VaR_Dur_Test1$rLL),nrow=4,byrow = TRUE)
Gas_VaR_Dur_Test_teststat5 <- matrix(c(Gas_SST_VaR_Dur_Test5$uLL,Gas_SST_VaR_Dur_Test5$rLL,
                                         Gas_SGED_VaR_Dur_Test5$uLL,Gas_SGED_VaR_Dur_Test5$rLL,
                                         Gas_JSU_VaR_Dur_Test5$uLL,Gas_JSU_VaR_Dur_Test5$rLL,
                                         Gas_Norm_VaR_Dur_Test5$uLL,Gas_Norm_VaR_Dur_Test5$rLL),nrow=4,byrow = TRUE)
Gas_VaR_Dur_Test_teststat10 <- matrix(c(Gas_SST_VaR_Dur_Test10$uLL,Gas_SST_VaR_Dur_Test10$rLL,
                                          Gas_SGED_VaR_Dur_Test10$uLL,Gas_SGED_VaR_Dur_Test10$rLL,
                                          Gas_JSU_VaR_Dur_Test10$uLL,Gas_JSU_VaR_Dur_Test10$rLL,
                                          Gas_Norm_VaR_Dur_Test10$uLL,Gas_Norm_VaR_Dur_Test10$rLL),nrow=4,byrow = TRUE)

## Round them up to 4 decimals

Heat_VaR_Dur_Test_pvalues1  <- round(Heat_VaR_Dur_Test_pvalues1,4)
Heat_VaR_Dur_Test_pvalues5  <- round(Heat_VaR_Dur_Test_pvalues5,4)
Heat_VaR_Dur_Test_pvalues10 <- round(Heat_VaR_Dur_Test_pvalues10,4)

Gas_VaR_Dur_Test_pvalues1  <- round(Gas_VaR_Dur_Test_pvalues1,4)
Gas_VaR_Dur_Test_pvalues5  <- round(Gas_VaR_Dur_Test_pvalues5,4)
Gas_VaR_Dur_Test_pvalues10 <- round(Gas_VaR_Dur_Test_pvalues10,4)

Heat_VaR_Dur_Test_teststat1  <- round(Heat_VaR_Dur_Test_teststat1,4)
Heat_VaR_Dur_Test_teststat5  <- round(Heat_VaR_Dur_Test_teststat5,4)
Heat_VaR_Dur_Test_teststat10 <- round(Heat_VaR_Dur_Test_teststat10,4)

Gas_VaR_Dur_Test_teststat1  <- round(Gas_VaR_Dur_Test_teststat1,4)
Gas_VaR_Dur_Test_teststat5  <- round(Gas_VaR_Dur_Test_teststat5,4)
Gas_VaR_Dur_Test_teststat10 <- round(Gas_VaR_Dur_Test_teststat10,4)

rownames(Heat_VaR_Dur_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
rownames(Heat_VaR_Dur_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
rownames(Heat_VaR_Dur_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")

rownames(Gas_VaR_Dur_Test_pvalues1)  <- c("SST","SGED","JSU","Norm")
rownames(Gas_VaR_Dur_Test_pvalues5)  <- c("SST","SGED","JSU","Norm")
rownames(Gas_VaR_Dur_Test_pvalues10)  <- c("SST","SGED","JSU","Norm")

rownames(Heat_VaR_Dur_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Dur_Test_teststat1)  <- c("Unrestricted","Restricted")
rownames(Heat_VaR_Dur_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Dur_Test_teststat5)  <- c("Unrestricted","Restricted")
rownames(Heat_VaR_Dur_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Heat_VaR_Dur_Test_teststat10)  <- c("Unrestricted","Restricted")

rownames(Gas_VaR_Dur_Test_teststat1)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Dur_Test_teststat1)  <- c("Unrestricted","Restricted")
rownames(Gas_VaR_Dur_Test_teststat5)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Dur_Test_teststat5)  <- c("Unrestricted","Restricted")
rownames(Gas_VaR_Dur_Test_teststat10)  <- c("SST","SGED","JSU","Norm")
colnames(Gas_VaR_Dur_Test_teststat10)  <- c("Unrestricted","Restricted")


## QL ratio and FZL ratio

QL_SGED_1  <- round(Gas_Test_SGED_1$Loss$Loss/Gas_Test_Norm_1$Loss$Loss,4)
QL_SGED_5  <- round(Gas_Test_SGED_5$Loss$Loss/Gas_Test_Norm_5$Loss$Loss,4)
QL_SGED_10 <- round(Gas_Test_SGED_10$Loss$Loss/Gas_Test_Norm_10$Loss$Loss,4)

QL_SST_1  <- round(Gas_Test_SST_1$Loss$Loss/Gas_Test_Norm_1$Loss$Loss,4)
QL_SST_5  <- round(Gas_Test_SST_5$Loss$Loss/Gas_Test_Norm_5$Loss$Loss,4)
QL_SST_10 <- round(Gas_Test_SST_10$Loss$Loss/Gas_Test_Norm_10$Loss$Loss,4)

QL_JSU_1  <- round(Gas_Test_JSU_1$Loss$Loss/Gas_Test_Norm_1$Loss$Loss,4)
QL_JSU_5  <- round(Gas_Test_JSU_5$Loss$Loss/Gas_Test_Norm_5$Loss$Loss,4)
QL_JSU_10 <- round(Gas_Test_JSU_10$Loss$Loss/Gas_Test_Norm_10$Loss$Loss,4)

QL <- matrix(c(QL_SST_1,QL_SST_5,QL_SST_10,QL_SGED_1,QL_SGED_5,QL_SGED_10,QL_JSU_1,QL_JSU_5,QL_JSU_10),nrow=3,byrow=TRUE)
rownames(QL) <- c("SST","SGED","JSU")
colnames(QL) <- c("1%","5%","10%")

AE_SST_1  <- round(Gas_Test_SST_1$AE,4)
AE_SST_5  <- round(Gas_Test_SST_5$AE,4)
AE_SST_10 <- round(Gas_Test_SST_10$AE,4)

AE_SGED_1  <- round(Gas_Test_SGED_1$AE,4)
AE_SGED_5  <- round(Gas_Test_SGED_5$AE,4)
AE_SGED_10 <- round(Gas_Test_SGED_10$AE,4)

AE_JSU_1  <- round(Gas_Test_JSU_1$AE,4)
AE_JSU_5  <- round(Gas_Test_JSU_5$AE,4)
AE_JSU_10 <- round(Gas_Test_JSU_10$AE,4)

AE_Norm_1  <- round(Gas_Test_Norm_1$AE,4)
AE_Norm_5  <- round(Gas_Test_Norm_5$AE,4)
AE_Norm_10 <- round(Gas_Test_Norm_10$AE,4)

AE <- matrix(c(AE_SST_1,AE_SST_5,AE_SST_10,AE_SGED_1,AE_SGED_5,AE_SGED_10,AE_JSU_1,AE_JSU_5,AE_JSU_10,AE_Norm_1,AE_Norm_5,AE_Norm_10),nrow=4,byrow=TRUE)
rownames(AE) <- c("SST","SGED","JSU","Norm")
colnames(AE) <- c("1%","5%","10%")

### Check the summary stats of the standardized residuals

Heat_Res8 <- Heat_Model8@fit$residuals/Heat_Model8@fit$sigma
Gas_Res8   <- Gas_Model8@fit$residuals/Gas_Model8@fit$sigma

Heat_Descr8 <- basicStats(Heat_Res8)
Gas_Descr8 <- basicStats(Gas_Res8)

### Fitting distributions on the standardized residuals 

Heat_Dist8 <- pearsonFitML(Heat_Res8)
Gas_Dist8   <- pearsonFitML(Gas_Res8)
Heat_JSU_param <- JohnsonFit(Heat_Res8,moment ="find")
Gas_JSU_param <- JohnsonFit(Gas_Res8,moment ="find")

### Store the parameters for different distributions  

Heat_P4_Param  <- Heat_Dist8[2:5]
Gas_P4_Param    <- Gas_Dist8[2:5]

### Checking the fit using KS tests

Heat_test1 <- ks.test(Heat_Res8,rpearson(length(Heat_Res8),params = Heat_Dist8))
Gas_test1 <- ks.test(Gas_Res8,rpearson(length(Gas_Res8),params = Gas_Dist8))


### Calculate the quantiles

Heat_Exceed_Pearson  <- array(data=0,dim=length(signific_levels))
Gas_Exceed_Pearson    <- array(data=0,dim=length(signific_levels))

for(j in 1:3)
{
  Z_Pearson_Heat <- qpearson(signific_levels[j],Heat_Dist8)
  Z_Pearson_Gas   <- qpearson(signific_levels[j],Gas_Dist8)
  
  mu_Heat   <- array(data=0,dim=length(Returns_Heat)-1000)
  VaR_Heat   <- array(data=0,dim=length(Returns_Heat)-1001)
  mu_Gas   <- array(data=0,dim=length(Returns_Gas)-1000)
  VaR_Gas   <- array(data=0,dim=length(Returns_Gas)-1001)
  
  for(i in 2:(length(Returns_Heat)-1000))
  {
    mu_Heat[i]    <- Heat_Model8@fit$coef[1]+Heat_Model8@fit$coef[2]*Returns_Heat[i-1]
    VaR_Heat[i-1] <- mu_Heat[i]+Heat_Model8@fit$sigma[i]*Z_Pearson_Heat
  }

  for(i in 2:(length(Returns_Gas)-1000))
  {
    mu_Gas[i] <- Gas_Model8@fit$coef[1]+Gas_Model8@fit$coef[2]*Returns_Gas[i-1]
    VaR_Gas[i-1] <- mu_Gas[i]+Gas_Model8@fit$sigma[i]*Z_Pearson_Gas
  }
  
  Heat_Violations <- VaRTest(signific_levels[j],Returns_Heat[1:(length(VaR_Heat))],VaR_Heat)
  Gas_Violations   <- VaRTest(signific_levels[j],Returns_Gas[1:(length(VaR_Gas))],VaR_Gas)
  Heat_Exceed_Pearson[j] <- Heat_Violations$actual.exceed
  Gas_Exceed_Pearson[j] <- Gas_Violations$actual.exceed
  
}


#### Out of sample VaR forecast  #####

#### Using Hold out sample

Holdout_Heat_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Gas_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Heat_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Gas_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Heat_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Gas_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Heat_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Gas_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Heat_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Gas_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Heat_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Gas_Pearson_pvalues) <- c("1%","5%","10%")

## Store the volatility and mu forecast
#signific_levels <- c(0.1,0.05,0.01)

Heat_Holdout_Mu <- as.data.frame(Heat_Roll8)[,'Mu']
Heat_Holdout_Sigma <- as.data.frame(Heat_Roll8)[,'Sigma']
Gas_Holdout_Mu <- as.data.frame(Gas_Roll8)[,'Mu']
Gas_Holdout_Sigma <- as.data.frame(Gas_Roll8)[,'Sigma']

for(j in 1:3)
{
 
  Z_Pearson_Heat <- qpearson(signific_levels[j],Heat_Dist8)
  Z_Pearson_Gas   <- qpearson(signific_levels[j],Gas_Dist8)

  Heat_Holdout_VaR <- Heat_Holdout_Mu + Heat_Holdout_Sigma*Z_Pearson_Heat
  Gas_Holdout_VaR   <- Gas_Holdout_Mu + Gas_Holdout_Sigma*Z_Pearson_Gas
  Heat_Violations <- BacktestVaR(Holdout_Returns_Heat[1:(length(Holdout_Returns_Heat)-1)],Heat_Holdout_VaR[-1],signific_levels[j])
  Gas_Violations   <- BacktestVaR(Holdout_Returns_Gas[1:(length(Holdout_Returns_Gas)-1)],Gas_Holdout_VaR[-1],signific_levels[j])
  Holdout_Heat_Pearson[j,] <- c(Heat_Violations$LRuc[1],Heat_Violations$LRcc[1],Heat_Violations$DQ$stat,Heat_Violations$AE,Heat_Violations$Loss$Loss)
  Holdout_Gas_Pearson[j,]   <- c(Gas_Violations$LRuc[1],Gas_Violations$LRcc[1],Gas_Violations$DQ$stat,Gas_Violations$AE,Gas_Violations$Loss$Loss)
  Holdout_Heat_Pearson_pvalues[j,] <- c(Heat_Violations$LRuc[2],Heat_Violations$LRcc[2],Heat_Violations$DQ$pvalue)
  Holdout_Gas_Pearson_pvalues[j,] <- c(Gas_Violations$LRuc[2],Gas_Violations$LRcc[2],Gas_Violations$DQ$pvalue)
}


#### Subsample Analysis ####

# ### JSU specification

Heat_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

Heat_Model_Sample_JSU1 <- ugarchfit(Heat_Spec_Sample_JSU,Heat_sub1,out.sample = 250)
Heat_Model_Sample_JSU2 <- ugarchfit(Heat_Spec_Sample_JSU,Heat_sub2,out.sample = 250)
Heat_Model_Sample_JSU3 <- ugarchfit(Heat_Spec_Sample_JSU,Heat_sub3,out.sample = 250)
Heat_Model_Sample_JSU4 <- ugarchfit(Heat_Spec_Sample_JSU,Heat_sub4,out.sample = 250)
Heat_Model_Sample_JSU5 <- ugarchfit(Heat_Spec_Sample_JSU,Heat_sub5,out.sample = 250)

Gas_Spec_Sample_JSU <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="jsu")

Gas_Model_Sample_JSU1 <- ugarchfit(Gas_Spec_Sample_JSU,Gas_sub1,out.sample = 250)
Gas_Model_Sample_JSU2 <- ugarchfit(Gas_Spec_Sample_JSU,Gas_sub2,out.sample = 250)
Gas_Model_Sample_JSU3 <- ugarchfit(Gas_Spec_Sample_JSU,Gas_sub3,out.sample = 250)
Gas_Model_Sample_JSU4 <- ugarchfit(Gas_Spec_Sample_JSU,Gas_sub4,out.sample = 250)
Gas_Model_Sample_JSU5 <- ugarchfit(Gas_Spec_Sample_JSU,Gas_sub5,out.sample = 250)


## Rollover Analysis

Heat_Roll_JSU_Sample1 <-  ugarchroll(Heat_Spec_Sample_JSU,Heat_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_JSU_Sample2 <-  ugarchroll(Heat_Spec_Sample_JSU,Heat_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_JSU_Sample3 <-  ugarchroll(Heat_Spec_Sample_JSU,Heat_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_JSU_Sample4 <-  ugarchroll(Heat_Spec_Sample_JSU,Heat_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_JSU_Sample5 <-  ugarchroll(Heat_Spec_Sample_JSU,Heat_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Gas_Roll_JSU_Sample1 <-  ugarchroll(Gas_Spec_Sample_JSU,Gas_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_JSU_Sample2 <-  ugarchroll(Gas_Spec_Sample_JSU,Gas_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_JSU_Sample3 <-  ugarchroll(Gas_Spec_Sample_JSU,Gas_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_JSU_Sample4 <-  ugarchroll(Gas_Spec_Sample_JSU,Gas_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_JSU_Sample5 <-  ugarchroll(Gas_Spec_Sample_JSU,Gas_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### SGED specification

Heat_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

Heat_Model_Sample_SGED1 <- ugarchfit(Heat_Spec_Sample_SGED,Heat_sub1,out.sample = 250)
Heat_Model_Sample_SGED2 <- ugarchfit(Heat_Spec_Sample_SGED,Heat_sub2,out.sample = 250)
Heat_Model_Sample_SGED3 <- ugarchfit(Heat_Spec_Sample_SGED,Heat_sub3,out.sample = 250)
Heat_Model_Sample_SGED4 <- ugarchfit(Heat_Spec_Sample_SGED,Heat_sub4,out.sample = 250)
Heat_Model_Sample_SGED5 <- ugarchfit(Heat_Spec_Sample_SGED,Heat_sub5,out.sample = 250)

Gas_Spec_Sample_SGED <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sged")

Gas_Model_Sample_SGED1 <- ugarchfit(Gas_Spec_Sample_SGED,Gas_sub1,out.sample = 250)
Gas_Model_Sample_SGED2 <- ugarchfit(Gas_Spec_Sample_SGED,Gas_sub2,out.sample = 250)
Gas_Model_Sample_SGED3 <- ugarchfit(Gas_Spec_Sample_SGED,Gas_sub3,out.sample = 250)
Gas_Model_Sample_SGED4 <- ugarchfit(Gas_Spec_Sample_SGED,Gas_sub4,out.sample = 250)
Gas_Model_Sample_SGED5 <- ugarchfit(Gas_Spec_Sample_SGED,Gas_sub5,out.sample = 250)


## Rollover Analysis

Heat_Roll_SGED_Sample1 <-  ugarchroll(Heat_Spec_Sample_SGED,Heat_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_SGED_Sample2 <-  ugarchroll(Heat_Spec_Sample_SGED,Heat_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_SGED_Sample3 <-  ugarchroll(Heat_Spec_Sample_SGED,Heat_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_SGED_Sample4 <-  ugarchroll(Heat_Spec_Sample_SGED,Heat_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_SGED_Sample5 <-  ugarchroll(Heat_Spec_Sample_SGED,Heat_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Gas_Roll_SGED_Sample1 <-  ugarchroll(Gas_Spec_Sample_SGED,Gas_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_SGED_Sample2 <-  ugarchroll(Gas_Spec_Sample_SGED,Gas_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_SGED_Sample3 <-  ugarchroll(Gas_Spec_Sample_SGED,Gas_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_SGED_Sample4 <-  ugarchroll(Gas_Spec_Sample_SGED,Gas_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_SGED_Sample5 <-  ugarchroll(Gas_Spec_Sample_SGED,Gas_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### SST Specification
# 
Heat_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")
# 
Heat_Model_Sample_SST1 <- ugarchfit(Heat_Spec_Sample_SST,Heat_sub1,out.sample = 250)
Heat_Model_Sample_SST2 <- ugarchfit(Heat_Spec_Sample_SST,Heat_sub2,out.sample = 250)
Heat_Model_Sample_SST3 <- ugarchfit(Heat_Spec_Sample_SST,Heat_sub3,out.sample = 250)
Heat_Model_Sample_SST4 <- ugarchfit(Heat_Spec_Sample_SST,Heat_sub4,out.sample = 250)
Heat_Model_Sample_SST5 <- ugarchfit(Heat_Spec_Sample_SST,Heat_sub5,out.sample = 250)

Gas_Spec_Sample_SST <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="sstd")

Gas_Model_Sample_SST1 <- ugarchfit(Gas_Spec_Sample_SST,Gas_sub1,out.sample = 250)
Gas_Model_Sample_SST2 <- ugarchfit(Gas_Spec_Sample_SST,Gas_sub2,out.sample = 250)
Gas_Model_Sample_SST3 <- ugarchfit(Gas_Spec_Sample_SST,Gas_sub3,out.sample = 250)
Gas_Model_Sample_SST4 <- ugarchfit(Gas_Spec_Sample_SST,Gas_sub4,out.sample = 250)
Gas_Model_Sample_SST5 <- ugarchfit(Gas_Spec_Sample_SST,Gas_sub5,out.sample = 250)


## Rollover Analysis

Heat_Roll_SST_Sample1 <-  ugarchroll(Heat_Spec_Sample_SST,Heat_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_SST_Sample2 <-  ugarchroll(Heat_Spec_Sample_SST,Heat_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_SST_Sample3 <-  ugarchroll(Heat_Spec_Sample_SST,Heat_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_SST_Sample4 <-  ugarchroll(Heat_Spec_Sample_SST,Heat_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_SST_Sample5 <-  ugarchroll(Heat_Spec_Sample_SST,Heat_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Gas_Roll_SST_Sample1 <-  ugarchroll(Gas_Spec_Sample_SST,Gas_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_SST_Sample2 <-  ugarchroll(Gas_Spec_Sample_SST,Gas_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_SST_Sample3 <-  ugarchroll(Gas_Spec_Sample_SST,Gas_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_SST_Sample4 <-  ugarchroll(Gas_Spec_Sample_SST,Gas_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_SST_Sample5 <-  ugarchroll(Gas_Spec_Sample_SST,Gas_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)


# ### Normal Specification

Heat_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Heat_Model_Sample_Normal1 <- ugarchfit(Heat_Spec_Sample_Normal,Heat_sub1,out.sample = 250)
Heat_Model_Sample_Normal2 <- ugarchfit(Heat_Spec_Sample_Normal,Heat_sub2,out.sample = 250)
Heat_Model_Sample_Normal3 <- ugarchfit(Heat_Spec_Sample_Normal,Heat_sub3,out.sample = 250)
Heat_Model_Sample_Normal4 <- ugarchfit(Heat_Spec_Sample_Normal,Heat_sub4,out.sample = 250)
Heat_Model_Sample_Normal5 <- ugarchfit(Heat_Spec_Sample_Normal,Heat_sub5,out.sample = 250)

Gas_Spec_Sample_Normal <- ugarchspec(variance.model=list(model="eGARCH",garchOrder=c(1,1)),mean.model = list(armaOrder=c(1,0)),distribution.model="norm")

Gas_Model_Sample_Normal1 <- ugarchfit(Gas_Spec_Sample_Normal,Gas_sub1,out.sample = 250)
Gas_Model_Sample_Normal2 <- ugarchfit(Gas_Spec_Sample_Normal,Gas_sub2,out.sample = 250)
Gas_Model_Sample_Normal3 <- ugarchfit(Gas_Spec_Sample_Normal,Gas_sub3,out.sample = 250)
Gas_Model_Sample_Normal4 <- ugarchfit(Gas_Spec_Sample_Normal,Gas_sub4,out.sample = 250)
Gas_Model_Sample_Normal5 <- ugarchfit(Gas_Spec_Sample_Normal,Gas_sub5,out.sample = 250)

Heat_Roll_Normal_Sample1 <-  ugarchroll(Heat_Spec_Sample_Normal,Heat_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_Normal_Sample2 <-  ugarchroll(Heat_Spec_Sample_Normal,Heat_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_Normal_Sample3 <-  ugarchroll(Heat_Spec_Sample_Normal,Heat_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_Normal_Sample4 <-  ugarchroll(Heat_Spec_Sample_Normal,Heat_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Heat_Roll_Normal_Sample5 <-  ugarchroll(Heat_Spec_Sample_Normal,Heat_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

Gas_Roll_Normal_Sample1 <-  ugarchroll(Gas_Spec_Sample_Normal,Gas_sub1,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_Normal_Sample2 <-  ugarchroll(Gas_Spec_Sample_Normal,Gas_sub2,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_Normal_Sample3 <-  ugarchroll(Gas_Spec_Sample_Normal,Gas_sub3,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_Normal_Sample4 <-  ugarchroll(Gas_Spec_Sample_Normal,Gas_sub4,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)
Gas_Roll_Normal_Sample5 <-  ugarchroll(Gas_Spec_Sample_Normal,Gas_sub5,forecast.length = 250,refit.every = 25,refit.window = "moving",solver = "hybrid",calculate.VaR = TRUE,VaR.alpha = c(0.01,0.05,0.10),keep.coef = TRUE)

### Sample Residuals Statistics

Heat_Res_Sample1 <- c(normalTest(Heat_Model_Sample_Normal1@fit$residuals/Heat_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(Heat_Model_Sample_Normal1@fit$residuals/Heat_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(Heat_Model_Sample_Normal1@fit$residuals/Heat_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
Heat_Res_Sample2 <- c(normalTest(Heat_Model_Sample_Normal2@fit$residuals/Heat_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(Heat_Model_Sample_Normal2@fit$residuals/Heat_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(Heat_Model_Sample_Normal2@fit$residuals/Heat_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
Heat_Res_Sample3 <- c(normalTest(Heat_Model_Sample_Normal3@fit$residuals/Heat_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(Heat_Model_Sample_Normal3@fit$residuals/Heat_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(Heat_Model_Sample_Normal3@fit$residuals/Heat_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
Heat_Res_Sample4 <- c(normalTest(Heat_Model_Sample_Normal4@fit$residuals/Heat_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(Heat_Model_Sample_Normal4@fit$residuals/Heat_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(Heat_Model_Sample_Normal4@fit$residuals/Heat_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
Heat_Res_Sample5 <- c(normalTest(Heat_Model_Sample_Normal5@fit$residuals/Heat_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(Heat_Model_Sample_Normal5@fit$residuals/Heat_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(Heat_Model_Sample_Normal5@fit$residuals/Heat_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)

Gas_Res_Sample1 <- c(normalTest(Gas_Model_Sample_Normal1@fit$residuals/Gas_Model_Sample_Normal1@fit$sigma,method="jb")@test$statistic,normalTest(Gas_Model_Sample_Normal1@fit$residuals/Gas_Model_Sample_Normal1@fit$sigma,method="sw")@test$statistic,normalTest(Gas_Model_Sample_Normal1@fit$residuals/Gas_Model_Sample_Normal1@fit$sigma,method="ks")@test$statistic)
Gas_Res_Sample2 <- c(normalTest(Gas_Model_Sample_Normal2@fit$residuals/Gas_Model_Sample_Normal2@fit$sigma,method="jb")@test$statistic,normalTest(Gas_Model_Sample_Normal2@fit$residuals/Gas_Model_Sample_Normal2@fit$sigma,method="sw")@test$statistic,normalTest(Gas_Model_Sample_Normal2@fit$residuals/Gas_Model_Sample_Normal2@fit$sigma,method="ks")@test$statistic)
Gas_Res_Sample3 <- c(normalTest(Gas_Model_Sample_Normal3@fit$residuals/Gas_Model_Sample_Normal3@fit$sigma,method="jb")@test$statistic,normalTest(Gas_Model_Sample_Normal3@fit$residuals/Gas_Model_Sample_Normal3@fit$sigma,method="sw")@test$statistic,normalTest(Gas_Model_Sample_Normal3@fit$residuals/Gas_Model_Sample_Normal3@fit$sigma,method="ks")@test$statistic)
Gas_Res_Sample4 <- c(normalTest(Gas_Model_Sample_Normal4@fit$residuals/Gas_Model_Sample_Normal4@fit$sigma,method="jb")@test$statistic,normalTest(Gas_Model_Sample_Normal4@fit$residuals/Gas_Model_Sample_Normal4@fit$sigma,method="sw")@test$statistic,normalTest(Gas_Model_Sample_Normal4@fit$residuals/Gas_Model_Sample_Normal4@fit$sigma,method="ks")@test$statistic)
Gas_Res_Sample5 <- c(normalTest(Gas_Model_Sample_Normal5@fit$residuals/Gas_Model_Sample_Normal5@fit$sigma,method="jb")@test$statistic,normalTest(Gas_Model_Sample_Normal5@fit$residuals/Gas_Model_Sample_Normal5@fit$sigma,method="sw")@test$statistic,normalTest(Gas_Model_Sample_Normal5@fit$residuals/Gas_Model_Sample_Normal5@fit$sigma,method="ks")@test$statistic)

Heat_Res_Sample_Test <- rbind(Heat_Res_Sample1,Heat_Res_Sample2,Heat_Res_Sample3,Heat_Res_Sample4,Heat_Res_Sample5)
Gas_Res_Sample_Test   <- rbind(Gas_Res_Sample1,Gas_Res_Sample2,Gas_Res_Sample3,Gas_Res_Sample4,Gas_Res_Sample5)
rownames(Heat_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")
colnames(Heat_Res_Sample_Test) <- c("JB","SW","KS") 
rownames(Gas_Res_Sample_Test) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")
colnames(Gas_Res_Sample_Test) <- c("JB","SW","KS") 

### Pearson Specification

Heat_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub1,out.sample = 250)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub1,out.sample = 250)@fit$sigma)
Heat_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub2)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub2,out.sample = 250)@fit$sigma)
Heat_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub3,out.sample = 250)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub3,out.sample = 250)@fit$sigma)
Heat_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub4,out.sample = 250)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub4,out.sample = 250)@fit$sigma)
Heat_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub5,out.sample = 250)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub5,out.sample = 250)@fit$sigma)

Heat_Param_Pearson_Subsample <- rbind(Heat_Param_Pearson_Sample1,Heat_Param_Pearson_Sample2,Heat_Param_Pearson_Sample3,Heat_Param_Pearson_Sample4,Heat_Param_Pearson_Sample5)
rownames(Heat_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")

Gas_Param_Pearson_Sample1 <- pearsonFitML(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub1,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub1,out.sample = 250)@fit$sigma)
Gas_Param_Pearson_Sample2 <- pearsonFitML(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub2)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub2,out.sample = 250)@fit$sigma)
Gas_Param_Pearson_Sample3 <- pearsonFitML(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub3,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub3,out.sample = 250)@fit$sigma)
Gas_Param_Pearson_Sample4 <- pearsonFitML(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub4,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub4,out.sample = 250)@fit$sigma)
Gas_Param_Pearson_Sample5 <- pearsonFitML(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub5,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub5,out.sample = 250)@fit$sigma)

Gas_Param_Pearson_Subsample <- rbind(Gas_Param_Pearson_Sample1,Gas_Param_Pearson_Sample2,Gas_Param_Pearson_Sample3,Gas_Param_Pearson_Sample4,Gas_Param_Pearson_Sample5)
rownames(Gas_Param_Pearson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")

### johnson SU parameter


Heat_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub1,out.sample = 250)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub1,out.sample = 250)@fit$sigma,moment = "find")
Heat_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub2)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub2)@fit$sigma,moment = "find")
Heat_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub3,out.sample = 250)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub3,out.sample = 250)@fit$sigma,moment = "find")
Heat_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub4,out.sample = 250)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub4,out.sample = 250)@fit$sigma,moment = "find")
Heat_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(Heat_Spec_Sample_Normal,Heat_sub5,out.sample = 250)@fit$residuals/ugarchfit(Heat_Spec_Sample_Normal,Heat_sub5,out.sample = 250)@fit$sigma,moment = "find")

Heat_Param_Johnson_Subsample <- rbind(Heat_Param_Johnson_Sample1,Heat_Param_Johnson_Sample2,Heat_Param_Johnson_Sample3,Heat_Param_Johnson_Sample4,Heat_Param_Johnson_Sample5)
rownames(Heat_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")

Gas_Param_Johnson_Sample1 <- JohnsonFit(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub1,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub1,out.sample = 250)@fit$sigma,moment = "find")
Gas_Param_Johnson_Sample2 <- JohnsonFit(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub2,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub2,out.sample = 250)@fit$sigma,moment = "find")
Gas_Param_Johnson_Sample3 <- JohnsonFit(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub3,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub3,out.sample = 250)@fit$sigma,moment = "find")
Gas_Param_Johnson_Sample4 <- JohnsonFit(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub4,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub4,out.sample = 250)@fit$sigma,moment = "find")
Gas_Param_Johnson_Sample5 <- JohnsonFit(ugarchfit(Gas_Spec_Sample_Normal,Gas_sub5,out.sample = 250)@fit$residuals/ugarchfit(Gas_Spec_Sample_Normal,Gas_sub5,out.sample = 250)@fit$sigma,moment = "find")

Gas_Param_Johnson_Subsample <- rbind(Gas_Param_Johnson_Sample1,Gas_Param_Johnson_Sample2,Gas_Param_Johnson_Sample3,Gas_Param_Johnson_Sample4,Gas_Param_Johnson_Sample5)
rownames(Gas_Param_Johnson_Subsample) <- c("Sample1","Sample2","Sample3","Sample4","Sample5")


#### Backtesting and Storing the results

Heat_Test_SST1_Sample1 <- BacktestVaR(Heat_Roll_SST_Sample1@forecast$VaR[,4],Heat_Roll_SST_Sample1@forecast$VaR[,1],0.01)
Heat_Test_SST5_Sample1 <- BacktestVaR(Heat_Roll_SST_Sample1@forecast$VaR[,4],Heat_Roll_SST_Sample1@forecast$VaR[,2],0.05)
Heat_Test_SST10_Sample1 <- BacktestVaR(Heat_Roll_SST_Sample1@forecast$VaR[,4],Heat_Roll_SST_Sample1@forecast$VaR[,3],0.1)

Gas_Test_SST1_Sample1 <- BacktestVaR(Gas_Roll_SST_Sample1@forecast$VaR[,4],Gas_Roll_SST_Sample1@forecast$VaR[,1],0.01)
Gas_Test_SST5_Sample1 <- BacktestVaR(Gas_Roll_SST_Sample1@forecast$VaR[,4],Gas_Roll_SST_Sample1@forecast$VaR[,2],0.05)
Gas_Test_SST10_Sample1 <- BacktestVaR(Gas_Roll_SST_Sample1@forecast$VaR[,4],Gas_Roll_SST_Sample1@forecast$VaR[,3],0.1)

Heat_Test_SGED1_Sample1 <- BacktestVaR(Heat_Roll_SGED_Sample1@forecast$VaR[,4],Heat_Roll_SGED_Sample1@forecast$VaR[,1],0.01)
Heat_Test_SGED5_Sample1 <- BacktestVaR(Heat_Roll_SGED_Sample1@forecast$VaR[,4],Heat_Roll_SGED_Sample1@forecast$VaR[,2],0.05)
Heat_Test_SGED10_Sample1 <- BacktestVaR(Heat_Roll_SGED_Sample1@forecast$VaR[,4],Heat_Roll_SGED_Sample1@forecast$VaR[,3],0.1)

Gas_Test_SGED1_Sample1 <- BacktestVaR(Gas_Roll_SGED_Sample1@forecast$VaR[,4],Gas_Roll_SGED_Sample1@forecast$VaR[,1],0.01)
Gas_Test_SGED5_Sample1 <- BacktestVaR(Gas_Roll_SGED_Sample1@forecast$VaR[,4],Gas_Roll_SGED_Sample1@forecast$VaR[,2],0.05)
Gas_Test_SGED10_Sample1 <- BacktestVaR(Gas_Roll_SGED_Sample1@forecast$VaR[,4],Gas_Roll_SGED_Sample1@forecast$VaR[,3],0.1)

Heat_Test_JSU1_Sample1 <- BacktestVaR(Heat_Roll_JSU_Sample1@forecast$VaR[,4],Heat_Roll_JSU_Sample1@forecast$VaR[,1],0.01)
Heat_Test_JSU5_Sample1 <- BacktestVaR(Heat_Roll_JSU_Sample1@forecast$VaR[,4],Heat_Roll_JSU_Sample1@forecast$VaR[,2],0.05)
Heat_Test_JSU10_Sample1 <- BacktestVaR(Heat_Roll_JSU_Sample1@forecast$VaR[,4],Heat_Roll_JSU_Sample1@forecast$VaR[,3],0.1)

Gas_Test_JSU1_Sample1 <- BacktestVaR(Gas_Roll_JSU_Sample1@forecast$VaR[,4],Gas_Roll_JSU_Sample1@forecast$VaR[,1],0.01)
Gas_Test_JSU5_Sample1 <- BacktestVaR(Gas_Roll_JSU_Sample1@forecast$VaR[,4],Gas_Roll_JSU_Sample1@forecast$VaR[,2],0.05)
Gas_Test_JSU10_Sample1 <- BacktestVaR(Gas_Roll_JSU_Sample1@forecast$VaR[,4],Gas_Roll_JSU_Sample1@forecast$VaR[,3],0.1)

Heat_Test_Normal1_Sample1 <- BacktestVaR(Heat_Roll_Normal_Sample1@forecast$VaR[,4],Heat_Roll_Normal_Sample1@forecast$VaR[,1],0.01)
Heat_Test_Normal5_Sample1 <- BacktestVaR(Heat_Roll_Normal_Sample1@forecast$VaR[,4],Heat_Roll_Normal_Sample1@forecast$VaR[,2],0.05)
Heat_Test_Normal10_Sample1 <- BacktestVaR(Heat_Roll_Normal_Sample1@forecast$VaR[,4],Heat_Roll_Normal_Sample1@forecast$VaR[,3],0.1)

Gas_Test_Normal1_Sample1 <- BacktestVaR(Gas_Roll_Normal_Sample1@forecast$VaR[,4],Gas_Roll_Normal_Sample1@forecast$VaR[,1],0.01)
Gas_Test_Normal5_Sample1 <- BacktestVaR(Gas_Roll_Normal_Sample1@forecast$VaR[,4],Gas_Roll_Normal_Sample1@forecast$VaR[,2],0.05)
Gas_Test_Normal10_Sample1 <- BacktestVaR(Gas_Roll_Normal_Sample1@forecast$VaR[,4],Gas_Roll_Normal_Sample1@forecast$VaR[,3],0.1)

##### Sample 2


Heat_Test_SST1_Sample2 <- BacktestVaR(Heat_Roll_SST_Sample2@forecast$VaR[,4],Heat_Roll_SST_Sample2@forecast$VaR[,1],0.01)
Heat_Test_SST5_Sample2 <- BacktestVaR(Heat_Roll_SST_Sample2@forecast$VaR[,4],Heat_Roll_SST_Sample2@forecast$VaR[,2],0.05)
Heat_Test_SST10_Sample2 <- BacktestVaR(Heat_Roll_SST_Sample2@forecast$VaR[,4],Heat_Roll_SST_Sample2@forecast$VaR[,3],0.1)

Gas_Test_SST1_Sample2 <- BacktestVaR(Gas_Roll_SST_Sample2@forecast$VaR[,4],Gas_Roll_SST_Sample2@forecast$VaR[,1],0.01)
Gas_Test_SST5_Sample2 <- BacktestVaR(Gas_Roll_SST_Sample2@forecast$VaR[,4],Gas_Roll_SST_Sample2@forecast$VaR[,2],0.05)
Gas_Test_SST10_Sample2 <- BacktestVaR(Gas_Roll_SST_Sample2@forecast$VaR[,4],Gas_Roll_SST_Sample2@forecast$VaR[,3],0.1)

Heat_Test_SGED1_Sample2 <- BacktestVaR(Heat_Roll_SGED_Sample2@forecast$VaR[,4],Heat_Roll_SGED_Sample2@forecast$VaR[,1],0.01)
Heat_Test_SGED5_Sample2 <- BacktestVaR(Heat_Roll_SGED_Sample2@forecast$VaR[,4],Heat_Roll_SGED_Sample2@forecast$VaR[,2],0.05)
Heat_Test_SGED10_Sample2 <- BacktestVaR(Heat_Roll_SGED_Sample2@forecast$VaR[,4],Heat_Roll_SGED_Sample2@forecast$VaR[,3],0.1)

Gas_Test_SGED1_Sample2 <- BacktestVaR(Gas_Roll_SGED_Sample2@forecast$VaR[,4],Gas_Roll_SGED_Sample2@forecast$VaR[,1],0.01)
Gas_Test_SGED5_Sample2 <- BacktestVaR(Gas_Roll_SGED_Sample2@forecast$VaR[,4],Gas_Roll_SGED_Sample2@forecast$VaR[,2],0.05)
Gas_Test_SGED10_Sample2 <- BacktestVaR(Gas_Roll_SGED_Sample2@forecast$VaR[,4],Gas_Roll_SGED_Sample2@forecast$VaR[,3],0.1)

Heat_Test_JSU1_Sample2 <- BacktestVaR(Heat_Roll_JSU_Sample2@forecast$VaR[,4],Heat_Roll_JSU_Sample2@forecast$VaR[,1],0.01)
Heat_Test_JSU5_Sample2 <- BacktestVaR(Heat_Roll_JSU_Sample2@forecast$VaR[,4],Heat_Roll_JSU_Sample2@forecast$VaR[,2],0.05)
Heat_Test_JSU10_Sample2 <- BacktestVaR(Heat_Roll_JSU_Sample2@forecast$VaR[,4],Heat_Roll_JSU_Sample2@forecast$VaR[,3],0.1)

Gas_Test_JSU1_Sample2 <- BacktestVaR(Gas_Roll_JSU_Sample2@forecast$VaR[,4],Gas_Roll_JSU_Sample2@forecast$VaR[,1],0.01)
Gas_Test_JSU5_Sample2 <- BacktestVaR(Gas_Roll_JSU_Sample2@forecast$VaR[,4],Gas_Roll_JSU_Sample2@forecast$VaR[,2],0.05)
Gas_Test_JSU10_Sample2 <- BacktestVaR(Gas_Roll_JSU_Sample2@forecast$VaR[,4],Gas_Roll_JSU_Sample2@forecast$VaR[,3],0.1)

Heat_Test_Normal1_Sample2 <- BacktestVaR(Heat_Roll_Normal_Sample2@forecast$VaR[,4],Heat_Roll_Normal_Sample2@forecast$VaR[,1],0.01)
Heat_Test_Normal5_Sample2 <- BacktestVaR(Heat_Roll_Normal_Sample2@forecast$VaR[,4],Heat_Roll_Normal_Sample2@forecast$VaR[,2],0.05)
Heat_Test_Normal10_Sample2 <- BacktestVaR(Heat_Roll_Normal_Sample2@forecast$VaR[,4],Heat_Roll_Normal_Sample2@forecast$VaR[,3],0.1)

Gas_Test_Normal1_Sample2 <- BacktestVaR(Gas_Roll_Normal_Sample2@forecast$VaR[,4],Gas_Roll_Normal_Sample2@forecast$VaR[,1],0.01)
Gas_Test_Normal5_Sample2 <- BacktestVaR(Gas_Roll_Normal_Sample2@forecast$VaR[,4],Gas_Roll_Normal_Sample2@forecast$VaR[,2],0.05)
Gas_Test_Normal10_Sample2 <- BacktestVaR(Gas_Roll_Normal_Sample2@forecast$VaR[,4],Gas_Roll_Normal_Sample2@forecast$VaR[,3],0.1)

### Sample3


Heat_Test_SST1_Sample3 <- BacktestVaR(Heat_Roll_SST_Sample3@forecast$VaR[,4],Heat_Roll_SST_Sample3@forecast$VaR[,1],0.01)
Heat_Test_SST5_Sample3 <- BacktestVaR(Heat_Roll_SST_Sample3@forecast$VaR[,4],Heat_Roll_SST_Sample3@forecast$VaR[,2],0.05)
Heat_Test_SST10_Sample3 <- BacktestVaR(Heat_Roll_SST_Sample3@forecast$VaR[,4],Heat_Roll_SST_Sample3@forecast$VaR[,3],0.1)

Gas_Test_SST1_Sample3 <- BacktestVaR(Gas_Roll_SST_Sample3@forecast$VaR[,4],Gas_Roll_SST_Sample3@forecast$VaR[,1],0.01)
Gas_Test_SST5_Sample3 <- BacktestVaR(Gas_Roll_SST_Sample3@forecast$VaR[,4],Gas_Roll_SST_Sample3@forecast$VaR[,2],0.05)
Gas_Test_SST10_Sample3 <- BacktestVaR(Gas_Roll_SST_Sample3@forecast$VaR[,4],Gas_Roll_SST_Sample3@forecast$VaR[,3],0.1)

Heat_Test_SGED1_Sample3 <- BacktestVaR(Heat_Roll_SGED_Sample3@forecast$VaR[,4],Heat_Roll_SGED_Sample3@forecast$VaR[,1],0.01)
Heat_Test_SGED5_Sample3 <- BacktestVaR(Heat_Roll_SGED_Sample3@forecast$VaR[,4],Heat_Roll_SGED_Sample3@forecast$VaR[,2],0.05)
Heat_Test_SGED10_Sample3 <- BacktestVaR(Heat_Roll_SGED_Sample3@forecast$VaR[,4],Heat_Roll_SGED_Sample3@forecast$VaR[,3],0.1)

Gas_Test_SGED1_Sample3 <- BacktestVaR(Gas_Roll_SGED_Sample3@forecast$VaR[,4],Gas_Roll_SGED_Sample3@forecast$VaR[,1],0.01)
Gas_Test_SGED5_Sample3 <- BacktestVaR(Gas_Roll_SGED_Sample3@forecast$VaR[,4],Gas_Roll_SGED_Sample3@forecast$VaR[,2],0.05)
Gas_Test_SGED10_Sample3 <- BacktestVaR(Gas_Roll_SGED_Sample3@forecast$VaR[,4],Gas_Roll_SGED_Sample3@forecast$VaR[,3],0.1)

Heat_Test_JSU1_Sample3 <- BacktestVaR(Heat_Roll_JSU_Sample3@forecast$VaR[,4],Heat_Roll_JSU_Sample3@forecast$VaR[,1],0.01)
Heat_Test_JSU5_Sample3 <- BacktestVaR(Heat_Roll_JSU_Sample3@forecast$VaR[,4],Heat_Roll_JSU_Sample3@forecast$VaR[,2],0.05)
Heat_Test_JSU10_Sample3 <- BacktestVaR(Heat_Roll_JSU_Sample3@forecast$VaR[,4],Heat_Roll_JSU_Sample3@forecast$VaR[,3],0.1)

Gas_Test_JSU1_Sample3 <- BacktestVaR(Gas_Roll_JSU_Sample3@forecast$VaR[,4],Gas_Roll_JSU_Sample3@forecast$VaR[,1],0.01)
Gas_Test_JSU5_Sample3 <- BacktestVaR(Gas_Roll_JSU_Sample3@forecast$VaR[,4],Gas_Roll_JSU_Sample3@forecast$VaR[,2],0.05)
Gas_Test_JSU10_Sample3 <- BacktestVaR(Gas_Roll_JSU_Sample3@forecast$VaR[,4],Gas_Roll_JSU_Sample3@forecast$VaR[,3],0.1)

Heat_Test_Normal1_Sample3 <- BacktestVaR(Heat_Roll_Normal_Sample3@forecast$VaR[,4],Heat_Roll_Normal_Sample3@forecast$VaR[,1],0.01)
Heat_Test_Normal5_Sample3 <- BacktestVaR(Heat_Roll_Normal_Sample3@forecast$VaR[,4],Heat_Roll_Normal_Sample3@forecast$VaR[,2],0.05)
Heat_Test_Normal10_Sample3 <- BacktestVaR(Heat_Roll_Normal_Sample3@forecast$VaR[,4],Heat_Roll_Normal_Sample3@forecast$VaR[,3],0.1)

Gas_Test_Normal1_Sample3 <- BacktestVaR(Gas_Roll_Normal_Sample3@forecast$VaR[,4],Gas_Roll_Normal_Sample3@forecast$VaR[,1],0.01)
Gas_Test_Normal5_Sample3 <- BacktestVaR(Gas_Roll_Normal_Sample3@forecast$VaR[,4],Gas_Roll_Normal_Sample3@forecast$VaR[,2],0.05)
Gas_Test_Normal10_Sample3 <- BacktestVaR(Gas_Roll_Normal_Sample3@forecast$VaR[,4],Gas_Roll_Normal_Sample3@forecast$VaR[,3],0.1)

### Sample 4


Heat_Test_SST1_Sample4 <- BacktestVaR(Heat_Roll_SST_Sample4@forecast$VaR[,4],Heat_Roll_SST_Sample4@forecast$VaR[,1],0.01)
Heat_Test_SST5_Sample4 <- BacktestVaR(Heat_Roll_SST_Sample4@forecast$VaR[,4],Heat_Roll_SST_Sample4@forecast$VaR[,2],0.05)
Heat_Test_SST10_Sample4 <- BacktestVaR(Heat_Roll_SST_Sample4@forecast$VaR[,4],Heat_Roll_SST_Sample4@forecast$VaR[,3],0.1)

Gas_Test_SST1_Sample4 <- BacktestVaR(Gas_Roll_SST_Sample4@forecast$VaR[,4],Gas_Roll_SST_Sample4@forecast$VaR[,1],0.01)
Gas_Test_SST5_Sample4 <- BacktestVaR(Gas_Roll_SST_Sample4@forecast$VaR[,4],Gas_Roll_SST_Sample4@forecast$VaR[,2],0.05)
Gas_Test_SST10_Sample4 <- BacktestVaR(Gas_Roll_SST_Sample4@forecast$VaR[,4],Gas_Roll_SST_Sample4@forecast$VaR[,3],0.1)

Heat_Test_SGED1_Sample4 <- BacktestVaR(Heat_Roll_SGED_Sample4@forecast$VaR[,4],Heat_Roll_SGED_Sample4@forecast$VaR[,1],0.01)
Heat_Test_SGED5_Sample4 <- BacktestVaR(Heat_Roll_SGED_Sample4@forecast$VaR[,4],Heat_Roll_SGED_Sample4@forecast$VaR[,2],0.05)
Heat_Test_SGED10_Sample4 <- BacktestVaR(Heat_Roll_SGED_Sample4@forecast$VaR[,4],Heat_Roll_SGED_Sample4@forecast$VaR[,3],0.1)

Gas_Test_SGED1_Sample4 <- BacktestVaR(Gas_Roll_SGED_Sample4@forecast$VaR[,4],Gas_Roll_SGED_Sample4@forecast$VaR[,1],0.01)
Gas_Test_SGED5_Sample4 <- BacktestVaR(Gas_Roll_SGED_Sample4@forecast$VaR[,4],Gas_Roll_SGED_Sample4@forecast$VaR[,2],0.05)
Gas_Test_SGED10_Sample4 <- BacktestVaR(Gas_Roll_SGED_Sample4@forecast$VaR[,4],Gas_Roll_SGED_Sample4@forecast$VaR[,3],0.1)

Heat_Test_JSU1_Sample4 <- BacktestVaR(Heat_Roll_JSU_Sample4@forecast$VaR[,4],Heat_Roll_JSU_Sample4@forecast$VaR[,1],0.01)
Heat_Test_JSU5_Sample4 <- BacktestVaR(Heat_Roll_JSU_Sample4@forecast$VaR[,4],Heat_Roll_JSU_Sample4@forecast$VaR[,2],0.05)
Heat_Test_JSU10_Sample4 <- BacktestVaR(Heat_Roll_JSU_Sample4@forecast$VaR[,4],Heat_Roll_JSU_Sample4@forecast$VaR[,3],0.1)

Gas_Test_JSU1_Sample4 <- BacktestVaR(Gas_Roll_JSU_Sample4@forecast$VaR[,4],Gas_Roll_JSU_Sample4@forecast$VaR[,1],0.01)
Gas_Test_JSU5_Sample4 <- BacktestVaR(Gas_Roll_JSU_Sample4@forecast$VaR[,4],Gas_Roll_JSU_Sample4@forecast$VaR[,2],0.05)
Gas_Test_JSU10_Sample4 <- BacktestVaR(Gas_Roll_JSU_Sample4@forecast$VaR[,4],Gas_Roll_JSU_Sample4@forecast$VaR[,3],0.1)

Heat_Test_Normal1_Sample4 <- BacktestVaR(Heat_Roll_Normal_Sample4@forecast$VaR[,4],Heat_Roll_Normal_Sample4@forecast$VaR[,1],0.01)
Heat_Test_Normal5_Sample4 <- BacktestVaR(Heat_Roll_Normal_Sample4@forecast$VaR[,4],Heat_Roll_Normal_Sample4@forecast$VaR[,2],0.05)
Heat_Test_Normal10_Sample4 <- BacktestVaR(Heat_Roll_Normal_Sample4@forecast$VaR[,4],Heat_Roll_Normal_Sample4@forecast$VaR[,3],0.1)

Gas_Test_Normal1_Sample4 <- BacktestVaR(Gas_Roll_Normal_Sample4@forecast$VaR[,4],Gas_Roll_Normal_Sample4@forecast$VaR[,1],0.01)
Gas_Test_Normal5_Sample4 <- BacktestVaR(Gas_Roll_Normal_Sample4@forecast$VaR[,4],Gas_Roll_Normal_Sample4@forecast$VaR[,2],0.05)
Gas_Test_Normal10_Sample4 <- BacktestVaR(Gas_Roll_Normal_Sample4@forecast$VaR[,4],Gas_Roll_Normal_Sample4@forecast$VaR[,3],0.1)

## Sample 5


Heat_Test_SST1_Sample5 <- BacktestVaR(Heat_Roll_SST_Sample5@forecast$VaR[,4],Heat_Roll_SST_Sample5@forecast$VaR[,1],0.01)
Heat_Test_SST5_Sample5 <- BacktestVaR(Heat_Roll_SST_Sample5@forecast$VaR[,4],Heat_Roll_SST_Sample5@forecast$VaR[,2],0.05)
Heat_Test_SST10_Sample5 <- BacktestVaR(Heat_Roll_SST_Sample5@forecast$VaR[,4],Heat_Roll_SST_Sample5@forecast$VaR[,3],0.1)

Gas_Test_SST1_Sample5 <- BacktestVaR(Gas_Roll_SST_Sample5@forecast$VaR[,4],Gas_Roll_SST_Sample5@forecast$VaR[,1],0.01)
Gas_Test_SST5_Sample5 <- BacktestVaR(Gas_Roll_SST_Sample5@forecast$VaR[,4],Gas_Roll_SST_Sample5@forecast$VaR[,2],0.05)
Gas_Test_SST10_Sample5 <- BacktestVaR(Gas_Roll_SST_Sample5@forecast$VaR[,4],Gas_Roll_SST_Sample5@forecast$VaR[,3],0.1)

Heat_Test_SGED1_Sample5 <- BacktestVaR(Heat_Roll_SGED_Sample5@forecast$VaR[,4],Heat_Roll_SGED_Sample5@forecast$VaR[,1],0.01)
Heat_Test_SGED5_Sample5 <- BacktestVaR(Heat_Roll_SGED_Sample5@forecast$VaR[,4],Heat_Roll_SGED_Sample5@forecast$VaR[,2],0.05)
Heat_Test_SGED10_Sample5 <- BacktestVaR(Heat_Roll_SGED_Sample5@forecast$VaR[,4],Heat_Roll_SGED_Sample5@forecast$VaR[,3],0.1)

Gas_Test_SGED1_Sample5 <- BacktestVaR(Gas_Roll_SGED_Sample5@forecast$VaR[,4],Gas_Roll_SGED_Sample5@forecast$VaR[,1],0.01)
Gas_Test_SGED5_Sample5 <- BacktestVaR(Gas_Roll_SGED_Sample5@forecast$VaR[,4],Gas_Roll_SGED_Sample5@forecast$VaR[,2],0.05)
Gas_Test_SGED10_Sample5 <- BacktestVaR(Gas_Roll_SGED_Sample5@forecast$VaR[,4],Gas_Roll_SGED_Sample5@forecast$VaR[,3],0.1)

Heat_Test_JSU1_Sample5 <- BacktestVaR(Heat_Roll_JSU_Sample5@forecast$VaR[,4],Heat_Roll_JSU_Sample5@forecast$VaR[,1],0.01)
Heat_Test_JSU5_Sample5 <- BacktestVaR(Heat_Roll_JSU_Sample5@forecast$VaR[,4],Heat_Roll_JSU_Sample5@forecast$VaR[,2],0.05)
Heat_Test_JSU10_Sample5 <- BacktestVaR(Heat_Roll_JSU_Sample5@forecast$VaR[,4],Heat_Roll_JSU_Sample5@forecast$VaR[,3],0.1)

Gas_Test_JSU1_Sample5 <- BacktestVaR(Gas_Roll_JSU_Sample5@forecast$VaR[,4],Gas_Roll_JSU_Sample5@forecast$VaR[,1],0.01)
Gas_Test_JSU5_Sample5 <- BacktestVaR(Gas_Roll_JSU_Sample5@forecast$VaR[,4],Gas_Roll_JSU_Sample5@forecast$VaR[,2],0.05)
Gas_Test_JSU10_Sample5 <- BacktestVaR(Gas_Roll_JSU_Sample5@forecast$VaR[,4],Gas_Roll_JSU_Sample5@forecast$VaR[,3],0.1)

Heat_Test_Normal1_Sample5 <- BacktestVaR(Heat_Roll_Normal_Sample5@forecast$VaR[,4],Heat_Roll_Normal_Sample5@forecast$VaR[,1],0.01)
Heat_Test_Normal5_Sample5 <- BacktestVaR(Heat_Roll_Normal_Sample5@forecast$VaR[,4],Heat_Roll_Normal_Sample5@forecast$VaR[,2],0.05)
Heat_Test_Normal10_Sample5 <- BacktestVaR(Heat_Roll_Normal_Sample5@forecast$VaR[,4],Heat_Roll_Normal_Sample5@forecast$VaR[,3],0.1)

Gas_Test_Normal1_Sample5 <- BacktestVaR(Gas_Roll_Normal_Sample5@forecast$VaR[,4],Gas_Roll_Normal_Sample5@forecast$VaR[,1],0.01)
Gas_Test_Normal5_Sample5 <- BacktestVaR(Gas_Roll_Normal_Sample5@forecast$VaR[,4],Gas_Roll_Normal_Sample5@forecast$VaR[,2],0.05)
Gas_Test_Normal10_Sample5 <- BacktestVaR(Gas_Roll_Normal_Sample5@forecast$VaR[,4],Gas_Roll_Normal_Sample5@forecast$VaR[,3],0.1)


Heat_Sample1 <- matrix(c(Heat_Test_Normal1_Sample1$LRuc[1],Heat_Test_SST1_Sample1$LRuc[1],Heat_Test_SGED1_Sample1$LRuc[1],Heat_Test_JSU1_Sample1$LRuc[1],
  Heat_Test_Normal1_Sample1$LRcc[1],Heat_Test_SST1_Sample1$LRcc[1],Heat_Test_SGED1_Sample1$LRcc[1],Heat_Test_JSU1_Sample1$LRcc[1],
  Heat_Test_Normal1_Sample1$DQ$stat,Heat_Test_SST1_Sample1$DQ$stat,Heat_Test_SGED1_Sample1$DQ$stat,Heat_Test_JSU1_Sample1$DQ$stat,
  Heat_Test_Normal1_Sample1$AE,Heat_Test_SST1_Sample1$AE,Heat_Test_SGED1_Sample1$AE,Heat_Test_JSU1_Sample1$AE,
  (Heat_Test_Normal1_Sample1$Loss$Loss/Heat_Test_Normal1_Sample1$Loss$Loss),(Heat_Test_SST1_Sample1$Loss$Loss/Heat_Test_Normal1_Sample1$Loss$Loss),
  (Heat_Test_SGED1_Sample1$Loss$Loss/Heat_Test_Normal1_Sample1$Loss$Loss),(Heat_Test_JSU1_Sample1$Loss$Loss/Heat_Test_Normal1_Sample1$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Heat_Sample1) <- c("Normal","SST","SGED","JSU")    
rownames(Heat_Sample1) <- c("UC","CC","DQ","AE","QL")    

Heat_Sample2 <- matrix(c(Heat_Test_Normal1_Sample2$LRuc[1],Heat_Test_SST1_Sample2$LRuc[1],Heat_Test_SGED1_Sample2$LRuc[1],Heat_Test_JSU1_Sample2$LRuc[1],
                          Heat_Test_Normal1_Sample2$LRcc[1],Heat_Test_SST1_Sample2$LRcc[1],Heat_Test_SGED1_Sample2$LRcc[1],Heat_Test_JSU1_Sample2$LRcc[1],
                          Heat_Test_Normal1_Sample2$DQ$stat,Heat_Test_SST1_Sample2$DQ$stat,Heat_Test_SGED1_Sample2$DQ$stat,Heat_Test_JSU1_Sample2$DQ$stat,
                          Heat_Test_Normal1_Sample2$AE,Heat_Test_SST1_Sample2$AE,Heat_Test_SGED1_Sample2$AE,Heat_Test_JSU1_Sample2$AE,
                          (Heat_Test_Normal1_Sample2$Loss$Loss/Heat_Test_Normal1_Sample2$Loss$Loss),(Heat_Test_SST1_Sample2$Loss$Loss/Heat_Test_Normal1_Sample2$Loss$Loss),
                          (Heat_Test_SGED1_Sample2$Loss$Loss/Heat_Test_Normal1_Sample2$Loss$Loss),(Heat_Test_JSU1_Sample2$Loss$Loss/Heat_Test_Normal1_Sample2$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Heat_Sample2) <- c("Normal","SST","SGED","JSU")    
rownames(Heat_Sample2) <- c("UC","CC","DQ","AE","QL")    

Heat_Sample3 <- matrix(c(Heat_Test_Normal1_Sample3$LRuc[1],Heat_Test_SST1_Sample3$LRuc[1],Heat_Test_SGED1_Sample3$LRuc[1],Heat_Test_JSU1_Sample3$LRuc[1],
                          Heat_Test_Normal1_Sample3$LRcc[1],Heat_Test_SST1_Sample3$LRcc[1],Heat_Test_SGED1_Sample3$LRcc[1],Heat_Test_JSU1_Sample3$LRcc[1],
                          Heat_Test_Normal1_Sample3$DQ$stat,Heat_Test_SST1_Sample3$DQ$stat,Heat_Test_SGED1_Sample3$DQ$stat,Heat_Test_JSU1_Sample3$DQ$stat,
                          Heat_Test_Normal1_Sample3$AE,Heat_Test_SST1_Sample3$AE,Heat_Test_SGED1_Sample3$AE,Heat_Test_JSU1_Sample3$AE,
                          (Heat_Test_Normal1_Sample3$Loss$Loss/Heat_Test_Normal1_Sample3$Loss$Loss),(Heat_Test_SST1_Sample3$Loss$Loss/Heat_Test_Normal1_Sample3$Loss$Loss),
                          (Heat_Test_SGED1_Sample3$Loss$Loss/Heat_Test_Normal1_Sample3$Loss$Loss),(Heat_Test_JSU1_Sample3$Loss$Loss/Heat_Test_Normal1_Sample3$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Heat_Sample3) <- c("Normal","SST","SGED","JSU")    
rownames(Heat_Sample3) <- c("UC","CC","DQ","AE","QL")    

Heat_Sample4 <- matrix(c(Heat_Test_Normal1_Sample4$LRuc[1],Heat_Test_SST1_Sample4$LRuc[1],Heat_Test_SGED1_Sample4$LRuc[1],Heat_Test_JSU1_Sample4$LRuc[1],
                          Heat_Test_Normal1_Sample4$LRcc[1],Heat_Test_SST1_Sample4$LRcc[1],Heat_Test_SGED1_Sample4$LRcc[1],Heat_Test_JSU1_Sample4$LRcc[1],
                          Heat_Test_Normal1_Sample4$DQ$stat,Heat_Test_SST1_Sample4$DQ$stat,Heat_Test_SGED1_Sample4$DQ$stat,Heat_Test_JSU1_Sample4$DQ$stat,
                          Heat_Test_Normal1_Sample4$AE,Heat_Test_SST1_Sample4$AE,Heat_Test_SGED1_Sample4$AE,Heat_Test_JSU1_Sample4$AE,
                          (Heat_Test_Normal1_Sample4$Loss$Loss/Heat_Test_Normal1_Sample4$Loss$Loss),(Heat_Test_SST1_Sample4$Loss$Loss/Heat_Test_Normal1_Sample4$Loss$Loss),
                          (Heat_Test_SGED1_Sample4$Loss$Loss/Heat_Test_Normal1_Sample4$Loss$Loss),(Heat_Test_JSU1_Sample4$Loss$Loss/Heat_Test_Normal1_Sample4$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Heat_Sample4) <- c("Normal","SST","SGED","JSU")    
rownames(Heat_Sample4) <- c("UC","CC","DQ","AE","QL")    

Heat_Sample5 <- matrix(c(Heat_Test_Normal1_Sample5$LRuc[1],Heat_Test_SST1_Sample5$LRuc[1],Heat_Test_SGED1_Sample5$LRuc[1],Heat_Test_JSU1_Sample5$LRuc[1],
                          Heat_Test_Normal1_Sample5$LRcc[1],Heat_Test_SST1_Sample5$LRcc[1],Heat_Test_SGED1_Sample5$LRcc[1],Heat_Test_JSU1_Sample5$LRcc[1],
                          Heat_Test_Normal1_Sample5$DQ$stat,Heat_Test_SST1_Sample5$DQ$stat,Heat_Test_SGED1_Sample5$DQ$stat,Heat_Test_JSU1_Sample5$DQ$stat,
                          Heat_Test_Normal1_Sample5$AE,Heat_Test_SST1_Sample5$AE,Heat_Test_SGED1_Sample5$AE,Heat_Test_JSU1_Sample5$AE,
                          (Heat_Test_Normal1_Sample5$Loss$Loss/Heat_Test_Normal1_Sample5$Loss$Loss),(Heat_Test_SST1_Sample5$Loss$Loss/Heat_Test_Normal1_Sample5$Loss$Loss),
                          (Heat_Test_SGED1_Sample5$Loss$Loss/Heat_Test_Normal1_Sample5$Loss$Loss),(Heat_Test_JSU1_Sample5$Loss$Loss/Heat_Test_Normal1_Sample5$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Heat_Sample5) <- c("Normal","SST","SGED","JSU")    
rownames(Heat_Sample5) <- c("UC","CC","DQ","AE","QL")    


Gas_Sample1 <- matrix(c(Gas_Test_Normal1_Sample1$LRuc[1],Gas_Test_SST1_Sample1$LRuc[1],Gas_Test_SGED1_Sample1$LRuc[1],Gas_Test_JSU1_Sample1$LRuc[1],
                          Gas_Test_Normal1_Sample1$LRcc[1],Gas_Test_SST1_Sample1$LRcc[1],Gas_Test_SGED1_Sample1$LRcc[1],Gas_Test_JSU1_Sample1$LRcc[1],
                          Gas_Test_Normal1_Sample1$DQ$stat,Gas_Test_SST1_Sample1$DQ$stat,Gas_Test_SGED1_Sample1$DQ$stat,Gas_Test_JSU1_Sample1$DQ$stat,
                          Gas_Test_Normal1_Sample1$AE,Gas_Test_SST1_Sample1$AE,Gas_Test_SGED1_Sample1$AE,Gas_Test_JSU1_Sample1$AE,
                          (Gas_Test_Normal1_Sample1$Loss$Loss/Gas_Test_Normal1_Sample1$Loss$Loss),(Gas_Test_SST1_Sample1$Loss$Loss/Gas_Test_Normal1_Sample1$Loss$Loss),
                          (Gas_Test_SGED1_Sample1$Loss$Loss/Gas_Test_Normal1_Sample1$Loss$Loss),(Gas_Test_JSU1_Sample1$Loss$Loss/Gas_Test_Normal1_Sample1$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Gas_Sample1) <- c("Normal","SST","SGED","JSU")    
rownames(Gas_Sample1) <- c("UC","CC","DQ","AE","QL")    

Gas_Sample2 <- matrix(c(Gas_Test_Normal1_Sample2$LRuc[1],Gas_Test_SST1_Sample2$LRuc[1],Gas_Test_SGED1_Sample2$LRuc[1],Gas_Test_JSU1_Sample2$LRuc[1],
                          Gas_Test_Normal1_Sample2$LRcc[1],Gas_Test_SST1_Sample2$LRcc[1],Gas_Test_SGED1_Sample2$LRcc[1],Gas_Test_JSU1_Sample2$LRcc[1],
                          Gas_Test_Normal1_Sample2$DQ$stat,Gas_Test_SST1_Sample2$DQ$stat,Gas_Test_SGED1_Sample2$DQ$stat,Gas_Test_JSU1_Sample2$DQ$stat,
                          Gas_Test_Normal1_Sample2$AE,Gas_Test_SST1_Sample2$AE,Gas_Test_SGED1_Sample2$AE,Gas_Test_JSU1_Sample2$AE,
                          (Gas_Test_Normal1_Sample2$Loss$Loss/Gas_Test_Normal1_Sample2$Loss$Loss),(Gas_Test_SST1_Sample2$Loss$Loss/Gas_Test_Normal1_Sample2$Loss$Loss),
                          (Gas_Test_SGED1_Sample2$Loss$Loss/Gas_Test_Normal1_Sample2$Loss$Loss),(Gas_Test_JSU1_Sample2$Loss$Loss/Gas_Test_Normal1_Sample2$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Gas_Sample2) <- c("Normal","SST","SGED","JSU")    
rownames(Gas_Sample2) <- c("UC","CC","DQ","AE","QL")    

Gas_Sample3 <- matrix(c(Gas_Test_Normal1_Sample3$LRuc[1],Gas_Test_SST1_Sample3$LRuc[1],Gas_Test_SGED1_Sample3$LRuc[1],Gas_Test_JSU1_Sample3$LRuc[1],
                          Gas_Test_Normal1_Sample3$LRcc[1],Gas_Test_SST1_Sample3$LRcc[1],Gas_Test_SGED1_Sample3$LRcc[1],Gas_Test_JSU1_Sample3$LRcc[1],
                          Gas_Test_Normal1_Sample3$DQ$stat,Gas_Test_SST1_Sample3$DQ$stat,Gas_Test_SGED1_Sample3$DQ$stat,Gas_Test_JSU1_Sample3$DQ$stat,
                          Gas_Test_Normal1_Sample3$AE,Gas_Test_SST1_Sample3$AE,Gas_Test_SGED1_Sample3$AE,Gas_Test_JSU1_Sample3$AE,
                          (Gas_Test_Normal1_Sample3$Loss$Loss/Gas_Test_Normal1_Sample3$Loss$Loss),(Gas_Test_SST1_Sample3$Loss$Loss/Gas_Test_Normal1_Sample3$Loss$Loss),
                          (Gas_Test_SGED1_Sample3$Loss$Loss/Gas_Test_Normal1_Sample3$Loss$Loss),(Gas_Test_JSU1_Sample3$Loss$Loss/Gas_Test_Normal1_Sample3$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Gas_Sample3) <- c("Normal","SST","SGED","JSU")    
rownames(Gas_Sample3) <- c("UC","CC","DQ","AE","QL")    

Gas_Sample4 <- matrix(c(Gas_Test_Normal1_Sample4$LRuc[1],Gas_Test_SST1_Sample4$LRuc[1],Gas_Test_SGED1_Sample4$LRuc[1],Gas_Test_JSU1_Sample4$LRuc[1],
                          Gas_Test_Normal1_Sample4$LRcc[1],Gas_Test_SST1_Sample4$LRcc[1],Gas_Test_SGED1_Sample4$LRcc[1],Gas_Test_JSU1_Sample4$LRcc[1],
                          Gas_Test_Normal1_Sample4$DQ$stat,Gas_Test_SST1_Sample4$DQ$stat,Gas_Test_SGED1_Sample4$DQ$stat,Gas_Test_JSU1_Sample4$DQ$stat,
                          Gas_Test_Normal1_Sample4$AE,Gas_Test_SST1_Sample4$AE,Gas_Test_SGED1_Sample4$AE,Gas_Test_JSU1_Sample4$AE,
                          (Gas_Test_Normal1_Sample4$Loss$Loss/Gas_Test_Normal1_Sample4$Loss$Loss),(Gas_Test_SST1_Sample4$Loss$Loss/Gas_Test_Normal1_Sample4$Loss$Loss),
                          (Gas_Test_SGED1_Sample4$Loss$Loss/Gas_Test_Normal1_Sample4$Loss$Loss),(Gas_Test_JSU1_Sample4$Loss$Loss/Gas_Test_Normal1_Sample4$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Gas_Sample4) <- c("Normal","SST","SGED","JSU")    
rownames(Gas_Sample4) <- c("UC","CC","DQ","AE","QL")    

Gas_Sample5 <- matrix(c(Gas_Test_Normal1_Sample5$LRuc[1],Gas_Test_SST1_Sample5$LRuc[1],Gas_Test_SGED1_Sample5$LRuc[1],Gas_Test_JSU1_Sample5$LRuc[1],
                          Gas_Test_Normal1_Sample5$LRcc[1],Gas_Test_SST1_Sample5$LRcc[1],Gas_Test_SGED1_Sample5$LRcc[1],Gas_Test_JSU1_Sample5$LRcc[1],
                          Gas_Test_Normal1_Sample5$DQ$stat,Gas_Test_SST1_Sample5$DQ$stat,Gas_Test_SGED1_Sample5$DQ$stat,Gas_Test_JSU1_Sample5$DQ$stat,
                          Gas_Test_Normal1_Sample5$AE,Gas_Test_SST1_Sample5$AE,Gas_Test_SGED1_Sample5$AE,Gas_Test_JSU1_Sample5$AE,
                          (Gas_Test_Normal1_Sample5$Loss$Loss/Gas_Test_Normal1_Sample5$Loss$Loss),(Gas_Test_SST1_Sample5$Loss$Loss/Gas_Test_Normal1_Sample5$Loss$Loss),
                          (Gas_Test_SGED1_Sample5$Loss$Loss/Gas_Test_Normal1_Sample5$Loss$Loss),(Gas_Test_JSU1_Sample5$Loss$Loss/Gas_Test_Normal1_Sample5$Loss$Loss)),nrow=5,byrow = TRUE)
colnames(Gas_Sample5) <- c("Normal","SST","SGED","JSU")    
rownames(Gas_Sample5) <- c("UC","CC","DQ","AE","QL")   

### Storing the p-values for subsamples 

# 1 percent
Heat_pvalues_1pc_Sample1 <- matrix(c(Heat_Test_Normal1_Sample1$LRuc[2],Heat_Test_SST1_Sample1$LRuc[2],Heat_Test_SGED1_Sample1$LRuc[2],Heat_Test_JSU1_Sample1$LRuc[2],
                          Heat_Test_Normal1_Sample1$LRcc[2],Heat_Test_SST1_Sample1$LRcc[2],Heat_Test_SGED1_Sample1$LRcc[2],Heat_Test_JSU1_Sample1$LRcc[2],
                          Heat_Test_Normal1_Sample1$DQ$pvalue,Heat_Test_SST1_Sample1$DQ$pvalue,Heat_Test_SGED1_Sample1$DQ$pvalue,Heat_Test_JSU1_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_1pc_Sample2 <- matrix(c(Heat_Test_Normal1_Sample2$LRuc[2],Heat_Test_SST1_Sample2$LRuc[2],Heat_Test_SGED1_Sample2$LRuc[2],Heat_Test_JSU1_Sample2$LRuc[2],
                                      Heat_Test_Normal1_Sample2$LRcc[2],Heat_Test_SST1_Sample2$LRcc[2],Heat_Test_SGED1_Sample2$LRcc[2],Heat_Test_JSU1_Sample2$LRcc[2],
                                      Heat_Test_Normal1_Sample2$DQ$pvalue,Heat_Test_SST1_Sample2$DQ$pvalue,Heat_Test_SGED1_Sample2$DQ$pvalue,Heat_Test_JSU1_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_1pc_Sample3 <- matrix(c(Heat_Test_Normal1_Sample3$LRuc[2],Heat_Test_SST1_Sample3$LRuc[2],Heat_Test_SGED1_Sample3$LRuc[2],Heat_Test_JSU1_Sample3$LRuc[2],
                                      Heat_Test_Normal1_Sample3$LRcc[2],Heat_Test_SST1_Sample3$LRcc[2],Heat_Test_SGED1_Sample3$LRcc[2],Heat_Test_JSU1_Sample3$LRcc[2],
                                      Heat_Test_Normal1_Sample3$DQ$pvalue,Heat_Test_SST1_Sample3$DQ$pvalue,Heat_Test_SGED1_Sample3$DQ$pvalue,Heat_Test_JSU1_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_1pc_Sample4 <- matrix(c(Heat_Test_Normal1_Sample4$LRuc[2],Heat_Test_SST1_Sample4$LRuc[2],Heat_Test_SGED1_Sample4$LRuc[2],Heat_Test_JSU1_Sample4$LRuc[2],
                                      Heat_Test_Normal1_Sample4$LRcc[2],Heat_Test_SST1_Sample4$LRcc[2],Heat_Test_SGED1_Sample4$LRcc[2],Heat_Test_JSU1_Sample4$LRcc[2],
                                      Heat_Test_Normal1_Sample4$DQ$pvalue,Heat_Test_SST1_Sample4$DQ$pvalue,Heat_Test_SGED1_Sample4$DQ$pvalue,Heat_Test_JSU1_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_1pc_Sample5 <- matrix(c(Heat_Test_Normal1_Sample5$LRuc[2],Heat_Test_SST1_Sample5$LRuc[2],Heat_Test_SGED1_Sample5$LRuc[2],Heat_Test_JSU1_Sample5$LRuc[2],
                                      Heat_Test_Normal1_Sample5$LRcc[2],Heat_Test_SST1_Sample5$LRcc[2],Heat_Test_SGED1_Sample5$LRcc[2],Heat_Test_JSU1_Sample5$LRcc[2],
                                      Heat_Test_Normal1_Sample5$DQ$pvalue,Heat_Test_SST1_Sample5$DQ$pvalue,Heat_Test_SGED1_Sample5$DQ$pvalue,Heat_Test_JSU1_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)


Gas_pvalues_1pc_Sample1 <- matrix(c(Gas_Test_Normal1_Sample1$LRuc[2],Gas_Test_SST1_Sample1$LRuc[2],Gas_Test_SGED1_Sample1$LRuc[2],Gas_Test_JSU1_Sample1$LRuc[2],
                                      Gas_Test_Normal1_Sample1$LRcc[2],Gas_Test_SST1_Sample1$LRcc[2],Gas_Test_SGED1_Sample1$LRcc[2],Gas_Test_JSU1_Sample1$LRcc[2],
                                      Gas_Test_Normal1_Sample1$DQ$pvalue,Gas_Test_SST1_Sample1$DQ$pvalue,Gas_Test_SGED1_Sample1$DQ$pvalue,Gas_Test_JSU1_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_1pc_Sample2 <- matrix(c(Gas_Test_Normal1_Sample2$LRuc[2],Gas_Test_SST1_Sample2$LRuc[2],Gas_Test_SGED1_Sample2$LRuc[2],Gas_Test_JSU1_Sample2$LRuc[2],
                                      Gas_Test_Normal1_Sample2$LRcc[2],Gas_Test_SST1_Sample2$LRcc[2],Gas_Test_SGED1_Sample2$LRcc[2],Gas_Test_JSU1_Sample2$LRcc[2],
                                      Gas_Test_Normal1_Sample2$DQ$pvalue,Gas_Test_SST1_Sample2$DQ$pvalue,Gas_Test_SGED1_Sample2$DQ$pvalue,Gas_Test_JSU1_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_1pc_Sample3 <- matrix(c(Gas_Test_Normal1_Sample3$LRuc[2],Gas_Test_SST1_Sample3$LRuc[2],Gas_Test_SGED1_Sample3$LRuc[2],Gas_Test_JSU1_Sample3$LRuc[2],
                                      Gas_Test_Normal1_Sample3$LRcc[2],Gas_Test_SST1_Sample3$LRcc[2],Gas_Test_SGED1_Sample3$LRcc[2],Gas_Test_JSU1_Sample3$LRcc[2],
                                      Gas_Test_Normal1_Sample3$DQ$pvalue,Gas_Test_SST1_Sample3$DQ$pvalue,Gas_Test_SGED1_Sample3$DQ$pvalue,Gas_Test_JSU1_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_1pc_Sample4 <- matrix(c(Gas_Test_Normal1_Sample4$LRuc[2],Gas_Test_SST1_Sample4$LRuc[2],Gas_Test_SGED1_Sample4$LRuc[2],Gas_Test_JSU1_Sample4$LRuc[2],
                                      Gas_Test_Normal1_Sample4$LRcc[2],Gas_Test_SST1_Sample4$LRcc[2],Gas_Test_SGED1_Sample4$LRcc[2],Gas_Test_JSU1_Sample4$LRcc[2],
                                      Gas_Test_Normal1_Sample4$DQ$pvalue,Gas_Test_SST1_Sample4$DQ$pvalue,Gas_Test_SGED1_Sample4$DQ$pvalue,Gas_Test_JSU1_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_1pc_Sample5 <- matrix(c(Gas_Test_Normal1_Sample5$LRuc[2],Gas_Test_SST1_Sample5$LRuc[2],Gas_Test_SGED1_Sample5$LRuc[2],Gas_Test_JSU1_Sample5$LRuc[2],
                                      Gas_Test_Normal1_Sample5$LRcc[2],Gas_Test_SST1_Sample5$LRcc[2],Gas_Test_SGED1_Sample5$LRcc[2],Gas_Test_JSU1_Sample5$LRcc[2],
                                      Gas_Test_Normal1_Sample5$DQ$pvalue,Gas_Test_SST1_Sample5$DQ$pvalue,Gas_Test_SGED1_Sample5$DQ$pvalue,Gas_Test_JSU1_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)

# 5 percent

Heat_pvalues_5pc_Sample1 <- matrix(c(Heat_Test_Normal5_Sample1$LRuc[2],Heat_Test_SST5_Sample1$LRuc[2],Heat_Test_SGED5_Sample1$LRuc[2],Heat_Test_JSU5_Sample1$LRuc[2],
                                      Heat_Test_Normal5_Sample1$LRcc[2],Heat_Test_SST5_Sample1$LRcc[2],Heat_Test_SGED5_Sample1$LRcc[2],Heat_Test_JSU5_Sample1$LRcc[2],
                                      Heat_Test_Normal5_Sample1$DQ$pvalue,Heat_Test_SST5_Sample1$DQ$pvalue,Heat_Test_SGED5_Sample1$DQ$pvalue,Heat_Test_JSU5_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_5pc_Sample2 <- matrix(c(Heat_Test_Normal5_Sample2$LRuc[2],Heat_Test_SST5_Sample2$LRuc[2],Heat_Test_SGED5_Sample2$LRuc[2],Heat_Test_JSU5_Sample2$LRuc[2],
                                      Heat_Test_Normal5_Sample2$LRcc[2],Heat_Test_SST5_Sample2$LRcc[2],Heat_Test_SGED5_Sample2$LRcc[2],Heat_Test_JSU5_Sample2$LRcc[2],
                                      Heat_Test_Normal5_Sample2$DQ$pvalue,Heat_Test_SST5_Sample2$DQ$pvalue,Heat_Test_SGED5_Sample2$DQ$pvalue,Heat_Test_JSU5_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_5pc_Sample3 <- matrix(c(Heat_Test_Normal5_Sample3$LRuc[2],Heat_Test_SST5_Sample3$LRuc[2],Heat_Test_SGED5_Sample3$LRuc[2],Heat_Test_JSU5_Sample3$LRuc[2],
                                      Heat_Test_Normal5_Sample3$LRcc[2],Heat_Test_SST5_Sample3$LRcc[2],Heat_Test_SGED5_Sample3$LRcc[2],Heat_Test_JSU5_Sample3$LRcc[2],
                                      Heat_Test_Normal5_Sample3$DQ$pvalue,Heat_Test_SST5_Sample3$DQ$pvalue,Heat_Test_SGED5_Sample3$DQ$pvalue,Heat_Test_JSU5_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_5pc_Sample4 <- matrix(c(Heat_Test_Normal5_Sample4$LRuc[2],Heat_Test_SST5_Sample4$LRuc[2],Heat_Test_SGED5_Sample4$LRuc[2],Heat_Test_JSU5_Sample4$LRuc[2],
                                      Heat_Test_Normal5_Sample4$LRcc[2],Heat_Test_SST5_Sample4$LRcc[2],Heat_Test_SGED5_Sample4$LRcc[2],Heat_Test_JSU5_Sample4$LRcc[2],
                                      Heat_Test_Normal5_Sample4$DQ$pvalue,Heat_Test_SST5_Sample4$DQ$pvalue,Heat_Test_SGED5_Sample4$DQ$pvalue,Heat_Test_JSU5_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_5pc_Sample5 <- matrix(c(Heat_Test_Normal5_Sample5$LRuc[2],Heat_Test_SST5_Sample5$LRuc[2],Heat_Test_SGED5_Sample5$LRuc[2],Heat_Test_JSU5_Sample5$LRuc[2],
                                      Heat_Test_Normal5_Sample5$LRcc[2],Heat_Test_SST5_Sample5$LRcc[2],Heat_Test_SGED5_Sample5$LRcc[2],Heat_Test_JSU5_Sample5$LRcc[2],
                                      Heat_Test_Normal5_Sample5$DQ$pvalue,Heat_Test_SST5_Sample5$DQ$pvalue,Heat_Test_SGED5_Sample5$DQ$pvalue,Heat_Test_JSU5_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)


Gas_pvalues_5pc_Sample1 <- matrix(c(Gas_Test_Normal5_Sample1$LRuc[2],Gas_Test_SST5_Sample1$LRuc[2],Gas_Test_SGED5_Sample1$LRuc[2],Gas_Test_JSU5_Sample1$LRuc[2],
                                    Gas_Test_Normal5_Sample1$LRcc[2],Gas_Test_SST5_Sample1$LRcc[2],Gas_Test_SGED5_Sample1$LRcc[2],Gas_Test_JSU5_Sample1$LRcc[2],
                                    Gas_Test_Normal5_Sample1$DQ$pvalue,Gas_Test_SST5_Sample1$DQ$pvalue,Gas_Test_SGED5_Sample1$DQ$pvalue,Gas_Test_JSU5_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_5pc_Sample2 <- matrix(c(Gas_Test_Normal5_Sample2$LRuc[2],Gas_Test_SST5_Sample2$LRuc[2],Gas_Test_SGED5_Sample2$LRuc[2],Gas_Test_JSU5_Sample2$LRuc[2],
                                    Gas_Test_Normal5_Sample2$LRcc[2],Gas_Test_SST5_Sample2$LRcc[2],Gas_Test_SGED5_Sample2$LRcc[2],Gas_Test_JSU5_Sample2$LRcc[2],
                                    Gas_Test_Normal5_Sample2$DQ$pvalue,Gas_Test_SST5_Sample2$DQ$pvalue,Gas_Test_SGED5_Sample2$DQ$pvalue,Gas_Test_JSU5_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_5pc_Sample3 <- matrix(c(Gas_Test_Normal5_Sample3$LRuc[2],Gas_Test_SST5_Sample3$LRuc[2],Gas_Test_SGED5_Sample3$LRuc[2],Gas_Test_JSU5_Sample3$LRuc[2],
                                    Gas_Test_Normal5_Sample3$LRcc[2],Gas_Test_SST5_Sample3$LRcc[2],Gas_Test_SGED5_Sample3$LRcc[2],Gas_Test_JSU5_Sample3$LRcc[2],
                                    Gas_Test_Normal5_Sample3$DQ$pvalue,Gas_Test_SST5_Sample3$DQ$pvalue,Gas_Test_SGED5_Sample3$DQ$pvalue,Gas_Test_JSU5_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_5pc_Sample4 <- matrix(c(Gas_Test_Normal5_Sample4$LRuc[2],Gas_Test_SST5_Sample4$LRuc[2],Gas_Test_SGED5_Sample4$LRuc[2],Gas_Test_JSU5_Sample4$LRuc[2],
                                    Gas_Test_Normal5_Sample4$LRcc[2],Gas_Test_SST5_Sample4$LRcc[2],Gas_Test_SGED5_Sample4$LRcc[2],Gas_Test_JSU5_Sample4$LRcc[2],
                                    Gas_Test_Normal5_Sample4$DQ$pvalue,Gas_Test_SST5_Sample4$DQ$pvalue,Gas_Test_SGED5_Sample4$DQ$pvalue,Gas_Test_JSU5_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_5pc_Sample5 <- matrix(c(Gas_Test_Normal5_Sample5$LRuc[2],Gas_Test_SST5_Sample5$LRuc[2],Gas_Test_SGED5_Sample5$LRuc[2],Gas_Test_JSU5_Sample5$LRuc[2],
                                    Gas_Test_Normal5_Sample5$LRcc[2],Gas_Test_SST5_Sample5$LRcc[2],Gas_Test_SGED5_Sample5$LRcc[2],Gas_Test_JSU5_Sample5$LRcc[2],
                                    Gas_Test_Normal5_Sample5$DQ$pvalue,Gas_Test_SST5_Sample5$DQ$pvalue,Gas_Test_SGED5_Sample5$DQ$pvalue,Gas_Test_JSU5_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)

# 10 percent

Heat_pvalues_10pc_Sample1 <- matrix(c(Heat_Test_Normal10_Sample1$LRuc[2],Heat_Test_SST10_Sample1$LRuc[2],Heat_Test_SGED10_Sample1$LRuc[2],Heat_Test_JSU10_Sample1$LRuc[2],
                                      Heat_Test_Normal10_Sample1$LRcc[2],Heat_Test_SST10_Sample1$LRcc[2],Heat_Test_SGED10_Sample1$LRcc[2],Heat_Test_JSU10_Sample1$LRcc[2],
                                      Heat_Test_Normal10_Sample1$DQ$pvalue,Heat_Test_SST10_Sample1$DQ$pvalue,Heat_Test_SGED10_Sample1$DQ$pvalue,Heat_Test_JSU10_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_10pc_Sample2 <- matrix(c(Heat_Test_Normal10_Sample2$LRuc[2],Heat_Test_SST10_Sample2$LRuc[2],Heat_Test_SGED10_Sample2$LRuc[2],Heat_Test_JSU10_Sample2$LRuc[2],
                                      Heat_Test_Normal10_Sample2$LRcc[2],Heat_Test_SST10_Sample2$LRcc[2],Heat_Test_SGED10_Sample2$LRcc[2],Heat_Test_JSU10_Sample2$LRcc[2],
                                      Heat_Test_Normal10_Sample2$DQ$pvalue,Heat_Test_SST10_Sample2$DQ$pvalue,Heat_Test_SGED10_Sample2$DQ$pvalue,Heat_Test_JSU10_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_10pc_Sample3 <- matrix(c(Heat_Test_Normal10_Sample3$LRuc[2],Heat_Test_SST10_Sample3$LRuc[2],Heat_Test_SGED10_Sample3$LRuc[2],Heat_Test_JSU10_Sample3$LRuc[2],
                                      Heat_Test_Normal10_Sample3$LRcc[2],Heat_Test_SST10_Sample3$LRcc[2],Heat_Test_SGED10_Sample3$LRcc[2],Heat_Test_JSU10_Sample3$LRcc[2],
                                      Heat_Test_Normal10_Sample3$DQ$pvalue,Heat_Test_SST10_Sample3$DQ$pvalue,Heat_Test_SGED10_Sample3$DQ$pvalue,Heat_Test_JSU10_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_10pc_Sample4 <- matrix(c(Heat_Test_Normal10_Sample4$LRuc[2],Heat_Test_SST10_Sample4$LRuc[2],Heat_Test_SGED10_Sample4$LRuc[2],Heat_Test_JSU10_Sample4$LRuc[2],
                                      Heat_Test_Normal10_Sample4$LRcc[2],Heat_Test_SST10_Sample4$LRcc[2],Heat_Test_SGED10_Sample4$LRcc[2],Heat_Test_JSU10_Sample4$LRcc[2],
                                      Heat_Test_Normal10_Sample4$DQ$pvalue,Heat_Test_SST10_Sample4$DQ$pvalue,Heat_Test_SGED10_Sample4$DQ$pvalue,Heat_Test_JSU10_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Heat_pvalues_10pc_Sample5 <- matrix(c(Heat_Test_Normal10_Sample5$LRuc[2],Heat_Test_SST10_Sample5$LRuc[2],Heat_Test_SGED10_Sample5$LRuc[2],Heat_Test_JSU10_Sample5$LRuc[2],
                                      Heat_Test_Normal10_Sample5$LRcc[2],Heat_Test_SST10_Sample5$LRcc[2],Heat_Test_SGED10_Sample5$LRcc[2],Heat_Test_JSU10_Sample5$LRcc[2],
                                      Heat_Test_Normal10_Sample5$DQ$pvalue,Heat_Test_SST10_Sample5$DQ$pvalue,Heat_Test_SGED10_Sample5$DQ$pvalue,Heat_Test_JSU10_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)


Gas_pvalues_10pc_Sample1 <- matrix(c(Gas_Test_Normal10_Sample1$LRuc[2],Gas_Test_SST10_Sample1$LRuc[2],Gas_Test_SGED10_Sample1$LRuc[2],Gas_Test_JSU10_Sample1$LRuc[2],
                                    Gas_Test_Normal10_Sample1$LRcc[2],Gas_Test_SST10_Sample1$LRcc[2],Gas_Test_SGED10_Sample1$LRcc[2],Gas_Test_JSU10_Sample1$LRcc[2],
                                    Gas_Test_Normal10_Sample1$DQ$pvalue,Gas_Test_SST10_Sample1$DQ$pvalue,Gas_Test_SGED10_Sample1$DQ$pvalue,Gas_Test_JSU10_Sample1$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_10pc_Sample2 <- matrix(c(Gas_Test_Normal10_Sample2$LRuc[2],Gas_Test_SST10_Sample2$LRuc[2],Gas_Test_SGED10_Sample2$LRuc[2],Gas_Test_JSU10_Sample2$LRuc[2],
                                    Gas_Test_Normal10_Sample2$LRcc[2],Gas_Test_SST10_Sample2$LRcc[2],Gas_Test_SGED10_Sample2$LRcc[2],Gas_Test_JSU10_Sample2$LRcc[2],
                                    Gas_Test_Normal10_Sample2$DQ$pvalue,Gas_Test_SST10_Sample2$DQ$pvalue,Gas_Test_SGED10_Sample2$DQ$pvalue,Gas_Test_JSU10_Sample2$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_10pc_Sample3 <- matrix(c(Gas_Test_Normal10_Sample3$LRuc[2],Gas_Test_SST10_Sample3$LRuc[2],Gas_Test_SGED10_Sample3$LRuc[2],Gas_Test_JSU10_Sample3$LRuc[2],
                                    Gas_Test_Normal10_Sample3$LRcc[2],Gas_Test_SST10_Sample3$LRcc[2],Gas_Test_SGED10_Sample3$LRcc[2],Gas_Test_JSU10_Sample3$LRcc[2],
                                    Gas_Test_Normal10_Sample3$DQ$pvalue,Gas_Test_SST10_Sample3$DQ$pvalue,Gas_Test_SGED10_Sample3$DQ$pvalue,Gas_Test_JSU10_Sample3$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_10pc_Sample4 <- matrix(c(Gas_Test_Normal10_Sample4$LRuc[2],Gas_Test_SST10_Sample4$LRuc[2],Gas_Test_SGED10_Sample4$LRuc[2],Gas_Test_JSU10_Sample4$LRuc[2],
                                    Gas_Test_Normal10_Sample4$LRcc[2],Gas_Test_SST10_Sample4$LRcc[2],Gas_Test_SGED10_Sample4$LRcc[2],Gas_Test_JSU10_Sample4$LRcc[2],
                                    Gas_Test_Normal10_Sample4$DQ$pvalue,Gas_Test_SST10_Sample4$DQ$pvalue,Gas_Test_SGED10_Sample4$DQ$pvalue,Gas_Test_JSU10_Sample4$DQ$pvalue),ncol=4,byrow = TRUE)

Gas_pvalues_10pc_Sample5 <- matrix(c(Gas_Test_Normal10_Sample5$LRuc[2],Gas_Test_SST10_Sample5$LRuc[2],Gas_Test_SGED10_Sample5$LRuc[2],Gas_Test_JSU10_Sample5$LRuc[2],
                                    Gas_Test_Normal10_Sample5$LRcc[2],Gas_Test_SST10_Sample5$LRcc[2],Gas_Test_SGED10_Sample5$LRcc[2],Gas_Test_JSU10_Sample5$LRcc[2],
                                    Gas_Test_Normal10_Sample5$DQ$pvalue,Gas_Test_SST10_Sample5$DQ$pvalue,Gas_Test_SGED10_Sample5$DQ$pvalue,Gas_Test_JSU10_Sample5$DQ$pvalue),ncol=4,byrow = TRUE)


Heat_pvalues_1pc   <- rbind(Heat_pvalues_1pc_Sample1,Heat_pvalues_1pc_Sample2,Heat_pvalues_1pc_Sample3,Heat_pvalues_1pc_Sample4,Heat_pvalues_1pc_Sample5)
Heat_pvalues_5pc   <- rbind(Heat_pvalues_5pc_Sample1,Heat_pvalues_5pc_Sample2,Heat_pvalues_5pc_Sample3,Heat_pvalues_5pc_Sample4,Heat_pvalues_5pc_Sample5)
Heat_pvalues_10pc  <- rbind(Heat_pvalues_10pc_Sample1,Heat_pvalues_10pc_Sample2,Heat_pvalues_10pc_Sample3,Heat_pvalues_10pc_Sample4,Heat_pvalues_10pc_Sample5)

Gas_pvalues_1pc   <- rbind(Gas_pvalues_1pc_Sample1,Gas_pvalues_1pc_Sample2,Gas_pvalues_1pc_Sample3,Gas_pvalues_1pc_Sample4,Gas_pvalues_1pc_Sample5)
Gas_pvalues_5pc   <- rbind(Gas_pvalues_5pc_Sample1,Gas_pvalues_5pc_Sample2,Gas_pvalues_5pc_Sample3,Gas_pvalues_5pc_Sample4,Gas_pvalues_5pc_Sample5)
Gas_pvalues_10pc  <- rbind(Gas_pvalues_10pc_Sample1,Gas_pvalues_10pc_Sample2,Gas_pvalues_10pc_Sample3,Gas_pvalues_10pc_Sample4,Gas_pvalues_10pc_Sample5)

rownames(Heat_pvalues_1pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Heat_pvalues_1pc) <- c("Normal","SST","SGED","JSU")

rownames(Heat_pvalues_5pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Heat_pvalues_5pc) <- c("Normal","SST","SGED","JSU")

rownames(Heat_pvalues_10pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Heat_pvalues_10pc) <- c("Normal","SST","SGED","JSU")

rownames(Gas_pvalues_1pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Gas_pvalues_1pc) <- c("Normal","SST","SGED","JSU")

rownames(Gas_pvalues_5pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Gas_pvalues_5pc) <- c("Normal","SST","SGED","JSU")

rownames(Gas_pvalues_10pc) <- c(rep(c("UC","CC","DQ"),5))
colnames(Gas_pvalues_10pc) <- c("Normal","SST","SGED","JSU")

#### Repeating for Pearson Distribution  #####

#signific_levels <- c(0.1,0.05,0.01)
## Store the volatility and mu forecast for Sample 1

Holdout_Sample1_Heat_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample1_Gas_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample1_Heat_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample1_Gas_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample1_Heat_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample1_Gas_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample1_Heat_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample1_Gas_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample1_Heat_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample1_Gas_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample1_Heat_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample1_Gas_Pearson_pvalues) <- c("1%","5%","10%")

Sample1_Heat_Holdout_Mu <- as.data.frame(Heat_Roll_Normal_Sample1)[,'Mu']
Sample1_Heat_Holdout_Sigma <- as.data.frame(Heat_Roll_Normal_Sample1)[,'Sigma']
Sample1_Gas_Holdout_Mu <- as.data.frame(Gas_Roll_Normal_Sample1)[,'Mu']
Sample1_Gas_Holdout_Sigma <- as.data.frame(Gas_Roll_Normal_Sample1)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample1_Heat <- qpearson(signific_levels[j],Heat_Param_Pearson_Sample1)
  Z_Pearson_Sample1_Gas   <- qpearson(signific_levels[j],Gas_Param_Pearson_Sample1)
  
  Sample1_Heat_Holdout_VaR <- Sample1_Heat_Holdout_Mu + Sample1_Heat_Holdout_Sigma*Z_Pearson_Sample1_Heat
  Sample1_Gas_Holdout_VaR   <- Sample1_Gas_Holdout_Mu + Sample1_Gas_Holdout_Sigma*Z_Pearson_Sample1_Gas
  Sample1_Heat_Violations <- BacktestVaR(Returns_Heat[(Heat_samples[1]-250+1):(Heat_samples[1]-1)],Sample1_Heat_Holdout_VaR[-1],signific_levels[j])
  Sample1_Gas_Violations   <- BacktestVaR(Returns_Gas[(Gas_samples[1]-250+1):(Gas_samples[1]-1)],Sample1_Gas_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample1_Heat_Pearson[j,] <- c(Sample1_Heat_Violations$LRuc[1],Sample1_Heat_Violations$LRcc[1],Sample1_Heat_Violations$DQ$stat,Sample1_Heat_Violations$AE,Sample1_Heat_Violations$Loss$Loss)
  Holdout_Sample1_Gas_Pearson[j,]   <- c(Sample1_Gas_Violations$LRuc[1],Sample1_Gas_Violations$LRcc[1],Sample1_Gas_Violations$DQ$stat,Sample1_Gas_Violations$AE,Sample1_Gas_Violations$Loss$Loss)
  Holdout_Sample1_Heat_Pearson_pvalues[j,] <- c(Sample1_Heat_Violations$LRuc[2],Sample1_Heat_Violations$LRcc[2],Sample1_Heat_Violations$DQ$pvalue)
  Holdout_Sample1_Gas_Pearson_pvalues[j,] <- c(Sample1_Gas_Violations$LRuc[2],Sample1_Gas_Violations$LRcc[2],Sample1_Gas_Violations$DQ$pvalue)
}


## Store the volatility and mu forecast for Sample 2

Holdout_Sample2_Heat_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample2_Gas_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample2_Heat_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample2_Gas_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample2_Heat_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample2_Gas_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample2_Heat_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample2_Gas_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample2_Heat_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample2_Gas_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample2_Heat_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample2_Gas_Pearson_pvalues) <- c("1%","5%","10%")

Sample2_Heat_Holdout_Mu <- as.data.frame(Heat_Roll_Normal_Sample2)[,'Mu']
Sample2_Heat_Holdout_Sigma <- as.data.frame(Heat_Roll_Normal_Sample2)[,'Sigma']
Sample2_Gas_Holdout_Mu <- as.data.frame(Gas_Roll_Normal_Sample2)[,'Mu']
Sample2_Gas_Holdout_Sigma <- as.data.frame(Gas_Roll_Normal_Sample2)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample2_Heat <- qpearson(signific_levels[j],Heat_Param_Pearson_Sample2)
  Z_Pearson_Sample2_Gas   <- qpearson(signific_levels[j],Gas_Param_Pearson_Sample2)
  
  Sample2_Heat_Holdout_VaR <- Sample2_Heat_Holdout_Mu + Sample2_Heat_Holdout_Sigma*Z_Pearson_Sample2_Heat
  Sample2_Gas_Holdout_VaR   <- Sample2_Gas_Holdout_Mu + Sample2_Gas_Holdout_Sigma*Z_Pearson_Sample2_Gas
  Sample2_Heat_Violations <- BacktestVaR(Returns_Heat[(Heat_samples[2]-250+1):(Heat_samples[2]-1)],Sample2_Heat_Holdout_VaR[-1],signific_levels[j])
  Sample2_Gas_Violations   <- BacktestVaR(Returns_Gas[(Gas_samples[2]-250+1):(Gas_samples[2]-1)],Sample2_Gas_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample2_Heat_Pearson[j,] <- c(Sample2_Heat_Violations$LRuc[1],Sample2_Heat_Violations$LRcc[1],Sample2_Heat_Violations$DQ$stat,Sample2_Heat_Violations$AE,Sample2_Heat_Violations$Loss$Loss)
  Holdout_Sample2_Gas_Pearson[j,]   <- c(Sample2_Gas_Violations$LRuc[1],Sample2_Gas_Violations$LRcc[1],Sample2_Gas_Violations$DQ$stat,Sample2_Gas_Violations$AE,Sample2_Gas_Violations$Loss$Loss)
  Holdout_Sample2_Heat_Pearson_pvalues[j,] <- c(Sample2_Heat_Violations$LRuc[2],Sample2_Heat_Violations$LRcc[2],Sample2_Heat_Violations$DQ$pvalue)
  Holdout_Sample2_Gas_Pearson_pvalues[j,] <- c(Sample2_Gas_Violations$LRuc[2],Sample2_Gas_Violations$LRcc[2],Sample2_Gas_Violations$DQ$pvalue)
}

## Store the volatility and mu forecast for Sample 3

Holdout_Sample3_Heat_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample3_Gas_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample3_Heat_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample3_Gas_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample3_Heat_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample3_Gas_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample3_Heat_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample3_Gas_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample3_Heat_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample3_Gas_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample3_Heat_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample3_Gas_Pearson_pvalues) <- c("1%","5%","10%")

Sample3_Heat_Holdout_Mu <- as.data.frame(Heat_Roll_Normal_Sample3)[,'Mu']
Sample3_Heat_Holdout_Sigma <- as.data.frame(Heat_Roll_Normal_Sample3)[,'Sigma']
Sample3_Gas_Holdout_Mu <- as.data.frame(Gas_Roll_Normal_Sample3)[,'Mu']
Sample3_Gas_Holdout_Sigma <- as.data.frame(Gas_Roll_Normal_Sample3)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample3_Heat <- qpearson(signific_levels[j],Heat_Param_Pearson_Sample3)
  Z_Pearson_Sample3_Gas   <- qpearson(signific_levels[j],Gas_Param_Pearson_Sample3)
  
  Sample3_Heat_Holdout_VaR <- Sample3_Heat_Holdout_Mu + Sample3_Heat_Holdout_Sigma*Z_Pearson_Sample3_Heat
  Sample3_Gas_Holdout_VaR   <- Sample3_Gas_Holdout_Mu + Sample3_Gas_Holdout_Sigma*Z_Pearson_Sample3_Gas
  Sample3_Heat_Violations <- BacktestVaR(Returns_Heat[(Heat_samples[3]-250+1):(Heat_samples[3]-1)],Sample3_Heat_Holdout_VaR[-1],signific_levels[j])
  Sample3_Gas_Violations   <- BacktestVaR(Returns_Gas[(Gas_samples[3]-250+1):(Gas_samples[3]-1)],Sample3_Gas_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample3_Heat_Pearson[j,] <- c(Sample3_Heat_Violations$LRuc[1],Sample3_Heat_Violations$LRcc[1],Sample3_Heat_Violations$DQ$stat,Sample3_Heat_Violations$AE,Sample3_Heat_Violations$Loss$Loss)
  Holdout_Sample3_Gas_Pearson[j,]   <- c(Sample3_Gas_Violations$LRuc[1],Sample3_Gas_Violations$LRcc[1],Sample3_Gas_Violations$DQ$stat,Sample3_Gas_Violations$AE,Sample3_Gas_Violations$Loss$Loss)
  Holdout_Sample3_Heat_Pearson_pvalues[j,] <- c(Sample3_Heat_Violations$LRuc[2],Sample3_Heat_Violations$LRcc[2],Sample3_Heat_Violations$DQ$pvalue)
  Holdout_Sample3_Gas_Pearson_pvalues[j,] <- c(Sample3_Gas_Violations$LRuc[2],Sample3_Gas_Violations$LRcc[2],Sample3_Gas_Violations$DQ$pvalue)
}

## Store the volatility and mu forecast for Sample 4

Holdout_Sample4_Heat_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample4_Gas_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample4_Heat_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample4_Gas_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample4_Heat_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample4_Gas_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample4_Heat_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample4_Gas_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample4_Heat_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample4_Gas_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample4_Heat_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample4_Gas_Pearson_pvalues) <- c("1%","5%","10%")

Sample4_Heat_Holdout_Mu <- as.data.frame(Heat_Roll_Normal_Sample4)[,'Mu']
Sample4_Heat_Holdout_Sigma <- as.data.frame(Heat_Roll_Normal_Sample4)[,'Sigma']
Sample4_Gas_Holdout_Mu <- as.data.frame(Gas_Roll_Normal_Sample4)[,'Mu']
Sample4_Gas_Holdout_Sigma <- as.data.frame(Gas_Roll_Normal_Sample4)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample4_Heat <- qpearson(signific_levels[j],Heat_Param_Pearson_Sample4)
  Z_Pearson_Sample4_Gas   <- qpearson(signific_levels[j],Gas_Param_Pearson_Sample4)
  
  Sample4_Heat_Holdout_VaR <- Sample4_Heat_Holdout_Mu + Sample4_Heat_Holdout_Sigma*Z_Pearson_Sample4_Heat
  Sample4_Gas_Holdout_VaR   <- Sample4_Gas_Holdout_Mu + Sample4_Gas_Holdout_Sigma*Z_Pearson_Sample4_Gas
  Sample4_Heat_Violations <- BacktestVaR(Returns_Heat[(Heat_samples[4]-250+1):(Heat_samples[4]-1)],Sample4_Heat_Holdout_VaR[-1],signific_levels[j])
  Sample4_Gas_Violations   <- BacktestVaR(Returns_Gas[(Gas_samples[4]-250+1):(Gas_samples[4]-1)],Sample4_Gas_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample4_Heat_Pearson[j,] <- c(Sample4_Heat_Violations$LRuc[1],Sample4_Heat_Violations$LRcc[1],Sample4_Heat_Violations$DQ$stat,Sample4_Heat_Violations$AE,Sample4_Heat_Violations$Loss$Loss)
  Holdout_Sample4_Gas_Pearson[j,]   <- c(Sample4_Gas_Violations$LRuc[1],Sample4_Gas_Violations$LRcc[1],Sample4_Gas_Violations$DQ$stat,Sample4_Gas_Violations$AE,Sample4_Gas_Violations$Loss$Loss)
  Holdout_Sample4_Heat_Pearson_pvalues[j,] <- c(Sample4_Heat_Violations$LRuc[2],Sample4_Heat_Violations$LRcc[2],Sample4_Heat_Violations$DQ$pvalue)
  Holdout_Sample4_Gas_Pearson_pvalues[j,] <- c(Sample4_Gas_Violations$LRuc[2],Sample4_Gas_Violations$LRcc[2],Sample4_Gas_Violations$DQ$pvalue)
}

## Store the volatility and mu forecast for Sample 5

Holdout_Sample5_Heat_Pearson  <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample5_Gas_Pearson    <- matrix(data=0,nrow=length(signific_levels),ncol=5)
Holdout_Sample5_Heat_Pearson_pvalues  <- matrix(data=0,nrow=length(signific_levels),ncol=3)
Holdout_Sample5_Gas_Pearson_pvalues    <- matrix(data=0,nrow=length(signific_levels),ncol=3)
colnames(Holdout_Sample5_Heat_Pearson) <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample5_Gas_Pearson)   <- c("UC","CC","DQ","AE","QL")
colnames(Holdout_Sample5_Heat_Pearson_pvalues) <- c("UC","CC","DQ")
colnames(Holdout_Sample5_Gas_Pearson_pvalues)   <- c("UC","CC","DQ")
rownames(Holdout_Sample5_Heat_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample5_Gas_Pearson) <- c("1%","5%","10%")
rownames(Holdout_Sample5_Heat_Pearson_pvalues) <- c("1%","5%","10%")
rownames(Holdout_Sample5_Gas_Pearson_pvalues) <- c("1%","5%","10%")

Sample5_Heat_Holdout_Mu <- as.data.frame(Heat_Roll_Normal_Sample5)[,'Mu']
Sample5_Heat_Holdout_Sigma <- as.data.frame(Heat_Roll_Normal_Sample5)[,'Sigma']
Sample5_Gas_Holdout_Mu <- as.data.frame(Gas_Roll_Normal_Sample5)[,'Mu']
Sample5_Gas_Holdout_Sigma <- as.data.frame(Gas_Roll_Normal_Sample5)[,'Sigma']

for(j in 1:3)
{
  
  Z_Pearson_Sample5_Heat <- qpearson(signific_levels[j],Heat_Param_Pearson_Sample5)
  Z_Pearson_Sample5_Gas   <- qpearson(signific_levels[j],Gas_Param_Pearson_Sample5)
  
  Sample5_Heat_Holdout_VaR <- Sample5_Heat_Holdout_Mu + Sample5_Heat_Holdout_Sigma*Z_Pearson_Sample5_Heat
  Sample5_Gas_Holdout_VaR   <- Sample5_Gas_Holdout_Mu + Sample5_Gas_Holdout_Sigma*Z_Pearson_Sample5_Gas
  Sample5_Heat_Violations <- BacktestVaR(Returns_Heat[(length(Returns_Heat)-250+1):(length(Returns_Heat)-1)],Sample5_Heat_Holdout_VaR[-1],signific_levels[j])
  Sample5_Gas_Violations   <- BacktestVaR(Returns_Gas[(length(Returns_Gas)-250+1):(length(Returns_Gas)-1)],Sample5_Gas_Holdout_VaR[-1],signific_levels[j])
  Holdout_Sample5_Heat_Pearson[j,] <- c(Sample5_Heat_Violations$LRuc[1],Sample5_Heat_Violations$LRcc[1],Sample5_Heat_Violations$DQ$stat,Sample5_Heat_Violations$AE,Sample5_Heat_Violations$Loss$Loss)
  Holdout_Sample5_Gas_Pearson[j,]   <- c(Sample5_Gas_Violations$LRuc[1],Sample5_Gas_Violations$LRcc[1],Sample5_Gas_Violations$DQ$stat,Sample5_Gas_Violations$AE,Sample5_Gas_Violations$Loss$Loss)
  Holdout_Sample5_Heat_Pearson_pvalues[j,] <- c(Sample5_Heat_Violations$LRuc[2],Sample5_Heat_Violations$LRcc[2],Sample5_Heat_Violations$DQ$pvalue)
  Holdout_Sample5_Gas_Pearson_pvalues[j,] <- c(Sample5_Gas_Violations$LRuc[2],Sample5_Gas_Violations$LRcc[2],Sample5_Gas_Violations$DQ$pvalue)
}

### Plotting the VaR for subsamples

Heat_Sample1_Dates <- format(Heat$Date[(Heat_samples[1]-249):Heat_samples[1]],"%b-%y")
Heat_Sample2_Dates <- format(Heat$Date[(Heat_samples[2]-249):Heat_samples[2]],"%b-%y")
Heat_Sample3_Dates <- format(Heat$Date[(Heat_samples[3]-249):Heat_samples[3]],"%b-%y")
Heat_Sample4_Dates <- format(Heat$Date[(Heat_samples[4]-249):Heat_samples[4]],"%b-%y")
Heat_Sample5_Dates <- format(Heat$Date[(length(Heat_xts)-249):length(Heat_xts)],"%b-%y")

Gas_Sample1_Dates <- format(Gas$Date[(Gas_samples[1]-249):Gas_samples[1]],"%b-%y")
Gas_Sample2_Dates <- format(Gas$Date[(Gas_samples[2]-249):Gas_samples[2]],"%b-%y")
Gas_Sample3_Dates <- format(Gas$Date[(Gas_samples[3]-249):Gas_samples[3]],"%b-%y")
Gas_Sample4_Dates <- format(Gas$Date[(Gas_samples[4]-249):Gas_samples[4]],"%b-%y")
Gas_Sample5_Dates <- format(Gas$Date[(length(Gas_xts)-249):length(Gas_xts)],"%b-%y")

par(mfrow=c(2,3))
plot(Gas_Roll_Normal_Sample1@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample1",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Gas_Roll_Normal_Sample1@forecast$VaR[,3],type="l",col="yellow")
lines(Gas_Roll_SST_Sample1@forecast$VaR[,3],type="l",col="green")
lines(Gas_Roll_SGED_Sample1@forecast$VaR[,3],type="l",col="gray")
lines(Gas_Roll_JSU_Sample1@forecast$VaR[,3],type="l",col="brown")
lines(Sample1_Gas_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Gas_Sample1_Dates[c(1,50,100,150,200,250)],las=2)


plot(Gas_Roll_JSU_Sample2@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample2",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Gas_Roll_Normal_Sample2@forecast$VaR[,3],type="l",col="yellow")
lines(Gas_Roll_SST_Sample2@forecast$VaR[,3],type="l",col="green")
lines(Gas_Roll_SGED_Sample2@forecast$VaR[,3],type="l",col="gray")
lines(Gas_Roll_JSU_Sample2@forecast$VaR[,3],type="l",col="brown")
lines(Sample2_Gas_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Gas_Sample2_Dates[c(1,50,100,150,200,250)],las=2)

plot(Gas_Roll_JSU_Sample3@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample3",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Gas_Roll_Normal_Sample3@forecast$VaR[,3],type="l",col="yellow")
lines(Gas_Roll_SST_Sample3@forecast$VaR[,3],type="l",col="green")
lines(Gas_Roll_SGED_Sample3@forecast$VaR[,3],type="l",col="gray")
lines(Gas_Roll_JSU_Sample3@forecast$VaR[,3],type="l",col="brown")
lines(Sample3_Gas_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Gas_Sample3_Dates[c(1,50,100,150,200,250)],las=2)

plot(Gas_Roll_JSU_Sample4@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample4",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Gas_Roll_Normal_Sample4@forecast$VaR[,3],type="l",col="yellow")
lines(Gas_Roll_SST_Sample4@forecast$VaR[,3],type="l",col="green")
lines(Gas_Roll_SGED_Sample4@forecast$VaR[,3],type="l",col="gray")
lines(Gas_Roll_JSU_Sample4@forecast$VaR[,3],type="l",col="brown")
lines(Sample4_Gas_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Gas_Sample4_Dates[c(1,50,100,150,200,250)],las=2)

plot(Gas_Roll_JSU_Sample5@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample5",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Gas_Roll_Normal_Sample5@forecast$VaR[,3],type="l",col="yellow")
lines(Gas_Roll_SST_Sample5@forecast$VaR[,3],type="l",col="green")
lines(Gas_Roll_SGED_Sample5@forecast$VaR[,3],type="l",col="gray")
lines(Gas_Roll_JSU_Sample5@forecast$VaR[,3],type="l",col="brown")
lines(Sample5_Gas_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Gas_Sample5_Dates[c(1,50,100,150,200,250)],las=2)

### For Heat

par(mfrow=c(2,3))
plot(Heat_Roll_Normal_Sample1@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample1",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Heat_Roll_Normal_Sample1@forecast$VaR[,3],type="l",col="yellow")
lines(Heat_Roll_SST_Sample1@forecast$VaR[,3],type="l",col="green")
lines(Heat_Roll_SGED_Sample1@forecast$VaR[,3],type="l",col="gray")
lines(Heat_Roll_JSU_Sample1@forecast$VaR[,3],type="l",col="brown")
lines(Sample1_Heat_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Heat_Sample1_Dates[c(1,50,100,150,200,250)],las=2)


plot(Heat_Roll_JSU_Sample2@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample2",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Heat_Roll_Normal_Sample2@forecast$VaR[,3],type="l",col="yellow")
lines(Heat_Roll_SST_Sample2@forecast$VaR[,3],type="l",col="green")
lines(Heat_Roll_SGED_Sample2@forecast$VaR[,3],type="l",col="gray")
lines(Heat_Roll_JSU_Sample2@forecast$VaR[,3],type="l",col="brown")
lines(Sample2_Heat_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Heat_Sample2_Dates[c(1,50,100,150,200,250)],las=2)

plot(Heat_Roll_JSU_Sample3@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample3",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Heat_Roll_Normal_Sample3@forecast$VaR[,3],type="l",col="yellow")
lines(Heat_Roll_SST_Sample3@forecast$VaR[,3],type="l",col="green")
lines(Heat_Roll_SGED_Sample3@forecast$VaR[,3],type="l",col="gray")
lines(Heat_Roll_JSU_Sample3@forecast$VaR[,3],type="l",col="brown")
lines(Sample3_Heat_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Heat_Sample3_Dates[c(1,50,100,150,200,250)],las=2)

plot(Heat_Roll_JSU_Sample4@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample4",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Heat_Roll_Normal_Sample4@forecast$VaR[,3],type="l",col="yellow")
lines(Heat_Roll_SST_Sample4@forecast$VaR[,3],type="l",col="green")
lines(Heat_Roll_SGED_Sample4@forecast$VaR[,3],type="l",col="gray")
lines(Heat_Roll_JSU_Sample4@forecast$VaR[,3],type="l",col="brown")
lines(Sample4_Heat_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Heat_Sample4_Dates[c(1,50,100,150,200,250)],las=2)

plot(Heat_Roll_JSU_Sample5@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample5",ylab="Returns/ VaR at 10%",xlab="Date",xaxt="n")
lines(Heat_Roll_Normal_Sample5@forecast$VaR[,3],type="l",col="yellow")
lines(Heat_Roll_SST_Sample5@forecast$VaR[,3],type="l",col="green")
lines(Heat_Roll_SGED_Sample5@forecast$VaR[,3],type="l",col="gray")
lines(Heat_Roll_JSU_Sample5@forecast$VaR[,3],type="l",col="brown")
lines(Sample5_Heat_Holdout_VaR[-1],type="l",col="black")
axis(1,at=c(0,50,100,150,200,250),labels=Heat_Sample5_Dates[c(1,50,100,150,200,250)],las=2)


### Incorporating Structural Breaks in a different way using dummies###
# 

par(mfrow=c(1,2))
plot(Heat_Roll8@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample for Heat",ylab="Returns/ VaR at 1%",xlab="",xaxt="n")
lines(Heat_Roll8@forecast$VaR[,1],type="l",col="yellow")
lines(Heat_Roll5@forecast$VaR[,1],type="l",col="green")
lines(Heat_Roll6@forecast$VaR[,1],type="l",col="gray")
lines(Heat_Roll7@forecast$VaR[,1],type="l",col="brown")
lines(Heat_Holdout_VaR[-1],type="l",col="black")
legend(x=1,y=0.14,legend=c("Returns","Normal","SST","SGED","JSU","Pearson"),col=c("blue","yellow","green","gray","brown","black"),cex=0.8,lty=1,bty="n",y.intersp = 0.7,ncol=2,x.intersp = 0.49)
axis(1,at=c(0,200,400,600,800,950),labels=format(Heat_Date[c(4100,4300,4500,4700,4900,5100)],"%b-%y"),las=2)

### Keep the plot width as 800 pixels and save as .tiff extension

plot(Gas_Roll8@forecast$VaR[,4],type="l",col="blue",main="Forecasting Sample for Gas",ylab="Returns/ VaR at 1%",xlab="",xaxt="n")
lines(Gas_Roll8@forecast$VaR[,1],type="l",col="yellow")
lines(Gas_Roll5@forecast$VaR[,1],type="l",col="green")
lines(Gas_Roll6@forecast$VaR[,1],type="l",col="gray")
lines(Gas_Roll7@forecast$VaR[,1],type="l",col="brown")
lines(Gas_Holdout_VaR[-1],type="l",col="black")
legend(x=300,y=0.61,legend=c("Returns","Normal","SST","SGED","JSU","Pearson"),col=c("blue","yellow","green","gray","brown","black"),cex=0.8,lty=1,bty="n",y.intersp = 0.7,x.intersp = 0.49,ncol=2)
axis(1,at=c(0,200,400,600,800,950),labels=format(Gas_Date[c(4001,4200,4400,4600,4800,5000)],"%b-%y"),las=2)

## Storing Pearson subsample results ###

Heat_subsample_1pc <- rbind(Holdout_Sample1_Heat_Pearson[1,],Holdout_Sample2_Heat_Pearson[1,],Holdout_Sample3_Heat_Pearson[1,],Holdout_Sample4_Heat_Pearson[1,],Holdout_Sample5_Heat_Pearson[1,])
Gas_subsample_1pc   <- rbind(Holdout_Sample1_Gas_Pearson[1,],Holdout_Sample2_Gas_Pearson[1,],Holdout_Sample3_Gas_Pearson[1,],Holdout_Sample4_Gas_Pearson[1,],Holdout_Sample5_Gas_Pearson[1,])

Heat_subsample_5pc <- rbind(Holdout_Sample1_Heat_Pearson[2,],Holdout_Sample2_Heat_Pearson[2,],Holdout_Sample3_Heat_Pearson[2,],Holdout_Sample4_Heat_Pearson[2,],Holdout_Sample5_Heat_Pearson[2,])
Gas_subsample_5pc   <- rbind(Holdout_Sample1_Gas_Pearson[2,],Holdout_Sample2_Gas_Pearson[2,],Holdout_Sample3_Gas_Pearson[2,],Holdout_Sample4_Gas_Pearson[2,],Holdout_Sample5_Gas_Pearson[2,])

Heat_subsample_10pc <- rbind(Holdout_Sample1_Heat_Pearson[3,],Holdout_Sample2_Heat_Pearson[3,],Holdout_Sample3_Heat_Pearson[3,],Holdout_Sample4_Heat_Pearson[3,],Holdout_Sample5_Heat_Pearson[3,])
Gas_subsample_10pc   <- rbind(Holdout_Sample1_Gas_Pearson[3,],Holdout_Sample2_Gas_Pearson[3,],Holdout_Sample3_Gas_Pearson[3,],Holdout_Sample4_Gas_Pearson[3,],Holdout_Sample5_Gas_Pearson[3,])

Heat_subsample_pvalues_1pc <- rbind(Holdout_Sample1_Heat_Pearson_pvalues[1,],Holdout_Sample2_Heat_Pearson_pvalues[1,],Holdout_Sample3_Heat_Pearson_pvalues[1,],Holdout_Sample4_Heat_Pearson_pvalues[1,],Holdout_Sample5_Heat_Pearson_pvalues[1,])
Gas_subsample_pvalues_1pc   <- rbind(Holdout_Sample1_Gas_Pearson_pvalues[1,],Holdout_Sample2_Gas_Pearson_pvalues[1,],Holdout_Sample3_Gas_Pearson_pvalues[1,],Holdout_Sample4_Gas_Pearson_pvalues[1,],Holdout_Sample5_Gas_Pearson_pvalues[1,])

Heat_subsample_pvalues_5pc <- rbind(Holdout_Sample1_Heat_Pearson_pvalues[2,],Holdout_Sample2_Heat_Pearson_pvalues[2,],Holdout_Sample3_Heat_Pearson_pvalues[2,],Holdout_Sample4_Heat_Pearson_pvalues[2,],Holdout_Sample5_Heat_Pearson_pvalues[2,])
Gas_subsample_pvalues_5pc   <- rbind(Holdout_Sample1_Gas_Pearson_pvalues[2,],Holdout_Sample2_Gas_Pearson_pvalues[2,],Holdout_Sample3_Gas_Pearson_pvalues[2,],Holdout_Sample4_Gas_Pearson_pvalues[2,],Holdout_Sample5_Gas_Pearson_pvalues[2,])

Heat_subsample_pvalues_10pc <- rbind(Holdout_Sample1_Heat_Pearson_pvalues[3,],Holdout_Sample2_Heat_Pearson_pvalues[3,],Holdout_Sample3_Heat_Pearson_pvalues[3,],Holdout_Sample4_Heat_Pearson_pvalues[3,],Holdout_Sample5_Heat_Pearson_pvalues[3,])
Gas_subsample_pvalues_10pc   <- rbind(Holdout_Sample1_Gas_Pearson_pvalues[3,],Holdout_Sample2_Gas_Pearson_pvalues[3,],Holdout_Sample3_Gas_Pearson_pvalues[3,],Holdout_Sample4_Gas_Pearson_pvalues[3,],Holdout_Sample5_Gas_Pearson_pvalues[3,])

#### VaRDuration Test for subsamples

Heat_Dur_Normal_Sample1 <- VaRDurTest(0.10,Heat_Roll_Normal_Sample1@forecast$VaR[,4],Heat_Roll_Normal_Sample1@forecast$VaR[,3])
Heat_Dur_Normal_Sample2 <- VaRDurTest(0.10,Heat_Roll_Normal_Sample2@forecast$VaR[,4],Heat_Roll_Normal_Sample2@forecast$VaR[,3])
Heat_Dur_Normal_Sample3 <- VaRDurTest(0.10,Heat_Roll_Normal_Sample3@forecast$VaR[,4],Heat_Roll_Normal_Sample3@forecast$VaR[,3])
Heat_Dur_Normal_Sample4 <- VaRDurTest(0.10,Heat_Roll_Normal_Sample4@forecast$VaR[,4],Heat_Roll_Normal_Sample4@forecast$VaR[,3])
Heat_Dur_Normal_Sample5 <- VaRDurTest(0.10,Heat_Roll_Normal_Sample5@forecast$VaR[,4],Heat_Roll_Normal_Sample5@forecast$VaR[,3])

Gas_Dur_Normal_Sample1 <- VaRDurTest(0.10,Gas_Roll_Normal_Sample1@forecast$VaR[,4],Gas_Roll_Normal_Sample1@forecast$VaR[,3])
Gas_Dur_Normal_Sample2 <- VaRDurTest(0.10,Gas_Roll_Normal_Sample2@forecast$VaR[,4],Gas_Roll_Normal_Sample2@forecast$VaR[,3])
Gas_Dur_Normal_Sample3 <- VaRDurTest(0.10,Gas_Roll_Normal_Sample3@forecast$VaR[,4],Gas_Roll_Normal_Sample3@forecast$VaR[,3])
Gas_Dur_Normal_Sample4 <- VaRDurTest(0.10,Gas_Roll_Normal_Sample4@forecast$VaR[,4],Gas_Roll_Normal_Sample4@forecast$VaR[,3])
Gas_Dur_Normal_Sample5 <- VaRDurTest(0.10,Gas_Roll_Normal_Sample5@forecast$VaR[,4],Gas_Roll_Normal_Sample5@forecast$VaR[,3])

Heat_Dur_SST_Sample1 <- VaRDurTest(0.10,Heat_Roll_SST_Sample1@forecast$VaR[,4],Heat_Roll_SST_Sample1@forecast$VaR[,3])
Heat_Dur_SST_Sample2 <- VaRDurTest(0.10,Heat_Roll_SST_Sample2@forecast$VaR[,4],Heat_Roll_SST_Sample2@forecast$VaR[,3])
Heat_Dur_SST_Sample3 <- VaRDurTest(0.10,Heat_Roll_SST_Sample3@forecast$VaR[,4],Heat_Roll_SST_Sample3@forecast$VaR[,3])
Heat_Dur_SST_Sample4 <- VaRDurTest(0.10,Heat_Roll_SST_Sample4@forecast$VaR[,4],Heat_Roll_SST_Sample4@forecast$VaR[,3])
Heat_Dur_SST_Sample5 <- VaRDurTest(0.10,Heat_Roll_SST_Sample5@forecast$VaR[,4],Heat_Roll_SST_Sample5@forecast$VaR[,3])

Gas_Dur_SST_Sample1 <- VaRDurTest(0.10,Gas_Roll_SST_Sample1@forecast$VaR[,4],Gas_Roll_SST_Sample1@forecast$VaR[,3])
Gas_Dur_SST_Sample2 <- VaRDurTest(0.10,Gas_Roll_SST_Sample2@forecast$VaR[,4],Gas_Roll_SST_Sample2@forecast$VaR[,3])
Gas_Dur_SST_Sample3 <- VaRDurTest(0.10,Gas_Roll_SST_Sample3@forecast$VaR[,4],Gas_Roll_SST_Sample3@forecast$VaR[,3])
Gas_Dur_SST_Sample4 <- VaRDurTest(0.10,Gas_Roll_SST_Sample4@forecast$VaR[,4],Gas_Roll_SST_Sample4@forecast$VaR[,3])
Gas_Dur_SST_Sample5 <- VaRDurTest(0.10,Gas_Roll_SST_Sample5@forecast$VaR[,4],Gas_Roll_SST_Sample5@forecast$VaR[,3])

Heat_Dur_SGED_Sample1 <- VaRDurTest(0.10,Heat_Roll_SGED_Sample1@forecast$VaR[,4],Heat_Roll_SGED_Sample1@forecast$VaR[,3])
Heat_Dur_SGED_Sample2 <- VaRDurTest(0.10,Heat_Roll_SGED_Sample2@forecast$VaR[,4],Heat_Roll_SGED_Sample2@forecast$VaR[,3])
Heat_Dur_SGED_Sample3 <- VaRDurTest(0.10,Heat_Roll_SGED_Sample3@forecast$VaR[,4],Heat_Roll_SGED_Sample3@forecast$VaR[,3])
Heat_Dur_SGED_Sample4 <- VaRDurTest(0.10,Heat_Roll_SGED_Sample4@forecast$VaR[,4],Heat_Roll_SGED_Sample4@forecast$VaR[,3])
Heat_Dur_SGED_Sample5 <- VaRDurTest(0.10,Heat_Roll_SGED_Sample5@forecast$VaR[,4],Heat_Roll_SGED_Sample5@forecast$VaR[,3])

Gas_Dur_SGED_Sample1 <- VaRDurTest(0.10,Gas_Roll_SGED_Sample1@forecast$VaR[,4],Gas_Roll_SGED_Sample1@forecast$VaR[,3])
Gas_Dur_SGED_Sample2 <- VaRDurTest(0.10,Gas_Roll_SGED_Sample2@forecast$VaR[,4],Gas_Roll_SGED_Sample2@forecast$VaR[,3])
Gas_Dur_SGED_Sample3 <- VaRDurTest(0.10,Gas_Roll_SGED_Sample3@forecast$VaR[,4],Gas_Roll_SGED_Sample3@forecast$VaR[,3])
Gas_Dur_SGED_Sample4 <- VaRDurTest(0.10,Gas_Roll_SGED_Sample4@forecast$VaR[,4],Gas_Roll_SGED_Sample4@forecast$VaR[,3])
Gas_Dur_SGED_Sample5 <- VaRDurTest(0.10,Gas_Roll_SGED_Sample5@forecast$VaR[,4],Gas_Roll_SGED_Sample5@forecast$VaR[,3])

Heat_Dur_JSU_Sample1 <- VaRDurTest(0.10,Heat_Roll_JSU_Sample1@forecast$VaR[,4],Heat_Roll_JSU_Sample1@forecast$VaR[,3])
Heat_Dur_JSU_Sample2 <- VaRDurTest(0.10,Heat_Roll_JSU_Sample2@forecast$VaR[,4],Heat_Roll_JSU_Sample2@forecast$VaR[,3])
Heat_Dur_JSU_Sample3 <- VaRDurTest(0.10,Heat_Roll_JSU_Sample3@forecast$VaR[,4],Heat_Roll_JSU_Sample3@forecast$VaR[,3])
Heat_Dur_JSU_Sample4 <- VaRDurTest(0.10,Heat_Roll_JSU_Sample4@forecast$VaR[,4],Heat_Roll_JSU_Sample4@forecast$VaR[,3])
Heat_Dur_JSU_Sample5 <- VaRDurTest(0.10,Heat_Roll_JSU_Sample5@forecast$VaR[,4],Heat_Roll_JSU_Sample5@forecast$VaR[,3])

Gas_Dur_JSU_Sample1 <- VaRDurTest(0.10,Gas_Roll_JSU_Sample1@forecast$VaR[,4],Gas_Roll_JSU_Sample1@forecast$VaR[,3])
Gas_Dur_JSU_Sample2 <- VaRDurTest(0.10,Gas_Roll_JSU_Sample2@forecast$VaR[,4],Gas_Roll_JSU_Sample2@forecast$VaR[,3])
Gas_Dur_JSU_Sample3 <- VaRDurTest(0.10,Gas_Roll_JSU_Sample3@forecast$VaR[,4],Gas_Roll_JSU_Sample3@forecast$VaR[,3])
Gas_Dur_JSU_Sample4 <- VaRDurTest(0.10,Gas_Roll_JSU_Sample4@forecast$VaR[,4],Gas_Roll_JSU_Sample4@forecast$VaR[,3])
Gas_Dur_JSU_Sample5 <- VaRDurTest(0.10,Gas_Roll_JSU_Sample5@forecast$VaR[,4],Gas_Roll_JSU_Sample5@forecast$VaR[,3])

## Store the results

Heat_Dur_Subsample_stat <- matrix(c(Heat_Dur_Normal_Sample1$rLL,
                                 Heat_Dur_SST_Sample1$rLL,
                                 Heat_Dur_SGED_Sample1$rLL,
                                 Heat_Dur_JSU_Sample1$rLL,
                                 Heat_Dur_Normal_Sample2$rLL,
                                 Heat_Dur_SST_Sample2$rLL,
                                 Heat_Dur_SGED_Sample2$rLL,
                                 Heat_Dur_JSU_Sample2$rLL,
                                 Heat_Dur_Normal_Sample3$rLL,
                                 Heat_Dur_SST_Sample3$rLL,
                                 Heat_Dur_SGED_Sample3$rLL,
                                 Heat_Dur_JSU_Sample3$rLL,
                                 Heat_Dur_Normal_Sample4$rLL,
                                 Heat_Dur_SST_Sample4$rLL,
                                 Heat_Dur_SGED_Sample4$rLL,
                                 Heat_Dur_JSU_Sample4$rLL,
                                 Heat_Dur_Normal_Sample5$rLL,
                                 Heat_Dur_SST_Sample5$rLL,
                                 Heat_Dur_SGED_Sample5$rLL,
                                 Heat_Dur_JSU_Sample5$rLL),nrow=4)

Gas_Dur_Subsample_stat <- matrix(c(Gas_Dur_Normal_Sample1$rLL,
                                 Gas_Dur_SST_Sample1$rLL,
                                 Gas_Dur_SGED_Sample1$rLL,
                                 Gas_Dur_JSU_Sample1$rLL,
                                 Gas_Dur_Normal_Sample2$rLL,
                                 Gas_Dur_SST_Sample2$rLL,
                                 Gas_Dur_SGED_Sample2$rLL,
                                 Gas_Dur_JSU_Sample2$rLL,
                                 Gas_Dur_Normal_Sample3$rLL,
                                 Gas_Dur_SST_Sample3$rLL,
                                 Gas_Dur_SGED_Sample3$rLL,
                                 Gas_Dur_JSU_Sample3$rLL,
                                 Gas_Dur_Normal_Sample4$rLL,
                                 Gas_Dur_SST_Sample4$rLL,
                                 Gas_Dur_SGED_Sample4$rLL,
                                 Gas_Dur_JSU_Sample4$rLL,
                                 Gas_Dur_Normal_Sample5$rLL,
                                 Gas_Dur_SST_Sample5$rLL,
                                 Gas_Dur_SGED_Sample5$rLL,
                                 Gas_Dur_JSU_Sample5$rLL),nrow=4)

Heat_Dur_Subsample_pvalues <- matrix(c(Heat_Dur_Normal_Sample1$LRp,
                                 Heat_Dur_SST_Sample1$LRp,
                                 Heat_Dur_SGED_Sample1$LRp,
                                 Heat_Dur_JSU_Sample1$LRp,
                                 Heat_Dur_Normal_Sample2$LRp,
                                 Heat_Dur_SST_Sample2$LRp,
                                 Heat_Dur_SGED_Sample2$LRp,
                                 Heat_Dur_JSU_Sample2$LRp,
                                 Heat_Dur_Normal_Sample3$LRp,
                                 Heat_Dur_SST_Sample3$LRp,
                                 Heat_Dur_SGED_Sample3$LRp,
                                 Heat_Dur_JSU_Sample3$LRp,
                                 Heat_Dur_Normal_Sample4$LRp,
                                 Heat_Dur_SST_Sample4$LRp,
                                 Heat_Dur_SGED_Sample4$LRp,
                                 Heat_Dur_JSU_Sample4$LRp,
                                 Heat_Dur_Normal_Sample5$LRp,
                                 Heat_Dur_SST_Sample5$LRp,
                                 Heat_Dur_SGED_Sample5$LRp,
                                 Heat_Dur_JSU_Sample5$LRp),nrow=4)

Gas_Dur_Subsample_pvalues <- matrix(c(Gas_Dur_Normal_Sample1$LRp,
                               Gas_Dur_SST_Sample1$LRp,
                               Gas_Dur_SGED_Sample1$LRp,
                               Gas_Dur_JSU_Sample1$LRp,
                               Gas_Dur_Normal_Sample2$LRp,
                               Gas_Dur_SST_Sample2$LRp,
                               Gas_Dur_SGED_Sample2$LRp,
                               Gas_Dur_JSU_Sample2$LRp,
                               Gas_Dur_Normal_Sample3$LRp,
                               Gas_Dur_SST_Sample3$LRp,
                               Gas_Dur_SGED_Sample3$LRp,
                               Gas_Dur_JSU_Sample3$LRp,
                               Gas_Dur_Normal_Sample4$LRp,
                               Gas_Dur_SST_Sample4$LRp,
                               Gas_Dur_SGED_Sample4$LRp,
                               Gas_Dur_JSU_Sample4$LRp,
                               Gas_Dur_Normal_Sample5$LRp,
                               Gas_Dur_SST_Sample5$LRp,
                               Gas_Dur_SGED_Sample5$LRp,
                               Gas_Dur_JSU_Sample5$LRp),nrow=4)


### Expected Shortfall Calculations ####

## Assuming Normal Distribution

Heat_index_normal1 <- which(Heat_Roll8@forecast$VaR[,1]>Heat_Roll8@forecast$VaR[,4])
Heat_VaR_values_normal1 <- Heat_Roll8@forecast$VaR[Heat_index_normal1,1]
Heat_return_values_normal1 <- Heat_Roll8@forecast$VaR[Heat_index_normal1,4]
Heat_ES_normal1 <- (Heat_return_values_normal1/Heat_VaR_values_normal1)
Heat_ES_Measure1_normal1 <- mean(Heat_ES_normal1)
Heat_ES_Measure2_normal1 <- max(Heat_ES_normal1)

Heat_index_normal5 <- which(Heat_Roll8@forecast$VaR[,2]>Heat_Roll8@forecast$VaR[,4])
Heat_VaR_values_normal5 <- Heat_Roll8@forecast$VaR[Heat_index_normal5,2]
Heat_return_values_normal5 <- Heat_Roll8@forecast$VaR[Heat_index_normal5,4]
Heat_ES_normal5 <- (Heat_return_values_normal5/Heat_VaR_values_normal5)
Heat_ES_Measure1_normal5 <- mean(Heat_ES_normal5)
Heat_ES_Measure2_normal5 <- max(Heat_ES_normal5)

Heat_index_normal10 <- which(Heat_Roll8@forecast$VaR[,3]>Heat_Roll8@forecast$VaR[,4])
Heat_VaR_values_normal10 <- Heat_Roll8@forecast$VaR[Heat_index_normal10,3]
Heat_return_values_normal10 <- Heat_Roll8@forecast$VaR[Heat_index_normal10,4]
Heat_ES_normal10 <- (Heat_return_values_normal10/Heat_VaR_values_normal10)
Heat_ES_Measure1_normal10 <- mean(Heat_ES_normal10)
Heat_ES_Measure2_normal10 <- max(Heat_ES_normal10)

Gas_index_normal1 <- which(Gas_Roll8@forecast$VaR[,1]>Gas_Roll8@forecast$VaR[,4])
Gas_VaR_values_normal1 <- Gas_Roll8@forecast$VaR[Gas_index_normal1,1]
Gas_return_values_normal1 <- Gas_Roll8@forecast$VaR[Gas_index_normal1,4]
Gas_ES_normal1 <- (Gas_return_values_normal1/Gas_VaR_values_normal1)
Gas_ES_Measure1_normal1 <- mean(Gas_ES_normal1)
Gas_ES_Measure2_normal1 <- max(Gas_ES_normal1)

Gas_index_normal5 <- which(Gas_Roll8@forecast$VaR[,2]>Gas_Roll8@forecast$VaR[,4])
Gas_VaR_values_normal5 <- Gas_Roll8@forecast$VaR[Gas_index_normal5,2]
Gas_return_values_normal5 <- Gas_Roll8@forecast$VaR[Gas_index_normal5,4]
Gas_ES_normal5 <- (Gas_return_values_normal5/Gas_VaR_values_normal5)
Gas_ES_Measure1_normal5 <- mean(Gas_ES_normal5)
Gas_ES_Measure2_normal5 <- max(Gas_ES_normal5)

Gas_index_normal10 <- which(Gas_Roll8@forecast$VaR[,3]>Gas_Roll8@forecast$VaR[,4])
Gas_VaR_values_normal10 <- Gas_Roll8@forecast$VaR[Gas_index_normal10,3]
Gas_return_values_normal10 <- Gas_Roll8@forecast$VaR[Gas_index_normal10,4]
Gas_ES_normal10 <- (Gas_return_values_normal10/Gas_VaR_values_normal10)
Gas_ES_Measure1_normal10 <- mean(Gas_ES_normal10)
Gas_ES_Measure2_normal10 <- max(Gas_ES_normal10)

## Assuming SST Distribution

Heat_index_SST1 <- which(Heat_Roll5@forecast$VaR[,1]>Heat_Roll5@forecast$VaR[,4])
Heat_VaR_values_SST1 <- Heat_Roll5@forecast$VaR[Heat_index_SST1,1]
Heat_return_values_SST1 <- Heat_Roll5@forecast$VaR[Heat_index_SST1,4]
Heat_ES_SST1 <- (Heat_return_values_SST1/Heat_VaR_values_SST1)
Heat_ES_Measure1_SST1 <- mean(Heat_ES_SST1)
Heat_ES_Measure2_SST1 <- max(Heat_ES_SST1)

Heat_index_SST5 <- which(Heat_Roll5@forecast$VaR[,2]>Heat_Roll5@forecast$VaR[,4])
Heat_VaR_values_SST5 <- Heat_Roll5@forecast$VaR[Heat_index_SST5,2]
Heat_return_values_SST5 <- Heat_Roll5@forecast$VaR[Heat_index_SST5,4]
Heat_ES_SST5 <- (Heat_return_values_SST5/Heat_VaR_values_SST5)
Heat_ES_Measure1_SST5 <- mean(Heat_ES_SST5)
Heat_ES_Measure2_SST5 <- max(Heat_ES_SST5)

Heat_index_SST10 <- which(Heat_Roll5@forecast$VaR[,3]>Heat_Roll5@forecast$VaR[,4])
Heat_VaR_values_SST10 <- Heat_Roll5@forecast$VaR[Heat_index_SST10,3]
Heat_return_values_SST10 <- Heat_Roll5@forecast$VaR[Heat_index_SST10,4]
Heat_ES_SST10 <- (Heat_return_values_SST10/Heat_VaR_values_SST10)
Heat_ES_Measure1_SST10 <- mean(Heat_ES_SST10)
Heat_ES_Measure2_SST10 <- max(Heat_ES_SST10)

Gas_index_SST1 <- which(Gas_Roll5@forecast$VaR[,1]>Gas_Roll5@forecast$VaR[,4])
Gas_VaR_values_SST1 <- Gas_Roll5@forecast$VaR[Gas_index_SST1,1]
Gas_return_values_SST1 <- Gas_Roll5@forecast$VaR[Gas_index_SST1,4]
Gas_ES_SST1 <- (Gas_return_values_SST1/Gas_VaR_values_SST1)
Gas_ES_Measure1_SST1 <- mean(Gas_ES_SST1)
Gas_ES_Measure2_SST1 <- max(Gas_ES_SST1)

Gas_index_SST5 <- which(Gas_Roll5@forecast$VaR[,2]>Gas_Roll5@forecast$VaR[,4])
Gas_VaR_values_SST5 <- Gas_Roll5@forecast$VaR[Gas_index_SST5,2]
Gas_return_values_SST5 <- Gas_Roll5@forecast$VaR[Gas_index_SST5,4]
Gas_ES_SST5 <- (Gas_return_values_SST5/Gas_VaR_values_SST5)
Gas_ES_Measure1_SST5 <- mean(Gas_ES_SST5)
Gas_ES_Measure2_SST5 <- max(Gas_ES_SST5)

Gas_index_SST10 <- which(Gas_Roll5@forecast$VaR[,3]>Gas_Roll5@forecast$VaR[,4])
Gas_VaR_values_SST10 <- Gas_Roll5@forecast$VaR[Gas_index_SST10,3]
Gas_return_values_SST10 <- Gas_Roll5@forecast$VaR[Gas_index_SST10,4]
Gas_ES_SST10 <- (Gas_return_values_SST10/Gas_VaR_values_SST10)
Gas_ES_Measure1_SST10 <- mean(Gas_ES_SST10)
Gas_ES_Measure2_SST10 <- max(Gas_ES_SST10)


## Assuming SGED Distribution

Heat_index_SGED1 <- which(Heat_Roll6@forecast$VaR[,1]>Heat_Roll6@forecast$VaR[,4])
Heat_VaR_values_SGED1 <- Heat_Roll6@forecast$VaR[Heat_index_SGED1,1]
Heat_return_values_SGED1 <- Heat_Roll6@forecast$VaR[Heat_index_SGED1,4]
Heat_ES_SGED1 <- (Heat_return_values_SGED1/Heat_VaR_values_SGED1)
Heat_ES_Measure1_SGED1 <- mean(Heat_ES_SGED1)
Heat_ES_Measure2_SGED1 <- max(Heat_ES_SGED1)

Heat_index_SGED5 <- which(Heat_Roll6@forecast$VaR[,2]>Heat_Roll6@forecast$VaR[,4])
Heat_VaR_values_SGED5 <- Heat_Roll6@forecast$VaR[Heat_index_SGED5,2]
Heat_return_values_SGED5 <- Heat_Roll6@forecast$VaR[Heat_index_SGED5,4]
Heat_ES_SGED5 <- (Heat_return_values_SGED5/Heat_VaR_values_SGED5)
Heat_ES_Measure1_SGED5 <- mean(Heat_ES_SGED5)
Heat_ES_Measure2_SGED5 <- max(Heat_ES_SGED5)

Heat_index_SGED10 <- which(Heat_Roll6@forecast$VaR[,3]>Heat_Roll6@forecast$VaR[,4])
Heat_VaR_values_SGED10 <- Heat_Roll6@forecast$VaR[Heat_index_SGED10,3]
Heat_return_values_SGED10 <- Heat_Roll6@forecast$VaR[Heat_index_SGED10,4]
Heat_ES_SGED10 <- (Heat_return_values_SGED10/Heat_VaR_values_SGED10)
Heat_ES_Measure1_SGED10 <- mean(Heat_ES_SGED10)
Heat_ES_Measure2_SGED10 <- max(Heat_ES_SGED10)

Gas_index_SGED1 <- which(Gas_Roll6@forecast$VaR[,1]>Gas_Roll6@forecast$VaR[,4])
Gas_VaR_values_SGED1 <- Gas_Roll6@forecast$VaR[Gas_index_SGED1,1]
Gas_return_values_SGED1 <- Gas_Roll6@forecast$VaR[Gas_index_SGED1,4]
Gas_ES_SGED1 <- (Gas_return_values_SGED1/Gas_VaR_values_SGED1)
Gas_ES_Measure1_SGED1 <- mean(Gas_ES_SGED1)
Gas_ES_Measure2_SGED1 <- max(Gas_ES_SGED1)

Gas_index_SGED5 <- which(Gas_Roll6@forecast$VaR[,2]>Gas_Roll6@forecast$VaR[,4])
Gas_VaR_values_SGED5 <- Gas_Roll6@forecast$VaR[Gas_index_SGED5,2]
Gas_return_values_SGED5 <- Gas_Roll6@forecast$VaR[Gas_index_SGED5,4]
Gas_ES_SGED5 <- (Gas_return_values_SGED5/Gas_VaR_values_SGED5)
Gas_ES_Measure1_SGED5 <- mean(Gas_ES_SGED5)
Gas_ES_Measure2_SGED5 <- max(Gas_ES_SGED5)

Gas_index_SGED10 <- which(Gas_Roll6@forecast$VaR[,3]>Gas_Roll6@forecast$VaR[,4])
Gas_VaR_values_SGED10 <- Gas_Roll6@forecast$VaR[Gas_index_SGED10,3]
Gas_return_values_SGED10 <- Gas_Roll6@forecast$VaR[Gas_index_SGED10,4]
Gas_ES_SGED10 <- (Gas_return_values_SGED10/Gas_VaR_values_SGED10)
Gas_ES_Measure1_SGED10 <- mean(Gas_ES_SGED10)
Gas_ES_Measure2_SGED10 <- max(Gas_ES_SGED10)


## Assuming JSU Distribution

Heat_index_JSU1 <- which(Heat_Roll7@forecast$VaR[,1]>Heat_Roll7@forecast$VaR[,4])
Heat_VaR_values_JSU1 <- Heat_Roll7@forecast$VaR[Heat_index_JSU1,1]
Heat_return_values_JSU1 <- Heat_Roll7@forecast$VaR[Heat_index_JSU1,4]
Heat_ES_JSU1 <- (Heat_return_values_JSU1/Heat_VaR_values_JSU1)
Heat_ES_Measure1_JSU1 <- mean(Heat_ES_JSU1)
Heat_ES_Measure2_JSU1 <- max(Heat_ES_JSU1)

Heat_index_JSU5 <- which(Heat_Roll7@forecast$VaR[,2]>Heat_Roll7@forecast$VaR[,4])
Heat_VaR_values_JSU5 <- Heat_Roll7@forecast$VaR[Heat_index_JSU5,2]
Heat_return_values_JSU5 <- Heat_Roll7@forecast$VaR[Heat_index_JSU5,4]
Heat_ES_JSU5 <- (Heat_return_values_JSU5/Heat_VaR_values_JSU5)
Heat_ES_Measure1_JSU5 <- mean(Heat_ES_JSU5)
Heat_ES_Measure2_JSU5 <- max(Heat_ES_JSU5)

Heat_index_JSU10 <- which(Heat_Roll7@forecast$VaR[,3]>Heat_Roll7@forecast$VaR[,4])
Heat_VaR_values_JSU10 <- Heat_Roll7@forecast$VaR[Heat_index_JSU10,3]
Heat_return_values_JSU10 <- Heat_Roll7@forecast$VaR[Heat_index_JSU10,4]
Heat_ES_JSU10 <- (Heat_return_values_JSU10/Heat_VaR_values_JSU10)
Heat_ES_Measure1_JSU10 <- mean(Heat_ES_JSU10)
Heat_ES_Measure2_JSU10 <- max(Heat_ES_JSU10)

Gas_index_JSU1 <- which(Gas_Roll7@forecast$VaR[,1]>Gas_Roll7@forecast$VaR[,4])
Gas_VaR_values_JSU1 <- Gas_Roll7@forecast$VaR[Gas_index_JSU1,1]
Gas_return_values_JSU1 <- Gas_Roll7@forecast$VaR[Gas_index_JSU1,4]
Gas_ES_JSU1 <- (Gas_return_values_JSU1/Gas_VaR_values_JSU1)
Gas_ES_Measure1_JSU1 <- mean(Gas_ES_JSU1)
Gas_ES_Measure2_JSU1 <- max(Gas_ES_JSU1)

Gas_index_JSU5 <- which(Gas_Roll7@forecast$VaR[,2]>Gas_Roll7@forecast$VaR[,4])
Gas_VaR_values_JSU5 <- Gas_Roll7@forecast$VaR[Gas_index_JSU5,2]
Gas_return_values_JSU5 <- Gas_Roll7@forecast$VaR[Gas_index_JSU5,4]
Gas_ES_JSU5 <- (Gas_return_values_JSU5/Gas_VaR_values_JSU5)
Gas_ES_Measure1_JSU5 <- mean(Gas_ES_JSU5)
Gas_ES_Measure2_JSU5 <- max(Gas_ES_JSU5)

Gas_index_JSU10 <- which(Gas_Roll7@forecast$VaR[,3]>Gas_Roll7@forecast$VaR[,4])
Gas_VaR_values_JSU10 <- Gas_Roll7@forecast$VaR[Gas_index_JSU10,3]
Gas_return_values_JSU10 <- Gas_Roll7@forecast$VaR[Gas_index_JSU10,4]
Gas_ES_JSU10 <- (Gas_return_values_JSU10/Gas_VaR_values_JSU10)
Gas_ES_Measure1_JSU10 <- mean(Gas_ES_JSU10)
Gas_ES_Measure2_JSU10 <- max(Gas_ES_JSU10)

### For Pearson

Heat_index_Pearson1 <- which(Heat_Holdout_VaR[-1]>Heat_Roll7@forecast$VaR[1:999,4])
Heat_VaR_values_Pearson1 <- Heat_Holdout_VaR[Heat_index_Pearson1]
Heat_return_values_Pearson1 <- Heat_Roll7@forecast$VaR[Heat_index_Pearson1,4]
Heat_ES_Pearson1 <- (Heat_return_values_Pearson1/Heat_VaR_values_Pearson1)
Heat_ES_Measure1_Pearson1 <- mean(Heat_ES_Pearson1)
Heat_ES_Measure2_Pearson1 <- max(Heat_ES_Pearson1)

Heat_index_Pearson5 <- which(Heat_Holdout_VaR[-1]>Heat_Roll7@forecast$VaR[1:999,4])
Heat_VaR_values_Pearson5 <- Heat_Holdout_VaR[Heat_index_Pearson5]
Heat_return_values_Pearson5 <- Heat_Roll7@forecast$VaR[Heat_index_Pearson5,4]
Heat_ES_Pearson5 <- (Heat_return_values_Pearson5/Heat_VaR_values_Pearson5)
Heat_ES_Measure1_Pearson5 <- mean(Heat_ES_Pearson5)
Heat_ES_Measure2_Pearson5 <- max(Heat_ES_Pearson5)

Heat_index_Pearson10 <- which(Heat_Holdout_VaR[-1]>Heat_Roll7@forecast$VaR[1:999,4])
Heat_VaR_values_Pearson10 <- Heat_Holdout_VaR[Heat_index_Pearson10]
Heat_return_values_Pearson10 <- Heat_Roll7@forecast$VaR[Heat_index_Pearson10,4]
Heat_ES_Pearson10 <- (Heat_return_values_Pearson10/Heat_VaR_values_Pearson10)
Heat_ES_Measure1_Pearson10 <- mean(Heat_ES_Pearson10)
Heat_ES_Measure2_Pearson10 <- max(Heat_ES_Pearson10)

Gas_index_Pearson1 <- which(Gas_Holdout_VaR[-1]>Gas_Roll7@forecast$VaR[1:999,4])
Gas_VaR_values_Pearson1 <- Gas_Holdout_VaR[Gas_index_Pearson1]
Gas_return_values_Pearson1 <- Gas_Roll7@forecast$VaR[Gas_index_Pearson1,4]
Gas_ES_Pearson1 <- (Gas_return_values_Pearson1/Gas_VaR_values_Pearson1)
Gas_ES_Measure1_Pearson1 <- mean(Gas_ES_Pearson1)
Gas_ES_Measure2_Pearson1 <- max(Gas_ES_Pearson1)

Gas_index_Pearson5 <- which(Gas_Holdout_VaR[-1]>Gas_Roll7@forecast$VaR[1:999,4])
Gas_VaR_values_Pearson5 <- Gas_Holdout_VaR[Gas_index_Pearson5]
Gas_return_values_Pearson5 <- Gas_Roll7@forecast$VaR[Gas_index_Pearson5,4]
Gas_ES_Pearson5 <- (Gas_return_values_Pearson5/Gas_VaR_values_Pearson5)
Gas_ES_Measure1_Pearson5 <- mean(Gas_ES_Pearson5)
Gas_ES_Measure2_Pearson5 <- max(Gas_ES_Pearson5)

Gas_index_Pearson10 <- which(Gas_Holdout_VaR[-1]>Gas_Roll7@forecast$VaR[1:999,4])
Gas_VaR_values_Pearson10 <- Gas_Holdout_VaR[Gas_index_Pearson10]
Gas_return_values_Pearson10 <- Gas_Roll7@forecast$VaR[Gas_index_Pearson10,4]
Gas_ES_Pearson10 <- (Gas_return_values_Pearson10/Gas_VaR_values_Pearson10)
Gas_ES_Measure1_Pearson10 <- mean(Gas_ES_Pearson10)
Gas_ES_Measure2_Pearson10 <- max(Gas_ES_Pearson10)

### Store for Pearson

Pearson_Full_ES_Measure1 <- matrix(c(Heat_ES_Measure1_Pearson1,Heat_ES_Measure1_Pearson5,Heat_ES_Measure1_Pearson10,
                                     Gas_ES_Measure1_Pearson1,Gas_ES_Measure1_Pearson5,Gas_ES_Measure1_Pearson10),nrow=3)

Pearson_Full_ES_Measure2 <- matrix(c(Heat_ES_Measure2_Pearson1,Heat_ES_Measure2_Pearson5,Heat_ES_Measure2_Pearson10,
                                     Gas_ES_Measure2_Pearson1,Gas_ES_Measure2_Pearson5,Gas_ES_Measure2_Pearson10),nrow=3)

rownames(Pearson_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(Pearson_Full_ES_Measure2) <- c("1%","5%","10%")
colnames(Pearson_Full_ES_Measure1) <- c("Heat","Gas")
colnames(Pearson_Full_ES_Measure2) <- c("Heat","Gas")

#### Store the values

Heat_Full_ES_Measure1 <- matrix(c(Heat_ES_Measure1_normal1,Heat_ES_Measure1_SST1,Heat_ES_Measure1_SGED1,Heat_ES_Measure1_JSU1,
                                  Heat_ES_Measure1_normal5,Heat_ES_Measure1_SST5,Heat_ES_Measure1_SGED5,Heat_ES_Measure1_JSU5,
                                  Heat_ES_Measure1_normal10,Heat_ES_Measure1_SST10,Heat_ES_Measure1_SGED10,Heat_ES_Measure1_JSU10),nrow=3,byrow=TRUE)
  
  
Gas_Full_ES_Measure1 <- matrix(c(Gas_ES_Measure1_normal1,Gas_ES_Measure1_SST1,Gas_ES_Measure1_SGED1,Gas_ES_Measure1_JSU1,
                                   Gas_ES_Measure1_normal5,Gas_ES_Measure1_SST5,Gas_ES_Measure1_SGED5,Gas_ES_Measure1_JSU5,
                                   Gas_ES_Measure1_normal10,Gas_ES_Measure1_SST10,Gas_ES_Measure1_SGED10,Gas_ES_Measure1_JSU10),nrow=3,byrow=TRUE)

Heat_Full_ES_Measure2 <- matrix(c(Heat_ES_Measure2_normal1,Heat_ES_Measure2_SST1,Heat_ES_Measure2_SGED1,Heat_ES_Measure2_JSU1,
                                   Heat_ES_Measure2_normal5,Heat_ES_Measure2_SST5,Heat_ES_Measure2_SGED5,Heat_ES_Measure2_JSU5,
                                   Heat_ES_Measure2_normal10,Heat_ES_Measure2_SST10,Heat_ES_Measure2_SGED10,Heat_ES_Measure2_JSU10),nrow=3,byrow=TRUE)

Gas_Full_ES_Measure2 <- matrix(c(Gas_ES_Measure2_normal1,Gas_ES_Measure2_SST1,Gas_ES_Measure2_SGED1,Gas_ES_Measure2_JSU1,
                                   Gas_ES_Measure2_normal5,Gas_ES_Measure2_SST5,Gas_ES_Measure2_SGED5,Gas_ES_Measure2_JSU5,
                                   Gas_ES_Measure2_normal10,Gas_ES_Measure2_SST10,Gas_ES_Measure2_SGED10,Gas_ES_Measure2_JSU10),nrow=3,byrow=TRUE)

rownames(Heat_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(Heat_Full_ES_Measure2) <- c("1%","5%","10%")
rownames(Gas_Full_ES_Measure1) <- c("1%","5%","10%")
rownames(Gas_Full_ES_Measure2) <- c("1%","5%","10%")

colnames(Heat_Full_ES_Measure1) <- c("Normal","SST","SGED","JSU")
colnames(Heat_Full_ES_Measure2) <- c("Normal","SST","SGED","JSU")
colnames(Gas_Full_ES_Measure1) <- c("Normal","SST","SGED","JSU")
colnames(Gas_Full_ES_Measure2) <- c("Normal","SST","SGED","JSU")

#### End of all the code/lines ######  
  