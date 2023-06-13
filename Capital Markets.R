rm(list=ls())
library(quantmod)
library(dplyr)
library(PerformanceAnalytics)
library(ggplot2)
library(xts)
library("fable") 
library("lubridate")
library("gridExtra")
library(tseries)
library(forecast)
library(rugarch)
library(CADFtest)
library(fGarch)
library(psych)
library(aTSA)
library(pdfetch)

mydata_xts = pdfetch_YAHOO(
  'MRK',
  fields = c("open", "high", "low", "close", "adjclose", "volume"),
  from = as.Date("2019-01-01"),
  to = as.Date("2022-12-31"),
  interval = "1d"
)

colnames(mydata_xts) <- c("MRK.Open"   ,  "MRK.High"   ,  "MRK.Low"   ,   "MRK.Close"  ,  "MRK.Adjusted","MRK.Volume")
lineChart(mydata_xts,theme = 'white', TA = c(addVo()), name = '', minor.ticks = FALSE)


#Data Frame
data <- cbind(
  Price = mydata_xts$MRK.Close,
  Return=CalculateReturns(mydata_xts$MRK.Close, method = 'log'))  #Calculating Returns and transform into log values
colnames(data) <- c('Price','Return')
head(data)
ggplot(data, aes(x = index(data), y = Return)) +
  geom_line(color = "blue", size = 1) +
  labs(x="",y = "Log Return") +
  theme_minimal()

#Distributions and statistics
describe(data, skew=TRUE,omit=TRUE)
histprice = ggplot(aes(Price), data=data) + geom_histogram(col='black',fill='lightblue', bins=50) + ggtitle('Close Price of MRK')
histreturn = ggplot(aes(Return), data=data) + geom_histogram(col='black',fill='lightblue',bins=50) + ggtitle('Log Return of MRK')
grid.arrange(histprice,histreturn, ncol = 2, nrow = 1)
jarque.bera.test(na.omit(data$Price)) #Normality Test -> not normal
jarque.bera.test(na.omit(data$Return)) #Normality Test -> not normal

#############################
#Stationarity and autocorrel#
#############################
(max.lag=round(sqrt(length(data$Price)))) # 32
CADFtest(data$Price, type= "trend", criterion= "BIC", max.lag.y=max.lag) # trend is stochastic
CADFtest(data$Price, type= "drift", criterion= "BIC", max.lag.y=max.lag) # data not stationary
CADFtest(data$Price, type= "none", criterion= "BIC", max.lag.y=max.lag) # data not stationary
CADFtest(diff(data$Price), type= "drift", criterion= "BIC", max.lag.y=max.lag) # data stationary now
CADFtest(diff(data$Price), type= "none", criterion= "BIC", max.lag.y=max.lag) # data stationary now

dprice = diff(data$Price)
plot(dprice) # looks stationary. what about white noise?

CADFtest(data$Return, type= "drift", criterion= "BIC", max.lag.y=max.lag) # log returns are stationary
CADFtest(data$Return, type= "none", criterion= "BIC", max.lag.y=max.lag) # log returns are stationary

#Charts
acfclose<- ggAcf(na.omit(data$Price), col='red',main='ACF of Close Price in levels')
pacfclose<- ggPacf(na.omit(data$Price),col='steelblue',main='PACF of Close Price in levels')
acfdclose<- ggAcf(na.omit(diff(data$Price)), col='red',main='ACF of Close Price in differences')
pacfdclose<- ggPacf(na.omit(diff(data$Price)),col='steelblue',main='PACF of Close Price in differences')
acfreturn<- ggAcf(na.omit(data$Return), col='red',main='ACF of Log Return')
pacfreturn<- ggPacf(na.omit(data$Return),col='steelblue',main='PACF of Log Return')
grid.arrange(acfclose, pacfclose,acfdclose,pacfdclose,acfreturn,pacfreturn, ncol = 2, nrow = 3)
Box.test(dprice$Price, lag = max.lag, type = "Ljung-Box") # Close prices have autocorrelation
Box.test(data$Return, lag = max.lag, type = "Ljung-Box") # have autocorrelation

#############
#Close Price#
#############
# We know that it's not stationary. Hence, we took first diff
# ACF and PACF suggest the same model specifications as for log returns
# so, we try AR(3), MA(3), ARIMA(3)
fit.close = auto.arima(data$Price)
checkresiduals(fit.close)

# Let's now perform the formal test for white noise-> the Q-test or Ljung-Box test.
Box.test(fit.close$residuals, lag=max.lag, type="Ljung-Box") # white noise
jarque.bera.test(fit.close$residuals) #Normality Test -> not normal
# looks close to normal but failed to pass jarque bera test
# Fit ARIMA model manually
fit_manual <- arima(data$Price, order = c(2,1,2))
summary(fit_manual)
arch.test(fit_manual) # hetero
forecast.close <- forecast::forecast(fit.close,h=20)
close_ts <- ts(data$Price["2019-01-02/"])
plot(close_ts, main = "", ylab = "Close Price",xlim=c(950,1050),ylim=c(100,120)) # MRK Close Price Forecast
lines(forecast.close$mean, col = "red")
lines(forecast.close$lower[,'95%'], col = "blue")
lines(forecast.close$upper[,'95%'], col = "blue")
legend("topright", legend=c("Prediction", "95% Prediction Bounds"), 
       col=c("red", "blue"), lty=1, cex=0.8)

#############
#Log Returns#
#############

fit.returns = auto.arima(data$Return)
fit.returns
checkresiduals(fit.returns)

# Let's now perform the formal test for white noise-> the Q-test or Ljung-Box test.
Box.test(fit.returns$residuals, lag=max.lag, type="Ljung-Box") # white noise
jarque.bera.test(fit.returns$residuals) #Normality Test -> not normal
# looks close to normal but failed to pass jarque bera test
# Fit ARIMA model manually
fit_manual_ret <- arima(data$Return, order = c(4,0,4),include.mean = FALSE)
summary(fit_manual_ret)
arch.test(fit_manual_ret) # hetero
# Forecasts
forecast.returns <- forecast::forecast(fit.returns,h=20)
return_ts <- ts(data$Return["2019-01-03/"])
plot(return_ts, main = "", ylab = "Log Returns",xlim=c(950,1040)) # MRK Log Returns Forecast
lines(forecast.returns$mean, col = "red")
lines(forecast.returns$lower[,'95%'], col = "blue")
lines(forecast.returns$upper[,'95%'], col = "blue")
legend("topright", legend=c("Prediction", "95% Prediction Bounds"), 
       col=c("red", "blue"), lty=1, cex=0.8)


