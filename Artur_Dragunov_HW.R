rm(list=ls())
library('dplyr')
library('forecast')
library(CADFtest)
library(vars)
library(urca)
###########################
#Time series specification#
###########################

# Source of data: https://www.fhfa.gov/DataTools/Downloads/Pages/House-Price-Index.aspx
# Data explanation https://www.fhfa.gov/Media/PublicAffairs/Pages/House-Price-Index-Frequently-Asked-Questions.aspx
US_HS_index = read.csv('HPI_master.csv') 

#filtered. Index for traditional housing of Texas is left.
Texas_HS = US_HS_index%>%filter(place_name=='Texas' & level == 'State' & hpi_flavor=='purchase-only')

all(complete.cases(Texas_HS)) # no missing observations


#time series defined. NSA stands for not seasonally adjusted
Texas_HS_ts <- ts(Texas_HS$index_nsa, frequency = 4,start = c(1991,1))
ts.plot(Texas_HS_ts) # not stationary data with a positive trend


#let's take logs and check whether the trend is deterministic or stochastic
Texas_HS_ts_log = log(Texas_HS_ts)
ts.plot(Texas_HS_ts_log)
n=dim(Texas_HS)[1]
max.lag<-round(sqrt(n))
CADFtest(Texas_HS_ts_log, type= "trend", criterion= "BIC", max.lag.y=max.lag) #The type = "trend" part of the command shows that we include
# both an intercept and a trend in the test equation.

#Stochastic trend (H0); deterministic trend (HA).
#p-value >0.05, so the trend is stochastic
#The model with only a constant and a trend is not valid. We keep on the analysis by taking log(Texas_HS) in first difference


#Let's take differences and recheck the chart
Texas_HS_ts_log_d = diff(log(Texas_HS_ts))
ts.plot(Texas_HS_ts_log_d) # we observe seasonality. The time series is still not stationary


#We check for seasonality
monthplot(Texas_HS_ts_log_d) 
#  The chart reports the seasonal plot of the log Texas Housing Index in first differences. There is a
#seasonal effect: the average Housing Index in the second quarter is higher than in the other quarters. It means that
# the real estates of Texas on average are sold for the higher prices during Q2 than for other quarters. Q4
#shows the lowest Housing index. The time series is not stationary.


#taking seasonal differences
Texas_HS_ts_log_ds <- diff(diff(log(Texas_HS_ts), lag = 4))


#We are checking now if the data are stationary
ts.plot(Texas_HS_ts_log_ds) # looks stationary. Let's confirm it with a test
nd = n-1-4 # initial nr - first difference and seasonal difference
max.lag<-round(sqrt(nd))
CADFtest(Texas_HS_ts_log_ds, type= "drift", criterion= "BIC", max.lag.y=max.lag)
#The null hypothesis states that there is a unit root. We reject it with 5% significance and say that time series is stationary


#Let's see if data are white noise by plotting correlogram and partial correlogram
par(mfrow=c(2,1))
acf(Texas_HS_ts_log_ds)
pacf(Texas_HS_ts_log_ds)
par(mfrow=c(1,1))

#There are many significant correlations and partial correlations. The time series is not a white noise.
#Both acf and pacf indicate that the SARMA model would work here because we have repetitions of auto correlations and partial autocorrelations at year 3

# We will try different model specifications and compare them: 
#1)SARIMA(2,1,0)(3,1,0), which means that we pick AR(2) which repeated 3 times. 1 in the middle indicate that we did first difference and seasonal difference
#2)SARIMA(0,1,3)(0,1,3), which means that we pick MA(3) which repeated 3 times.
#3)SARIMA(1,1,1)(3,1,3) ARIMA combination of both MA and AR models repeated 3 times
#4)SARIMA(4,1,0)(3,1,0), which means that we pick AR(4) which repeated 3 times.

##################
#Model definition#
##################

#Model 1
fit_sarma1 = arima(Texas_HS_ts_log,order=c(2,1,0),seasonal = c(3,1,0))
fit_sarma1 # both ar2 and sar3 are significant. AR2 CI: 0.248+-0.0922*2 doesn't have 0 in the interval. SAR3: -0.224+-0.1095*2 doesn't have 0 in the interval
AIC(fit_sarma1) #AIC of the model is -849.94
l=length(Texas_HS_ts)
AIC(fit_sarma1,k=log(l)) # BIC is -832.92

#Model 2
fit_sarma2 = arima(Texas_HS_ts_log,order=c(0,1,3),seasonal = c(0,1,3))
fit_sarma2 # both ma3 and sar3 are significant. MA3 CI: 0.631+-0.113*2 doesn't have 0 in the interval. SAR3: -0.320+-0.096*2 doesn't have 0 in the interval
AIC(fit_sarma2)#AIC of the model is -860.25
AIC(fit_sarma2,k=log(l)) # -840.40

#Model 3
fit_sarma3 = arima(Texas_HS_ts_log,order=c(1,1,1),seasonal = c(3,1,3))
fit_sarma3 # ar1, ma1 and sar3 are significant. BUT sma3 is not significant.
AIC(fit_sarma3)#AIC of the model is -858.40
AIC(fit_sarma3,k=log(l)) # -832.87

#Model 4
fit_sarma4 = arima(Texas_HS_ts_log,order=c(4,1,0),seasonal = c(3,1,0))
fit_sarma4 # ar4 and sar3 are significant. 
AIC(fit_sarma4) #AIC of the model is -855.43
AIC(fit_sarma4,k=log(l)) # -832.73


##################
#Model validation#
##################
#Model 1-4 correlograms of residuals
par(mfrow=c(2,2))
acf(fit_sarma1$res,plot=T,main="residual correlogram Model 1") # OK in practice. There're autocorrelations at lags 3 and 16 near borders.
acf(fit_sarma2$res,plot=T,main="residual correlogram Model 2")# OK
acf(fit_sarma3$res,plot=T,main="residual correlogram Model 3") # OK
acf(fit_sarma4$res,plot=T,main="residual correlogram Model 4") # OK. There's an autocorrelation near the border at lag 16, but it's not significant.
par(mfrow=c(1,1))
#all acf demonstrate white noises for residuals. 

# Let's now perform the formal test for white noise-> the Q-test or Ljung-Box test. Lag 10 is enough as one cycle is 4 (quarterly data):
Box.test(fit_sarma1$res, lag=10, type="Ljung-Box")
Box.test(fit_sarma2$res, lag=10, type="Ljung-Box")
Box.test(fit_sarma3$res, lag=10, type="Ljung-Box")
Box.test(fit_sarma4$res, lag=10, type="Ljung-Box")
# H0: autocorrelation=0. We don't reject H0 for all the tests and say that residuals are indeed a white noise for all 4 models.

#############
#Forecasting#
#############

#Model 1 
myforecastSARMA1<-predict(fit_sarma1,n.ahead=6)
expectedSARMA1=myforecastSARMA1$pred # 
lower_sarma1=myforecastSARMA1$pred-qnorm(0.975)*myforecastSARMA1$se
upper_sarma1=myforecastSARMA1$pred+qnorm(0.975)*myforecastSARMA1$se
#Model 2
myforecastSARMA2<-predict(fit_sarma2,n.ahead=6)
expectedSARMA2=myforecastSARMA2$pred
lower_sarma2=myforecastSARMA2$pred-qnorm(0.975)*myforecastSARMA2$se
upper_sarma2=myforecastSARMA2$pred+qnorm(0.975)*myforecastSARMA2$se
#Model 3
myforecastSARMA3<-predict(fit_sarma3,n.ahead=6)
expectedSARMA3=myforecastSARMA3$pred
lower_sarma3=myforecastSARMA3$pred-qnorm(0.975)*myforecastSARMA3$se
upper_sarma3=myforecastSARMA3$pred+qnorm(0.975)*myforecastSARMA3$se
#Model 4
myforecastSARMA4<-predict(fit_sarma4,n.ahead=6)
expectedSARMA4=myforecastSARMA4$pred
lower_sarma4=myforecastSARMA4$pred-qnorm(0.975)*myforecastSARMA4$se
upper_sarma4=myforecastSARMA4$pred+qnorm(0.975)*myforecastSARMA4$se


'Now we plot the forecasted values, together with the prediction interval.'
par(mfrow=c(2,2))
plot.ts(Texas_HS_ts_log,xlim=c(2020,2024),ylim=c(5.5,7.5),ylab="",xlab="",main='Model 1')
#visualize prediction and intervals for Model 1
lines(expectedSARMA1,col="red")
lines(lower_sarma1,col="blue")
lines(upper_sarma1,col="blue")

#visualize prediction and intervals for Model 2
plot.ts(Texas_HS_ts_log,xlim=c(2020,2024),ylim=c(5.5,7.5),ylab="",xlab="",main='Model 2')
lines(expectedSARMA2,col="red")
lines(lower_sarma2,col="blue")
lines(upper_sarma2,col="blue")

#visualize prediction and intervals for Model 3
plot.ts(Texas_HS_ts_log,xlim=c(2020,2024),ylim=c(5.5,7.5),ylab="",xlab="",main='Model 3')
lines(expectedSARMA3,col="red")
lines(lower_sarma3,col="blue")
lines(upper_sarma3,col="blue")

#visualize prediction and intervals for Model 4
plot.ts(Texas_HS_ts_log,xlim=c(2020,2024),ylim=c(5.5,7.5),ylab="",xlab="",main='Model 4')
lines(expectedSARMA4,col="red")
lines(lower_sarma4,col="blue")
lines(upper_sarma4,col="blue")
par(mfrow=c(1,1))
#all four models have similar forecasts and prediction intervals. For all the models, these intervals become wider with increasing h, which is expected to be

' Now we will compare the forecast performance of the models using absolute value loss. 
We will compute the Mean Absolute Forecast Error for both models and test if the forecast performance of the models using absolute value loss significantly
different with a Diebold-Mariaono test.'

#Now we use an expanding-window approach to forecast values of log(Texas_HS_ts) for h-periods.
#We use h = 1.

y<-Texas_HS_ts_log # 126 obs in total
S=round(0.75*length(y)) # we set out sample = 0.75 length of y. Now we are going to predict each value until we get to the last one. We calculate average mistakes and compare them
h=1 # doing prediction for the next value

#model 1
error1.h=c() # created an empty sequence
for (i in S:(length(y)-h)) # not the whole length but range starting from S until the one before the last value
{
  mymodel.sub<-arima(y[1:i], order = c(2,1,0),seasonal=c(3,1,0)) # each time we estimate a new model increasing data set by 1 element. 
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h] # y_predicted 
  error1.h<-c(error1.h,y[i+h]-predict.h) # saved all previous errors + the new one: y_observed-y_predicted
}
#model 2
error2.h<-c() 
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(0,1,3),seasonal=c(0,1,3))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h<-c(error2.h,y[i+h]-predict.h)
}
#model 3
error3.h<-c() 
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(1,1,1),seasonal=c(3,1,3), method="ML") # for model 3, without method='ML' specification it gives an error. See link: https://stackoverflow.com/questions/7233288/non-stationary-seasonal-ar-part-from-css-error-in-r
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error3.h<-c(error3.h,y[i+h]-predict.h)
}
#model 4
error4.h<-c() 
for (i in S:(length(y)-h))
{
  mymodel.sub<-arima(y[1:i], order = c(4,1,0),seasonal=c(3,1,0))
  predict.h<-predict(mymodel.sub,n.ahead=h)$pred[h]
  error4.h<-c(error4.h,y[i+h]-predict.h)
}


(MAE1 <- mean(abs(error1.h)))
(MAE2 <- mean(abs(error2.h)))
(MAE3 <- mean(abs(error3.h)))
(MAE4 <- mean(abs(error4.h)))

# as there're 4 models, we can create 6 unique pairs: 4C2, and test whether Mean absolute errors are different or equal. H0: error1=error2
dm.test(error1.h,error2.h,h=h,power=1) # MAE1 and MAE2. 
#We reject H0 and state that the forecast performance of the two models, using the absolute value loss, is significantly different. 
#Performance of the first model is better (MAE1 is significantly lower than MAE2)

dm.test(error1.h,error3.h,h=h,power=1) # MAE1 and MAE3
# We don't reject H0 and conclude that the forecast performance of the two models, using the absolute value loss, is not significantly different

dm.test(error1.h,error4.h,h=h,power=1) # MAE1 and MAE4. We don't reject H0. Forecast performance is not statistically different based on MAE
dm.test(error2.h,error3.h,h=h,power=1) # MAE2 and MAE3. We reject H0. Model 3 has a better forecast performance based on the MAE, than Model 2
dm.test(error2.h,error4.h,h=h,power=1) # MAE2 and MAE4. We reject H0. Model 4 has a better forecast performance based on the MAE, than Model 2
dm.test(error3.h,error4.h,h=h,power=1) # MAE3 and MAE4. We don't reject H0. Forecast performance is not statistically different based on MAE


############################
#Compare results in a table#
############################


df_frame=data.frame(Model=c('Model 1: SARIMA(2,1,0)(3,1,0)','Model 2: SARIMA(0,1,3)(0,1,3)','Model 3: SARIMA(1,1,1)(3,1,3)',
                            'Model 4: SARIMA(4,1,0)(3,1,0)'),
                    Highest_order_terms_significant = c('YES: AR(2), SAR(3)','YES: MA(3), SAR(3)','NO: SMA(3) not signif',
                                                        'YES: AR(4), SAR(3)'),
                    Correlogram_residuals = c('OK in practicec','OK','OK','OK in practice'),
                    Box_test = c('White noise','White noise','White noise','White noise'),
                    AIC = c(AIC(fit_sarma1),AIC(fit_sarma2),AIC(fit_sarma3),AIC(fit_sarma4)),
                    BIC = c(AIC(fit_sarma1,k=log(l)),AIC(fit_sarma2,k=log(l)),AIC(fit_sarma3,k=log(l)),AIC(fit_sarma4,k=log(l))),
                    MAE = c(MAE1,MAE2,MAE3,MAE4),
                    Nr_of_param=c(5,6,8,7))


df_frame # prepared table

# Based on the table, Model 2 has the lowest AIC and BIC. But it has the highest MAE which is significantly different than other MAE results, see dm.test outputs.
# So, Model 2 is not a winner. Model 3 doesn't have all the highest order terms significant. This model is also the most complicated as it has 8 parameters.
# So, Model 3 is not a winner. Model 1 and Model 4 have in practice good correlograms and their forecasting performance is not statistically different. These models
# also have almost the same BIC performances (-833 if we round to the nearest integer). BUT Model 1 has the worst AIC even though it has the lowest # of parameters.
# So, out of these 4 models, we would pick Model 4 as a winner.


########################################Multivariate time series analysis###############################################

###########################
#Time series specification#
###########################

#For this part let's conduct a multivariate time series analysis on three time series. 
#Texas House Price Index was introduced before. 30-year fixed rate mortgage average in the US and M2 - money stock measure as a proxy for Money Supply in the US
mortgage_rate = read.csv('MORTGAGE30US.csv') # Quarterly (end of period) time series of mortgage percentage points. Not seasonally adjusted. Period 1991Q1-2022Q2
M2_rate = read.csv('M2SL.csv') # Quarterly (end of period) seasonally adjusted billions of dollars. Period 1991Q1-2022Q2
#Both time series were retrieved from Federal Reserve Economic Data https://fred.stlouisfed.org/series/MORTGAGE30US; https://fred.stlouisfed.org/series/M2SL#0
# A research hypothesis is that these three time series should be interconnected. M2 and mortgage rates should have a significant relationship with Texas House Price Index.
# We already know that log(Texas_HS_ts) should be first and seasonally differenced to be stationary. Let's observe if other time series are stationary.
mortgage_rate_ts =  ts(mortgage_rate$MORTGAGE30US, frequency = 4,start = c(1991,1))
M2_rate_ts_log = ts(log(M2_rate$M2SL), frequency = 4,start = c(1991,1))
par(mfrow=c(2,1))
ts.plot(mortgage_rate_ts) # negative trend with an increase in 2022
ts.plot(M2_rate_ts_log) # positive trend. Both time series are not stationary.
par(mfrow=c(1,1))
# let's test if trends are deterministic or stochastic

CADFtest(mortgage_rate_ts, type= "trend", criterion= "BIC", max.lag.y=max.lag) #trend is stochastic
CADFtest(M2_rate_ts_log, type= "trend", criterion= "BIC", max.lag.y=max.lag) #trend is stochastic

#we take first difference and check if data is I(1)
M2_rate_ts_log_d = diff(M2_rate_ts_log)
mortgage_rate_ts_d = diff(mortgage_rate_ts)
ts.plot(mortgage_rate_ts_d)
ts.plot(M2_rate_ts_log_d)
CADFtest(mortgage_rate_ts_d, type= "drift", criterion= "BIC", max.lag.y=max.lag) #stationary
CADFtest(M2_rate_ts_log_d, type= "drift", criterion= "BIC", max.lag.y=max.lag) #stationary
#both time series are I(1)

######
#VECM#
######

#We will try first a vector error correction model to test for cointegration (long-term relationship). If there's no such, then we will proceed with VAR model on differenced data

#We will use Johansen's approach to test for cointegration among log(M2), log(Texas_HS) and mortgage_rate.
'We will use automatic lag selection to determine the optimal lag length needed in the Johansen
test equation by running a VAR in levels'

df_ts = data.frame(M2_rate_ts_log,Texas_HS_ts_log,mortgage_rate_ts)
VARselect(df_ts,lag.max=10,type="const")
#The order of the VAR model in levels selected by Schwarz Information Criterion is 6 (the lowest SIC). Such a big order might be due to the fact that the recent quarters have big fluctuations
#for real estate, M2 and mortgage rates
"Now we test for cointegration using Johansen's trace test statistic. We will use the ca.jo function from the urca package:"

trace_test<-ca.jo(df_ts,type="trace",K=6,ecdet="const",spec="transitory")
summary(trace_test)

'For r=0 and r=1, the test statistics is larger than then the critical value (52.02 > 34.91; 21.19>19.96), thus there
are two cointegrating relations. For r=2, the test statistics is smaller than then the critical value (1.61 < 9.24).
Hence, log(M2), log(Texas_HS) and log(mortgage_rate) are cointegrated and we will proceed with VECM model'
maxeigen_test<-ca.jo(df_ts,type="eigen",K=6,ecdet="const",spec="transitory")
summary(maxeigen_test) # same results with maximum eigenvalue test statistic

#####################################
#Model specification and forecasting#
#####################################
fit_vecm1<-cajorls(trace_test,r=2)
fit_vecm1
summary(fit_vecm1$rlm)
#We see that not all Error Correcting Terms are negative, which is against the theoretical expectations. BUT it may be that these coefficients are not statistically significant.

# $beta: estimated cointegration vectors. log(M2_rate)=-0.289*mortage_rate+10.058 
#Log(Texas_HS) = -0.166*mortgage_rate+5.630

# Now we use the VECM to forecast the values of log(M2), log(Texas_HS) and mortgage_rate for the next 6 quarters. 
#We first have to retransform the VECM to a VAR specification in levels:

fit_var<-vec2var(trace_test,r=2)
myforecast<-predict(fit_var,n.ahead=6)

'A plot of the forecast for log(M2) together with the 95% prediction interval is below.'
par(mfrow=c(3,2))
ts.plot(M2_rate_ts_log)
M2_rate_ts_log_forecast<-ts(myforecast$fcst$M2_rate_ts_log[,1],frequency=4,start=c(2022,3))
M2_rate_ts_log_lower<-ts(myforecast$fcst$M2_rate_ts_log[,2],frequency=4,start=c(2022,3))
M2_rate_ts_log_upper<-ts(myforecast$fcst$M2_rate_ts_log[,3],frequency=4,start=c(2022,3))
ts.plot(M2_rate_ts_log_forecast,M2_rate_ts_log_lower,M2_rate_ts_log_upper,col=c("black","red","red"))
title(main = "6-step-ahead forecast of log(M2)",cex.main=1) # based on the forecast, money supply will continue smoothly growing

'A plot of the forecast for log(Texas_HS) together with the 95% prediction interval is below.'
plot.ts(Texas_HS_ts_log)
Texas_HS_ts_log_forecast<-ts(myforecast$fcst$Texas_HS_ts_log[,1],frequency=4,start=c(2022,3))
Texas_HS_ts_log_lower<-ts(myforecast$fcst$Texas_HS_ts_log[,2],frequency=4,start=c(2022,3))
Texas_HS_ts_log_upper<-ts(myforecast$fcst$Texas_HS_ts_log[,3],frequency=4,start=c(2022,3))
ts.plot(Texas_HS_ts_log_forecast,Texas_HS_ts_log_lower,Texas_HS_ts_log_upper,col=c("black","red","red"))
title(main = "6-step-ahead forecast of log(Texas_HS)",cex.main=1) # Texas House index should continue growing slowly till the end of the year but then lie on a plateau in 2023, based on the forecast.
# The result is similar to SARIMA predictions.

'A plot of the forecast for mortgage_rate together with the 95% prediction interval is below.'
plot.ts(mortgage_rate_ts)
mortgage_rate_ts_forecast<-ts(myforecast$fcst$mortgage_rate_ts[,1],frequency=4,start=c(2022,3))
mortgage_rate_ts_lower<-ts(myforecast$fcst$mortgage_rate_ts[,2],frequency=4,start=c(2022,3))
mortgage_rate_ts_upper<-ts(myforecast$fcst$mortgage_rate_ts[,3],frequency=4,start=c(2022,3))
ts.plot(mortgage_rate_ts_forecast,mortgage_rate_ts_lower,mortgage_rate_ts_upper,col=c("black","red","red"))
title(main = "6-step-ahead forecast of mortgage rate",cex.main=1)
#Forecast demonstrates that the mortgage rates should decrease and move back to its downtrend. However, the forecast goes against the expectations of the Federal Reserve,
#which intends to keep increasing the interest rates at least till the end of 2023 https://www.bankrate.com/banking/federal-reserve/when-will-the-fed-stop-raising-rates/ 
