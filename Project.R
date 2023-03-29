#Packages used:
#tidyverse, urca, forecast, tseries, TSstudio, Quandl, GGplot2, moments

library(FinTS)
library(TSA)
library(PerformanceAnalytics)
library(tseries)
library('quantmod')
library(purrr)
library(tidyr)
library(tseries)
library(ggplot2)
library(TSstudio)
library(moments)

#Importing data

sdate <- as.Date("2000-01-01")
edate <- as.Date("2022-03-31")
ss_stock <- getSymbols('PEP',from=sdate,to=edate,auto.assign = F)#we get daily data

# Typically use previous value for NA
no.na <- which(is.na(ss_stock[,4]))  # no for NA

ss_stock1<- to.weekly(ss_stock)
head(ss_stock1)
ss_stock1[no.na,4] <- ss_stock1[no.na - 1,4]
ss_price <- ss_stock1[,4]
head(ss_price) #Final closing price data

#Price plot

ts_plot(ss_price, title = "PepsiCo price over last 22 years",
        Xtitle = "Time",
        Ytitle = "Price in USD",
        slider = TRUE)

#ADF test on closing price

adf.test(ss_price, alternative = "stationary")

head(ss_price)

#ACF of closing price
myacf = acf(ss_price)
myacf
dfd = data.frame(Column1 = c(myacf$lag), Column2 = c(myacf$acf))
dfd


#Returns calculation

pep_ret <- diff(log(ss_price))
pep1 <-na.omit(pep_ret)
ts_plot(pep1)
head(pep1)
head(ss_price)

#ADF test for returns

adf.test(pep1, alternative = "stationary")


#Distribution of returns


kurtosis(pep1$ss_stock.Close)

qqnorm(pep1, main = "QQ plot of Returns")
qqline(pep1)

boxplot(pep1, horizontal = FALSE, main = "Boxplot of Returns")

shapiro.test(as.vector(pep1))


# VaR computation

VaR(pep1, p = 0.95, method = "gaussian")
VaR(pep1, p = 0.99, method = "gaussian")
VaR(pep1, p = 0.95, method = "historical")
VaR(pep1, p = 0.99, method = "historical")
VaR(pep1, p = 0.95, method = "modified")
VaR(pep1, p = 0.99, method = "modified")

#Determining possible candidate model

acf(pep1, main = "ACF of Pepsico Log returns")
pacf(pep1, main = "PACF of Pepsico Log returns")

#Fitting ARIMA model

ARIMAfitfinal = auto.arima(pep1, approximation = FALSE, trace=TRUE)
ARIMAfitfinal
myarima1 <- arima(pep1, order = c(2,0,2))
myarima1
myarima2 <- arima(pep1, order = c(0,0,0))
myarima2
myarima3 <- arima(pep1, order = c(1,0,0))
myarima3
myarima4 <- arima(pep1, order = c(0,0,1))
myarima4
myarima5 <- arima(pep1, order = c(1,0,2))
myarima5
myarima6 <- arima(pep1, order = c(0,0,2))
myarima6
myarima7 <- arima(pep1, order = c(1,0,1))
myarima7
myarima8 <- arima(pep1, order = c(2,0,1))
myarima8
myarima9 <- arima(pep1, order = c(2,0,0))
myarima9

#White noise check


acf(ARIMAfitfinal$residuals, main = "ACF of ARIMA(1,0,1) Residuals")
qqnorm(ARIMAfitfinal$residuals, main = "QQ plot of ARIMA (1,0,1) Residuals")
qqline(ARIMAfitfinal$residuals, main = "QQ plot of ARIMA (1,0,1) Residuals")

shapiro.test(ARIMAfitfinal$residuals)



#Box cox transformation

library(MASS)

min(pep1)
pep2 = pep1+.21
head(pep2)


lambda <- BoxCox.lambda(pep2)
lambda
new_ts <- BoxCox(pep2,lambda)
new_ts
BoxCox.lambda(new_ts)
hist(new_ts, breaks = 100)
max(new_ts)

#https://stackoverflow.com/questions/26617587/finding-optimal-lambda-for-box-cox-transform-in-r



#Outlier removal 

#1st method

pepfinal = tsclean(pep1) 
hist(pepfinal)
qqnorm(pepfinal);qqline(pepfinal)
boxplot(pepfinal)

#2nd method

library(ggstatsplot)

Q <- quantile(pep1, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(pep1)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range


pep3<- subset(pep1, pep1 > low & pep1 < up)
head(pep3)

boxplot(pep3)
hist(pep3, breaks = 50)
qqnorm(pep3)
qqline(pep3)

nrow(pep1)
nrow(pep3)

library(forecast)


#Refitting ARIMA on data without outliers


ARIMAfitfinal1 = auto.arima(pep3, approximation = FALSE, trace=TRUE)
ARIMAfitfinal1

checkresiduals(ARIMAfitfinal1)
tsdiag(ARIMAfitfinal1)
arimaorder(ARIMAfitfinal1)


Box.test(ARIMAfitfinal1$residuals, type="Ljung-Box")

#Confidence intervals

confint(ARIMAfitfinal1, level = 0.95)


#Out of sample forecasts

return_forecast1 <- forecast(ARIMAfitfinal1,h=50)
plot(return_forecast1)
return_forecast1



#GARCH modelling


modelresiduals = ARIMAfitfinal1$residuals
resisquare<- modelresiduals^2
head(resisquare)



# lagrange multiplier test

ArchTest(resisquare)


# Fitting garch

library(rugarch)
pepsi_garch <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), mean.model = list(armaOrder = c(1,0,1)), distribution.model = "norm")
pepsi_garch_fit <- ugarchfit(spec = pepsi_garch, data = pep3)
pepsi_garch_fit
str(pepsi_garch_fit)
plot(pepsi_garch_fit@fit$residuals, main = "Residuals plot after fitting GARCH")




garchforecast = ugarchforecast(pepsi_garch_fit,data = pep3, n.ahead = 50)


print(sigma(garchforecast))
plot(sigma(garchforecast), main ="Volatility forecast of next 50 weeks")

#White noise check

jarque.bera.test(pepsi_garch_fit@fit$residuals)

Box.test(pepsi_garch_fit@fit$residuals^2, type="Ljung-Box")
