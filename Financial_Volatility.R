#### Packages, wd, general setup ####
# Load the required packages 
library(xts)
library(quantmod)
library(vrtest)
library(urca)
library(ggplot2)
library(rugarch)
library(readxl)
library(vars)
library(tsDyn)
library(zoo)
library(lubridate)
library(highfrequency)
library(forecast)
library(strucchange)
library(grid)
library(gtable)
library(highfrequency)
library(quantmod)
library(TTR)
library(qqplotr)
library(ggpubr)
library(moments)



skewness(NEE_returns)
kurtosis(NEE_returns)

# Get the directory path of the current code
PATH <- dirname(rstudioapi::getSourceEditorContext()$path)
# Set the working directory to where the code is 
setwd(PATH)

# Loading the crude oil prices
WTI <- read_excel("C:/Users/marco/Documents/HSG/Financial Volatility/WTI.xlsx")
WTI <- xts(WTI[,-1], order.by=WTI$Date)
WTI <- diff(log(WTI))[-1]



#### XOM ####
# Loading the daily closing prices of XOM from yahoo finance
XOM_daily <- getSymbols("XOM", src="yahoo", from="2016-01-01", to="2023-05-13", auto.assign=FALSE)
XOM_daily <- XOM_daily$XOM.Close

# Converting to log returns
XOM_returns <- diff(log(XOM_daily))[-1]

# Shapiro-Wilkins normality test 
shapiro.test(as.numeric(XOM_returns)) # H0 is rejected, time-series is not normally distributed
# Variance-ratio test for log of XOM
Auto.VR(XOM_returns) # H0 accepted Variance is constant, ARCH/GARCH warranted

# ADF test 
T <- length(XOM_returns)
lag.max <- as.integer(12*(T/100)^0.25)
summary(ur.df(XOM_returns, type="trend", lags=lag.max, selectlags="BIC"))
#KPSS test
summary(ur.kpss(XOM_returns, type="tau", lags="long")) # Stationary 

# ACF and PACF 
acf(XOM_returns)
pacf(XOM_returns)
# Plotting the time-series
plot(XOM_returns)

# Testing for white-noise with Ljung-Box test
Box.test(XOM_returns, type = "Ljung-Box") # White noise

# GARCH models
acf(XOM_returns^2)
pacf(XOM_returns^2)

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE), 
                    distribution.model = "sstd")
XOM_SGARCH <- ugarchfit(spec1, data = XOM_returns)
show(XOM_SGARCH)


# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(5, 2), include.mean = TRUE), 
                   distribution.model = "sstd")
XOM_GJR <- ugarchfit(spec2, data = XOM_returns)
show(XOM_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE), 
                    distribution.model = "sstd")
XOM_EGARCH <- ugarchfit(spec3, data = XOM_returns)
show(XOM_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE), 
                    distribution.model = "sstd")
XOM_TGARCH <- ugarchfit(spec4, data = XOM_returns)
show(XOM_TGARCH)

# Comparing the models
# Likelihood
likelihood(XOM_SGARCH)
likelihood(XOM_GJR)
likelihood(XOM_EGARCH)
likelihood(XOM_TGARCH)

# AIC/BIC
infocriteria(XOM_SGARCH)
infocriteria(XOM_GJR)
infocriteria(XOM_EGARCH)
infocriteria(XOM_TGARCH)


#### CVX ####
# Loading the daily closing prices of CVX from yahoo finance
CVX_daily <- getSymbols("CVX", src="yahoo", from="2016-01-01", to="2023-05-13", auto.assign=FALSE)
CVX_daily <- CVX_daily$CVX.Close

# Converting to log returns
CVX_returns <- diff(log(CVX_daily))[-1]

# Shapiro-Wilkins normality test 
shapiro.test(as.numeric(CVX_returns)) # H0 is rejected, time-series is not normally distributed
# Variance-ratio test for log of CVX
Auto.VR(CVX_returns) # H0 rejected Variance is not constant, ARCH/GARCH warranted

# ADF test 
T <- length(CVX_returns)
lag.max <- as.integer(12*(T/100)^0.25)
summary(ur.df(CVX_returns, type="trend", lags=lag.max, selectlags="BIC"))
#KPSS test
summary(ur.kpss(CVX_returns, type="tau", lags="long")) # not Stationary 

# ACF and PACF 
acf(CVX_returns)
pacf(CVX_returns)
# Plotting the time-series
plot(CVX_returns)

# Testing for white-noise with Ljung-Box test
Box.test(CVX_returns, type = "Ljung-Box") # not White noise

# Finding the optimal ARIMA parameters
auto.arima(CVX_returns) # ARIMA(5,0,2)

# ARIMA(5,1,2):
CVX_ARIMA <- Arima(CVX_returns, order=c(5,0,2))
coeftest(CVX_ARIMA)
# modulus
abs(polyroot(c(1, -CVX_ARIMA$model$phi)))
abs(polyroot(c(1, CVX_ARIMA$model$theta)))
checkresiduals(CVX_ARIMA, lag=30)
plot(CVX_ARIMA)
tsdiag(CVX_ARIMA)
acf(CVX_ARIMA$residuals)
pacf(CVX_ARIMA$residuals)
BIC(CVX_ARIMA)
plot(CVX_ARIMA$residuals)


# GARCH models
acf(CVX_returns^2)
pacf(CVX_returns^2)

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE), 
                    distribution.model = "sstd")
CVX_SGARCH <- ugarchfit(spec1, data = CVX_returns)
show(CVX_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE), 
                    distribution.model = "sstd")
CVX_GJR <- ugarchfit(spec2, data = CVX_returns)
show(CVX_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE), 
                    distribution.model = "sstd")
CVX_EGARCH <- ugarchfit(spec3, data = CVX_returns)
show(CVX_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE), 
                    distribution.model = "sstd")
CVX_TGARCH <- ugarchfit(spec4, data = CVX_returns)
show(CVX_TGARCH)

# Comparing the models
# Likelihood
likelihood(CVX_SGARCH)
likelihood(CVX_GJR)
likelihood(CVX_EGARCH)
likelihood(CVX_TGARCH)

# AIC/BIC
infocriteria(CVX_SGARCH)
infocriteria(CVX_GJR)
infocriteria(CVX_EGARCH)
infocriteria(CVX_TGARCH)


#### NEE ####
# Loading the daily closing prices of NEE from yahoo finance
NEE_daily <- getSymbols("NEE", src="yahoo", from="2016-01-01", to="2023-05-13", auto.assign=FALSE)
NEE_daily <- NEE_daily$NEE.Close

# Converting to log returns
NEE_returns <- diff(log(NEE_daily))[-1]

# Shapiro-Wilkins normality test 
shapiro.test(as.numeric(NEE_returns)) # H0 is rejected, time-series is not normally distributed
# Variance-ratio test for log of NEE
Auto.VR(NEE_returns) # H0 rejected Variance is not constant, ARCH/GARCH warranted

# ADF test 
T <- length(NEE_returns)
lag.max <- as.integer(12*(T/100)^0.25)
summary(ur.df(NEE_returns, type="trend", lags=lag.max, selectlags="BIC"))
#KPSS test
summary(ur.kpss(NEE_returns, type="tau", lags="long")) # not Stationary 

# ACF and PACF 
acf(NEE_returns)
pacf(NEE_returns)
# Plotting the time-series
plot(NEE_returns)

# Testing for white-noise with Ljung-Box test
Box.test(NEE_returns, type = "Ljung-Box") # not White noise

# Finding the optimal ARIMA parameters
auto.arima(NEE_returns) # ARIMA(3,0,3)

# ARIMA(3,1,3):
NEE_ARIMA <- Arima(NEE_returns, order=c(3,0,3))
coeftest(NEE_ARIMA)
# modulus
abs(polyroot(c(1, -NEE_ARIMA$model$phi)))
abs(polyroot(c(1, NEE_ARIMA$model$theta)))
checkresiduals(NEE_ARIMA, lag=30)
plot(NEE_ARIMA)
tsdiag(NEE_ARIMA)
acf(NEE_ARIMA$residuals)
pacf(NEE_ARIMA$residuals)
BIC(NEE_ARIMA)
plot(NEE_ARIMA$residuals)


# GARCH models
acf(NEE_returns^2)
pacf(NEE_returns^2)

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE), 
                    distribution.model = "sstd")
NEE_SGARCH <- ugarchfit(spec1, data = NEE_returns)
show(NEE_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE), 
                    distribution.model = "sstd")
NEE_GJR <- ugarchfit(spec2, data = NEE_returns)
show(NEE_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE), 
                    distribution.model = "sstd")
NEE_EGARCH <- ugarchfit(spec3, data = NEE_returns)
show(NEE_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE), 
                    distribution.model = "sstd")
NEE_TGARCH <- ugarchfit(spec4, data = NEE_returns)
show(NEE_TGARCH)
plot(residuals(NEE_TGARCH))

# Comparing the models
# Likelihood
likelihood(NEE_SGARCH)
likelihood(NEE_GJR)
likelihood(NEE_EGARCH)
likelihood(NEE_TGARCH)

# AIC/BIC
infocriteria(NEE_SGARCH)
infocriteria(NEE_GJR)
infocriteria(NEE_EGARCH)
infocriteria(NEE_TGARCH)

#### Structural break ####
# XOM
# Use the breakpoints function to detect structural breaks
#time <- 1:length(XOM_daily)
#bp <- breakpoints(XOM_daily ~ time, breaks=1)
#print(bp) # Breakpoint observation: 1050 (2020.03.05)
#XOM_daily[1050]


# Create the structural change dummy
XOM_breakpoint <- as.Date("2020-03-05")
XOM_dummy0 <- as.numeric(index(XOM_returns) > XOM_breakpoint)
# Bind the dummy variable as an external regressor
XOM_dummy <- as.matrix(cbind(XOM_dummy0))

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = XOM_dummy), 
                    distribution.model = "sstd")
XOM_SGARCH <- ugarchfit(spec1, data = XOM_returns)
show(XOM_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = XOM_dummy), 
                    distribution.model = "sstd")
XOM_GJR <- ugarchfit(spec2, data = XOM_returns)
show(XOM_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = XOM_dummy), 
                    distribution.model = "sstd")
XOM_EGARCH <- ugarchfit(spec3, data = XOM_returns)
show(XOM_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = XOM_dummy), 
                    distribution.model = "sstd")
XOM_TGARCH <- ugarchfit(spec4, data = XOM_returns)
show(XOM_TGARCH)

# Comparing the models
# Likelihood
likelihood(XOM_SGARCH)
likelihood(XOM_GJR)
likelihood(XOM_EGARCH)
likelihood(XOM_TGARCH)

# AIC/BIC
infocriteria(XOM_SGARCH)
infocriteria(XOM_GJR)
infocriteria(XOM_EGARCH)
infocriteria(XOM_TGARCH)




# CVX
#time <- 1:length(CVX_daily)
#bp <- breakpoints(CVX_daily ~ time, breaks=1)
#print(bp) # Breakpoint observation: 1044 (2020.02.26)
#CVX_daily[1044]

# Create the structural change dummy
CVX_breakpoint <- as.Date("2020-02-26")
CVX_dummy0 <- as.numeric(index(CVX_returns) > CVX_breakpoint)
# Bind the dummy variable as an external regressor
CVX_dummy <- as.matrix(cbind(CVX_dummy0))

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = CVX_dummy), 
                    distribution.model = "sstd")
CVX_SGARCH <- ugarchfit(spec1, data = CVX_returns)
show(CVX_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = CVX_dummy), 
                    distribution.model = "sstd")
CVX_GJR <- ugarchfit(spec2, data = CVX_returns)
show(CVX_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = CVX_dummy), 
                    distribution.model = "sstd")
CVX_EGARCH <- ugarchfit(spec3, data = CVX_returns)
show(CVX_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = CVX_dummy), 
                    distribution.model = "sstd")
CVX_TGARCH <- ugarchfit(spec4, data = CVX_returns)
show(CVX_TGARCH)

# Comparing the models
# Likelihood
likelihood(CVX_SGARCH)
likelihood(CVX_GJR)
likelihood(CVX_EGARCH)
likelihood(CVX_TGARCH)

# AIC/BIC
infocriteria(CVX_SGARCH)
infocriteria(CVX_GJR)
infocriteria(CVX_EGARCH)
infocriteria(CVX_TGARCH)



# NEE
#time <- 1:length(NEE_daily)
#bp <- breakpoints(NEE_daily ~ time, breaks=1)
#print(bp) # Breakpoint observation: 1142 (2020.07.16)
#NEE_daily[1142]

# Create the structural change dummy
NEE_breakpoint <- as.Date("2020-07-16")
NEE_dummy0 <- as.numeric(index(NEE_returns) > NEE_breakpoint)
# Bind the dummy variable as an external regressor
NEE_dummy <- as.matrix(cbind(NEE_dummy0))

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = NEE_dummy), 
                    distribution.model = "sstd")
NEE_SGARCH <- ugarchfit(spec1, data = NEE_returns)
show(NEE_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = NEE_dummy), 
                    distribution.model = "sstd")
NEE_GJR <- ugarchfit(spec2, data = NEE_returns)
show(NEE_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = NEE_dummy), 
                    distribution.model = "sstd")
NEE_EGARCH <- ugarchfit(spec3, data = NEE_returns)
show(NEE_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = NEE_dummy), 
                    distribution.model = "sstd")
NEE_TGARCH <- ugarchfit(spec4, data = NEE_returns)
show(NEE_TGARCH)
plot(residuals(NEE_TGARCH))

# Comparing the models
# Likelihood
likelihood(NEE_SGARCH)
likelihood(NEE_GJR)
likelihood(NEE_EGARCH)
likelihood(NEE_TGARCH)

# AIC/BIC
infocriteria(NEE_SGARCH)
infocriteria(NEE_GJR)
infocriteria(NEE_EGARCH)
infocriteria(NEE_TGARCH)


#### Exogenous WTI ####
external_matrix <- as.matrix(cbind(WTI))

# XOM
# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
XOM_SGARCH <- ugarchfit(spec1, data = XOM_returns)
show(XOM_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
XOM_GJR <- ugarchfit(spec2, data = XOM_returns)
show(XOM_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
XOM_EGARCH <- ugarchfit(spec3, data = XOM_returns)
show(XOM_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
XOM_TGARCH <- ugarchfit(spec4, data = XOM_returns)
show(XOM_TGARCH)

# Comparing the models
# Likelihood
likelihood(XOM_SGARCH)
likelihood(XOM_GJR)
likelihood(XOM_EGARCH)
likelihood(XOM_TGARCH)

# AIC/BIC
infocriteria(XOM_SGARCH)
infocriteria(XOM_GJR)
infocriteria(XOM_EGARCH)
infocriteria(XOM_TGARCH)

# CVX
# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
CVX_SGARCH <- ugarchfit(spec1, data = CVX_returns)
show(CVX_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
CVX_GJR <- ugarchfit(spec2, data = CVX_returns)
show(CVX_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
CVX_EGARCH <- ugarchfit(spec3, data = CVX_returns)
show(CVX_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
CVX_TGARCH <- ugarchfit(spec4, data = CVX_returns)
show(CVX_TGARCH)

# Comparing the models
# Likelihood
likelihood(CVX_SGARCH)
likelihood(CVX_GJR)
likelihood(CVX_EGARCH)
likelihood(CVX_TGARCH)

# AIC/BIC
infocriteria(CVX_SGARCH)
infocriteria(CVX_GJR)
infocriteria(CVX_EGARCH)
infocriteria(CVX_TGARCH)

# NEE
# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
NEE_SGARCH <- ugarchfit(spec1, data = NEE_returns)
show(NEE_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
NEE_GJR <- ugarchfit(spec2, data = NEE_returns, solver = "hybrid")
show(NEE_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
NEE_EGARCH <- ugarchfit(spec3, data = NEE_returns)
show(NEE_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = external_matrix), 
                    distribution.model = "sstd")
NEE_TGARCH <- ugarchfit(spec4, data = NEE_returns)
show(NEE_TGARCH)
plot(residuals(NEE_TGARCH))

# Comparing the models
# Likelihood
likelihood(NEE_SGARCH)
likelihood(NEE_GJR)
likelihood(NEE_EGARCH)
likelihood(NEE_TGARCH)

# AIC/BIC
infocriteria(NEE_SGARCH)
infocriteria(NEE_GJR)
infocriteria(NEE_EGARCH)
infocriteria(NEE_TGARCH)




#### Structural break + exogenous WTI #### 
# XOM
XOM_WTI <- as.matrix(cbind(WTI,XOM_dummy0))

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = XOM_WTI), 
                    distribution.model = "sstd")
XOM_SGARCH <- ugarchfit(spec1, data = XOM_returns)
show(XOM_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = XOM_WTI), 
                    distribution.model = "sstd")
XOM_GJR <- ugarchfit(spec2, data = XOM_returns)
show(XOM_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = XOM_WTI), 
                    distribution.model = "sstd")
XOM_EGARCH <- ugarchfit(spec3, data = XOM_returns)
show(XOM_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(0, 0), include.mean = TRUE,
                                      external.regressors = XOM_WTI), 
                    distribution.model = "sstd")
XOM_TGARCH <- ugarchfit(spec4, data = XOM_returns)
show(XOM_TGARCH)

# Comparing the models
# Likelihood
likelihood(XOM_SGARCH)
likelihood(XOM_GJR)
likelihood(XOM_EGARCH)
likelihood(XOM_TGARCH)

# AIC/BIC
infocriteria(XOM_SGARCH)
infocriteria(XOM_GJR)
infocriteria(XOM_EGARCH)
infocriteria(XOM_TGARCH)




# CVX
# Bind the dummy variable as an external regressor
CVX_WTI <- as.matrix(cbind(WTI, CVX_dummy0))

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = CVX_WTI), 
                    distribution.model = "sstd")
CVX_SGARCH <- ugarchfit(spec1, data = CVX_returns)
show(CVX_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = CVX_WTI), 
                    distribution.model = "sstd")
CVX_GJR <- ugarchfit(spec2, data = CVX_returns)
show(CVX_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = CVX_WTI), 
                    distribution.model = "sstd")
CVX_EGARCH <- ugarchfit(spec3, data = CVX_returns)
show(CVX_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(5, 2), include.mean = TRUE,
                                      external.regressors = CVX_WTI), 
                    distribution.model = "sstd")
CVX_TGARCH <- ugarchfit(spec4, data = CVX_returns)
show(CVX_TGARCH)

# Comparing the models
# Likelihood
likelihood(CVX_SGARCH)
likelihood(CVX_GJR)
likelihood(CVX_EGARCH)
likelihood(CVX_TGARCH)

# AIC/BIC
infocriteria(CVX_SGARCH)
infocriteria(CVX_GJR)
infocriteria(CVX_EGARCH)
infocriteria(CVX_TGARCH)



# NEE
NEE_WTI <- as.matrix(cbind(WTI, NEE_dummy0))

# Standard GARCH
spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = NEE_WTI), 
                    distribution.model = "sstd")
NEE_SGARCH <- ugarchfit(spec1, data = NEE_returns)
show(NEE_SGARCH)

# GJR-GARCH
spec2 <- ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = NEE_WTI), 
                    distribution.model = "sstd")
NEE_GJR <- ugarchfit(spec2, data = NEE_returns)
show(NEE_GJR)

# E-GARCH
spec3 <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = NEE_WTI), 
                    distribution.model = "sstd")
NEE_EGARCH <- ugarchfit(spec3, data = NEE_returns)
show(NEE_EGARCH)

# T-GARCH
spec4 <- ugarchspec(variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1, 1)),
                    mean.model = list(armaOrder = c(3, 3), include.mean = TRUE,
                                      external.regressors = NEE_WTI), 
                    distribution.model = "sstd")
NEE_TGARCH <- ugarchfit(spec4, data = NEE_returns)
show(NEE_TGARCH)
plot(residuals(NEE_TGARCH))

# Comparing the models
# Likelihood
likelihood(NEE_SGARCH)
likelihood(NEE_GJR)
likelihood(NEE_EGARCH)
likelihood(NEE_TGARCH)

# AIC/BIC
infocriteria(NEE_SGARCH)
infocriteria(NEE_GJR)
infocriteria(NEE_EGARCH)
infocriteria(NEE_TGARCH)

##### Rolling window forecast ####
rolling_window_forecast <- function(spec, returns_data, start_date, end_date) {
  
  # Identify the indices for the training and testing periods
  start_training_index <- min(which(index(returns_data) >= as.Date(start_date)))
  end_training_index <- max(which(index(returns_data) <= as.Date(end_date)))
  
  # Number of predictions
  n <- length(returns_data) - end_training_index
  
  # Preallocate vector to store predictions
  predicted_values <- vector(mode="numeric", length=n)
  
  # Rolling window forecast
  for(i in 1:n) {
    
    # Create training set
    training_set <- returns_data[1:(end_training_index + i - 1), ]
    
    # Fit GARCH model
    fit <- ugarchfit(spec, data = training_set)
    
    # Forecast one step ahead
    forecast <- ugarchforecast(fit, n.ahead = 1)
    
    # Store predicted value
    predicted_values[i] <- forecast@forecast$seriesFor[1]
  }
  # Calculate the out-of-sample MSE
  actual_values <- as.numeric(returns_data[(end_training_index + 1):(end_training_index + n), ])
  mse <- mean((predicted_values - actual_values)^2)
  
  # Return the predictions and the MSE
  return(list(predicted_values = predicted_values, mse = mse))
}

# Call the function with the GARCH specification and the returns data
predicted_values <- rolling_window_forecast(spec1, XOM_returns, "2016-01-05", "2022-05-11")

# Print the predictions
print(predicted_values$predicted_values)
print(predicted_values$mse)

# DM test
y_real <- as.data.frame(tail(XOM_returns, 252))
xom_normal_pred <- predicted_values$predicted_values
xom_gjr_pred <- predicted_values$predicted_values
dm_test <- dm.test(y_real, xom_normal_pred, xom_gjr_pred, alternative = "two.sided")
summary(dm_test)

#### test ####
# Define start and end dates
start_date <- as.Date("2016-01-05")
end_date <- as.Date("2022-05-11")

# Extract the portion of the data for training
train_data <- CVX_returns[index(CVX_returns) >= start_date & index(CVX_returns) <= end_date]

# Initialize a list to store the results
results <- list()

train_and_forecast <- function(spec, data, start_date, end_date) {
  # Identify the indices for the training and testing periods
  start_training_index <- min(which(index(data) >= as.Date(start_date)))
  end_training_index <- max(which(index(data) <= as.Date(end_date)))
  
  # Number of predictions
  n <- length(data) - end_training_index
  
  # Preallocate vector to store predictions
  predicted_values <- vector(mode = "numeric", length = n)
  
  # Extract ARMA order from spec
  armaOrder <- spec@mean.model$armaOrder
  
  # Rolling window forecast
  for (i in 1:n) {
    # Create training set
    training_set <- data[1:(end_training_index + i - 1)]
    
    # Update spec with the same ARMA order
    spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                       mean.model = list(armaOrder = armaOrder, include.mean = TRUE), 
                       distribution.model = "sstd")
    
    # Fit GARCH model
    fit <- ugarchfit(spec, data = training_set)
    
    # Forecast one step ahead
    forecast <- ugarchforecast(fit, n.ahead = 1)
    
    # Store predicted value
    predicted_values[i] <- forecast@forecast$seriesFor[1]
  }
  
  # Calculate the out-of-sample MSE
  actual_values <- data[(end_training_index + 1):(end_training_index + n)]
  mse <- calculate_MSE(actual_values, predicted_values)
  
  # Return the predictions and the MSE
  return(list(predicted_values = predicted_values, mse = mse))
}

# Function to calculate Mean Squared Error (MSE)
calculate_MSE <- function(actual, forecast) {
  mse <- mean((actual - forecast)^2)
  return(mse)
}


forecast_results <- train_and_forecast(spec1, train_data, as.character(start_date), as.character(end_date))
mse <- forecast_results$mse

#### VAR/VECM ####
# Creating the dataframe for VAR/VECM testing
XOM_vol <- as.data.frame(sigma(XOM_EGARCH))
CVX_vol <- as.data.frame(sigma(CVX_EGARCH))
NEE_vol <- as.data.frame(sigma(NEE_TGARCH))
vardata <- cbind(XOM_vol, CVX_vol, NEE_vol)
colnames(vardata) <- c("XOM_vol", "CVX_vol", "NEE_vol")

# Testing for cointegration
jotest <- ca.jo(vardata, type="eigen")
summary(jotest)
jotest2 <- ca.jo(vardata, type = "trace")
summary(jotest2)
# r=2, there is cointegration

# VECM model lag select
vardata <- as.data.frame(vardata)
lag_select <- VARselect(vardata, lag.max = 10, type = "const")
lag_select$selection
# lag 8 or 5

# VECM model
WTI <- as.data.frame(WTI)
model1 <- VECM(vardata[,1:3], lag = 8, include = "const", estim = "ML", r = 2, exogen = WTI)
summary(model1)

# Testing the VECM model
vecmresid <- as.data.frame(resid(model1))
lapply(vecmresid, function(i)Box.test(i,lag=2,type = "Ljung-Box"))

# IVF
irf_model <- irf(model1, impulse = "WTI", response = c("XOM_vol", "CVX_vol", "NEE_vol"), n.ahead = 10)

model1VAR <- vec2var(jotest, r = 2)
plot(irf(model1VAR, impulse = "XOM_vol", n.ahead = 30, ortho = F))

# Variance decomposition
fevd_total <- vars:: fevd(model1, n.ahead = 2)
plot(vars:: fevd(model1))

#### Multivariate GARCH ####

#### HAR ####
# XOM
XOM_intraday <- read_excel("C:/Users/Marco/Documents/HSG/Financial Volatility/XOM.xlsx", col_types = c("date", "numeric"))
rownames(XOM_intraday) <- XOM_intraday$Date
XOM_intraday <- xts(XOM_intraday[, -1], order.by = XOM_intraday$Date)
XOM_intraday <- diff(log(XOM_intraday))[-1]



#### Plots ####
# Plotting WTI prices
WTI2 <- read_excel("C:/Users/marco/Documents/HSG/Financial Volatility/WTI.xlsx")
WTI2$Date <- as.Date(WTI2$Date, format = "%Y/%m/%d")

ggplot(data = WTI2, aes(x = Date, y = Close)) +
  geom_line(color = "black", size = 0.5) +  
  labs(
    y = "Price (USD)", 
    x ="",
    title = "WTI crude oil daily closing prices",
    subtitle = "2016.01.04-2023.05.12",
    caption = "Source: Yahoo Finance"
  ) +
  theme_bw() + 
  theme(
    plot.title = element_text(size= 14, hjust = 0.5, face ="bold"), 
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold") 
  )

# Plot the time series with ggplot
XOM_daily_df <- data.frame(time=index(XOM_daily), XOM_value=coredata(XOM_daily))
CVX_daily_df <- data.frame(CVX_value=coredata(CVX_daily))
NEE_daily_df <- data.frame(NEE_value=coredata(NEE_daily))

# Combine them into a single data frame
data <- cbind(XOM_daily_df, CVX_daily_df, NEE_daily_df)

# Plotting the three stock prices
ggplot(data, aes(x = time)) +
  geom_line(aes(y = XOM.Close, colour = "XOM"), size = 0.5) +
  geom_line(aes(y = CVX.Close, colour = "CVX"), size = 0.5) +
  geom_line(aes(y = NEE.Close, colour = "NEE"), size = 0.5) +
  scale_colour_manual("", 
                      values = c("XOM" = "#F01523", 
                                 "CVX" = "#0054A4", 
                                 "NEE" = "#78C239")) +
  labs(title = "Daily closing prices of the three stocks",
       subtitle = "2016.01.04-2023.05.12",
       caption = "Source: Yahoo Finance",
       y = "Price (USD)",
       x = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, hjust = 0.5 ,face = "bold"),
        plot.subtitle = element_text(size=10, hjust=0.5),
        axis.title.y = element_text(size=10, hjust = 0.5, face='bold'),
        axis.text.x = element_text(angle=45, hjust = 1))

plot(NEE_returns)

# Plot the log returns
# Convert XTS objects to data frames and add a new column for the stock name
XOM_df <- data.frame(time = index(XOM_returns), XOM_value = coredata(XOM_returns))
CVX_df <- data.frame(CVX_value = coredata(CVX_returns))
NEE_df <- data.frame(NEE_value = coredata(NEE_returns))

data <- cbind(XOM_df, CVX_df, NEE_df)

# Plot each data frame
XOM_plot <- ggplot(data, aes(x = time, y = XOM.Close)) +
  geom_line(colour = "#F01523", size = 0.3) +
  labs(x = "", y = "XOM Log Returns") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 10, hjust = 0.5, face = 'bold'))

CVX_plot <- ggplot(data, aes(x = time, y = CVX.Close)) +
  geom_line(colour = "#0054A4", size = 0.3) +
  labs(x = "", y = "CVX Log Returns") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 10, hjust = 0.5, face = 'bold'))

NEE_plot <- ggplot(data, aes(x = time, y = NEE.Close)) +
  geom_line(colour = "#78C239", size = 0.3) +
  labs(x = "Time", y = "NEE Log Returns", caption="Source: Yahoo Finance") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 10, hjust = 0.5, face = 'bold'))

# Arrange plots
combined_plot <- ggarrange(XOM_plot, CVX_plot, NEE_plot, nrow = 3)

# Add title
title <- expression(atop(bold("Logarithmic returns for XOM, CVX, and NEE"), scriptstyle("2016.01.05 - 2023.05.12")))
annotate_figure(combined_plot, top = text_grob(title, color = "black", face = "bold", size = 12))


# Plotting the 3 GARCH volatilities
# Create the dataframe

volatility_df <- as.data.frame(sigma(XOM_GJR)) 
volatility_df2 <- as.data.frame(sigma(CVX_EGARCH)) 
volatility_df3 <- as.data.frame(sigma(NEE_EGARCH)) 

VOL_FINAL <- cbind(volatility_df, volatility_df2, volatility_df3)

# Reset the index and add it as a column
VOL_FINAL$Index <- rownames(VOL_FINAL)

# Rename the columns in VOL_FINAL
colnames(VOL_FINAL) <- c("XOM_VOL", "CVX_VOL", "NEE_VOL", "time")
VOL_FINAL$time <- as.Date(VOL_FINAL$time)

XOM_VOL_PLOT <- ggplot(data=VOL_FINAL, aes(time,XOM_VOL)) +
  geom_line(aes(y=XOM_VOL), colour="#F01523", size=0.3)+
  ylab("XOM Volatility")+xlab("")+
  scale_x_date(date_breaks = "6 month" ,date_labels = "%Y %b")+
  theme_bw()+
  theme(plot.title = element_text(size = 10, hjust = 0.5 ,face = "bold"),
        plot.subtitle = element_text(size=8, hjust=0.5),
        axis.title.y = element_text(size=8, hjust = 0.5, face='bold'),
        axis.text.x = element_text(angle=45, hjust = 1))

CVX_VOL_PLOT <- ggplot(data=VOL_FINAL, aes(time,XOM_VOL)) +
  geom_line(aes(y=XOM_VOL), colour="#0054A4", size=0.3)+
  ylab("CVX Volatility")+xlab("")+
  scale_x_date(date_breaks = "6 month" ,date_labels = "%Y %b")+
  theme_bw()+
  theme(plot.title = element_text(size = 10, hjust = 0.5 ,face = "bold"),
        plot.subtitle = element_text(size=8, hjust=0.5),
        axis.title.y = element_text(size=8, hjust = 0.5, face='bold'),
        axis.text.x = element_text(angle=45, hjust = 1))

NEE_VOL_PLOT <- ggplot(data=VOL_FINAL, aes(time,NEE_VOL)) +
  geom_line(aes(y=XOM_VOL), colour="#78C239", size=0.3)+
  ylab("NEE Volatility")+xlab("")+
  scale_x_date(date_breaks = "6 month" ,date_labels = "%Y %b")+
  theme_bw()+
  theme(plot.title = element_text(size = 10, hjust = 0.5 ,face = "bold"),
        plot.subtitle = element_text(size=8, hjust=0.5),
        axis.title.y = element_text(size=8, hjust = 0.5, face='bold'),
        axis.text.x = element_text(angle=45, hjust = 1))

threevol <- ggarrange(XOM_VOL_PLOT, CVX_VOL_PLOT, NEE_VOL_PLOT, nrow=3)
title <- expression(atop(bold("Volatility of XOM, CVX and NEE"), scriptstyle("2016.01.05 - 2023.05.12")))
annotate_figure(threevol, top=text_grob(title, color="black", face="bold", size = 14))  
