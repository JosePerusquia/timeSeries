####################################################################
# Time series examples and basic models
# Author: Jose Antonio Perusquia Cortes
# Afil: Facultad de Ciencias - UNAM
# Module: Time series 
####################################################################

####################################################################
# Libraries
library(astsa)        # Version 2.3
library(ggplot2)      # Version 3.5.2
library(ggthemes)     # Version 5.1.0
library(dplyr)        # Version 1.1.4
library(tidyquant)    # Version 1.0.11
library(tseries)      # Version 0.10-58
library(here)         # Version 1.0.1
library(TSA)          # Version 1.3.1
####################################################################

####################################################################
# Functions
source(here('plot_ts.R'))
####################################################################

####################################################################
# Johnson & Johnson Quarterly Earnings
data(jj)
plot_TS(jj)
####################################################################

####################################################################
# Global warming
data(gtemp_both)
plot_TS(gtemp_both)
####################################################################

####################################################################
# Speech
data(speech)
plot_TS(speech)
####################################################################

####################################################################
# Financial time series of Apple
getSymbols("AAPL", from = '2024-09-21',to = "2025-09-21",
           warnings = FALSE,auto.assign = TRUE)
AAPL=as.ts(AAPL$AAPL.Close)
plot_TS(AAPL)
####################################################################

####################################################################
#El Niño and Fish Population
data(soi)
plot_TS(soi)

data(rec)
plot_TS(rec)
####################################################################

####################################################################
# Magnetic resonance
data(fmri1)
fmri1_df=as.data.frame(fmri1)

ggplot(data=fmri1_df)+
  geom_line(aes(x=time,y=cort1),col='darkred')+
  geom_line(aes(x=time,y=cort2),col='darkblue')+
  geom_line(aes(x=time,y=cort3),col='darkgreen')+
  geom_line(aes(x=time,y=cort4),col='darkorange')+
  theme_minimal()+
  labs(x='',y='')
####################################################################

####################################################################
# Gaussian white noise
n = 1000
set.seed(3141)
Zt = as.ts(rnorm(n))
plot_TS(Zt)

# acf
acf_zt = acf(Zt,lag.max = 30,plot=F)
lag = acf_zt$lag
val = acf_zt$acf
u = rep(1.96/sqrt(n),30)
l = -u

acf_zt = data.frame(lag,val,u,l)  
plot_acf(acf_zt)

# pacf
pacf_zt = pacf(Zt,lag.max = 30,plot=F)
lag = pacf_zt$lag
val = pacf_zt$acf
u = rep(1.96/sqrt(n),30)
l = -u

pacf_zt = data.frame(lag,val,u,l)  
plot_pacf(pacf_zt)
####################################################################

####################################################################
# MA(1)
set.seed(314159)
ma = arima.sim(n=n,list(order = c(0,0,1),ma=c(.5),sd=4))
plot_TS(ma)

# ACF
acf_ma = acf(ma,lag.max = 30,plot=F)
lag = acf_ma$lag
val = acf_ma$acf
u = rep(1.96/sqrt(n),30)
l = -u

acf_ma = data.frame(lag,val,u,l)  
plot_acf(acf_ma)

# PACF
pacf_ma = pacf(ma,lag.max = 30,plot=F)
lag = pacf_ma$lag
val = pacf_ma$acf
u = rep(1.96/sqrt(n),30)
l = -u

pacf_ma = data.frame(lag,val,u,l)  
plot_pacf(pacf_ma)
####################################################################

####################################################################
# AR(1)
set.seed(31415)
ar = arima.sim(n=n,list(order = c(1,0,0),ar=c(.9),sd=4))
plot_TS(ar)

# ACF
acf_ar = acf(ar,lag.max = 30,plot=F)
lag = acf_ar$lag
val = acf_ar$acf
u = rep(1.96/sqrt(n),30)
l = -u
acf_ar = data.frame(lag,val,u,l)  

plot_acf(acf_ar)

# PACF
pacf_ar = pacf(ar,lag.max = 30,plot=F)
lag = pacf_ar$lag
val = pacf_ar$acf
u = rep(1.96/sqrt(n),30)
l = -u
pacf_ar = data.frame(lag,val,u,l)  

plot_pacf(pacf_ar)
####################################################################