####################################################################
# GARCH Model
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
library(forecast)     # Version 8.24.0    
library(TSA)          # Version 1.3.1
library(rugarch)      # Version 4.3.3
library(zoo)          # Version 1.8-12
####################################################################

####################################################################
# Functions
source(here('plot_ts.R'))
####################################################################

####################################################################
# Financial time series of Apple
getSymbols("AAPL", from = '2022-01-01',to = "2025-11-20",
           warnings = FALSE,auto.assign = TRUE)

t = time(AAPL)
AAPLCl=data.frame(t=t,x=AAPL$AAPL.Close)
names(AAPLCl)=c('t','x')
n=length(AAPLCl$x)

ggplot(data=AAPLCl,aes(x=t,y=x))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# Calculate the log-returns
Zt = log(AAPLCl$x[-1]/AAPLCl$x[-n])
n_Zt = length(Zt)

# Plot the log-returns
Zt=data.frame(t=t[-1],y=Zt)
ggplot(data=Zt,aes(x=t,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# Daily, monthly and annaual volatility
sd_d=sd(Zt$y);sd_d
sd_m=sqrt(22)*sd_d;sd_m
sd_y=sqrt(252)*sd_d;sd_y

# Rolling volatility
rolling_vol = rollapply(Zt$y,width = 21,FUN = sd,fill = NA,
                        align = "right"
)

annual_vol = rolling_vol * sqrt(252)
RollingVol = data.frame(t=t[-1],annual_vol)

ggplot(data=RollingVol,aes(x=t,y=annual_vol))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# ACF of the log-returns
acf_apple = acf(Zt$y,lag.max = 30,plot=F)
lag = acf_apple$lag
val = acf_apple$acf
u = rep(1.96/sqrt(n_Zt),30)
l = -u
acf_apple = data.frame(lag,val,u,l)  
plot_acf(acf_apple)

# ACF of the squared process (there is correlation)
acf_apple = acf(Zt$y^2,lag.max = 30,plot=F)
lag = acf_apple$lag
val = acf_apple$acf
u = rep(1.96/sqrt(n_Zt),30)
l = -u
acf_apple = data.frame(lag,val,u,l)  
plot_acf(acf_apple)
####################################################################

####################################################################
# Fit a GARCH (1,1) model

# Specify the model
mod = ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
    mean.model     = list(armaOrder = c(0,0), include.mean = TRUE),
    distribution.model = "norm"
)

# Fit the model
fit = ugarchfit(spec = mod, data = Zt$y);fit

# Plot the estimated volatility
fit_plot = data.frame(t=t[-1],x=fit@fit$sigma,y=abs(Zt$y))

ggplot(data=fit_plot,aes(x=t,y=x))+
  geom_line(aes(x=t,y=y),alpha=1,col='lightgray')+
  geom_line(col='darkblue')+
  theme_minimal()+
  labs(y=expression(sigma[t]),x='')

# Make a 20 day prediction
k=20
fc = ugarchforecast(fit, n.ahead = k)
sigma_fc = sigma(fc)

# Data frame and plot
sigmaOld = data.frame(t=t[-1],X=fit@fit$sigma)

tnew=t[n]+c(1:k)
sigmaNew = data.frame(t=tnew,x=sigma_fc)
names(sigmaNew)=c('t','x')

plot_preds(sigmaOld[c(950:974),],sigmaNew,CI=F)
####################################################################
