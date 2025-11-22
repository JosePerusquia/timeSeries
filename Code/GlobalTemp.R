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
library(tseries)      # Version 0.10-58
library(here)         # Version 1.0.1
library(forecast)     # Version 8.24.0    
library(TSA)          # Version 1.3.1
####################################################################

####################################################################
# Functions
source(here('plot_ts.R'))
####################################################################

####################################################################
# Global temperature
data(gtemp_both)
n=length(gtemp_both)
plot_TS(gtemp_both)
####################################################################

####################################################################
# Auto arima function
auto.arima(gtemp_both)
mod = Arima(gtemp_both, order = c(0,1,1));mod

# Residuals
df_res=data.frame(x=c(1:n),y=mod$residuals)
ggplot(data=df_res,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# ACF
acf_res = acf(df_res$y,lag.max=30,plot=F)
lag = acf_res$lag[,1,]
val = acf_res$acf[,1,]
u = rep(1.96/sqrt(n),30)
l = -u
acf_df = data.frame(lag,val,u,l)  
plot_acf(acf_df)

# PACF
pacf_res = pacf(df_res$y,lag.max = 30,plot=F)
lag = pacf_res$lag[,1,]
val = pacf_res$acf[,1,]
u = rep(1.96/sqrt(n),30)
l = -u
pacf_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_df)

# Ljung-Box test
Box.test(df_res$y,lag=12,type='Ljung-Box')

# Prediction for Xt
t=time(gtemp_both)
df_Xt=data.frame(t=t,X=gtemp_both)

k=12
preds=forecast(mod,h=k)
tnew=t[n]+c(1:k)
df_preds=data.frame(t=tnew,x=preds$mean,
                    l=preds$lower[,2],u=preds$upper[,2])

plot_preds(df_Xt,df_preds)
####################################################################

####################################################################
# ARIMA
gtemp_diff = diff(gtemp_both,lag=1)
n_diff=length(gtemp_diff)
plot_TS(gtemp_diff)

# ACF
acf_diff = acf(gtemp_diff,lag.max = 30,plot=F)
lag = acf_diff$lag
val = acf_diff$acf
u = rep(1.96/sqrt(n_diff),30)
l = -u
acf_diff = data.frame(lag,val,u,l)  
plot_acf(acf_diff)

# PACF
pacf_diff = pacf(gtemp_diff,lag.max = 30,plot=F)
lag = pacf_diff$lag
val = pacf_diff$acf
u = rep(1.96/sqrt(n_diff),30)
l = -u
pacf_diff = data.frame(lag,val,u,l)  
plot_pacf(pacf_diff) 

# ARIMA model
mod = Arima(gtemp_both, order = c(1,1,1));mod

# Residuals
df_res=data.frame(x=c(1:n),y=mod$residuals)
ggplot(data=df_res,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# ACF
acf_res = acf(df_res$y,lag.max=30,plot=F)
lag = acf_res$lag[,1,]
val = acf_res$acf[,1,]
u = rep(1.96/sqrt(n),30)
l = -u
acf_df = data.frame(lag,val,u,l)  
plot_acf(acf_df)

# PACF
pacf_res = pacf(df_res$y,lag.max = 30,plot=F)
lag = pacf_res$lag[,1,]
val = pacf_res$acf[,1,]
u = rep(1.96/sqrt(n),30)
l = -u
pacf_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_df)

# Ljung-Box test
Box.test(df_res$y,lag=12,type='Ljung-Box')

# Prediction for Xt
t=time(gtemp_both)
df_Xt=data.frame(t=t,X=gtemp_both)

k=12
preds=forecast(mod,h=k)
tnew=t[n]+c(1:k)
df_preds=data.frame(t=tnew,x=preds$mean,
                    l=preds$lower[,2],u=preds$upper[,2])

plot_preds(df_Xt,df_preds)
####################################################################

####################################################################
# SARIMA

# Differentiate at lag 12
gtemp_diff = diff(gtemp_diff,lag=12)
n_diff=length(gtemp_diff)
plot_TS(gtemp_diff)

# acf
acf_diff = acf(gtemp_diff,lag.max = 30,plot=F)
lag = acf_diff$lag
val = acf_diff$acf
u = rep(1.96/sqrt(n_diff),30)
l = -u
acf_diff = data.frame(lag,val,u,l)  
plot_acf(acf_diff)

# pacf
pacf_diff = pacf(gtemp_diff,lag.max = 30,plot=F)
lag = pacf_diff$lag
val = pacf_diff$acf
u = rep(1.96/sqrt(n_diff),30)
l = -u
pacf_diff = data.frame(lag,val,u,l)  
plot_pacf(pacf_diff) 

# SARIMA
mod = Arima(gtemp_both, order = c(3,1,1), 
            seasonal = list(order=c(2,1,1),period=12));mod

# Residuales
df_res=data.frame(x=c(1:n),y=mod$residuals)
ggplot(data=df_res,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# ACF
acf_res = acf(df_res$y,lag.max=30,plot=F)
lag = acf_res$lag[,1,]
val = acf_res$acf[,1,]
u = rep(1.96/sqrt(n),30)
l = -u
acf_df = data.frame(lag,val,u,l)  
plot_acf(acf_df)

# PACF
pacf_res = pacf(df_res$y,lag.max = 30,plot=F)
lag = pacf_res$lag[,1,]
val = pacf_res$acf[,1,]
u = rep(1.96/sqrt(n),30)
l = -u
pacf_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_df)

# Ljung-Box test
Box.test(df_res$y,lag=12,type='Ljung-Box')

# Prediction for Xt
k=12
preds=forecast(mod,h=k)
df_preds=data.frame(t=tnew,x=preds$mean,
                    l=preds$lower[,2],u=preds$upper[,2])

plot_preds(df_Xt,df_preds)
####################################################################

####################################################################
# Exponential smoothing state space model
gtemp_ETS = ets(gtemp_both);gtemp_ETS

# Residuals
res=gtemp_ETS$residuals

# ACF
acf_exp = acf(res,lag.max = 30,plot=F)
lag = acf_exp$lag
val = acf_exp$acf
u = rep(1.96/sqrt(n),30)
l = -u
acf_exp = data.frame(lag,val,u,l)  
plot_acf(acf_exp)

# PACF
pacf_exp = pacf(res,lag.max = 30,plot=F)
lag = pacf_exp$lag
val = pacf_exp$acf
u = rep(1.96/sqrt(n),30)
l = -u
pacf_exp = data.frame(lag,val,u,l)  
plot_pacf(pacf_exp)

# Ljung-Box test
Box.test(res,lag=12,type='Ljung-Box')

# Prediction for Xt
k=12
preds=forecast(gtemp_ETS,h=k)
df_preds=data.frame(t=tnew,x=preds$mean,
                    l=preds$lower[,2],u=preds$upper[,2])

plot_preds(df_Xt,df_preds)
####################################################################
