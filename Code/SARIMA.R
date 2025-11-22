####################################################################
# SARIMA models 
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
library(forecast)     # Version 8.24.0    
library(TSA)          # Version 1.3.1
library(here)         # Version 1.0.1
####################################################################

####################################################################
# Functions
source(here('plot_ts.R'))
####################################################################

####################################################################
# SARIMA (1,0,1)x(0,1,0)_12
n=200
h=n/4

set.seed(31415)
Xt = sarima.sim(d=0,ar=.3,ma=.2,D=1,S=12,n=n)
plot_TS(Xt)

# ACF
acf_Xt = acf(Xt,lag.max=h,plot=F)
lag = acf_Xt$lag[,1,]
val = acf_Xt$acf[,1,]
u = rep(1.96/sqrt(n),h)
l = -u
acf_df = data.frame(lag,val,u,l)  
plot_acf(acf_df)

# PACF
pacf_Xt = pacf(Xt,lag.max = h,plot=F)
lag = pacf_Xt$lag[,1,]
val = pacf_Xt$acf[,1,]
u = rep(1.96/sqrt(n),h)
l = -u
pacf_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_df)

# Differentiate the process 
Yt = diff(Xt,12)
plot_TS(Yt)
n_Yt=length(Yt)
h_Yt=n_Yt/4

# ACF
acf_Yt = acf(Yt,lag.max=h_Yt,plot=F)
lag = acf_Yt$lag[,1,]
val = acf_Yt$acf[,1,]
u = rep(1.96/sqrt(n_Yt),h_Yt)
l = -u
acf_df = data.frame(lag,val,u,l)  
plot_acf(acf_df)

# PACF
pacf_Yt = pacf(Yt,lag.max = h_Yt,plot=F)
lag = pacf_Yt$lag[,1,]
val = pacf_Yt$acf[,1,]
u = rep(1.96/sqrt(n_Yt),h_Yt)
l = -u
pacf_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_df)

# ARMA(1,1) for Yt
mod1 = Arima(Yt,order=c(1,0,1));mod1

# Residuals of Yt
df_res=data.frame(x=c(1:n_Yt),y=mod1$residuals)
ggplot(data=df_res,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

#ACF
acf_res = acf(df_res$y,lag.max=h,plot=F)
lag = acf_res$lag[,1,]
val = acf_res$acf[,1,]
u = rep(1.96/sqrt(n_Yt),h)
l = -u
acf_df = data.frame(lag,val,u,l)  
plot_acf(acf_df)

# PACF
pacf_res = pacf(df_res$y,lag.max = h,plot=F)
lag = pacf_res$lag[,1,]
val = pacf_res$acf[,1,]
u = rep(1.96/sqrt(n_Yt),h)
l = -u
pacf_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_df)

# Ljung-Box test
Box.test(df_res$y,lag=12,type='Ljung-Box')

# Prediction for Yt
df_Yt=data.frame(t=c(1:n_Yt),X=Yt)

k=10
preds=forecast(mod1,h=k)
df_preds=data.frame(t=c((n_Yt+1):(n_Yt+k)),x=preds$mean,
                    l=preds$lower[,2],u=preds$upper[,2])

plot_preds(df_Yt,df_preds)

# Model for Xt
mod =Arima(Xt,order = c(1,0,1), seasonal = list(order=c(0,1,0),
                                                  period=12));mod

# Residuals of Xt
df_res=data.frame(x=c(1:n),y=mod$residuals)
ggplot(data=df_res,aes(x=x,y=y))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

#ACF
acf_res = acf(df_res$y,lag.max=h,plot=F)
lag = acf_res$lag[,1,]
val = acf_res$acf[,1,]
u = rep(1.96/sqrt(n),h)
l = -u
acf_df = data.frame(lag,val,u,l)  
plot_acf(acf_df)

# PACF
pacf_res = pacf(df_res$y,lag.max = h,plot=F)
lag = pacf_res$lag[,1,]
val = pacf_res$acf[,1,]
u = rep(1.96/sqrt(n),h)
l = -u
pacf_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_df)

# Ljung-Box test
Box.test(df_res$y,lag=12,type='Ljung-Box')

# Prediction for Xt
df_Xt=data.frame(t=c(1:n),X=Xt)

k=10
preds=forecast(mod,h=k)
df_preds=data.frame(t=c((n+1):(n+k)),x=preds$mean,
                    l=preds$lower[,2],u=preds$upper[,2])

plot_preds(df_Xt,df_preds,CI=F)
plot_preds(df_Xt,df_preds)+
  coord_cartesian(x=c(180,210))
####################################################################
