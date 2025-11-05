####################################################################
# ARMA models 
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
####################################################################

####################################################################
# Functions
plot_ARMA = function(mod){
  t=time(mod)
  vals=as.data.frame(mod)               
  vals$t=t
  
  p=ggplot(data=vals,aes(x=t,y=x))+
    geom_line()+
    theme_minimal()+
    labs(x='',y='')
  return(p)
}

plot_acf =function(df_acf){
  p=ggplot(data=df_acf)+
    geom_point(aes(x=lag,y=val),size=1)+
    geom_segment(x=lag,xend=lag,y=0,yend=val)+
    geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
    geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
    geom_hline(yintercept=0)+
    theme_minimal()+
    labs(x=expression(h),y=expression(ACF))
  return(p)
}

plot_pacf =function(df_pacf){
  p=ggplot(data=df_pacf)+
    geom_point(aes(x=lag,y=val),size=1)+
    geom_segment(x=lag,xend=lag,y=0,yend=val)+
    geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
    geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
    geom_hline(yintercept=0)+
    theme_minimal()+
    labs(x=expression(h),y=expression(PACF))
  return(p)
}
####################################################################

####################################################################
# MA(1)
n=100
h=n/4

set.seed(314159)
ma = arima.sim(n=n,list(order = c(0,0,1),ma=c(.9)),sd=sqrt(5))
plot_ARMA(ma)

# ACF with confidence bands of white noise
acf_ma = acf(ma,lag.max = h,plot=F)
lag = acf_ma$lag[,1,]
val = acf_ma$acf[,1,]

u = rep(1.96/sqrt(n),h+1)
l = -u

acf_ma_df = data.frame(lag,val,u,l)  
plot_acf(acf_ma_df)

# ACF with confidence bands of moving average
w = numeric(h+1)
w[1]=NA
w[2]=1-3*(val[2]^2)+4*(val[2]^4)
w[3:(h+1)]=1+2*(val[2]^2)

u = val+rep(1.96*sqrt(w)/sqrt(n))
l = val-rep(1.96*sqrt(w)/sqrt(n))

acf_ma_df = data.frame(lag,val,u,l)  
plot_acf(acf_ma_df)

# PACF with confidence bands of white noise
pacf_ma = pacf(ma,lag.max = h,plot=F)
lag = pacf_ma$lag[,1,]
val = pacf_ma$acf[,1,]

u = rep(1.96/sqrt(n),h)
l = -u

pacf_ma_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_ma_df)
####################################################################

####################################################################
# AR(1)
set.seed(314159)
ar = arima.sim(n=n,list(order = c(1,0,0),ar=c(.9)),sd=sqrt(10))
plot_ARMA(ar)

# ACF with white noise confidence bands
acf_ar = acf(ar,lag.max = h,plot=F)
lag = acf_ar$lag[,1,]
val = acf_ar$acf[,1,]
u = rep(1.96/sqrt(n),h+1)
l = -u

acf_ar_df = data.frame(lag,val,u,l)  
plot_acf(acf_ar_df)

# ACF with true autoregressive confidence bands
phi=.9
w = numeric(h+1)
w[1]=NA
for(i in 1:h){
  aux1=(1-(phi^(2*i)))
  aux2=(1+(phi^2))
  aux3=(1-(phi^2))^(-1)
  aux4=2*i*(phi^(2*i))
  w[i+1]=aux1*aux2*aux3-aux4
}

u = val+rep(1.96*sqrt(w)/sqrt(n))
l = val-rep(1.96*sqrt(w)/sqrt(n))

acf_ar_df = data.frame(lag,val,u,l)  
plot_acf(acf_ar_df)

# PACF with confidence bands of white noise
pacf_ar = pacf(ar,lag.max = h,plot=F)
lag = pacf_ar$lag[,1,]
val = pacf_ar$acf[,1,]

u = rep(1.96/sqrt(n),h)
l = -u

pacf_ar_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_ar_df)
####################################################################

####################################################################
# MA(3)
set.seed(314159)
ma = arima.sim(n=n,list(order = c(0,0,3),ma=c(.9,1.5,2)),sd=sqrt(5))
plot_ARMA(ma)

# ACF with confidence bands of white noise
acf_ma = acf(ma,lag.max = h,plot=F)
lag = acf_ma$lag[,1,]
val = acf_ma$acf[,1,]

u = rep(1.96/sqrt(n),h+1)
l = -u

acf_ma_df = data.frame(lag,val,u,l)  
plot_acf(acf_ma_df)

# PACF with confidence bands of white noise
pacf_ma = pacf(ma,lag.max = h,plot=F)
lag = pacf_ma$lag[,1,]
val = pacf_ma$acf[,1,]

u = rep(1.96/sqrt(n),h)
l = -u

pacf_ma_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_ma_df)
####################################################################

####################################################################
# AR(2)
set.seed(314159)
ar = arima.sim(n=n,list(order = c(2,0,0),ar=c(-.7,-.5)),sd=sqrt(10))
plot_ARMA(ar)

# ACF with white noise confidence bands
acf_ar = acf(ar,lag.max = h,plot=F)
lag = acf_ar$lag[,1,]
val = acf_ar$acf[,1,]
u = rep(1.96/sqrt(n),h+1)
l = -u

acf_ar_df = data.frame(lag,val,u,l)  
plot_acf(acf_ar_df)

# PACF with confidence bands of white noise
pacf_ar = pacf(ar,lag.max = h,plot=F)
lag = pacf_ar$lag[,1,]
val = pacf_ar$acf[,1,]

u = rep(1.96/sqrt(n),h)
l = -u

pacf_ar_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_ar_df)
####################################################################

####################################################################
# ARMA (1,1)
set.seed(3141)
Xt = arima.sim(n=n,list(order = c(1,0,1),ar=c(.5),ma=c(.4)),
               sd=sqrt(5))
plot_ARMA(Xt)

# ACF with white noise confidence bands
acf_arma = acf(Xt,lag.max = h,plot=F)
lag = acf_arma$lag[,1,]
val = acf_arma$acf[,1,]
u = rep(1.96/sqrt(n),h+1)
l = -u

acf_arma_df = data.frame(lag,val,u,l)  
plot_acf(acf_arma_df)

# PACF with confidence bands of white noise
pacf_arma = pacf(Xt,lag.max = h,plot=F)
lag = pacf_arma$lag[,1,]
val = pacf_arma$acf[,1,]

u = rep(1.96/sqrt(n),h)
l = -u

pacf_arma_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_arma_df)
####################################################################

####################################################################
# ARMA (2,1)
set.seed(3141)
Xt = arima.sim(n=n,list(order = c(2,0,1),ar=c(.75,-.5625),ma=c(1.25)),
               sd=sqrt(5))
plot_ARMA(Xt)

# ACF with white noise confidence bands
acf_arma = acf(Xt,lag.max = h,plot=F)
lag = acf_arma$lag[,1,]
val = acf_arma$acf[,1,]
u = rep(1.96/sqrt(n),h+1)
l = -u

acf_arma_df = data.frame(lag,val,u,l)  
plot_acf(acf_arma_df)

# PACF with confidence bands of white noise
pacf_arma = pacf(Xt,lag.max = h,plot=F)
lag = pacf_arma$lag[,1,]
val = pacf_arma$acf[,1,]

u = rep(1.96/sqrt(n),h)
l = -u

pacf_arma_df = data.frame(lag,val,u,l)  
plot_pacf(pacf_arma_df)
####################################################################
