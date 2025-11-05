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
library(forecast)     # Version 8.24.0    
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
# Air passenger data
air = read.csv(here("../Datos/AirPassengers.csv"))
air = air%>%
  rename("x" = X.Passengers)
air$t = c(1:144)

ggplot(data=air,aes(x=t,y=x,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()
####################################################################

####################################################################
# Moving average q=1
air$mt = stats::filter(air$x,rep(1/3,3),sides=2)

ggplot(data=air,aes(x=t,y=x,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_line(aes(x=t,y=mt),col='red')+
  geom_point(size=.5)+
  theme_minimal()

air$yt = air$x-air$mt
ggplot(data=air,aes(x=t,y=yt,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()
####################################################################

####################################################################
# Moving average with q = 6
air$mt = stats::filter(air$x,rep(1/13,13),sides=2)

ggplot(data=air,aes(x=t,y=x,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_line(aes(x=t,y=mt),col='red')+
  geom_point(size=.5)+
  theme_minimal()

air$yt = air$x-air$mt
ggplot(data=air,aes(x=t,y=yt,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()
####################################################################

####################################################################
# Exponential smoothing M,Ad,M
airExp = ets(AirPassengers)
air$exp = airExp$fitted
air$resExp = airExp$residuals

ggplot(data=air,aes(x=t,y=x,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_line(aes(x=t,y=exp),col='red')+
  geom_point(size=.5)+
  theme_minimal()

ggplot(data=air,aes(x=t,y=resExp,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()

# acf
acf_exp = acf(air$resExp,lag.max = 30,plot=F)
lag = acf_exp$lag
val = acf_exp$acf
u = rep(1.96/sqrt(144),31)
l = -u

acf_exp = data.frame(lag,val,u,l)  
plot_acf(acf_exp)

# pacf
pacf_exp = pacf(air$resExp,lag.max = 30,plot=F)
lag = pacf_exp$lag
val = pacf_exp$acf
u = rep(1.96/sqrt(144),30)
l = -u

pacf_exp = data.frame(lag,val,u,l)  
plot_pacf(pacf_exp)

# Ljung-Box test
Box.test(air$resExp,type='Ljung-Box')
####################################################################

####################################################################
# Lag operator 
air$xDif=c(NA,diff(air$x,lag=1))

ggplot(data=air,aes(x=t,y=xDif,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()

# Log on the data and lag operator at 12 and 1
air$log.x = log(air$x)
airDiff12 = diff(air$log.x,lag=12)
airDiff1 = diff(airDiff12,lag=1)
t=c(1:length(airDiff1))

airDiff = data.frame(t,airDiff1)
ggplot(data=airDiff,aes(x=t,y=airDiff1,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()

# acf
acf_diff = acf(airDiff1,lag.max = 30,plot=F)
lag = acf_diff$lag
val = acf_diff$acf
u = rep(1.96/sqrt(131),31)
l = -u

acf_diff = data.frame(lag,val,u,l)  
plot_acf(acf_diff)

# pacf
pacf_diff = pacf(airDiff1,lag.max = 30,plot=F)
lag = pacf_diff$lag
val = pacf_diff$acf
u = rep(1.96/sqrt(131),30)
l = -u

pacf_diff = data.frame(lag,val,u,l)  
plot_pacf(pacf_diff)

# Ljung-Box test
Box.test(airDiff1,type='Ljung-Box')
####################################################################

####################################################################
# ARMA (1,1) for the lagged differences
res1=arima(airDiff$airDiff1,order=c(1,0,1))
summary(res1)

#Residuales
residuals=data.frame(t,y=res1$residuals)
ggplot(data=residuals,aes(x=t,y=y,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()

# acf
acf_diff = acf(residuals$y,lag.max = 30,plot=F)
lag = acf_diff$lag
val = acf_diff$acf
u = rep(1.96/sqrt(131),31)
l = -u

acf_diff = data.frame(lag,val,u,l)  
plot_acf(acf_diff)

# pacf
pacf_diff = pacf(residuals$y,lag.max = 30,plot=F)
lag = pacf_diff$lag
val = pacf_diff$acf
u = rep(1.96/sqrt(131),30)
l = -u

pacf_diff = data.frame(lag,val,u,l)  
plot_pacf(pacf_diff)

# Ljung-Box test
Box.test(residuals$y,type='Ljung-Box')
####################################################################

####################################################################
# ARMA (12,12) for the lagged differences
res1=arima(airDiff$airDiff1,order=c(12,0,12))
summary(res1)

#Residuales
residuals=data.frame(t,y=res1$residuals)
ggplot(data=residuals,aes(x=t,y=y,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()

# acf
acf_diff = acf(residuals$y,lag.max = 30,plot=F)
lag = acf_diff$lag
val = acf_diff$acf
u = rep(1.96/sqrt(131),31)
l = -u

acf_diff = data.frame(lag,val,u,l)  
plot_acf(acf_diff)

# pacf
pacf_diff = pacf(residuals$y,lag.max = 30,plot=F)
lag = pacf_diff$lag
val = pacf_diff$acf
u = rep(1.96/sqrt(131),30)
l = -u

pacf_diff = data.frame(lag,val,u,l)  
plot_pacf(pacf_diff)

# Ljung-Box test
Box.test(residuals$y,type='Ljung-Box')

# Prediction
preds=forecast(res1,h=12)
df_preds=data.frame(t=c(132:143),x=preds$mean,
                    l=preds$lower[,2],u=preds$upper[,2])

ggplot(data=airDiff,aes(x=t,y=airDiff1,group = 1))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  geom_line(data=df_preds,aes(x=t,y=x),col='lightblue')+
  geom_line(data=df_preds,aes(x=t,y=l),col='skyblue4',
            linetype=2)+
  geom_line(data=df_preds,aes(x=t,y=u),col='skyblue4',
            linetype=2)+
  theme_minimal()
####################################################################

