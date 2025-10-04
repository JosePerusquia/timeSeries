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
####################################################################

####################################################################
# Johnson & Johnson Quarterly Earnings
data(jj)
jj_df = as.data.frame(jj)
jj_df$t = time(jj)

ggplot(data=jj_df,aes(x=t,y=x))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()
####################################################################

####################################################################
# Global warming
data(gtemp_both)
globtemp = as.data.frame(gtemp_both)
globtemp$t = time(gtemp_both)

ggplot(data=globtemp,aes(x=t,y=x))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()
####################################################################

####################################################################
# Speech
data(speech)
speech_df = as.data.frame(speech)
speech_df$t = time(speech)

ggplot(data=speech_df,aes(x=t,y=x))+
  labs(x='',y='')+
  geom_line(col='black')+
  geom_point(size=.5)+
  theme_minimal()
####################################################################

####################################################################
# Financial time series of Apple
getSymbols("AAPL", from = '2024-09-21',
           to = "2025-09-21",warnings = FALSE,
           auto.assign = TRUE)

AAPL_df=as.data.frame(AAPL)
AAPL_df$t=as.Date(time(AAPL))
ggplot(data=AAPL_df,aes(x=t,y=AAPL.Close))+
  geom_line()+
  geom_point(size=.5)+
  theme_minimal()+
  labs(x='',y='')
####################################################################

####################################################################
#El Niño and Fish Population
data(soi)
soi_df=as.data.frame(soi)
soi_df$t=time(soi)

ggplot(data=soi_df,aes(x=t,y=x))+
  geom_line()+
  geom_point(size=.5)+
  theme_minimal()+
  labs(x='',y='')


data(rec)
rec_df=as.data.frame(rec)
rec_df$t=time(rec)

ggplot(data=rec_df,aes(x=t,y=x))+
  geom_line()+
  geom_point(size=.5)+
  theme_minimal()+
  labs(x='',y='')
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
# Basic time series models

# Gaussian white noise
t = c(1:1000)

set.seed(31415)
Zt = rnorm(1000)

GWN = data.frame(t,Zt)
ggplot(data=GWN,aes(x=t,y=Zt))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# acf
acf_zt = acf(Zt,lag.max = 30,plot=F)
lag = acf_zt$lag
val = acf_zt$acf
u = rep(1.96/sqrt(1000),31)
l = -u

acf_zt = data.frame(lag,val,u,l)  

ggplot(data=acf_zt)+
  geom_point(aes(x=lag,y=val),size=.5)+
  geom_segment(x=lag,xend=lag,y=0,yend=val)+
  geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
  geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
  geom_hline(yintercept=0)+
  theme_minimal()+
  labs(x=expression(h),y=expression(ACF))

# Moving average of order 1
set.seed(314159)
ma = arima.sim(n=1000,list(order = c(0,0,1),
                           ma=c(5),sd=sqrt(2)))

t=time(ma)
vals=as.data.frame(ma)               
vals$t=t

ggplot(data=vals,aes(x=t,y=x))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# acf
acf_ma = acf(ma,lag.max = 30,plot=F)

lag = acf_ma$lag
val = acf_ma$acf
u = rep(1.96/sqrt(1000),31)
l = -u

acf_ma = data.frame(lag,val,u,l)  

ggplot(data=acf_ma)+
  geom_point(aes(x=lag,y=val),size=.5)+
  geom_segment(x=lag,xend=lag,y=0,yend=val)+
  geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
  geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
  geom_hline(yintercept=0)+
  theme_minimal()+
  labs(x=expression(h),y=expression(ACF))

# autoregressive of order 1
set.seed(31415)
ar = arima.sim(n=1000,list(order = c(1,0,0),
                           ar=c(.9),sd=sqrt(4)))

t=time(ar)
vals=as.data.frame(ar)               
vals$t=t

ggplot(data=vals,aes(x=t,y=x))+
  geom_line()+
  theme_minimal()+
  labs(x='',y='')

# acf
acf_ar = acf(ar,lag.max = 50,plot=F)

lag = acf_ar$lag
val = acf_ar$acf
u = rep(1.96/sqrt(1000),51)
l = -u

acf_ar = data.frame(lag,val,u,l)  
ggplot(data=acf_ar)+
  geom_point(aes(x=lag,y=val),size=.5)+
  geom_segment(x=lag,xend=lag,y=0,yend=val)+
  geom_line(aes(x=lag,y=u),linetype=2,col='blue')+
  geom_line(aes(x=lag,y=l),linetype=2,col='blue')+
  geom_hline(yintercept=0)+
  theme_minimal()+
  labs(x=expression(h),y=expression(ACF))
####################################################################


