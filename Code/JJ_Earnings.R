####################################################################
# Johnson and Johnson earnings
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
# Johnson & Johnson Quarterly Earnings
data(jj)
n=length(jj)
plot_TS(jj)

# Apply logarithms
jj_log = log(jj)
plot_TS(jj_log)

# Differentiate one lag
jj_diff = diff(jj_log,lag=1)
n_diff=length(jj_diff)
plot_TS(jj_diff)

# Periodogram
period_jj=periodogram(jj_diff,plot=F)
df_period=data.frame(freq=period_jj$freq,spec=period_jj$spec)
plot_periodogram(df_period)

# Dominating cycles
head(order(period_jj$spec,decreasing=T))
1/period_jj$freq[45]
1/period_jj$freq[23]
1/period_jj$freq[44]
1/period_jj$freq[21]

# Differentiate at lag 4
jj_diff = diff(jj_diff,lag=4)
n_diff=length(jj_diff)

# ACF
acf_diff = acf(jj_diff,lag.max = 30,plot=F)
lag = acf_diff$lag
val = acf_diff$acf
u = rep(1.96/sqrt(n_diff),30)
l = -u
acf_diff = data.frame(lag,val,u,l)  
plot_acf(acf_diff)

# PACF
pacf_diff = pacf(jj_diff,lag.max = 30,plot=F)
lag = pacf_diff$lag
val = pacf_diff$acf
u = rep(1.96/sqrt(n_diff),30)
l = -u
pacf_diff = data.frame(lag,val,u,l)  
plot_pacf(pacf_diff) 

# SARIMA
mod = Arima(jj_log,order=c(1,1,1),seasonal=list(order=c(1,1,1),
                                                period=4));mod

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
Box.test(df_res$y,lag=2,type='Ljung-Box')

# Prediction for the log(JJ) series
t=time(jj)
df_Xt=data.frame(t=t,X=jj_log)

k=8
preds=forecast(mod,h=k)
tnew=c(1981.00,1981.25,1981.50,1981.75,1982.00,1982.25,
       1982.50,1982.75)

df_preds=data.frame(t=tnew,x=preds$mean,
                    l=preds$lower[,2],u=preds$upper[,2])

plot_preds(df_Xt,df_preds)

# Prediction for the original series}
df_Xt=data.frame(t=t,X=jj)
df_preds=data.frame(t=tnew,x=exp(preds$mean),
                    l=exp(preds$lower[,2]),u=exp(preds$upper[,2]))

plot_preds(df_Xt,df_preds)
####################################################################
