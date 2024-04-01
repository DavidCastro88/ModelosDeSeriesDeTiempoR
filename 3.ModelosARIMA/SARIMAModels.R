#LIBRERIAS Y PAQUETES ----------------------------------------------------------
rm(list=ls(all=TRUE))
library(forecast)
library(FitAR)
library(car)
library(carData)
library(TSA)
library(fANCOVA)
library("ggplot2")
library("showtext")
library(pdR)
library(timsac);library(lmtest) 
library(uroot)
windowsFonts(A = windowsFont("Times New Roman"))
source(file.choose()) 
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funcion-regexponencial.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/Funciones-BP.LB.test-pruebaDW1.R")
source("https://raw.githubusercontent.com/NelfiGonzalez/Funciones-de-Usuario-Estadistica-III/main/
Funciones-Criterios.Informacion-Calidad.Intervalos.R")

#Leer los datos
Datos16=read.table(file.choose(),header=T,sep=";",skip=9,dec=",",colClasses=c(rep("NULL",8),"numeric",rep("NULL",10)))
serie=ts(Datos16,freq=12,start=c(2003,1))
#GRAFICAS DE LA SERIE
win.graph(width=5,height=5,pointsize=8) 
plot(serie,xlab = "Time",ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
#Serie recortada y ACF
win.graph(width=5,height=5,pointsize=8) 
plot(yt,xlab = "Time",ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
win.graph(width=5,height=5,pointsize=8) 
acf(as.numeric(yt),main="",lag.max=36,ci.type="ma",ci.col=2,lwd=2,family = "A",cex.axis = 1)
title(main="ACF Yt",cex.main=1)
#Diferencia regular tendencia
difReg=diff(yt)
win.graph(width=5,height=5,pointsize=8) 
plot(difReg,xlab = "Time",ylab='???Yt',family = "A",cex.lab=20,cex.axis = 0.5)
abline(h=mean(difReg),lwd=1)


win.graph(width=5,height=5,pointsize=8) 
acf(as.numeric(difReg),main="",lag.max=36,ci.type="ma",ci.col=2,lwd=2,family = "A",cex.axis = 1)
title(main="ACF diferencia regular",cex.main=1)
abline(v=c(12,24,36),lty=2,col=2)
#Diferencia estacional
difEst=diff(yt,lag=12)
win.graph(width=5,height=5,pointsize=8) 
plot(difEst,xlab = "Time",ylab='???12Yt',family = "A",cex.lab=20,cex.axis = 0.5)
abline(h=mean(difEst),lwd=1)

win.graph(width=5,height=5,pointsize=8) 
acf(as.numeric(difEst),main="",lag.max=36,ci.type="ma",ci.col=2,lwd=2,family = "A",cex.axis = 1)
title(main="ACF diferencia estacional",cex.main=1)
abline(v=c(12,24,36),lty=2,col=2)
#Diferencia Regular y Estacional
difRegEst=diff(diff(yt,12),1)
win.graph(width=5,height=5,pointsize=8) 
plot(difRegEst,xlab = "Time",ylab='??????12Yt',family = "A",cex.lab=20,cex.axis = 0.5)
abline(h=mean(difRegEst),lwd=1)
win.graph(width=5,height=5,pointsize=8)
acf(as.numeric(difRegEst),main="",lag.max=36,ci.type="ma",ci.col=2,lwd=2,family = "A",cex.axis = 1)
title(main="ACF diferencia regular y estacional (d=D=1)",cex.main=1)
abline(v=c(12,24,36),lty=2,col=2)
win.graph(width=3.7,height=2.8)
pacf(as.numeric(difRegEst),lag.max=36,lwd=3,main="",cex.lab=0.5,cex.axis=0.5,family = "A")
title(main="PACF diferencia regular y estacional (d=D=1)",cex.main=0.5)
abline(v=c(12,24,36),lty=2,col=2)


#Test Hegy sobre serie recortada
HEGY.test(wts=yt,itsd=c(0,0,c(0)),selectlags=list(mode="aic", Pmax=12))$stats
n=length(serie)-12
t=1:n


#Serie recortada
yt=ts(serie[t],freq=12,start=c(2003,1))
#Definiendo variables para pron?sticos
tnuevo=(n+1):length(serie) #?ndice de tiempo en los pron?sticos
ytnuevo=ts(serie[tnuevo],freq=12,start=c(2018,7))

#Auto.arima()
auto.arima(yt,ic="aic",seasonal.test="ocsb")
auto.arima(yt,ic="aic",seasonal.test="ch")
auto.arima(yt,ic="aic",seasonal.test="seas")
auto.arima(yt,ic="bic",seasonal.test="ocsb")
auto.arima(yt,ic="bic",seasonal.test="ch")
auto.arima(yt,ic="bic",seasonal.test="seas")

#armasubsets()
win.graph(width=9,height=5,pointsize=8)
plot(armasubsets(difRegEst,nar=24,nma=24,y.name='Y',ar.method='ols'))


#-----MODELO 1:ARIMA(4,1,0)(0,1,1)[12]------------------------------------------------------
modelo1=Arima(yt,order=c(4,1,0),seasonal=list(order=c(0,1,1)),method="ML")
k1=length(modelo1$coef[modelo1$coef!=0]);k1 #Calcular k el total de par?ametros del modelo
coeftest(modelo1)
summary(modelo1)

#-----MODELO 2 ARIMA(0,1,1)(0,1,1)[12] ------------------------------------------------------
modelo2=Arima(yt,order=c(0,1,1),seasonal=list(order=c(0,1,1)),method="ML")
k2=length(modelo2$coef[modelo2$coef!=0]);k2 
coeftest(modelo2)
summary(modelo2)

#-----MODELO 3 armasubsets() ARIMA (4,1,1)(0,1,1)[12] ------------------------------------------------------
modelo3=Arima(yt,order=c(4,1,1),seasonal=list(order=c(0,1,1)),fixed=c(rep(0,3),NA,NA,NA),method="ML")
k3=length(modelo3$coef[modelo3$coef!=0]);k3
coeftest(modelo3)
summary(modelo3)

#-----MODELO 4 armasubsets() ARIMA (4,1,11)(0,1,1)[12] ------------------------------------------------------
modelo4=Arima(yt,order=c(4,1,11),seasonal=list(order=c(0,1,1)),fixed=c(rep(0,3),NA,NA,rep(0,5),NA,rep(0,3),NA,NA),method="ML")
k4=length(modelo4$coef[modelo4$coef!=0]);k4
coeftest(modelo4)
summary(modelo4)

#------------------------------------------------------------------
#------------------------------------------------------------------
#INDICADORES-PRUEBAS
#AIC-BIC
Criteriosmodelo1=exp.crit.inf.resid(residuales=residuals(modelo1),n.par=k1); Criteriosmodelo1
Criteriosmodelo2=exp.crit.inf.resid(residuales=residuals(modelo2),n.par=k2); Criteriosmodelo2 
Criteriosmodelo3=exp.crit.inf.resid(residuales=residuals(modelo3),n.par=k3); Criteriosmodelo3
Criteriosmodelo4=exp.crit.inf.resid(residuales=residuals(modelo4),n.par=k4); Criteriosmodelo4
#Ljun-Box
BP.LB.test(residuals(modelo1),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo2),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo3),maxlag=36,type="Ljung")
BP.LB.test(residuals(modelo4),maxlag=36,type="Ljung")
#Normalidad
shapiro.test(residuals(modelo1))
shapiro.test(residuals(modelo2))
shapiro.test(residuals(modelo3))
shapiro.test(residuals(modelo4))
#Pronostico
predmod1=ts(as.data.frame(forecast(modelo1,h=12,level=95)),freq=12,start=c(2018,7));predmod1
ytpron1=predmod1[,1] #Pron?ostico puntual
predmod2=ts(as.data.frame(forecast(modelo2,h=12,level=95)),freq=12,start=c(2018,7));predmod2
ytpron2=predmod2[,1] #Pron?ostico puntual
predmod3=ts(as.data.frame(forecast(modelo3,h=12,level=95)),freq=12,start=c(2018,7));predmod3
ytpron3=predmod3[,1] #Pron?ostico puntual
predmod4=ts(as.data.frame(forecast(modelo4,h=12,level=95)),freq=12,start=c(2018,7));predmod4
ytpron4=predmod4[,1] #Pron?ostico puntual
#Medidas de pronostico
accuracy(ytpron1,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmod1[,2],LSP=predmod1[,3]) 
accuracy(ytpron2,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmod2[,2],LSP=predmod2[,3]) 
accuracy(ytpron3,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmod3[,2],LSP=predmod3[,3])
accuracy(ytpron4,ytnuevo)
amplitud.cobertura(real=ytnuevo,LIP=predmod4[,2],LSP=predmod4[,3]) 
#------------------------------------------------------------------
#------------------------------------------------------------------
#GRAFICAS
#Ajustes
#Modelo 1
win.graph(width=5,height=5,pointsize=8)
ythat1=modelo1$fitted #Este objeto ya tiene fechas
plot(serie,ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythat1,col=2)
legend("topleft",legend=c("Datos","Ajuste Modelo 1"),lty=1,col=1:2,lwd=2, cex=0.8)
#Modelo 2
win.graph(width=5,height=5,pointsize=8)
ythat2=modelo2$fitted #Este objeto ya tiene fechas
plot(serie,ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythat2,col=2)
legend("topleft",legend=c("Datos","Ajuste Modelo 2"),lty=1,col=1:2,lwd=2, cex=0.8)
#Modelo 3
win.graph(width=5,height=5,pointsize=8)
ythat3=modelo3$fitted #este objeto ya queda con las fechas de la serie
plot(serie,ylab='Serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythat3,col=2)
legend("topleft",legend=c("Datos","Ajuste Modelo 3"),lty=1,col=1:2,lwd=2, cex=0.8)
#Modelo 4
win.graph(width=5,height=5,pointsize=8)
ythat4=modelo4$fitted #este objeto ya queda con las fechas de la serie
plot(serie,ylab='serie',family = "A",cex.lab=20,cex.axis = 0.5)
lines(ythat4,col=2)
legend("topleft",legend=c("Datos","Ajuste modelo 4"),lty=1,col=1:2,lwd=2, cex=0.8)
#Graficas de ACF
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(modelo1)),ci.type="ma",lag.max=36,main="ACF Modelo 1",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(modelo2)),ci.type="ma",lag.max=36,main="ACF Modelo 2",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(modelo3)),ci.type="ma",lag.max=36,main="ACF Modelo 3",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
acf(as.numeric(residuals(modelo4)),ci.type="ma",lag.max=36,main="ACF Modelo 4",ci.col=2,family = "A")
#Graficas de PACF 
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(modelo1)),lag.max=36,main="PACF Modelo 1",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(modelo2)),lag.max=36,main="PACF Modelo 2",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(modelo3)),lag.max=36,main="PACF Modelo 3",ci.col=2,family="A")
win.graph(width=5,height=5,pointsize=10)
pacf(as.numeric(residuals(modelo4)),lag.max=36,main="PACF Modelo 4",ci.col=2,family = "A")
#Graficas de Normalidad
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(modelo1),main="Gr?fico de normalidad residuos de ajuste Modelo 1",family="A")
qqline(residuals(modelo1),col=2) 
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(modelo2),main="Gr?fico de normalidad residuos de ajuste Modelo 2",family="A")
qqline(residuals(modelo2),col=2) 
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(modelo3),main="Gr?fico de normalidad residuos de ajuste Modelo 3",family = "A")
qqline(residuals(modelo3),col=2)
win.graph(width=5,height=5,pointsize=10)
qqnorm(residuals(modelo4),main="Gr?fico de normalidad residuos de ajuste Modelo 4",family = "A")
qqline(residuals(modelo4),col=2) 
#Graficas de Residuos vs Tiempo
win.graph(width=5,height=5,pointsize=10)
plot(residuals(modelo1),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo1$sigma2),2*sqrt(modelo1$sigma2)),lty=2)
win.graph(width=5,height=5,pointsize=10)
plot(residuals(modelo2),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo2$sigma2),2*sqrt(modelo2$sigma2)),lty=2)
win.graph(width=5,height=5,pointsize=10)
plot(residuals(modelo3),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo3$sigma2),2*sqrt(modelo3$sigma2)),lty=2)
win.graph(width=5,height=5,pointsize=10)
plot(residuals(modelo4),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo4$sigma2),2*sqrt(modelo4$sigma2)),lty=2)
#Graficas de Residuos vs valores ajustados 
win.graph(width=5,height=5,pointsize=10)
plot(as.numeric(modelo1$fitted),residuals(modelo1),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo1$sigma2),2*sqrt(modelo1$sigma2)),lty=2) 
win.graph(width=5,height=5,pointsize=10)
plot(as.numeric(modelo2$fitted),residuals(modelo2),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo2$sigma2),2*sqrt(modelo2$sigma2)),lty=2)
win.graph(width=5,height=5,pointsize=10)
plot(as.numeric(modelo3$fitted),residuals(modelo3),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo3$sigma2),2*sqrt(modelo3$sigma2)),lty=2) 
win.graph(width=5,height=5,pointsize=10)
plot(as.numeric(modelo4$fitted),residuals(modelo4),family = "A");abline(h=0)
abline(h=c(-2*sqrt(modelo4$sigma2),2*sqrt(modelo4$sigma2)),lty=2)
#Comparaci?on gr?afica de los pron?osticos y la serie real en periodos ex-post
win.graph(width=5,height=5,pointsize=10)
p=cbind(ytnuevo,ytpron1,ytpron2,ytpron3,ytpron4)
plot(main="Valores reales y pron?sticos ?ltimos 12 meses",family="A",p,plot.type="single",type="b",pch=c(19,1:4),col=c("black","blue","orange","darkgreen","darkred"),lty=1,lwd=2,ylab="Ytf",xaxt="n")
axis(1,at=time(p),labels=c("Jul-18","Ago-18","Sep-18","Oct-18","Nov-18","Dic-18","Ene-19","Feb-19","Mar-19","Abr-19","May-19","Jun-19"),family="A")
legend("topleft",legend=c("Serie Original","Pron?stico Modelo 1","Pron?stico Modelo 2","Pron?stico Modelo 3","Pron?stico Modelo 4"),
       pch=c(19,1:4),lty=c(1:4),col=c("black","blue","orange","darkgreen","darkred"),lwd=2,cex=0.8)